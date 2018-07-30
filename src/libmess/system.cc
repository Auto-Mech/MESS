#include <iostream>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <cstdarg>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <dlfcn.h>

#include <string>
#include <vector>

#include "system.hh"
#include "error.hh"
#include "array.hh"

//#include "tmatrix.hh"

// remove files from directory

#include <ftw.h>

int rm_file (const char* fname, const struct stat* fstat, int flag)
{
  std::remove(fname);
  return 0;
}

void System::clean_dir (const char* dname)
{
    ftw(dname, rm_file, 100); // 100 descriptors is allowed
}

// copy files
int System::file_copy(const char* old_fname, const char* new_fname)
{
    std::ifstream from(old_fname);
  if(!from)
    return 1;
  std::ofstream   to(new_fname);
  if (!to)
    return 2;

  char ch;
  while(from.get(ch))
    if (!to.put(ch))
      return -1;
  return 0;
}

/************************************************************************
 **************         External call      ******************************
 ************************************************************************/

// list of arguments should be terminated by 0 pointer
int System::call_exe (const char* exename ...) 
{
    const char funame [] = "System::call_exe: ";

    char* argv [7]; //maximum number of parameters = 5
    Array<char> xn(int(strlen(exename) + 1));

    va_list ap;
    va_start(ap, exename);
    char** argp = argv;
    *argp++ = basename(strcpy(xn, exename));
    while (*argp++ = va_arg(ap, char*)) {}
    va_end(ap);

    std::cout.flush();
    int fd, exe_stat;
    switch(fork()) {      
	case -1: //error

	    std::cerr << funame << "fork failed\n";
	    return 1000;

	case 0: //child

	    // redirect standard output
	    fd = open("std.out", O_WRONLY | O_CREAT | O_TRUNC, 
		      S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
	    dup2(fd, 1);
	    close(fd);
	    // redirect standard error
	    fd = open("std.err", O_WRONLY | O_CREAT | O_TRUNC,
		      S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );
	    dup2(fd, 2);
	    close(fd);
	    // substitute another image
	    exe_stat = execvp(exename, argv);
	    std::cerr << funame << "running <" << exename << "> with arguments ";
	    argp = argv;
	    while(*argp)
		std::cerr << "<" << *argp++ << ">";
	    std::cerr << "failed with the exit status " << exe_stat << std::endl;
	    exit(1);

	default: //parent

	    int child_stat;
	    wait(&child_stat);
	    return WEXITSTATUS(child_stat);

    }
}

/*************************************************************************
 *                      Pipes
 *************************************************************************/

System::Pipe_base::Pipe_base() throw(Error::General)
{
    const char funame [] = "System::Pipe_base::Pipe_base: ";

    if (pipe(fd)) //create a pipe
	switch(errno) {
	    case EMFILE:
		std::cerr << "pipe: too many file descriptors are in use\n";
		throw Error::Init();
	    case ENFILE:
		std::cerr << funame << "system file table is full\n";
		throw Error::Init();
	    case EFAULT:
		std::cerr << funame << "file descriptor is not valid\n";
		throw Error::Init();
	    default:
		std::cerr << funame << "unknown error\n";
		throw Error::Init();
	}
}

// it seems that a new gcc standard does not support attaching streams to file descriptors
//Pipe::Pipe () : pin(fd[0]), pout(fd[1])
System::Pipe::Pipe ()
{
  pin.precision (14);
  pout.precision(14);
}

/*********************************************************************************************
 *                                     Dynamic Libraries                                     *
 ********************************************************************************************/

void System::DynLib::open (const std::string& lib) throw(Error::General) 
{    
    const char funame [] = "System::DynLib::_open: ";

    if(_lib == lib)
	return;

    _delete_ref();

    _lib.clear();
    _count = 0;
    _handle = 0;

    if(!lib.size())
	return;

    _handle = dlopen(lib.c_str(), RTLD_NOW);
    if(!_handle) {
	std::cerr << funame << "error encountered while opening " 
		  << lib << " library:\n";
	const char* messg = dlerror();
	if(messg)
	    std::cerr <<messg << "\n";
	throw Error::Init();
    }
    
    _lib = lib;
    _count = new int(1);

}

void System::DynLib::_delete_ref ()
{
    const char funame [] = "System::DynLib::_delete_ref: ";
    
    if(!_count)
	return;
    
    if(!(--(*_count))) {
	delete _count;
	_count  = 0;
	if(dlclose(_handle)) {
	    std::cerr << funame << "error occurred while closing " 
		      << _lib << " library:\n";
	    const char* messg = dlerror();
	    if(messg)
		std::cerr << messg << "\n";
	}		
	_lib.clear();
	_handle = 0;

    }
}

void* System::DynLib::member (const std::string& sym) throw(Error::General)
{
    const char funame [] = "System::DynLib::member: ";

    if(!_handle) {
	std::cerr << funame << "library has not been opened\n";
	throw Error::Init();
    }

    dlerror();
    void* res = dlsym(_handle, sym.c_str());
    const char* messg = dlerror();
    if(messg) {
	std::cerr << funame << "failed to get " << sym << " symbol from the " 
		  << _lib  << " library: " << messg << "\n";
	throw Error::Init();
    }
    return res;
}

