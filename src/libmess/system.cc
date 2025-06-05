/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

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
//
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

int System::make_dir (const std::string& dir, const std::string& funame)
{
  std::string token = funame;
  token += dir + ": mkdir: ";

  struct stat wstat;

  if(mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH ))
    //
    switch(errno) {
      //
    case EEXIST:
      //
      if(stat(dir.c_str(), &wstat)) {
	//
	std::cerr << token << strerror(errno) << "\n";
	return 1;
      }

      if(!S_ISDIR(wstat.st_mode)) {
	//
	std::cerr << token << "not a directory\n";
	return 1;
      }
      
      if(!(wstat.st_mode & S_IRUSR) || !(wstat.st_mode & S_IWUSR) || !(wstat.st_mode & S_IXUSR)) {
	std::cerr << token << "wrong permissions\n";
	return 1;
      }

      return 0;

    case ENOENT:
      std::cerr << token << "the parent directory does not exist\n";
      return 1;

    default:
      std::cerr << token << strerror(errno) << "\n";
      return 1;
    }

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

System::Pipe_base::Pipe_base() 
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

/*********************************************************
 *                   Semaphores 
 *********************************************************/

#if defined(_POSIX_C_SOURCE) && !defined(_DARWIN_C_SOURCE)

union semun
{
  int val;                          // value for SETVAL
  struct semid_ds *buf;             // buffer for IPC_STAT & IPC_SET
  unsigned short int *array;        // array for GETALL & SETALL
  struct seminfo *__buf;            // buffer for IPC_INFO
};

#endif

System::Semaphore::Semaphore(int n) 
   : key(IPC_PRIVATE), num(n), creator(getpid())
{
    const char funame [] = "System::Semaphore::Semaphore(int): ";

   // create
   while((id = semget(++key, num, 0666 | IPC_CREAT | IPC_EXCL)) == -1)
   {} 

   // initialize
   semun su;
   Array<unsigned short> init_val (n);
   for (int i = 0; i < n; ++i)
      init_val[i] = 1;
   su.array = init_val;

   if (semctl(id, 0, SETALL, su) == -1)
      {
	 semctl(id, 0, IPC_RMID);
	 std::cerr << funame << "couldn't initialize\n";
	 throw Error::Init();
      }
}

System::Semaphore::Semaphore(key_t k, int n) 
   : key(k), num(n), creator(0)
{
    const char funame [] = "System::Semaphore::Semaphore(key_t, int): ";
    
    if((id = semget(k, n, 0666)) == -1) {
       std::cerr << funame << "couldn't open existing semaphore set\n";
       throw Error::Init();
    }
}

System::Semaphore::~Semaphore()
{
   if (getpid() == creator)
      semctl(id, 0, IPC_RMID);
} 

void System::Semaphore::busy (int n) const 
{
    const char funame [] = "System::Semaphore::busy: ";

    if (n >= num) {
	std::cerr << funame << "wrong semaphore number\n";
	throw Error::Range();
    }

   struct sembuf sb;
   sb.sem_num = n;
   sb.sem_op = -1;
   sb.sem_flg = SEM_UNDO;

   if (semop(id, &sb, 1) == -1) {
       std::cerr << funame << "failed\n";
       throw Error::Run();
   }
}

void System::Semaphore::free (int n) const 
{
    const char funame [] = "System::Semaphore::free: ";
    if (n >= num) {
	std::cerr << funame << "wrong semaphore number\n";
	throw Error::Range();
    }

    struct sembuf sb;
    sb.sem_num = n;
    sb.sem_op = 1;
    sb.sem_flg = SEM_UNDO;

    if (semop(id, &sb, 1) == -1) {
       std::cerr << funame << "failed\n";
       throw Error::Run();
    }
}


