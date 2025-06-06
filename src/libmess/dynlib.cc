#include "dynlib.hh"

#include <dlfcn.h>

/*********************************************************************************************
 *                                     Dynamic Libraries                                     *
 ********************************************************************************************/

void DynLib::open (const std::string& lib)  
{    
    const char funame [] = "DynLib::_open: ";

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

void DynLib::_delete_ref ()
{
    const char funame [] = "DynLib::_delete_ref: ";
    
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

void* DynLib::member (const std::string& sym) 
{
    const char funame [] = "DynLib::member: ";

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
