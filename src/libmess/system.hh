

#ifndef SYSTEM_HH
#define SYSTEM_HH

#include <sys/types.h>
#include <fstream>
#include <iostream>
#include "error.hh"
#include "io.hh"

namespace System {

    // delete all files in the directory
    void clean_dir (const char*);

    // file copy
    int file_copy(const char*, const char*);

    /**********************************************************************
     * call to external executable (substitute for system());             *
     * list of arguments in call_exe() should be terminated by 0 pointer  *
     *********************************************************************/

    int call_exe (const char* ...);

    /*************************************************************************
     *                                 Pipes                                 *
     *************************************************************************/

    class Pipe_base
    {
	int fd [2];

	Pipe_base (const Pipe_base&);
	Pipe_base& operator= (const Pipe_base&);

	Pipe_base () throw(Error::General);

	friend class Pipe;
    };

    class Pipe : public Pipe_base
    {
	Pipe (const Pipe&);
	Pipe& operator= (const Pipe&);

    public:

	std::ifstream pin;
	std::ofstream pout;
	Pipe();
    };

    /*********************************************************************************************
     *                                     Dynamic Libraries                                     *
     ********************************************************************************************/

    class DynLib : public IO::Read {
	std::string _lib;
	void* _handle;
	int*  _count;

	void _delete_ref ();
	void _create_ref (const DynLib& dl);

    public:
	void open (const std::string&) throw(Error::General);

	DynLib () : _handle(0), _count(0) {}
	explicit DynLib (const std::string& lib) throw(Error::General) 
	    : _handle(0), _count(0) { open(lib); }

	DynLib (const DynLib& dl) { _create_ref(dl); }
	DynLib& operator= (const DynLib& dl) { _delete_ref(); _create_ref(dl); return *this; }

	virtual ~DynLib () { _delete_ref(); }
	virtual void read (std::istream&) throw(Error::General);
	
	bool isopen () const;
	void* member (const std::string&) throw(Error::General);
    };

    inline bool DynLib::isopen () const
    {
	if(_handle)
	    return true;
	return false;
    }

    inline void DynLib::read (std::istream& from) throw(Error::General)
    {
	const char funame [] = "System::DynLib::read: ";
    
	std::string lib;
	from >> lib;
	if(!from) {
	    std::cerr << funame << "input stream is corrupted\n";
	    throw Error::Form();
	}
	open(lib);
    }

    inline void DynLib::_create_ref (const DynLib& dl)
    { 
	_handle = dl._handle;
	_count = dl._count;
	_lib = dl._lib;
	if(_count)
	    ++(*_count); 
    }
 
}// System

#endif

