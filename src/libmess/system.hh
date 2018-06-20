/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

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

    /*************************************************************************************
     *                                       Semaphores                                  *
     ************************************************************************************/

    class Semaphore
    {
	key_t key;   // semaphore key
	int   id;    // semaphore id
	int   num;   // number of semaphores in the set
	pid_t creator; // creator of semaphore

	Semaphore(const Semaphore&);
	Semaphore& operator= (const Semaphore&);

    public:

	explicit Semaphore (int) throw(Error::General); // creates set of n semaphores & initilizes them 
	Semaphore (key_t, int)   throw(Error::General); // initializes existing semaphore set
	~Semaphore();
  
	key_t get_key () const {return key;}
	void busy (int) const throw(Error::General); // raise n-th semaphore (P, wait)
	void free (int) const throw(Error::General); // free n-th semaphore  (V, signal)
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

