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

	Pipe_base () ;

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

	explicit Semaphore (int) ; // creates set of n semaphores & initilizes them 
	Semaphore (key_t, int)   ; // initializes existing semaphore set
	~Semaphore();
  
	key_t get_key () const {return key;}
	void busy (int) const ; // raise n-th semaphore (P, wait)
	void free (int) const ; // free n-th semaphore  (V, signal)
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
	void open (const std::string&) ;

	DynLib () : _handle(0), _count(0) {}
	explicit DynLib (const std::string& lib)  
	    : _handle(0), _count(0) { open(lib); }

	DynLib (const DynLib& dl) { _create_ref(dl); }
	DynLib& operator= (const DynLib& dl) { _delete_ref(); _create_ref(dl); return *this; }

	virtual ~DynLib () { _delete_ref(); }
	virtual void read (std::istream&) ;
	
	bool isopen () const;
	void* member (const std::string&) ;
    };

    inline bool DynLib::isopen () const
    {
	if(_handle)
	    return true;
	return false;
    }

    inline void DynLib::read (std::istream& from) 
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

