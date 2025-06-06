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

#ifndef DYNLIB_HH
#define DYNLIB_HH

//#include <sys/types.h>
#include <fstream>
#include <iostream>
#include "error.hh"
#include "io.hh"

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
 
#endif
