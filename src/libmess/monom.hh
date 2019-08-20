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

#ifndef MONOM_HH
#define MONOM_HH

#include "error.hh"
#include "lapack.hh"

#include <vector>
#include <map>

/************************************************************************************************
 ******************************************** MONOMIAL ******************************************
 ************************************************************************************************/

class Monom {
  int _rank;
  int _size;

  std::vector<std::vector<int> > _index_multi_map;
  std::map<std::vector<int>, int>     _multi_index_map;

  class iterator : private std::vector<int> {
    bool _end;

  public:
    iterator () : _end(true) {}
    iterator (const std::vector<int>&) ;

    void operator++ ();
    void operator++ (int) { operator++(); }
    bool end () const { return _end; }

    const std::vector<int>& operator* () const { return *this; }
  };
  
  iterator _begin;

public:
  Monom (int rn, int sz) ; 
  
  int rank () const { return _rank; }
  int size () const { return _size; }
  int linear_size () const { return _index_multi_map.size(); }

  const std::vector<int>& operator() (int i)                    const { return _index_multi_map[i]; }
  Lapack::Matrix          operator() (const Lapack::Matrix&)    const ;

};

#endif
