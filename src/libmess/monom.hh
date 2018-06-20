/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
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
    iterator (const std::vector<int>&) throw(Error::General);

    void operator++ ();
    void operator++ (int) { operator++(); }
    bool end () const { return _end; }

    const std::vector<int>& operator* () const { return *this; }
  };
  
  iterator _begin;

public:
  Monom (int rn, int sz) throw(Error::General); 
  
  int rank () const { return _rank; }
  int size () const { return _size; }
  int linear_size () const { return _index_multi_map.size(); }

  const std::vector<int>& operator() (int i)                    const { return _index_multi_map[i]; }
  Lapack::Matrix          operator() (const Lapack::Matrix&)    const throw(Error::General);

};

#endif
