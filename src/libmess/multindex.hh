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

#ifndef MULTINDEX_HH
#define MULTINDEX_HH

#include "error.hh"
#include "io.hh"

#include <vector>
#include <set>
#include <iostream>
#include <map>

// symmetric tensor index
//
class SymIndex : private std::vector<int> {
  //
  int _range;
    
  void _isinit () const;

  SymIndex ();

public:

  explicit SymIndex(int r) : _range(r) { _isinit(); }

  void operator++ ();
  void operator++ (int) { operator++(); }

  int range () const { return _range; }

  int size () const { return std::vector<int>::size(); }

  int operator() (int i) const { return (*this)[i]; }

  operator std::multiset<int> () const;

  operator std::map<int, int> () const;
};

inline void SymIndex::_isinit () const
{
  const char funame [] = "Thermo::SymIndex::_isinit: ";
    
  if(_range <= 0) {
    //
    ErrOut err_out;

    err_out << funame << "not initialized properly";
  }
}

// generic fixed size multi-index 
//
class GenIndex : private std::vector<int> {
  //
  int _range;

  bool _fin;
    
  void _isinit () const;

  GenIndex ();

public:

  GenIndex(int s, int r) : std::vector<int>(s), _range(r), _fin(false) { _isinit(); }

  void operator++ ();
  
  void operator++ (int) { operator++(); }

  int range () const { return _range; }

  int size () const { return std::vector<int>::size(); }

  bool fin () const { return _fin; }

  int operator() (int i) const { return (*this)[i]; }

  operator std::multiset<int> () const;
  //
  operator std::map<int, int> () const;
};

inline void GenIndex::_isinit () const
{
  const char funame [] = "Thermo::GenIndex::_isinit: ";
    
  if(_range <= 0) {
    //
    ErrOut err_out;

    err_out << funame << "not initialized properly";
  }
}

//
// numerical derivative generator
//
class NumDer {
  //
  // configurational space dimensionality
  //
  int _size;
    
  void _check ()  const;

public:
  //
  // derivative signature
  //
  typedef std::map<int, int> der_t;

  // configuration map
  //
  typedef std::map<std::vector<int>, int> cmap_t;

  // function-on-the-grid data
  //
  typedef std::map<std::vector<int>, double> fdata_t;

  NumDer () : _size(0) {}
    
  void resize (int s) { _size = s; _check(); }

  explicit NumDer (int s) { resize(s); }

  int size () const { return _size; }

  // converts derivative signature to configuration map
  //
  cmap_t convert (const der_t& der) const;

  // numerical derivative from the function-on-the-grid data
  //
  double operator() (const der_t& der, const fdata_t&  fdata, const std::vector<double>& step) const;

  static int order (const der_t& der);

private:

  void _assert (const der_t&) const;
};

inline void NumDer::_check () const
{
  const char funame [] = "Thermo::NumDer::_check: ";

  if(_size <= 0) {
    //
    ErrOut err_out;
      
    err_out << funame << "dimension out of range: " << _size;
  }
}

/************************************************************************************************
 ************************************** MULTI-DIMENSIONAL INDEX *********************************
 ************************************************************************************************/

// several groups of ordered indices of the same size
//
// index range: 0 <= index <= index_max

class MultiIndex : private std::vector<int> {
  bool _end;

  std::vector<int> _group_rank;  // ordered group rank
  std::vector<int> _group_size;  // ordered group size
  std::vector<int> _group_end;   // ordered group end index

  std::vector<int> _group_map;    

  template <typename V>
  int _compare (const V&) const ;

  void _clear () { clear(); _group_rank.clear(); _group_size.clear(); 
    _group_end.clear(); _group_map.clear(); _end = true; }

public:
  enum {
    ORDERED = 1
  };

  void set (const std::vector<int>&, const std::vector<int>&) ;

  MultiIndex () : _end(true) {}

  MultiIndex(int, int, int =0);
  explicit MultiIndex(const std::vector<int>&);
  MultiIndex(const std::vector<int>&, int);
  MultiIndex(const std::vector<int>&, const std::vector<int>&);

  bool end()  const { return _end; }

  int       size ()      const { return std::vector<int>::size(); }
  int operator[] (int i) const { return std::vector<int>::operator[](i); }

  const std::vector<int>&  base () const { return *this; }
  std::vector<std::multiset<int> > convert_multiset () const;

  void operator++() ;
  void operator++(int)  { operator++(); }

  int cumulative_index () const;

  template<class T> bool operator== (const T& t) const { return _compare(t) == 0; }
  template<class T> bool operator!= (const T& t) const { return _compare(t) != 0; }
  template<class T> bool operator<  (const T& t) const { return _compare(t) <  0; }
  template<class T> bool operator>  (const T& t) const { return _compare(t) >  0; }
  template<class T> bool operator<= (const T& t) const { return _compare(t) <= 0; }
  template<class T> bool operator>= (const T& t) const { return _compare(t) >= 0; }

  void reset () { if(!size()) return; _end = false; int itemp = size(); clear(); resize(itemp, 0); }
};

inline MultiIndex::MultiIndex (int d, int m, int flags) 
{
  if(flags & ORDERED)
    set(std::vector<int>(1, d), std::vector<int>(1, m));
  else
    set(std::vector<int>(d, 1), std::vector<int>(d, m));
}

inline MultiIndex::MultiIndex (const std::vector<int>& m) 
{
  set(std::vector<int>(m.size(), 1), m);
}

inline MultiIndex::MultiIndex (const std::vector<int>& g, int m) 
{
  set(g, std::vector<int>(g.size(), m));
}

inline MultiIndex::MultiIndex (const std::vector<int>& g, const std::vector<int>& m) 
{
  set(g, m);
}

template<class T>
int MultiIndex::_compare (const T& m) const 
{
  const char funame [] = "MultiIndex::_compare: ";
    
  if(m.size() != size()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Range();
  }

  typename T::const_iterator mit = m.begin();
  for(const_iterator it = begin(); it != std::vector<int>::end(); ++mit, ++it)
    if(*it < *mit)
      return -1;
    else if(*it > *mit)
      return 1;

  return 0;
}

  // convert multi-dimensional index to one-dimensional index
template<class T>
int multi2one (const T& pr, int n) 
{
  const char funame [] = "multi2one: ";
  
  if(n <= 1) {
    std::cerr << funame << "wrong base\n";
    throw Error::Range();
  }

  int res = 0;
  int base = 1;
  for(typename T::const_iterator mit = pr.begin(); mit != pr.end(); ++mit, base *= n) {
    if(*mit < 0 || *mit >= n) {
      std::cerr << funame << " out of range\n";
      throw Error::Range();
    }
    res += *mit * base;
  }
  return res;
}

/************************************************************************************************
 **************************** MULTI-DIMENSIONAL INDEX CONVERTER *********************************
 ************************************************************************************************/

class MultiIndexConvert : private std::vector<int> {
  long _linear_size;

public:
  void resize (const std::vector<int>&) ;
  void resize (int d, int s) { std::vector<int> m(d, s); resize(m); }

  MultiIndexConvert          () : _linear_size(0)       {}
  MultiIndexConvert          (int d, int s)             { resize(d, s); }
  explicit MultiIndexConvert (const std::vector<int> s) { resize(s); }

  int  size (int i) const { return (*this)[i];     }
  int  rank ()      const { return std::vector<int>::size(); }
  long size ()      const { return _linear_size; }

  std::vector<int> operator() (long)                     const ;
  long             operator() (const std::vector<int>&) const ;
  
  long  conjugate (long) const;
};

#endif
