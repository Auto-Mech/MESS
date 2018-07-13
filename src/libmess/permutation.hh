

#ifndef PERMUTATION_HH
#define PERMUTATION_HH

#include <set>
#include <vector>
#include <iostream>

#include "error.hh"

/************************************************************************************************
 ************************************ INDEX PERMUTATION *****************************************
 ************************************************************************************************/

class Permutation : private std::vector<int> {
  bool _end;
  int _compare (const Permutation&) const;
  
public:
  enum { NOCHECK = 1 };
    
  Permutation () : _end(true) {}

  // identical permutation
  explicit Permutation (int);
  // explicit permutation description
  explicit Permutation (const std::vector<int>&, int = 0);
  // permutation from orbit
  Permutation (const std::set<std::vector<int> >&, int);
  // two elements permutation
  Permutation (int, int, int);

  int                     size ()      const { return std::vector<int>::size(); }
  int               operator[] (int i) const { return std::vector<int>::operator[](i); }
  const std::vector<int>& base ()      const { return *this; }

  operator bool () { return size(); }

  // elments orbits
  std::set<std::vector<int> > orbit () const;
  std::vector<int>        orbit_max () const;

  // group operations
  template<typename V> 
  V operator() (const V&) const;

  Permutation operator* (const Permutation&) const;
  Permutation invert () const;

  // index operations
  bool operator== (const Permutation& p) const { return _compare(p) == 0; }
  bool operator!= (const Permutation& p) const { return _compare(p) != 0; }

  bool operator<  (const Permutation& p) const { return _compare(p) < 0; }
  bool operator>  (const Permutation& p) const { return _compare(p) > 0; }

  bool operator<= (const Permutation& p) const { return _compare(p) <= 0; }
  bool operator>= (const Permutation& p) const { return _compare(p) >= 0; }

  void operator++ ();
  void operator++ (int) { operator++(); }
  bool end () const { return _end; }
  void clear () { _end = false; }

  static std::set<std::vector<int> > read_orbit(std::istream&);
  static std::vector<int>            read_cycle(std::istream&);

  static int factorial (int);
};

template<typename V> 
V Permutation::operator() (const V& v) const
{
  Exception::Base funame = "Permutation::operator(): ";

  if(v.size() != size())
    throw funame << "dimensions mismatch\n";

  V res(size());
  for(int i = 0; i < size(); ++i)
    res[(*this)[i]] = v[i];

  return res;
}

std::set<Permutation> permutation_group (const std::set<Permutation>&);
std::set<Permutation> permutation_group (int);

/************************************************************************************************
 ***************************** ITERATION OVER MULTIPLE PERMUTATIONS *****************************
 ************************************************************************************************/

class MultiPerm : private std::vector<Permutation> {
  bool _end;
  
public:
  MultiPerm (const std::vector<int>&) throw(Error::General);

  const Permutation& operator[] (int i) const { return std::vector<Permutation>::operator[](i); }
  int                      size ()      const { return std::vector<Permutation>::size(); }
  bool                      end ()      const { return _end; }

  void operator++ ();
  void operator++ (int) { return operator++(); }
};

#endif
