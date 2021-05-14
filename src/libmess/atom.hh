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

#ifndef ATOM_HH
#define ATOM_HH

#include <iomanip>
#include <map>
#include <set>

#include "d3.hh"
#include "units.hh"
#include "permutation.hh"
#include "symmetry.hh"
#include "lapack.hh"

/*************************** Atom description ***************************/

class AtomBase 
{
  int          _num; // atomic number 
  int         _isot; // isotope
  double      _mass; // atom mass
  //std::string _name; // atom name

  struct Data {
    std::string              name;
    int                   isotope;
    int                   valence;
    std::map<int, double>    mass;
  };

  class DataBase : private std::map<int, Data> {
    std::map<std::string, int> _name_num_map;
    
  public:
    DataBase ();

    int              number (const std::string&) const ;
    const std::string& name (int)                const ;
    int     default_isotope (int)                const ;
    int             valence (int)                const ;
    double             mass (int, int)           const ;
  };

  static DataBase _data;

  void _set () { _mass = _data.mass(_num, _isot); }

public: 
  enum { // atomic number
    DUMMY      =  0,
    HYDROGEN   =  1,
    HELIUM     =  2,
    CARBON     =  6,
    NITROGEN   =  7,
    OXYGEN     =  8,
    FLUORINE   =  9,
    SODIUM     = 11,
    SILICON    = 14,
    PHOSPHORUS = 15,
    SULFUR     = 16,
    CHLORINE   = 17,
    TITANIUM   = 22,
    BROMINE    = 35,
    URANIUM    = 92
  };

  void set         (int)                     ;
  void set         (int, int)                ;
  void set         (const std::string&)      ;
  void set         (const std::string&, int) ;
  void set_isotope (int)                     ;

  AtomBase () : _num(DUMMY), _isot(0), _mass(0.) {}
  explicit AtomBase (const std::string& s)         { set(s);    }
  AtomBase          (const std::string& s, int i)  { set(s, i); }
    
  int       number () const { return _num; }
  int      isotope () const { return _isot; }
  double      mass () const { return _mass; }
  const char* name () const { return _data.name(_num).c_str(); }
  int      valence () const { return _data.valence(_num); }
};

inline void AtomBase::set (int n)  
{
  _num  = n; 
  _isot = _data.default_isotope(_num);

  _set();
}

inline void AtomBase::set (int n, int i)  
{
  _num  = n; 
  _isot = i;

  _set();
}

inline void AtomBase::set (const std::string& s)  
{
  _num  = _data.number(s); 
  _isot = _data.default_isotope(_num);

  _set();
}

inline void AtomBase::set (const std::string& s, int i) 
{
  _num  = _data.number(s);
  _isot = i;

  _set();
}

inline void AtomBase::set_isotope (int i) 
{
  _isot = i; 

  _set();
}

inline bool operator< (const AtomBase& a, const AtomBase& b)
{ 
  return a.number() < b.number() || a.number() == b.number() && a.isotope() < b.isotope();
}

inline bool operator> (const AtomBase& a, const AtomBase& b)
{ 
  return a.number() > b.number() || a.number() == b.number() && a.isotope() > b.isotope();
}

inline bool operator== (const AtomBase& a, const AtomBase& b)
{ 
  return a.number() == b.number() &&  a.isotope() == b.isotope();
}

inline bool operator!= (const AtomBase& a, const AtomBase& b)
{ 
  return !(a == b);
}

inline bool operator<= (const AtomBase& a, const AtomBase& b)
{
  return a < b || a == b;
}

inline bool operator>= (const AtomBase& a, const AtomBase& b)
{
  return a > b || a == b;
}

class Atom : public AtomBase, public D3::Vector
{
  void _read(std::istream&) ;

public:
  Atom () {}
  Atom          (const std::string& s, int i) : AtomBase(s, i) {}
  Atom          (const AtomBase& b)           : AtomBase(b)    {}
  explicit Atom (const std::string& s)        : AtomBase(s)    {}
  explicit Atom (std::istream& from)  { _read(from); }

  Atom& operator= (const D3::Vector& v) { D3::Vector::operator=(v); return *this; }
  
  friend std::istream& operator>> (std::istream&, Atom&) ;
};

std::ostream& operator<< (std::ostream& out , const Atom& a);

inline std::istream& operator>> (std::istream& from, Atom& a)  
{ 
  a._read(from); 
  return from; 
}

enum {IGNORE_ISOTOPE = 1};

bool are_equal (double, double, double);
bool are_equal (const Atom&, const Atom&, double, int =0);
bool are_equal (const std::vector<Atom>&, const std::vector<Atom>&, double, int =0);

Permutation is_symmetric (const std::vector<Atom>&, const Symmetry::SpaceElement&, double, int =0);

std::set<Permutation> identical_atoms_permutation_symmetry_group (const Lapack::SymmetricMatrix& dist, double tol) ;

std::set<Permutation> permutation_symmetry_group (const std::vector<Atom>&, double, int = 0) ;

std::pair<int, int> symmetry_number (const std::vector<Atom>&, double, int = 0) ;

Symmetry::SpaceGroup spatial_symmetry_group (std::vector<Atom>&, double, int = 0) ;

#endif
