/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2019, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

// different atomic coordinates representations (zmatrix, cartesian, ...)
//

#ifndef COORD_HH
#define COORD_HH

#include <iostream>
#include <set>

#include "array.hh"
#include "d3.hh"

namespace Coord {

  // internal coordinate: distance, angle, or dihedral
  //
  double internal (const std::vector<D3::Vector>&);
  
  // cartesian coordinates
  //
  class Cartesian : public Array<double> {
    //
    void _assert (int) const;

    int _atom_size;
    
  public:
    //
    Cartesian () : _atom_size(0) {}
    
    explicit Cartesian (int s) : Array<double>(s), _atom_size(s / 3) { _assert(s); };

    void resize (int s) { _assert(s); _atom_size = s / 3; Array<double>::resize(s); }

    int atom_size () const { return _atom_size; }

    double*       atom_pos (int a)       { return *this + 3 * a; }
    
    const double* atom_pos (int a) const { return *this + 3 * a; }
  };

  // z-matrix coordinates
  //
  class ZMat : private std::vector<double> {
    //
    int _atom_size;

    void _assert(int, int) const;
    
  public:
    //
    enum {DISTANCE, ANGLE, DIHEDRAL};

    ZMat () : _atom_size(0) {}

    ZMat (int s) : std::vector<double>(3 * s), _atom_size(s) {}

    // number of atoms
    //
    int atom_size () const { return _atom_size; }

    // element access
    //
    double   operator() (int v, int a) const { _assert(v, a); return (*this)[v + 3 * a]; }

    double&  operator() (int v, int a)       { _assert(v, a); return (*this)[v + 3 * a]; }
  };

  // z-matrix signature
  //
  class ZBase : private std::vector<int> {

    void _assert (int, int) const;

    void _assert () const;
    
  public:
    //
    ZBase () {}

    void init (std::istream&);

    void init (const std::vector<int>& z) { *this = z; _assert(); }

    ZBase (std::istream& from) { init(from); }

    ZBase (const std::vector<int>& z) : std::vector<int>(z) { _assert(); }
    
    bool isinit () const { return size(); }

    int atom_size () const { return size() / 3; }

    int operator () (int v, int a) const { _assert(v, a); return (*this)[v + 3 * a]; }
    
    // converters
    //
    Cartesian operator() (const ZMat&)      const;

    ZMat      operator() (const Cartesian&) const;
  };

  // abstract function in cartesian space
  //
  class CartFun {
    //
    // dependent atoms pool
    //
    std::set<int> _dep_pool;

    void _assert (int, int = 0) const;
    
    bool _depend (int) const;
    
  protected:
    //
    bool add_dep (int a) { _assert(a); return _dep_pool.insert(a).second; }
    
  public:
    //
    virtual double evaluate (const Cartesian&) const = 0;

    // first derivative
    //
    double grad (const Cartesian&, int) const;

    // second derivative
    //
    double hess (const Cartesian&, int, int) const;
    
    double grad_length (const Cartesian&) const;
    
    static double step;
  };

  inline void CartFun::_assert (int i, int imax) const
  {
    if(imax > 0 && i >= imax || i < 0) {
      //
      std::cerr << "Coord::CartFun::_assert: out of range: " << i << ", " << imax << "\n";

      throw Error::Range();
    }
  }

  inline bool CartFun::_depend (int i) const
  {
    if(_dep_pool.size() && _dep_pool.find(i / 3) == _dep_pool.end())
      //
      return false;

    return true;
  }
}

#endif
