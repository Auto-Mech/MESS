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
#include "atom.hh"
#include "lapack.hh"

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

  /********************************************************************************************
   ****************************************** Z-MATRIX ****************************************
   ********************************************************************************************/
  
  enum {DISTANCE, ANGLE, DIHEDRAL};

  class ZRecord : private std::pair<AtomBase, std::vector<std::pair<int, std::string> > > {
    //
  public:
    //
    AtomBase&       atom ()       { return first; }

    const AtomBase& atom () const { return first; }

    int size () const { return second.size(); }

    void resize (int s) { second.resize(s); }

    // reference atom index
    //
    int  ref (int i) const { return second[i].first; }
    
    int& ref (int i)       { return second[i].first; }

    const std::string& var (int i) const { return second[i].second; }
    
    std::string&       var (int i)       { return second[i].second; }
  };

  class ZMat : public std::vector<ZRecord> {
    //
    // variable name to index map
    //
    typedef std::map<std::string, std::pair<int, int> > _vmap_t;
    
    _vmap_t _var_to_ind;

    const std::pair<int, int>& _pos (const std::string&) const;

  public:
    //
    ZMat () {}
    
    bool isinit () const { return size(); }

    void init (std::istream&);
    
    ZMat (std::istream& from) { init(from); }

    // name to index convertion
    //
    int v2i (const std::string& var) const { return _pos(var).first + 3 * _pos(var).second; }

    // internal variable signature
    //
    std::vector<int> sign (const std::string&) const;

    // is variable used
    //
    bool isvar (const std::string& var) const { if(_var_to_ind.find(var) != _var_to_ind.end()) return true; return false; }

    int type (const std::string& var) const { return _pos(var).first; }

    int record_number (const std::string& var) const { return _pos(var).second; }

    void cm_shift (Cartesian&) const;

    double conv_factor (const std::string& v) const {if(type(v) == DISTANCE) return Phys_const::angstrom; return M_PI / 180.; }
  };
  
  inline const std::pair<int, int>& ZMat::_pos (const std::string& var) const
  {
    _vmap_t::const_iterator vit = _var_to_ind.find(var);
  
    if(vit != _var_to_ind.end())
      //
      return vit->second;

    std::cerr << "Coord::ZMat::_pos: unknown varable name: " << var << "\n";

    throw Error::Init();
  }

  /***************************************************************************************************
   ****************************************** Z-MATRIX DATA ******************************************
   ***************************************************************************************************/
  
  class ZData : private std::vector<double> {
    //
    const ZMat& _base;

    void _assert(int, int) const;
    
  public:
    //
    ZData (const ZMat& z) : std::vector<double>(3 * z.size()), _base(z) {}

    const ZMat& zmat () const { return _base; }
    
    // number of atoms
    //
    int atom_size () const { return _base.size(); }

    // element access
    //
    double  operator[] (const std::string& v) const { return std::vector<double>::operator[](_base.v2i(v)); }
    
    double& operator[] (const std::string& v)       { return std::vector<double>::operator[](_base.v2i(v)); }

    double   operator() (int v, int a) const { _assert(v, a); return at(v + 3 * a); }

    double&  operator() (int v, int a)       { _assert(v, a); return at(v + 3 * a); }

    int type (const std::string& v) const { return _base.type(v); }

    double conv_factor (const std::string& v) const { return _base.conv_factor(v); }
    
    // conversion
    //
    operator Cartesian () const;

    void import (const Cartesian&);

    Lapack::Vector inertia_moments () const;

    Lapack::SymmetricMatrix mobility_matrix () const;
  };

  /********************************************************************************************
   ****************************** FUNCTION IN CARTESIAN SPACE *********************************
   ********************************************************************************************/
  
  // abstract function
  //
  class CartFun {
    //
    // dependent atoms pool
    //
    std::set<int> _dep_pool;

    void _assert (int, int = 0) const;
    
    bool _depend (int) const;

    void _adjust_dihedral (double&) const;

  protected:
    //
    bool add_dep (int a) { _assert(a); return _dep_pool.insert(a).second; }

    int _type;
    
  public:
    //
    CartFun () : _type(-1) {}

    virtual double evaluate (const Cartesian&) const = 0;

    virtual int atom_size () const = 0;

    int type () const { return _type; }

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
  
  // internal coordinate: distance, angle, or dihedral
  //
  class Internal : public CartFun, private std::vector<int> {
    //
    int _atom_size;
    
    void _assert();
    
  public:
    //
    Internal(int, std::istream&);

    Internal (int as, const std::vector<int>& v) : std::vector<int>(v), _atom_size(as) { _assert(); }

    int atom_size () const { return _atom_size; }
    
    double evaluate (const Cartesian&) const;
  };
}

#endif
