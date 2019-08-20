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

#ifndef SYMMETRY_HH
#define SYMMETRY_HH

#include "d3.hh"
#include "permutation.hh"

namespace Symmetry {

  /***************************************************************************************************
   ************************************** ABSTRACT GROUP CLASS ***************************************
   ***************************************************************************************************/

  // unity corresponds to 0
  class GroupBase : private std::vector<Permutation> {
    bool                _init;
    std::vector<int> _inverse;

    void _isinit () const ;

  public:
    void init (const std::vector<Permutation>&) ;

    GroupBase () : _init(false) {}
    GroupBase (const std::vector<Permutation>& p) { init(p); }

    int    size ()             const { _isinit(); return std::vector<Permutation>::size(); }    
    int product (int e, int g) const { _isinit(); return (*this)[e][g]; }
    int inverse (int e)        const { _isinit(); return   _inverse[e]; }
    
    bool operator== (const GroupBase& g) const { 
      return std::operator==((const std::vector<Permutation>&)(*this), (const std::vector<Permutation>&)g);
    }
    bool operator!= (const GroupBase& g) const { return !operator==(g); }

    const std::vector<Permutation>& multiplication_table () const { return *this; }
  };

  /************************************************************************************************
   ******************************** SPATIAL SYMMETRY OPERATIONS ***********************************
   ************************************************************************************************/

  // spatial symmetry element: rotation and inversion 
  class SpaceElement : private Quaternion {
    bool       _inversion;
    D3::Matrix _rmatrix;

  public:
    SpaceElement () : Quaternion(1.), _inversion(false), _rmatrix(1.) {}

    SpaceElement (const Quaternion& q, bool i, int = 0) ;

    SpaceElement (const D3::Matrix& m, int = 0)          ;

    operator const Quaternion& () const { return *this; }

    const double& operator[] (int i) const { return Quaternion::operator[](i); }

    SpaceElement operator* (const SpaceElement&) const;
    
    bool inversion () const { return _inversion; }

    int     compare (const SpaceElement& e) const;
    bool operator== (const SpaceElement& e) const { if(compare(e) == 1) return true; return false; }
    bool operator!= (const SpaceElement& e) const { return !operator==(e); }

    const D3::Matrix& rmatrix () const { return _rmatrix; }

    // symmetry application to 3D vector
    void apply (const double*, double*) const;

    // symmetry application to the molecule
    template <typename V>
    void apply(const std::vector<V>&, std::vector<V>&) const ;
  };

  template <typename V>
  void SpaceElement::apply(const std::vector<V>& m1, std::vector<V>& m2) const 
  {
    const char funame [] = "Symmetry::SpaceElement::apply: ";

    if(m1.size() != m2.size()) {
      std::cerr << funame << "sizes mismatch\n";
      throw Error::Range();
    }
    
    typename std::vector<V>::const_iterator cat = m1.begin();
    for(typename std::vector<V>::iterator bat = m2.begin(); bat != m2.end(); ++bat, ++cat)
      apply((const double*)(*cat), (double*)(*bat));
  }

  inline int compare (const SpaceElement& e, const SpaceElement& g) { return e.compare(g); }

  /******************************************** SPECIFIC SYMMETRY ELEMENTS *****************************************/

  // unity space symmetry element
  inline SpaceElement unity () 
  {
    return SpaceElement(Quaternion(1.), 0, Quaternion::NOCHECK);
  }

  // rotation about reference frame axis by 180 degrees  
  inline SpaceElement rotation (int i)
  {
    const char funame [] = "Symmetry::rotation: ";

    if(i < 0 || i > 2) {
      std::cerr << funame << "index out of range\n";
      throw Error::Range();
    }

    Quaternion q;
    q[i + 1] = 1.;
    return SpaceElement(q, 0, Quaternion::NOCHECK);
  }

  // rotation about reference frame axis by arbitrary angle
  inline SpaceElement rotation (int i, double angle)
  {
    const char funame [] = "Symmetry::rotation: ";

    if(i < 0 || i > 2) {
      std::cerr << funame << "index out of range\n";
      throw Error::Range();
    }

    Quaternion q;
    q[0]     = std::cos(angle / 2.);
    q[i + 1] = std::sin(angle / 2.);

    return SpaceElement(q, 0, Quaternion::NOCHECK);
  }

  inline SpaceElement rotation (int i, int order)
  {
    const char funame [] = "Symmetry::rotation: ";

    if(order < 2) {
      std::cerr << funame << "order out of range\n";
      throw Error::Range();
    }
    return rotation(i, 2. * M_PI / (double)order);
  }

  // rotation about arbitrary axis by 180 degrees
  inline SpaceElement rotation (const double* n)
  {
    Quaternion q;
    for(double* it = q.begin() + 1; it != q.end(); ++it, ++n)
      *it = *n;
    
    return SpaceElement(q, false);    
  }

  // rotation about arbitrary axis by arbitrary angle
  inline SpaceElement rotation (const double* n, double angle) 
  {
    Quaternion q;
    q[0] = std::cos(angle / 2.);

    double dtemp = std::sin(angle / 2.);
    for(double* it = q.begin() + 1; it != q.end(); ++it, ++n)
      *it = *n * dtemp;
    
    return SpaceElement(q, false);    
  }

  inline SpaceElement rotation (const double* n, int order)
  {
    return rotation(n, 2. * M_PI / (double)order);
  }

  // inversion 
  inline SpaceElement inversion () 
  {
    return SpaceElement(Quaternion(1.), 1, Quaternion::NOCHECK);
  }

  // reflection of the reference frame plane
  inline SpaceElement reflection (int i)
  {
    const char funame [] = "Symmetry::reflection: ";

    if(i < 0 || i > 2) {
      std::cerr << funame << "index out of range\n";
      throw Error::Range();
    }

    Quaternion q;
    q[i + 1] = 1.;
    return SpaceElement(q, 1, Quaternion::NOCHECK);
  }

  // reflection of the arbitrary plane 
  inline SpaceElement reflection (const double* n)
  {
    Quaternion q;
    for(double* it = q.begin() + 1; it != q.end(); ++it, ++n)
      *it = *n;

    return SpaceElement(q, true);    
  }

  /************************************************************************************************
   *********************************** SPATIAL SYMMETRY GROUP *************************************
   ************************************************************************************************/

  class SpaceGroup : public GroupBase, private std::vector<SpaceElement> { 
    std::string _name;
    Array_2<int> _sign;
    
    typedef std::vector<SpaceElement>::iterator             iterator;
    typedef std::vector<SpaceElement>::const_iterator const_iterator;

    iterator       begin ()       { return std::vector<SpaceElement>::begin(); }
    const_iterator begin () const { return std::vector<SpaceElement>::begin(); }
    iterator         end ()       { return std::vector<SpaceElement>::end();   }
    const_iterator   end () const { return std::vector<SpaceElement>::end();   }

  public:
    SpaceGroup (const std::string&, const std::vector<SpaceElement>&, int =0) ;

    const std::string& name () const { return _name; }

    bool operator== (const SpaceGroup& g) const;
    bool operator!= (const SpaceGroup& g) const { return !operator==(g); }

    int                              sign (int i, int j) const { return _sign(i, j); }
    const std::vector<SpaceElement>& base ()             const { return *this; }
    const SpaceElement&        operator[] (int g)        const { return base()[g]; }
    int                              size ()             const { return GroupBase::size(); }

  };

  inline bool SpaceGroup::operator== (const SpaceGroup& g) const 
  { 
    return GroupBase::operator==(g) && 
      std::operator==((const std::vector<SpaceElement>&)(*this), (const std::vector<SpaceElement>&)g);
  }
  

}//symmetry

#endif
