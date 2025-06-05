/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2025, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef DYNAMIC_HH
#define DYNAMIC_HH

#include "structure.hh"
#include "array.hh"
#include "logical.hh"
#include "lapack.hh"

#include <iostream>
#include <string>

namespace Dynamic {

  /************************************************************************************
   *********************** DYNAMICAL VARIABLES DERIVATIVES ****************************
   ************************************************************************************/

  void set_dvd (const D3::Vector* torque, const double* dv, double* dvd);

  /**************************************************************************************
   ***************************** GENERALIZED COORDINATES ********************************
   **************************************************************************************/

  class Coordinates : private Array<double>
  {
    mutable D3::Matrix _mfo [2];                // rotational matrix for a nonlinear fragment
    
    mutable std::vector<bool> _need_mfo_update; // delayed mfo update flags

    void _update_mfo (int frag) const;

    double _imm (int, int, int) const;

  public:
    //
    static int size () { return Structure::pos_size(); }

    static double ang_len_tol; // maximal allowed angular vector length deviation from unity

    void get (const double*, double* = 0); // copy from
    void put (double*) const; // copy to

    // constructors
    //
    Coordinates () : Array<double>(size(), 0.), _need_mfo_update(2, true) {}
    
    explicit Coordinates (const double* dc) : Array<double>(size()), _need_mfo_update(2, true) { get(dc); }

    // cm-to-cm vector, 1->2
    //
    const double* orb_pos  ()  const { return begin() + Structure::orb_pos(); }
    double*       orb_pos  ()        { return begin() + Structure::orb_pos(); }

    double orb_len () const { return ::vlength(orb_pos(), 3); }

    // orientational variables
    //
    const double*       ang_pos (int frag) const;
    double*       write_ang_pos (int frag);

    double get_ang_pos (int frag, const double* pos);

    double normalize (int frag);// normalize fragment angular vector

    void mf2lf (int frag, const double* mf, double* lf)  const;
    void lf2mf (int frag, const double* lf, double* mf)  const;

    void set_rel_pos (int frag, std::vector<D3::Vector>& rel_pos) const;

    Lapack::SymmetricMatrix imm () const;

    static double min_atom_dist;
    bool are_atoms_close () const;

    void print_geom (std::ostream&, const std::string&) const;

    void                        dipole (int frag,                  double* d) const;
    void     quadrupole_vector_product (int frag, const double* v, double* q) const;
    void polarizability_vector_product (int frag, const double* v, double* p) const;
  };

  inline const double* Coordinates::ang_pos (int f) const 
  {
    const char funame [] = "Dynamic::Coordinates::ang_pos: ";
    
    if(Structure::type(f) == Molecule::MONOATOMIC) {
      //
      std::cerr << funame << "should not be called for monoatomic molecule\n";

      throw Error::Logic();
    }

    return begin() + Structure::ang_pos(f);
  }

  inline double* Coordinates::write_ang_pos  (int f) 
  {
    const char funame [] = "Dynamic::Coordinates::write_ang_pos: ";
    
    switch(Structure::type(f)) {
      //
    case Molecule::MONOATOMIC:
      //
      std::cerr << funame << "should not be called for monoatomic molecule\n";

      throw Error::Logic();
      //
    case Molecule::NONLINEAR:
      //
      _need_mfo_update[f] = true;
    }
    
    return begin() + Structure::ang_pos(f);
  }

  inline void Coordinates::_update_mfo (int f) const 
  {
    if(_need_mfo_update[f] == false)
      //
      return;

    if(Structure::type(f) == Molecule::NONLINEAR)
      //
      quat2mat(ang_pos(f), _mfo[f]);

    _need_mfo_update[f] = false;
  }

  /*************************************************************************
   ************************ GENERALIZED MOMENTA ****************************
   *************************************************************************/

  class Momenta : public Array<double>
  {
    double* _orb_vel;
    double* _ang_vel [2];

    void _init(); // set pointers

  public:
    //
    static int size () { return Structure::vel_size(); }

    Momenta () : Array<double>(size()) { _init(); }
    Momenta (const double* dv) : Array<double>(size()) { Array<double>::operator=(dv);  _init(); }
    Momenta (const Momenta& m) : Array<double>((const Array<double>&)m) { _init(); }
    Momenta& operator= (const Momenta& m) { Array<double>::operator=((const Array<double>&)m); return *this; }

    void get (const double* dv) { Array<double>::operator=(dv); }
    void put (double*) const;

    // cm-to-cm velocity
    //
    double*       orb_vel ()       { return _orb_vel; }
    const double* orb_vel () const { return _orb_vel; }

    // rotational frequencies
    //
    double*       ang_vel (int frag);
    const double* ang_vel (int frag) const;

    double orbital_kinetic_energy        () const;
    double total_kinetic_energy          () const;
    double fragment_rotational_energy (int) const;
  };

  inline const double* Momenta::ang_vel (int frag) const 
  { 

    const char funame [] = "Dynamic::Momenta::ang_vel: ";

    if(Structure::type(frag) == Molecule::MONOATOMIC) {
      //
      std::cerr << funame << "should not be called for monoatomic molecule\n";

      throw Error::Logic();
    }

    return _ang_vel[frag];
  }

  inline double* Momenta::ang_vel (int frag) 
  { 
    const char funame [] = "Dynamic::Momenta::ang_vel: ";

    if(Structure::type(frag) == Molecule::MONOATOMIC) {
      //
      std::cerr << funame << "should not be called for monoatomic molecule\n";

      throw Error::Logic();
    }

    return _ang_vel[frag];
  }

  /*************************************************************************
   *************************** Dynamic Variables ***************************
   *************************************************************************/

  class Vars : public Coordinates, public Momenta
  {
  public:
    //
    static int size () { return Structure::dv_size(); }

    void get (const double* dv, double* =0);
    void put (double* dv) const { Coordinates::put(dv); Momenta::put(dv + Coordinates::size()); }

    Vars () {}
    Vars (const double* dv) { get(dv); }
    Vars (const Coordinates& dc) : Coordinates(dc) {}

    void    total_angular_momentum (double*)      const;
    void  orbital_angular_momentum (double*)      const;
    void fragment_angular_momentum (int, double*) const;

    D3::Vector    total_angular_momentum  ()    const;
    D3::Vector  orbital_angular_momentum  ()    const;
    D3::Vector fragment_angular_momentum  (int) const;

    double radial_kinetic_energy         ()    const;
    double total_angular_momentum_length ()    const;
    double angular_momentum_k_projection (int) const;
    double angular_momentum_m_projection (int) const;
  };

  inline void Vars::get (const double* dv, double* len) 
  {
    Coordinates::get(dv, len); 

    Momenta::get(dv + Coordinates::size()); 

    for(int f = 0; f < 2; ++f)
      //
      if(Structure::type(f) == Molecule::LINEAR)
	//
	::orthogonalize(ang_vel(f), ang_pos(f), 3);
  }

  inline D3::Vector Vars::total_angular_momentum () const 
  { 
    D3::Vector res; 
    total_angular_momentum(res); 
    return res;
  }

  inline D3::Vector Vars::orbital_angular_momentum () const 
  { 
    D3::Vector res; 
    orbital_angular_momentum(res); 
    return res;
  }

  inline D3::Vector Vars::fragment_angular_momentum (int frag) const 
  { 
    D3::Vector res; 
    fragment_angular_momentum(frag, res); 
    return res;
  }

  /*************************************************************************
   ******************************** CLASSIFIER *****************************
   *************************************************************************/

  // abstract class to classify the configuration
  class Classifier {
  public:
    virtual int  classify (const Coordinates&)  const =0;
    virtual ~Classifier () {}
  };

  /*************************************************************************
   ******************************* CONDITION *******************************
   *************************************************************************/

  // abstract logical condition
  //
  class Condition {
  public:
    //
    virtual bool test (const Coordinates&) const =0;

    virtual ~Condition () {}

    virtual void print (std::ostream&, const std::string&) const {}
  };

  typedef ConstSharedPointer<Condition> CCP;

  // configuration region which is not sampled 
  //
  extern CCP exclude_region;

  //distance condition
  //
  class DistanceCondition : public Condition {
    //
    double _dist;

  public:
    //
    DistanceCondition (double d) : _dist(d) {}

    ~DistanceCondition () {}

    bool test (const Coordinates&) const;
  };

  inline bool DistanceCondition::test (const Coordinates& dc) const 
  {
    if(::vlength(dc.orb_pos(), 3) > _dist)
      //
      return true;

    return false;
  }

  // binary condition
  //
  class BinaryCondition : public Condition {
    //
    Logical::BinExpr::Oper _op;
    CCP _c1, _c2;

    BinaryCondition(CCP c1, CCP c2, Logical::BinExpr::Oper o) : _op(o), _c1(c1) , _c2(c2) {}

  public:
    //
    bool test (const Coordinates&) const;

    ~BinaryCondition () {}

    friend CCP operator& (CCP c1, CCP c2); 
    friend CCP operator| (CCP c1, CCP c2); 
  };

  inline CCP operator& (CCP c1, CCP c2) { return CCP(new BinaryCondition(c1, c2, Logical::BinExpr::AND)); }
  inline CCP operator| (CCP c1, CCP c2) { return CCP(new BinaryCondition(c1, c2, Logical::BinExpr::OR)); }

  inline CCP& operator&= (CCP& c1, CCP c2) { c1 = c1 & c2; return c1; }
  inline CCP& operator|= (CCP& c1, CCP c2) { c1 = c1 | c2; return c1; }


  inline bool BinaryCondition::test (const Coordinates& dc) const
  {
    const char funame [] = "Dynamic::BinaryCondition::test: ";

    switch(_op) {
    //
    case Logical::BinExpr::AND: 
      //
      return  _c1->test(dc) && _c2->test(dc);

    case Logical::BinExpr::OR: 
      //
      return  _c1->test(dc) || _c2->test(dc);

    default:
      //
      std::cerr << funame << "unknown operation\n";

      throw Error::Logic();
    }
  }

  // negation condition
  //
  class UnaryCondition : public Condition
  {
    CCP _c;

    UnaryCondition(CCP c) : _c(c) {}

  public:
    //
    bool test (const Coordinates& dc) const { return !_c->test(dc); }

    ~UnaryCondition () {}

    friend CCP negate (CCP c);
  };

  inline CCP negate (CCP c) { return CCP(new UnaryCondition(c)); }

  /*************************************************************************
   ********************* CONFIGURATIONAL OBSERVABLE ************************
   *************************************************************************/

  class ObservBase {
    //
  public:
    //
    virtual double evaluate (const Coordinates&) const =0;
  };

  // atom-atom distance, atom-atom-atom angle, or atom-atom-atom-atom dihedral
  //
  class ZmatObserv : public ObservBase {
    //
    std::vector<int> _atom;

    static D3::Vector _vector(int, int, const Coordinates&, const std::vector<D3::Vector> []);

  public:
    //
    explicit ZmatObserv (const std::vector<int>&);
    
    double evaluate (const Coordinates&) const;
  };

  class ObservCondition : public Condition {
    ConstSharedPointer<ObservBase> _observ;
    std::pair<double, double>      _range;
    int                            _mode;

  public:
    enum {LOWER_BOUND, 
	  UPPER_BOUND, 
	  RANGE_BOUND
    };

    ObservCondition (ConstSharedPointer<ObservBase> p, double limit, int m);

    ObservCondition (ConstSharedPointer<ObservBase> p, double lower, double upper);

    bool test (const Coordinates&) const;
  };


  /*************************************************************************
   ********************************** WATCH ********************************
   *************************************************************************/

  // watch if some condition has changed along the trajectory
  class Watch
  {
    CCP _cond;
    bool _start;
    bool _changed;

  public:
    Watch (CCP c, const Coordinates& dv) : _cond(c), _start(c->test(dv)), _changed(false) {}
    Watch (CCP c, bool s) : _cond(c), _start(s), _changed(false) {}

    void test (const Coordinates&);
    bool is_changed () const { return _changed; }
  };

  inline void Watch::test (const Coordinates& c) 
  { 
    if(_changed) 
      return;
    if(_cond->test(c) != _start)
      _changed = true;
  }

  /*************************************************************************
   ********************************** ACTION *******************************
   *************************************************************************/

  class Action
  {
  public:
    virtual void execute (const Vars&)  =0;
  };

}// Dynamic

#endif
