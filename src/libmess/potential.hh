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

#ifndef POTENTIAL_HH
#define POTENTIAL_HH

#include "dynamic.hh"
#include "io.hh"
#include "system.hh"
#include "dynlib.hh"
#include "slatec.hh"

#include <omp.h>

namespace Potential {

  // potential exceptions
  //
  class Exception {};

  // potential types
  //
  enum {
    TN,            // Torque-Numeric 
    CL,            // Charge - linear molecule multipole potential
    CN,            // Charge - non-linear molecule multipole potential
    DD,            // Dipole -dipole multipole potential
    MULTIPOLE,     // General long-range multipole potential
    HARMONIC       // harmonic expansion
  };

  // abstract class base for potential objects
  //
  class Base { 
    //      
    mutable std::vector<long> _count;

  protected:
    //
    void _increment_count () const { _count[omp_get_thread_num()]++; }

  public:
    //
    Base ();

    virtual ~Base () {}

    virtual double operator() (const Dynamic::Coordinates&, D3::Vector*) const =0;

    virtual int type () const =0;

    long reset_count () const { long res = 0; for(int i = 0; i < _count.size(); ++i) { res += _count[i]; _count[i] = 0; } return res; }
  };

  inline Base::Base ()
  {
    const char funame [] = "Potential::Base::Base: ";

#pragma omp parallel
#pragma omp master
    {
      IO::log << IO::log_offset << funame << "number of threads = " << omp_get_num_threads() << "\n\n";
 
      _count.resize(omp_get_num_threads());
    }
  }

  // potential object wrapper
  //
  class Wrap : public IO::Read 
  {
    ConstSharedPointer<Base> _fun;
	
  public:
    //
    Wrap () {}
    Wrap (ConstSharedPointer<Base> pf) : _fun(pf) {}
    ~Wrap () {}

    void read (std::istream&);

    operator  bool () const { return  _fun; }
    bool operator! () const { return !_fun; }

    void   isinit () const ;

    double operator() (const Dynamic::Coordinates&, D3::Vector* =0) const;
    int type () const;
    long reset_count () const;
  };

  inline void Wrap::isinit () const 
  {
    const char funame [] = "Potential::wrap::isinit: ";

    if(!_fun) {
      //
      std::cerr << funame << "not initialized\n";

      throw Error::Init();
    }
  }

  inline double Wrap::operator() (const Dynamic::Coordinates& dc, D3::Vector* force) const 
  {
    isinit();
    return (*_fun)(dc, force);
  }

  inline int Wrap::type () const 
  {
    isinit();
    return _fun->type();
  }
  
  inline long Wrap::reset_count () const 
  {
    isinit();
    return _fun->reset_count();
  }
  
  // low energy condition
  //
  class Condition : public Dynamic::Condition 
  {
    Potential::Wrap _pot;
    double _emin;

  public:
    //
    Condition (Potential::Wrap pot, double e);

    bool test (const Dynamic::Coordinates&) const;

    ~Condition () {}
  };

  inline Condition::Condition (Potential::Wrap pot, double e) : _pot(pot), _emin(e)
  {
    const char funame [] = "Potential::Condition::Condition: ";

    if(!_pot) {
      //
      std::cerr << funame << "potential not initialized\n";

      throw Error::Init();
    }
  }

  inline bool Condition::test (const Dynamic::Coordinates& dc) const
  {
    if(_pot(dc) < _emin) {
      //
      return true;
    }
    else
      //
      return false;
  }

  /********************************************************************************************
   ********************** TORQUE BY NUMERICAL DIFFERENTIATION FACILITY ************************
   ********************************************************************************************/
  
  // base bare energy (no torque) potential class
  //
  class BareEnergy {
    //
  public:
    //
    virtual double evaluate (const Dynamic::Coordinates&) const =0;

    virtual ~BareEnergy () {}
  };

  // Analytic sjk-style bare energy potential
  //
  extern "C" {
    //
    typedef double (*ener_t) (const double* r, const double* rpar, const int* ipar, int& ifail);
    typedef double (*sjk_t)  (const double* r, const double* rpar, const int* ipar);
    typedef void   (*init_t) (const char* data_file);
  }

  class SJK : public BareEnergy {
    //
    enum {SJK_FORM};
    
    int             _pot_form;
    
    int _format () const { return _pot_form; }
    
    DynLib	    _pot_libr;
    void*           _pot_ener;
    init_t          _pot_init;
    Array<double>   _pot_rpar;
    Array<int>      _pot_ipar;

    DynLib	   _corr_libr;
    void*          _corr_ener;
    init_t         _corr_init;
    Array<double>  _corr_rpar;
    Array<int>     _corr_ipar;

    // no copies
    //
    SJK (const SJK&);
    SJK& operator= (const SJK&);
    
  public:
    //
    SJK (std::istream&);
    ~SJK () {}

    double evaluate (const Dynamic::Coordinates&) const;
  };

  class TorqueNumeric : public Base
  {
    ConstSharedPointer<BareEnergy> _bare_ener;

    // cartesian & angular increments for numerical differentiation
    //
    double _dist_incr, _ang_incr;

    // get torques by numerical differentiation
    //
    void _set_torque (const Dynamic::Coordinates&, D3::Vector*) const;

  public:
    //
    TorqueNumeric (std::istream&);
    ~TorqueNumeric () {}
    
    double operator() (const Dynamic::Coordinates&, D3::Vector*) const;

    int type () const { return TN; }
  };

  /********************************************************************************************
   ************************************** YG POTENTIAL ****************************************
   ********************************************************************************************/
  namespace YG {
    //
    inline void assert ()
    {
      const char funame [] = "Potential::YG::assert: ";

      if(Structure::type(0) == Molecule::NONLINEAR)
	//
	return;

      std::cerr << funame << "the first frament of the potential should be non-linear\n";
      
      throw Error::Init();
    }

    // potential base class
    //
    class Base : public BareEnergy {
      //
    protected:
      //
      typedef std::map<std::vector<int>, Slatec::Spline> coef_t;

      coef_t _coef;

      double _dist_min, _dist_max;
      
    public:
      //
      int size () const { return _coef.size(); }
      
      virtual double basis_function (const std::vector<int>&, const Array<double>&) const =0;

      virtual std::vector<int> index_range () const =0;
      
      void read_coefficients (const std::string&);
      
      double evaluate (const Dynamic::Coordinates&) const;

      virtual ~Base () {}
    };

    // Anchor fragment with Cs symmetry + Atom
    //
    class Cs_Atom : public Base {
      //
      int _x;
      
    public:
      //
      double basis_function (const std::vector<int>&, const Array<double>&) const;

      std::vector<int> index_range () const { return std::vector<int>(2, -1); }

      Cs_Atom (int);
      
      ~Cs_Atom () {}
    };

    // Anchor fragment with C1 symmetry + Atom
    //
    class C1_Atom : public Base {
      //
      int _x;
      
    public:
      //
      double basis_function (const std::vector<int>&, const Array<double>&) const;

      std::vector<int> index_range () const { std::vector<int> res(3, -1); res[2] = 2; return res; }

      C1_Atom (int);
      
      ~C1_Atom () {}
    };

    extern ConstSharedPointer<Base> pot;
    
    void init (std::istream&);
    //
  }// YG namespace

  /********************************************************************************
   **************************** MULTIPOLE POTENTIALS ******************************
   ********************************************************************************/
  
  // Multipole potential for a charge(1) and a linear(2) molecule
  //
  class ChargeLinear : public Base
  {
    double                     _charge;
    double                     _dipole;
    double                 _quadrupole;  // occording to Landau & Lifshitz definition
    double   _isotropic_polarizability;
    double _anisotropic_polarizability;

    enum _mode_t {POT_VALUE, DIST_DERIV, ANGLE_DERIV};
    double _pot (double dist, double cos_theta, _mode_t mode) const;

    // no copies
    //
    ChargeLinear (const ChargeLinear&); 
    ChargeLinear& operator= (const ChargeLinear&);

  public:
    //
    ChargeLinear (std::istream&);
    ~ChargeLinear () {}

    double operator() (const Dynamic::Coordinates&, D3::Vector*) const ;
    int type () const { return CL; }
  };

  // Multipole potential for a charge(1) and a nonlinear(2) molecule
  //
  class ChargeNonlinear : public Base
  {
    double             _charge;
    D3::Vector         _dipole; // dx, dy, dz
    D3::Matrix     _quadrupole; // qxx, qxy, qxz, qyy, qyz
    D3::Matrix _polarizability; // pxx, pxy, pxz, pyy, pyz, pzz

    // no copies
    //
    ChargeNonlinear (const ChargeNonlinear&);
    ChargeNonlinear& operator= (const ChargeNonlinear&);

  public:
    //
    ChargeNonlinear (std::istream&) ;
    ~ChargeNonlinear () {}

    double operator() (const Dynamic::Coordinates&, D3::Vector*) const ;
    int type () const { return CN; }
  };


  // Dipole-dipole potential
  //
  class DipoleDipole : public Base
  {
    D3::Vector _dipole [2];
    double _polarizability [2] ; // isotropic polarizability is assumed
    double _dispersion;  // dispersion interaction

    // no copies
    //
    DipoleDipole (const DipoleDipole&);
    DipoleDipole& operator= (const DipoleDipole&);

  public:
    //
    DipoleDipole (std::istream&);
    ~DipoleDipole () {}

    double operator() (const Dynamic::Coordinates& dc, D3::Vector*) const;
    int type () const { return DD; }
  };


  // Multipole potential
  //
  class Multipole : public Base
  {
    double _dispersion;  // dispersion interaction

    // no copies
    //
    Multipole (const Multipole&);
    Multipole& operator= (const Multipole&);

  public:
    //
    Multipole (std::istream&);
    ~Multipole () {}

    double operator() (const Dynamic::Coordinates& dc, D3::Vector*) const;
    int type () const { return MULTIPOLE; }
  };

}// Potential namespace

extern Potential::Wrap default_pot;

#endif
