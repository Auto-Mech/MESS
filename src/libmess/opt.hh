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

#ifndef OPT_HH
#define OPT_HH

#include "coord.hh"
#include "shared.hh"
#include "harding.hh"
#include "linpack.hh"
#include "lapack.hh"

// molecular geometry optimization

namespace Opt {

  using namespace Coord;

  /********************* POTENTIALS ************************/
  
  // wrapper for Harding potential
  //
  class HardPot : public CartFun, private Harding {
    //
  public:
    //
    HardPot(std::istream& from) : Harding(from) {}
    
    double evaluate (const Cartesian&) const;

    int atom_size () const { return Harding::atom_size(); }
  };
    
  /********************* OPTIMIZATION METHODS ************************/
  
  // constrained optimization in cartesian configurational space
  //
  class CcOpt {
    //
    // potential
    //
    ConstSharedPointer<CartFun> _pot;

    // set of constrains
    //
    typedef std::set<ConstSharedPointer<CartFun> > _con_t;
    
    _con_t _constrain;

    // energy gradient tolerance
    //
    double _grad_tol;

    // initial time step
    //
    double _time_step;

    // integration time
    //
    double _time_length;

    // (local) independence tolerance of constrains
    //
    double _ci_tol;

    // gradient calculation count
    //
    static int _grad_count;

    // end of integration exception
    //
    class _Fin {};
    
  public:
    //
    mutable int use_constrain;
    
    CcOpt () : _grad_tol(1.e-5), _time_step(0.1), _time_length(1.e5), _ci_tol(1.e-5), use_constrain(1) {}

    void set (std::istream&);
    
    // potential
    //
    ConstSharedPointer<CartFun> pot () const { return _pot; }

    void set_pot (ConstSharedPointer<CartFun> p);
    
    void add_constrain (ConstSharedPointer<CartFun> c);

    // RHS of differential equation calculator
    //
    void operator () (const Cartesian&, Cartesian&, double = 0.) const;

    // find a minimum
    //
    void execute (Cartesian&) const;
  };

  inline void CcOpt::add_constrain (ConstSharedPointer<CartFun> c)
  {
    const char funame [] = "Opt::CcOpt::add_constrain: ";

    if(!c) {
      //
      std::cerr << funame << "empty pointer\n";

      throw Error::Range();
    }
      
    if(!_constrain.insert(c).second) {
      //
      std::cerr << funame << "duplicated constrain\n";
    }
  }

  inline void CcOpt::set_pot (ConstSharedPointer<CartFun> p)
  {
    const char funame [] = "Opt::CcOpt::set_pot: ";

    if(!p) {
      //
      std::cerr << funame << "empty pointer\n";

      throw Error::Range();
    }
      
    if(_pot) {
      //
      std::cerr << funame << "already initialized\n";

      throw Error::Init();
    }

    _pot = p;
  }

  /********************************************************************************************************
   *********************** OPTIMIZATION AND IMPORTANCE SAMPLING IN Z-MATRIX COORDINATES********************
   ********************************************************************************************************/

  class ZOpt {
    //
  public:
    //
    typedef std::map<std::string, std::pair<double, double> > mode_t;

  private:
    //
    // conserved (non-fluxional) modes and their sampling limits
    //
    mode_t _con_modes;
    
    // z-matrix coordinates of the constrained minimum
    //
    ZData _zmin;

    // energy at the constrained minimum
    //
    double _ener_min;

    // all modes generalized mass factor (including inertia moments factor)  at the constrained minimum
    //
    double _mass_factor;
    
    // non-fluxional (conserved) modes force constant matrix eigenvalues square roots at the constrained minimum
    //
    Lapack::Vector _fc_eval_sqrt;

    // non-fluxional (conserved) modes force constant matrix eigenvectors at the constrained minimum
    //
    Lapack::Matrix _fc_evec;

    // potential surface for non-fluxional (conserved) modes
    //
    double _con_pot (const std::vector<double>&) const;
    
    // potential gradient for non-fluxional (conserved) modes
    //
    Lapack::Vector          _con_grad () const;

    // hessian for non-fluxional (conserved) modes
    //
    Lapack::SymmetricMatrix _con_hess () const;

  public:
    //
    ZOpt (const ZData&, const mode_t&);
 
    // partition function unharmonic correction by importance sampling 
    //
    double anharmonic_correction (double, int, double* =0, int* =0) const;

    double mass_factor () const { return _mass_factor; }

    // non-fluxional modes force constant factor
    //
    double fc_factor () const { return product(_fc_eval_sqrt); }

    double ener_min () const { return _ener_min; }

    // z-matrix coordinates of the constrained minimum
    //
    const ZData& zmin () const { return _zmin; }

    // potential
    //
    static ConstSharedPointer<CartFun> pot;

    // gradient tolerance
    //
    static double grad_tol;

    // maximal optimization step
    //
    static double max_opt_step;

    // maximal number of optimization iterations
    //
    static int max_opt_count;

    // diffirentiation step
    //
    static double diff_step;

    // no convergence exception
    //
    class NoConv {};
  };
}
#endif
