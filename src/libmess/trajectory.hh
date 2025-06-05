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

#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH

#include "dynamic.hh"
#include "slatec.hh"
#include "potential.hh"

#include <setjmp.h>

enum {BACKWARD = 0, FORWARD = 1};

namespace Trajectory {
  //
  class Exception {};
  class PotentialFailure : public Exception {};
  class       RunFailure : public Exception {};

  // dynamic variables derivatives calculator
  //
  extern "C" void set_dvd (const double& time, const double* dv, double* dvd, void*, void* par);

  extern std::vector<double> traj_rel_tol; // trajectory relative tolerance
  extern std::vector<double> traj_abs_tol; // trajectory absolute tolerance

  // default integration step
  //
  extern double default_step;

  // container to transfer data between dvd calculator and the calling function
  //
  struct DvdPar {
    //
    const Potential::Wrap& pot;
    double ener; // energy
    jmp_buf jmp; // jump buffer for for long jumps

    DvdPar (const Potential::Wrap& p) : pot(p) {}
  };

  class Out {
    //
    std::ostream* _out;

    void _assert () const
    {
      if(!_out) {
	//
	std::cerr << "Trajectory::Out:_assert: not initialized\n";

	 throw Error::Init();
      }
    }

  public:
    //
    Out () : _out(0) {}

    Out (std::ostream& out) : _out(&out) {}

    operator bool () { return _out; }

    template<typename T>
    Out& operator<< (const T& t) { _assert(); *_out << t; return *this; }

    void init (std::ostream& out) 
    { 
      if(_out) { 
	//
	std::cerr << "Trajectory::Out::init: already initialized\n";

	throw Error::Init();
      }
      _out = &out;
    }

    int precision (int p) { _assert(); return _out->precision(p); }
  };

  class Propagator : public Dynamic::Vars, public Slatec::AdamSolver
  {
    int    _dir;
    int    _spec;
    double _time;
    double _ener;
    double _step;
    
    double     _start_ener;
    D3::Vector _start_amom;

    const Potential::Wrap& _pot;

    Out _out;

    int _old_precision;

  public:
    //
    void set_step (double s) { _step = s; }

    double time         () const { return _time; }
    int    direction    () const { return _dir; }

    double current_energy () const { return total_kinetic_energy() + _pot(*this); }
    double   final_energy () const { return _ener; }

    int    species      () const { return _spec; }

    void setout (std::ostream&, int = 3);

    ~Propagator () { if(_out) _out.precision(_old_precision); }

    Propagator(const Potential::Wrap& pot, const Dynamic::Vars& dv, int dir) 
	//
      : _pot(pot), Dynamic::Vars(dv), Slatec::AdamSolver(Structure::dv_size(), set_dvd),
	//
	_time(0.0), _dir(dir), _step(default_step)
      {
	rel_tol = traj_rel_tol;
	abs_tol = traj_abs_tol;

	_start_ener = current_energy();
	_start_amom = total_angular_momentum();
      }

    Propagator(const Potential::Wrap& pot, const Dynamic::Vars& dv, int dir, std::ostream& out, int = 3);

    Propagator(const Potential::Wrap& pot, const Dynamic::Vars& dv, const Slatec::AdamSolver& as, int dir) 
      //
      : _pot(pot), Dynamic::Vars(dv), Slatec::AdamSolver(as), _time(0.0), _dir(dir), _step(default_step)
      {
	rel_tol = traj_rel_tol;
	abs_tol = traj_abs_tol;

	_start_ener = current_energy();
	_start_amom = total_angular_momentum();
      }

    Propagator(const Potential::Wrap& pot, const Dynamic::Vars& dv, int dir, double rt, double at) 
      //
      : _pot(pot), Dynamic::Vars(dv), Slatec::AdamSolver(Structure::dv_size(), set_dvd, rt, at), 
	//
	_time(0.0), _dir(dir), _step(default_step) 
      {
	_start_ener = current_energy();
	_start_amom = total_angular_momentum();
      }

    void run (const Dynamic::CCP& stop, const Dynamic::Classifier& sort, Dynamic::CCP = Dynamic::CCP());
  };

  inline void Propagator::setout (std::ostream& out, int prec) { _out.init(out); _old_precision = _out.precision(prec); }

  inline Propagator::Propagator(const Potential::Wrap& pot, const Dynamic::Vars& dv, int dir, std::ostream& out, int prec) 
	//
      : _pot(pot), Dynamic::Vars(dv), Slatec::AdamSolver(Structure::dv_size(), set_dvd),
	//
	_time(0.0), _dir(dir), _step(default_step), _out(out)
    {
      rel_tol = traj_rel_tol;
      abs_tol = traj_abs_tol;

      _start_ener = current_energy();
      _start_amom = total_angular_momentum();

      _old_precision = _out.precision(prec);
    }
}

#endif
