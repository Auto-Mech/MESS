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

#include "trajectory.hh"
#include "units.hh"

namespace Trajectory {
  double Propagator::step = 100;
  //Flags Propagator::flags;
  //Dynamic::CCP fail_condition;
}

// dvd calculator
extern "C" void Trajectory::set_dvd (const double& time, const double* dv, double* dvd, 
				     void*, void* par) 
{
  DvdPar& dvd_par = *static_cast<DvdPar*>(par);

  // Comulative force acting on a second fragment
  // and Torques on each fragment: This can be either 
  // laboratory frame torque for a linear fragment or
  // molecular  frame torque for a nonlinear fragment
  D3::Vector torque [3]; 
  Dynamic::Coordinates dc(dv);

  try {
      // energies, forces and torques
      dvd_par.ener = dvd_par.pot(dc, torque);
  }
  catch (Error::General) {
    longjmp(dvd_par.jmp, 1);
  }
  // dynamic variables derivatives
  Dynamic::set_dvd(torque, dv, dvd);

  // something else
  // ...
}

void Trajectory::Propagator::run (Dynamic::CCP stop, const Dynamic::Classifier& sort) 
{
  const char funame [] = "Trajectory::Propagator::run: ";

  double dtemp;
  int itemp;
  D3::Vector vtemp;
	
  // Long jump facility to transfer failure information 
  // between the dvd calculator and calling function
  DvdPar dvd_par;
  dvd_par.pot = _pot;

  if(setjmp(dvd_par.jmp)) {
    //std::cerr << funame << "potential calculation failed\n"; 
    throw PotentialFailure();
  }

  // dynamic variables
  //
  ::Array<double> dv(size());

  put(dv);
  Mode mode = RESTART;

  double timeout;
  int adjust_count = 0;
  while(1) {// main cycle

    switch(_dir) {
    case FORWARD:
      timeout = _time + step;
      break;
    case BACKWARD:
      timeout = _time - step;
      break;
    default:
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }

    try {
      AdamSolver::run(_time, dv, timeout, static_cast<void*>(&dvd_par), mode);
    }
    catch(Error::General) {
      throw RunFailure();
    }

    // get dynamical variables data and normalize
    get(dv);

    // checking orthogonality of angular velocity to the molecular axis 
    // for linear fragments  and normalization of the angular vectors
    // for all nonatomic fragments

    bool need_adjustment = false;
    for(int frag = 0; frag < 2; ++frag) {// fragment cycle
      if(Structure::fragment(frag).type() == Molecule::MONOATOMIC)
	continue;

      if(length(frag) > 2. || length(frag) < 0.5) {
	std::cerr << funame << "WARNING: length of " << frag << "-th fragment is not normalized\n";
	need_adjustment = true;
      }

      if(Structure::fragment(frag).type() == Molecule::LINEAR) {
	// velocity projection
	dtemp = vdot(ang_pos(frag), ang_vel(frag), 3);
	dtemp = dtemp > 0. ? dtemp : -dtemp;

	// checking if velocity projection on the molecular axis does exceed the calculation error
	itemp = Structure::pos_size() + Structure::ang_vel(frag);
	double avl = vlength(ang_vel(frag), 3);
	if(dtemp > abs_tol[itemp] + rel_tol[itemp] * avl) {
	  std::cerr << funame << "WARNING: angular velocity of the " << frag 
		    << "-th fragment is not orthogonal, angle = "
		    << dtemp / avl / length(frag) << " rad, adjusting " << ++adjust_count << " time\n";

	  orthogonalize(ang_vel(frag), ang_pos(frag), 3);
	  need_adjustment = true;
	}
      }
      // ...
    }// fragment cycle

    if(need_adjustment) { 
      put(dv);
      mode = RESTART;
    }
    else
      mode = CONTINUE;

    /*	
    // do some output with the flags
    // ...
	
    // run watch tests
    for(int i = 0; i < watch.size(); ++i)
    watch[i].test(*this);

    // execute registered actions
    for(int i = 0; i < act.size(); ++i)
    act[i]->execute(*this);
    */

    // stop condition
    if(stop->test(*this)) {
      _spec = sort.classify(*this);
      _ener = total_kinetic_energy() + _pot(*this);
      return;
    }

    // exclude region
    if(Dynamic::exclude_region && Dynamic::exclude_region->test(*this))
      throw ExcludeRegionHit();
  }
}
