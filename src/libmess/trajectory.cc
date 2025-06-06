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
  //
  // default integration step
  //
  double default_step = 100.;

  std::vector<double> traj_rel_tol; // trajectory relative tolerance
  std::vector<double> traj_abs_tol; // trajectory absolute tolerance

}

// dvd calculator
//
extern "C" void Trajectory::set_dvd (const double& time, const double* dv, double* dvd, void*, void* par) 
{
  DvdPar& dvd_par = *static_cast<DvdPar*>(par);

 /****************************************************
  * Comulative force acting on a second fragment     *
  * and Torques on each fragment: This can be either *
  * laboratory frame torque for a linear fragment or *
  * molecular  frame torque for a nonlinear fragment *
  ****************************************************/

  D3::Vector torque[3]; 
  Dynamic::Coordinates dc(dv);

  //dvd_par.ener = dvd_par.pot(dc, torque);

  // long jump facility seems to interfere with OpenMP parallelization
  //
  try {
    //
    // energies, forces and torques

    dvd_par.ener = dvd_par.pot(dc, torque);
  }
  catch (Error::General) {
    //
    longjmp(dvd_par.jmp, 1);
  }

  // dynamic variables derivatives
  //
  Dynamic::set_dvd(torque, dv, dvd);
}

// propagate trajectory
//
void Trajectory::Propagator::run (const Dynamic::CCP& stop, const Dynamic::Classifier& sort, Dynamic::CCP special_stop) 
{
  const char funame [] = "Trajectory::Propagator::run: ";

  //std::cerr << funame << "relative tolerance = " << rel_tol[0] << "\n";

  double dtemp;
  int itemp;
  D3::Vector vtemp;
	
  // Long jump facility to transfer failure information  between the dvd calculator and calling function
  //
  DvdPar dvd_par(_pot);

  // longjmp facility seems to inerfere with OpenMP parallelization
  //
  if(setjmp(dvd_par.jmp)) {
    //
    std::cerr << funame << "potential calculation failed\n"; 

    throw PotentialFailure();
  }

  if(_out) {
    //
    _out << std::setw(15) << "time, au" 
	//
	 << std::setw(15) << "dist, bohr"
	//
	 << std::setw(15) << "ener, kcal"
	//
	 << std::setw(15) << "amom, au"
	//
	 << "\n";
  }

  // dynamic variables
  //
  ::Array<double> dv(size());

  put(dv);
  Mode mode = RESTART;

  double timeout;
  int adjust_count = 0;

  // main cycle
  //
  while(1) {
    //
    switch(_dir) {
      //
    case FORWARD:
      //
      timeout = _time + _step;

      break;
      //
    case BACKWARD:
      //
      timeout = _time - _step;

      break;
      //
    default:
      //
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }

    try {
      //
      AdamSolver::run(_time, dv, timeout, static_cast<void*>(&dvd_par), mode);
    }
    catch(Error::General) {
      //
      throw RunFailure();
    }

    // get dynamical variables data and normalize
    //
    double ang_len [2];
    get(dv, ang_len);

    /********************************************************************
     * checking orthogonality of angular velocity to the molecular axis *
     * for linear fragments  and normalization of the angular vectors   *
     * for all nonatomic fragments                                      *
    **********************************************************************/

    bool need_adjustment = false;

    for(int frag = 0; frag < 2; ++frag) {
      //
      if(Structure::type(frag) == Molecule::MONOATOMIC)
	//
	continue;

      if(ang_len[frag] > ang_len_tol || ang_len[frag] < 1. / ang_len_tol) {
	//
	std::cerr << funame << "WARNING: " << frag << "-th fragment: angular vector not normalized: " << ang_len[frag] << "\n";

	need_adjustment = true;
      }

      if(Structure::type(frag) == Molecule::LINEAR) {
	//
	// velocity projection to molecular axis
	//
	dtemp = vdot(ang_pos(frag), ang_vel(frag), 3);
	dtemp = dtemp > 0. ? dtemp : -dtemp;

	double avl = ::vlength(ang_vel(frag), 3);

	// checking if velocity projection on the molecular axis does exceed the calculation error
	//
	itemp = Structure::pos_size() + Structure::ang_vel(frag);

	if(dtemp > abs_tol[itemp] + rel_tol[itemp] * avl) {
	  //
	  std::cerr << funame << "WARNING: " << frag << "-th fragment: angular velocity is not orthogonal to molecular axis, angle = "
	    //
		    << dtemp / avl << " rad, adjusting " << ++adjust_count << " time\n";

	  orthogonalize(ang_vel(frag), ang_pos(frag), 3);

	  need_adjustment = true;
	}
      }
    }

    if(need_adjustment) { 
      //
      put(dv);

      mode = RESTART;
    }
    else
      //
      mode = CONTINUE;

    // intermediate output
    //
    if(_out) {
      //
      _out << std::setw(15) << _time
	//
	   << std::setw(15) << orb_len()
	//
	   << std::setw(15) << (current_energy() - _start_ener) / Phys_const::kcal
	//
	   << std::setw(15) << (total_angular_momentum() - _start_amom).vlength()
	//
	   << "\n";
    }

    // stop condition
    //
    if(special_stop && !special_stop->test(*this) || stop->test(*this) || Dynamic::exclude_region && Dynamic::exclude_region->test(*this)) {
      //
      // remove colinear velocity component
      //
      for(int f = 0; f < 2; ++f)
	//
	if(Structure::type(f) == Molecule::LINEAR)
	  //
	  orthogonalize(ang_vel(f), ang_pos(f), 3);

      _spec = sort.classify(*this);

      _ener = current_energy();

      return;
    }
  }
}
