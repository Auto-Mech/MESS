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

#include "slatec.h"

/****************************************************************************************************
 *                   Differential equations solver by the Adams-Bashforth-Moulton                   *
 *                    Predictor-Corrector formulas of orders one through twelve                     *
 ****************************************************************************************************/

const int Slatec::AdamSolver::_liw;

Slatec::AdamSolver::AdamSolver (int s, dde_f d, double rt, double at)   
  : _deriv(d), _size(s),  _rwork(130 + 21 * s),  rel_tol(s, rt), abs_tol(s, at)
{
  const char funame [] = "AdamSolver::AdamSolver: ";

  if(!s) {
    std::cerr << funame << "zero dimension\n";
    throw Error::Logic();
  }

  if(!_deriv) {
    std::cerr << funame << "the function to calculate derivatives should be provided\n";
    throw Error::Logic();
  }

  _info [0] = 0;     // start new trajectory
  _info [1] = 1;     // rel_tol and abs_tol are vectors
  _info [2] = 0;     // the solution only at the end point
  _info [3] = 0;     // the integration can go beyond final point

  _lrw = _rwork.size();
} 
    
void Slatec::AdamSolver::run(double& x, double* y, double xout,
			     void* param, Mode mode) 
{
  const char funame [] = "AdamSolver::run: ";

  const int count_max = 5;// number of attempts to run an integrator

#ifdef DEBUG
  if(rel_tol.size() != _size || abs_tol.size() != _size) {
      std::cerr << funame << "error tolerance: wrong size\n";
      throw Error::Logic();
  }

  for(int i = 0; i < _size; ++i)
    if(rel_tol[i] < 0. || abs_tol[i] <= 0.) {
      std::cerr << funame << "error tolerance: negative\n";
      throw Error::Logic();
    }
#endif 
	  
  if(mode == RESTART)
    _info[0] = 0;

  int_t idid;
  for(int count = 0; count < count_max; ++count) {
      if(count)
	  std::cerr << funame << count + 1 << "-th etempt ... ";
    ddeabm_(_deriv, _size, x, y, xout, _info, rel_tol, abs_tol, idid, 
	    _rwork, _lrw, _iwork, _liw, 0, param);

    _info[0] = 1;

    // check for errors;
    if(idid < 0) {
	if(count)
	    std::cerr << "failure:\n";
      switch(idid) {
      case -1:
	std::cerr << funame << "DDEABM has attempted 500 steps:"
	  " relaxing the error tolerance by 2 and trying again\n";
	rel_tol *= 2.;
	abs_tol *= 2.;
	break;

      case -2:
	std::cerr << funame << "the error tolerances have been increased by DDEABM;"
	  " trying again\n";
	break;
	    
      case -3:
	std::cerr << funame << "one of the absolute error tolerances is zero?:\n";
	for(int i = 0; i < _size; ++i)
	  std::cerr << "component = " << i 
		    << ", absolute error tolerance = " << abs_tol[i] << "\n";
	throw Error::Logic();

      case -4: 
	std::cerr << funame << "the problem appears to be stiff:"
	  "one may try using DDEBDF integrator instead of DDEABM;"
	  " meanwhile relaxing the error tolerance by 2 and trying again\n";
	rel_tol *= 2.;
	abs_tol *= 2.;
	break;

      default:
	std::cerr << funame << "integrator failed; error code idid = " 
		  << idid << "\n";
	throw Error::Run();
      }
    }
    else {
	if(count)
	    std::cerr << "success!!!\n";
	break;
    }
  }
  if(idid < 0) {
      std::cerr << funame << "maximum number of steps have been reached;"
      "inegrator failed with the code idid = " << idid << "\n";
    throw Error::Run();
  }
}


/**********************************************************************************
 *                                  Spline fitting                                *
 *********************************************************************************/

void Slatec::Spline::init (const double* x, const double* y, 
		       int_t s) 
{
    const char funame [] = "Slatec::Spline::init: ";

    if(_size) {
      std::cerr << funame << "already initialized\n";
      throw Error::Init();
    }
    
    // check dimensions
    if(s < 2) {
	std::cerr << funame << "wrong array size: " << s << "\n";
	throw Error::Init();
    }

    _size = s;

    // check order
    const double* const end = x + _size;
    for(const double* p = x; p != end; ++p) {
	if(p != x && *p <= _xmax) {
	    std::cerr << funame << "x values are not in increasing order at i = " 
		      << p - x << "\n";
	    throw Error::Init();
	}
	_xmax = *p;
    }

    _xmin = *x;
    _ymin = *y;
    _ymax = *(y + _size - 1);

    // initialization
    _inbv = 1;
    _dim = _size + 2;
    _kn.resize(_dim + 4);
    _bc.resize(_dim);
    Array<double> work(5 * _dim);

    int_t n, k;
    dbint4_(x, y, _size, 2, 2, 0., 0., 1, 
	     _kn, _bc, n, k, work);

    if (n != _dim || k != 4) {
	std::cerr << funame << "slatec dbint4 failed\n";
	throw Error::Init();
    }
}

double Slatec::Spline::operator()(double x, int_t drv) const 
{
    const char funame [] = "Slatec::Spline::fit: ";

    if(!_size) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    static double work[12];
    if(x < _xmin || x > _xmax) {
	std::cerr << funame << " x is out of range: xmin = " 
		  << _xmin   << ", x = " << x << ", xmax = " 
		  << _xmax << "\n";
	throw Error::Range();
    }

    return dbvalu_(_kn, _bc, _dim, 4, drv, x, _inbv, work);
}
