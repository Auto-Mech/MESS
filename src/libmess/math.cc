/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2021, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#include "math.hh"
#include "array.hh"

#include <cmath>
#include <vector>
#include <set>
#include <map>

void Math::NewtonRaphsonSearch::find (double& x) const
{
  const char funame [] = "Math::NewtonRaphsonSearch::find: ";

  static const int imax = 100;// maximum number of iterations

  double dtemp;
  for(int i = 0; i < imax; ++i) {
    dtemp = (*this)(x, 1);
    if(dtemp == 0.) {
      std::cerr << funame << "zero derivative\n";
      throw Error::Range();
    }

    dtemp = (*this)(x, 0) / dtemp;
    x = x - dtemp;
    if(dtemp < tol && dtemp > -tol)
      return;
  }
  std::cerr << funame << "maximum number of iterations has been reached\n";
  throw Error::Run();
}

// bisection zero search
void Math::BisectionSearch::find (double& xmin, double& xmax) const
{
  const char funame [] = "Math::BisectionSearch::find: ";

  double dtemp, x, y;

  double ymin = (*this)(xmin);
  if(ymin == 0.) {
    xmax = xmin;
    return;
  }
  double ymax = (*this)(xmax);
  if(ymax == 0.) {
    xmin = xmax;
    return;
  }
  if(ymin > 0. && ymax > 0. || ymin < 0. && ymax < 0.) {
    std::cerr << funame << "function has the same sign on both sides of the search interval\n";
    throw Error::Range();
  }

  if(xmin > xmax) {
    std::swap(xmin, xmax);
    std::swap(ymin, ymax);
  }
  double sign = ymin < ymax ? 1. : -1.;
  while(xmax - xmin > tol) {
    x = (xmin + xmax) / 2.;
    y = (*this)(x);
    if(y == 0.) {
      xmin = x;
      xmax = x;
      return;
    }

    if(y * sign > 0.) {
      xmax = x;
      ymax = y;
    }
    else {
      xmin = x;
      ymin = y;
    }
  }
}

double Math::ZeroSearch::find (double x1, double x2) const 
{
  static const char funame [] = " Math::ZeroSearch::find: ";

  if(tol <= 0.) {
    std::cerr << funame << "tolerance was not initialized\n";
    throw Error::Init();
  }

  double dtemp;

  double v1 = (*this)(x1);
  double v2 = (*this)(x2);

  if(v1 == 0.)
    return x1;

  if(v2 == 0.)
    return x2;

  if(v1 < 0. && v2 < 0. || v1 > 0. && v2 > 0.) {
    std::cerr << funame << "function values are of the same sign on both ends, no zero\n";
    throw Math::Exception();
  }

  double x, v;

  while(1) {
    // gradient search
    x = (x1 * v2 - x2 * v1) / (v2 - v1);
    v = (*this)(x);

    dtemp = (x2 - x1) * v;
    if(dtemp > -tol && dtemp < tol) {
      return x;
    }
   

    if(v > 0. && v1 > 0. || v < 0. && v1 < 0.) {
      x1 = x;
      v1 = v;
    }
    else {
      x2 = x;
      v2 = v;
    }

    // bisection search
    x = (x2 + x1) / 2.;
    v = (*this)(x);

    dtemp = (x2 - x1) * v;
    if(dtemp > -tol && dtemp < tol) {
      return x;
    }
   

    if(v > 0. && v1 > 0. || v < 0. && v1 < 0.) {
      x1 = x;
      v1 = v;
    }
    else {
      x2 = x;
      v2 = v;
    }
  }
}

double Math::GradientSearch::find (double x) const 
{
  static const char funame [] = "Math::GradientSearch::find: ";

  static const int step_max = 100;

  double v, s, dtemp;

  for(int i = 0; i < step_max; ++i) {
    v = (*this)(x, 1);
    s = (*this)(x, 2);
    dtemp = v * v / s;
    if( dtemp < tol && -dtemp < tol) {
      return x;
    }
    x -= v / s;
  }

  std::cerr << funame << "did not converged\n";
  throw Math::Exception();
}

double Math::parabola_minimum (double x1, double y1, double x2, double y2, double x3, double y3, double* yp)
  
{
  static const char funame [] = "Math::find_parabola_minimum: ";

  if(x1 == x2 || x2 == x3) {
    std::cerr << funame << "arguments are the same\n";
    throw Error::Run();
  }

  double u1 = (x2 + x1) / 2.;
  double v1 = (y2 - y1) / (x2 - x1);

  double u2 = (x3 + x2) / 2.;
  double v2 = (y3 - y2) / (x3 - x2);

  if(v1 == v2) {
    std::cerr << funame << "points are on the same line\n";
    throw Exception();
  }

  double x = (v2 * u1 - v1 * u2) / (v2 - v1);

  double y = 
    y1 * (x - x2) * (x - x3) / (x1 - x2) / (x1 - x3) +
    y2 * (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3) +
    y3 * (x - x1) * (x - x2) / (x3 - x1) / (x3 - x2);

  if(yp)
    *yp = y;

  return x;

}

double Math::MinimumSearch::find (double x1, double step, double* xp) 
{
  static const char funame [] = "Math::MinimumSearch::find: ";
  static const double eps = 1. + 1.e-10;

  double y1, x0, y0, x2, y2, x3, y3;
  
  y1 = _fun(x1);
   
  x2 = x1 - step;
  y2 = _fun(x2);
    
  if(y2 < y1) {
    while (1) {

      x0 = x2 - step;
      y0 = _fun(x0);

      if(y0 >= y2)
	break;
	  
      x1 = x2;
      y1 = y2;

      x2 = x0;
      y2 = y0;
    }
    std::swap(x0, x2);
    std::swap(y0, y2);
  }
  else {
    while(1) {

      x0 = x1 + step;
      y0 = _fun(x0);

      if(y0 >= y1)
	break;

      x2 = x1;
      y2 = y1;

      x1 = x0;
      y1 = y0;
    }
    std::swap(x0, x1);
    std::swap(y0, y1);
  }


  while(1) {// x2 < x0 < x1, y1 >= y0 <= y2 situation

    if(x1 - x2 <= _xtol) {

#ifdef DEBUG
      std::cout << funame << "x-converged\n";
#endif

      if(xp)
	*xp = x0;
      return y0;
    }
      
    // find a minimum of the parabola going throw x1, y1, x, y, and x2, y2
    try{
      x3 = Math::parabola_minimum(x1, y1, x0, y0, x2, y2, &y3);
    }
    catch (Math::Exception) {
      if(xp)
	*xp = x0;
      return y0;
    }

    if(y3 > (y0 > 0. ? y0 / eps : y0 * eps)) {
      std::cerr << funame << "WARNING: numerically marginal\n";
      if(xp)
	*xp = x0;
      return y0;
    }

    if(x3 > (x1 > 0. ? x1 / eps : x1 * eps)) {
      std::cerr << funame << "WARNING: numerically marginal\n";
      if(xp)
	*xp = x1;
      return y1;

    }

    if(x3 < (x2 > 0. ? x2 * eps : x2 / eps)) {
      std::cerr << funame << "WARNING: numerically marginal\n";
      if(xp)
	*xp = x2;
      return y2;
    }
      
    if(y0 - y3 <= _ytol) {

#ifdef DEBUG
      std::cout << funame << "y-converged\n";
#endif

      if(xp)
	*xp = x3;
      return y3;
    }
      
    y3 = _fun(x3);

    if(y3 < y0) {// x3, y3 new minimum

      if(x3 < x0) {
	x1 = x0;
	y1 = y0;
      }
      else {
	x2 = x0;
	y2 = y0;
      }

      x0 = x3;
      y0 = y3;

    }
    else {// x3, y3 new border
      
      if(x3 < x0) {
	x2 = x3;
	y2 = y3;
      }
      else {
	x1 = x3;
	y1 = y3;
      }
    }
  }// x2 < x0 < x1, y1 >= y0 <= y2 situation
}

void Math::Spline::init (const std::map<double, double>& data, int t) throw(Error::General)
{
  const char funame [] = "Math::Spline::init:";
  
  int    itemp;
  double dtemp;

  _type = t;

  _arg_min = data.begin()->first;

  _arg_max = data.rbegin()->first;

  _arg_range = _arg_max - _arg_min;
  
  const gsl_interp_type* spline_type;

  if(_type & STEFFEN) {
    //
    spline_type = gsl_interp_steffen;
  }
  else if((_type & AKIMA) && (_type & PERIODIC)) {
    //
    spline_type = gsl_interp_akima_periodic;
  }
  else if((_type & AKIMA) && !(_type & PERIODIC)) {
    //
    spline_type = gsl_interp_akima;
  }
  else if(!(_type & AKIMA) && (_type & PERIODIC)) {
    //
    spline_type = gsl_interp_cspline_periodic;
  }
  else if(!(_type & AKIMA) && !(_type & PERIODIC)) {
    //
    spline_type = gsl_interp_cspline;
  }

  itemp = gsl_interp_type_min_size(spline_type);
  
  if(data.size() < itemp) {
    //
    std::cerr << funame << "number of data points, " << data.size() << ", out of range: " << itemp << "\n";

    throw Error::Range();
  }

  if(_type & PERIODIC && data.begin()->second != data.rbegin()->second) {
    //
    std::cerr << funame << "end point values mismatch: " << data.begin()->second << ", " << data.rbegin()->second << "\n";

    throw Error::Range();
  }

  _spline = gsl_spline_alloc(spline_type, data.size());

  Array<double> x((int)data.size());
  
  Array<double> y((int)data.size());

  itemp = 0;
  
  for(std::map<double, double>::const_iterator cit = data.begin(); cit != data.end(); ++cit, ++itemp) {
    //
    x[itemp] = cit->first;

    y[itemp] = cit->second;
  }

  gsl_spline_init(_spline, x, y, data.size());
}

double Math::Spline::operator() (double arg, int der) const throw(Error::General)
{
  const char funame [] = "Math::Spline::operator(): ";

  double dtemp;
  int    itemp;

  if(!_spline) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }
  
  if(arg > _arg_max || arg < _arg_min) {
    //
    if(!(_type & PERIODIC)) {
      //
      std::cerr << funame << "argument value, " << arg << ", out of range: [" << _arg_min << ", " << _arg_max << "]\n";

      throw Error::Range();
    }

    arg -= _arg_min;

    arg -= std::floor(arg / _arg_range) * _arg_range;

    arg += _arg_min;

    if(arg < _arg_min) {
      //
      std::cerr << funame << "WARNING: argument, " << arg << ", out of range, ["
		<< _arg_min << ", " << _arg_max << "]: rounding up\n";
      
      arg = _arg_min;
    }
    else if(arg > _arg_max) {
      //
      std::cerr << funame << "WARNING: argument, " << arg << ", out of range, ["
		<< _arg_min << ", " << _arg_max << "]: rounding down\n";
      
      arg = _arg_max;
    }
  }

  switch(der) {
    //
  case 0:
    //
    return gsl_spline_eval(_spline, arg, _accel.value());

  case 1:
    //
    return gsl_spline_eval_deriv(_spline, arg, _accel.value());

  case 2:
    //
    return gsl_spline_eval_deriv2(_spline, arg, _accel.value());

  default:
    //
    std::cerr << funame << "derivative index out of range: " << der << "\n";

    throw Error::Range();
  }
}

Math::DirectFFT::DirectFFT (int n) : _size(n)
{
  const char funame [] = "Math::DirectFFT::DirectFFT: ";

  if(n < 1) {
    //
    std::cerr << funame << "size out of range: " << n << "\n";

    throw Error::Range();
  }

  _table = gsl_fft_real_wavetable_alloc(n);

  _work  = gsl_fft_real_workspace_alloc(n);
}

Math::DirectFFT::~DirectFFT ()
{
  gsl_fft_real_workspace_free(_work);
  
  gsl_fft_real_wavetable_free(_table);
}
  
void Math::DirectFFT::operator() (Array<double>& data) const
{
  const char funame [] = "Math::DirectFFT::operator(): ";

  if(data.size() != size()) {
    //
    std::cerr << funame << "sizes mismatch: " << data.size() << ", " << size() << "\n";

    throw Error::Range();
  }
  
  gsl_fft_real_transform(data, 1, data.size(), _table, _work);
}

void Math::DirectFFT::normalize (Array<double>& data)
{
  double dtemp;
  
  if(!data.size())
    //
    return;
  
  data[0] /= data.size();

  dtemp = data.size() / 2.;
  
  for(int i = 1; i < data.size(); ++i)
    //
    if(i % 2) {
      //
      data[i] /= dtemp;
    }
    else
      //
      data[i] /= -dtemp;

  if(!(data.size() % 2))
    //
    data.back() /= 2.;
}
