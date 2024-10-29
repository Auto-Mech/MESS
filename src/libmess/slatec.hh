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

#ifndef SLATEC_HH
#define SLATEC_HH

#include "error.hh"
#include "array.hh"

#include <map>

namespace Slatec {

  typedef int int_t;// this is apparently the right type for fortran int

  extern "C" {
    typedef double (*arg1_f) (const double&);
    typedef void   (*dde_f) (const double& x, const double* y, 
			     double* dy, void* rp, void* ip);
  }

  /********************************************************************************
   *       Differential equations solver by the Adams-Bashforth-Moulton           *
   *       Predictor-Corrector formulas of orders one through twelve              *
   *******************************************************************************/

  class AdamSolver
  {
    dde_f _deriv;
    const int_t _size;
    int_t _info [15];
    Array<double> _rwork;
    int_t _lrw;

    static const int_t _liw = 51;
    int_t _iwork [_liw];

  public:
    enum Mode {RESTART, CONTINUE};

    Array<double> rel_tol;
    Array<double> abs_tol;

    AdamSolver(int, dde_f, double = -1., double = -1.) ;
    void run(double& x, double* y, double xout, void* param = 0, Mode mode =CONTINUE) ; 
  };

  /********************************************************************************
   *                              Spline fitting                                  *
   *******************************************************************************/

  class Spline
  {
    int_t _size;

    int_t         _dim;
    Array<double> _kn;   // spline knots
    Array<double> _bc;   // spline coefficients
    mutable int_t _inbv; // slatec internal parameter

    double _xmax;
    double _xmin;
    double _ymin;
    double _ymax;

    //Spline (const Spline&);
    //Spline& operator= (const Spline&);

  public:

    Spline () : _size(0) {}
    
    void init (const double*, const double*, int_t);

    void init (const std::map<double, double>&);
    
    Spline (const double* x, const double* y, int_t n)  : _size(0) { init(x, y, n); }

    Spline (const std::map<double, double>& m)  : _size(0) { init(m); }

    int_t size () const { return _size; }

    bool isinit () const { return _size; }

    double arg_min () const { return _xmin; }
    double arg_max () const { return _xmax; }
    double fun_min () const { return _ymin; }
    double fun_max () const { return _ymax; }
    
    ~Spline () {}

    // evaluate i-th derivative at x
    double operator() (double, int_t =0) const ; 
  };

} // Slatec

#endif
