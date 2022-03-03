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

#ifndef MATH_HH
#define MATH_HH

#include "error.hh"
#include "array.hh"

#include <iostream>
#include <map>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>

namespace Math {

  class Exception {};

  // Newton-Raphson search method
  //
  class NewtonRaphsonSearch {
    //
  public:
    //
    double tol;

    virtual double operator() (double, int) const =0;
    
    void find (double& guess) const;
  };

  // Bisection search method
  //
  class BisectionSearch {
    //
  public:
    //
    double tol;

    virtual double operator() (double) const =0;
    
    void find (double& x1, double& x2) const;
  };

  // zero seach algorithm for finding function minimum
  //
  class ZeroSearch {
  public:
    double tol; // composite tolerance: (x2 - x1) * y < tol


    virtual double operator() (double) const  =0;
    double find (double x1, double x2) const ;
    
    ZeroSearch (double t = -1.);
  };
  inline ZeroSearch::ZeroSearch (double t) : tol(t) {} 

  class GradientSearch {
    //
  public:
    //
    // composite tolerance: (x2 - x1) * y < tol
    //
    double tol;

    // function and its first derivative
    //
    virtual double operator() (double x, int n) const  =0;
    
    double find (double x) const ;
    
    GradientSearch (double t = -1.);
  };
  
  inline GradientSearch::GradientSearch (double t) : tol(t) {} 

  // find minimum of the parabola passing throw three points
  //
  double parabola_minimum (double x1, double y1, double x2, double y2, double x3, double y3, double* yp =0);

  class MinimumSearch {
    //
    double _xtol;
    
    double _ytol;
    
  public:
    //
    class Exception {
      //
    public:
      //
      double arg;
      
      Exception (double x) : arg(x) {}
    };

    virtual double _fun (double x)  = 0;
    
    double find (double xguess, double step, double* xp =0) ;
    
    MinimumSearch (double xt, double yt) : _xtol(xt), _ytol(yt) {}
  };

  // GNU scientific library interpolation wrapper
  //
  class Spline {
    //
    class Accel {
      //
      gsl_interp_accel* _value;

      Accel (const Accel&);

      Accel& operator= (const Accel&);
      
    public:
      //
      Accel () : _value(gsl_interp_accel_alloc()) { }

      ~Accel () { gsl_interp_accel_free(_value); }

      gsl_interp_accel* value() const { return _value; }
    };

    Accel              _accel;
    
    gsl_spline*        _spline;

    int                _type;

    double             _arg_min;

    double             _arg_max;

    double             _arg_range;
    
    // no copies
    //
    Spline (const Spline&);

    Spline& operator= (const Spline&);
    
  public:
    //
    enum {AKIMA = 1, PERIODIC = 2, STEFFEN = 4};

    Spline () : _spline(0) {}
    
    void init (const std::map<double, double>& data, int type =0);

    explicit Spline (const std::map<double, double>& data, int type =0);

    ~Spline () { if(_spline) gsl_spline_free(_spline); }

    bool isinit () const { return _spline; };

    double operator() (double, int =0) const;
  };

  inline Spline::Spline (const std::map<double, double>& data, int type) : _spline(0) { init(data, type); }

  // GSL mixed radix FFT
  //
  class DirectFFT {
    //
    int _size;

    gsl_fft_real_wavetable* _table;

    gsl_fft_real_workspace*  _work;
    
    // no copies
    //
    DirectFFT (const DirectFFT&);

    DirectFFT& operator= (const DirectFFT&);
    
  public:
    //
    DirectFFT (int);

    ~DirectFFT ();

    void operator() (Array<double>& data) const;

    static void normalize (Array<double>& data);

    int size () const { return _size; }
  };
}

#endif
