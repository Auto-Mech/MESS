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

#ifndef MATH_HH
#define MATH_HH

#include "error.hh"

#include <iostream>
#include <map>

namespace Math {

  class Exception {};

  // Newton-Raphson search method
  class NewtonRaphsonSearch {
  public:
    double tol;

    virtual double operator() (double, int) const =0;
    void find (double& guess) const;
  };

  // Bisection search method
  class BisectionSearch {
  public:
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
  public:
    double tol; // composite tolerance: (x2 - x1) * y < tol

    virtual double operator() (double x, int n) const  =0; // function and its first derivative
    double find (double x) const ;
    
    GradientSearch (double t = -1.);
  };
  inline GradientSearch::GradientSearch (double t) : tol(t) {} 

  // find minimum of the parabola passing throw three points
  double parabola_minimum (double x1, double y1, double x2, double y2, double x3, double y3, double* yp =0)
    ;

  class MinimumSearch {
    double _xtol;
    double _ytol;
    
  public:
    class Exception {
    public:
      double arg;
      Exception (double x) : arg(x) {}
    };

    virtual double _fun (double x)  = 0;
    double find (double xguess, double step, double* xp =0) ;
    
    MinimumSearch (double xt, double yt) : _xtol(xt), _ytol(yt) {}
  };

}

#endif
