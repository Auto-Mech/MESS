

#ifndef SLATEC_HH
#define SLATEC_HH

#include "error.hh"
#include "array.hh"

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

    AdamSolver(int, dde_f, double = -1., double = -1.) throw(Error::General);
    void run(double& x, double* y, double xout, void* param = 0, Mode mode =CONTINUE) throw(Error::General); 
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
    void init (const double*, const double*, int_t) throw(Error::General);
    Spline (const double* x, const double* y, int_t n) throw(Error::General) : _size(0) { init(x, y, n); }

    int_t size () const { return _size; }

    double arg_min () const { return _xmin; }
    double arg_max () const { return _xmax; }
    double fun_min () const { return _ymin; }
    double fun_max () const { return _ymax; }
    
    ~Spline () {}

    // evaluate i-th derivative at x
    double operator() (double, int_t =0) const throw(Error::General); 
  };

} // Slatec

#endif
