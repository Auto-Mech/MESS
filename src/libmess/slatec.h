#ifndef SLATEC_H
#define SLATEC_H

#include "slatec.hh"

extern "C" 
{

  /********************************************************************************
   *                                 Grid integrator                              *
   *******************************************************************************/
  void davint_ (const double* x, const double* y, const Slatec::int_t& n, 
		const double& xlo, const double& xup, double& ans, 
		Slatec::int_t& ierr); 
  
  /*******************************************************************************
   *                                   Spline fit                                *
   ******************************************************************************/
  void  dbint4_ (const double* x, const double* y, const Slatec::int_t& ndata, 
		 const Slatec::int_t& ibcl, const Slatec::int_t& ibcr, 
		 const double& fbcl, const double& fbcr, const Slatec::int_t& kntopt, 
		 double* t, double* bcoef, Slatec::int_t& n, Slatec::int_t& k, double* w);

  double dbvalu_ (const double* t, const double* a, const Slatec::int_t& n, 
		  const Slatec::int_t& k, const Slatec::int_t& ideriv, const double& x, 
		  Slatec::int_t& inbv, double* work);

  /********************************************************************************
   *       Differential equations solver by the Adams-Bashforth-Moulton           *
   *       Predictor-Corrector formulas of orders one through twelve              *
   *******************************************************************************/
  void ddeabm_(Slatec::dde_f DF, const Slatec::int_t& NEQ, double& T, double* Y, 
	       const double& TOUT, const Slatec::int_t* INFO, double* RTOL, 
	       double* ATOL, Slatec::int_t& IDID, double* RWORK, const Slatec::int_t& LRW, 
	       Slatec::int_t* IWORK, const Slatec::int_t& LIW, void* RPAR, void* IPAR);

  /*********************************************************************************
   *                  Differential equations solver by Runge-Kutta method          *
   ********************************************************************************/
  void dderkf_(Slatec::dde_f df, const Slatec::int_t& neq, double& t, double* y, 
	       const double& tout, const Slatec::int_t* info, double& rtol, 
	       double& atol, Slatec::int_t& idid, double*  rwork, const Slatec::int_t& lrw,
	       Slatec::int_t* iwork, const Slatec::int_t& liw, void* rpar, void* ipar);

  /**********************************************************************************
   *                     Machine specific integer constants                         *
   *********************************************************************************/
  Slatec::int_t i1mach_(const Slatec::int_t&);

  /**********************************************************************************
   *           Machine specific double precision constants (see also dlamch)        *
   *********************************************************************************/
  double d1mach_(const Slatec::int_t&);

  /**********************************************************************************
   *         Symmetric system of linear equations solver (together with dsifa)      *
   *********************************************************************************/
  void dsisl_(double* a, const Slatec::int_t& lda, const Slatec::int_t& n, 
	      Slatec::int_t* kpvt, double* b);

  /***********************************************************************************
   *          Symmetric system of linear equations solver (with dsisl)               *
   **********************************************************************************/
  void dsifa_(double* a, const Slatec::int_t& lda, const Slatec::int_t& n, 
	      Slatec::int_t* kpvt, Slatec::int_t& info);

  /***********************************************************************************
   *              General system of linear equations solver (with dsisl)             *
   **********************************************************************************/
  void dgefs_(double* a,const Slatec::int_t& lda,const Slatec::int_t& n, double* v, 
	      const Slatec::int_t& itask, Slatec::int_t& ind, double* work,
	      Slatec::int_t* iwork);

  /***********************************************************************************
   *                                  Gaus integrator                                *
   **********************************************************************************/
  void dgaus8_(Slatec::arg1_f fun, const double& a, const double& b, 
	       double& err, double& ans, Slatec::int_t& ierr);

  /************************************************************************************
   *                                Elliptic integral                                 *
   ***********************************************************************************/
  double drf_(const double& X,const double& Y,const double& Z,Slatec::int_t& IERR);

  /************************************************************************************
   *                            Bessel of first kind, zero order                      *
   ***********************************************************************************/
  double dbesj0_(const double&);

  /************************************************************************************
   *                        Bessel of first kind, first order                         *
   ***********************************************************************************/
  double dbesj1_(const double&);
  
  /************************************************************************************
   *                  Modified bessel of first kind, zero order                       *
   ***********************************************************************************/
  double dbesi0_(const double&);

  /************************************************************************************
   *              Modified bessel of first kind, zero order (predexponent)            *
   ***********************************************************************************/
  double dbsi0e_(const double&);

  /************************************************************************************
   *                     Modified bessel of third kind, zero order                    *
   ***********************************************************************************/
  double dbesk0_(const double&);

  /************************************************************************************
   *          Modified bessel of third kind, zero order (predexponent)                *
   ***********************************************************************************/
  double dbsk0e_(const double&);

  /************************************************************************************
   *                       Search for a zero of a function                            *
   ***********************************************************************************/

  void dfzero_(Slatec::arg1_f fun, double& b, double& c, const double& guess, 
	       const double& rel_tol, const double& abs_tol, Slatec::int_t& iflag);


}

#endif
