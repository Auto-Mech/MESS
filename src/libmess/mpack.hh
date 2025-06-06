/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2016, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef MPACK_HH
#define MPACK_HH

#include "io.hh"
#include "lapack.hh"
#include "lapack.h"
#include "linpack.hh"

#if defined(WITH_MPACK) || defined(WITH_MPLAPACK)

#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <gmpxx.h>

#endif

#ifdef WITH_MPACK

# define MPACK_INT mpackint

#include <mpack/mpack_config.h>

#include <mpack/dd_complex.h>
#include <mpack/mlapack_dd.h>

#include <mpack/qd_complex.h>
#include <mpack/mlapack_qd.h>

#include <mpack/mpreal.h>
#include <mpack/mpcomplex.h>
#include <mpack/mlapack_mpfr.h>

#include <mpack/mpc_class.h>
#include <mpack/mlapack_gmp.h>

#ifdef WITH_FLOAT128
#include <mpack/mlapack___float128.h>
#endif

// with mplapack
//
#elif defined WITH_MPLAPACK

# define MPACK_INT mplapackint

#include <mplapack/mplapack_config.h>

#include <mplapack/dd_complex.h>
#include <mplapack/mplapack_dd.h>

#include <mplapack/qd_complex.h>
#include <mplapack/mplapack_qd.h>

#include <mplapack/mpreal.h>
#include <mplapack/mpcomplex.h>
#include <mplapack/mplapack_mpfr.h>

#include <mplapack/mpc_class.h>
#include <mplapack/mplapack_gmp.h>

#ifdef WITH_FLOAT128
#include <mplapack/mplapack__Float128.h>
#include <mplapack/mplapack__Float64x.h>
#endif

// with double
//
#else

# define MPACK_INT Lapack::int_t

#endif

namespace Mpack {
  //
  // float type
  //
  enum {DOUBLE, DD, QD, MPFR, GMP, FLOAT128, FLOAT64X};

  extern int mp_type;

  template<typename A>
  //
  void Rsyev (const char *jobz, const char *uplo, MPACK_INT n, A * a, MPACK_INT lda,
	      //
	      A * w, A * work, MPACK_INT lwork, MPACK_INT *info)
  {
#if defined WITH_MPACK

    ::Rsyev(jobz, uplo, n, a, lda, w, work, lwork, info);

#elif defined WITH_MPLAPACK

    ::Rsyev(jobz, uplo, n, a, lda, w, work, lwork, *info);

#endif
  }

  template<typename A>
  Lapack::Vector eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec = 0); 

  template<typename A>
  Lapack::Vector eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec) 
  {
    const char funame [] = "Mpack::eigenvalues: ";

    if(!mat.isinit()) {
      //
      ErrOut err_out;

      err_out << funame << "not initialized";
    }

    const char* jobz = "N";
    //
    if(evec)
      //
      jobz = "V";

    MPACK_INT lwork, info;

    Array<A> a(mat.size() * mat.size());

    Array<A> eval(mat.size());

    for(int i = 0; i < mat.size(); ++i)
      //
      for(int j = i; j < mat.size(); ++j)
	//
	a[i + j * mat.size()] = mat(i, j);

    //work space query
    //
    lwork = -1;

    Array<A> work(1);

    Rsyev<A>(jobz, "U", mat.size(), (A*)a, mat.size(), (A*)eval, (A*)work, lwork, &info);

    // matrix diagonalization
    //
    //lwork = (int)cast2double(work[0]);

    lwork = (long)work[0];

    work.resize(std::max((MPACK_INT)1, lwork));

    Rsyev<A>(jobz, "U", mat.size(), (A*)a, mat.size(), (A*)eval, (A*)work, lwork, &info);

    if(info < 0) {
      //
      ErrOut err_out;

      err_out << funame << "Rsyev: " << -info << "-th argument has an illegal value";
    }
    else if(info > 0) {
      //
      ErrOut err_out;

      err_out << funame << "Rsyev: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero";
    }

    Lapack::Vector res(mat.size());
  
    for(int i = 0; i < mat.size(); ++i)
      //
      res[i] = (double)eval[i];

    if(evec) {
      //
      for(int i = 0; i < mat.size(); ++i) {
	//
	for(int j = 0; j < mat.size(); ++j) {
	  //
	  (*evec)(i, j) = (double)a[i + j * mat.size()];
	}
      }
    }

    return res;
  }

  template<typename T>
  //
  void eigenvalues(const Matrix<T>& m, Vector<T>& eval, Matrix<T>* =0);
  
  template<typename T>
  //
  void eigenvalues(const Matrix<T>& m, Vector<T>& eval, Matrix<T>* p) {
    //
    const char funame [] = "Mpack::Matrix::eigenvalues: ";

    if(!m.isinit()) {
      //
      std::cerr << funame << "not initialized\n";

      throw Error::Init();
    }
    
    const char* jobz = "N";

    Matrix<T> a;

    Matrix<T>* evec;
    
    if(p) {
      //
      jobz = "V";

      evec = p;
    }
    else
      //
      evec = &a;

    evec->resize(m.size());

    eval.resize(m.size());
    
    (*evec) = m;
    
    MPACK_INT lwork, info;

    //work space query
    //
    lwork = -1;

    Array<T> work(1);

    Rsyev<T>(jobz, "U", m.size(), (T*)(*evec), m.size(), (T*)eval, (T*)work, lwork, &info);

    lwork = (long)work[0];

    work.resize(std::max((MPACK_INT)1, lwork));

    Rsyev<T>(jobz, "U", m.size(), (T*)(*evec), m.size(), (T*)eval, (T*)work, lwork, &info);

    if(info < 0) {
      //
      IO::log << IO::log_offset << funame << "Rsyev: " << -info << "-th argument has an illegal value\n";

      throw Error::Logic();
    }
    else if(info > 0) {
      //
      IO::log << IO::log_offset << funame << "Rsyev: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero";

      throw Error::Run();
    }
  }
  
  Lapack::Vector eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec = 0);

  inline Lapack::Vector eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec)
  {
#if defined(WITH_MPACK) || defined(WITH_MPLAPACK)
      
      if(mp_type == DD) {
        //
        return eigenvalues<dd_real>(mat, evec);
      }
      else if(mp_type == QD) {
        //
        return eigenvalues<qd_real>(mat, evec);
      }
      else if(mp_type == MPFR) {
        //
        return eigenvalues<mpreal>(mat, evec);
      }
      else if(mp_type == GMP) {
        //
        return eigenvalues<mpf_class>(mat, evec);
      }

#ifdef WITH_FLOAT128
#ifdef WITH_MPACK
        
      else if(mp_type == FLOAT128) {
        //
        return eigenvalues<__float128>(mat, evec);
      }
#else
      else if(mp_type == FLOAT128) {
        //
        return eigenvalues<_Float128>(mat, evec);
      }
      else if(mp_type == FLOAT64X) {
        //
        return eigenvalues<_Float64x>(mat, evec);
      }
#endif
#endif
      else {
        //
        return  mat.eigenvalues(evec);
      }
#else
      return  mat.eigenvalues(evec);
#endif
  }
  //
}// Mpack namespace

#endif
