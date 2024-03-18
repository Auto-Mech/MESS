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

#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <gmpxx.h>

#ifdef WITH_MPACK

# define INT mpackint

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

#include <mpack/mlapack___float128.h>

// with mplapack
//
#else

# define INT mplapackint

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

#include <mplapack/mplapack__Float128.h>
#include <mplapack/mplapack__Float64x.h>

#endif

namespace Mpack {

  template<typename A>
  //
  void Rsyev (const char *jobz, const char *uplo, INT n, A * a, INT lda,
	      //
	      A * w, A * work, INT lwork, INT *info)
  {
#ifdef WITH_MPACK
    
    ::Rsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
    
#else
    
    ::Rsyev(jobz, uplo, n, a, lda, w, work, lwork, *info);

#endif
 
  }

  template<typename A>
  //
  Lapack::Vector eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec) 
  {
    const char funame [] = "Mpack::dd_eigenvalues: ";

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

    INT lwork, info;

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

    work.resize(std::max((INT)1, lwork));

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

}// Mpack namespace

#endif
