/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpack_debug.h,v 1.55 2010/08/18 10:28:48 nakatamaho Exp $
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#ifndef _MPACK_DEBUG_H_
#define _MPACK_DEBUG_H_

#ifndef ___MPACK_DEBUG_BUILD___
#define ___MPACK_DEBUG_BUILD___
#endif

extern int blas_errno;		//for BLAS xerbla dispatch

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::complex;

#include <mpreal.h>
#include <mpcomplex.h>
using namespace mpfr;

#define EPSILON    1e-13
#define EPSILON6    1e-6
#define EPSILON7    1e-7
#define EPSILON8    1e-8
#define EPSILON9    1e-9
#define EPSILON10  1e-10
#define EPSILON11  1e-11
#define EPSILON12  1e-12
#define EPSILON13  1e-13
#define EPSILON14  1e-14
#define TRUE     1
#define FALSE    0
//#define GMP_P_FORMAT  "%+36.32Fe"
//#define MPFR_P_FORMAT "%+36.32Re"
#define GMP_P_FORMAT  "%+68.64Fe"
#define MPFR_P_FORMAT "%+68.64Re"
//double
//#define GMP_P_FORMAT  "%+21.16Fe"
//#define MPFR_P_FORMAT "%+21.16Re"
#define P_FORMAT      "%+21.16e"
#define QD_PRECISION  64
#define DD_PRECISION  32

#if defined __MPACK_BUILD_DEBUG_CPP__
#define _MPACK_DEBUG_EXTERN_
#else
#define _MPACK_DEBUG_EXTERN_ extern
#endif

#if defined ___MPACK_BUILD_WITH_GMP___ || ___MPACK_BUILD_WITH_QD___ || ___MPACK_BUILD_WITH_DD___ || ___MPACK_BUILD_WITH_DOUBLE___
#include <mblas_mpfr.h>
#include <mlapack_mpfr.h>
#endif

_MPACK_DEBUG_EXTERN_ gmp_randstate_t uniformrandomstate_mpfr;
#if defined ___MPACK_BUILD_WITH_GMP___
_MPACK_DEBUG_EXTERN_ gmp_randclass *uniformrandomstate_gmp;
#endif

#if defined ___MPACK_BUILD_WITH_MPFR___
typedef int INTEGER_REF;
typedef double REAL_REF;
typedef complex < double >COMPLEX_REF;
#else
typedef mpackint INTEGER_REF;
typedef mpreal REAL_REF;
typedef mpcomplex COMPLEX_REF;
#endif

mpreal mpf_randomnumber(mpreal);
mpcomplex mpc_randomnumber(mpcomplex);
double mpf_randomnumber(double);
complex < double >mpc_randomnumber(complex < double >);

void printnum(mpreal rtmp);
void printnum(mpcomplex ctmp);
void printnum(complex < double >ctmp);
void printnum(double rtmp);
void set_random_number(double &a, mpreal & b);
void set_random_number(complex < double >&a, mpcomplex & b);
void set_random_number(INTEGER_REF & a, INTEGER & b);

void set_random_number1to2(double &a, mpreal & b);
void set_random_number1to2(complex < double >&a, mpcomplex & b);
void set_random_number1to2(INTEGER_REF & a, INTEGER & b);

REAL_REF infnorm(COMPLEX_REF * vec_ref, COMPLEX * vec, int len, int inc);
REAL_REF infnorm(REAL_REF * vec_ref, REAL * vec, int len, int inc);
INTEGER_REF infnorm(INTEGER_REF * vec_ref, INTEGER * vec, int len, int inc);
REAL_REF reldiff(REAL_REF * vec_ref, REAL * vec, int len, int inc);
REAL_REF infnorm_mat(REAL_REF * vec_ref, REAL * vec, int n, int m, int lda);

inline void set_random_vector(INTEGER_REF * vec_ref, INTEGER * vec, int len, int n)
{
    if (len <= 0)
	return;
    for (int i = 0; i < len; i++) {
	vec_ref[i] = int (drand48() * (n - 1)) + 1;
	vec[i] = vec_ref[i];
    }
}

inline int veclen(int n, int inc)
{
    return std::max(1, (1 + (abs(n) - 1) * abs(inc)));
}

inline int vecplen(int n)
{
    if (n <= 0)
	n = 1;
    return n * (n + 1) / 2;
}

inline int matlen(int lda, int n)
{
    return std::max(1, abs(lda) * abs(n));
}

#if defined ___MPACK_BUILD_WITH_GMP___
void printnum(mpf_class rtmp);
void printnum(mpc_class ctmp);
mpf_class mpf_randomnumber(mpf_class);
mpc_class mpc_randomnumber(mpc_class);
void set_random_number(mpreal & a, mpf_class & b);
void set_random_number(mpcomplex & a, mpc_class & b);
void set_random_number1to2(mpreal & a, mpf_class & b);
void set_random_number1to2(mpcomplex & a, mpc_class & b);
#endif
#if defined ___MPACK_BUILD_WITH_QD___
void printnum(qd_real rtmp);
void printnum(qd_complex rtmp);
qd_real mpf_randomnumber(qd_real);
qd_complex mpc_randomnumber(qd_complex);
void set_random_number(mpreal & a, qd_real & b);
void set_random_number(mpcomplex & a, qd_complex & b);
void set_random_number1to2(mpreal & a, qd_real & b);
void set_random_number1to2(mpcomplex & a, qd_complex & b);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
void printnum(dd_real rtmp);
void printnum(dd_complex rtmp);
dd_real mpf_randomnumber(dd_real);
dd_complex mpc_randomnumber(dd_complex);
void set_random_number(mpreal & a, dd_real & b);
void set_random_number(mpcomplex & a, dd_complex & b);
void set_random_number1to2(mpreal & a, dd_real & b);
void set_random_number1to2(mpcomplex & a, dd_complex & b);
#endif
#if defined ___MPACK_BUILD_WITH_DOUBLE___
void set_random_number(mpreal & a, double &b);
void set_random_number(mpcomplex & a, std::complex < double >&b);
void set_random_number1to2(mpreal & a, double &b);
void set_random_number1to2(mpcomplex & a, std::complex < double >&b);
#endif

template < class X > void printmat(int N, int M, X * A, int LDA)
{
    X tmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    tmp = A[i + j * LDA];
	    printnum(tmp);
	    if (j < M - 1)
		printf(", ");
	}
	if (i < N - 1)
	    printf("]; ");
	else
	    printf("] ");
    }
    printf("]");
}

template < class X > void printvec(X * A, int len)
{
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
	tmp = A[i];
	printnum(tmp);
	if (i < len - 1)
	    printf(", ");
    }
    printf("]");
}

template < class X_REF, class X > void set_random_vector(X_REF * vec_ref, X * vec, int len)
{
    if (len <= 0)
	return;
    for (int i = 0; i < len; i++) {
	set_random_number(vec_ref[i], vec[i]);
    }
}

template < class X_REF, class X > void set_random_psdmat(X_REF * p_ref, X * p, int ldp, int n)
{
    X_REF *tmpmat_ref = new X_REF[matlen(ldp, n)];
    X_REF mtmp_ref;
    X *tmpmat = new X[matlen(ldp, n)];
    X mtmp;

    for (int j = 0; j < n; j++) {
	for (int i = 0; i < n; i++) {
	    set_random_number1to2(tmpmat_ref[i + j * ldp], tmpmat[i + j * ldp]);
	}
    }
    for (int j = 0; j < n; j++) {
	for (int i = 0; i < n; i++) {
	    mtmp_ref = 0.0;
	    mtmp = 0.0;
	    for (int k = 0; k < n; k++) {
		mtmp_ref = mtmp_ref + tmpmat_ref[i + k * ldp] * tmpmat_ref[j + k * ldp];
		mtmp = mtmp + tmpmat[i + k * ldp] * tmpmat[j + k * ldp];
	    }
	    p_ref[i + j * ldp] = mtmp_ref;
	    p[i + j * ldp] = mtmp;
	}
    }
    delete[]tmpmat;
    delete[]tmpmat_ref;
}

template < class X_REF, class X > void set_random_symmmat_cond(X_REF * p_ref, X * p, int ldp, int n, int cond)
{
    //all calculations should be done with mpfr.
    mpreal *tmpmat1_mpreal = new mpreal[matlen(ldp, n)];
    mpreal *tmpmat2_mpreal = new mpreal[matlen(ldp, n)];
    mpreal *tmpmat3_mpreal = new mpreal[matlen(ldp, n)];
    mpreal rtmp;

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    tmpmat1_mpreal[i + j * ldp] = 0.0;
	}
    }
    for (int i = 0; i < n; i++) {
	rtmp = mpreal(cond) * mpreal(2.0 * (i + 1) - n) / mpreal(2.0 * n);
	tmpmat1_mpreal[i + i * ldp] = pow((mpreal) 10.0, rtmp);
    }
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    tmpmat2_mpreal[i + j * ldp] = (mpf_randomnumber(rtmp) + (mpreal) 1.0) / (mpreal) 2.0;
	}
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    rtmp = 0.0;
	    for (int k = 0; k < n; k++) {
		for (int l = 0; l < n; l++) {
		    rtmp = rtmp + tmpmat2_mpreal[i + k * ldp] * tmpmat1_mpreal[k + l * ldp] * tmpmat2_mpreal[j + l * ldp];
		}
	    }
	    tmpmat3_mpreal[i + j * ldp] = rtmp;
	}
    }

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
#if defined (___MPACK_BUILD_WITH_MPFR___)
	    p[i + j * ldp] = (tmpmat3_mpreal[i + j * ldp]);
#elif defined (___MPACK_BUILD_WITH_GMP___)
	    p[i + j * ldp] = cast2mpf_class(tmpmat3_mpreal[i + j * ldp]);
#elif defined (___MPACK_BUILD_WITH_QD___)
	    p[i + j * ldp] = cast2qd_real(tmpmat3_mpreal[i + j * ldp]);
#elif defined (___MPACK_BUILD_WITH_DD___)
	    p[i + j * ldp] = cast2dd_real(tmpmat3_mpreal[i + j * ldp]);
#elif defined (___MPACK_BUILD_WITH_DOUBLE___)
	    p[i + j * ldp] = (double) (tmpmat3_mpreal[i + j * ldp]);
#endif
	    p_ref[i + j * ldp] = tmpmat3_mpreal[i + j * ldp];
	}
    }

    delete[]tmpmat1_mpreal;
    delete[]tmpmat2_mpreal;
    delete[]tmpmat3_mpreal;
}

template < class X_REF, class X > void set_hilbertmat(X_REF * p_ref, X * p, int ldp, int n)
{
    X_REF mtmp_ref, One_ref;
    X mtmp, One;
    One_ref = 1.0; One = 1.0;

    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtmp_ref = (i + 1) + (j + 1) - 1;
	    mtmp = (i + 1) + (j + 1) - 1;
	    p_ref[i + j * ldp] = One_ref / mtmp_ref;
	    p[i + j * ldp] = One / mtmp;
	}
    }
}

#endif
