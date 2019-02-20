/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
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
#include <mblas.h>
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_M -2
#define MAX_M 4
#define MIN_N -2
#define MAX_N 4
#define MIN_KL -2
#define MAX_KL 4
#define MIN_KU -2
#define MAX_KU 4
#define MIN_LDA -2
#define MAX_LDA 4
#define MIN_INCX -2
#define MAX_INCX 2
#define MIN_INCY -2
#define MAX_INCY 2
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Cgbmv_test3(const char *trans, COMPLEX_REF alpha_ref, COMPLEX_REF beta_ref, COMPLEX alpha, COMPLEX beta)
{
    int errorflag = FALSE;
    int mpack_errno1, mpack_errno2;
    for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	    for (int lda = MIN_LDA; lda < MAX_LDA; lda++) {
		for (int ku = MIN_KU; ku < MAX_KU; ku++) {
		    for (int kl = MIN_KL; kl < MAX_KL; kl++) {
			for (int n = MIN_N; n < MAX_N; n++) {
			    for (int m = MIN_M; m < MAX_M; m++) {
#if defined VERBOSE_TEST
				printf
				    ("#m is %d, n is %d, kl is %d, ku is %d, lda is %d, incx is %d, incy is %d, trans is %s.\n",
				     m, n, kl, ku, lda, incx, incy, trans);
#endif
				COMPLEX *x;
				COMPLEX *y;
				COMPLEX *A;
				COMPLEX_REF *x_ref;
				COMPLEX_REF *y_ref;
				COMPLEX_REF *A_ref;
				REAL_REF diff;

				if (Mlsame(trans, "N")) {
				    x_ref = new COMPLEX_REF[veclen(n, incx)];
				    y_ref = new COMPLEX_REF[veclen(m, incy)];
				    A_ref = new COMPLEX_REF[matlen(lda, n)];
				    x = new COMPLEX[veclen(n, incx)];
				    y = new COMPLEX[veclen(m, incy)];
				    A = new COMPLEX[matlen(lda, n)];
				} else {
				    x_ref = new COMPLEX_REF[veclen(m, incx)];
				    y_ref = new COMPLEX_REF[veclen(n, incy)];
				    A_ref = new COMPLEX_REF[matlen(lda, n)];
				    x = new COMPLEX[veclen(m, incx)];
				    y = new COMPLEX[veclen(n, incy)];
				    A = new COMPLEX[matlen(lda, n)];
				}

				for (int iter = 0; iter < MAX_ITER; iter++) {
				    if (Mlsame(trans, "N")) {
					set_random_vector(x_ref, x, veclen(n, incx));
					set_random_vector(y_ref, y, veclen(m, incy));
					set_random_vector(A_ref, A, matlen(lda, n));
				    } else {
					set_random_vector(x_ref, x, veclen(m, incx));
					set_random_vector(y_ref, y, veclen(n, incy));
					set_random_vector(A_ref, A, matlen(lda, n));
				    }
				    mpack_errno = 0; blas_errno = 0;
#if defined ___MPACK_BUILD_WITH_MPFR___
				    zgbmv_f77(trans, &m, &n, &kl, &ku, &alpha_ref, A_ref, &lda, x_ref, &incx, &beta_ref, y_ref,
					      &incy);
				    mpack_errno1 = blas_errno;
#else
				    Cgbmv(trans, m, n, kl, ku, alpha_ref, A_ref, lda, x_ref, incx, beta_ref, y_ref, incy);
				    mpack_errno1 = mpack_errno;
#endif
				    Cgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
				    mpack_errno2 = mpack_errno;

#if defined VERBOSE_TEST
				    printf("errno: mpack %d, ref  %d\n", mpack_errno1, mpack_errno2);
#endif
				    if (mpack_errno1 != mpack_errno2) {
					printf("error in Mxerbla!!\n");
					exit(1);
				    }
				    if (Mlsame(trans, "N")) {
					diff = infnorm(y_ref, y, veclen(m, incy), 1);
				    } else {
					diff = infnorm(y_ref, y, veclen(n, incy), 1);
				    }
				    if (diff > EPSILON) {
					printf("error: "); printnum(diff); printf("\n");
					errorflag = TRUE;
					exit(-1);
				    }
				    if (maxdiff < diff)
					maxdiff = diff;
				    printf("max error: "); printnum(maxdiff); printf("\n");
				}
				delete[]x;
				delete[]y;
				delete[]A;
				delete[]x_ref;
				delete[]y_ref;
				delete[]A_ref;
			    }
			}
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Cgbmv test failed...\n");
	exit(1);
    } else {
	printf("Cgbmv test passed...\n");
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Cgbmv_test2(const char *trans)
{
    COMPLEX_REF alpha_ref;
    COMPLEX_REF beta_ref;
    COMPLEX alpha;
    COMPLEX beta;

//a=0, b=*;
    alpha_ref = 0.0;
    alpha = 0.0;
    set_random_number(beta_ref, beta);
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0;
    beta = 0.0;
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0;
    beta = 1.0;
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0;  beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0;
    alpha = 1.0;
    set_random_number(beta_ref, beta);
    Cgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);
}

void Cgbmv_test()
{
    Cgbmv_test2("N");
    Cgbmv_test2("C");
    Cgbmv_test2("T");
}

int main(int argc, char *argv[])
{
    Cgbmv_test();
    printf("Cgbmv test passed...\n");
    return (0);
}
