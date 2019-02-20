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
#define MAX_M 10
#define MIN_N -2
#define MAX_N 10
#define MIN_KL -2
#define MAX_KL 10
#define MIN_KU -2
#define MAX_KU 10
#define MIN_LDA -2
#define MAX_LDA 10
#define MIN_INCX -3
#define MAX_INCX 3
#define MIN_INCY -3
#define MAX_INCY 3
#define MAX_ITER 1

REAL_REF maxdiff = 0.0;

void Rgbmv_test3(const char *trans, REAL_REF alpha_ref, REAL_REF beta_ref, REAL alpha, REAL beta)
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
				REAL_REF *x_ref;
				REAL_REF *y_ref;
				REAL_REF *A_ref;
				REAL_REF diff;
				REAL *x;
				REAL *y;
				REAL *A;

				if (Mlsame(trans, "N")) {
				    x_ref = new REAL_REF[veclen(n, incx)];
				    y_ref = new REAL_REF[veclen(m, incy)];
				    A_ref = new REAL_REF[matlen(lda, n)];
				    x = new REAL[veclen(n, incx)];
				    y = new REAL[veclen(m, incy)];
				    A = new REAL[matlen(lda, n)];
				} else {
				    x_ref = new REAL_REF[veclen(m, incx)];
				    y_ref = new REAL_REF[veclen(n, incy)];
				    A_ref = new REAL_REF[matlen(lda, n)];
				    x = new REAL[veclen(m, incx)];
				    y = new REAL[veclen(n, incy)];
				    A = new REAL[matlen(lda, n)];
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
				    dgbmv_f77(trans, &m, &n, &kl, &ku, &alpha_ref, A_ref, &lda, x_ref, &incx, &beta_ref, y_ref, &incy);
				    mpack_errno1 = blas_errno;
#else
				    Rgbmv(trans, m, n, kl, ku, alpha_ref, A_ref, lda, x_ref, incx, beta_ref, y_ref, incy);
				    mpack_errno1 = mpack_errno;
#endif
				    Rgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
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
	printf("Rgbmv test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Rgbmv_test2(const char *trans)
{
    REAL_REF alpha_ref;
    REAL_REF beta_ref;
    REAL alpha;
    REAL beta;

//a=0, b=*;
    alpha_ref = 0.0;
    alpha = 0.0;
    set_random_number(beta_ref, beta);
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0;
    beta = 0.0;
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0;
    beta = 1.0;
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0;  beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0;
    alpha = 1.0;
    set_random_number(beta_ref, beta);
    Rgbmv_test3(trans, alpha_ref, beta_ref, alpha, beta);
}

void Rgbmv_test()
{
    Rgbmv_test2("N");
    Rgbmv_test2("T");
    Rgbmv_test2("C");
}

int main(int argc, char *argv[])
{
    Rgbmv_test();
    printf("Rgbmv test passed...\n");
    return (0);
}
