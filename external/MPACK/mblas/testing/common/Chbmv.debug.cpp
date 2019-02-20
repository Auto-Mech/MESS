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

#define MIN_N 1
#define MAX_N 10
#define MAX_K 10
#define MIN_LDA 1
#define MAX_LDA 10
#define MIN_INCX -2
#define MAX_INCX 3
#define MIN_INCY -2
#define MAX_INCY 2
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Chbmv_test3(const char *uplo, COMPLEX_REF alpha_ref, COMPLEX_REF beta_ref, COMPLEX alpha, COMPLEX beta)
{
    int errorflag = FALSE;
    int mpack_errno1, mpack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	    for (int n = MIN_N; n < MAX_N; n++) {
		for (int k = 0; k < MAX_K; k++) {
		    for (int lda = k + 1; lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
			printf("#k is %d, n is %d, lda is %d incx %d incy %d uplo %s\n", k, n, lda, incx, incy, uplo);
#endif
			COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
			COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
			COMPLEX_REF *y_ref = new COMPLEX_REF[veclen(n, incy)];
			COMPLEX *A = new COMPLEX[matlen(lda, n)];
			COMPLEX *x = new COMPLEX[veclen(n, incx)];
			COMPLEX *y = new COMPLEX[veclen(n, incy)];

			for (int iter = 0; iter < MAX_ITER; iter++) {
			    set_random_vector(A_ref, A, matlen(lda, n));
			    set_random_vector(x_ref, x, veclen(n, incx));
			    set_random_vector(y_ref, y, veclen(n, incy));

			    mpack_errno = 0;  blas_errno = 0;
#if defined ___MPACK_BUILD_WITH_MPFR___
			    zhbmv_f77(uplo, &n, &k, &alpha_ref, A_ref, &lda, x_ref, &incx, &beta_ref, y_ref, &incy);
			    mpack_errno1 = blas_errno;
#else
			    Chbmv(uplo, n, k, alpha_ref, A_ref, lda, x_ref, incx, beta_ref, y_ref, incy);
			    mpack_errno1 = mpack_errno;
#endif
			    Chbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
			    mpack_errno2 = mpack_errno;
#if defined VERBOSE_TEST
			    printf("errno: mpack %d, ref %d\n", mpack_errno1, mpack_errno2);
#endif
			    if (mpack_errno1 != mpack_errno2) {
				printf("error in Mxerbla!!\n");
				exit(1);
			    }
			    REAL_REF diff = infnorm(y_ref, y, veclen(n, incy), 1);
			    if (diff > EPSILON) {
				printf("error: "); printnum(diff); printf("\n");
				errorflag = TRUE;
			    }
			    if (maxdiff < diff)
				maxdiff = diff;
			    printf("max error: "); printnum(maxdiff); printf("\n");
			}
			delete[]A_ref;
			delete[]x_ref;
			delete[]y_ref;
			delete[]x;
			delete[]y;
			delete[]A;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Chbmv test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Chbmv_test2(const char *uplo)
{
    COMPLEX_REF alpha_ref, beta_ref;
    COMPLEX alpha, beta;

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0; beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0; beta = 0.0;
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0; beta = 1.0;
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=*;
    alpha_ref = 0.0; alpha = 0.0;
    set_random_number(beta_ref, beta);
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0; alpha = 1.0;
    set_random_number(beta_ref, beta);
    Chbmv_test3(uplo, alpha_ref, beta_ref, alpha, beta);
}

void Chbmv_test()
{
    Chbmv_test2("U");
    Chbmv_test2("L");
}

int main(int argc, char *argv[])
{
    Chbmv_test();
    printf("Chbmv test passed...\n");
    return (0);
}
