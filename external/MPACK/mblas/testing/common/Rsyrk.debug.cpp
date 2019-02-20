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

#define MIN_K -2
#define MAX_K 8
#define MIN_N -2
#define MAX_N 8
#define MIN_LDA -2
#define MAX_LDA 8
#define MIN_LDC -2
#define MAX_LDC 8
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rsyrk_test3(const char *uplo, const char *trans, REAL_REF alpha_ref, REAL_REF beta_ref, REAL alpha, REAL beta)
{
    int errorflag = FALSE;
    int mpack_errno1, mpack_errno2;
    for (int n = MIN_N; n < MAX_N; n++) {
	for (int k = MIN_K; k < MAX_K; k++) {
	    int minlda;
	    if (Mlsame(trans, "N"))
		minlda = max(1, n);
	    else
		minlda = max(1, k);
	    for (int lda = minlda; lda < MAX_LDA; lda++) {
		for (int ldc = max(1, n); ldc < MAX_LDC; ldc++) {
#if defined VERBOSE_TEST
		    printf("#n is %d, k is %d, lda is %d, ldc is %d ", n, k, lda, ldc);
		    printf("uplo is %s, trans is %s \n", uplo, trans);
#endif
		    REAL_REF *A_ref;
		    REAL_REF *C_ref;
		    REAL *A;
		    REAL *C;

		    int ka;
		    if (Mlsame(trans, "N"))
			ka = k;
		    else
			ka = n;

		    A_ref = new REAL_REF[matlen(lda, ka)];
		    C_ref = new REAL_REF[matlen(ldc, n)];
		    A = new REAL[matlen(lda, ka)];
		    C = new REAL[matlen(ldc, n)];

		    for (int iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(A_ref, A, matlen(lda, ka));
			set_random_vector(C_ref, C, matlen(ldc, n));

			mpack_errno = 0; blas_errno = 0;
#if defined ___MPACK_BUILD_WITH_MPFR___
			dsyrk_f77(uplo, trans, &n, &k, &alpha_ref, A_ref, &lda, &beta_ref, C_ref, &ldc);
			mpack_errno1 = blas_errno;
#else
			Rsyrk(uplo, trans, n, k, alpha_ref, A_ref, lda, beta_ref, C_ref, ldc);
			mpack_errno1 = mpack_errno;
#endif
			Rsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
			mpack_errno2 = mpack_errno;
#if defined VERBOSE_TEST
			printf("errno: mpack %d, ref %d\n", mpack_errno1, mpack_errno2);
#endif
			if (mpack_errno1 != mpack_errno2) {
			    printf("error in Mxerbla!!\n");
			    exit(1);
			}

			REAL_REF diff = infnorm(C_ref, C, matlen(ldc, n), 1);
			if (diff > EPSILON) {
			    printf("error: "); printnum(diff);  printf("\n");
			    errorflag = TRUE;
			    exit(1);
			}
			if (maxdiff < diff)
			    maxdiff = diff;
			printf("max error: "); printnum(maxdiff); printf("\n");
		    }
		    delete[]C_ref;
		    delete[]A_ref;
		    delete[]C;
		    delete[]A;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rsyrk test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Rsyrk_test2(const char *side, const char *uplo)
{
    REAL_REF alpha_ref, beta_ref;
    REAL alpha, beta;

//a=*, b=*;
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0;
    beta = 1.0;
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0;
    alpha = 1.0;
    set_random_number(beta_ref, beta);
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0;
    beta = 0.0;
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=*;
    alpha_ref = 0.0;
    alpha = 0.0;
    set_random_number(beta_ref, beta);
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0; beta = 0.0;
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0; beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Rsyrk_test3(side, uplo, alpha_ref, beta_ref, alpha, beta);
}

void Rsyrk_test()
{
    Rsyrk_test2("U", "N");
    Rsyrk_test2("L", "N");
    Rsyrk_test2("U", "C");
    Rsyrk_test2("L", "C");
}

int main(int argc, char *argv[])
{
    Rsyrk_test();
    printf("Rsyrk test passed...\n");
    return (0);
}
