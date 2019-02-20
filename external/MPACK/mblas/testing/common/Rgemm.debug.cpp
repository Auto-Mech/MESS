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
#define MIN_K -2
#define MAX_K 10
#define MIN_LDA -2
#define MAX_LDA 10
#define MIN_LDB -2
#define MAX_LDB 10
#define MIN_LDC -2
#define MAX_LDC 10
#define MAX_ITER 3

REAL_REF maxdiff = 0.0;

void Rgemm_test3(const char *transa, const char *transb, REAL_REF alpha_ref, REAL_REF beta_ref, REAL alpha, REAL beta)
{
    int errorflag = FALSE;
    int mpack_errno1, mpack_errno2;
    for (int k = MIN_K; k < MAX_K; k++) {
	for (int n = MIN_N; n < MAX_N; n++) {
	    for (int m = MIN_M; m < MAX_M; m++) {

		int minlda, minldb;

		if (Mlsame(transa, "N"))
		    minlda = max(1, m);
		else
		    minlda = max(1, k);

		if (Mlsame(transb, "N"))
		    minldb = max(1, k);
		else
		    minldb = max(1, n);

		for (int lda = minlda; lda < MAX_LDA; lda++) {
		    for (int ldb = minldb; ldb < MAX_LDB; ldb++) {
			for (int ldc = max(1, m); ldc < MAX_LDC; ldc++) {
#if defined VERBOSE_TEST
			    printf("#k is %d, n is %d, m is %d, lda is %d, ldb is %d, ldc is %d\n", k, n, m, lda, ldb, ldc);
			    printf("#transa is %s, transb is %s\n", transa, transb);
#endif
			    REAL_REF *A_ref;
			    REAL_REF *B_ref;
			    REAL_REF *C_ref;
			    REAL_REF diff;
			    REAL *A;
			    REAL *B;
			    REAL *C;

			    int ka, kb;

			    if (Mlsame(transa, "N"))
				ka = k;
			    else
				ka = m;

			    if (Mlsame(transb, "N"))
				kb = n;
			    else
				kb = k;

			    A_ref = new REAL_REF[matlen(lda, ka)];
			    B_ref = new REAL_REF[matlen(ldb, kb)];
			    C_ref = new REAL_REF[matlen(ldc, n)];
			    A = new REAL[matlen(lda, ka)];
			    B = new REAL[matlen(ldb, kb)];
			    C = new REAL[matlen(ldc, n)];

			    for (int iter = 0; iter < MAX_ITER; iter++) {
				set_random_vector(A_ref, A, matlen(lda, ka));
				set_random_vector(B_ref, B, matlen(ldb, kb));
				set_random_vector(C_ref, C, matlen(ldc, n));

				mpack_errno = 0; blas_errno = 0;
#if defined ___MPACK_BUILD_WITH_MPFR___
				dgemm_f77(transa, transb, &m, &n, &k, &alpha_ref, A_ref, &lda, B_ref, &ldb, &beta_ref, C_ref, &ldc);
				mpack_errno1 = blas_errno;
#else
				Rgemm(transa, transb, m, n, k, alpha_ref, A_ref, lda, B_ref, ldb, beta_ref, C_ref, ldc);
				mpack_errno1 = mpack_errno;
#endif
				Rgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
				mpack_errno2 = mpack_errno;

#if defined VERBOSE_TEST
				printf("errno: mpack %d, ref %d\n", mpack_errno1, mpack_errno2);
#endif
				if (mpack_errno1 != mpack_errno2) {
				    printf("error in Mxerbla!!\n");
				    exit(1);
				}
				diff = infnorm(C_ref, C, matlen(ldc, n), 1);
				if (diff > EPSILON) {
				    printf("error: "); printnum(diff); printf("\n");
				    errorflag = TRUE;
				}
				if (maxdiff < diff)
				    maxdiff = diff;
				printf("max error: "); printnum(maxdiff); printf("\n");
			    }
			    delete[]C_ref;
			    delete[]B_ref;
			    delete[]A_ref;
			    delete[]C;
			    delete[]B;
			    delete[]A;
			}
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rgemm test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Rgemm_test2(const char *uplo, const char *trans)
{
    REAL_REF alpha_ref;
    REAL_REF beta_ref;
    REAL alpha;
    REAL beta;

//a=0, b=*;
    alpha_ref = 0.0;
    alpha = 0.0;
    set_random_number(beta_ref, beta);
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=0;
    set_random_number(alpha_ref, alpha);
    beta_ref = 0.0;
    beta = 0.0;
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=*, b=1;
    set_random_number(alpha_ref, alpha);
    beta_ref = 1.0;
    beta = 1.0;
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//alpha=*, beta=*
    set_random_number(alpha_ref, alpha);
    set_random_number(beta_ref, beta);
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=0;
    alpha_ref = 0.0; beta_ref = 0.0;
    alpha = 0.0;  beta = 0.0;
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=0, b=1;
    alpha_ref = 0.0; beta_ref = 1.0;
    alpha = 0.0; beta = 1.0;
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=0;
    alpha_ref = 1.0; beta_ref = 0.0;
    alpha = 1.0; beta = 0.0;
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=1;
    alpha_ref = 1.0; beta_ref = 1.0;
    alpha = 1.0; beta = 1.0;
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);

//a=1, b=*;
    alpha_ref = 1.0;
    alpha = 1.0;
    set_random_number(beta_ref, beta);
    Rgemm_test3(uplo, trans, alpha_ref, beta_ref, alpha, beta);
}

void Rgemm_test()
{
    Rgemm_test2("N", "N");
    Rgemm_test2("N", "T");
    Rgemm_test2("N", "C");
    Rgemm_test2("T", "N");
    Rgemm_test2("T", "C");
    Rgemm_test2("T", "T");
    Rgemm_test2("C", "N");
    Rgemm_test2("C", "T");
    Rgemm_test2("C", "C");
}

int main(int argc, char *argv[])
{
    Rgemm_test();
    printf("Rgemm test passed...\n");
    return (0);
}
