/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Ctrtrs.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
#include <mlapack.h>
#include <mpack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N 1
#define MAX_N 10		//shold not be so large
#define MIN_LDA 1
#define MAX_LDA 10		//shold not be so large
#define MIN_LDB 1
#define MAX_LDB 10		//shold not be so large
#define MAX_NRHS 10		//shold not be so large
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Ctrtrs_test2(const char *uplo, const char *trans, const char *diag)
{
    int errorflag = FALSE;
    int iter;
    INTEGER_REF n, lda, ldb, nrhs;
    REAL_REF diff;
    INTEGER_REF info_ref;
    INTEGER info;

    for (n = MIN_N; n <= MAX_N; n++) {
	for (lda = max(1, (int)n); lda <= MAX_LDA; lda++) {
	    for (ldb = max(1, (int)n); ldb <= MAX_LDB; ldb++) {
		for (nrhs = 1; nrhs <= MAX_NRHS; nrhs++) {
#if defined VERBOSE_TEST
		    printf("# uplo %s, trans %s, diag %s, n %d, lda %d, ldb %d, nrhs %d\n", uplo, trans, diag, (int) n, (int) lda, (int) ldb, (int) nrhs);
#endif
		    COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
		    COMPLEX_REF *B_ref = new COMPLEX_REF[matlen(ldb, nrhs)];

		    COMPLEX *A = new COMPLEX[matlen(lda, n)];
		    COMPLEX *B = new COMPLEX[matlen(ldb, nrhs)];

		    for (iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(A_ref, A, matlen(lda, n));
			set_random_vector(B_ref, B, matlen(ldb, nrhs));
#if defined ___MPACK_BUILD_WITH_MPFR___
			ztrtrs_f77(uplo, trans, diag, &n, &nrhs, A_ref, &lda, B_ref, &ldb, &info_ref);
#else
			Ctrtrs(uplo, trans, diag, n, nrhs, A_ref, lda, B_ref, ldb, &info_ref);
#endif
			Ctrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, &info);

			if (info != info_ref) {
			    printf("info differ! %d, %d\n", (int) info, (int) info_ref);
			    errorflag = TRUE;
			}
			diff = infnorm(B_ref, B, matlen(ldb, nrhs), 1);

			if (diff > EPSILON7) {
			    printf("error: "); printnum(diff); printf("\n");
			    errorflag = TRUE;
			    exit(1);
			}
			if (maxdiff < diff)
			    maxdiff = diff;
			printf("max error: "); printnum(maxdiff); printf("\n");
		    }
		    delete[]A_ref;
		    delete[]B_ref;
		    delete[]A;
		    delete[]B;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Ctrtrs test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Ctrtrs_test()
{
    Ctrtrs_test2("U", "N", "N");
    Ctrtrs_test2("U", "N", "U");
    Ctrtrs_test2("U", "T", "N");
    Ctrtrs_test2("U", "T", "U");
    Ctrtrs_test2("U", "C", "N");
    Ctrtrs_test2("U", "C", "U");
    Ctrtrs_test2("L", "N", "N");
    Ctrtrs_test2("L", "N", "U");
    Ctrtrs_test2("L", "T", "N");
    Ctrtrs_test2("L", "T", "U");
    Ctrtrs_test2("L", "C", "N");
    Ctrtrs_test2("L", "C", "U");
}

int main(int argc, char *argv[])
{
    Ctrtrs_test();
    printf("Ctrtrs test passed...\n");
    return (0);
}
