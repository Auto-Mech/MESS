/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rpotrf.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N      0
#define MAX_N     30		//should not be so large
#define MAX_LDA   30		//should not be so large
#define MAX_ITER  30

REAL_REF maxdiff = 0.0;

void Rpotrf_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
	    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
	    REAL *A = new REAL[matlen(lda, n)];
#if defined VERBOSE_TEST
	    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
#endif
	    j = 0;
	    while (j < MAX_ITER) {
//general (not necessary psd) case
		set_random_vector(A_ref, A, matlen(lda, n));
#if defined ___MPACK_BUILD_WITH_MPFR___
		dpotf2_f77(uplo, &n, A_ref, &lda, &info_ref);
#else
		Rpotrf(uplo, n, A_ref, lda, &info_ref);
#endif
		Rpotrf(uplo, n, A, lda, &info);

		if (info < 0) {
		    printf("info %d error\n", -(int) info);
		    errorflag = TRUE;
                    exit(1);
		}
		if (info > 0) {
		    printf("non psd matrix in %d-th (not an error)\n", (int) info);
		}
		if (info_ref != info) {
		    printf("info error! %d, %d\n", (int)info_ref, (int)info);
		    errorflag = TRUE;
                    exit(1);
		}
                diff = infnorm(A_ref, A, matlen(lda, n), 1);
		if (diff > EPSILON8) {
		    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
		    printf("error1: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
                    exit(1);
		}
	        if (maxdiff < diff)
		    maxdiff = diff;
	        printf("max error: "); printnum(maxdiff); printf("\n");
//psd case
		set_random_psdmat(A_ref, A, lda, n);
#if defined ___MPACK_BUILD_WITH_MPFR___
		dpotrf_f77(uplo, &n, A_ref, &lda, &info_ref);
#else
		Rpotrf(uplo, n, A_ref, lda, &info_ref);
#endif
		Rpotrf(uplo, n, A, lda, &info);

		if (info < 0) {
		    printf("info %d error\n", -(int) info);
		    errorflag = TRUE;
                    exit(1);
		}
		if (info > 0) {
		    printf("non psd matrix in %d-th (not an error)\n", (int) info);
		}
		if (info_ref != info) {
		    printf("info error! %d, %d\n", (int) info_ref, (int)info);
		    errorflag = TRUE;
                    exit(1);
		}
                diff = infnorm(A_ref, A, matlen(lda, n), 1);
		if (diff > EPSILON8) {
		    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
		    printf("error2: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
                    exit(1);
		}
	        if (maxdiff < diff)
		    maxdiff = diff;
	        printf("max error: "); printnum(maxdiff); printf("\n");

		j++;
	    }
	    delete[]A;
	    delete[]A_ref;
	}
	if (errorflag == TRUE) {
	    printf("Rpotrf test failed...\n");
	    exit(1);
	}
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Rpotrf_test(void)
{
    Rpotrf_test2("U");
    Rpotrf_test2("L");
}

int main(int argc, char *argv[])
{
    Rpotrf_test();
    printf("Rpotrf test passed...\n");
    return (0);
}
