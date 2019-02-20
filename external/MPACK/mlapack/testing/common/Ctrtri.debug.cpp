/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Ctrtri.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 0
#define MAX_N 10		//shold not be so large
#define MIN_LDA 0
#define MAX_LDA 10		//shold not be so large
#define MAX_ITER 10

REAL_REF maxdiff = 0.0;

void Ctrtri_test2(const char *uplo, const char *diag)
{
    int errorflag = FALSE;
    int iter;
    int n, lda;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;

    for (n = MIN_N; n <= MAX_N; n++) {
	for (lda = max(1, n); lda <= MAX_LDA; lda++) {
#if defined VERBOSE_TEST
	    printf("# uplo %s, diag %s, n %d, lda %d\n", uplo, diag, n, lda);
#endif
	    COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
	    COMPLEX *A = new COMPLEX[matlen(lda, n)];

	    for (iter = 0; iter < MAX_ITER; iter++) {
		set_random_vector(A_ref, A, matlen(lda, n));
#if defined ___MPACK_BUILD_WITH_MPFR___
		ztrtri_f77(uplo, diag, &n, A_ref, &lda, &info_ref);
#else
		Ctrtri(uplo, diag, n, A_ref, lda, &info_ref);
#endif
		Ctrtri(uplo, diag, n, A, lda, &info);

		if (info_ref != info) {
		    printf("info differ! %d, %d\n", (int) info_ref, (int) info);
		    errorflag = TRUE;
		}

		diff = infnorm(A_ref, A, matlen(lda, n), 1);
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
	    delete[]A;
	}
    }
    if (errorflag == TRUE) {
	printf("Ctrtri test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Ctrtri_test()
{
    Ctrtri_test2("U", "N");
    Ctrtri_test2("U", "U");
    Ctrtri_test2("L", "N");
    Ctrtri_test2("L", "U");
}

int main(int argc, char *argv[])
{
    Ctrtri_test();
    printf("Ctrtri test passed...\n");
    return (0);
}
