/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rlaswp.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_INCX  -2
#define MAX_INCX  2
#define MIN_LDA   5
#define MAX_LDA   10
#define MIN_N     5
#define MAX_N     10
#define MAX_ITER  100

void Rlaswp_test()
{
    int errorflag = FALSE;
    int i, iter;
    REAL_REF diff;

    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = MIN_N; n <= MAX_N; n++) {
	    for (int lda = MIN_LDA; lda <= MAX_LDA; lda++) {
		for (int k1 = lda - 1; k1 <= lda; k1++) {
		    for (int k2 = k1 + 1; k2 <= lda; k2++) {
#if defined VERBOSE_TEST
			printf("# n:%d, lda: %d, incx:%d k1:%d k2: %d\n", n, lda, incx, k1, k2);
#endif
			REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
			INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(k2, incx)];

			REAL *A = new REAL[matlen(lda, n)];
			INTEGER *ipiv = new INTEGER[veclen(k2, incx)];

			for (iter = 0; iter < MAX_ITER; iter++) {
			    set_random_vector(A_ref, A, matlen(lda, n));
			    set_random_vector(ipiv_ref, ipiv, veclen(k2, incx), lda);
#if defined ___MPACK_BUILD_WITH_MPFR___
			    dlaswp_f77(&n, A_ref, &lda, &k1, &k2, ipiv_ref, &incx);
#else
			    Rlaswp(n, A_ref, lda, k1, k2, ipiv_ref, incx);
#endif
			    Rlaswp(n, A, lda, k1, k2, ipiv, incx);
			    diff = infnorm(A_ref, A, matlen(lda, n), 1);
			    if (diff > EPSILON) {
				for (i = 0; i < veclen(k2, incx); i++) {
				    printf("%d %d\n", (int) ipiv_ref[i], (int)ipiv[i]);
				}
		                printf("error: "); printnum(diff); printf("\n");
				errorflag = TRUE;
				exit(1);
			    }
			}
			delete[]ipiv_ref;
			delete[]ipiv;
			delete[]A_ref;
			delete[]A;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rlaswp test failed...\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    Rlaswp_test();
    printf("Rlaswp test passed...\n");
    return (0);
}
