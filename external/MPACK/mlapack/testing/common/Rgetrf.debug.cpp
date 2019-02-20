/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgetrf.debug.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N     0
#define MAX_N     30
#define MIN_M     0
#define MAX_M     30
#define MAX_LDA   30
#define MAX_ITER  3

REAL_REF maxdiff = 0.0;

void Rgetrf_test()
{
    int errorflag = FALSE;
    int j = 0, i;
    REAL_REF diff;
    INTEGER_REF info_ref;
    INTEGER info;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int m = MIN_M; m < MAX_M; m++) {
	    for (int lda = max(m, 1); lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
		printf("# n:%d m:%d lda:%d\n", n, m, lda);
#endif
		REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
		INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(min(m, n), 1)];

		REAL *A = new REAL[matlen(lda, n)];
		INTEGER *ipiv = new INTEGER[veclen(min(m, n), 1)];

		j = 0;
		while (j < MAX_ITER) {
		    set_random_vector(A_ref, A, matlen(lda, n));
		    ipiv_ref[0] = 0;
		    ipiv[0] = 0;
		    for (i = 0; i < min(m, n); i++) {
			ipiv_ref[i] = 0;
			ipiv[i] = 0;
		    }
#if defined ___MPACK_BUILD_WITH_MPFR___
		    dgetrf_f77(&m, &n, A_ref, &lda, ipiv_ref, &info_ref);
#else
		    Rgetrf(m, n, A_ref, lda, ipiv_ref, &info_ref);
#endif
		    Rgetrf(m, n, A, lda, ipiv, &info);

		    diff = infnorm(A_ref, A, matlen(lda, n), 1);
		    if (diff > EPSILON) {
			printf("error: "); printnum(diff); printf("\n");
			errorflag = TRUE;
			exit(1);
		    }
		    if (maxdiff < diff)
			maxdiff = diff;
		    printf("max error: "); printnum(maxdiff); printf("\n");

		    int idiff = infnorm(ipiv_ref, ipiv, veclen(min(m, n), 1), 1);
//          for(i=0;i<min(m,n);i++){ printf("%d %d\n",ipiv[i],ipivd[i]);} 
		    if (idiff > 0) {
			printf("error pivoting %d!!\n", (int) idiff);
			errorflag = TRUE;
		    }
		    j++;
		}
		delete[]ipiv_ref;
		delete[]ipiv;
		delete[]A_ref;
		delete[]A;
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rgetrf test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Rgetrf_test();
    printf("Rgetrf test passed...\n");
    return (0);
}
