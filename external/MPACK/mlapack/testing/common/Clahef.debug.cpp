/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clahef.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 */
#include <mblas.h>
#include <mlapack.h>
#include <mpack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_N      2
#define MAX_N     10
#define MAX_NB    10
#define MAX_LDA    9
#define MAX_LDW    9
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Clahef_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    REAL_REF diff;
    INTEGER_REF info_ref, kb = 0;
    INTEGER info;

    for (int n = MIN_N; n <= MAX_N; n++) {
	for (int nb = 2; nb < MAX_NB; nb++) {
	    for (int k = 0; k < 4; k++) {
		if (k == 0) { kb = nb - 1; }
		if (k == 1) { kb = nb; }
		if (k == 2) {
		  if (n <= nb) { kb = n; } else { kb = nb - 1; }
		}
		if (k == 3) {
		    if (n <= nb) { kb = n; } else { kb = nb; }
		}

		for (int lda = max(n, 1); lda <= MAX_LDA; lda++) {
		    for (int ldw = max(n, 1); ldw <= MAX_LDW; ldw++) {

			COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
			COMPLEX_REF *W_ref = new COMPLEX_REF[matlen(ldw, nb)];
			INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(n, 1)];

			COMPLEX *A = new COMPLEX[matlen(lda, n)];
			COMPLEX *W = new COMPLEX[matlen(ldw, nb)];
			INTEGER *ipiv = new INTEGER[veclen(n, 1)];
#if defined VERBOSE_TEST
			printf("#uplo %s n:%d nb:%d kb:%d ldw:%d\n", uplo, (int)n, (int)nb, (int)kb, (int)ldw);
#endif
			j = 0;
			while (j < MAX_ITER) {
			    set_random_vector(A_ref, A, matlen(lda, n));
#if defined ___MPACK_BUILD_WITH_MPFR___
			    zlahef_f77(uplo, &n, &nb, &kb, A_ref, &lda, ipiv_ref, W_ref, &ldw, &info_ref);
#else
			    Clahef(uplo, n, nb, kb, A_ref, lda, ipiv_ref, W_ref, ldw, &info_ref);
#endif
			    Clahef(uplo, n, nb, kb, A, lda, ipiv, W, ldw, &info);

			    if (info < 0) {
				printf("info %d error\n", -(int) info);
			    }
			    if (info != info_ref) {
				printf("info differ! %d, %d\n", (int) info, (int) info_ref);
				errorflag = TRUE;
			    }
			    diff = infnorm(A_ref, A, matlen(lda, n), 1);
			    if (diff > EPSILON) {
				printf("error: "); printnum(diff); printf("\n");
				errorflag = TRUE;
				exit(1);
			    }
			    if (maxdiff < diff)
				maxdiff = diff;
			    printf("max error: "); printnum(maxdiff); printf("\n");
			    j++;
			}
			delete[]ipiv_ref;
			delete[]W_ref;
			delete[]A_ref;
			delete[]ipiv;
			delete[]W;
			delete[]A;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Clahef test failed...\n");
	exit(-1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Clahef_test()
{
    Clahef_test2("L");
    Clahef_test2("U");
}

int main(int argc, char *argv[])
{
    Clahef_test();
    printf("Clahef test passed...\n");
    return (0);
}
