/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rsytd2.debug.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MAX_N     50
#define MAX_LDA   50
#define MAX_ITER  1

void Rsytd2_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    REAL_REF diff;
    INTEGER info;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
	    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
	    REAL_REF *d_ref = new REAL_REF[veclen(n, 1)];
	    REAL_REF *e_ref = new REAL_REF[veclen(n - 1, 1)];
	    REAL_REF *tau_ref = new REAL_REF[veclen(n - 1, 1)];

	    REAL *A = new REAL[matlen(lda, n)];
	    REAL *d = new REAL[veclen(n, 1)];
	    REAL *e = new REAL[veclen(n - 1, 1)];
	    REAL *tau = new REAL[veclen(n - 1, 1)];
#if defined VERBOSE_TEST
	    printf("#uplo %s, n:%d lda %d\n", uplo, n, lda);
#endif
	    j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(A_ref, A, matlen(lda, n));
		set_random_vector(d_ref, d, veclen(n, 1));
		set_random_vector(e_ref, e, veclen(n - 1, 1));
		set_random_vector(tau_ref, tau, veclen(n - 1, 1));
#if defined ___MPACK_BUILD_WITH_MPFR___
		dsytd2_f77(uplo, &n, A_ref, &lda, d_ref, e_ref, tau_ref, &info_ref);
#else
		Rsytd2(uplo, n, A_ref, lda, d_ref, e_ref, tau_ref, &info_ref);
#endif
		Rsytd2(uplo, n, A, lda, d, e, tau, &info);

		if (info < 0) {
		    printf("info %d error\n", -(int) info);
		    errorflag = TRUE;
		}
		if (info_ref != info) {
		    printf("info differ! %d, %d\n", (int) info_ref, (int)info);
		    errorflag = TRUE;
		}
		diff = infnorm(A_ref, A, matlen(lda, n), 1);
		if (diff > EPSILON) {
		    printf("Ret mat A error !!"); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		diff = infnorm(d_ref, d, veclen(n, 1), 1);
		if (diff > EPSILON) {
		    printf("diagonal part of tridiagonal mat error!!"); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		diff = infnorm(e_ref, e, veclen(n - 1, 1), 1);
		if (diff > EPSILON) {
		    printf("off-diagonal part of tridiagonal mat error!!"); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		j++;
	    }
	    delete[]tau_ref;
	    delete[]e_ref;
	    delete[]d_ref;
	    delete[]A_ref;
	    delete[]tau;
	    delete[]e;
	    delete[]d;
	    delete[]A;
	}
	if (errorflag == TRUE) {
	    printf("Rsytd2 test failed...\n");
	    exit(1);
	}
    }
}

void Rsytd2_test(void)
{
    Rsytd2_test2("U");
    Rsytd2_test2("L");
}

int main(int argc, char *argv[])
{
    Rsytd2_test();
    printf("Rsytd2 test passed...\n");
    return (0);
}
