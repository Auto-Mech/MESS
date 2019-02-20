/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Chetrd.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N     1
#define MAX_N     12
#define MAX_LDA   12
#define MAX_ITER  2

REAL_REF maxdiff = 0.0;

void Chetrd_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    INTEGER info, lwork;
    REAL_REF diff;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {

	    COMPLEX_REF *A_ref = new COMPLEX_REF[matlen(lda, n)];
	    REAL_REF *d_ref = new REAL_REF[veclen(n, 1)];
	    REAL_REF *e_ref = new REAL_REF[veclen(n - 1, 1)];
	    COMPLEX_REF *tau_ref = new COMPLEX_REF[veclen(n - 1, 1)];

	    COMPLEX *A = new COMPLEX[matlen(lda, n)];
	    REAL *d = new REAL[veclen(n, 1)];
	    REAL *e = new REAL[veclen(n - 1, 1)];
	    COMPLEX *tau = new COMPLEX[veclen(n - 1, 1)];

#if defined VERBOSE_TEST
	    printf("#uplo %s, n:%d lda %d\n", uplo, n, lda);
#endif
	    lwork_ref = -1;
	    lwork = -1;
	    COMPLEX_REF *work_ref = new COMPLEX_REF[1];
	    COMPLEX *work = new COMPLEX[1];
#if defined ___MPACK_BUILD_WITH_MPFR___
	    zhetrd_f77(uplo, &n, A_ref, &lda, d_ref, e_ref, tau_ref, work_ref, &lwork_ref, &info_ref);
#else
	    Chetrd(uplo, n, A_ref, lda, d_ref, e_ref, tau_ref, work_ref, lwork_ref, &info_ref);
#endif
	    Chetrd(uplo, n, A, lda, d, e, tau, work, lwork, &info);
	    lwork_ref = (int) cast2double(work_ref[0].real());
	    lwork = (int) cast2double(work[0].real());
#if defined VERBOSE_TEST
	    printf("optimized worksize by Chetrd %d : by zhetrd %d.\n", (int) lwork, (int)lwork_ref);
#endif
#ifdef DUMMY
//comparison of workspace is nonsense...
	    if (worksize != worksized)
		printf("error in worksize\n");
#endif
	    delete[]work;
	    delete[]work_ref;
	    work_ref = new COMPLEX_REF[max(1, (int)lwork_ref)];
	    work = new COMPLEX[max(1, (int) lwork)];
	    j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(A_ref, A, matlen(lda, n));
		set_random_vector(d_ref, d, veclen(n, 1));
		set_random_vector(e_ref, e, veclen(n - 1, 1));
		set_random_vector(tau_ref, tau, veclen(n - 1, 1));
		set_random_vector(work_ref, work, veclen(lwork, 1));

#if defined ___MPACK_BUILD_WITH_MPFR___
		zhetrd_f77(uplo, &n, A_ref, &lda, d_ref, e_ref, tau_ref, work_ref, &lwork_ref, &info_ref);
#else
		Chetrd(uplo, n, A_ref, lda, d_ref, e_ref, tau_ref, work_ref, lwork_ref, &info_ref);
#endif
		Chetrd(uplo, n, A, lda, d, e, tau, work, lwork, &info);

		if (info != info_ref) {
		    printf("info differ! %d, %d\n", (int) info, (int) info_ref);
		    errorflag = TRUE;
		}
		if (info < 0) {
		    continue;
		}
		diff = infnorm(A_ref, A, matlen(lda, n), 1);
		if (diff > EPSILON) {
		    printf("error in A: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		    exit(1);
		}
		if (maxdiff < diff)
		    maxdiff = diff;
		diff = infnorm(d_ref, d, veclen(n, 1), 1);
		if (diff > EPSILON) {
		    printf("error in d: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		    exit(1);
		}
		if (maxdiff < diff)
		    maxdiff = diff;

		diff = infnorm(e_ref, e, veclen(n - 1, 1), 1);
		if (diff > EPSILON) {
		    printf("error in e: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		    exit(1);
		}
		if (maxdiff < diff)
		    maxdiff = diff;
		j++;
	    }
	    delete[]work;
	    delete[]work_ref;
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
	    printf("Chetrd test failed...\n");
	    exit(1);
	}
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Chetrd_test(void)
{
    Chetrd_test2("U");
    Chetrd_test2("L");
}

int main(int argc, char *argv[])
{
    Chetrd_test();
    printf("Chetrd test passed...\n");
    return (0);
}
