/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rorgqr.debug.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N 10
#define MAX_N 30
#define MIN_M 10
#define MAX_M 30
#define MIN_K 10
#define MAX_K 30
#define MIN_LDA 10
#define MAX_LDA 30
#define MAX_ITER 30

REAL_REF maxdiff = 0.0;

void Rorgqr_test()
{
    int errorflag = FALSE;
    int iter;
    int m, n, k;
    int lda;
    REAL_REF diff;
    INTEGER_REF info_ref, worksize_ref, lwork;
    INTEGER info, worksize;

    for (m = MIN_M; m <= MAX_M; m++) {
	for (n = MIN_N; n <= m; n++) {
	    for (k = MIN_K; k <= n; k++) {
		for (lda = max(1, m); lda <= MAX_LDA; lda++) {
#if defined VERBOSE_TEST
		    printf("# m %d n %d k %d lda %d\n", m, n, k, lda);
#endif
		    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
		    REAL_REF *tau_ref = new REAL_REF[veclen(k, 1)];
		    REAL_REF *work_ref = new REAL_REF[veclen(n, 1) * 1024];

		    REAL *A = new REAL[matlen(lda, n)];
		    REAL *tau = new REAL[veclen(k, 1)];
		    REAL *work = new REAL[veclen(n, 1) * 1024];

		    for (iter = 0; iter < MAX_ITER; iter++) {
			set_random_vector(A_ref, A, matlen(lda, n));
			set_random_vector(tau_ref, tau, veclen(k, 1));
			set_random_vector(work_ref, work, veclen(n, 1) * 1024);
// these workspace query might not be the same value.
			lwork = -1;
#if defined ___MPACK_BUILD_WITH_MPFR___
			dorgqr_f77(&m, &n, &k, A_ref, &lda, tau_ref, work_ref, &lwork, &info_ref);
#else
			Rorgqr(m, n, k, A_ref, lda, tau_ref, work_ref, lwork, &info_ref);
#endif
			Rorgqr(m, n, k, A, lda, tau, work, lwork, &info);
			worksize_ref = (INTEGER_REF) work_ref[0];
			worksize = (INTEGER) cast2double(work[0]);
#if defined VERBOSE_TEST
			printf("optimized worksize by dorgqr %d : by Rorgqr %d.\n", (int)worksize_ref, (int)worksize);
#endif
#ifdef DUMMY
//comparison of workspace is nonsense...
			if (worksize != worksize_ref)
			    printf("error in worksize\n");
#endif
			lwork = worksize;
#if defined ___MPACK_BUILD_WITH_MPFR___
			dorgqr_f77(&m, &n, &k, A_ref, &lda, tau_ref, work_ref, &lwork, &info_ref);
#else
			Rorgqr(m, n, k, A_ref, lda, tau_ref, work_ref, lwork, &info_ref);
#endif
			Rorgqr(m, n, k, A, lda, tau, work, lwork, &info);

			diff = infnorm(A_ref, A, matlen(lda, n), 1);
		        if (diff > EPSILON) {
		            printf("error in A: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
                            exit(1);
		        }
			diff = infnorm(tau_ref, tau, veclen(k, 1), 1);
		        if (diff > EPSILON) {
		            printf("error in t: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
                            exit(1);
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;
#ifdef DUMMY
//comparison of workspace is nonsense...
			diff = infnorm(work_ref, work, veclen(n, 1), 1);
		        if (diff > EPSILON) {
		            printf("error in t: "); printnum(diff); printf("\n");
		            errorflag = TRUE;
                            exit(1);
		        }
	                if (maxdiff < diff)
		            maxdiff = diff;
#endif
	                printf("max error: "); printnum(maxdiff); printf("\n");
		    }
		    delete[]tau_ref;
		    delete[]work_ref;
		    delete[]A_ref;
		    delete[]tau;
		    delete[]work;
		    delete[]A;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rorgqr test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Rorgqr_test();
    printf("Rorgqr test passed...\n");
    return (0);
}
