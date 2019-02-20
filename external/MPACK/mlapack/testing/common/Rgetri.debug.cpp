/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgetri.debug.cpp,v 1.8 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_N      1
#define MAX_N     30
#define MAX_LDA   30
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Rgetri_test()
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref, lwork_ref;
    REAL_REF diff;
    INTEGER info, lwork;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
	    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
	    INTEGER_REF *ipiv_ref = new INTEGER_REF[veclen(n, 1)];

	    REAL *A = new REAL[matlen(lda, n)];
	    INTEGER *ipiv = new INTEGER[veclen(n, 1)];
#if defined VERBOSE_TEST
	    printf("#n:%d lda %d\n", n, lda);
#endif
//these workspace query might not be the same value.
	    lwork_ref = -1;
	    lwork = -1;
	    REAL_REF *work_ref = new REAL_REF[1];
	    REAL *work = new REAL[1];
#if defined ___MPACK_BUILD_WITH_MPFR___
	    dgetri_f77(&n, A_ref, &lda, ipiv_ref, work_ref, &lwork_ref, &info_ref);
#else
	    Rgetri(n, A_ref, lda, ipiv_ref, work_ref, lwork_ref, &info_ref);
#endif
	    Rgetri(n, A, lda, ipiv, work, lwork, &info);

	    lwork_ref = (int) cast2double (work_ref[0]);
	    lwork = (int) cast2double (work[0]);
#if defined VERBOSE_TEST
	    printf("optimized worksize by Rgetri %d : by dgetri %d.\n", (int) lwork_ref, (int) lwork);
#endif
#ifdef DUMMY
//comparison of workspace is nonsense...
	    if (worksize != worksize_ref)
		printf("error in worksize\n");
#endif
	    delete[]work;
	    delete[]work_ref;
	    work_ref = new REAL_REF[max(1, (int)lwork_ref)];
	    work = new REAL[max(1, (int) lwork)];
	    j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(A_ref, A, matlen(lda, n));
		set_random_vector(work_ref, work, veclen(lwork, 1));
#if defined ___MPACK_BUILD_WITH_MPFR___
		dgetrf_f77(&n, &n, A_ref, &lda, ipiv_ref, &info_ref);
		dgetri_f77(&n, A_ref, &lda, ipiv_ref, work_ref, &lwork_ref, &info_ref);
#else
		Rgetrf(n, n, A_ref, lda, ipiv_ref, &info_ref);
		Rgetri(n, A_ref, lda, ipiv_ref, work_ref, lwork_ref, &info_ref);
#endif
		Rgetrf(n, n, A, lda, ipiv, &info);
		Rgetri(n, A, lda, ipiv, work, lwork, &info);

		if (info < 0) {
		    printf("info %d error\n", -(int) info);
		}
		if (info_ref != info) {
		    printf("info differ! %d, %d\n", (int) info_ref, (int) info);
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
	    delete[]work;
	    delete[]work_ref;
	    delete[]ipiv_ref;
	    delete[]A_ref;
	    delete[]ipiv;
	    delete[]A;
	}
	if (errorflag == TRUE) {
	    printf("Rgetri test failed...\n");
	    exit(-1);
	}
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Rgetri_test();
    printf("Rgetri test passed...\n");
    return (0);
}
