/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rposv.debug.cpp,v 1.2 2010/08/07 05:50:10 nakatamaho Exp $
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
#define MIN_NRHS   1
#define MAX_N      8		//should not be so large
#define MAX_LDA    8		//should not be so large
#define MAX_LDB    8		//should not be so large
#define MAX_NRHS   8		//should not be so large
#define MAX_ITER   3

REAL_REF maxdiff = 0.0;

void Rposv_test2(const char *uplo)
{
    int errorflag = FALSE;
    int j = 0;
    INTEGER_REF info_ref;
    INTEGER info;
    REAL_REF diff;
    REAL_REF rtmp_ref;
    REAL rtmp;

    for (int n = MIN_N; n < MAX_N; n++) {
	for (int nrhs = MIN_NRHS; nrhs < n; nrhs++) {
	    for (int lda = max(n, 1); lda < MAX_LDA; lda++) {
		for (int ldb = max(n, 1); ldb < MAX_LDB; ldb++) {
		    REAL_REF *A_ref = new REAL_REF[matlen(lda, n)];
		    REAL_REF *B_ref = new REAL_REF[matlen(ldb, nrhs)];
		    REAL_REF *C_ref = new REAL_REF[matlen(ldb, nrhs)];
		    REAL *A = new REAL[matlen(lda, n)];
		    REAL *B = new REAL[matlen(ldb, nrhs)];
		    REAL *C = new REAL[matlen(ldb, nrhs)];
#if defined VERBOSE_TEST
		    printf("n:%d lda %d, ldb %d, nrhs %d, uplo %s\n", n, lda, ldb, nrhs, uplo);
#endif
		    j = 0;
		    while (j < MAX_ITER) {
			set_random_psdmat(A_ref, A, lda, n);
			set_random_vector(C_ref, C, matlen(ldb, nrhs));
			set_random_vector(B_ref, B, matlen(ldb, nrhs)); //clear B matrix
//made an answer first so that not distrubed by the numerical errors.
		        for (int p = 0; p < n; p++ ) {
			  for (int q = 0; q < nrhs; q++ ) {
			    rtmp_ref = 0.0; rtmp = 0.0;
			    for (int r = 0; r < n; r++ ) {
			      rtmp_ref += A_ref [ p + r * lda ] * C_ref [ r + q * ldb ];
			      rtmp += A [ p + r * lda ] * C [ r + q * ldb ];
			    }
			    B_ref [ p + q * ldb ] = rtmp_ref;
			    B [ p + q * ldb ] = rtmp;
			  }
			}
#if defined ___MPACK_BUILD_WITH_MPFR___
			dposv_f77(uplo, &n, &nrhs, A_ref, &lda, B_ref, &ldb, &info_ref);
#else
			Rposv(uplo, n, nrhs, A_ref, lda, B_ref, ldb, &info_ref);
#endif
			Rposv(uplo, n, nrhs, A, lda, B, ldb, &info);
			if (info < 0) {
			    printf("info %d error\n", -(int) info);
			    errorflag = TRUE;
			    exit(1);
			}
			if (info_ref > 0) {
			    printf("not a psd matrix in %d-th (not an error) %d\n", (int) info_ref, (int) info);
			    continue;
			}
			if (info_ref != info) {
			    printf("info error! %d, %d\n", (int) info_ref, (int) info); errorflag = TRUE; exit(1);
			}
			diff = infnorm_mat(A_ref, A, n, n, lda);
			if (diff > EPSILON7) {
			    printf("n:%d lda %d, uplo %s\n", n, lda, uplo);
			    printf("error2: "); printnum(diff); printf("\n");
			    printf("A_ref = "); printmat (n, n, A_ref, lda); printf("\n");
			    printf("A     = "); printmat (n, n, A, lda); printf("\n");
			    errorflag = TRUE;
			    exit(1);
			}
			diff = infnorm(B_ref, B, matlen(ldb, nrhs), 1);
			if (diff > EPSILON7) {
			    printf("n:%d ldb %d, uplo %s\n", n, ldb, uplo);
			    printf("error2: "); printnum(diff); printf("\n");
			    printf("B_ref = "); printmat (n, nrhs, B_ref, ldb); printf("\n");
			    printf("B     = "); printmat (n, nrhs, B, ldb); printf("\n");
			    errorflag = TRUE;
			    exit(1);
			}
			if (maxdiff < diff) maxdiff = diff;
			printf("max error: "); printnum(maxdiff); printf("\n");
			j++;
		    }
		    delete[]A;
		    delete[]B;
		    delete[]C;
		    delete[]A_ref;
		    delete[]B_ref;
		    delete[]C_ref;
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rposv test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

void Rposv_test(void)
{
    Rposv_test2("U");
    Rposv_test2("L");
}

int main(int argc, char *argv[])
{
    Rposv_test();
    printf("Rposv test passed...\n");
    return (0);
}
