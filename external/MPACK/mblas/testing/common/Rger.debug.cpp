/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
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
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MIN_INCX  -3
#define MAX_INCX  3
#define MIN_INCY  -3
#define MAX_INCY  3
#define MIN_N     1
#define MAX_N     12
#define MIN_M     1
#define MAX_M     12
#define MIN_LDA   1
#define MAX_LDA   12
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

void Rger_test()
{
    int errorflag = FALSE;
    int mpack_errno1, mpack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int incy = MIN_INCY; incy < MAX_INCY; incy++) {
	    for (int n = MIN_N; n < MAX_N; n++) {
		for (int m = MIN_M; m < MAX_M; m++) {
		    for (int lda = m; lda < MAX_LDA; lda++) {
#if defined VERBOSE_TEST
			printf("#m is %d, n is %d, lda is %d, incx is %d, incy is %d\n", m, n, lda, incx, incy);
#endif
			REAL_REF *A_ref;
			REAL_REF *x_ref;
			REAL_REF *y_ref;
			REAL *A;
			REAL *x;
			REAL *y;

			A_ref = new REAL_REF[matlen(lda, n)];
			x_ref = new REAL_REF[veclen(m, incx)];
			y_ref = new REAL_REF[veclen(n, incy)];
			A = new REAL[matlen(lda, n)];
			x = new REAL[veclen(m, incx)];
			y = new REAL[veclen(n, incy)];

			REAL_REF alpha_ref;
			REAL alpha;

			for (int i = 0; i < MAX_ITER; i++) {
			    set_random_vector(A_ref, A, matlen(lda, n));
			    set_random_vector(x_ref, x, veclen(m, incx));
			    set_random_vector(y_ref, y, veclen(n, incy));
			    set_random_number(alpha_ref, alpha);
			    mpack_errno = 0; blas_errno = 0;
#if defined ___MPACK_BUILD_WITH_MPFR___
			    dger_f77(&m, &n, &alpha_ref, x_ref, &incx, y_ref, &incy, A_ref, &lda);
			    mpack_errno1 = blas_errno;
#else
			    Rger(m, n, alpha_ref, x_ref, incx, y_ref, incy, A_ref, lda);
			    mpack_errno1 = mpack_errno;
#endif
			    Rger(m, n, alpha, x, incx, y, incy, A, lda);
			    mpack_errno2 = mpack_errno;
#if defined VERBOSE_TEST
			    printf("errno: mpack %d, ref %d\n", mpack_errno1, mpack_errno2);
#endif
			    if (mpack_errno1 != mpack_errno2) {
				printf("error in Mxerbla!!\n");
				exit(1);
			    }
			    REAL_REF diff = infnorm(A_ref, A, matlen(lda, n), 1);
			    if (diff > EPSILON) {
				printf("error: "); printnum(diff); printf("\n");
				errorflag = TRUE;
			    }
			    if (maxdiff < diff)
				maxdiff = diff;
			    printf("max error: "); printnum(maxdiff); printf("\n");
			}
			delete[]y;
			delete[]x;
			delete[]A;
			delete[]y_ref;
			delete[]x_ref;
			delete[]A_ref;
		    }
		}
	    }
	}
    }
    if (errorflag == TRUE) {
	printf("Rger test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Rger_test();
    printf("Rger test passed...\n");
    return (0);
}
