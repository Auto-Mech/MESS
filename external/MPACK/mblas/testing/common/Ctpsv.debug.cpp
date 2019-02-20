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

#define MIN_N -1
#define MAX_N 10
#define MIN_INCX -3
#define MAX_INCX 3
#define MAX_ITER 5

REAL_REF maxdiff = 0.0;

void Ctpsv_test2(const char *uplo, const char *trans, const char *diag)
{
    int errorflag = FALSE;
    int mpack_errno1, mpack_errno2;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = MIN_N; n < MAX_N; n++) {
#if defined VERBOSE_TEST
	    printf("#n is %d, incx is %d ", n, incx);
	    printf("uplo is %s trans is %s, diag is %s \n", uplo, trans, diag);
#endif
	    COMPLEX_REF *AP_ref = new COMPLEX_REF[vecplen(n)];
	    COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
	    COMPLEX *AP = new COMPLEX[vecplen(n)];
	    COMPLEX *x = new COMPLEX[veclen(n, incx)];

	    for (int iter = 0; iter < MAX_ITER; iter++) {
		set_random_vector(AP_ref, AP, vecplen(n));
		set_random_vector(x_ref, x, veclen(n, incx));

		mpack_errno = 0; blas_errno = 0;
#if defined ___MPACK_BUILD_WITH_MPFR___
		ztpsv_f77(uplo, trans, diag, &n, AP_ref, x_ref, &incx);
		mpack_errno1 = blas_errno;
#else
		Ctpsv(uplo, trans, diag, n, AP_ref, x_ref, incx);
		mpack_errno1 = mpack_errno;
#endif
		Ctpsv(uplo, trans, diag, n, AP, x, incx);
		mpack_errno2 = mpack_errno;

#if defined VERBOSE_TEST
		printf("errno: mpack %d, ref %d\n", mpack_errno1, mpack_errno2);
#endif
		if (mpack_errno1 != mpack_errno2) {
		    printf("error in Mxerbla!!\n");
		    exit(1);
		}
		REAL_REF diff = infnorm(x_ref, x, veclen(n, incx), 1);
		if (diff > EPSILON10) {
		    printf("error: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
		if (maxdiff < diff)
		    maxdiff = diff;
		printf("max error: "); printnum(maxdiff); printf("\n");
	    }
	    delete[]AP_ref;
	    delete[]AP;
	    delete[]x_ref;
	    delete[]x;
	}
    }
    if (errorflag == TRUE) {
	printf("Ctpsv test failed...\n");
	exit(1);
    }
    printf("max error: ");
    printnum(maxdiff);
    printf("\n");
}

void Ctpsv_test()
{
    Ctpsv_test2("U", "N", "U");
    Ctpsv_test2("U", "N", "N");
    Ctpsv_test2("U", "T", "U");
    Ctpsv_test2("U", "T", "N");
    Ctpsv_test2("U", "C", "U");
    Ctpsv_test2("U", "C", "N");

    Ctpsv_test2("L", "N", "U");
    Ctpsv_test2("L", "N", "N");
    Ctpsv_test2("L", "T", "U");
    Ctpsv_test2("L", "T", "N");
    Ctpsv_test2("L", "C", "U");
    Ctpsv_test2("L", "C", "N");
}

int main(int argc, char *argv[])
{
    Ctpsv_test();
    printf("Ctpsv test passed...\n");
    return (0);
}
