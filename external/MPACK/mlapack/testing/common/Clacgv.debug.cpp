/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Clacgv.debug.cpp,v 1.9 2010/08/07 05:50:10 nakatamaho Exp $
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

#define MIN_INCX -10
#define MAX_INCX  10
#define MAX_N     100
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Clacgv_test()
{
    int errorflag = FALSE;
    REAL_REF diff;
    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = 0; n < MAX_N; n++) {
#if defined VERBOSE_TEST
	    printf("# n:%d incx:%d\n", n, incx);
#endif
	    COMPLEX_REF *x_ref = new COMPLEX_REF[veclen(n, incx)];
	    COMPLEX *x = new COMPLEX[veclen(n, incx)];
	    int j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(x_ref, x, veclen(n, incx));
#if defined ___MPACK_BUILD_WITH_MPFR___
		zlacgv_f77(&n, x_ref, &incx);
#else
		Clacgv(n, x_ref, incx);
#endif
		Clacgv(n, x, incx);

		diff = infnorm(x_ref, x, veclen(n, incx), 1);
		if (diff > EPSILON) {
		    printf("error: "); printnum(diff); printf("\n");
		    errorflag = TRUE;
		}
	        if (maxdiff < diff)
		    maxdiff = diff;
	        printf("max error: "); printnum(maxdiff); printf("\n");
		j++;
	    }
	    delete[]x_ref;
	    delete[]x;
	}
    }
    if (errorflag == TRUE) {
	printf("Clacgv test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Clacgv_test();
    printf("Clacgv test passed...\n");
    return (0);
}
