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

#define MIN_INCX -10
#define MAX_INCX  10
#define MAX_N     100
#define MAX_ITER  10

REAL_REF maxdiff = 0.0;

void Rasum_test()
{
    int errorflag = FALSE;
    REAL_REF dtemp;
    REAL rtemp;

    for (int incx = MIN_INCX; incx <= MAX_INCX; incx++) {
	for (int n = 0; n < MAX_N; n++) {
#if defined VERBOSE_TEST
	    printf("# n:%d incx:%d\n", n, incx);
#endif
	    REAL_REF *x_ref = new REAL_REF[veclen(n, incx)];
	    REAL *x = new REAL[veclen(n, incx)];

	    int j = 0;
	    while (j < MAX_ITER) {
		set_random_vector(x_ref, x, veclen(n, incx));
#if defined ___MPACK_BUILD_WITH_MPFR___
		dtemp = dasum_f77(&n, x_ref, &incx);
#else
		dtemp = Rasum(n, x_ref, incx);
#endif
		rtemp = Rasum(n, x, incx);

		REAL_REF diff = rtemp - dtemp;
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
    if (errorflag != FALSE) {
	printf("Rasum test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Rasum_test();
    printf("Rasum test passed...\n");
}
