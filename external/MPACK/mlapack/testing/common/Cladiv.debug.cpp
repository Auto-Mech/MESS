/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Cladiv.debug.cpp,v 1.4 2010/08/07 05:50:10 nakatamaho Exp $
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

REAL_REF maxdiff = 0.0;

void Cladiv_test()
{
    int errorflag = FALSE;
    REAL_REF diff;
    COMPLEX_REF x_ref, y_ref, ret_ref;
    COMPLEX x, y, ret;

    int count = 100;
    while (count--) {
	set_random_number(x_ref, x);
	set_random_number(y_ref, y);

#if defined ___MPACK_BUILD_WITH_MPFR___
	ret_ref = zladiv_f77(&x_ref, &y_ref);
#else
	ret_ref = Cladiv(x_ref, y_ref);
#endif
	ret = Cladiv(x, y);

#if defined VERBOSE_TEST
        cout << "x_ref  " ; printnum (x_ref) ; cout << endl;
        cout << "y_ref  " ; printnum (y_ref) ; cout << endl;
        cout << "ret_ref" ; printnum (ret_ref) ; cout << endl;

        cout << "x      " ; printnum (x) ; cout << endl;
        cout << "y      " ; printnum (y) ; cout << endl;
        cout << "ret    " ; printnum (ret) ; cout << endl;
#endif
        diff = abs (ret_ref - ret);
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("error1: "); printnum(diff); printf("\n");
	}
        if (maxdiff < diff)
	    maxdiff = diff;
        printf("max error: "); printnum(maxdiff); printf("\n");
    }
    if (errorflag == TRUE) {
	printf("Cladiv test failed...\n");
	exit(1);
    }
    printf("max error: "); printnum(maxdiff); printf("\n");
}

int main(int argc, char *argv[])
{
    Cladiv_test();
    printf("Cladiv test passed...\n");
    return (0);
}
