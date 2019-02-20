/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rdot_dd.cpp,v 1.3 2010/08/07 05:50:09 nakatamaho Exp $
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

#define ___MPACK_BUILD_WITH_DD___
#include <stdio.h>
#include <mblas_dd.h>
#include <mpack_benchmark.h>

#define MAXLOOP 10

int main()
{
    mpackint n;
    mpackint incx = 1, incy = 1;
    dd_real alpha, dummy;
    double elapsedtime, t1, t2;

    for (n = 10000000; n < 10500000; n = n + 100) {
	printf("n: %d\n", (int) n);
	dd_real *x = new dd_real[n];
	dd_real *y = new dd_real[n];
	for (int i = 0; i < n; i++) {
	    x[i] = mpf_randomnumber(dummy);
	    y[i] = mpf_randomnumber(dummy);
	}
	t1 = gettime();
	for (int p = 0; p < MAXLOOP; p++) {
	    dummy = Rdot(n, x, incx, y, incy);
	}
	t2 = gettime();
	elapsedtime = (t2 - t1);
	printf("time : %f\n", elapsedtime);
	printf("Mflops : %f\n", 2.0 * (double) n * MAXLOOP / elapsedtime / (1000 * 1000));
	delete[]y;
	delete[]x;
    }
}
