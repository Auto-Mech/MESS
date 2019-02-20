/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: dgemm.cpp,v 1.5 2010/08/19 01:29:39 nakatamaho Exp $
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

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#define F77_FUNC(name,NAME) name ## _
#include <blas.h>
#include <mpack_benchmark.h>

#define MAXLOOP 3

int main()
{
    int n;
    double alpha, beta, dummy;
    double elapsedtime, t1, t2;

    for (n = 3000; n < 10000; n = n + 100) {
	printf("n: %d\n", (int) n);
	double *A = new double[n * n];
	double *B = new double[n * n];
	double *C = new double[n * n];
	for (int i = 0; i < n; i++) {
	    for (int j = 0; j < n; j++) {
		A[i * n + j] = mpf_randomnumber(dummy);
		B[i * n + j] = mpf_randomnumber(dummy);
		C[i * n + j] = mpf_randomnumber(dummy);
	    }
	}
	t1 = gettime();
	for (int p = 0; p < MAXLOOP; p++) {
	    dgemm_f77("n", "n", &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
	}
	t2 = gettime();
	elapsedtime = (t2 - t1);
	printf("time : %f\n", elapsedtime);
	printf("Mflops : %f\n",
	       (2.0 * (double) n * (double) n * (double) n +
		2.0 * (double) n * (double) n) * MAXLOOP / elapsedtime / (1000 * 1000));
	delete[]C;
	delete[]B;
	delete[]A;
    }
}
