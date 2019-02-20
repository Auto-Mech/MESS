/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cspr.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

/*
Based on http://www.netlib.org/blas/dspr.f
Rspr performs the symmetric rank 1 operation
 A := alpha*x*x' + A,
where alpha is a real scalar, x is an n element vector and A is an
n by n symmetric matrix, supplied in packed form.
*/

#include <mblas.h>

void Cspr(const char *uplo, INTEGER n, COMPLEX alpha, COMPLEX * x, INTEGER incx, COMPLEX * AP)
{
    INTEGER ix, jx, kx, kk;
    COMPLEX temp;
    REAL Zero = 0.0;

//Test the input parameters.
    INTEGER info = 0;

    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L"))
	info = 1;
    else if (n < 0)
	info = 2;
    else if (incx == 0)
	info = 5;
    if (info != 0) {
	Mxerbla("Cspr  ", info);
	return;
    }
    //quick return if possible.
    if ((n == 0) || (alpha == Zero))
	return;

    //set the start point in x if the increment is not unity.
    if (incx <= 0)
	kx = -(n - 1) * incx;
    else
	kx = 0;

    //start the operations. in this version the elements of the array ap
    //are accessed sequentially with one pass through ap.
    kk = 0;
    if (Mlsame(uplo, "U")) {
	//Form  A  when upper triangle is stored in AP.
	jx = kx;
	for (INTEGER j = 0; j < n; j++) {
	    if (x[jx] != Zero) {
		temp = alpha * x[jx];
		ix = kx;
		for (INTEGER k = kk; k < kk + j + 1; k++) {
		    AP[k] = AP[k] + x[ix] * temp;
		    ix = ix + incx;
		}
	    }
	    jx = jx + incx;
	    kk = kk + j + 1;
	}
    } else {
	//Form  A  when lower triangle is stored in AP.
	jx = kx;
	for (INTEGER j = 0; j < n; j++) {
	    if (x[jx] != Zero) {
		temp = alpha * x[jx];
		ix = jx;
		for (INTEGER k = kk; k < kk + n - j; k++) {
		    AP[k] = AP[k] + x[ix] * temp;
		    ix = ix + incx;
		}
	    }
	    jx = jx + incx;
	    kk = kk + n - j;
	}
    }
}
