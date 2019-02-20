/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarzb.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>

void
Clarzb(const char *side, const char *trans, const char *direct,
       const char *storev, INTEGER m, INTEGER n, INTEGER k, INTEGER l, COMPLEX * v, INTEGER ldv, COMPLEX * t, INTEGER ldt, COMPLEX * c, INTEGER ldc, COMPLEX * work, INTEGER ldwork)
{
    INTEGER i, j, info;
    char transt;
    REAL One = 1.0;

    if (m <= 0 || n <= 0) {
	return;
    }
//Check for currently supported options
    info = 0;
    if (!Mlsame(direct, "B")) {
	info = -3;
    } else if (!Mlsame(storev, "R")) {
	info = -4;
    }
    if (info != 0) {
	Mxerbla("Clarzb", -(info));
	return;
    }

    if (Mlsame(trans, "N")) {
	transt = 'C';
    } else {
	transt = 'N';
    }

    if (Mlsame(side, "L")) {
//Form  H * C  or  H' * C
//W( 1:n, 1:k ) = conjg( C( 1:k, 1:n )' )
	for (j = 0; j < k; j++) {
	    Ccopy(n, &c[j + ldc], ldc, &work[j * ldwork + 1], 1);
	}
//W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
//                conjg( C( m-l+1:m, 1:n )' ) * V( 1:k, 1:l )'
	if (l > 0) {
	    Cgemm("Transpose", "Conjugate transpose", n, k, l, (COMPLEX) One, &c[m - l + 1 + ldc], ldc, &v[0], ldv, (COMPLEX) One, &work[0], ldwork);
	}
//W( 1:n, 1:k ) = W( 1:n, 1:k ) * T'  or  W( 1:m, 1:k ) * T
	Ctrmm("Right", "Lower", &transt, "Non-unit", n, k, One, &t[0], ldt, &work[0], ldwork);
//C( 1:k, 1:n ) = C( 1:k, 1:n ) - conjg( W( 1:n, 1:k )' )
	for (j = 0; j < n; j++) {
	    for (i = 0; i < k; i++) {
		c[i + j * ldc] = c[i + j * ldc] - work[j + i * ldwork];
	    }
	}
//C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
//            conjg( V( 1:k, 1:l )' ) * conjg( W( 1:n, 1:k )' )
	if (l > 0) {
	    Cgemm("Transpose", "Transpose", l, n, k, (COMPLEX) - One, &v[0], ldv, &work[0], ldwork, (COMPLEX) One, &c[m - l + 1 + ldc], ldc);
	}
    } else if (Mlsame(side, "R")) {
//Form  C * H  or  C * H'
//W( 1:m, 1:k ) = C( 1:m, 1:k )
	for (j = 0; j < k; j++) {
	    Ccopy(m, &c[j * ldc + 1], 1, &work[j * ldwork + 1], 1);
	}
//W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
//                C( 1:m, n-l+1:n ) * conjg( V( 1:k, 1:l )' )
	if (l > 0) {
	    Cgemm("No transpose", "Transpose", m, k, l, (COMPLEX) One, &c[(n - l + 1) * ldc + 1], ldc, &v[0], ldv, (COMPLEX) One, &work[0], ldwork);
	}
//W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
//                W( 1:m, 1:k ) * conjg( T' )
	for (j = 0; j < k; j++) {
	    Clacgv(k - j + 1, &t[j + j * ldt], 1);
	}
	Ctrmm("Right", "Lower", trans, "Non-unit", m, k, (COMPLEX) One, &t[0], ldt, &work[0], ldwork);
	for (j = 0; j < k; j++) {
	    Clacgv(k - j + 1, &t[j + j * ldt], 1);
	}
//C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )
	for (j = 0; j < k; j++) {
	    for (i = 0; i < m; i++) {
		c[i + j * ldc] = c[i + j * ldc] - work[i + j * ldwork];
	    }
	}
//C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
//                    W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )
	for (j = 0; j < l; j++) {
	    Clacgv(k, &v[j * ldv + 1], 1);
	}
	if (l > 0) {
	    Cgemm("No transpose", "No transpose", m, l, k, (COMPLEX) - One, &work[0], ldwork, &v[0], ldv, (COMPLEX) One, &c[(n - l + 1) * ldc + 1], ldc);
	}
	for (j = 0; j < l; j++) {
	    Clacgv(k, &v[j * ldv + 1], 1);
	}
    }
    return;
}
