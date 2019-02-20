/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rptcon.cpp,v 1.7 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rptcon(INTEGER n, REAL * d, REAL * e, REAL anorm, REAL * rcond, REAL * work, INTEGER * info)
{
    INTEGER i, ix;
    REAL ainvnm;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    if (n < 0) {
	*info = -1;
    } else if (anorm < Zero) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rptcon", -(*info));
	return;
    }
//Quick return if possible
    *rcond = Zero;
    if (n == 0) {
	*rcond = One;
	return;
    } else if (anorm == Zero) {
	return;
    }
// Check that D(1:N) is positive.
    for (i = 0; i < n; i++) {
	if (d[i] <= Zero) {
	    return;
	}

    }

//Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
// m(i,j) =  abs(A(i,j)), i = j,
// m(i,j) = -abs(A(i,j)), i !=j,
//and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.
//Solve M(L) * x = e.
    work[0] = One;
    for (i = 1; i < n; i++) {
	work[i] = work[i - 1] * abs(e[i - 1]) + One;
    }
//Solve D * M(L)' * x = b.
    work[n] /= d[n];
    for (i = n - 2; i >= 0; i--) {
	work[i] = work[i] / d[i] + work[i + 1] * abs(e[i]);

    }
//Compute AINVNM = max(x(i)), 1<=i<=n.

    ix = iRamax(n, &work[0], 1);
    ainvnm = abs(work[ix]);

//Compute the reciprocal condition number.
    if (ainvnm != Zero) {
	*rcond = One / ainvnm / anorm;
    }
    return;
}
