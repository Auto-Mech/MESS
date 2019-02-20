/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rormqr.cpp,v 1.12 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define LDT   (NBMAX+1)
#define NBMAX 64
void
Rormqr(const char *side, const char *trans, INTEGER m, INTEGER n, INTEGER k, REAL * A, INTEGER lda, REAL * tau, REAL * c, INTEGER ldc, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i;
    REAL t[LDT * NBMAX];
    char ch[3];
    INTEGER i1, i2, i3, ib, ic, jc, nb, mi = 0, ni = 0, nq, nw, iws;
    INTEGER left;
    INTEGER nbmin, iinfo;
    INTEGER notran;
    INTEGER ldwork, lwkopt;
    INTEGER lquery;
    REAL One = 1.0;

//Test the input arguments
    *info = 0;
    left = Mlsame(side, "L");
    notran = Mlsame(trans, "N");
    lquery = (lwork == -1);

//NQ is the order of Q and NW is the minimum dimension of WORK
    if (left) {
	nq = m;
	nw = n;
    } else {
	nq = n;
	nw = m;
    }
    if (!left && !Mlsame(side, "R")) {
	*info = -1;
    } else if (!notran && !Mlsame(trans, "T")) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (k < 0 || k > nq) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, nq)) {
	*info = -7;
    } else if (ldc < max((INTEGER) 1, m)) {
	*info = -10;
    } else if (lwork < max((INTEGER) 1, nw) && !lquery) {
	*info = -12;
    }

    if (*info == 0) {
//Determine the block size.  NB may be at most NBMAX, where NBMAX
//is used to define the local array T.

	ch[0] = (*side);
	ch[1] = (*trans);
	ch[2] = '\0';
	nb = min((INTEGER) NBMAX, iMlaenv(1, "Rormqr", ch, m, n, k, -1));
	lwkopt = max((INTEGER) 1, nw) * nb;
	work[0] = lwkopt;
    }

    if (*info != 0) {
	Mxerbla("Rormqr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0 || k == 0) {
	work[0] = One;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < k) {
	iws = nw * nb;
	if (lwork < iws) {
	    nb = lwork / ldwork;

	    ch[0] = (*side);
	    ch[1] = (*trans);
	    ch[2] = '\0';
	    nbmin = max((INTEGER) 2, iMlaenv(2, "Rormqr", ch, m, n, k, -1));
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= k) {
//Use unblocked code
	Rorm2r(side, trans, m, n, k, &A[0], lda, &tau[0], &c[0], ldc, &work[0], &iinfo);
    } else {

//Use blocked code
	if ((left && !notran) || (!left && notran)) {
	    i1 = 1;
	    i2 = k;
	    i3 = nb;
	} else {
	    i1 = (k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}
	if (left) {
	    ni = n;
	    jc = 1;
	} else {
	    mi = m;
	    ic = 1;
	}
	for (i = i1; i < i2; i += i3) {
	    ib = min(nb, k - i + 1);

//Form the triangular factor of the block reflector
//H = H(i) H(i+1) . . . H(i+ib-1)
	    Rlarft("Forward", "Columnwise", nq - i + 1, ib, &A[i + i * lda], lda, &tau[i], t, LDT);
	    if (left) {
//H or H' is applied to C(i:m,1:n)
		mi = m - i + 1;
		ic = i;
	    } else {
//H or H' is applied to C(1:m,i:n)
		ni = n - i + 1;
		jc = i;
	    }
//Apply H or H'
	    Rlarfb(side, trans, "Forward", "Columnwise", mi, ni, ib, &A[i + i * lda], lda, t, LDT, &c[ic + jc * ldc], ldc, &work[0], ldwork);
	}
    }
    work[0] = (REAL) double (lwkopt);
    return;
}
