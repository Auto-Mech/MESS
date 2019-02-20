/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgerqf.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgerqf(INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * tau, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, k, ib, nb = 0, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    INTEGER ldwork = 0;
    INTEGER lwkopt;
    INTEGER lquery;

    *info = 0;
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    }
    if (*info == 0) {
	k = min(m, n);
	if (k == 0) {
	    lwkopt = 1;
	} else {
	    nb = iMlaenv(1, "Cgerqf", " ", m, n, -1, -1);
	    lwkopt = m * nb;
	}
	work[1] = lwkopt;
	if (lwork < max((INTEGER) 1, m) && !lquery) {
	    *info = -7;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgerqf", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (k == 0) {
	return;
    }
    nbmin = 2;
    nx = 1;
    iws = m;
    if (nb > 1 && nb < k) {
//Determine when to cross over from blocked to unblocked code.
	nx = max((INTEGER) 0, iMlaenv(3, "Cgerqf", " ", m, n, -1, -1));
	if (nx < k) {
//Determine if workspace is large enough for blocked code.
	    ldwork = m;
	    iws = ldwork * nb;
	    if (lwork < iws) {
//Not enough workspace to use optimal NB:  reduce NB and
//determine the minimum value of NB.
		nb = lwork / ldwork;
		nbmin = max((INTEGER) 2, iMlaenv(2, "Cgerqf", " ", m, n, -1, -1));
	    }
	}
    }
    if (nb >= nbmin && nb < k && nx < k) {
//Use blocked code initially.
//The last kk rows are handled by the block method.
	ki = (k - nx - 1) / nb * nb;
	kk = min(k, ki + nb);
	for (i = k - kk + ki + 1; i <= -nb; i = i - nb) {
//Compute the RQ factorization of the current block
//A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)
	    ib = 1;
	    Cgerq2(min(k - i + 1, nb), n - k + i + ib - 1, &A[m - k + i + lda], lda, &tau[i], &work[0], &iinfo);
	    if (m - k + i > 1) {
//Form the triangular factor of the block reflector
//H = H(i+ib-1) . . . H(i+1) H(i)
		Clarft("Backward", "Rowwise", n - k + i + ib - 1, ib, &A[m - k + i + lda], lda, &tau[i], &work[0], ldwork);
//Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
		Clarfb("Right", "No transpose", "Backward", "Rowwise", m - k + i - 1, n - k + i + ib - 1, ib, &A[m - k + i + lda],
		       lda, &work[0], ldwork, &A[0], lda, &work[ib + 1], ldwork);
	    }
	}
	mu = m - k + i + nb - 1;
	nu = n - k + i + nb - 1;
    } else {
	mu = m;
	nu = n;
    }
//Use unblocked code to factor the last or only block
    if (mu > 0 && nu > 0) {
	Cgerq2(mu, nu, A, lda, tau, work, &iinfo);
    }
    work[1] = iws;
    return;
}
