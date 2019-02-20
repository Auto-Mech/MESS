/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsygvx.cpp,v 1.3 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsygvx(INTEGER itype, const char *jobz, const char *range, const char *uplo, INTEGER n, REAL * A, INTEGER lda, REAL * B,
	    INTEGER ldb, REAL vl, REAL vu, INTEGER il, INTEGER iu, REAL abstol, INTEGER * m, REAL * w, REAL * z, INTEGER ldz,
	    REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * ifail, INTEGER * info)
{
    INTEGER nb;
    char trans;
    INTEGER upper, wantz, alleig, indeig, valeig;
    INTEGER lwkmin;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL One = 1.0;

//Test the input parameters.
    upper = Mlsame(uplo, "U");
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    lquery = lwork == -1;

    *info = 0;
    if (itype < 1 || itype > 3) {
	*info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -2;
    } else if (!(alleig || valeig || indeig)) {
	*info = -3;
    } else if (!(upper || Mlsame(uplo, "L"))) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else {
	if (valeig) {
	    if (n > 0 && vu <= vl) {
		*info = -11;
	    }
	} else if (indeig) {
	    if (il < 1 || il > max((INTEGER) 1, n)) {
		*info = -12;
	    } else if (iu < min(n, il) || iu > n) {
		*info = -13;
	    }
	}
    }
    if (*info == 0) {
	if (ldz < 1 || (wantz && ldz < n)) {
	    *info = -18;
	}
    }
    if (*info == 0) {
	lwkmin = max((INTEGER) 1, n * 8);
	nb = iMlaenv(1, "Rsytrd", uplo, n, -1, -1, -1);
	lwkopt = max(lwkmin, (nb + 3) * n);
	work[1] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -20;
	}
    }
    if (*info != 0) {
	Mxerbla("Rsygvx", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    m = 0;
    if (n == 0) {
	return;
    }
//Form a Cholesky factorization of B.
    Rpotrf(uplo, n, &B[0], ldb, info);
    if (*info != 0) {
	*info = n + *info;
	return;
    }
//Transform problem to standard eigenvalue problem and solve.
    Rsygst(itype, uplo, n, &A[0], lda, &B[0], ldb, info);
    Rsyevx(jobz, range, uplo, n, &A[0], lda, vl, vu, il, iu, abstol, m, &w[1], &z[0], ldz, &work[0], lwork, &iwork[1], &ifail[1], info);
    if (wantz) {
//Backtransform eigenvectors to the original problem.
	if (*info > 0) {
	    *m = *info - 1;
	}
	if (itype == 1 || itype == 2) {
//For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
//backtransform eigenvectors: x = inv(L)'*y or inv(U)*y
	    if (upper) {
		trans = 'N';
	    } else {
		trans = 'T';
	    }
	    Rtrsm("Left", uplo, &trans, "Non-unit", n, *m, One, &B[0], ldb, &z[0], ldz);
	} else if (itype == 3) {
//For B*A*x=(lambda)*x;
//backtransform eigenvectors: x = L*y or U'*y
	    if (upper) {
		trans = 'T';
	    } else {
		trans = 'N';
	    }
	    Rtrmm("Left", uplo, &trans, "Non-unit", n, *m, One, &B[0], ldb, &z[0], ldz);
	}
    }
//Set WORK(1) to optimal workspace size.
    work[1] = lwkopt;
    return;
}
