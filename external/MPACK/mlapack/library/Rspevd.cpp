/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rspevd.cpp,v 1.6 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rspevd(const char *jobz, const char *uplo, INTEGER n, REAL * ap, REAL * w, REAL * z, INTEGER ldz, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    REAL eps;
    INTEGER inde;
    REAL anrm, rmin, rmax;
    REAL sigma = 0.0;
    INTEGER iinfo, lwmin;
    INTEGER wantz;
    INTEGER iscale;
    REAL safmin;
    REAL bignum;
    INTEGER indtau;
    INTEGER indwrk, liwmin;
    INTEGER llwork;
    REAL smlnum;
    INTEGER lquery;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    lquery = lwork == -1 || liwork == -1;
    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(Mlsame(uplo, "U") || Mlsame(uplo, "L"))) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -7;
    }
    if (*info == 0) {
	if (n <= 1) {
	    liwmin = 1;
	    lwmin = 1;
	} else {
	    if (wantz) {
		liwmin = n * 5 + 3;
		lwmin = n * 6 + 1 + n * n;
	    } else {
		liwmin = 1;
		lwmin = n * 2;
	    }
	}
	iwork[1] = liwmin;
	work[1] = lwmin;

	if (lwork < lwmin && !lquery) {
	    *info = -9;
	} else if (liwork < liwmin && !lquery) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	Mxerbla("Rspevd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
	w[1] = ap[1];
	if (wantz) {
	    z[ldz + 1] = One;
	}
	return;
    }
//Get machine constants.
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
    smlnum = safmin / eps;
    bignum = One / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
//Scale matrix to allowable range, if necessary.
    anrm = Rlansp("M", uplo, n, &ap[1], &work[0]);
    iscale = 0;
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	Rscal(n * (n + 1) / 2, sigma, &ap[1], 1);
    }
//Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
    inde = 1;
    indtau = inde + n;
    Rsptrd(uplo, n, &ap[1], &w[1], &work[inde], &work[indtau], &iinfo);
//For eigenvalues only, call DSTERF.  For eigenvectors, first call
//DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
//tridiagonal matrix, then call DOPMTR to multiply it by the
//Householder transformations represented in AP.
    if (!wantz) {
	Rsterf(n, &w[1], &work[inde], info);
    } else {
	indwrk = indtau + n;
	llwork = lwork - indwrk + 1;
	Rstedc("I", n, &w[1], &work[inde], &z[0], ldz, &work[indwrk]
	       , llwork, &iwork[1], liwork, info);
	Ropmtr("L", uplo, "N", n, n, &ap[1], &work[indtau], &z[0], ldz, &work[indwrk], &iinfo);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
    if (iscale == 1) {
	Rscal(n, One / sigma, &w[1], 1);
    }
    work[1] = lwmin;
    iwork[1] = liwmin;
    return;
}
