/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsyevx.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rsyevx(const char *jobz, const char *range, const char *uplo, INTEGER n,
	    REAL * A, INTEGER lda, REAL vl, REAL vu, INTEGER il, INTEGER iu, REAL abstol, INTEGER * m, REAL * w,
	    REAL * z, INTEGER ldz, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * ifail, INTEGER * info)
{
    INTEGER i, j, nb, jj;
    REAL eps, vll = 0.0, vuu = 0.0, tmp1;
    INTEGER indd, inde;
    REAL anrm;
    INTEGER imax;
    REAL rmin, rmax;
    INTEGER test;
    INTEGER itmp1, indee;
    REAL sigma = 0.0;
    INTEGER iinfo;
    char order;
    INTEGER lower, wantz;
    INTEGER alleig, indeig;
    INTEGER iscale, indibl;
    INTEGER valeig;
    REAL safmin;
    REAL abstll, bignum;
    INTEGER indtau, indisp, indiwo, indwkn;
    INTEGER llwrkn, llwork, nsplit;
    REAL smlnum;
    INTEGER lwkopt = 0, lwkmin, indwrk;
    INTEGER lquery;
    REAL mtemp1, mtemp2;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
    lower = Mlsame(uplo, "L");
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    lquery = lwork == -1;

    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(alleig || valeig || indeig)) {
	*info = -2;
    } else if (!(lower || Mlsame(uplo, "U"))) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else {
	if (valeig) {
	    if (n > 0 && vu <= vl) {
		*info = -8;
	    }
	} else if (indeig) {
	    if (il < 1 || il > max((INTEGER) 1, n)) {
		*info = -9;
	    } else if (iu < min(n, il) || iu > n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
	if (ldz < 1 || (wantz && ldz < n)) {
	    *info = -15;
	}
    }
    if (*info == 0) {
	if (n <= 1) {
	    lwkmin = 1;
	    work[1] = lwkmin;
	} else {
	    lwkmin = n * 8;
	    nb = iMlaenv(1, "Rsytrd", uplo, n, -1, -1, -1);
	    nb = max(nb, iMlaenv(1, "Rormtr", uplo, n, -1, -1, -1));
	    lwkopt = max(lwkmin, (nb + 3) * n);
	    work[1] = lwkopt;
	}
	if (lwork < lwkmin && !lquery) {
	    *info = -17;
	}
    }
    if (*info != 0) {
	Mxerbla("Rsyevx", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    (*m) = 0;
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (alleig || indeig) {
	    (*m) = 1;
	    w[1] = A[lda + 1];
	} else {
	    if (vl < A[lda + 1] && vu >= A[lda + 1]) {
		(*m) = 1;
		w[1] = A[lda + 1];
	    }
	}
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
    mtemp1 = sqrt(bignum), mtemp2 = One / sqrt(sqrt(safmin));
    rmax = min(mtemp1, mtemp2);
//Scale matrix to allowable range, if necessary.
    iscale = 0;
    abstll = abstol;
    if (valeig) {
	vll = vl;
	vuu = vu;
    }
    anrm = Rlansy("M", uplo, n, &A[0], lda, &work[0]);
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	if (lower) {
	    for (j = 0; j < n; j++) {
		Rscal(n - j + 1, sigma, &A[j + j * lda], 1);
	    }
	} else {
	    for (j = 0; j < n; j++) {
		Rscal(j, sigma, &A[j * lda], 1);
	    }
	}
	if (abstol > Zero) {
	    abstll = abstol * sigma;
	}
	if (valeig) {
	    vll = vl * sigma;
	    vuu = vu * sigma;
	}
    }
//Call DSYTRD to reduce symmetric matrix to tridiagonal form.
    indtau = 1;
    inde = indtau + n;
    indd = inde + n;
    indwrk = indd + n;
    llwork = lwork - indwrk + 1;
    Rsytrd(uplo, n, &A[0], lda, &work[indd], &work[inde], &work[indtau], &work[indwrk], llwork, &iinfo);
//If all eigenvalues are desired and ABSTOL is less than or equal to
//zero, then call DSTERF or DORGTR and SSTEQR.  If this fails for
//some eigenvalue, then try DSTEBZ.
    test = MFALSE;
    if (indeig) {
	if (il == 1 && iu == n) {
	    test = MTRUE;
	}
    }
    if ((alleig || test) && abstol <= Zero) {
	Rcopy(n, &work[indd], 1, &w[1], 1);
	indee = indwrk + (n * 2);
	if (!wantz) {
	    Rcopy(n - 1, &work[inde], 1, &work[indee], 1);
	    Rsterf(n, &w[1], &work[indee], info);
	} else {
	    Rlacpy("A", n, n, &A[0], lda, &z[0], ldz);
	    Rorgtr(uplo, n, &z[0], ldz, &work[indtau], &work[indwrk], llwork, &iinfo);
	    Rcopy(n - 1, &work[inde], 1, &work[indee], 1);
	    Rsteqr(jobz, n, &w[1], &work[indee], &z[0], ldz, &work[indwrk], info);
	    if (*info == 0) {
		for (i = 0; i < n; i++) {
		    ifail[i] = 0;
		}
	    }
	}
	if (*info == 0) {
	    (*m) = n;
	    goto L40;
	}
	*info = 0;
    }
//Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
    if (wantz) {
	order = 'B';
    } else {
	order = 'E';
    }
    indibl = 0;
    indisp = indibl + n;
    indiwo = indisp + n;
    Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &work[indd], &work[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk], &iwork[indiwo], info);
    if (wantz) {
	Rstein(n, &work[indd], &work[inde], *m, &w[1], &iwork[indibl], &iwork[indisp], &z[0], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], info);
//Apply orthogonal matrix used in reduction to tridiagonal
//form to eigenvectors returned by DSTEIN.
	indwkn = inde;
	llwrkn = lwork - indwkn + 1;
	Rormtr("L", uplo, "N", n, *m, &A[0], lda, &work[indtau], &z[0], ldz, &work[indwkn], llwrkn, &iinfo);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
  L40:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = (*m);
	} else {
	    imax = *info - 1;
	}
	Rscal(imax, One / sigma, &w[1], 1);
    }
//If eigenvalues are not in order, then sort them, along with
//eigenvectors.
    if (wantz) {
	for (j = 0; j < (*m) - 1; j++) {
	    i = 0;
	    tmp1 = w[j];
	    for (jj = j + 1; jj <= (*m); jj++) {
		if (w[jj] < tmp1) {
		    i = jj;
		    tmp1 = w[jj];
		}
	    }
	    if (i != 0) {
		itmp1 = iwork[indibl + i - 1];
		w[i] = w[j];
		iwork[indibl + i - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		Rswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
		if (*info != 0) {
		    itmp1 = ifail[i];
		    ifail[i] = ifail[j];
		    ifail[j] = itmp1;
		}
	    }

	}
    }
//Set WORK(1) to optimal workspace size.
    work[1] = lwkopt;
    return;
}
