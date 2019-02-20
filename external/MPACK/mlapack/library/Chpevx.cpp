/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chpevx.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chpevx(const char *jobz, const char *range, const char *uplo, INTEGER n,
	    COMPLEX * ap, REAL vl, REAL vu, INTEGER il,
	    INTEGER iu, REAL abstol, INTEGER * m, REAL * w, COMPLEX * z, INTEGER ldz, COMPLEX * work, REAL * rwork, INTEGER * iwork, INTEGER * ifail, INTEGER * info)
{
    INTEGER i, j, jj;
    REAL eps, vll, vuu, tmp1;
    INTEGER indd, inde;
    REAL anrm;
    INTEGER imax;
    REAL rmin, rmax;
    INTEGER test;
    INTEGER itmp1, indee;
    REAL sigma = 0.0;
    INTEGER iinfo;
    char order;
    INTEGER wantz;
    INTEGER alleig, indeig;
    INTEGER iscale, indibl;
    INTEGER valeig;
    REAL safmin;
    REAL abstll, bignum;
    INTEGER indiwk, indisp, indtau;
    INTEGER indrwk, indwrk, nsplit;
    REAL smlnum;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");

    *info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(alleig || valeig || indeig)) {
	*info = -2;
    } else if (!(Mlsame(uplo, "L") || Mlsame(uplo, "U"))) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else {
	if (valeig) {
	    if (n > 0 && vu <= vl) {
		*info = -7;
	    }
	} else if (indeig) {
	    if (il < 1 || il > max((INTEGER) 1, n)) {
		*info = -8;
	    } else if (iu < min(n, il) || iu > n) {
		*info = -9;
	    }
	}
    }
    if (*info == 0) {
	if (ldz < 1 || (wantz && ldz < n)) {
	    *info = -14;
	}
    }
    if (*info != 0) {
	Mxerbla("Chpevx", -(*info));
	return;
    }
//Quick return if possible
    m = 0;
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = ap[1].real();
	} else {
	    if (vl < ap[1].real() && vu >= ap[1].real()) {
		*m = 1;
		w[1] = ap[1].real();
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
    } else {
	vll = Zero;
	vuu = Zero;
    }
    anrm = Clanhp("M", uplo, n, &ap[1], &rwork[1]);
    if (anrm > Zero && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {
	CRscal(n * (n + 1) / 2, sigma, &ap[1], 1);
	if (abstol > Zero) {
	    abstll = abstol * sigma;
	}
	if (valeig) {
	    vll = vl * sigma;
	    vuu = vu * sigma;
	}
    }
//Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.
    indd = 1;
    inde = indd + n;
    indrwk = inde + n;
    indtau = 1;
    indwrk = indtau + n;
    Chptrd(uplo, n, &ap[1], &rwork[indd], &rwork[inde], &work[indtau], &iinfo);
//If all eigenvalues are desired and ABSTOL is less than or equal
//to zero, then call DSTERF or ZUPGTR and ZSTEQR.  If this fails
//for some eigenvalue, then try DSTEBZ.
    test = MFALSE;
    if (indeig) {
	if (il == 1 && iu == n) {
	    test = MTRUE;
	}
    }
    if ((alleig || test) && abstol <= Zero) {
	Rcopy(n, &rwork[indd], 1, &w[1], 1);
	indee = indrwk + (n * 2);
	if (!wantz) {
	    Rcopy(n - 1, &rwork[inde], 1, &rwork[indee], 1);
	    Rsterf(n, &w[1], &rwork[indee], info);
	} else {
	    Cupgtr(uplo, n, &ap[1], &work[indtau], &z[0], ldz, &work[indwrk], &iinfo);
	    Rcopy(n - 1, &rwork[inde], 1, &rwork[indee], 1);
	    Csteqr(jobz, n, &w[1], &rwork[indee], &z[0], ldz, &rwork[indrwk], info);
	    if (*info == 0) {
		for (i = 0; i < n; i++) {
		    ifail[i] = 0;
		}
	    }
	}
	if (*info == 0) {
	    *m = n;
	    goto L20;
	}
	*info = 0;
    }
//Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
    if (wantz) {
	order = 'B';
    } else {
	order = 'E';
    }
    indibl = 0;

    indisp = indibl + n;
    indiwk = indisp + n;
    Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &rwork[indd], &rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &rwork[indrwk], &iwork[indiwk], info);
    if (wantz) {
	Cstein(n, &rwork[indd], &rwork[inde], *m, &w[1], &iwork[indibl], &iwork[indisp], &z[0], ldz, &rwork[indrwk], &iwork[indiwk], &ifail[1], info);
//Apply unitary matrix used in reduction to tridiagonal
//form to eigenvectors returned by ZSTEIN.
	indwrk = indtau + n;
	Cupmtr("L", uplo, "N", n, *m, &ap[1], &work[indtau], &z[0], ldz, &work[indwrk], info);
    }
//If matrix was scaled, then rescale eigenvalues appropriately.
  L20:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	Rscal(imax, One / sigma, &w[1], 1);
    }
//If eigenvalues are not in order, then sort them, along with
//eigenvectors.
    if (wantz) {
	for (j = 0; j < *m - 1; j++) {
	    i = 0;
	    tmp1 = w[j];
	    for (jj = j + 1; jj <= *m; jj++) {
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
		Cswap(n, &z[i * ldz + 1], 1, &z[j * ldz + 1], 1);
		if (*info != 0) {
		    itmp1 = ifail[i];
		    ifail[i] = ifail[j];
		    ifail[j] = itmp1;
		}
	    }
	}
    }
    return;
}
