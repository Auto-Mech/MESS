/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgees.cpp,v 1.15 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rgees(const char *jobvs, const char *sort, LFP select, INTEGER n, REAL * A,
      INTEGER lda, INTEGER * sdim, REAL * wr, REAL * wi, REAL * vs, INTEGER ldvs, REAL * work, INTEGER lwork, LOGICAL * bwork, INTEGER * info)
{
    INTEGER i;
    REAL s;
    INTEGER i1, i2, ip, ihi, ilo;
    REAL dum[1], eps, sep;
    INTEGER ibal;
    REAL anrm;
    INTEGER idum[1], ierr, itau, iwrk, inxt, icond, ieval;
    INTEGER cursl;
    INTEGER lst2sl, scalea;
    REAL cscale = 0.0;
    REAL bignum;
    INTEGER lastsl;
    INTEGER minwrk, maxwrk;
    REAL smlnum;
    INTEGER hswork;
    INTEGER wantst, lquery, wantvs;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    wantvs = Mlsame(jobvs, "V");
    wantst = Mlsame(sort, "S");
    if (!wantvs && !Mlsame(jobvs, "N")) {
	*info = -1;
    } else if (!wantst && !Mlsame(sort, "N")) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldvs < 1 || (wantvs && ldvs < n)) {
	*info = -11;
    }
//Compute workspace
//(Note: Comments in the code beginning "Workspace:" describe the
//minimal amount of workspace needed at that point in the code,
//as well as the preferred amount for good performance.
//NB refers to the optimal block size for the immediately
//following subroutine, as returned by ILAENV.
//HSWORK refers to the workspace preferred by DHSEQR, as
//calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
//the worst case.)
    if (*info == 0) {
	if (n == 0) {
	    minwrk = 0;
	    maxwrk = 0;
	} else {
	    maxwrk = (n * 2) + n * iMlaenv(1, "Rgehrd", " ", n, 1, n, 0);
	    minwrk = n * 3;
	    Rhseqr("S", jobvs, n, 1, n, A, lda, wr, wi, vs, ldvs, work, -1, &ieval);
	    hswork = (INTEGER) (cast2double(work[1]));
	    if (!wantvs) {
		maxwrk = max(maxwrk, n + hswork);
	    } else {
		maxwrk = max(maxwrk, (n * 2) + (n - 1) * iMlaenv(1, "Rorghr", " ", n, 1, n, -1));
		maxwrk = max(maxwrk, n + hswork);
	    }
	}
	work[1] = (REAL) double (maxwrk);
	if (lwork < minwrk && !lquery) {
	    *info = -13;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgees ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	*sdim = 0;
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = Rlamch("S");
    bignum = One / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Rlange("M", n, n, A, lda, dum);
    scalea = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	scalea = MTRUE;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = MTRUE;
	cscale = bignum;
    }
    if (scalea) {
	Rlascl("G", 0, 0, anrm, cscale, n, n, &A[0], lda, &ierr);
    }
// Permute the matrix to make it more nearly triangular
//(Workspace: need N)
    ibal = 0;
    Rgebal("P", n, A, lda, &ilo, &ihi, &work[ibal], &ierr);
//Reduce to upper Hessenberg form
//(Workspace: need 3*N, prefer 2*N+N*NB)
    itau = n + ibal;
    iwrk = n + itau;
    Rgehrd(n, ilo, ihi, A, lda, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    if (wantvs) {
//Copy Householder vectors to VS
	Rlacpy("L", n, n, A, lda, vs, ldvs);
//Generate orthogonal matrix in VS
//(Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
	Rorghr(n, ilo, ihi, vs, ldvs, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    }
    *sdim = 0;
//Perform QR iteration, accumulating Schur vectors in VS if desired
//(Workspace: need N+1, prefer N+HSWORK (see comments) )
    iwrk = itau;
    Rhseqr("S", jobvs, n, ilo, ihi, A, lda, wr, wi, vs, ldvs, &work[iwrk], lwork - iwrk + 1, &ieval);
    if (ieval > 0) {
	*info = ieval;
    }
//Sort eigenvalues if desired
    if (wantst && *info == 0) {
	if (scalea) {
	    Rlascl("G", 0, 0, cscale, anrm, n, 1, &wr[1], n, &ierr);
	    Rlascl("G", 0, 0, cscale, anrm, n, 1, &wi[1], n, &ierr);
	}
	for (i = 0; i < n; i++) {
	    bwork[i] = (*select) (&wr[i], &wi[i]);
	}
//Reorder eigenvalues and transform Schur vectors
//(Workspace: none needed)
	Rtrsen("N", jobvs, bwork, n, A, lda, vs, ldvs, wr, wi, *sdim, &s, &sep, &work[iwrk], lwork - iwrk + 1, idum, 1, &icond);
	if (icond > 0) {
	    *info = n + icond;
	}
    }
    if (wantvs) {
//Undo balancing
//(Workspace: need N)
	Rgebak("P", "R", n, ilo, ihi, &work[ibal], n, vs, ldvs, &ierr);
    }
    if (scalea) {
//Undo scaling for the Schur form of A
	Rlascl("H", 0, 0, cscale, anrm, n, n, A, lda, &ierr);
	Rcopy(n, A, lda + 1, wr, 1);
	if (cscale == smlnum) {
//If scaling back towards underflow, adjust WI if an
//offdiagonal element of a 2-by-2 block in the Schur form
//underflows.
	    if (ieval > 0) {
		i1 = ieval + 1;
		i2 = ihi - 1;
		Rlascl("G", 0, 0, cscale, anrm, ilo - 1, 1, &wi[1], max(ilo - 1, (INTEGER) 1), &ierr);
	    } else if (wantst) {
		i1 = 1;
		i2 = n - 1;
	    } else {
		i1 = ilo;
		i2 = ihi - 1;
	    }
	    inxt = i1 - 1;
	    for (i = i2; i <= i1; i++) {
		if (i < inxt) {
		    goto L20;
		}
		if (wi[i] == Zero) {
		    inxt = i + 1;
		} else {
		    if (A[i + 1 + i * lda] == Zero) {
			wi[i] = Zero;
			wi[i + 1] = Zero;
		    } else if (A[i + 1 + i * lda] != Zero && A[i + (i + 1) * lda] == Zero) {
			wi[i] = Zero;
			wi[i + 1] = Zero;
			if (i > 1) {
			    Rswap(i - 1, &A[i * lda], 1, &A[(i + 1) * lda], 1);
			}
			if (n > i + 1) {
			    Rswap(n - i - 1, &A[i + (i + 2) * lda], lda, &A[i + 1 + (i + 2) * lda], lda);
			}
			if (wantvs) {
			    Rswap(n, &vs[i * ldvs + 1], 1, &vs[(i + 1) * ldvs + 1], 1);
			}
			A[i + (i + 1) * lda] = A[i + 1 + i * lda];
			A[i + 1 + i * lda] = Zero;
		    }
		    inxt = i + 2;
		}
	      L20:
		;
	    }
	}
//Undo scaling for the imaginary part of the eigenvalues
	Rlascl("G", 0, 0, cscale, anrm, n - ieval, 1, &wi[ieval + 1], max(n - ieval, (INTEGER) 1), &ierr);
    }
    if (wantst && *info == 0) {
//Check if reordering successful
	lastsl = MTRUE;
	lst2sl = MTRUE;
	*sdim = 0;
	ip = 0;
	for (i = 0; i < n; i++) {
	    cursl = (*select) (&wr[i], &wi[i]);
	    if (wi[i] == Zero) {
		if (cursl) {
		    ++(*sdim);
		}
		ip = 0;
		if (cursl && !lastsl) {
		    *info = n + 2;
		}
	    } else {
		if (ip == 1) {
//Last eigenvalue of conjugate pair
		    cursl = cursl || lastsl;
		    lastsl = cursl;
		    if (cursl) {
			*sdim += 2;
		    }
		    ip = -1;
		    if (cursl && !lst2sl) {
			*info = n + 2;
		    }
		} else {
//First eigenvalue of conjugate pair
		    ip = 1;
		}
	    }
	    lst2sl = lastsl;
	    lastsl = cursl;
	}
    }
    work[1] = (REAL) double (maxwrk);
    return;
}
