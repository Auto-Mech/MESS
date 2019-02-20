/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlatrs.cpp,v 1.15 2010/08/19 01:20:10 nakatamaho Exp $ 
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

#include <stdio.h> //for printf, this routine will be tesed and then printf will be removed. 

void Rlatrs(const char *uplo, const char *trans, const char *diag, const char *normin, INTEGER n, REAL * A, INTEGER lda, REAL * x,
	    REAL * scale, REAL * cnorm, INTEGER * info)
{
    LOGICAL notran, nounit, upper;
    INTEGER i, imax, j, jfirst, jinc, jlast;
    REAL bignum, grow, rec, smlnum, sumj, tjj, tjjs = 0.0;
    REAL tmax, tscal, uscal, xbnd, xj, xmax;
    REAL Zero = 0.0, Half = 0.5, One = 1.0;
    REAL mtemp1, mtemp2, mtemp3;

    *info = 0;
    upper = Mlsame(uplo, "U");
    notran = Mlsame(trans, "N");
    nounit = Mlsame(diag, "N");
//Test the input parameters.
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (!nounit && !Mlsame(diag, "U")) {
	*info = -3;
    } else if (!Mlsame(normin, "Y") && !Mlsame(normin, "N")) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    }
    if (*info != 0) {
	Mxerbla("Rlatrs", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
//Determine machine dependent parameters to control overflow.
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = One / smlnum;
    *scale = One;
    if (Mlsame(normin, "N")) {
//Compute the 1-norm of each column, not including the diagonal.
	if (upper) {
//A is upper triangular.
	    for (j = 1; j <= n; j++) {
		cnorm[j - 1] = Rasum(j - 1, &A[(j - 1) * lda], 1);
	    }
	} else {
//A is lower triangular.
	    for (j = 1; j <= n - 1; j++) {
		cnorm[j - 1] = Rasum(n - j, &A[j + (j - 1) * lda], 1);
	    }
	    cnorm[n - 1] = Zero;
	}
    }
//Scale the column norms by TSCAL if the maximum element in CNORM is
//greater than BIGNUM.
    imax = iRamax(n, cnorm, 1);
    tmax = cnorm[imax - 1];
    if (tmax <= bignum) {
	tscal = One;
    } else {
	tscal = One / (smlnum * tmax);
	Rscal(n, tscal, cnorm, 1);
    }
//Compute a bound on the computed solution vector to see if the
//Level 2 BLAS routine DTRSV can be used.
    j = iRamax(n, &x[0], 1);
    xmax = abs(x[j - 1]);
    xbnd = xmax;
    if (notran) {
// Compute the growth in A * x = b.
	if (upper) {
	    jfirst = n;
	    jlast = 1;
	    jinc = -1;
	} else {
	    jfirst = 1;
	    jlast = n;
	    jinc = 1;
	}
	if (tscal != One) {
	    grow = Zero;
	    goto L50;
	}
	if (nounit) {
//A is non-unit triangular.
//Compute GROW = 1/G(j) and XBND = 1/M(j).
//Initially, G(0) = max{x(i), i=1,...,n}.
	    grow = One / max(xbnd, smlnum);
	    xbnd = grow;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum)
		    goto L50;
//M(j) = G(j-1) / abs(A(j,j))
		tjj = abs(A[(j - 1) + (j - 1) * lda]);
		mtemp1 = xbnd;
		mtemp2 = min(One, tjj) * grow;
		xbnd = min(mtemp1, mtemp2);
		if (tjj + cnorm[j - 1] >= smlnum) {
//G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
		    grow = grow * tjj / (tjj + cnorm[j - 1]);
		} else {
//G(j) could overflow, set GROW to Zero
		    grow = Zero;
		}
	    }
	    grow = xbnd;
	} else {
//A is unit triangular.
//Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
	    mtemp1 = One, mtemp2 = One / max(xbnd, smlnum);
	    grow = min(mtemp1, mtemp2);
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum)
		    goto L50;
//G(j) = G(j-1)*( 1 + CNORM(j) )
		grow = grow * One / (cnorm[j - 1] + One);
	    }
	}
      L50:;
    } else {
//Compute the growth in A' * x = b.
	if (upper) {
	    jfirst = 1;
	    jlast = n;
	    jinc = 1;
	} else {
	    jfirst = n;
	    jlast = 1;
	    jinc = -1;
	}
	if (tscal != One) {
	    grow = Zero;
	    goto L80;
	}
	if (nounit) {
//A is non-unit triangular.
//Compute GROW = 1/G(j) and XBND = 1/M(j).
//Initially, M(0) = max{x(i), i=1,...,n}.
	    grow = One / max(xbnd, smlnum);
	    xbnd = grow;
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum) {
		    goto L80;
		}
//G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
		xj = cnorm[j - 1] + One;
		mtemp1 = grow, mtemp2 = xbnd / xj;
		grow = min(mtemp1, mtemp2);
//M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
		tjj = abs(A[(j - 1) + (j - 1) * lda]);
		if (xj > tjj) {
		    xbnd = xbnd * tjj / xj;
		}
	    }
	    grow = min(grow, xbnd);
	} else {
//A is unit triangular.
//Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
	    mtemp1 = One, mtemp2 = One / max(xbnd, smlnum);
	    grow = min(mtemp1, mtemp2);
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {
//Exit the loop if the growth factor is too small.
		if (grow <= smlnum) {
		    goto L80;
		}
//G(j) = ( 1 + CNORM(j) )*G(j-1)
		xj = cnorm[j - 1] + One;
		grow = grow / xj;
	    }
	}
      L80:;
    }
    if (grow * tscal > smlnum) {
//Use the Level 2 BLAS solve if the reciprocal of the bound on
//elements of X is not too small.
	Rtrsv(uplo, trans, diag, n, A, lda, x, 1);
    } else {
        printf("Rlatrs: not yet tested very much\n"); exit(-1);
//Use a Level 1 BLAS solve, scaling intermediate results.
	if (xmax > bignum) {
//Scale X so that its components are less than or equal to
//BIGNUM in absolute value.
	    *scale = bignum / xmax;
	    Rscal(n, *scale, &x[0], 1);
	    xmax = bignum;
	}
	if (notran) {
//Solve A * x = b
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {
//Compute x(j) = b(j) / A(j,j), scaling x if necessary.
		xj = abs(x[j - 1]);
		if (nounit) {
		    tjjs = A[(j - 1) + (j - 1) * lda] * tscal;
		} else {
		    tjjs = tscal;
		    if (tscal == One) {
			goto L100;
		    }
		}
		tjj = abs(tjjs);
		if (tjj > smlnum) {
//abs(A(j,j)) > SMLNUM:
		    if (tjj < One) {
			if (xj > tjj * bignum) {
//Scale x by 1/b(j).
			    rec = One / xj;
			    Rscal(n, rec, x, 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    x[j - 1] = x[j - 1] / tjjs;
		    xj = abs(x[j - 1]);
		} else if (tjj > Zero) {
//0 < abs(A(j,j)) <= SMLNUM:
		    if (xj > tjj * bignum) {
//Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
//to avoid overflow when dividing by A(j,j).
			rec = tjj * bignum / xj;
			if (cnorm[j - 1] > One) {
//Scale by 1/CNORM(j) to avoid overflow when
//multiplying x(j) times column j.
			    rec = rec / cnorm[j - 1];
			}
			Rscal(n, rec, &x[0], 1);
			*scale = (*scale) * rec;
			xmax = xmax * rec;
		    }
		    x[j - 1] = x[j - 1] / tjjs;
		    xj = abs(x[j - 1]);
		} else {
//A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//scale = 0, and compute a solution to A*x = Zero
		    for (i = 0; i < n; i++) {
			x[i] = Zero;
		    }
		    x[j - 1] = One;
		    xj = One;
		    *scale = Zero;
		    xmax = Zero;
		}
	      L100:
//Scale x if necessary to avoid overflow when adding a
//multiple of column j of A.
		if (xj > One) {
		    rec = One / xj;
		    if (cnorm[j - 1] > (bignum - xmax) * rec) {
//Scale x by 1/(2*abs(x(j))).
			rec = rec * Half;
			Rscal(n, rec, &x[0], 1);
			*scale = *scale * rec;
		    }
		} else if (xj * cnorm[j - 1] > bignum - xmax) {
//Scale x by 1/Two
		    Rscal(n, Half, &x[0], 1);
		    *scale = *scale * Half;
		}
		if (upper) {
		    if (j > 1) {
//Compute the update
//       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
			Raxpy(j - 1, -x[j - 1] * tscal, &A[(j - 1) * lda], 1, x, 1);
			i = iRamax(j - 1, &x[0], 1);
			xmax = abs(x[i-1]);
		    }
		} else {
		    if (j < n) {
//Compute the update
//x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
			Raxpy(n - j, -x[j - 1] * tscal, &A[j + (j - 1) * lda], 1, &x[j], 1);
			i = j + iRamax(n - j, &x[j], 1);
			xmax = abs(x[i-1]);
		    }
		}
	    }
	} else {
//Solve A' * x = b
	    for (j = jfirst; jinc < 0 ? j >= jlast : j <= jlast; j += jinc) {
//Compute x(j) = b(j) - sum A(k,j)*x(k).
//             k<>j
		xj = abs(x[j - 1]);
		uscal = tscal;
		rec = One / max(xmax, One);
		if (cnorm[j - 1] > (bignum - xj) * rec) {
//If x(j) could overflow, scale x by 1/(2*XMAX).
		    rec = rec * Half;
		    if (nounit) {
			tjjs = A[(j - 1) + (j - 1) * lda] * tscal;
		    } else {
			tjjs = tscal;
		    }
		    tjj = abs(tjjs);
		    if (tjj > One) {
//Divide by A(j,j) when scaling x if A(j,j) > One
			mtemp1 = One, mtemp2 = rec * tjj;
			rec = min(mtemp1, mtemp2);
			uscal = uscal / tjjs;
		    }
		    if (rec < One) {
			Rscal(n, rec, &x[0], 1);
			*scale = *scale * rec;
			xmax = xmax * rec;
		    }
		}
		sumj = Zero;
		if (uscal == One) {
//If the scaling needed for A in the dot product is 1,
//call DDOT to perform the dot product.
		    if (upper) {
			sumj = Rdot(j - 1, &A[(j - 1) * lda], 1, &x[0], 1);
		    } else if (j < n) {
			sumj = Rdot(n - j, &A[j + (j - 1) * lda], 1, &x[j], 1);
		    }
		} else {
//Otherwise, use in-line code for the dot product.
		    if (upper) {
			for (i = 1; i <= j - 1; i++) {
			    sumj = sumj + A[(i - 1) + (j - 1) * lda] * uscal * x[i - 1];
			}
		    } else if (j < n) {
			for (i = j + 1; i <= n; i++) {
			    sumj = sumj + A[(i - 1) + (j - 1) * lda] * uscal * x[i - 1];
			}
		    }
		}
		if (uscal == tscal) {
//Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
//was not used to scale the dotproduct. */
		    x[j - 1] = x[j - 1] - sumj;
		    xj = abs(x[j - 1]);
		    if (nounit) {
			tjjs = A[(j - 1) + (j - 1) * lda] * tscal;
		    } else {
			tjjs = tscal;
			if (tscal == One) {
			    goto L150;
			}
		    }
//Compute x(j) = x(j) / A(j,j), scaling if necessary.
		    tjj = abs(tjjs);
		    if (tjj > smlnum) {
//abs(A(j,j)) > SMLNUM:
			if (tjj < One) {
			    if (xj > tjj * bignum) {
//Scale X by 1/abs(x(j)).
				rec = One / xj;
				Rscal(n, rec, &x[0], 1);
				*scale = *scale * rec;
				xmax = xmax * rec;
			    }
			}
			x[j - 1] = x[j - 1] / tjjs;
		    } else if (tjj > Zero) {
//0 < abs(A(j,j)) <= SMLNUM:
			if (xj > tjj * bignum) {
//Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
			    rec = tjj * bignum / xj;
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
			x[j - 1] = x[j - 1] / tjjs;
		    } else {
//A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
//scale = 0, and compute a solution to A'*x = Zero
			for (i = 0; i < n; i++) {
			    x[i] = Zero;
			}
			x[j - 1] = One;
			*scale = Zero;
			xmax = Zero;
		    }
		  L150:
		    ;
		} else {
//Compute x(j) := x(j) / A(j,j)  - sumj if the dot
//product has already been divided by 1/A(j,j).
		    x[j - 1] = x[j - 1] / tjjs - sumj;
		}
		mtemp2 = xmax, mtemp3 = abs(x[j-1]);
		xmax = max(mtemp2, mtemp3);
	    }
	}
	*scale = *scale / tscal;
    }
//Scale the column norms by 1/TSCAL for return.
    if (tscal != One) {
	Rscal(n, One / tscal, cnorm, 1);
    }
    return;
}
