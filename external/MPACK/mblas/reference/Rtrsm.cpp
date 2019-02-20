/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rtrsm.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
 *
 * $Id: Rtrsm.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $

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

/*
Based on http://www.netlib.org/blas/dtrsm.f
Rtrsm solves one of the matrix equations
 op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
where alpha is a scalar, X and B are m by n matrices, A is a unit, or
non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   op( A ) = A   or   op( A ) = A'.
The matrix X is overwritten on B.
*/

#include <mblas.h>

void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER m, INTEGER n, REAL alpha, REAL * A,
	   INTEGER lda, REAL * B, INTEGER ldb)
{
    INTEGER i, info, j, k, lside, nrowa, nounit, upper;
    REAL Zero = 0.0, One = 1.0;
    REAL temp;

//test the input parameters.
    lside = Mlsame(side, "L");
    if (lside)
	nrowa = m;
    else
	nrowa = n;

    nounit = Mlsame(diag, "N");
    upper = Mlsame(uplo, "U");

    info = 0;
    if ((!lside) && (!Mlsame(side, "R")))
	info = 1;
    else if ((!upper) && (!Mlsame(uplo, "L")))
	info = 2;
    else if ((!Mlsame(transa, "N")) && (!Mlsame(transa, "T")) && (!Mlsame(transa, "C")))
	info = 3;
    else if ((!Mlsame(diag, "U")) && (!Mlsame(diag, "N")))
	info = 4;
    else if (m < 0)
	info = 5;
    else if (n < 0)
	info = 6;
    else if (lda < max((INTEGER) 1, nrowa))
	info = 9;
    else if (ldb < max((INTEGER) 1, m))
	info = 11;
    if (info != 0) {
	Mxerbla("Rtrsm ", info);
	return;
    }
//quick return if possible.
    if (m == 0 || n == 0)
	return;

//and when alpha==zero.
    if (alpha == Zero) {
	for (j = 0; j < n; j++) {
	    for (i = 0; i < m; i++) {
		B[i + j * ldb] = Zero;
	    }
	}
	return;
    }
//start the operations.
    if (lside) {
	if (Mlsame(transa, "N")) {
//Form B := alpha*inv(A)*B.
	    if (upper) {
		for (j = 0; j < n; j++) {
		    if (alpha != One) {
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = alpha * B[i + j * ldb];
			}
		    }
		    for (k = m - 1; k >= 0; k--) {
			if (B[k + j * ldb] != Zero) {
			    if (nounit)
				B[k + j * ldb] = B[k + j * ldb] / A[k + k * lda];
			    for (i = 0; i < k; i++) {
				B[i + j * ldb] = B[i + j * ldb] - B[k + j * ldb] * A[i + k * lda];
			    }
			}
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    if (alpha != One) {
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = alpha * B[i + j * ldb];
			}
		    }
		    for (k = 0; k < m; k++) {
			if (B[k + j * ldb] != Zero) {
			    if (nounit)
				B[k + j * ldb] = B[k + j * ldb] / A[k + k * lda];
			    for (i = k + 1; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] - B[k + j * ldb] * A[i + k * lda];
			    }
			}
		    }
		}
	    }
	} else {
	    //Form B := alpha*inv(A')*B.
	    if (upper) {
		for (j = 0; j < n; j++) {
		    for (i = 0; i < m; i++) {
			temp = alpha * B[i + j * ldb];
			for (k = 0; k < i; k++) {
			    temp = temp - A[k + i * lda] * B[k + j * ldb];
			}
			if (nounit)
			    temp = temp / A[i + i * lda];
			B[i + j * ldb] = temp;
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (i = m - 1; i >= 0; i--) {
			temp = alpha * B[i + j * ldb];
			for (k = i + 1; k < m; k++) {
			    temp = temp - A[k + i * lda] * B[k + j * ldb];
			}
			if (nounit)
			    temp = temp / A[i + i * lda];
			B[i + j * ldb] = temp;
		    }
		}
	    }
	}
    } else {
	if (Mlsame(transa, "N")) {
//Form B := alpha*B*inv(A).
	    if (upper) {
		for (j = 0; j < n; j++) {
		    if (alpha != One) {
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = alpha * B[i + j * ldb];
			}
		    }
		    for (k = 0; k < j; k++) {
			if (A[k + j * lda] != Zero) {
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] - A[k + j * lda] * B[i + k * ldb];
			    }
			}
		    }
		    if (nounit) {
			temp = One / A[j + j * lda];
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = temp * B[i + j * ldb];
			}
		    }
		}
	    } else {
		for (j = n - 1; j >= 0; j--) {
		    if (alpha != One) {
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = alpha * B[i + j * ldb];
			}
		    }
		    for (k = j + 1; k < n; k++) {
			if (A[k + j * lda] != Zero) {
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] - A[k + j * lda] * B[i + k * ldb];
			    }
			}
		    }
		    if (nounit) {
			temp = One / A[j + j * lda];
			for (i = 0; i < m; i++) {
			    B[i + j * ldb] = temp * B[i + j * ldb];
			}
		    }
		}
	    }
	} else {
//Form  B := alpha*B*inv(A').
	    if (upper) {
		for (k = n - 1; k >= 0; k--) {
		    if (nounit) {
			temp = One / A[k + k * lda];
			for (i = 0; i < m; i++) {
			    B[i + k * ldb] = temp * B[i + k * ldb];
			}
		    }
		    for (j = 0; j < k; j++) {
			if (A[j + k * lda] != Zero) {
			    temp = A[j + k * lda];
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] - temp * B[i + k * ldb];
			    }
			}
		    }
		    if (alpha != One) {
			for (i = 0; i < m; i++) {
			    B[i + k * ldb] = alpha * B[i + k * ldb];
			}
		    }
		}
	    } else {
		for (k = 0; k < n; k++) {
		    if (nounit) {
			temp = One / A[k + k * lda];
			for (i = 0; i < m; i++) {
			    B[i + k * ldb] = temp * B[i + k * ldb];
			}
		    }
		    for (j = k + 1; j < n; j++) {
			if (A[j + k * lda] != Zero) {
			    temp = A[j + k * lda];
			    for (i = 0; i < m; i++) {
				B[i + j * ldb] = B[i + j * ldb] - temp * B[i + k * ldb];
			    }
			}
		    }
		    if (alpha != One) {
			for (i = 0; i < m; i++) {
			    B[i + k * ldb] = alpha * B[i + k * ldb];
			}
		    }
		}
	    }
	}
    }
    return;
}
