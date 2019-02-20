/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rpotrf.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rpotrf(const char *uplo, INTEGER n, REAL * A, INTEGER lda, INTEGER * info)
{
    INTEGER upper;
    INTEGER j, jb, nb;
    REAL One = 1.0;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rpotrf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Determine the block size for this environment.
    nb = iMlaenv(1, "Rpotrf", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
//Use unblocked code.
	Rpotf2(uplo, n, A, lda, info);
    } else {
//Use blocked code.
	if (upper) {
//Compute the Cholesky factorization A = U'*U.
	    for (j = 1; j <= n; j = j + nb) {
//Update and factorize the current diagonal block and test
//for non-positive-definiteness.
		jb = min(nb, n - j + 1);
		Rsyrk("Upper", "Transpose", jb, j - 1, -One, &A[0 + (j - 1) * lda], lda, One, &A[(j - 1) + (j - 1) * lda], lda);
		Rpotf2("Upper", jb, &A[(j - 1) + (j - 1) * lda], lda, info);
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= n) {
//Compute the current block row.
		    Rgemm("Transpose", "No transpose", jb, n - j - jb + 1,
			  j - 1, -One, &A[0 + (j - 1) * lda], lda, &A[0 + (j + jb - 1) * lda], lda, One, &A[(j - 1) + (j + jb - 1) * lda], lda);
		    Rtrsm("Left", "Upper", "Transpose", "Non-unit", jb, n - j - jb + 1, One, &A[(j - 1) + (j - 1) * lda], lda, &A[(j - 1) + (j + jb - 1) * lda], lda);
		}
	    }
	} else {
//Compute the Cholesky factorization A = L*L'.
	    for (j = 1; j <= n; j = j + nb) {
//Update and factorize the current diagonal block and test
//for non-positive-definiteness.
		jb = min(nb, n - j + 1);
		Rsyrk("Lower", "No transpose", jb, j - 1, -One, &A[(j - 1) + 0 * lda], lda, One, &A[(j - 1) + (j - 1) * lda], lda);
		Rpotf2("Lower", jb, &A[(j - 1) + (j - 1) * lda], lda, info);
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= n) {
//Compute the current block column.
		    Rgemm("No transpose", "Transpose", n - j - jb + 1, jb,
			  j - 1, -One, &A[(j + jb - 1) + 0 * lda], lda, &A[(j - 1) + 0 * lda], lda, One, &A[(j + jb - 1) + (j - 1) * lda], lda);
		    Rtrsm("Right", "Lower", "Transpose", "Non-unit", n - j - jb + 1, jb, One, &A[(j - 1) + (j - 1) * lda], lda, &A[(j + jb - 1) + (j - 1) * lda], lda);
		}
	    }
	}
    }
    return;
  L30:
    *info = *info + j - 1;
    return;
}
