/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rsyrk.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Rsyrk.cpp,v 1.6 2010/08/07 05:50:10 nakatamaho Exp $

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
Rsyrk performs one of the symmetric rank k operations
 C := alpha*A*A' + beta*C,
or
 C := alpha*A'*A + beta*C,
where  alpha and beta  are scalars, C is an n by n symmetric matrix
and A is an n by k matrix in the first case and a k by n matrix
in the second case.
*/

#include <mblas.h>

void Rsyrk(const char *uplo, const char *trans, INTEGER n, INTEGER k, REAL alpha, REAL * A, INTEGER lda, REAL beta, REAL * C,
	   INTEGER ldc)
{
    INTEGER i, j, l, nrowa, upper, info;
    REAL Zero = 0.0, One = 1.0;
    REAL temp;

//Test the input parameters.
    if (Mlsame(trans, "N"))
	nrowa = n;
    else
	nrowa = k;
    upper = Mlsame(uplo, "U");

    info = 0;
    if ((!upper) && (!Mlsame(uplo, "L")))
	info = 1;
    else if ((!Mlsame(trans, "N")) && (!Mlsame(trans, "T")) && (!Mlsame(trans, "C")))
	info = 2;
    else if (n < 0)
	info = 3;
    else if (k < 0)
	info = 4;
    else if (lda < max((INTEGER) 1, nrowa))
	info = 7;
    else if (ldc < max((INTEGER) 1, n))
	info = 10;
    if (info != 0) {
	Mxerbla("Rsyrk ", info);
	return;
    }
//quick return if possible.
    if ((n == 0) || (((alpha == Zero) || (k == 0)) && (beta == One)))
	return;

//and when alpha==Zero.
    if (alpha == Zero) {
	if (upper) {
	    if (beta == Zero) {
		for (j = 0; j < n; j++) {
		    for (i = 0; i <= j; i++) {
			C[i + j * ldc] = Zero;
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (i = 0; i <= j; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
	    }
	} else {
	    if (beta == Zero) {
		for (j = 0; j < n; j++) {
		    for (i = j; i < n; i++) {
			C[i + j * ldc] = Zero;
		    }
		}
	    } else {
		for (j = 0; j < n; j++) {
		    for (i = j; i < n; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
	    }
	}
	return;
    }
//start the operations.
    if (Mlsame(trans, "N")) {
//Form C := alpha*A*A' + beta*C.
	if (upper) {
	    for (j = 0; j < n; j++) {
		if (beta == Zero) {
		    for (i = 0; i <= j; i++) {
			C[i + j * ldc] = Zero;
		    }
		} else if (beta != One) {
		    for (i = 0; i <= j; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
		for (l = 0; l < k; l++) {
		    if (A[j + l * lda] != Zero) {
			temp = alpha * A[j + l * lda];
			for (i = 0; i <= j; i++) {
			    C[i + j * ldc] = C[i + j * ldc] + temp * A[i + l * lda];
			}
		    }
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		if (beta == Zero) {
		    for (i = j; i < n; i++) {
			C[i + j * ldc] = Zero;
		    }
		} else if (beta != One) {
		    for (i = j; i < n; i++) {
			C[i + j * ldc] = beta * C[i + j * ldc];
		    }
		}
		for (l = 0; l < k; l++) {
		    if (A[j + l * lda] != Zero) {
			temp = alpha * A[j + l * lda];
			for (i = j; i < n; i++) {
			    C[i + j * ldc] = C[i + j * ldc] + temp * A[i + l * lda];
			}
		    }
		}
	    }
	}
    } else {
//Form C := alpha*A'*A + beta*C.
	if (upper) {
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= j; i++) {
		    temp = Zero;
		    for (l = 0; l < k; l++) {
			temp = temp + A[l + i * lda] * A[l + j * lda];
		    }
		    if (beta == Zero) {
			C[i + j * ldc] = alpha * temp;
		    } else {
			C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
		    }
		}
	    }
	} else {
	    for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
		    temp = Zero;
		    for (l = 0; l < k; l++) {
			temp = temp + A[l + i * lda] * A[l + j * lda];
		    }
		    if (beta == Zero)
			C[i + j * ldc] = alpha * temp;
		    else
			C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
		}
	    }
	}
    }
    return;
}
