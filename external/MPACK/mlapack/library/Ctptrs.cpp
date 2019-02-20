/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctptrs.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctptrs(const char *uplo, const char *trans, const char *diag, INTEGER n, INTEGER nrhs, COMPLEX * ap, COMPLEX * B, INTEGER ldb, INTEGER * info)
{
    INTEGER j, jc;
    INTEGER upper;
    INTEGER nounit;
    REAL Zero = 0.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    nounit = Mlsame(diag, "N");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "T")
	       && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (!nounit && !Mlsame(diag, "U")) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (nrhs < 0) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    }
    if (*info != 0) {

	Mxerbla("CTPTRS", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Check for singularity.
    if (nounit) {
	if (upper) {
	    jc = 1;
	    for (*info = 1; *info <= n; ++(*info)) {
		if (ap[jc + *info - 1] == Zero) {
		    return;
		}
		jc = jc + *info;
	    }
	} else {
	    jc = 1;
	    for (*info = 1; *info <= n; ++(*info)) {
		if (ap[jc] == Zero) {
		    return;
		}
		jc = jc + n - *info + 1;
	    }
	}
    }
    *info = 0;
//Solve  A * x = b,  A**T * x = b,  or  A**H * x = b.
    for (j = 0; j < nrhs; j++) {
	Ctpsv(uplo, trans, diag, n, ap, &B[j * ldb + 1], 1);
    }
    return;
}
