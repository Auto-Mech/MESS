/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Mutils.debug.cpp,v 1.13 2010/08/07 05:50:10 nakatamaho Exp $
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
 * $Id: Mutils.debug.cpp,v 1.13 2010/08/07 05:50:10 nakatamaho Exp $

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
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#include <math.h>

#define ITERATION 100

#define MAXSIZE 10
#define MAXCOND 16

int Mutils_test_highcond()
{
  for (int n = 2 ; n < MAXSIZE;  n++ ) {
    REAL_REF *A_ref = new REAL_REF [n * n];
    REAL_REF *w_ref = new REAL_REF[veclen(n, 1)];
    REAL_REF *work_ref = new REAL_REF[veclen(3*n-1, 1)];
    REAL *A = new REAL [n * n];
    REAL *w = new REAL [veclen(n, 1)];
    REAL *work = new REAL[veclen(3*n-1, 1)];

    for (int c = 1 ; c < MAXCOND;  c++ ) {
      set_random_symmmat_cond(A_ref, A, n, n, c); 
      printf("approx cond: %d\n",c);
      printf("A="); printmat(n, n, A, n); printf("\n");printf("\n");
    }

    delete []work;
    delete []w;
    delete []A;
    delete []work_ref;
    delete []w_ref;
    delete []A_ref;
  }
}

//ln2 test
int Mutils_test_log2()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    printf("approx ln2 test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	a_ref = abs(a_ref);
        a = abs(a);

	b_ref = log2(a_ref);
	b = log2(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "ln2_ref= "; printnum(b_ref); cout << endl;
	cout << "ln2a=    "; printnum(b); cout << endl;
//	cout << "residue=ln2a-ln2_ref\n";
//      cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//       cout << endl;
#endif
        diff = abs( b_ref - b);
        printf("diff     "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error in log2\n");
	    exit(1);
	}
    }
    return errorflag;
}

int Mutils_test_log10()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    printf("approx log10 test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	a_ref = abs(a_ref);
	b_ref = log10(a_ref);
        a = abs(a);
	b = log10(a);
#if defined VERBOSE_TEST
	cout << "a_ref=      "; printnum(a_ref); cout << endl;
	cout << "a=          "; printnum(a); cout << endl;
	cout << "log10a_ref= "; printnum(b_ref);cout << endl;
	cout << "log10a=     "; printnum(b); cout << endl;
//	cout << "residue=log10a-log10a_ref" << endl;
//	cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//      cout << endl;
#endif
        diff = abs (b_ref - b);
        printf("diff        "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error in log10\n");
	    exit(1);
	}
    }
    return errorflag;
}

//ln test
int Mutils_test_log()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    printf("approx log test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	a_ref = abs(a_ref);
	a = abs(a);

	b_ref = log(a_ref);
	b = log(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "lna_ref= "; printnum(b_ref); cout << endl;
	cout << "lna=     "; printnum(b); cout << endl;
//	cout << "residue=lna-lna_ref" << endl;
//	cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//	cout << endl;
#endif
        diff = abs (b_ref - b);
        printf("diff     "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error log\n");
	    exit(1);
	}
    }
    return errorflag;
}

int Mutils_test_sign()
{
    int errorflag = TRUE;
    REAL mtemp;

    printf("sign test\n");
    mtemp = sign((REAL)1.0, (REAL)-1.0);
    if (mtemp != -1.0) {
	printf("sign error\n");
	errorflag = FALSE;
    }
    mtemp = sign((REAL)-1.0, (REAL)1.0);
    if (mtemp != 1.0) {
	printf("sign error\n");
	errorflag = FALSE;
    }
    mtemp = sign((REAL)-1.0, (REAL)0.0);
    if (mtemp != 1.0) {
	printf("sign error\n");
	errorflag = FALSE;
    }

    mtemp = sign((REAL)1.0, (REAL)0.0);
    if (mtemp != 1.0) {
	printf("sign error\n");
	errorflag = FALSE;
    }

    mtemp = sign((REAL)0.0, (REAL)0.0);
    if (mtemp != 0.0) {
	printf("sign error\n");
	errorflag = FALSE;
    }

    if (errorflag == FALSE) {
	printf("Error in sign\n");
	exit(1);
    }
    return errorflag;
}

int Mutils_test_pow()
{
    int errorflag = TRUE;
    REAL_REF x_ref, y_ref, z_ref, diff;
    REAL x, y, z;

    printf("approx pow test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(x_ref, x);
	set_random_number(y_ref, y);
	x_ref = abs(x_ref);
	y_ref = abs(y_ref);
#if defined ___MPACK_BUILD_WITH_MPFR___
        z_ref = mpfr::pow(x_ref, y_ref); //somehow cannot avoid
#else
        z_ref = pow(x_ref, y_ref);
#endif
	x = abs(x);
	y = abs(y);
#if defined ___MPACK_BUILD_WITH_DOUBLE___
        z = std::pow(x, y);  //somehow cannot avoid
#else
	z = pow(x, y);
#endif

#if defined VERBOSE_TEST
	cout << "x_ref=   "; printnum(x_ref); cout << endl;
	cout << "y_ref=   "; printnum(y_ref); cout << endl;
	cout << "z_ref=   "; printnum(z_ref); cout << endl;
	cout << "x=       "; printnum(x); cout << endl;
	cout << "y=       "; printnum(y); cout << endl;
	cout << "z=       "; printnum(z); cout << endl;
//	cout << "residue=z-z_ref "; cout << endl;
//	cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//	cout << endl;
#endif
        diff = abs (z_ref - z);
        printf("diff     "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error in pow\n");
	    exit(1);
	}
    }
    return errorflag;
}

//sin test
int Mutils_test_sin()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    printf("approx sin test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	b_ref = sin(a_ref);
	b = sin(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "sin_ref= "; printnum(b_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "sin=     "; printnum(b); cout << endl;
	cout << endl;
//	cout << "residue=sre-s" << endl;
//	cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//	cout << endl;
#endif
        diff = abs( b_ref - b);
        printf("diff     "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error in sin\n");
	    exit(1);
	}
    }
    return errorflag;
}

//cos test
int Mutils_test_cos()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    printf("approx cos test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	b_ref = cos(a_ref);
	b = cos(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "cos_ref= "; printnum(b_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "cos=     "; printnum(b); cout << endl;
//	cout << "residue=s-sd" << endl;
//	cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//	cout << endl;
#endif
        diff = abs( b_ref - b);
        printf("diff     "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error in cos\n");
	    exit(1);
	}
    }
    return errorflag;
}

//exp test
int Mutils_test_exp()
{
    int errorflag = TRUE;
    REAL_REF a_ref, b_ref, diff;
    REAL a, b;

    printf("approx exp test\n");
    for (int i = 0; i < ITERATION; i++) {
	set_random_number(a_ref, a);
	b_ref = exp(a_ref);
	b = exp(a);
#if defined VERBOSE_TEST
	cout << "a_ref=   "; printnum(a_ref); cout << endl;
	cout << "exp_ref= "; printnum(b_ref); cout << endl;
	cout << "a=       "; printnum(a); cout << endl;
	cout << "exp=     "; printnum(b); cout << endl;
//	cout << "residue=s-sd" << endl;
//	cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//	cout << endl;
#endif
        diff = abs( b_ref - b);
        printf("diff     "); printnum(diff); printf("\n\n");
	if (diff > EPSILON) {
	    errorflag = TRUE;
	    printf("Error in exp\n");
	    exit(1);
	}
    }
    return errorflag;
}

int Mutils_test_pi()
{
    int errorflag = TRUE;
    REAL_REF p_ref, diff, dummy_ref = 0.0;
    REAL p, dummy = 0.0;

    printf("approx constant pi test\n");
    p_ref = pi(dummy_ref);
    p = pi(dummy);
#if defined VERBOSE_TEST
    cout << "p=       "; printnum(p); cout << endl;
    cout << "p_ref=   "; printnum(p_ref); cout << endl;
    cout << "residue=p-p_ref" << endl;
//    cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
//    cout << endl;
#endif
    diff = abs( p_ref - p);
    printf("diff     "); printnum(diff); printf("\n\n");
    if (diff > EPSILON) {
	errorflag = TRUE;
	printf("Error in pi\n");
	exit(1);
    }
    return errorflag;
}

void Mutils_test()
{
    Mutils_test_highcond();
    Mutils_test_log2();
    Mutils_test_log10();
    Mutils_test_log();
    Mutils_test_sign();
    Mutils_test_pow();
    Mutils_test_sin();
    Mutils_test_cos();
    Mutils_test_exp();
    Mutils_test_pi();
}

int main(int argc, char *argv[])
{
    Mutils_test();
    printf("Mutils test passed succesfully\n");
    return (0);
}
