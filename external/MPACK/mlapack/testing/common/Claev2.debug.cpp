/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Claev2.debug.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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

#include <mblas.h>
#include <mlapack.h>
#include <mpack_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define ITERATION 100

int errorflag = FALSE;
REAL_REF maxdiff = 0.0;

void errorcheck(REAL_REF rt1_ref, REAL_REF rt2_ref, REAL_REF cs1_ref, COMPLEX_REF sn1_ref, REAL rt1, REAL rt2, REAL cs1, COMPLEX sn1)
{
    REAL_REF diff;
    diff = abs(rt1_ref - rt1);
    printf("diff1    ");  printnum(diff); printf("\n");
    if (diff > EPSILON) {
	errorflag = TRUE;
	printf("Error1\n");
    }
    if (maxdiff < diff) maxdiff = diff;

    diff = abs(rt2_ref - rt2);
    printf("diff2    "); printnum(diff); printf("\n");
    if (diff > EPSILON) {
	errorflag = TRUE;
	printf("Error2\n");
    }
    if (maxdiff < diff) maxdiff = diff;

    diff = abs(cs1_ref - cs1);
    printf("diff3    "); printnum(diff); printf("\n");
    if (diff > EPSILON) {
	errorflag = TRUE;
	printf("Error3\n");
    }
    if (maxdiff < diff) maxdiff = diff;

    diff = abs(sn1_ref - sn1);
    printf("diff4    "); printnum(diff); printf("\n");
    if (diff > EPSILON) {
	errorflag = TRUE;
	printf("Error4\n");
    }
    if (maxdiff < diff) maxdiff = diff;
    printf("maxdiff  ");  printnum(maxdiff); printf("\n");
}

void Claev2_test()
{
    COMPLEX_REF a_ref, b_ref, c_ref, sn1_ref;
    REAL_REF rt1_ref, rt2_ref, cs1_ref;
    COMPLEX a, b, c, sn1;
    REAL rt1, rt2, cs1;

    int count = 100;
    while (count--) {
#if defined VERBOSE_TEST
	printf("Claev2: general random case\n");
#endif
	set_random_number(a_ref, a);
	set_random_number(b_ref, b);
	set_random_number(c_ref, c);
#if defined ___MPACK_BUILD_WITH_MPFR___
	zlaev2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#else
	Claev2(a_ref, b_ref, c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#endif
	Claev2(a, b, c, &rt1, &rt2, &cs1, &sn1);
	errorcheck(rt1_ref, rt2_ref, cs1_ref, sn1_ref, rt1, rt2, cs1, sn1);

/* checking adf = ab case (|a-c| = 2|b|) */
#if defined VERBOSE_TEST
	printf("Claev2: |a-c| = 2|b| case\n");
#endif
	set_random_number(a_ref, a);
	set_random_number(c_ref, c);
	b_ref = (a_ref - c_ref) / 2.0;
	b = (a - c) / 2.0;
#if defined ___MPACK_BUILD_WITH_MPFR___
	zlaev2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#else
	Claev2(a_ref, b_ref, c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#endif
	Claev2(a, b, c, &rt1, &rt2, &cs1, &sn1);
	errorcheck(rt1_ref, rt2_ref, cs1_ref, sn1_ref, rt1, rt2, cs1, sn1);

/* checking sm = 0 case */
#if defined VERBOSE_TEST
	printf("Claev2: sm = 0 case\n");
#endif
	set_random_number(a_ref, a);
	set_random_number(b_ref, b);
	c_ref = -a_ref;
	c = -a;
#if defined ___MPACK_BUILD_WITH_MPFR___
	zlaev2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#else
	Claev2(a_ref, b_ref, c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#endif
	Claev2(a, b, c, &rt1, &rt2, &cs1, &sn1);
	errorcheck(rt1_ref, rt2_ref, cs1_ref, sn1_ref, rt1, rt2, cs1, sn1);

/*zero eigenvalue case */
#if defined VERBOSE_TEST
	printf("Claev2: zero eigenvalue case\n");
#endif
	set_random_number(a_ref, a);
	b_ref = c_ref = 0.0;
	b = c = 0.0;
#if defined ___MPACK_BUILD_WITH_MPFR___
	zlaev2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#else
	Claev2(a_ref, b_ref, c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#endif
	Claev2(a, b, c, &rt1, &rt2, &cs1, &sn1);
	errorcheck(rt1_ref, rt2_ref, cs1_ref, sn1_ref, rt1, rt2, cs1, sn1);

/*zero matrix case */
#if defined VERBOSE_TEST
	printf("Claev2: zero matrix case\n");
#endif
	a_ref = b_ref = c_ref = 0.0;
	a = b = c = 0.0;
#if defined ___MPACK_BUILD_WITH_MPFR___
	zlaev2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#else
	Claev2(a_ref, b_ref, c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#endif
	Claev2(a, b, c, &rt1, &rt2, &cs1, &sn1);
	errorcheck(rt1_ref, rt2_ref, cs1_ref, sn1_ref, rt1, rt2, cs1, sn1);

/*Identity matrix case */
#if defined VERBOSE_TEST
	printf("Claev2: identity matrix case\n");
#endif
	a_ref = c_ref = 1.0;
	b_ref = 0.0;
	a = c = 1.0;
	b = 0.0;
#if defined ___MPACK_BUILD_WITH_MPFR___
	zlaev2_f77(&a_ref, &b_ref, &c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#else
	Claev2(a_ref, b_ref, c_ref, &rt1_ref, &rt2_ref, &cs1_ref, &sn1_ref);
#endif
	Claev2(a, b, c, &rt1, &rt2, &cs1, &sn1);
	errorcheck(rt1_ref, rt2_ref, cs1_ref, sn1_ref, rt1, rt2, cs1, sn1);
    }
}

int main(int argc, char *argv[])
{
    Claev2_test();
    printf("Claev2 test passed...\n");
    return (0);
}
