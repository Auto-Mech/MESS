/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mutils_mpfr.h,v 1.9 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_MPFR_H_
#define _MUTILS_MPFR_H_

mpreal pi(mpreal dummy);
mpreal sign(mpreal a, mpreal b);
double cast2double(mpreal a);
long nint(mpreal a);
mpcomplex Real2Complex(mpreal a, mpreal b);
mpreal Cabs1(mpcomplex zdum);

//implementation of sign transfer function.
inline mpreal sign(mpreal a, mpreal b)
{
    mpreal mtmp;
    mtmp = abs(a);
    if (b < 0.0) {
	mtmp = -mtmp;
    }
    return mtmp;
}

inline double cast2double(mpreal a)
{
    double tmp;
    tmp = a;
    return a;
}

inline long nint(mpreal a)
{
    long i;
    mpreal tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = tmp; //cast to long
    return i;
}

inline mpcomplex Real2Complex(mpreal a, mpreal b)
{
    mpcomplex ret(a,b);
    return ret;
}

inline mpreal Cabs1(mpcomplex zdum)
{
    mpreal ret;
    ret = abs(zdum.real()) + abs(zdum.imag());
    return ret;
}

#endif
