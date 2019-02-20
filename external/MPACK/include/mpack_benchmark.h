/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpack_benchmark.h,v 1.4 2010/08/07 03:15:46 nakatamaho Exp $
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

#include <sys/types.h>
#include <sys/time.h>
#if !defined  _WIN32
#include <sys/resource.h>
#endif

#if defined ___MPACK_BUILD_WITH_DD___
dd_real mpf_randomnumber(dd_real dummy)
{
    dd_real mtmp;
    mtmp = ddrand();            //uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___MPACK_BUILD_WITH_QD___
qd_real mpf_randomnumber(qd_real dummy)
{
    qd_real mtmp;
    mtmp = qdrand();            //uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
mpf_class mpf_randomnumber(mpf_class dummy)
{
    mpf_class mtmp;
 
    mtmp = uniformrandomstate_gmp->get_f();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___MPACK_BUILD_WITH_MPFR___
mpreal mpf_randomnumber(mpreal dummy)
{
    mpreal mtmp;

    mtmp = urandomb(uniformrandomstate_mpfr);
    mtmp = 2.0 * mtmp - 1.0;

    return mtmp;
}
#endif

double mpf_randomnumber(double dummy)
{
#if defined _WIN32 //XXX
    double mtmp = (double)rand();
#else
    double mtmp = drand48();
#endif
    return mtmp;
}

#if !defined _WIN32
unsigned long microseconds(void)
{
    rusage  t;
    timeval tv;
    getrusage( RUSAGE_SELF, &t );
    tv = t.ru_utime;
    return ((unsigned long)tv.tv_sec)*1000000 + tv.tv_usec;
}
#endif

inline double gettime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

