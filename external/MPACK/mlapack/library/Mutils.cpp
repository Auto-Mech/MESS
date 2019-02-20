/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Mutils.cpp,v 1.14 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define __MUTILS_CPP__
#if defined ___MPACK_BUILD_WITH_GMP___
REAL log2(REAL x)
{
    double d;
    double ln2_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln2_app = (double) exp + log10(d) / log10(2);
    return ln2_app;
}
#endif

#if defined ___MPACK_BUILD_WITH_QD___
REAL log2(REAL x)
{
    return log10(x) / (qd_real::_log2 / qd_real::_log10);
}
#endif

#if defined ___MPACK_BUILD_WITH_DD___
REAL log2(REAL x)
{
    return log10(x) / (dd_real::_log2 / dd_real::_log10);
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
REAL log(REAL x)
{
    double d;
    double ln_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln_app = (double) exp *log(2.0) + log(d);
    return ln_app;
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
REAL log10(REAL x)
{
    double d;
    double ln10_app;
    signed long int exp;

    d = mpf_get_d_2exp(&exp, x.get_mpf_t());
    ln10_app = (double) exp *log10(2.0) + log10(d);
    return ln10_app;
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
REAL pow(REAL x, REAL y)
{
    REAL mtemp1, mtemp2;
    mtemp1 = y * log(x);
    mtemp2 = exp(mtemp1);
    return mtemp2;
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
REAL cos(REAL x)
{
    REAL mtemp1;
    mtemp1 = cos(x.get_d());
    return mtemp1;
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
REAL sin(REAL x)
{
    REAL mtemp1;
    mtemp1 = sin(x.get_d());
    return mtemp1;
}
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
REAL exp(REAL x)
{
    REAL mtemp1;
    mtemp1 = exp(x.get_d());
    return mtemp1;
}
#endif

COMPLEX exp(COMPLEX x)
{
    REAL ex;
    REAL c;
    REAL s;
    COMPLEX ans;
#if defined ___MPACK_BUILD_WITH_MPFR___
    ex = mpfr::exp(x.real());
    c = mpfr::cos(x.imag());
    s = mpfr::sin(x.imag());
#else
    ex = exp(x.real());
    c = cos(x.imag());
    s = sin(x.imag());
#endif
    ans.real() = ex * c;
    ans.imag() = ex * s;
    return ans;
}

REAL pi(REAL dummy)
{
#if defined ___MPACK_BUILD_WITH_GMP___
    REAL mtemp1;
    mtemp1 = M_PI;
    return mtemp1;
#endif
#if defined ___MPACK_BUILD_WITH_QD___
    return qd_real::_pi;
#endif
#if defined ___MPACK_BUILD_WITH_DD___
    return dd_real::_pi;
#endif
#if defined ___MPACK_BUILD_WITH_MPFR___
    mpfr_free_cache();
    return const_pi();
#endif
#if defined ___MPACK_BUILD_WITH_DOUBLE___
    return M_PI;
#endif
}

#if defined ___MPACK_BUILD_WITH_MPFR___
mpreal __mpfr_sqrt(mpreal z)
{
    mpreal tmp;
    tmp = sqrt(z);
    return tmp;
}

mpcomplex __mpfr_sqrt(mpcomplex z)
{
    mpcomplex tmp;
    tmp = sqrt(z);
    return tmp;
}

mpreal __mpfr_log(mpreal x)
{
    return mpfr::log(x);
}

mpreal __mpfr_pow(mpreal x, mpreal y)
{
    return pow(x, y);
}

mpreal __mpfr_sin(mpreal x)
{
    return mpfr::sin(x);
}

mpreal __mpfr_cos(mpreal x)
{
    return mpfr::cos(x);
}

#endif
