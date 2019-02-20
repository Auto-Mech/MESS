/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
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
#include <complex>
#include <mblas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MAX_ITER 10

void mp_sub_test2()
{
    COMPLEX_REF Ctemp1r, Ctemp2r, dummy;
    COMPLEX Ctemp1;
    REAL_REF diff;

    set_random_number(Ctemp1r, Ctemp1);

    cout << "C1R  = "; printnum(Ctemp1r); cout << endl;
    cout << "C1   = "; printnum(Ctemp1);  cout << endl;

    Ctemp2r = Ctemp1 - Ctemp1r;
    diff = abs(Ctemp2r);
    cout << "diff = "; printnum(diff); cout << endl;
}

void mp_sub_test1()
{
    REAL_REF Rtemp1r, Rtemp2r, Rtemp3r;
    REAL Rtemp1;

    printf("Substitution test\n");
    set_random_number(Rtemp1r, Rtemp1);
    printnum(Rtemp1r); printf("\n");
    printnum(Rtemp1);  printf("\n");
    Rtemp2r = Rtemp1 - Rtemp1r;
    printnum(Rtemp2r); printf("\n");
    Rtemp3r = Rtemp1r - Rtemp1;
    printnum(Rtemp3r); printf("\n");
    printf("Substitution test done\n");

    printf("Subtraction test\n");
    Rtemp1r = 1.0;
    Rtemp1 = 2.0;
    Rtemp2r = Rtemp1r - Rtemp1;    
    printnum(Rtemp2r);  printf("\n");

    Rtemp2r = Rtemp1 - Rtemp1r;
    printnum(Rtemp2r);  printf("\n");
    printf("Subtraction test done\n");
}

int main(int argc, char *argv[])
{
    mp_sub_test1();
    mp_sub_test2();
    return (0);
}
