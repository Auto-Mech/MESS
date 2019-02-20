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
#include <mblas.h>
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

void Mlsame_test()
{
//  char a="A";
//  char b="A";
    int errorflag = FALSE;

    if (Mlsame("A", "A"))
	printf("same letter ok\n");
    else {
	printf("same letter NG\n");
	errorflag = TRUE;
    }
    if (Mlsame("A", "a"))
	printf("Uppercase/lowercase ok\n");
    else {
	printf("Uppercase/lowercase NG\n");
	errorflag = TRUE;
    }

    if (Mlsame("A", "aho"))
	printf("Only looks for the first char ok\n");
    else {
	printf("Only looks for the first char NG\n");
	errorflag = TRUE;
    }

    if (Mlsame("Inv", "inv"))
	printf("Only looks for the first char ok\n");
    else {
	printf("Only looks for the first char NG\n");
	errorflag = TRUE;
    }

    if (errorflag == TRUE) {
	printf("Mlsame test failed...\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    Mlsame_test();
    printf("Mlsame test passed...\n");
    return (0);
}
