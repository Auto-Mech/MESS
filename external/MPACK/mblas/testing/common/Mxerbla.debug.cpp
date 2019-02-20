#include <mblas.h>
#include <blas.h>
#include <mpack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

void Mxerbla_test()
{
    Mxerbla("Fasum", 10);
    Mxerbla("Maho", 100);
}

int main(int argc, char *argv[])
{
    Mxerbla_test();
    printf("Mxerbla test passed...\n");
    return (0);
}
