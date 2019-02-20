// Get inverse of Hilbert matrix via Rgetri, using double-double
// Very difficult to calculate the inverse due to very large
// condition number.
// For details, see http://en.wikipedia.org/wiki/Hilbert_matrix
// This file is freely usable.
// written by Nakata Maho, 2009/9/23.

#include <mblas_dd.h>
#include <mlapack_dd.h>
#include <stdio.h>

//Matlab/Octave format
void printmat(int N, int M, dd_real * A, int LDA)
{
    dd_real mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%3.1f", mtmp.x[0]);
	    if (j < M - 1)
		printf(", ");
	}
	if (i < N - 1)
	    printf("]; ");
	else
	    printf("] ");
    }
    printf("]");
}

void inv_hilbert_matrix(int n)
{
    mpackint lwork, info;
    dd_real *A = new dd_real[n * n];	//A is overwritten
    dd_real *Aorg = new dd_real[n * n];
    dd_real *C = new dd_real[n * n];
    mpackint *ipiv = new mpackint[n];
    dd_real One = 1.0, Zero = 0.0, mtmp;

//setting A matrix
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtmp = (i + 1) + (j + 1) - 1;
	    A[i + j * n] = One / mtmp;
	    Aorg[i + j * n] = One / mtmp;
	}
    }

    printf("A = ");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    dd_real *work = new dd_real[1];

    Rgetri(n, A, n, ipiv, work, lwork, &info);
    lwork = int (work[0].x[0]);
    delete[]work;
    work = new dd_real[std::max(1, (int) lwork)];
//inverse matrix
    Rgetrf(n, n, A, n, ipiv, &info);
    Rgetri(n, A, n, ipiv, work, lwork, &info);

    printf("invA = ");
    printmat(n, n, A, n);
    printf("\n");

    One = 1.0, Zero = 0.0;
    Rgemm("N", "N", n, n, n, One, Aorg, n, A, n, Zero, C, n);
    printf("A * invA =");
    printmat(n, n, C, n);
    printf("\n");

    printf("Infnorm(A*invA)=");
    mtmp = 0.0;
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    if (mtmp < abs(C[i + j * n] - ((i == j) ? 1.0 : 0.0)))
		mtmp = abs(C[i + j * n] - ((i == j) ? 1.0 : 0.0));
	}
    }
    printf("%.16e\n", mtmp.x[0]);
    delete[]A;
    delete[]Aorg;
    delete[]C;
    delete[]ipiv;
    delete[]work;
}

int main()
{
    for (int n = 1; n < 21; n++) {
	printf("# inversion of Hilbert matrix of order n=%d\n", n);
	inv_hilbert_matrix(n);
    }
}
