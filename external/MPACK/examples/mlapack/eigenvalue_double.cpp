// Get eigenvalues and eigenvecs
// of Matrix A via Rsyev, using double.
// This file is freely usable.
// written by Nakata Maho, 2009/9/24.

#include <mblas_double.h>
#include <mlapack_double.h>
#include <stdio.h>

//Matlab/Octave format
void printmat(int N, int M, double *A, int LDA)
{
    double mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%8.6e", mtmp);
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

int main()
{
    mpackint n = 3;
    mpackint lwork, info;

    double *A = new double[n * n];
    double *w = new double[n];

//setting A matrix
    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 2;    A[0 + 2 * n] = 3;
    A[1 + 0 * n] = 2;    A[1 + 1 * n] = 5;    A[1 + 2 * n] = 4;
    A[2 + 0 * n] = 3;    A[2 + 1 * n] = 4;    A[2 + 2 * n] = 6;

    printf("A =");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    double *work = new double[1];

    Rsyev("V", "U", n, A, n, w, work, lwork, &info);
    lwork = (int) work[0];
    delete[]work;
    work = new double[std::max((mpackint) 1, lwork)];
//inverse matrix
    Rsyev("V", "U", n, A, n, w, work, lwork, &info);
//print out some results.
    printf("#eigenvalues \n");
    printf("w =");
    printmat(n, 1, w, 1);
    printf("\n");
    printf("#eigenvecs \n");
    printf("U =");
    printmat(n, n, A, n);
    printf("\n");
    printf("#you can check eigenvalues using octave/Matlab by:\n");
    printf("eig(A)\n");
    printf("#you can check eigenvectors using octave/Matlab by:\n");
    printf("U'*A*U\n");

    delete[]work;
    delete[]w;
    delete[]A;
}
