#ifndef BLAS_H
#define BLAS_H

typedef struct { float r, i; } complex;
typedef struct { double r, i; } doublecomplex;

#if defined __cplusplus
extern "C" {
#endif

int
dspmv_(const char& uplo, const Lapack::int_t& N,
       const double& alpha,
       const double* ap, const double* X, const Lapack::int_t& incX, 
       const double& beta, 
       double* Y, const Lapack::int_t& incY);

int
dgemv_(const char& trans, const Lapack::int_t& M, const Lapack::int_t& N,
       const double& alpha, const double* A, const Lapack::int_t& lda,
       const double* X, const Lapack::int_t& incX, const double& beta,
       double* Y, const Lapack::int_t& incY);

int
zgemm_(const char& transA, const char& transB, 
       const Lapack::int_t& M, const Lapack::int_t& N, const Lapack::int_t& K,
       const Lapack::complex& alpha,
       const Lapack::complex* A, const Lapack::int_t& lda,
       const Lapack::complex* B, const Lapack::int_t& ldb,
       const Lapack::complex& beta,
       Lapack::complex* C, const Lapack::int_t& ldc);

int
dgemm_(const char& transA, const char& transB, 
       const Lapack::int_t& M, const Lapack::int_t& N, const Lapack::int_t& K,
       const double& alpha,
       const double* A, const Lapack::int_t& lda,
       const double* B, const Lapack::int_t& ldb,
       const double& beta,
       double* C, const Lapack::int_t& ldc);

float 
sdot_(Lapack::int_t* N, 
      float* X, Lapack::int_t* incX, 
      float* Y, Lapack::int_t* incY);

double
ddot_(Lapack::int_t* N, 
      double* X, Lapack::int_t* incX, 
      double* Y, Lapack::int_t* incY);

void 
cdotu_(complex* retval,
       Lapack::int_t* N, 
       complex* X, Lapack::int_t* incX, 
       complex* Y, Lapack::int_t* incY);

void
cdotc_(complex* retval,
       Lapack::int_t* N, 
       complex* X, Lapack::int_t* incX, 
       complex* Y, Lapack::int_t* incY);

void
zdotu_(doublecomplex* retval,
       Lapack::int_t* N, 
       doublecomplex* X, Lapack::int_t* incX, 
       doublecomplex* Y, Lapack::int_t* incY);

void
zdotc_(doublecomplex* retval,
       Lapack::int_t* N, 
       doublecomplex* X, Lapack::int_t* incX, 
       doublecomplex* Y, Lapack::int_t* incY);

float 
snrm2_(Lapack::int_t* N, 
       float* X, Lapack::int_t* incX);

float
sasum_(Lapack::int_t* N, 
       float* X, Lapack::int_t* incX);

double
dnrm2_(Lapack::int_t* N, 
       double* X, Lapack::int_t* incX);

double
dasum_(Lapack::int_t* N, 
       double* X, Lapack::int_t* incX);

float 
scnrm2_(Lapack::int_t* N, 
        complex* X, Lapack::int_t* incX);

float
scasum_(Lapack::int_t* N, 
        complex* X, Lapack::int_t* incX);

double 
dznrm2_(Lapack::int_t* N, 
        doublecomplex* X, Lapack::int_t* incX);

double
dzasum_(Lapack::int_t* N, 
        doublecomplex* X, Lapack::int_t* incX);

int
isamax_(Lapack::int_t* N,
        float* X, Lapack::int_t* incX);

int
idamax_(Lapack::int_t* N,
        double* X, Lapack::int_t* incX);

int
icamax_(Lapack::int_t* N,
        complex* X, Lapack::int_t* incX);

int
izamax_(Lapack::int_t* N,
        doublecomplex* X, Lapack::int_t* incX);

int
sswap_(Lapack::int_t* N,
       float* X, Lapack::int_t* incX,
       float* Y, Lapack::int_t* incY);

int
scopy_(Lapack::int_t* N,
       float* X, Lapack::int_t* incX,
       float* Y, Lapack::int_t* incY);

int
saxpy_(Lapack::int_t* N,
       float* alpha,
       float* X, Lapack::int_t* incX,
       float* Y, Lapack::int_t* incY);

int
dswap_(Lapack::int_t* N,
       double* X, Lapack::int_t* incX,
       double* Y, Lapack::int_t* incY);

int
dcopy_(Lapack::int_t* N,
       double* X, Lapack::int_t* incX,
       double* Y, Lapack::int_t* incY);

int
daxpy_(Lapack::int_t* N,
       double* alpha,
       double* X, Lapack::int_t* incX,
       double* Y, Lapack::int_t* incY);

int
cswap_(Lapack::int_t* N,
       complex* X, Lapack::int_t* incX,
       complex* Y, Lapack::int_t* incY);

int
ccopy_(Lapack::int_t* N,
       complex* X, Lapack::int_t* incX,
       complex* Y, Lapack::int_t* incY);

int
caxpy_(Lapack::int_t* N,
      complex* alpha,
      complex* X, Lapack::int_t* incX,
      complex* Y, Lapack::int_t* incY);

int
zswap_(Lapack::int_t* N,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY);

int
zcopy_(Lapack::int_t* N,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY);

int
zaxpy_(Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY);

int
srotg_(float* a, float* b, float* c, float* s);

int
srot_(Lapack::int_t* N,
      float* X, Lapack::int_t* incX,
      float* Y, Lapack::int_t* incY,
      float* c, float* s);

int
drotg_(double* a, double* b, double* c, double* s);

int
drot_(Lapack::int_t* N,
      double* X, Lapack::int_t* incX,
      double* Y, Lapack::int_t* incY,
      double* c, double* s);

int
sscal_(Lapack::int_t* N,
       float* alpha,
       float* X, Lapack::int_t* incX);

int
dscal_(Lapack::int_t* N,
       double* alpha,
       double* X, Lapack::int_t* incX);

int
cscal_(Lapack::int_t* N,
       complex* alpha,
       complex* X, Lapack::int_t* incX);

int
zscal_(Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* X, Lapack::int_t* incX);

int
csscal_(Lapack::int_t* N,
        float* alpha,
        complex* X, Lapack::int_t* incX);

int
zdscal_(Lapack::int_t* N,
        double* alpha,
        doublecomplex* X, Lapack::int_t* incX);

int
sgemv_(char* trans, Lapack::int_t* M, Lapack::int_t* N,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* X, Lapack::int_t* incX,
       float* beta,
       float* Y, Lapack::int_t* incY);

int
sgbmv_(char *trans, Lapack::int_t *M, Lapack::int_t *N, Lapack::int_t *KL, Lapack::int_t *KU, 
       float *alpha, 
       float *A, Lapack::int_t *lda, 
       float *X, Lapack::int_t *incX, 
       float *beta, 
       float *Y, Lapack::int_t *incY);

int 
strmv_(char* uplo, char *trans, char* diag, Lapack::int_t *N,  
       float *A, Lapack::int_t *lda, 
       float *X, Lapack::int_t *incX);

int
stbmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       float* A, Lapack::int_t* lda,
       float* X, Lapack::int_t* incX);

int
stpmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       float* Ap, 
       float* X, Lapack::int_t* incX);

int
strsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N,
       float* A, Lapack::int_t* lda,
       float* X, Lapack::int_t* incX);

int
stbsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       float* A, Lapack::int_t* lda, 
       float* X, Lapack::int_t* incX);

int
stpsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       float* Ap, 
       float* X, Lapack::int_t* incX);

int 
dgbmv_(char *trans, Lapack::int_t *M, Lapack::int_t *N, Lapack::int_t *KL, Lapack::int_t *KU, 
       double *alpha, 
       double *A, Lapack::int_t *lda, 
       double *X, Lapack::int_t *incX, 
       double *beta, 
       double *Y, Lapack::int_t *incY);

int 
dtrmv_(char* uplo, char *trans, char* diag, Lapack::int_t *N,  
       double *A, Lapack::int_t *lda, 
       double *X, Lapack::int_t *incX);

int
dtbmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       double* A, Lapack::int_t* lda,
       double* X, Lapack::int_t* incX);

int
dtpmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       double* Ap, 
       double* X, Lapack::int_t* incX);

int
dtrsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N,
       double* A, Lapack::int_t* lda,
       double* X, Lapack::int_t* incX);

int
dtbsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       double* A, Lapack::int_t* lda, 
       double* X, Lapack::int_t* incX);

int
dtpsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       double* Ap, 
       double* X, Lapack::int_t* incX);

int
cgemv_(char* trans, Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* X, Lapack::int_t* incX,
       complex* beta,
       complex* Y, Lapack::int_t* incY);

int 
cgbmv_(char *trans, Lapack::int_t *M, Lapack::int_t *N, Lapack::int_t *KL, Lapack::int_t *KU, 
       complex *alpha, 
       complex *A, Lapack::int_t *lda, 
       complex *X, Lapack::int_t *incX, 
       complex *beta, 
       complex *Y, Lapack::int_t *incY);

int 
ctrmv_(char* uplo, char *trans, char* diag, Lapack::int_t *N,  
       complex *A, Lapack::int_t *lda, 
       complex *X, Lapack::int_t *incX);

int
ctbmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       complex* A, Lapack::int_t* lda,
       complex* X, Lapack::int_t* incX);

int
ctpmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       complex* Ap, 
       complex* X, Lapack::int_t* incX);

int
ctrsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N,
       complex* A, Lapack::int_t* lda,
       complex* X, Lapack::int_t* incX);

int
ctbsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       complex* A, Lapack::int_t* lda, 
       complex* X, Lapack::int_t* incX);

int
ctpsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       complex* Ap, 
       complex* X, Lapack::int_t* incX);

int
zgemv_(char* trans, Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* beta,
       doublecomplex* Y, Lapack::int_t* incY);

int 
zgbmv_(char *trans, Lapack::int_t *M, Lapack::int_t *N, Lapack::int_t *KL, Lapack::int_t *KU, 
       doublecomplex *alpha, 
       doublecomplex *A, Lapack::int_t *lda, 
       doublecomplex *X, Lapack::int_t *incX, 
       doublecomplex *beta, 
       doublecomplex *Y, Lapack::int_t *incY);

int 
ztrmv_(char* uplo, char *trans, char* diag, Lapack::int_t *N,  
       doublecomplex *A, Lapack::int_t *lda, 
       doublecomplex *X, Lapack::int_t *incX);

int
ztbmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* X, Lapack::int_t* incX);

 void  
ztpmv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
      doublecomplex* Ap, 
      doublecomplex* X, Lapack::int_t* incX);

int
ztrsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* X, Lapack::int_t* incX);

int
ztbsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, Lapack::int_t* K,
       doublecomplex* A, Lapack::int_t* lda, 
       doublecomplex* X, Lapack::int_t* incX);

int
ztpsv_(char* uplo, char* trans, char* diag, Lapack::int_t* N, 
       doublecomplex* Ap, 
       doublecomplex* X, Lapack::int_t* incX);

int
ssymv_(char* uplo, Lapack::int_t* N,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* X, Lapack::int_t* incX,
       float* beta,
       float* Y, Lapack::int_t* incY);

int 
ssbmv_(char* uplo, Lapack::int_t* N, Lapack::int_t* K,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* X, Lapack::int_t* incX,
       float* beta,
       float* Y, Lapack::int_t* incY);

int
sspmv_(char* uplo, Lapack::int_t* N,
       float* alpha,
       float* Ap,
       float* X, Lapack::int_t* incX,
       float* beta,
       float* Y, Lapack::int_t* incY);

int
sger_(Lapack::int_t* M, Lapack::int_t* N,
      float* alpha,
      float* X, Lapack::int_t* incX,
      float* Y, Lapack::int_t* incY,
      float* A, Lapack::int_t* lda);

int
ssyr_(char* uplo, Lapack::int_t* N,
      float* alpha,
      float* X, Lapack::int_t* incX,
      float* A, Lapack::int_t* lda);

int
sspr_(char* uplo, Lapack::int_t* N,
      float* alpha,
      float* X, Lapack::int_t* incX,
      float* Ap);

int
ssyr2_(char* uplo, Lapack::int_t* N,
       float* alpha,
       float* X, Lapack::int_t* incX,
       float* Y, Lapack::int_t* incY,
       float* A, Lapack::int_t* lda);

int
sspr2_(char* uplo, Lapack::int_t* N,
       float* alpha, 
       float* X, Lapack::int_t* incX,
       float* Y, Lapack::int_t* incY,
       float* A);

int
dsymv_(char* uplo, Lapack::int_t* N,
       double* alpha,
       double* A, Lapack::int_t* lda,
       double* X, Lapack::int_t* incX,
       double* beta,
       double* Y, Lapack::int_t* incY);

int 
dsbmv_(char* uplo, Lapack::int_t* N, Lapack::int_t* K,
       double* alpha,
       double* A, Lapack::int_t* lda,
       double* X, Lapack::int_t* incX,
       double* beta,
       double* Y, Lapack::int_t* incY);

int
dger_(Lapack::int_t* M, Lapack::int_t* N,
      double* alpha,
      double* X, Lapack::int_t* incX,
      double* Y, Lapack::int_t* incY,
      double* A, Lapack::int_t* lda);

int
dsyr_(char* uplo, Lapack::int_t* N,
      double* alpha,
      double* X, Lapack::int_t* incX,
      double* A, Lapack::int_t* lda);

int
dspr_(char* uplo, Lapack::int_t* N,
      double* alpha,
      double* X, Lapack::int_t* incX,
      double* Ap);

int
dsyr2_(char* uplo, Lapack::int_t* N,
       double* alpha,
       double* X, Lapack::int_t* incX,
       double* Y, Lapack::int_t* incY,
       double* A, Lapack::int_t* lda);

int
dspr2_(char* uplo, Lapack::int_t* N,
       double* alpha, 
       double* X, Lapack::int_t* incX,
       double* Y, Lapack::int_t* incY,
       double* A);

int
chemv_(char* uplo, Lapack::int_t* N,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* X, Lapack::int_t* incX,
       complex* beta,
       complex* Y, Lapack::int_t* incY);

int
chbmv_(char* uplo, Lapack::int_t* N, Lapack::int_t* K,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* X, Lapack::int_t* incX,
       complex* beta,
       complex* Y, Lapack::int_t* incY);

int
chpmv_(char* uplo, Lapack::int_t* N, 
       complex* alpha,
       complex* Ap, 
       complex* X, Lapack::int_t* incX,
       complex* beta,
       complex* Y, Lapack::int_t* incY);

int
cgeru_(Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* X, Lapack::int_t* incX,
       complex* Y, Lapack::int_t* incY,
       complex* A, Lapack::int_t* lda);

int
cgerc_(Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* X, Lapack::int_t* incX,
       complex* Y, Lapack::int_t* incY,
       complex* A, Lapack::int_t* lda);

int
cher_(char* uplo, Lapack::int_t* N,
      float* alpha,
      complex* X, Lapack::int_t* incX,
      complex* A, Lapack::int_t* lda);

int
chpr_(char* uplo, Lapack::int_t* N,
      float* alpha,
      complex* X, Lapack::int_t* incX,
      complex* Ap);

int
cher2_(char* uplo, Lapack::int_t* N,
       complex* alpha,
       complex* X, Lapack::int_t* incX,
       complex* Y, Lapack::int_t* incY,
       complex* A, Lapack::int_t* lda);

int
chpr2_(char* uplo, Lapack::int_t* N,
       complex* alpha,
       complex* X, Lapack::int_t* incX,
       complex* Y, Lapack::int_t* incY,
       complex* Ap);

int
zhemv_(char* uplo, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* beta,
       doublecomplex* Y, Lapack::int_t* incY);

int
zhbmv_(char* uplo, Lapack::int_t* N, Lapack::int_t* K,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* beta,
       doublecomplex* Y, Lapack::int_t* incY);

int
zhpmv_(char* uplo, Lapack::int_t* N, 
       doublecomplex* alpha,
       doublecomplex* Ap, 
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* beta,
       doublecomplex* Y, Lapack::int_t* incY);

int
zgeru_(Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY,
       doublecomplex* A, Lapack::int_t* lda);

int
zgerc_(Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY,
       doublecomplex* A, Lapack::int_t* lda);

int
zher_(char* uplo, Lapack::int_t* N,
      double* alpha,
      doublecomplex* X, Lapack::int_t* incX,
      doublecomplex* A, Lapack::int_t* lda);

int
zhpr_(char* uplo, Lapack::int_t* N,
      double* alpha,
      doublecomplex* X, Lapack::int_t* incX,
      doublecomplex* Ap);

int
zher2_(char* uplo, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY,
       doublecomplex* A, Lapack::int_t* lda);

int
zhpr2_(char* uplo, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* X, Lapack::int_t* incX,
       doublecomplex* Y, Lapack::int_t* incY,
       doublecomplex* Ap);

int
sgemm_(char* transA, char* transB, Lapack::int_t* M, Lapack::int_t* N, Lapack::int_t* K,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* B, Lapack::int_t* ldb,
       float* beta,
       float* C, Lapack::int_t* ldc);

int
ssymm_(char* side, char* uplo, Lapack::int_t* M, Lapack::int_t* N,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* B, Lapack::int_t* ldb,
       float* beta,
       float* C, Lapack::int_t* ldc);

int
ssyrk_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* beta,
       float* C, Lapack::int_t* ldc);

int
ssyr2k_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
        float* alpha,
        float* A, Lapack::int_t* lda,
        float* B, Lapack::int_t* ldb,
        float* beta,
        float* C, Lapack::int_t* ldc);

int
strmm_(char* side, char* uplo, char* trans, char* diag, 
       Lapack::int_t* M, Lapack::int_t* N,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* B, Lapack::int_t* ldb);

int 
strsm_(char* side, char* uplo, char* trans, char* diag,
       Lapack::int_t* M, Lapack::int_t* N,
       float* alpha,
       float* A, Lapack::int_t* lda,
       float* B, Lapack::int_t* ldb);


int
dsymm_(char* side, char* uplo, Lapack::int_t* M, Lapack::int_t* N,
       double* alpha,
       double* A, Lapack::int_t* lda,
       double* B, Lapack::int_t* ldb,
       double* beta,
       double* C, Lapack::int_t* ldc);

int
dsyrk_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
       double* alpha,
       double* A, Lapack::int_t* lda,
       double* beta,
       double* C, Lapack::int_t* ldc);

int
dsyr2k_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
        double* alpha,
        double* A, Lapack::int_t* lda,
        double* B, Lapack::int_t* ldb,
        double* beta,
        double* C, Lapack::int_t* ldc);

int
dtrmm_(char* side, char* uplo, char* trans, char* diag, 
       Lapack::int_t* M, Lapack::int_t* N,
       double* alpha,
       double* A, Lapack::int_t* lda,
       double* B, Lapack::int_t* ldb);

int 
dtrsm_(char* side, char* uplo, char* trans, char* diag,
       Lapack::int_t* M, Lapack::int_t* N,
       double* alpha,
       double* A, Lapack::int_t* lda,
       double* B, Lapack::int_t* ldb);

int
cgemm_(char* transA, char* transB, Lapack::int_t* M, Lapack::int_t* N, Lapack::int_t* K,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* B, Lapack::int_t* ldb,
       complex* beta,
       complex* C, Lapack::int_t* ldc);

int
csymm_(char* side, char* uplo, Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* B, Lapack::int_t* ldb,
       complex* beta,
       complex* C, Lapack::int_t* ldc);

int
csyrk_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* beta,
       complex* C, Lapack::int_t* ldc);

int
csyr2k_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
        complex* alpha,
        complex* A, Lapack::int_t* lda,
        complex* B, Lapack::int_t* ldb,
        complex* beta,
        complex* C, Lapack::int_t* ldc);

int
ctrmm_(char* side, char* uplo, char* trans, char* diag, 
       Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* B, Lapack::int_t* ldb);

int 
ctrsm_(char* side, char* uplo, char* trans, char* diag,
       Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* B, Lapack::int_t* ldb);


int
zsymm_(char* side, char* uplo, Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* B, Lapack::int_t* ldb,
       doublecomplex* beta,
       doublecomplex* C, Lapack::int_t* ldc);

int
zsyrk_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* beta,
       doublecomplex* C, Lapack::int_t* ldc);

int
zsyr2k_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
        doublecomplex* alpha,
        doublecomplex* A, Lapack::int_t* lda,
        doublecomplex* B, Lapack::int_t* ldb,
        doublecomplex* beta,
        doublecomplex* C, Lapack::int_t* ldc);

int
ztrmm_(char* side, char* uplo, char* trans, char* diag, 
       Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* B, Lapack::int_t* ldb);

int 
ztrsm_(char* side, char* uplo, char* trans, char* diag,
       Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* B, Lapack::int_t* ldb);

int
chemm_(char* side, char* uplo, Lapack::int_t* M, Lapack::int_t* N,
       complex* alpha,
       complex* A, Lapack::int_t* lda,
       complex* B, Lapack::int_t* ldb,
       complex* beta,
       complex* C, Lapack::int_t* ldc);

int
cherk_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
       float* alpha,
       complex* A, Lapack::int_t* lda,
       float* beta,
       complex* C, Lapack::int_t* ldc);

int
cher2k_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
        complex* alpha,
        complex* A, Lapack::int_t* lda,
        complex* B, Lapack::int_t* ldb,
        float* beta,
        complex* C, Lapack::int_t* ldc);

int
zhemm_(char* side, char* uplo, Lapack::int_t* M, Lapack::int_t* N,
       doublecomplex* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       doublecomplex* B, Lapack::int_t* ldb,
       doublecomplex* beta,
       doublecomplex* C, Lapack::int_t* ldc);

int
zherk_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
       double* alpha,
       doublecomplex* A, Lapack::int_t* lda,
       double* beta,
       doublecomplex* C, Lapack::int_t* ldc);

int
zher2k_(char* uplo, char* trans, Lapack::int_t* N, Lapack::int_t* K,
        doublecomplex* alpha,
        doublecomplex* A, Lapack::int_t* lda,
        doublecomplex* B, Lapack::int_t* ldb,
        double* beta,
        doublecomplex* C, Lapack::int_t* ldc);

#if defined __cplusplus
}
#endif

#endif /* __BLAS_H */

