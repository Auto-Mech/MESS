#ifndef __LAPACK_H
#define __LAPACK_H

#if defined __cplusplus
extern "C" {
#endif

  int dspev_(const char& jobz, const char& uplo, const Lapack::int_t& n, double* ap,
	     double* w, double* z, const Lapack::int_t& ldz, double* work, Lapack::int_t& info);

  int zhpev_(const char& jobz, const char& uplo, const Lapack::int_t& n, Lapack::complex* ap,
	     double* w, Lapack::complex* z, const Lapack::int_t& ldz, Lapack::complex* work, 
	     double* rwork, Lapack::int_t& info);

  int dpptrf_(const char& uplo, const Lapack::int_t& n, double* ap, Lapack::int_t& info);

  int dpptri_(const char& uplo, const Lapack::int_t& n, double* ap, Lapack::int_t& info);

  int dpptrs_(const char& uplo, const Lapack::int_t& n, const Lapack::int_t& nrhs, 
	      const double* ap, double* b, const Lapack::int_t& ldb, Lapack::int_t& info);

  int dspgvd_(const Lapack::int_t& itype, const char& job, const char& uplo, const Lapack::int_t& n, 
	      double* a, double* b, double* w, double* z, const Lapack::int_t& ldz, 
	      double* work, const Lapack::int_t& lwork, Lapack::int_t* iwork, const Lapack::int_t& liwork, 
	      Lapack::int_t& info);

  // simultaneous diagonalization of two hermitian matrices, B should be positively defined
  //
  int zhpgvd_(const Lapack::int_t& itype, // calculation type: 1 - A * x = lamda * B * x; 2 - A * B * x = lamda * x; 3 - B * A * x = lamda * x
	      //
	      const char&          job,   // job type: V - eigenvalues and eigenvectors; N - eigenvalues only
	      //
	      const char&          uplo,  // matrix packing type: U - upper triangle packing; L - low triangle packing
	      //
	      const Lapack::int_t& n,     // matrices order
	      //
	      Lapack::complex*     a,     // A matrix
	      //
	      Lapack::complex*     b,     // B matrix
	      //
	      double*              w,     // eigenvalues
	      //
	      Lapack::complex*     z,     // eigenvectors
	      //
	      const Lapack::int_t& ldz,   // leading dimeinsion of matrix z
	      //
	      Lapack::complex*     work,  // complex workspace
	      //
	      const Lapack::int_t& lwork, // complex workspace size
	      //
	      double*              rwork, // real workspace
	      //
	      const Lapack::int_t& lrwork,// real workspace size
	      //
	      Lapack::int_t*       iwork, // integer workspace
	      //
	      const Lapack::int_t& liwork,// integer workspace size
	      //
	      Lapack::int_t&       info   // exit status
	      );

  int dsbevd_(const char& job, const char& uplo, const Lapack::int_t& n, const Lapack::int_t& kd, 
	      double* ab, const Lapack::int_t& ldab, double* w, double* z, const Lapack::int_t& ldz,
	      double* work, const Lapack::int_t& lwork, Lapack::int_t* iwork, const Lapack::int_t& liwork,
	      Lapack::int_t& info);
  
  int dgesv_(const Lapack::int_t& n, const Lapack::int_t& nrhs, double* a, const Lapack::int_t& lda,
	     Lapack::int_t* ipiv, double* b, const Lapack::int_t& ldb, Lapack::int_t& info);

  int dsyev_(const char& jobz, const char& uplo, const Lapack::int_t& n, double* a, const Lapack::int_t& lda, 
	     double* w, double* work, const Lapack::int_t& lwork, Lapack::int_t& info);

  int dsyevd_(const char& jobz, const char& uplo, const Lapack::int_t& n, double* a, const Lapack::int_t& lda, 
	     double* w, double* work, const Lapack::int_t& lwork, Lapack::int_t* iwork, const Lapack::int_t& liwork, 
	     Lapack::int_t& info);

  int dspsv_(const char& uplo, const Lapack::int_t& n, const Lapack::int_t& nrhs, double* ap,
	     Lapack::int_t* ipiv, double* b, const Lapack::int_t& ldb, Lapack::int_t& info);

  int dppsv_(const char& uplo, const Lapack::int_t& n, const Lapack::int_t& nrhs, double* ap,
	     double* b, const Lapack::int_t& ldb, Lapack::int_t& info);

  int dgetrf_(const Lapack::int_t& m, const Lapack::int_t& n, double* a, const Lapack::int_t& lda,
	      Lapack::int_t* ipiv, Lapack::int_t& info);

  int dgetri_(const Lapack::int_t& n, double* a, const Lapack::int_t& lda, const Lapack::int_t* ipiv,
	      double* work, const Lapack::int_t& lwork, Lapack::int_t& info);

  int dgetrs_(const char& trans, const Lapack::int_t& n, const Lapack::int_t& nrhs,
	      const double* a, const Lapack::int_t& lda, const Lapack::int_t* ipiv,
	      double* b, const Lapack::int_t& ldb, Lapack::int_t& info);

  int dsptrf_(const char& uplo, const Lapack::int_t& n, double* ap, Lapack::int_t* ipiv, 
	      Lapack::int_t& info);

  int dsptri_(const char& uplo, const Lapack::int_t& n, double* ap, const Lapack::int_t* ipiv,
	      double* work, Lapack::int_t& info);

  int dsptrs_(const char& uplo, const Lapack::int_t& n, const Lapack::int_t& nrhs, const double* ap,
	      const Lapack::int_t* ipiv, double* b, const Lapack::int_t& ldb, Lapack::int_t& info);

  // solving system of linear equations via singular value decomposition
  //
  void dgelsd_(const Lapack::int_t& m,
	       const Lapack::int_t& n,
	       const Lapack::int_t& nrhs,
	       double*              a,
	       const Lapack::int_t& lda,
	       double*              b,
	       const Lapack::int_t& ldb,
	       double*              s,
	       const double&        rcond,
	       Lapack::int_t&       rank,
	       double*              work,
	       const Lapack::int_t& lwork,
	       Lapack::int_t*       iwork,
	       Lapack::int_t&       info);

  // singular value decomposition by divide-and-conquer method
  //
  void dgesdd_ (const char&          jobz,
		const Lapack::int_t& m,
		const Lapack::int_t& n,
		double*              a,
		const Lapack::int_t& lda,
		double*              s,
		double*              u,
		const Lapack::int_t& ldu,
		double*              vt,
		const Lapack::int_t& ldvt,
		double*              work,
		const Lapack::int_t& lwork,
		Lapack::int_t*       iwork,
		Lapack::int_t&       info);

  // singular value decomposition
  //
  void dgesvd_ (const char&          jobu,
		const char&          jobvt,
		const Lapack::int_t& m,
		const Lapack::int_t& n,
		double*              a,
		const Lapack::int_t& lda,
		double*              s,
		double*              u,
		const Lapack::int_t& ldu,
		double*              vt,
		const Lapack::int_t& ldvt,
		double*              work,
		const Lapack::int_t& lwork,
		Lapack::int_t&       info);

#if defined __cplusplus
}
#endif

#endif /* __LAPACK_H */

