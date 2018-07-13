#include "lapack.hh"
#include "error.hh"
#include "io.hh"
#include "blas.h"
#include "lapack.h"

#include <iostream>

/****************************************************************
 ************************** Vector ******************************
 ****************************************************************/

double Lapack::Vector::operator* (const Vector& v) const throw(Error::General)
{
  const char funame [] = "Lapack::Vector::operator*: ";

  if(!isinit() || !v.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(v.size() != size()) {
    std::cerr << funame << "the dimensions are different:"
	      << " size1 = " << size()
	      << " size2 = " << v.size()
	      << std::endl;
    throw Error::Range();
  }
  
  const double* p  = *this;
  const double* p1 = v;
  double res = 0.;
  for(int_t i = 0; i < size(); ++i)
    res += *p++ * *p1++;
  return res;
}

double Lapack::Vector::vdot () const
{
  const char funame [] = "Lapack::Vector::vdot: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  double res = 0.;
  for(const double* p = begin(); p != end(); ++p)
    res += *p * *p;
  return res;
}

Lapack::Vector Lapack::Vector::operator* (const Matrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::Vector::operator*: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
  
  if(m.size1() != size()) {
    std::cerr << funame << "the dimensions are different\n";
    throw Error::Range();
  }
  
  Vector res(m.size2());
  dgemv_('T', size(), m.size2(), 1.,  m, size(), *this, 1, 0., res, 1);
  return res;
}

Lapack::Vector Lapack::Vector::operator* (const SymmetricMatrix& m)
  const throw(Error::General)
{
  const char funame [] = "Lapack::Vector::operator*: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
  
  if(size() != m.size()) {
    std::cerr << funame << "the dimensions are different:"
	      << " vector_size = " << size()
	      << " matrix_size = " << m.size()
	      << std::endl;
    throw Error::Range();
  }
  
  Vector res(size());
  dspmv_('U', size(), 1., m,  *this, 1, 0., res, 1);
  return res;
}

/****************************************************************
 ************************** Matrix ******************************
 ****************************************************************/

Lapack::Matrix Lapack::Matrix::transpose () const
{
  const char funame [] = "Lapack::Matrix::transpose: ";
  
  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
  
  Matrix res(size2(), size1());
  for(int_t i = 0; i < size1(); ++i) {
    res.column(i) = row(i);
  }
  return res;
}

Lapack::Matrix Lapack::Matrix::operator* (const Matrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size2() != m.size1()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Range();
  }

  Matrix res(size1(), m.size2());
  dgemm_('N', 'N', size1(), m.size2(), size2(), 1., 
	 *this, size1(), m, m.size1(), 0., res, size1());

  return res;
}

Lapack::Matrix Lapack::Matrix::operator* (const SymmetricMatrix& m)
  const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
  
  if(size2() != m.size()) {
    std::cerr << funame << "dimensions mismatch:"
	      << " mat1_size2 = " << size2()
	      << " mat2_size = " << m.size()
	      << "\n";
    throw Error::Range();
  }

  Matrix res(size1(), size2());
  for(int_t i = 0; i < size1(); ++i)
    dspmv_('U', size2(),  1., m, *this + i, size1(), 0., res + i, size1());
  return res;
}

Lapack::Vector Lapack::Matrix::operator* (const double* v) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::operator*: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  Vector res(size1());
  dgemv_('N', size1(), size2(), 1.,  *this, size1(), v, 1, 0., res, 1);
  return res;
}

Lapack::Vector Lapack::operator* (const double* v, const Matrix& m) throw(Error::General)
{
  const char funame [] = "Lapack::operator*: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  Vector res(m.size2());
  dgemv_('T', m.size1(), m.size2(), 1., m, m.size1(), v, 1, 0., res, 1);
  //for(int_t i = 0; i < m.size2(); ++i)
  //res[i] = m.column(i) * v;
  return res;
}

// inverse matrix
Lapack::Matrix Lapack::Matrix::invert () const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::invert: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size1() != size2()) {
    std::cerr << funame << "not square\n";
    throw Error::Logic();
  }

  Matrix lu = copy();
  Array<int_t> ipiv (size1());
  Matrix res(size1(), size1());
  res = 0.;
  res.diagonal() = 1.;

  int_t info = 0;
  dgesv_(size1(), size1(), lu, size1(), ipiv, res, size1(), info);

#ifdef DEBUG

  double dtemp;

  // matrix inversion check
  Matrix unit_mat = *this * res;
  double maxel = 0.;
  for(int_t i = 0; i < unit_mat.size1(); ++i)
    for(int_t j = 0; j < unit_mat.size2(); ++j)
      if(i != j) {
	dtemp = unit_mat(i, j);
	if(dtemp < 0.)
	  dtemp = -dtemp;
	if(dtemp > maxel)
	  maxel = dtemp;
      }
  IO::log << IO::log_offset << funame 
	  << "matrix inversion check: maximal nondiagonal element = "
	  << maxel << "\n";

#endif
      
 if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dgesv: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgesv: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe  factorization has  been  completed,  but the factor U is\n"
      "\texactly singular, so the solution could not be computed.\n";
    throw Error::Math();
  }
}

// solve linear equations
Lapack::Vector Lapack::Matrix::invert (const Vector& v) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::invert: ";

  if(!isinit() || !v.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size1() != size2()) {
    std::cerr << funame << "not square\n";
    throw Error::Logic();
  }

  if(size1() != v.size()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Logic();
  }

  Matrix lu = copy();
  Array<int_t> ipiv(size1());
  Vector res = v.copy();

  int_t info = 0;
  dgesv_(size1(), 1, lu, size1(), ipiv, res, size1(), info);

 if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dgesv: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgesv: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe  factorization has  been  completed,  but the factor U is\n"
      "\texactly singular, so the solution could not be computed.\n";
    throw Error::Math();
  }
}

// solve linear equations
Lapack::Matrix Lapack::Matrix::invert (const Matrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::invert: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size1() != size2()) {
    std::cerr << funame << "not square\n";
    throw Error::Logic();
  }

  if(size1() != m.size1()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Logic();
  }

  Matrix lu = copy();
  Array<int_t> ipiv(size1());
  Matrix res = m.copy();

  int_t info = 0;
  dgesv_(size1(), m.size2(), lu, size1(), ipiv, res, size1(), info);

 if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dgesv: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgesv: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe  factorization has  been  completed,  but the factor U is\n"
      "\texactly singular, so the solution could not be computed.\n";
    throw Error::Math();
  }
}

  // matrix row
Slice<double> Lapack::Matrix::row (int_t i) throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::row: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i < 0 || i >= size1()) {
    std::cerr << funame << "out of range: index = " << i << "\n";
    throw Error::Range();
  }

  return Slice<double>(*this + i, size2(), size1());
}

ConstSlice<double> Lapack::Matrix::row (int_t i) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::row: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i < 0 || i >= size1()) {
    std::cerr << funame << "out of range: index = " << i << "\n";
    throw Error::Range();
  }

  return ConstSlice<double>(*this + i, size2(), size1());
}

  // matrix column
Slice<double> Lapack::Matrix::column (int_t i) throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::column: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i < 0 || i >= size2()) {
    std::cerr << funame << "out of range: index = " << i << "\n";
    throw Error::Range();
  }

  return Slice<double>(*this + i * size1(), size1());
}

ConstSlice<double> Lapack::Matrix::column (int_t i) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::column: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i < 0 || i >= size2()) {
    std::cerr << funame << "out of range: index = " << i << "\n";
    throw Error::Range();
  }

  return ConstSlice<double>(*this + i * size1(), size1());
}

  // matrix diagonal
Slice<double> Lapack::Matrix::diagonal (int_t i) throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::diagonal: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i <= -size1() || i >= size2()) {
    std::cerr << funame << "out of range: index = " << i << "\n";
    throw Error::Range();
  }

  int_t sz;
  if(i > 0) {// upper diagonal
    sz = size1() < size2() - i ? size1() : size2() - i;
    return Slice<double>(*this + i * size1(), sz, size1() + 1);
  }
  else {// lower diagonal
    sz = size1() + i < size2() ? size1() + i : size2();
    return Slice<double>(*this - i, sz, size1() + 1);
  }
}

ConstSlice<double> Lapack::Matrix::diagonal (int_t i) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::diagonal: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i <= -size1() || i >= size2()) {
    std::cerr << funame << "out of range: index = " << i << "\n";
    throw Error::Range();
  }

  int_t sz;
  if(i > 0) {// upper diagonal
    sz = size1() < size2() - i ? size1() : size2() - i;
    return ConstSlice<double>(*this + i * size1(), sz, size1() + 1);
  }
  else {// lower diagonal
    sz = size1() + i < size2() ? size1() + i : size2();
    return ConstSlice<double>(*this - i, sz, size1() + 1);
  }
}

/****************************************************************
 ******************* Band Symmetric Matrix **********************
 ****************************************************************/

double  Lapack::BandMatrix::operator() (int_t i, int_t j) const
{
  const char funame [] = "Lapack::BandMatrix::operator(): ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i > j)
    std::swap(i, j);

  if(j - i < band_size())
    return Matrix::operator()(band_size() - 1 + i - j, j);

  return 0.;
}

double& Lapack::BandMatrix::operator() (int_t i, int_t j) throw(Error::General)
{
  const char funame [] = "Lapack::BandMatrix::operator(): ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(i > j)
    std::swap(i, j);
 
  if(j - i < band_size())
    return Matrix::operator()(band_size() - 1 + i - j, j);

  std::cerr << funame << "out of range\n";
  throw Error::Range();

}

Lapack::Vector Lapack::BandMatrix::eigenvalues (Matrix* evec)
  const throw(Error::General)
{
  const char funame [] = "Lapack::BandMatrix::eigenvalues: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  char job;
  double* z;
  int_t lwork = -1;
  int_t liwork = -1;
  if(!evec) {
    job = 'N';
    z = 0;
    lwork = 2 * size();
    liwork = 1;
  }
  else {
    job = 'V';
    evec->resize(size());
    z = *evec;
    lwork = 1 + 5 * size() + 2 * size() * size();
    liwork = 3 + 5 * size();
  }

  Vector work(lwork);
  Array<int_t> iwork(liwork);

  BandMatrix cp = copy();
  Vector res(size());
  int_t      info = 0;
  dsbevd_(job, 'U', size(), band_size() - 1, cp, band_size(), res, z, size(),
	  work, lwork, iwork, liwork, info);
 
  if(!info)
    return res;

  else if(info < 0) {
    std::cerr << funame << "dsbevd: " << -info 
	      << "-th argument has an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dsbevd: " << info 
	      << "-th off-diagonal  elements  of an "
      "intermediate tridiagonal form did not converge to zero\n";
    throw Error::Math();
  }
}

/****************************************************************
 ********************* Symmetric Matrix *************************
 ****************************************************************/

Lapack::SymmetricMatrix::SymmetricMatrix 
(const Lapack::Matrix& m, char uplo) throw(Error::General)
  : RefArr<double>(m.size1() * (m.size1() + 1) / 2), _size(new int_t(m.size1()))
{
  const char funame [] = "Lapack::SymmetricMatrix::SymmetricMatrix: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(m.size1() != m.size2()) {
    std::cerr << funame << "not square:"
              << " size1 = " << m.size1()
              << " size2 = " << m.size2()
              << "\n";
    throw Error::Range();
  }

  switch(uplo) {
  case 'U':
    for(int_t i = 0; i < size(); ++i) {
      Slice<double> res(*this + i * (i + 1) / 2, i + 1);
      ConstSlice<double> col(m + i * size(), i + 1);
      res = col;
    }
    return;

  case 'L':
    for(int_t i = 0; i < size(); ++i) {
      Slice<double> res(*this + i * (i + 1) / 2, i + 1);
      ConstSlice<double> row(m + i, i + 1, size());
      res = row;
    }
    return;
  default:
    std::cerr << funame << "wrong uplo: " << uplo << std::endl;
    throw Error::Init();
  }
}

Lapack::Vector Lapack::SymmetricMatrix::operator* (const Vector& v) const throw(Error::General)
{
  const char funame [] = "Lapack::SymmetricMatrix::operator*: ";

  if(!isinit() || !v.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != v.size()) {
    std::cerr << funame << "the dimensions are different:"
	      << " vector_size = " << v.size()
	      << " matrix_size = " << size()
	      << "\n";
    throw Error::Range();
  }
  
  Vector res(size());
  dspmv_('U', size(),  1.,  *this, v, 1, 0., res, 1);
  return res;
}

Lapack::Vector Lapack::SymmetricMatrix::operator* (const double* v) const
{
  const char funame [] = "Lapack::SymmetricMatrix::operator*: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  Vector res(size());
  dspmv_('U', size(),  1., *this, v, 1, 0., res, 1);
  return res;
}

Lapack::Matrix Lapack::SymmetricMatrix::operator* (const Matrix& m)
  const throw(Error::General)
{
  const char funame [] = "Lapack::SymmetricMatrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != m.size1()) {
    std::cerr << funame << "the dimensions are different:"
	      << " mat1_size = " << size()
	      << " mat2_size1 = " << m.size1()
	      << std::endl;
    throw Error::Range();
  }

  Matrix res(size(), m.size2());
  for(int_t j = 0; j < m.size2(); ++j)
    dspmv_('U', size(),  1.,  *this, m + j * size(), 1, 0., res + j * size(), 1);  
  return res;
}

Lapack::Matrix Lapack::SymmetricMatrix::operator* (const SymmetricMatrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::SymmetricMatrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != m.size()) {
    std::cerr << funame << "the dimensions are different:"
	      << " mat1_size = " << size()
	      << " mat2_size = " << m.size()
	      << std::endl;
    throw Error::Range();
  }

  Matrix res(size());
  for(int_t i = 0; i < size(); ++i)
    for(int_t j = 0; j < size(); ++j) {
      double dtemp = 0.;
      const double* dp1 = *this + i * (i + 1) / 2;
      const double* dp2 = m     + j * (j + 1) / 2;

      for(int_t k = 0; k < size(); ++k) {
	dtemp += *dp1 * *dp2;
	// next first multiplier element
	if(k < i)
	  dp1++;
	else
	  dp1 += k + 1;
	// next second multiplier element
	if(k < j)
	  dp2++;
	else
	  dp2 += k + 1;
      }
      res(i,j) = dtemp;
    }
  return res;
}

Lapack::Vector Lapack::SymmetricMatrix::eigenvalues (Matrix* evec)
  const throw(Error::General)
{
  const char funame [] = "Lapack::SymmetricMatrix::eigenvalues: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  SymmetricMatrix sm = copy();
  Vector res(size());
  Vector work(3 * size());
  int_t info = 0;
  if(!evec) {
    dspev_('N', 'U', size(), sm, res, 0, size(), work, info);
  }
  else {
    evec->resize(size());
    dspev_('V', 'U', size(), sm, res, *evec, size(), work, info);
  }
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dspev: " << -info 
	      << "-th argument has an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dspev: " << info 
	      << "-th off-diagonal  elements  of an "
      "intermediate tridiagonal form did not converge to zero\n";
    throw Error::Math();
  }
}

Lapack::SymmetricMatrix Lapack::SymmetricMatrix::invert () const throw(Error::General)
{
  const char funame [] = "Lapack::SymmetricMatrix::invert: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  SymmetricMatrix lu = copy();

  Array<int_t> ipiv (size());
  Matrix res(size());
  res = 0.;
  res.diagonal() = 1.;

  int_t info = 0;
  dspsv_('U', size(), size(), lu, ipiv, res, size(), info);

  if(!info)
    return SymmetricMatrix(res);
  else if(info < 0) {
    std::cerr << funame << "dspsv: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dspsv: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe  factorization has  been  completed,  but the factor U is\n"
      "\texactly singular, so the solution could not be computed.\n";
    throw Error::Math();
  }
}

Lapack::SymmetricMatrix Lapack::SymmetricMatrix::positive_invert ()
  const throw(Error::General)
{
  const char funame [] = "Lapack::SymmetricMatrix::positive_invert: ";

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  SymmetricMatrix lu = copy();

  Matrix res(size());
  res = 0.;
  res.diagonal() = 1.;

  int_t info = 0;
  dppsv_('U', size(), size(), lu, res, size(), info);

  if(!info)
    return SymmetricMatrix(res);
  else if(info < 0) {
    std::cerr << funame << "dppsv: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dppsv: "
	      << "the leading minor of the " << info <<  "-th order of  A is\n"
      "\tnot positive definite, so the factorization could not be completed,\n"
      "\tand the solution has not been computed.\n";
    throw Error::Math();
  }
}

Lapack::Vector Lapack::diagonalize(SymmetricMatrix a0, SymmetricMatrix b0, Matrix* evec) throw(Error::General)
{
  const char funame [] = "Lapack::diagonalize: ";
  const char lapack_funame [] = "DSPGVD: ";

  if(!a0.isinit() || !b0.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(a0.size() != b0.size()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Range();
  }

  SymmetricMatrix a = a0.copy();
  SymmetricMatrix b = b0.copy();
  Vector res(a.size());

  char job = 'N';
  int_t lwork = 2 * a.size();
  int_t liwork = 1;
  double* z = 0;
  if(evec) {
    evec->resize(a.size());
    z = *evec;
    job = 'V';
    lwork = 1 + 6 * a.size() + 2 * a.size() * a.size();
    liwork = 3 + 5 * a.size();
  }

  Array<double> work(lwork);
  Array<int_t> iwork(liwork);
  int_t info = 0;

  dspgvd_(1, job, 'U', a.size(), a, b , res, z, a.size(), work, lwork, iwork, liwork, info);

  if(!info)
    return res;

  if(info < 0) {
    std::cerr << funame << lapack_funame << -info 
	      << "-th argument has an illegal value\n";
    throw Error::Range();
  }

  if(info <= a.size() ){
    std::cerr << funame << lapack_funame << info 
	      << " off-diagonal  elements  of an intermediate tridiagonal form did not converge to zero\n";
    throw Error::Math();
  }

  if(info <= 2 * a.size() ){
    std::cerr << funame << lapack_funame << " the leading minor of order " << info - a.size()
	      << " of B is not positively definite\n";
    throw Error::Range();
  }

  std::cerr << funame << lapack_funame <<  "unknown error\n";
  throw Error::Run();
}

Lapack::Vector Lapack::diagonalize (const HermitianMatrix& a0, const HermitianMatrix& b0, ComplexMatrix* evec) throw(Error::General)
{
  const char funame [] = "Lapack::diagonalize: ";

  const char lapack_funame [] = "ZHPGVD: ";

  if(!a0.size() || !b0.size()) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }

  if(a0.size() != b0.size()) {
    //
    std::cerr << funame << "dimensions mismatch\n";

    throw Error::Range();
  }

  HermitianMatrix a = a0;
  HermitianMatrix b = b0;

  Vector res(a.size());

  char   job = 'N';

  complex* z = 0;

  Array<complex> work(5);
  Array<double> rwork(5);
  Array<int_t>  iwork(5);

  int_t  lwork = -1;
  int_t lrwork = -1;
  int_t liwork = -1;


  if(evec) {
    //
    evec->resize(a.size());

    z = *evec;

    job = 'V';
  }

  int_t info = 1;

  zhpgvd_(1, job, 'U', a.size(), a, b , res, z, a.size(), work, lwork, rwork, lrwork, iwork, liwork, info);

  if(info < 0) {
    //
    std::cerr << funame << lapack_funame << -info 
	      << "-th argument has an illegal value\n";

    throw Error::Range();
  }

  if(info) {
    //
    std::cerr << funame << lapack_funame
	      << "failed with info = " << info << "\n";

    throw Error::Logic();
  }

  lwork  = (int_t)work[0].real();

  lrwork = (int_t)rwork[0];

  liwork = iwork[0];
  
  /*
  std::cout << funame 
	    << "  lwork = " << std::setw(6) <<  lwork 
	    << " lrwork = " << std::setw(6) << lrwork 
	    << " liwork = " << std::setw(6) << liwork 
	    << std::endl;
  */
    
  work.resize(lwork);

  rwork.resize(lrwork);

  iwork.resize(liwork);

  zhpgvd_(1, job, 'U', a.size(), a, b , res, z, a.size(), work, lwork, rwork, lrwork, iwork, liwork, info);

  if(!info)
    //
    return res;

  if(info < 0) {
    //
    std::cerr << funame << lapack_funame << -info 
	      << "-th argument has an illegal value\n";

    throw Error::Logic();
  }

  if(info <= a.size() ) {
    //
    std::cerr << funame << lapack_funame << info 
	      << " off-diagonal  elements  of an intermediate tridiagonal form did not converge to zero\n";

    throw Error::Math();
  }

  if(info <= 2 * a.size() ){
    //
    std::cerr << funame << lapack_funame << " the leading minor of order " << info - a.size()
	      << " of B is not positively definite\n";

    throw Error::Math();
  }

  std::cerr << funame << lapack_funame <<  "failed with info = " << info << "\n";

  throw Error::Run();
}

/****************************************************************
 ********************** LU Factorization ************************
 ****************************************************************/

// calculates a parity of the permutation
//
Lapack::int_t Lapack::parity (RefArr<int_t> perm)
{
  if(perm.size() < 2)
    return 0;

  int_t s = perm.size() - 1;
  int_t last = perm[s];

  for(int_t i = 0; i < s; ++i)
    if(perm[i] > last)
      --perm[i];

  perm.resize(s);
  return s - last + parity(perm);
}

double Lapack::LU::det () const
{
  const char funame [] = "Lapack::LU::det: ";

  double res = 1.;
  for(int_t i = 0; i < size(); ++i)
    res *= (*this)(i, i);

  int_t sign = 1;
  for(int_t i = 0; i < size(); ++i)
    if(_ipiv[i] != i + 1)
      sign = -sign;

  // std::cout << funame << "permutation sign = " << sign << "\n";

  if(sign > 0)
    return  res;
  else
    return -res;
}

Lapack::LU::LU (const Matrix& m) throw(Error::General)
  : Matrix(m.copy())
{
  const char funame [] = "Lapack::LU::LU: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  _ipiv.resize(m.size1());

  if(m.size1() != m.size2()) {
    std::cerr << funame << "not square\n";
    throw Error::Init();
  }

  int_t info = 0;
  dgetrf_(size(), size(), *this, size(), _ipiv, info);

  if(!info)
    return;
  if(info < 0) {
    std::cerr << funame << "dgetrf: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgetrf: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe  factorization has  been  completed,  but the factor U is\n"
      "\texactly singular, and division  by zero will occur if it is used\n"
      "\tto solve a system of equations.\n";
    throw Error::Math();
  }
}

Lapack::Matrix Lapack::LU::invert() const throw(Error::General)
{
  const char funame [] = "Lapack::LU::invert: ";

  Matrix res = copy();

  int_t info = 0;
  Array<double> work(size());
  dgetri_(size(), res, size(), _ipiv, work, size(), info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dgetri: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgetri: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe matrix is singular and its inverse could not be computed.\n";
    throw Error::Math();
  }
}

Lapack::Vector Lapack::LU::invert(const Vector& v) const throw(Error::General)
{
  const char funame [] = "Lapack::LU::invert: ";

  if(!v.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != v.size()) {
    std::cerr << funame << "dimensions mismatch:"
	      << " matrix size = " << size()
	      << " vector size = " << v.size()
	      << "\n";
    throw Error::Range();
  }

  Vector res = v.copy();

  int_t info = 0;
  dgetrs_('N', size(), 1, *this, size(), _ipiv, res, size(), info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dgetrs: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgetrs: unknown error code "<< info << std::endl;
    throw Error::Math();
  }
}

Lapack::Matrix Lapack::LU::invert(const Matrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::LU::invert: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != m.size1()) {
    std::cerr << funame << "dimensions mismatch:"
	      << " matrix size = " << size()
	      << " vector size = " << m.size()
	      << "\n";
    throw Error::Range();
  }

  Matrix res = m.copy();

  int_t info = 0;
  dgetrs_('N', size(), m.size2(), *this, size(), _ipiv, res, size(), info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dgetrs: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dgetrs: unknown error code "<< info << std::endl;
    throw Error::Math();
  }
}

/****************************************************************
 ******** LU factorization for symmetric packed matrix **********
 ****************************************************************/

Lapack::SymLU::SymLU (const SymmetricMatrix& m) throw(Error::General)
  : SymmetricMatrix(m.copy())
{
  const char funame [] = "Lapack::SymLU::SymLU: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  _ipiv.resize(m.size());

  int_t info = 0;
  dsptrf_('U', size(), *this, _ipiv, info);

  if(!info)
    return;
  if(info < 0) {
    std::cerr << funame << "dsptrf: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dsptrf: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe  factorization has  been  completed,  but the factor U is\n"
      "\texactly singular, and division  by zero will occur if it is used\n"
      "\tto solve a system of equations.\n";
    throw Error::Math();
  }
}

Lapack::SymmetricMatrix Lapack::SymLU::invert() const throw(Error::General)
{
  const char funame [] = "Lapack::SymLU::invert: ";

  SymmetricMatrix res = copy();

  int_t info = 0;
  Array<double> work(size());
  dsptri_('U', size(), res, _ipiv, work, info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dsptri: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dsptri: U("<< info <<"," << info 
	      << ") is exactly zero.\n"
      "\tThe matrix is singular and its inverse could not be computed.\n";
    throw Error::Math();
  }
}

Lapack::Vector Lapack::SymLU::invert(const Vector& v) const throw(Error::General)
{
  const char funame [] = "Lapack::SymLU::invert: ";

  if(!v.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != v.size()) {
    std::cerr << funame << "dimensions mismatch:"
	      << " matrix size = " << size()
	      << " vector size = " << v.size()
	      << "\n";
    throw Error::Range();
  }

  Vector res = v.copy();

  int_t info = 0;
  dsptrs_('U', size(), 1, *this, _ipiv, res, size(), info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dsptrs: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dsptrs: unknown error code "<< info << std::endl;
    throw Error::Math();
  }
}

Lapack::Matrix Lapack::SymLU::invert(const Matrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::SymLU::invert: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != m.size1()) {
    std::cerr << funame << "dimensions mismatch:"
	      << " matrix size = " << size()
	      << " rhs size = " << m.size1()
	      << "\n";
    throw Error::Range();
  }

  Matrix res = m.copy();

  int_t info = 0;
  dsptrs_('U', size(), m.size2(), *this, _ipiv, res, size(), info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dsptrs: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dsptrs: unknown error code "<< info << std::endl;
    throw Error::Math();
  }
}

/****************************************************************
 ******************* Cholesky Factorization *********************
 ****************************************************************/

Lapack::Cholesky::Cholesky (const SymmetricMatrix& m) throw(Error::General)
  : SymmetricMatrix(m.copy())
{
  const char funame [] = "Lapack::Cholesky::Cholesky: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  int_t info;
  dpptrf_('U', size(), *this, info);

  if(!info)
    return;
  else if(info < 0) {
    std::cerr << funame << "dpptrf: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dpptrf: the leading minor of the " 
	      << info <<  "-th order of  A is not\n"
      "\tpositive definite, and the factorization could not be completed.\n";
    throw Error::Math();
  }
}

Lapack::SymmetricMatrix Lapack::Cholesky::invert() const throw(Error::General)
{
  const char funame [] = "Lapack::Cholesky::invert: ";

  SymmetricMatrix res = copy();

  int_t info;
  dpptri_('U', size(), res, info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dpptri: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dpptri: U("<< info <<"," << info 
	      << ") is exactly zero and\n"
      "\tthe inverse could not be computed.\n";
    throw Error::Math();
  }
}

Lapack::Vector Lapack::Cholesky::invert(const Vector& v) const throw(Error::General)
{
  const char funame [] = "Lapack::Cholesky::invert: ";

  if(!v.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != v.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " matrix size = " << size()
	      << " vector size = " << v.size()
	      << "\n";
    throw Error::Range();
  }

  Vector res = v.copy();

  int_t info;
  dpptrs_('U', size(), 1, *this, res, size(), info);

  // error codes
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "dpptrs: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dpptrs: unknown error code "<< info << std::endl;
    throw Error::Math();
  }
}

Lapack::Matrix Lapack::Cholesky::invert(const Matrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::Cholesky::invert: ";

  if(!m.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(size() != m.size1()) {
    std::cerr << funame << "dimensions are different:"
	      << " matrix size = " << size()
	      << " rigt-hand size = " << m.size1()
	      << "\n";
    throw Error::Range();
  }

  Matrix res = m.copy();

  int_t info;
  dpptrs_('U', size(), res.size2(), *this, res, size(), info);

  // error codes
  if(!info)
    return res;

  if(info < 0) {
    std::cerr << funame << "dpptrs: " << -info 
	      << "-th argument had an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "dpptrs: unknown error code "<< info << std::endl;
    throw Error::Math();
  }
}

double Lapack::Cholesky::det_sqrt ()
{
  if(!size())
    return 0.;

  double* dp = *this;
  double res = *dp;
  for(int_t i = 2; i <= size(); ++i) {
    dp += i;
    res *= *dp;
  }
  return res;
}

/****************************************************************
 *********************** Complex Matrix *************************
 ****************************************************************/

Lapack::ComplexMatrix Lapack::ComplexMatrix::operator* (const ComplexMatrix& m) const throw(Error::General)
{
  const char funame [] = "Lapack::Matrix::operator*: ";

  if(size2() != m.size1()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Range();
  }

  if(!size1())
    return ComplexMatrix();
    
  ComplexMatrix res(size1(), m.size2());
  zgemm_('N', 'N', size1(), m.size2(), size2(), 1., *this, size1(), m, m.size1(), 0., res, size1());

  return res;
}

/****************************************************************
 ********************* Hermitian Matrix *************************
 ****************************************************************/

Lapack::HermitianMatrix::HermitianMatrix 
(const Lapack::ComplexMatrix& m, char uplo) throw(Error::General)
  : Array<complex>(m.size1() * (m.size1() + 1) / 2), _size(m.size1())
{
  const char funame [] = "Lapack::HermitianMatrix::HermitianMatrix: ";

  if(!m.size1()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(m.size1() != m.size2()) {
    std::cerr << funame << "dimensions mismatch:"
              << " size1 = " << m.size1()
              << " size2 = " << m.size2()
              << "\n";
    throw Error::Range();
  }

  switch(uplo) {
  case 'U':
    for(int_t i = 0; i < size(); ++i) {
      Slice<complex> res(*this + i * (i + 1) / 2, i + 1);
      ConstSlice<complex> col(m + i * size(), i + 1);
      res = col;
    }
    return;

  case 'L':
    for(int_t i = 0; i < size(); ++i) {
      Slice<complex> res(*this + i * (i + 1) / 2, i + 1);
      ConstSlice<complex> row(m + i, i + 1, size());
      res = row;
    }
    return;
  default:
    std::cerr << funame << "wrong uplo: " << uplo << std::endl;
    throw Error::Init();
  }
}

Lapack::Vector Lapack::HermitianMatrix::eigenvalues (ComplexMatrix* evec)
  const throw(Error::General)
{
  const char funame [] = "Lapack::HermitianMatrix::eigenvalues: ";

  if(!size()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  HermitianMatrix cp = *this;
  Vector res(size());
  Array<complex> work(2 * size() - 1);
  Array<double> rwork(3 * size() - 2);
  int_t info = 0;
  
  if(!evec) {
    zhpev_('N', 'U', size(), cp, res,     0, size(), work, rwork, info);
  }
  else {
    evec->resize(size());
    zhpev_('V', 'U', size(), cp, res, *evec, size(), work, rwork, info);
  }
  if(!info)
    return res;
  else if(info < 0) {
    std::cerr << funame << "zhpev: " << -info 
	      << "-th argument has an illegal value\n";
    throw Error::Range();
  }
  else {
    std::cerr << funame << "zhpev: " << info 
	      << "-th off-diagonal  elements  of an "
      "intermediate tridiagonal form did not converge to zero\n";
    throw Error::Math();
  }
}

/**********************************************************************************************
 ************************************** COMPLEX FOURIER TRANSFORM *****************************
 **********************************************************************************************/

Lapack::ComplexVector Lapack::fourier_transform (const Vector& fun, const MultiIndexConvert& mi)
{
  const char funame [] = "Lapack::fourier_transform: ";

  if(!fun.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Range();
  }
  
  if(mi.size() != fun.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

  double dtemp;
  int itemp;

  const double nfac = (double)mi.size();

  ComplexVector res(mi.size());

#pragma omp parallel for default(shared) private(itemp, dtemp)  schedule(static)

  for(int g = 0; g < mi.size(); ++g) {

    itemp = mi.conjugate(g);
    if(itemp < g)
      continue;

    std::vector<int> gv = mi(g);
    double re = 0.;
    double im = 0.;
    for(int h = 0; h < mi.size(); ++h) {
      std::vector<int> hv = mi(h);
      dtemp = 0.;
      for(int i = 0; i < mi.rank(); ++i)
	dtemp += double(gv[i] * hv[i]) / (double) mi.size(i);
      dtemp *= 2. * M_PI;
      re += std::cos(dtemp) * fun[h];
      im += std::sin(dtemp) * fun[h];
    }
    re /= nfac;
    im /= nfac;

    if(g != itemp)
      res[g] = complex(re, -im);
    else 
      res[g] = re;
  }

#pragma omp parallel for default(shared) private(itemp, dtemp)  schedule(static)

  for(int g = 0; g < mi.size(); ++g) {
    itemp = mi.conjugate(g);
    if(itemp < g)
      res[g] = std::conj(res[itemp]);
  }

  return res;
}

Lapack::ComplexVector Lapack::fourier_transform (const ComplexVector& fun, const MultiIndexConvert& mi)
{
  const char funame [] = "Lapack::fourier_transform: ";

  if(!fun.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Range();
  }
  
  if(mi.size() != fun.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

  ComplexVector res(mi.size());

#pragma omp parallel for default(shared) schedule(static)

  for(int g = 0; g < mi.size(); ++g) {

    double dtemp;
    int itemp;

    std::vector<int> gv = mi(g);
    complex cval = 0.;
    for(int h = 0; h < mi.size(); ++h) {
      std::vector<int> hv = mi(h);
      dtemp = 0.;
      for(int i = 0; i < mi.rank(); ++i)
	dtemp += double(gv[i] * hv[i]) / (double) mi.size(i);
      dtemp *= 2. * M_PI;
      cval += complex(std::cos(dtemp), std::sin(dtemp)) * fun[h];
    }
    res[g] = cval;
  }

  return res;
}

Lapack::ComplexVector Lapack::fourier_transform (const std::map<int, complex>& fun, const MultiIndexConvert& mi)
{
  const char funame [] = "Lapack::fourier_transform: ";

  ComplexVector res(mi.size());

#pragma omp parallel for default(shared)  schedule(static)

  for(int g = 0; g < mi.size(); ++g) {
    std::vector<int> gv = mi(g);

    double dtemp;
    complex cval = 0.;
    for(std::map<int, complex>::const_iterator it = fun.begin(); it != fun.end(); ++it) {
      std::vector<int> hv;

      try {
	hv = mi(it->first);
      }
      catch(Error::General) {
	std::cerr << funame << "index out of range\n";
	throw;
      }

      dtemp = 0.;
      for(int i = 0; i < mi.rank(); ++i)
	dtemp += double(gv[i] * hv[i]) / (double) mi.size(i);
      dtemp *= 2. * M_PI;
      cval += complex(std::cos(dtemp), std::sin(dtemp)) * it->second;
    }
    res[g] = cval;
  }

  return res;
}
