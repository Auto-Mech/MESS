/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2021, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/


#include "mplapack_dd.hh"

#ifdef WITH_MPLAPACK

#include "error.hh"
#include "io.hh"

//#include <mplapack/mpblas_dd.h>
#include <mplapack/mplapack_dd.h>

#include "linpack.hh"

#include <iostream>

#pragma omp declare reduction (+: dd_real: omp_out += omp_in) initializer(omp_priv = 0.)
	
/****************************************************************
 ************************** Vector ******************************
 ****************************************************************/

dd_real Mpack_dd::Vector::operator* (Vector v) const 
{
  const char funame [] = "Mpack_dd::Vector::operator*: ";

  if(!isinit() || !v.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", "<< v.size() << "\n";
    
    throw Error::Range();
  }

  return *this * (const dd_real*)v;
}

dd_real Mpack_dd::Vector::operator* (ConstSlice<dd_real> v) const 
{
  const char funame [] = "Mpack_dd::Vector::operator*: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", "<< v.size() << "\n";
    
    throw Error::Range();
  }
  
  dd_real res = 0.;
  
#pragma omp parallel for default(shared) reduction(+: res) schedule(static)
  //
  for(int i = 0; i < size(); ++i)
    //
    res += (*this)[i] * v[i];
  
  return res;
}

dd_real Mpack_dd::Vector::operator* (const dd_real* v) const 
{
  const char funame [] = "Mpack_dd::Vector::operator*: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  dd_real res = 0.;
  
#pragma omp parallel for default(shared) reduction(+: res) schedule(static)
  //
  for(int i = 0; i < size(); ++i)
    //
    res += (*this)[i] * v[i];
  
  return res;
}

dd_real Mpack_dd::Vector::vdot () const
{
  const char funame [] = "Mpack_dd::Vector::vdot: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  return *this * (const dd_real*)(*this);
}

Mpack_dd::Vector Mpack_dd::Vector::operator* (const Matrix& m) const
{
  const char funame [] = "Mpack_dd::Vector::operator*: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }
  
  if(size() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size1() << "\n";
    
    throw Error::Range();
  }
  
  //Rgemv("T", size(), m.size2(), 1.,  m, size(), *this, 1, 0., res, 1);

  return (const dd_real*)(*this) * m;
}

Mpack_dd::Vector Mpack_dd::Vector::operator* (const SymmetricMatrix& m) const
{
  const char funame [] = "Mpack_dd::Vector::operator*: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }
  
  if(size() != m.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size() << "\n";
    
    throw Error::Range();
  }
  
  //Rspmv("U", size(), 1., m,  *this, 1, 0., res, 1);

  return (const dd_real*)(*this) * m;
}

/****************************************************************
 ************************** Matrix ******************************
 ****************************************************************/

// copy by value
//
Mpack_dd::Matrix Mpack_dd::Matrix::copy () const
{
  const char funame [] = "Mpack_dd::Matrix::copy: ";

  _assert();

  Matrix res(size1(), size2());

#pragma omp parallel for default(shared)  schedule(static)
  //
  for(int i = 0; i < RefArr<dd_real>::size(); ++i)
    //
    res._at(i) = _at(i);

  return res;
}
  
Mpack_dd::Matrix::Matrix (const SymmetricMatrix& m)
{
  const char funame [] = "Mpack_dd::Matrix::Matrix: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }

  resize(m.size());

  const int n = m.size() * m.size();

#pragma omp parallel for default(shared) schedule(static)

  for(int k = 0; k < n; ++k) {
    //
    int i = k / m.size();

    int j = k % m.size();

    (*this)(i, j) = m(i, j);
  }
}

Mpack_dd::Vector Mpack_dd::Matrix::eigenvalues (Matrix* evec) const
{
  const char funame [] = "Mpack_dd::Matrix::eigenvalues: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }

  if(size1() != size2()) {
    //
    std::cerr << funame << "not square: "<< size1() << ", " << size2() << "\n";

    throw Error::Init();
  }

  const char* jobz = "N";
  //
  if(evec)
    //
    jobz = "V";

  Matrix sm = copy();

  Vector res(size());

  int_t lwork, info = 0;

  //work space query
  //
  lwork = -1;

  Array<dd_real> work(1);

  Rsyev(jobz, "U", size(), sm, size(), res, work, lwork, info);

  // matrix diagonalization
  //
  lwork = std::max((int_t)1, (int_t)work[0]._hi());

  work.resize((int)lwork);

  Rsyev(jobz, "U", size(), sm, size(), res, work, lwork, info);
    
  if(evec)
    //
    *evec = sm;

  if(!info)
    //
    return res;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rsyev: " << -info << "-th argument has illegal value\n";

    throw Error::Logic();
  }
  
  std::cerr << funame << "Rsyev: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero\n";

  throw Error::Lapack();
}

Mpack_dd::Matrix Mpack_dd::Matrix::transpose () const
{
  const char funame [] = "Mpack_dd::Matrix::transpose: ";
  
  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }
  
  Matrix res(size2(), size1());

  const int n = size1() * size2();

#pragma omp parallel for default(shared)  schedule(static)

  for(int k = 0; k < n; ++k) {
    //
    int i = k / size2();

    int j = k % size2();

    res(j, i) = (*this)(i, j);
  }
  
  return res;
}

Mpack_dd::Matrix Mpack_dd::Matrix::operator* (const Matrix& m) const
{
  const char funame [] = "Mpack_dd::Matrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size2() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch\n";
    
    throw Error::Range();
  }

  Matrix res(size1(), m.size2());
  
  const int n = size1() * m.size2();

#pragma omp parallel for default(shared)  schedule(static)

  for(int k = 0; k < n; ++k) {
    //
    int i = k / m.size2();

    int j = k % m.size2();

    res(i, j) = vdot(row(i), m.column(j));
  }
  
  //Rgemm("N", "N", size1(), m.size2(), size2(), 1., *this, size1(), m, m.size1(), 0., res, size1());

  return res;
}

Mpack_dd::Matrix Mpack_dd::Matrix::operator* (const SymmetricMatrix& m) const
{
  const char funame [] = "Mpack_dd::Matrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }
  
  if(size2() != m.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size2() << ", " << m.size() << "\n";
    
    throw Error::Range();
  }

  Matrix res(size1(), size2());
  
  const int n = size1() * size2();

#pragma omp parallel for default(shared)  schedule(static)

  for(int k = 0; k < n; ++k) {
    //
    int i = k / size2();

    int j = k % size2();

    dd_real dtemp = 0.;

    for(int l = 0; l < size2(); ++l)
      //
      dtemp += (*this)(i, l) * m(l, j);
    
    res(i, j) = dtemp;
  }
  
  //Rspmv("U", size2(),  1., sm, m + i, size1(), 0., res + i, size1());

  return res;
}

Mpack_dd::Vector Mpack_dd::Matrix::operator* (const dd_real* v) const
{
  const char funame [] = "Mpack_dd::Matrix::operator*: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  Vector res(size1());
  
#pragma omp parallel for default(shared)  schedule(static)

  for(int i = 0; i < size1(); ++i)
    //
    res[i] = vdot(row(i), v);
    
  //Rgemv("N", size1(), size2(), 1.,  *this, size1(), v, 1, 0., res, 1);
  
  return res;
}

Mpack_dd::Vector Mpack_dd::operator* (const dd_real* v, const Matrix& m)
{
  const char funame [] = "Mpack_dd::operator*: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  Vector res(m.size2());
  
#pragma omp parallel for default(shared)  schedule(static)

  for(int i = 0; i < m.size2(); ++i)
    //
    res[i] = vdot(m.column(i), v);
    
  //Rgemv("T", m.size1(), m.size2(), 1., m, m.size1(), v, 1, 0., res, 1);
  
  return res;
}

// inverse matrix
//
Mpack_dd::Matrix Mpack_dd::Matrix::invert () const
{
  const char funame [] = "Mpack_dd::Matrix::invert: ";

  _assert();

  if(size1() != size2()) {
    //
    std::cerr << funame << "not square\n";
    
    throw Error::Logic();
  }

  Matrix lu = copy();
  
  Array<int_t> ipiv ((int)size1());
  
  Matrix res(size1(), size1());
  
  res = 0.;
  
  res.diagonal() = (dd_real)1.;

  int_t info = 0;
  
  Rgesv(size1(), size1(), lu, size1(), ipiv, res, size1(), info);

 if(!info)
   //
   return res;

 // error codes
 //
 if(info < 0) {
   //
    std::cerr << funame << "Rgesv: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
 }

 std::cerr << funame << "Rgesv: U("<< info <<"," << info << ") is exactly zero.\n";
 
 throw Error::Run();
}

// solve linear equations
//
Mpack_dd::Vector Mpack_dd::Matrix::invert (Vector v) const
{
  const char funame [] = "Mpack_dd::Matrix::invert: ";

  if(!isinit() || !v.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << v.size() << "\n";
    
    throw Error::Logic();
  }

  Matrix lu = copy();
  
  Array<int_t> ipiv((int)size1());
  
  Vector res = v.copy();

  int_t info = 0;
  
  Rgesv(size1(), 1, lu, size1(), ipiv, res, size1(), info);

  if(!info)
    //
    return res;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rgesv: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rgesv: U("<< info <<"," << info << ") is exactly zero.\n";

  throw Error::Run();
}

// solve linear equations
//
Mpack_dd::Matrix Mpack_dd::Matrix::invert (Matrix m) const
{
  const char funame [] = "Mpack_dd::Matrix::invert: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size1() << "\n";
    
    throw Error::Logic();
  }

  Matrix lu = copy();
  
  Array<int_t> ipiv((int)size1());
  
  Matrix res = m.copy();

  int_t info = 0;
  
  Rgesv(size1(), m.size2(), lu, size1(), ipiv, res, size1(), info);

 if(!info)
   //
   return res;

 // error codes
 //
 if(info < 0) {
   //
   std::cerr << funame << "Rgesv: " << -info << "-th argument had an illegal value\n";
   
   throw Error::Range();
 }

 std::cerr << funame << "Rgesv: U("<< info <<"," << info << ") is exactly zero.\n";

 throw Error::Math();
}

// matrix row
//
Slice<dd_real> Mpack_dd::Matrix::row (int_t i) 
{
  const char funame [] = "Mpack_dd::Matrix::row: ";

  _assert();

  if(i < 0 || i >= size1()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";
    
    throw Error::Range();
  }

  return Slice<dd_real>(*this + i, size2(), size1());
}

ConstSlice<dd_real> Mpack_dd::Matrix::row (int_t i) const
{
  const char funame [] = "Mpack_dd::Matrix::row: ";

  _assert();

  if(i < 0 || i >= size1()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";
    
    throw Error::Range();
  }

  return ConstSlice<dd_real>(*this + i, size2(), size1());
}

// matrix column
//
Slice<dd_real> Mpack_dd::Matrix::column (int_t i) 
{
  const char funame [] = "Mpack_dd::Matrix::column: ";

  _assert();

  if(i < 0 || i >= size2()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";
    
    throw Error::Range();
  }

  return Slice<dd_real>(*this + i * size1(), size1());
}

ConstSlice<dd_real> Mpack_dd::Matrix::column (int_t i) const
{
  const char funame [] = "Mpack_dd::Matrix::column: ";

  _assert();

  if(i < 0 || i >= size2()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";
    
    throw Error::Range();
  }

  return ConstSlice<dd_real>(*this + i * size1(), size1());
}

// matrix diagonal
//
Slice<dd_real> Mpack_dd::Matrix::diagonal (int_t i) 
{
  const char funame [] = "Mpack_dd::Matrix::diagonal: ";

  _assert();

  if(i <= -size1() || i >= size2()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";
    
    throw Error::Range();
  }

  int_t sz;
  
  if(i > 0) {// upper diagonal
    //
    sz = size1() < size2() - i ? size1() : size2() - i;
    
    return Slice<dd_real>(*this + i * size1(), sz, size1() + 1);
  }
  else {// lower diagonal
    //
    sz = size1() + i < size2() ? size1() + i : size2();
    
    return Slice<dd_real>(*this - i, sz, size1() + 1);
  }
}

ConstSlice<dd_real> Mpack_dd::Matrix::diagonal (int_t i) const
{
  const char funame [] = "Mpack_dd::Matrix::diagonal: ";

  _assert();

  if(i <= -size1() || i >= size2()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";
    
    throw Error::Range();
  }

  int_t sz;
  
  if(i > 0) {// upper diagonal
    //
    sz = size1() < size2() - i ? size1() : size2() - i;
    
    return ConstSlice<dd_real>(*this + i * size1(), sz, size1() + 1);
  }
  else {// lower diagonal
    //
    sz = size1() + i < size2() ? size1() + i : size2();
    
    return ConstSlice<dd_real>(*this - i, sz, size1() + 1);
  }
}

/****************************************************************
 ********************* Symmetric Matrix *************************
 ****************************************************************/

Mpack_dd::SymmetricMatrix Mpack_dd::SymmetricMatrix::copy () const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::copy: ";

  _assert();

  SymmetricMatrix res(size());

#pragma omp parallel for default(shared)  schedule(static)
  //
  for(int i = 0; i < RefArr<dd_real>::size(); ++i)
    //
    res._at(i) = _at(i);

  return res;
}
  
Mpack_dd::SymmetricMatrix::SymmetricMatrix (Matrix m, char uplo)
  //
  : RefArr<dd_real>(m.size() * (m.size() + 1) / 2), _size(new int_t(m.size()))
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::SymmetricMatrix: ";

  const int n = size() * size();
  
#pragma omp parallel for default(shared)  schedule(dynamic)
  //
  for(int k = 0; k < n; ++k) {
    //
    int i = k / size();

    int j = k % size();

    if(i > j)
      //
      continue;
      

    switch(uplo) {
      //
    case 'U':
      //
      (*this)(i, j) = m(i, j);

      break;
      //
    case 'L':
      //
      (*this)(i, j) = m(j, i);

      break;
      //
    default:
      //
      std::cerr << funame << "wrong uplo: " << uplo << std::endl;
    
      throw Error::Init();
    }
  }
}

Mpack_dd::Vector Mpack_dd::SymmetricMatrix::operator* (Vector v) const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator*: ";

  if(!isinit() || !v.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << v.size() << "\n";
    
    throw Error::Range();
  }

  return *this * (const dd_real*)v;
}

Mpack_dd::Vector Mpack_dd::SymmetricMatrix::operator* (ConstSlice<dd_real> v) const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator*: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << v.size() << "\n";
    
    throw Error::Range();
  }
  
  Vector res(size());
  
#pragma omp parallel for default(shared)  schedule(static)

  for(int i = 0; i < size(); ++i) {
    //
    dd_real dtemp = 0.;

    for(int j = 0; j < size(); ++j)
      //
      dtemp += (*this)(i, j) * v[j];

    res[i] = dtemp;
  }
  
  //Rspmv("U", size(),  1.,  *this, v.begin(), v.stride(), 0., res, 1);
  
  return res;
}

Mpack_dd::Vector Mpack_dd::SymmetricMatrix::operator* (const dd_real* v) const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator*: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  Vector res(size());

#pragma omp parallel for default(shared)  schedule(static)

  for(int i = 0; i < size(); ++i) {
    //
    dd_real dtemp = 0.;

    for(int j = 0; j < size(); ++j)
      //
      dtemp += (*this)(i, j) * v[j];

    res[i] = dtemp;
  }
  
  //Rspmv("U", size(),  1., *this, v, 1, 0., res, 1);
  
  return res;
}

Mpack_dd::Matrix Mpack_dd::SymmetricMatrix::operator* (const Matrix& m) const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size1() << "\n";

    throw Error::Range();
  }

  Matrix res(m.size1(), m.size2());
  
  const int n = m.size1() * m.size2();

#pragma omp parallel for default(shared)  schedule(static)

  for(int k = 0; k < n; ++k) {
    //
    int i = k / m.size2();

    int j = k % m.size2();

    dd_real dtemp = 0.;

    for(int l = 0; l < size(); ++l)
      //
      dtemp += (*this)(i, l) * m(l, j);
    
    res(i, j) = dtemp;
  }
  
  //Rspmv("U", size(),  1.,  *this, m + j * size(), 1, 0., res + j * size(), 1);
  
  return res;
}

Mpack_dd::Matrix Mpack_dd::SymmetricMatrix::operator* (const SymmetricMatrix& m) const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator*: ";

  if(!isinit() || !m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != m.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size() << "\n";
    
    throw Error::Range();
  }

  Matrix res(size());
  
  const int n = size() * size();

#pragma omp parallel for default(shared)  schedule(static)
  //
  for(int k = 0; k < n; ++k) {
    //
    int i = k / size();

    int j = k % size();

    dd_real dtemp = 0.;

    for(int l = 0; l < size(); ++l)
      //
      dtemp += (*this)(i, l) * m(l, j);
    
    res(i, j) = dtemp;
  }
  
  return res;
}

Mpack_dd::SymmetricMatrix Mpack_dd::SymmetricMatrix::operator= (dd_real d)
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator=: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  RefArr<dd_real>::operator=(0.);

#pragma omp parallel for default(shared)  schedule(static)
  //
  for(int i = 0; i < size(); ++i)
    //
    (*this)(i, i) = d;
  
  return *this;
}

Mpack_dd::SymmetricMatrix Mpack_dd::SymmetricMatrix::operator+= (dd_real d)
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator+=: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
      
    throw Error::Init();
  }

#pragma omp parallel for default(shared)  schedule(static)
  //
  for(int i = 0; i < size(); ++i)
    //
    (*this)(i, i) += d;
  
  return *this;
}

Mpack_dd::SymmetricMatrix Mpack_dd::SymmetricMatrix::operator-= (dd_real d)
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::operator-=: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
      
    throw Error::Init();
  }

#pragma omp parallel for default(shared)  schedule(static)
  //
  for(int i = 0; i < size(); ++i)
    //
    (*this)(i, i) -= d;
  
  return *this;
}

Mpack_dd::Vector Mpack_dd::SymmetricMatrix::eigenvalues (Matrix* evec) const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::eigenvalues: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }

  const char* jobz = "N";
  //
  if(evec)
    //
    jobz = "V";

  Matrix sm(*this);

  Vector res(size());

  int_t lwork, info = 0;

  //work space query
  //
  lwork = -1;

  Array<dd_real> work(1);

  Rsyev(jobz, "U", size(), sm, size(), res, work, lwork, info);

  // matrix diagonalization
  //
  lwork = std::max((int_t)1, (int_t)work[0]._hi());

  work.resize(lwork);

  Rsyev(jobz, "U", size(), sm, size(), res, work, lwork, info);

  if(evec)
    //
    *evec = sm;

  if(!info)
    //
    return res;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rsyev: " << -info << "-th argument has an illegal value\n";

    throw Error::Logic();
  }
  
  std::cerr << funame << "Rsyev: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero\n";
    
  throw Error::Run();
}

Mpack_dd::SymmetricMatrix Mpack_dd::SymmetricMatrix::invert () const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::invert: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  SymmetricMatrix lu = copy();

  Array<int_t> ipiv ((int)size());
  
  Matrix res(size());
  
  res = 0.;
  
  res.diagonal() = (dd_real)1.;

  int_t info = 0;
  
  Rspsv("U", size(), size(), lu, ipiv, res, size(), info);

  if(!info)
    //
    return SymmetricMatrix(res);

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "dspsv: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "dspsv: U("<< info <<"," << info << ") is exactly zero\n";
    
  throw Error::Run();
}

Mpack_dd::SymmetricMatrix Mpack_dd::SymmetricMatrix::positive_invert () const
{
  const char funame [] = "Mpack_dd::SymmetricMatrix::positive_invert: ";

  if(!isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  SymmetricMatrix lu = copy();

  Matrix res(size());
  
  res = 0.;
  res.diagonal() = (dd_real)1.;

  int_t info = 0;
  
  Rppsv("U", size(), size(), lu, res, size(), info);

  if(!info)
    //
    return SymmetricMatrix(res);

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rppsv: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }
  
  std::cerr << funame << "Rppsv: " << "the leading minor of the " << info <<  "-th order of  A is not positive definite\n";
  
  throw Error::Run();
}

Mpack_dd::Vector Mpack_dd::diagonalize(SymmetricMatrix a0, SymmetricMatrix b0, Matrix* evec) 
{
  const char funame [] = "Mpack_dd::diagonalize: ";
  
  const char lapack_funame [] = "Rspgvd: ";

  if(!a0.isinit() || !b0.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(a0.size() != b0.size()) {
    //
    std::cerr << funame << "dimensions mismatch\n";
    
    throw Error::Range();
  }

  SymmetricMatrix a = a0.copy();
  
  SymmetricMatrix b = b0.copy();
  
  Vector res(a.size());

  const char* job = "N";
  
  int_t lwork = 2 * a.size();
  
  int_t liwork = 1;
  
  dd_real* z = 0;
  
  if(evec) {
    //
    evec->resize(a.size());
    
    z = *evec;
    
    job = "V";
    
    lwork = 1 + 6 * a.size() + 2 * a.size() * a.size();
    
    liwork = 3 + 5 * a.size();
  }

  Array<dd_real> work((int)lwork);
  
  Array<int_t> iwork((int)liwork);
  
  int_t info = 0;

  Rspgvd(1, job, "U", a.size(), a, b , res, z, a.size(), work, lwork, iwork, liwork, info);

  if(!info)
    //
    return res;

  // error codes
  //
  if(info < 0) {
    std::cerr << funame << lapack_funame << -info << "-th argument has an illegal value\n";
    
    throw Error::Range();
  }

  if(info <= a.size() ){
    //
    std::cerr << funame << lapack_funame << info
      //
	      << " off-diagonal  elements  of an intermediate tridiagonal form did not converge to zero\n";
    
    throw Error::Run();
  }

  if(info <= 2 * a.size() ){
    //
    std::cerr << funame << lapack_funame << " the leading minor of order " << info - a.size()
      //
	      << " of B is not positively definite\n";
    
    throw Error::Run();
  }

  std::cerr << funame << lapack_funame <<  "unknown error\n";
  
  throw Error::Run();
}

/****************************************************************
 ******************* Band Symmetric Matrix **********************
 ****************************************************************/

std::pair<Mpack_dd::int_t, Mpack_dd::int_t> Mpack_dd::BandMatrix::_convert (int_t i, int_t j) const
{
  const char funame [] = "Mpack_dd::BandMatrix::operator(): ";

  _assert();
  
  if(i > j)
    //
    std::swap(i, j);
 
  if(j - i < band_size())
    //
    return std::make_pair(band_size() - 1 + i - j, j);

  std::cerr << funame << "out of range\n";
  
  throw Error::Range();

}

Mpack_dd::Vector Mpack_dd::BandMatrix::eigenvalues (Matrix* evec) const
{
  const char funame [] = "Mpack_dd::BandMatrix::eigenvalues: ";

  _assert();

  const char* job;
  
  dd_real* z;
  
  int_t lwork = -1;
  
  int_t liwork = -1;
  
  if(!evec) {
    //
    job = "N";
    
    z = 0;
    
    lwork = 2 * size();
    
    liwork = 1;
  }
  else {
    //
    job = "V";
    
    evec->resize(size());
    
    z = *evec;
    
    lwork = 1 + 5 * size() + 2 * size() * size();
    
    liwork = 3 + 5 * size();
  }

  Vector work(lwork);
  
  Array<int_t> iwork((int)liwork);

  Matrix a = copy();
  
  Vector res(size());
  
  int_t info = 0;
  
  Rsbevd(job, "U", size(), band_size() - 1, a, band_size(), res, z, size(), work, lwork, iwork, liwork, info);
 
  if(!info)
    //
    return res;

  if(info < 0) {
    //
    std::cerr << funame << "Rsbevd: " << -info << "-th argument has an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rsbevd: " << info << "-th off-diagonal  elements  of an  intermediate tridiagonal form did not converge to zero\n";
  
  throw Error::Run();
}

/****************************************************************
 ********************** LU Factorization ************************
 ****************************************************************/

dd_real Mpack_dd::LU::det () const
{
  const char funame [] = "Mpack_dd::LU::det: ";

  dd_real res = 1.;
  
  for(int_t i = 0; i < size(); ++i)
    //
    res *= (*this)(i, i);

  int_t sign = 1;
  
  for(int i = 0; i < size(); ++i)
    //
    if(_ipiv[i] != i + 1)
      //
      sign = -sign;

  // std::cout << funame << "permutation sign = " << sign << "\n";

  if(sign > 0) {
    //
    return  res;
  }
  else
    //
    return -res;
}

Mpack_dd::LU::LU (Matrix m) : Matrix(m.copy())
{
  const char funame [] = "Mpack_dd::LU::LU: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  _ipiv.resize(m.size());

  int_t info = 0;
  
  Rgetrf(size(), size(), *this, size(), _ipiv, info);

  if(!info)
    //
    return;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rgetrf: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }
  
  std::cerr << funame << "Rgetrf: U("<< info <<"," << info << ") is exactly zero.\n";
  
  throw Error::Run();
}

Mpack_dd::Matrix Mpack_dd::LU::invert() const
{
  const char funame [] = "Mpack_dd::LU::invert: ";

  Matrix res = copy();

  int_t info = 0;
  
  Array<dd_real> work((int)size());
  
  Rgetri(size(), res, size(), _ipiv, work, size(), info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rgetri: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }
  
  std::cerr << funame << "Rgetri: U("<< info <<"," << info << ") is exactly zero.\n";

  throw Error::Run();
}

Mpack_dd::Vector Mpack_dd::LU::invert(Vector v) const
{
  const char funame [] = "Mpack_dd::LU::invert: ";

  if(!v.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << v.size() << "\n";
    
    throw Error::Range();
  }

  Matrix a = copy();
  
  Vector res = v.copy();

  int_t info = 0;
  
  Rgetrs("N", size(), 1, a, size(), _ipiv, res, size(), info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rgetrs: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }
  
  std::cerr << funame << "Rgetrs: unknown error code "<< info << "\n";
  
  throw Error::Run();
}

Mpack_dd::Matrix Mpack_dd::LU::invert(Matrix m) const
{
  const char funame [] = "Mpack_dd::LU::invert: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size() << "\n";
    
    throw Error::Range();
  }

  Matrix a = copy();
  
  Matrix res = m.copy();

  int_t info = 0;
  
  Rgetrs("N", size(), m.size2(), a, size(), _ipiv, res, size(), info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rgetrs: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }
  
  std::cerr << funame << "Rgetrs: unknown error code "<< info << "\n";
  
  throw Error::Math();
}

/****************************************************************
 ******** LU factorization for symmetric packed matrix **********
 ****************************************************************/

Mpack_dd::SymLU::SymLU (SymmetricMatrix m) : SymmetricMatrix(m.copy())
{
  const char funame [] = "Mpack_dd::SymLU::SymLU: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  _ipiv.resize(m.size());

  int_t info = 0;
  
  Rsptrf("U", size(), *this, _ipiv, info);

  if(!info)
    //
    return;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rsptrf: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rsptrf: U("<< info <<"," << info << ") is exactly zero.\n";
  
  throw Error::Run();
}

dd_real Mpack_dd::SymLU::det () const
{
  const char funame [] = "Mpack_dd::SymLU::det: ";

  dd_real res = 1.;
  
  for(int i = 0; i < size(); ++i)
    //
    if(_ipiv[i] > 0) {
      //
      res *= (*this)(i, i);
    }
    else if(i && _ipiv[i] < 0 && _ipiv[i] == _ipiv[i - 1]) {
      //
      res *= (*this)(i, i) * (*this)(i - 1, i - 1) - (*this)(i - 1, i) * (*this)(i - 1, i);
    }

  return res;
}

Mpack_dd::SymmetricMatrix Mpack_dd::SymLU::invert() const
{
  const char funame [] = "Mpack_dd::SymLU::invert: ";

  SymmetricMatrix res = copy();

  int_t info = 0;
  
  Array<dd_real> work((int)size());
  
  Rsptri("U", size(), res, _ipiv, work, info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rsptri: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }
  
  std::cerr << funame << "Rsptri: U("<< info <<"," << info << ") is exactly zero.\n";

  throw Error::Run();
}

Mpack_dd::Vector Mpack_dd::SymLU::invert(Vector v) const
{
  const char funame [] = "Mpack_dd::SymLU::invert: ";

  if(!v.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << v.size() << "\n";
    
    throw Error::Range();
  }

  SymmetricMatrix a = copy();
  
  Vector res = v.copy();

  int_t info = 0;
  
  Rsptrs("U", size(), 1, a, _ipiv, res, size(), info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rsptrs: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rsptrs: unknown error code "<< info << "\n";

  throw Error::Run();
}

Mpack_dd::Matrix Mpack_dd::SymLU::invert(Matrix m) const
{
  const char funame [] = "Mpack_dd::SymLU::invert: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size1() << "\n";
    
    throw Error::Range();
  }

  SymmetricMatrix a = copy();
  
  Matrix res = m.copy();

  int_t info = 0;
  
  Rsptrs("U", size(), m.size2(), a, _ipiv, res, size(), info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rsptrs: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rsptrs: unknown error code "<< info << "\n";

  throw Error::Run();
}

/****************************************************************
 ******************* Cholesky Factorization *********************
 ****************************************************************/

Mpack_dd::Cholesky::Cholesky (SymmetricMatrix m)
  //
  : SymmetricMatrix(m.copy())
{
  const char funame [] = "Mpack_dd::Cholesky::Cholesky: ";

  int_t info;
  
  Rpptrf("U", size(), *this, info);

  if(!info)
    //
    return;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rpptrf: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rpptrf: the leading minor of the " << info <<  "-th order of  A is not\n"
    //
	    << "\tpositive definite, and the factorization could not be completed.\n";
  
  throw Error::Run();
}

Mpack_dd::SymmetricMatrix Mpack_dd::Cholesky::invert() const
{
  const char funame [] = "Mpack_dd::Cholesky::invert: ";

  SymmetricMatrix res = copy();

  int_t info;

  Rpptri("U", size(), res, info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rpptri: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rpptri: U("<< info <<"," << info << ") is exactly zero\n";

  throw Error::Math();
}

Mpack_dd::Vector Mpack_dd::Cholesky::invert(Vector v) const
{
  const char funame [] = "Mpack_dd::Cholesky::invert: ";

  if(!v.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != v.size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << v.size() << "\n";
    
    throw Error::Range();
  }

  SymmetricMatrix a = copy();
  
  Vector res = v.copy();

  int_t info;
  
  Rpptrs("U", size(), 1, a, res, size(), info);

  if(!info)
    //
    return res;
  
  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rpptrs: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "Rpptrs: unknown error code "<< info << std::endl;

  throw Error::Run();
}

Mpack_dd::Matrix Mpack_dd::Cholesky::invert(Matrix m) const
{
  const char funame [] = "Mpack_dd::Cholesky::invert: ";

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  if(size() != m.size1()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size() << ", " << m.size1() << "\n";
    
    throw Error::Range();
  }

  SymmetricMatrix a = copy();
  
  Matrix res = m.copy();

  int_t info;
  
  Rpptrs("U", size(), res.size2(), a, res, size(), info);

  if(!info)
    //
    return res;

  // error codes
  //
  if(info < 0) {
    //
    std::cerr << funame << "Rpptrs: " << -info << "-th argument had an illegal value\n";
    
    throw Error::Range();
  }

  std::cerr << funame << "dpptrs: unknown error code "<< info << "\n";

  throw Error::Math();
}

dd_real Mpack_dd::Cholesky::det_sqrt () const
{
  if(!size())
    //
    return 0.;

  const dd_real* dp = *this;
  
  dd_real res = *dp;
  
  for(int_t i = 2; i <= size(); ++i) {
    //
    dp += i;
    
    res *= *dp;
  }
  
  return res;
}

#endif
