/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef LAPACK_HH
#define LAPACK_HH

#include "array.hh"
#include "shared.hh"
#include "multindex.hh"

#include <complex>
#include <map>

namespace Lapack {

  typedef int int_t;// index type in lapack & blas libraries

  typedef std::complex<double> complex;

  class Matrix; // general matrix
  class SymmetricMatrix;// packed symmetric matrix

  /****************************************************************
   ************************** Vector ******************************
   ****************************************************************/

  class Vector : private RefArr<double> {
   explicit Vector(const RefArr<double>& a) : RefArr<double>(a) {}

  public:
    typedef       double*       iterator;
    typedef const double* const_iterator;
    typedef       double      value_type;

    Vector () {}
    explicit Vector (int_t s)            : RefArr<double>(s)    {}
    Vector          (int_t s, double p)  : RefArr<double>(s, p) {}

    bool isinit () const                       { return RefArr<double>::isinit(); }
    int_t  size   () const  { return RefArr<double>::size(); }

    double*       begin ()       { return  RefArr<double>::begin(); }
    double*         end ()       { return  RefArr<double>::end(); }
    const double* begin () const { return  RefArr<double>::begin(); }
    const double*   end () const { return  RefArr<double>::end(); }

    operator       double* ()       { return RefArr<double>::operator       double*(); }
    operator const double* () const { return RefArr<double>::operator const double*(); }

    double&       operator[] (int_t i)       { return RefArr<double>::operator[](i); }
    const double& operator[] (int_t i) const { return RefArr<double>::operator[](i); }

    double&       back ()       { return RefArr<double>::back(); }
    const double& back () const { return RefArr<double>::back(); }

    double&       front ()       { return RefArr<double>::front(); }
    const double& front () const { return RefArr<double>::front(); }

    void resize (int_t s) { RefArr<double>::resize(s); }

    Vector copy () const {return Vector(RefArr<double>::copy()); }

    // scalar product
    //
    double operator* (const Vector&)      const;
    double operator* (const double*)      const;
    double operator* (ConstSlice<double>) const;
    double vdot      ()                   const;
    
    // block operations
    Vector operator* (const Matrix&) const ;
    Vector operator* (const SymmetricMatrix&) const ;

    Vector operator+ (const Vector&) const ;
    Vector operator- (const Vector&) const ;

    const Vector& operator+= (const Vector&) ;
    const Vector& operator-= (const Vector&) ;

    const Vector& operator=  (const double*);

    const Vector& operator-        ();
    const Vector& operator=  (double);
    const Vector& operator*= (double);
    const Vector& operator/= (double);
  };

  inline double operator* (const double* p, Vector v) { return v * p; }
  
  inline double operator* (ConstSlice<double> p, Vector v) { return v * p; }
  
  inline Vector Vector::operator+ (const Vector& v) const 
  {
    return Vector(RefArr<double>::operator+(v));
  }

  inline Vector Vector::operator- (const Vector& v) const 
  {
    return Vector(RefArr<double>::operator-(v));
  }

  inline const Vector& Vector::operator+= (const Vector& m) 
  {
    RefArr<double>::operator+=(m);
    return *this;
  }

  inline const Vector& Vector::operator-= (const Vector& m) 
  {
    RefArr<double>::operator-=(m);
    return *this;
  }

  inline const Vector& Vector::operator-  ()
  {
    RefArr<double>::operator-();
    return *this;
  }

  inline const Vector& Vector::operator= (double d)
  {
    RefArr<double>::operator=(d);
    return *this;
  }

  inline const Vector& Vector::operator= (const double* p)
  {
    RefArr<double>::operator=(p);
    return *this;
  }

  inline const Vector& Vector::operator*= (double d)
  {
    RefArr<double>::operator*=(d);
    return *this;
  }

  inline const Vector& Vector::operator/= (double d)
  {
    RefArr<double>::operator/=(d);
    return *this;
  }


  /****************************************************************
   ************************** Matrix ******************************
   ****************************************************************/

  // Fortran style indexing
  class Matrix : private RefArr<double> {
    SharedPointer<int_t> _size1;
    SharedPointer<int_t> _size2;

    void _check_dim   (const Matrix&) const ;
    void _check_index (int_t, int_t)  const ;

  protected:
    Matrix (const Matrix&, int_t);// copy constructor by value

  public:
    void resize (int_t, int_t) ;
    void resize (int_t s)       { resize(s, s); }

    bool isinit () const { return _size1; }

    Matrix () {}
    explicit Matrix  (int_t s1)  { resize(s1);     }
    Matrix (int_t s1, int_t s2)  { resize(s1, s2); }

    operator       double* ()       { return RefArr<double>::operator       double*(); }
    operator const double* () const { return RefArr<double>::operator const double*(); }

    Matrix copy () const { return Matrix(*this, 0); }

    const double&  operator() (int_t, int_t) const ;
    double&        operator() (int_t, int_t)       ;

    int_t size1 () const ;
    int_t size2 () const ;
    int_t size  () const ; // square matrix size

    // matrix-matrix multiplication
    Matrix operator* (const Matrix&)          const ;
    Matrix operator* (const SymmetricMatrix&) const ;

    // matrix-vector multiplication
    Vector operator* (const Vector&) const ;
    Vector operator* (const double*) const ;

    Matrix operator+ (const Matrix&) const ;
    Matrix operator- (const Matrix&) const ;

    Matrix& operator+= (const Matrix&) ;
    Matrix& operator-= (const Matrix&) ;

    Matrix& operator-        ();
    Matrix& operator=  (double);
    Matrix& operator*= (double);
    Matrix& operator/= (double);

    Slice<double> row      (int_t)    ;
    Slice<double> column   (int_t)    ;
    Slice<double> diagonal (int_t =0) ;

    ConstSlice<double> row      (int_t)    const ;
    ConstSlice<double> column   (int_t)    const ;
    ConstSlice<double> diagonal (int_t =0) const ;

    Matrix transpose () const;

    Matrix  invert ()              const;
    Vector  invert (const Vector&) const;
    Matrix  invert (const Matrix&) const;

    // orthogonalize column-wise
    //
    double orthogonalize (double = -1.);

    // null space of the matrix
    //
    Matrix kernel () const;
  };
  
  Vector operator* (const double*, const Matrix&) ;

  inline int_t Matrix::size1 () const 
  { 
    const char funame [] = "Lapack::Matrix::size1: ";

    if(_size1) 
      return *_size1; 

    std::cerr << funame << "not initialized\n";
    throw Error::Init();

  }

  inline int_t Matrix::size2 () const 
  { 
    const char funame [] = "Lapack::Matrix::size2: ";

    if(_size2) 
      return *_size2; 

    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  inline int_t Matrix::size () const  
  {
    const char funame [] = "Lapack::Matrix::size: ";

    if(size1() != size2()) {
      std::cerr << funame << "not square\n";
      throw Error::Logic();
    }

    return size1();
  }

  inline void Matrix::resize(int_t s1, int_t s2) 
  {
    const char funame [] = "Lapack::Matrix::resize: ";

    if(s1 < 0 || s2 < 0) {
      std::cerr << funame << "negative dimensions: s1 = " << s1 << " s2 = " << s2 << "\n";
      throw Error::Range();
    }

    if(!s1 || !s2)
      s1 = s2 = 0;

    if(!isinit()) {
      _size1.init(new int_t(s1));
      _size2.init(new int_t(s2));
    }
    else {
      *_size1 = s1; 
      *_size2 = s2;
    } 

    RefArr<double>::resize(s1 * s2);
  }

  // copy constructor by value
  inline Matrix::Matrix (const Matrix& m, int_t)
  {
    if(m.isinit()) {
      _size1.init(new int_t(m.size1())); 
      _size2.init(new int_t(m.size2()));
      RefArr<double>::operator=(m.RefArr<double>::copy());
    }
  }

  inline void Matrix::_check_index (int_t i1, int_t i2) const 
  {
    const char funame [] = "Lapack::Matrix::_check_index: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    if(i1 < 0 || i1 >= size1() || i2 < 0 || i2 >= size2()) {
      std::cerr << funame << "out of range: i1 = " << i1 << " i2 = " << i2 << "\n";;
      throw Error::Range();
    }
  }

  inline const double& Matrix::operator() (int_t i1, int_t i2) const 
  {
    _check_index(i1, i2);
    return RefArr<double>::operator[](i1 + size1() * i2);
  }

  inline double& Matrix::operator() (int_t i1, int_t i2) 
  {
    _check_index(i1, i2);
    return RefArr<double>::operator[](i1 + size1() * i2);
  }

  inline void Matrix::_check_dim (const Matrix& m) const 
  {
    const char funame [] = "Lapack::Matrix::_check_dim: ";
    
    if(!isinit() || !m.isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    if(m.size1() != size1() || m.size2() != size2()) {
      std::cerr << funame << "dimensions mismatch:" 
		<< " mat1_size1 = " << size1() 
		<< " mat1_size2 = " << size2() 
		<< " mat2_size1 = " << m.size1() 
		<< " mat2_size2 = " << m.size2()
		<< "\n";
      throw Error::Range();
    }
  }

  inline Matrix Matrix::operator+ (const Matrix& m) const 
  {
    Matrix res(*this, 0);
    res += m;
    return res;
  }

  inline Matrix Matrix::operator- (const Matrix& m) const 
  {
    Matrix res(*this, 0);
    res -= m;
    return res;
  }

  inline Matrix& Matrix::operator+= (const Matrix& m) 
  {
    _check_dim(m);
    RefArr<double>::operator+=(m);
    return *this;
  }

  inline Matrix& Matrix::operator-= (const Matrix& m) 
  {
    _check_dim(m);
    RefArr<double>::operator-=(m);
    return *this;
  }

  // matrix-vector product
  inline Vector Matrix::operator* (const Vector& v) const 
  {
    const char funame [] = "Lapack::Matrix::operator*: ";

    if(!isinit() || !v.isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    if(size2() != v.size()) {
      std::cerr << funame << "the dimensions are different\n";
      throw Error::Range();
    }

    return *this * (const double*)v;
  }

  // arithmetic operations with double
  inline Matrix& Matrix::operator-  ()
  {
    const char funame [] = "Lapack::Matrix::operator-: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator-();
    return *this;
  }

  inline Matrix& Matrix::operator= (double d)
  {
    const char funame [] = "Lapack::Matrix::operator=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator=(d);
    return *this;
  }

  inline Matrix& Matrix::operator*= (double d)
  {
    const char funame [] = "Lapack::Matrix::operator*=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator*=(d);
    return *this;
  }

  inline Matrix& Matrix::operator/= (double d)
  {
    const char funame [] = "Lapack::Matrix::operator/=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator/=(d);
    return *this;
  }

  /****************************************************************
   ********************* Symmetric Matrix *************************
   ****************************************************************/

  // packed symmetric matrix with upper triangle reference
  class SymmetricMatrix : private RefArr<double> {
    SharedPointer<int_t> _size;
    explicit SymmetricMatrix (const SymmetricMatrix&, int_t); // copy constructor by value

  public:
    void resize (int_t) ;

    bool isinit () const { return _size; }

    SymmetricMatrix () {}
    explicit SymmetricMatrix (int_t s)  { resize(s); }
    explicit SymmetricMatrix (const Matrix&, char ='U')       ;

    operator       double* ()       { return RefArr<double>::operator       double*(); }
    operator const double* () const { return RefArr<double>::operator const double*(); }

    SymmetricMatrix copy () const { return SymmetricMatrix(*this, 0); }

    int_t size () const ;

    // referencing as an upper triangle column-wise
    const double& operator() (int_t, int_t) const ;
    double&       operator() (int_t, int_t)       ;

    Vector operator* (const Vector&)      const;
    Vector operator* (const double*)      const;
    Vector operator* (ConstSlice<double>) const;

    Matrix operator* (const Matrix&)          const ;
    Matrix operator* (const SymmetricMatrix&) const ;
    
    SymmetricMatrix operator+ (const SymmetricMatrix&) const ;
    SymmetricMatrix operator- (const SymmetricMatrix&) const ;

    SymmetricMatrix operator+= (const SymmetricMatrix&) ;
    SymmetricMatrix operator-= (const SymmetricMatrix&) ;

    SymmetricMatrix operator-  ()      ;
    SymmetricMatrix operator=  (double);
    SymmetricMatrix operator+= (double);
    SymmetricMatrix operator-= (double);
    SymmetricMatrix operator*= (double);
    SymmetricMatrix operator/= (double);

    Vector    eigenvalues (Matrix* =0) const ;

    SymmetricMatrix invert ()             const ;
    SymmetricMatrix positive_invert ()    const ;
  };

  inline Vector operator* (ConstSlice<double> v, SymmetricMatrix m) { return m * v; }
  
  inline Vector operator* (const double* v, SymmetricMatrix m) { return m * v; }
  
  inline void SymmetricMatrix::resize(int_t s) 
  {
    const char funame [] = "Lapack::SymmetricMatrix::resize: ";

    if(s < 0) {
      std::cerr << funame << "negative size: s = " 
		<< s << "\n";
      throw Error::Range();
    }

    if(!_size)
      _size.init(new int_t(s));
    else
      *_size = s;
    
    RefArr<double>::resize(s*(s+1)/2);
  }

  inline int_t SymmetricMatrix::size () const 
  { 
    const char funame [] = "Lapack::SymmetricMatrix::size: ";

    if(_size) 
      return *_size; 

    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  // copy constructor by value
  inline SymmetricMatrix::SymmetricMatrix (const SymmetricMatrix& m, int_t)
  {
    if(m.isinit()) {
      _size.init(new int_t(m.size()));
      RefArr<double>::operator=(m.RefArr<double>::copy());
    }
  }

  inline double& SymmetricMatrix::operator() (int_t i1, int_t i2) 
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator(): ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    if(i1 > i2) 
      std::swap(i1, i2);

    if(i1 < 0 || i2 >= size()) {
      std::cerr << funame << "out of range: i1 = " << i1 
		<< " i2 = " << i2 << std::endl;
      throw Error::Range();
    }

    return RefArr<double>::operator[](i1 + i2 * (i2 + 1) / 2);
  }

  inline const double& SymmetricMatrix::operator() (int_t i1, int_t i2) const 
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator(): ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    if(i1 > i2) 
      std::swap(i1, i2);

    if(i1 < 0 || i2 >= size()) {
      std::cerr << funame << "out of range: i1 = " << i1 
		<< " i2 = " << i2 << std::endl;
      throw Error::Range();
    }

    return RefArr<double>::operator[](i1 + i2 * (i2 + 1) / 2);
  }

  inline SymmetricMatrix SymmetricMatrix::operator+ (const SymmetricMatrix& m)
    const 
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator+: ";

    if(!isinit() || !m.isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    SymmetricMatrix res(*this, 0);
    res += m;
    return res;
  }

  inline SymmetricMatrix SymmetricMatrix::operator- (const SymmetricMatrix& m)
    const 
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator-: ";

    if(!isinit() || !m.isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    SymmetricMatrix res(*this, 0);
    res -= m;
    return res;
  }

  inline SymmetricMatrix SymmetricMatrix::operator+= (const SymmetricMatrix& m) 
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator+=: ";

    if(!isinit() || !m.isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator+=(m);
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator-= (const SymmetricMatrix& m) 
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator-=: ";

    if(!isinit() || !m.isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator-=(m);
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator-  ()
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator-: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator-();
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator= (double d)
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator=(0.);
    double* p = *this;
    int_t n = size() + 2;
    for(int_t i = 2; i < n; ++i) {
      *p = d;
      p += i;
    }
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator+= (double d)
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator+=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    double* p = *this;
    int_t n = size() + 2;
    for(int_t i = 2; i < n; ++i) {
      *p += d;
      p += i;
    }
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator-= (double d)
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator-=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    double* p = *this;
    int_t n = size() + 2;
    for(int_t i = 2; i < n; ++i) {
      *p -= d;
      p += i;
    }
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator*= (double d)
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator*=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator*=(d);
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator/= (double d)
  {
    const char funame [] = "Lapack::SymmetricMatrix::operator/=: ";

    if(!isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    RefArr<double>::operator/=(d);
    return *this;
  }

  // Generalized eigenvalue problem
  //
  Vector diagonalize(SymmetricMatrix, SymmetricMatrix, Matrix* = 0) ;

  /****************************************************************
   ******************* Band Symmetric Matrix **********************
   ****************************************************************/

  class BandMatrix : private Matrix {
    void _check_size () const ;

    BandMatrix (const BandMatrix& m, int_t) : Matrix(m, 0) { } // copy construction by value

  public:
    BandMatrix () {}
    BandMatrix   (int_t s, int_t b)  : Matrix(b, s) { _check_size(); }
    void  resize (int_t s, int_t b) { Matrix::resize(b, s); _check_size(); }

    bool isinit () const { return Matrix::isinit(); }

    BandMatrix copy () const { return BandMatrix(*this, 0); }

    int_t size      () const  { return size2(); }
    int_t band_size () const  { return size1(); }

    double  operator() (int_t, int_t) const;
    double& operator() (int_t, int_t) ;

    BandMatrix& operator= (double d) { Matrix::operator=(d); return *this; }

    Vector eigenvalues (Matrix* =0) const ;
  };

  inline void BandMatrix::_check_size () const 
  {
    const char funame [] = "Lapack::BandMatrix::_check_size(): ";
    
    if(band_size() <= size())
      return;

    std::cerr << funame << "band size bigger than the matrix size\n";
    throw Error::Range();
  }

  /****************************************************************
   ********************** LU Factorization ************************
   ****************************************************************/
 
  int_t parity (RefArr<int_t>);

  class LU : private Matrix {
    RefArr<int_t> _ipiv;

  public:
    explicit LU (const Matrix&) ;
    
    int_t size ()            const { return Matrix::size1(); }
    RefArr<int_t> ipiv ()  const { return _ipiv.copy(); }
    double det ()          const;

    Matrix invert ()              const ; // inverse matrix
    Vector invert (const Vector&) const ; // solve linear equations
    Matrix invert (const Matrix&) const ; // solve linear equations
  };

  /****************************************************************
   ******* LU Factorization for symmetric packed matrices *********
   ****************************************************************/
  class SymLU : private SymmetricMatrix {
    RefArr<int_t> _ipiv;

  public:
    explicit SymLU (const SymmetricMatrix&) ;
    int_t size () const { return SymmetricMatrix::size(); }
    double det () const;

    SymmetricMatrix invert ()              const ; // inverse matrix
    Vector          invert (const Vector&) const ; // solve linear equations
    Matrix          invert (const Matrix&) const ; // solve linear equations
  };

  /****************************************************************
   ******************* Cholesky Factorization *********************
   ****************************************************************/

  // for positively defined matrices
  class Cholesky : private SymmetricMatrix {

  public:
    explicit Cholesky (const SymmetricMatrix&) ;
    int_t size () const { return SymmetricMatrix::size(); }

    SymmetricMatrix invert () const ; // inverse matrix
    Vector invert (const Vector&) const ; // solve linear equations
    Matrix invert (const Matrix&) const ; // solve linear equations
    double det_sqrt ();
  };

  /****************************************************************
   ************************ Complex Matrix ************************
   ****************************************************************/

  class ComplexMatrix : public Array<complex> {
    int_t _size1;
    int_t _size2;

    void _index_check (int_t, int_t) const ;

  public:
    void resize (int_t, int_t) ;
    void resize (int_t s)       { resize(s, s); }

    ComplexMatrix () : _size1(0), _size2(0) {}
    explicit ComplexMatrix  (int_t s1)  : _size1(0), _size2(0) { resize(s1);     }
    ComplexMatrix (int_t s1, int_t s2)  : _size1(0), _size2(0) { resize(s1, s2); }

    const complex&  operator() (int_t, int_t) const ;
    complex&        operator() (int_t, int_t)       ;

    int_t size1 () const { return _size1; }
    int_t size2 () const { return _size2; }
    int_t size  () const ;

    // matrix-matrix multiplication
    ComplexMatrix operator* (const ComplexMatrix&) const ;
  };

  inline void ComplexMatrix::resize(int_t s1, int_t s2) 
  {
    const char funame [] = "Lapack::ComplexMatrix::resize: ";

    if(s1 <= 0 || s2 <=0) {
      std::cerr << funame << "wrong dimensions: s1 = " << s1 << " s2 = " << s2 << "\n";
      throw Error::Range();
    }

    _size1 = s1;
    _size2 = s2;

    Array<complex>::resize(s1*s2);
  }

  inline int_t ComplexMatrix::size () const  
  {
    const char funame [] = "Lapack::ComplexMatrix::size: ";

    if(size1() != size2()) {
      std::cerr << funame << "not square\n";
      throw Error::Logic();
    }

    return size1();
  }

  inline void ComplexMatrix::_index_check (int_t i1, int_t i2) const 
  {
    const char funame [] = "Lapack::ComplexMatrix::_index_check: ";

    if(i1 < 0 || i1 >= size1() || i2 < 0 || i2 >= size2()) {
      std::cerr << funame << "out of range: i1 = " << i1 << " i2 = " << i2 << "\n";
      throw Error::Range();
    }
  }

  inline const complex& ComplexMatrix::operator() (int_t i1, int_t i2) const 
  {
    _index_check(i1, i2);
    return Array<complex>::operator[](i1 + size1() * i2);
  }

  inline complex& ComplexMatrix::operator() (int_t i1, int_t i2) 
  {
    _index_check(i1, i2);
    return Array<complex>::operator[](i1 + size1() * i2);
  }

  /******************************************************************
   ************************ Hermitian Matrix ************************
   ******************************************************************/

  class HermitianMatrix : public Array<complex> {
    int_t _size;
    
  public:
    void resize (int_t) ;

    HermitianMatrix () : _size(0) {}
    explicit HermitianMatrix (int_t s1)  : _size(0) { resize(s1); } 
    explicit HermitianMatrix (const ComplexMatrix&, char = 'U') ; 

    const complex&  operator() (int_t, int_t) const ;
    complex&        operator() (int_t, int_t)       ;

    void operator= (complex v) { Array<complex>::operator=(v); }

    int_t size () const { return _size; }
    Vector eigenvalues (ComplexMatrix* =0) const ;
  };

  inline void HermitianMatrix::resize(int_t s) 
  {
    const char funame [] = "Lapack::HermitianMatrix::resize: ";

    if(s <= 0) {
      std::cerr << funame << "wrong size: s = " 
		<< s << std::endl;
      throw Error::Range();
    }

    _size = s;

    Array<complex>::resize(s * (s + 1) / 2);
  }

  inline complex& HermitianMatrix::operator() (int_t i1, int_t i2) 
  {
    const char funame [] = "Lapack::HermitianMatrix::operator(): ";

    if(i1 > i2) 
      std::swap(i1, i2);

    if(i1 < 0 || i2 >= size()) {
      std::cerr << funame << "out of range: i1 = " << i1 
		<< " i2 = " << i2 << std::endl;
      throw Error::Range();
    }

    return Array<complex>::operator[](i1 + i2 * (i2 + 1) / 2);
  }

  inline const complex& HermitianMatrix::operator() (int_t i1, int_t i2) const 
  {
    const char funame [] = "Lapack::HermitianMatrix::operator(): ";

    if(i1 > i2) 
      std::swap(i1, i2);

    if(i1 < 0 || i2 >= size()) {
      std::cerr << funame << "out of range: i1 = " << i1 
		<< " i2 = " << i2 << std::endl;
      throw Error::Range();
    }

    return Array<complex>::operator[](i1 + i2 * (i2 + 1) / 2);
  }

  // Generalized eigenvalue problem
  //
  Vector diagonalize(const HermitianMatrix&, const HermitianMatrix&, ComplexMatrix* = 0) ;

  /************************************************************************
   ************************** COMPLEX VECTOR ******************************
   ************************************************************************/

  class ComplexVector : public RefArr<complex> {
  public:
    ComplexVector () {}
    explicit ComplexVector (int_t s) : RefArr<complex>(s) {}
  };

  ComplexVector fourier_transform (const Vector&,                         const MultiIndexConvert&);
  ComplexVector fourier_transform (const ComplexVector&,                  const MultiIndexConvert&);
  ComplexVector fourier_transform (const std::map<int, Lapack::complex>&, const MultiIndexConvert&);

  // find a complimentary orthogonal basis set to the non-orthogonal vector set
  //
  double orthogonalize (Matrix basis, int vsize);

  // solve linear equations by svd-decomposition
  //
  Vector svd_solve (Matrix a, Vector b, double* residue = 0,double pres = -1., double (*weight)(double) = 0);

  // singular value decomposition, A = U * S * V**T
  //
  void svd (Matrix a, Vector s, Matrix u, Matrix v);
  //
}// namespace Lapack

#endif
