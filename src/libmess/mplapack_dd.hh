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

#ifdef WITH_MPLAPACK

#ifndef MPACK_DD_HH
#define MPACK_DD_HH

#include "array.hh"
#include "shared.hh"
#include "multindex.hh"


#include <qd/dd_real.h>
#include <mplapack/mplapack_config.h>

namespace Mpack_dd {

  typedef mplapackint int_t;// index type in lapack & blas libraries

  class Matrix; // general matrix
  class SymmetricMatrix;// packed symmetric matrix

  /****************************************************************
   ************************** Vector ******************************
   ****************************************************************/

  class Vector : private RefArr<dd_real> {
    //
    explicit Vector(const RefArr<dd_real>& a) : RefArr<dd_real>(a) {}
    
    void _assert () const;
    
  public:
    //
    typedef       dd_real*       iterator;
    typedef const dd_real* const_iterator;
    typedef       dd_real      value_type;

    Vector () {}
    explicit Vector (int_t s)            : RefArr<dd_real>(s)    {}
    Vector          (int_t s, dd_real p)  : RefArr<dd_real>(s, p) {}

    bool isinit () const { return RefArr<dd_real>::isinit(); }
    int_t  size () const { return RefArr<dd_real>::size(); }

    dd_real*       begin ()       { return  RefArr<dd_real>::begin(); }
    dd_real*         end ()       { return  RefArr<dd_real>::end(); }
    const dd_real* begin () const { return  RefArr<dd_real>::begin(); }
    const dd_real*   end () const { return  RefArr<dd_real>::end(); }

    operator       dd_real* ()       { return RefArr<dd_real>::operator       dd_real*(); }
    operator const dd_real* () const { return RefArr<dd_real>::operator const dd_real*(); }

    dd_real&       operator[] (int_t i)       { return RefArr<dd_real>::operator[](i); }
    const dd_real& operator[] (int_t i) const { return RefArr<dd_real>::operator[](i); }

    dd_real&       back ()       { return RefArr<dd_real>::back(); }
    const dd_real& back () const { return RefArr<dd_real>::back(); }

    dd_real&       front ()       { return RefArr<dd_real>::front(); }
    const dd_real& front () const { return RefArr<dd_real>::front(); }

    void resize (int_t s) { RefArr<dd_real>::resize(s); }

    Vector copy () const {return Vector(RefArr<dd_real>::copy()); }

    // scalar product
    //
    dd_real operator* (Vector)              const;
    dd_real operator* (const dd_real*)      const;
    dd_real operator* (ConstSlice<dd_real>) const;
    dd_real vdot      ()                    const;
    
    // block operations
    //
    Vector operator* (const Matrix&)          const;
    Vector operator* (const SymmetricMatrix&) const;

    Vector operator+ (Vector v) const { return Vector(RefArr<dd_real>::operator+(v)); }
    Vector operator- (Vector v) const { return Vector(RefArr<dd_real>::operator-(v)); }

    Vector operator+= (Vector m) { RefArr<dd_real>::operator+=(m); return *this; }
    Vector operator-= (Vector m) { RefArr<dd_real>::operator-=(m); return *this; }

    Vector operator-  ()          { RefArr<dd_real>::operator-  (); return *this; }
    Vector operator=  (dd_real d) { RefArr<dd_real>::operator= (d); return *this; }
    Vector operator*= (dd_real d) { RefArr<dd_real>::operator*=(d); return *this; }
    Vector operator/= (dd_real d) { RefArr<dd_real>::operator/=(d); return *this; }

    Vector operator= (const dd_real* p) { RefArr<dd_real>::operator=(p); return *this; }
  };

  inline dd_real operator* (const dd_real*      p, Vector v) { return v * p; }
  
  inline dd_real operator* (ConstSlice<dd_real> p, Vector v) { return v * p; }

  inline void Vector::_assert () const
  {
    const char funame [] = "Mpack_dd::Vector::_assert: ";

    if(isinit())
      //
      return;

    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }
      
  /****************************************************************
   ************************** Matrix ******************************
   ****************************************************************/

  // Fortran style indexing
  //
  class Matrix : private RefArr<dd_real> {
    //
    SharedPointer<int_t> _size1;
    SharedPointer<int_t> _size2;

    void _check_dim   (const Matrix&) const;

    void _assert ()         const;
    void _assert (int, int) const;

    dd_real&       _at (int i)       { return RefArr<dd_real>::operator[](i); }
    const dd_real& _at (int i) const { return RefArr<dd_real>::operator[](i); }

  public:
    //
    void resize (int_t, int_t);
    
    void resize (int_t s)  { resize(s, s); }

    bool isinit () const { return _size1; }

    Matrix () {}
    explicit Matrix  (int_t s1)  { resize(s1);     }
    Matrix (int_t s1, int_t s2)  { resize(s1, s2); }

    explicit Matrix (const SymmetricMatrix&);

    operator       dd_real* ()       { return RefArr<dd_real>::operator       dd_real*(); }

    operator const dd_real* () const { return RefArr<dd_real>::operator const  dd_real*(); }

    Matrix copy () const;

    dd_real&       operator() (int_t i1, int_t i2)       { _assert(i1, i2); return _at(i1 + size1() * i2); }

    const dd_real& operator() (int_t i1, int_t i2) const { _assert(i1, i2); return _at(i1 + size1() * i2); }

    int_t size1 () const { _assert(); return *_size1; }
    int_t size2 () const { _assert(); return *_size2; }
    
    int_t size  () const; // square matrix size

    // matrix-matrix multiplication
    //
    Matrix operator* (const          Matrix&) const;
    Matrix operator* (const SymmetricMatrix&) const;

    // matrix-vector multiplication
    //
    Vector operator* (Vector)         const;
    Vector operator* (const dd_real*) const;

    Matrix operator+= (Matrix m) { _check_dim(m); RefArr<dd_real>::operator+=(m); return *this; }
    Matrix operator-= (Matrix m) { _check_dim(m); RefArr<dd_real>::operator-=(m); return *this; }

    Matrix operator+ (Matrix m) const { Matrix res = copy(); res += m; return res; }
    Matrix operator- (Matrix m) const { Matrix res = copy(); res -= m; return res; }

    // arithmetic operations with dd_real
    //
    Matrix operator-  ()          { _assert(); RefArr<dd_real>::operator-  (); return *this; }
    Matrix operator=  (dd_real d) { _assert(); RefArr<dd_real>::operator= (d); return *this; }
    Matrix operator*= (dd_real d) { _assert(); RefArr<dd_real>::operator*=(d); return *this; }
    Matrix operator/= (dd_real d) { _assert(); RefArr<dd_real>::operator/=(d); return *this; }

    Slice<dd_real> row      (int_t);
    Slice<dd_real> column   (int_t);
    Slice<dd_real> diagonal (int_t =0);

    ConstSlice<dd_real> row      (int_t)    const;
    ConstSlice<dd_real> column   (int_t)    const;
    ConstSlice<dd_real> diagonal (int_t =0) const;

    Matrix transpose () const;

    Matrix  invert ()       const;
    Vector  invert (Vector) const;
    Matrix  invert (Matrix) const;

    // orthogonalize column-wise
    //
    dd_real orthogonalize (dd_real = -1.);

    // null space of the matrix
    //
    Matrix kernel () const;

    Vector eigenvalues (Matrix* =0) const;
    //
  };// Matrix
  
  Vector operator* (const dd_real*, const Matrix&);

  inline void Matrix::_assert () const
  {
    const char funame [] = "Mpack_dd::Matrix::_assert: ";

    if(isinit())
      //
      return;

    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }
  
  inline void Matrix::_assert (int i1, int i2) const 
  {
    const char funame [] = "Mpack_dd::Matrix::_assert: ";

    _assert();
    
    if(i1 < 0 || i1 >= size1() || i2 < 0 || i2 >= size2()) {
      //
      std::cerr << funame << "out of range: " << i1 << ", " << i2 << "\n";
      
      throw Error::Range();
    }
  }

  inline void Matrix::_check_dim (const Matrix& m) const 
  {
    const char funame [] = "Mpack_dd::Matrix::_check_dim: ";

    _assert();

    m._assert();
    
    if(m.size1() != size1() || m.size2() != size2()) {
      //
      std::cerr << funame << "dimensions mismatch: " << m.size1() << ", " << m.size2() << "\n";
      
      throw Error::Range();
    }
  }

  inline int_t Matrix::size () const
  {
    const char funame [] = "Mpack_dd::Matrix::size: ";

    _assert();
    
    if(size1() != size2()) {
      //
      std::cerr << funame << "not square: " << size1() << ", " << size2() << "\n";
      
      throw Error::Logic();
    }

    return size1();
  }

  inline void Matrix::resize(int_t s1, int_t s2) 
  {
    const char funame [] = "Mpack_dd::Matrix::resize: ";

    if(s1 <= 0 || s2 <= 0) {
      //
      std::cerr << funame << "out of range: " << s1 << ", " << s2 << "\n";
      
      throw Error::Range();
    }

    if(!isinit()) {
      //
      _size1.init(new int_t(s1));
      
      _size2.init(new int_t(s2));
    }
    else {
      //
      *_size1 = s1;
      
      *_size2 = s2;
    } 

    RefArr<dd_real>::resize(s1 * s2);
  }

  // matrix-vector product
  //
  inline Vector Matrix::operator* (Vector v) const
  {
    const char funame [] = "Mpack_dd::Matrix::operator*: ";

    _assert();

    if(size2() != v.size()) {
      //
      std::cerr << funame << "dimensions mismatch: " << v.size() << "\n";
      
      throw Error::Range();
    }

    return *this * (const dd_real*)v;
  }

  /****************************************************************
   ********************* Symmetric Matrix *************************
   ****************************************************************/

  // packed symmetric matrix with upper triangle reference
  //
  class SymmetricMatrix : private RefArr<dd_real> {
    //
    SharedPointer<int_t> _size;

    void _assert ()    const;
    void _assert (int) const;

    dd_real&       _at(int i)       { return RefArr<dd_real>::operator[](i); }
    const dd_real& _at(int i) const { return RefArr<dd_real>::operator[](i); }
    
  public:
    //
    void resize (int_t);

    bool isinit () const { return _size; }

    SymmetricMatrix () {}
    
    explicit SymmetricMatrix (int_t s)  { resize(s); }
    
    explicit SymmetricMatrix (Matrix, char ='U');

    operator       dd_real* ()       { _assert(); return RefArr<dd_real>::operator       dd_real*(); }
    operator const dd_real* () const { _assert(); return RefArr<dd_real>::operator const dd_real*(); }

    // copy by value
    //
    SymmetricMatrix copy () const;

    int_t size () const { _assert(); return *_size; }

    // referencing as an upper triangle column-wise
    //
    dd_real&       operator() (int_t i1, int_t i2)       { _assert(i1); _assert(i2); if(i1 > i2) std::swap(i1, i2); return _at(i1 + i2 * (i2 + 1) / 2); }
    const dd_real& operator() (int_t i1, int_t i2) const { _assert(i1); _assert(i2); if(i1 > i2) std::swap(i1, i2); return _at(i1 + i2 * (i2 + 1) / 2); }

    Vector operator* (Vector)              const;
    Vector operator* (const dd_real*)      const;
    Vector operator* (ConstSlice<dd_real>) const;

    Matrix operator* (const          Matrix&) const;
    Matrix operator* (const SymmetricMatrix&) const;
    
    SymmetricMatrix operator+= (SymmetricMatrix m)       { _assert(); m._assert(); RefArr<dd_real>::operator+=(m); return *this; }
    SymmetricMatrix operator-= (SymmetricMatrix m)       { _assert(); m._assert(); RefArr<dd_real>::operator-=(m); return *this; }

    SymmetricMatrix operator+  (SymmetricMatrix m) const { _assert(); m._assert(); SymmetricMatrix res = copy(); res += m; return res; }
    SymmetricMatrix operator-  (SymmetricMatrix m) const { _assert(); m._assert(); SymmetricMatrix res = copy(); res -= m; return res; }

    SymmetricMatrix operator = (dd_real);
    SymmetricMatrix operator+= (dd_real);
    SymmetricMatrix operator-= (dd_real);
    
    SymmetricMatrix operator - (dd_real d) { _assert(); RefArr<dd_real>::operator-  (); return *this; }
    SymmetricMatrix operator*= (dd_real d) { _assert(); RefArr<dd_real>::operator*=(d); return *this; }
    SymmetricMatrix operator/= (dd_real d) { _assert(); RefArr<dd_real>::operator/=(d); return *this; }

    Vector eigenvalues (Matrix* =0) const;

    SymmetricMatrix          invert () const;
    SymmetricMatrix positive_invert () const;
  };

  inline Vector operator* (ConstSlice<dd_real> v, const SymmetricMatrix& m) { return m * v; }
  inline Vector operator* (const      dd_real* v, const SymmetricMatrix& m) { return m * v; }
  
  inline void SymmetricMatrix::_assert () const
  {
    const char funame [] = "Mpack_dd::SymmetricMatrix::_assert: ";

    if(!isinit()) {
      //
      std::cerr << funame << "not initialized\n";
      
      throw Error::Init();
    }
  }
  
  inline void SymmetricMatrix::_assert (int i) const
  {
    const char funame [] = "Mpack_dd::SymmetricMatrix::_assert: ";

    _assert();

    if(i < 0 || i >= size()) {
      //
      std::cerr << funame << "out of range: " << i << "\n";

      throw Error::Range();
    }
  }
  
  inline void SymmetricMatrix::resize(int_t s) 
  {
    const char funame [] = "Mpack_dd::SymmetricMatrix::resize: ";

    if(s <= 0) {
      //
      std::cerr << funame << "out of range: " << s << "\n";
      
      throw Error::Range();
    }

    if(!_size) {
      //
      _size.init(new int_t(s));
    }
    else
      //
      *_size = s;
    
    RefArr<dd_real>::resize(s * (s + 1) / 2);
  }

  // Generalized eigenvalue problem
  //
  Vector diagonalize(SymmetricMatrix, SymmetricMatrix, Matrix* = 0) ;

  /****************************************************************
   ******************* Band Symmetric Matrix **********************
   ****************************************************************/

  class BandMatrix : private Matrix {
    //
    void _assert () const;

    std::pair<int_t, int_t> _convert (int_t, int_t) const;
    
  public:
    //
    BandMatrix () {}

    void resize (int_t s, int_t b) { Matrix::resize(b, s); _assert(); }

    BandMatrix (int_t s, int_t b) { resize(s, b); }
    
    bool isinit () const { return Matrix::isinit(); }

    int_t      size () const  { return size2(); }
    int_t band_size () const  { return size1(); }

    dd_real&       operator() (int_t i, int_t j)       { std::pair<int_t, int_t> p = _convert(i, j); return Matrix::operator()(p.first, p.second); }
    const dd_real& operator() (int_t i, int_t j) const { std::pair<int_t, int_t> p = _convert(i, j); return Matrix::operator()(p.first, p.second); }
    
    BandMatrix& operator= (dd_real d) { Matrix::operator=(d); return *this; }

    Vector eigenvalues (Matrix* =0) const;
  };

  inline void BandMatrix::_assert ()  const
  {
    const char funame [] = "Mpack_dd::BandMatrix::_assert: ";
    
    if(!isinit()) {
      //
      std::cerr << funame << "not initialized\n";
    
      throw Error::Init();
    }

    if(band_size() <= size())
      //
      return;

    std::cerr << funame << "out of range: " << band_size() << ", " << size() << "\n";
    
    throw Error::Range();
  }

  /****************************************************************
   ********************** LU Factorization ************************
   ****************************************************************/
 
  int_t parity (RefArr<int_t>);

  class LU : private Matrix {
    //
    mutable ::RefArr<int_t> _ipiv;

  public:
    //
    explicit LU (Matrix);
    
    int_t    size () const { return Matrix::size1(); }
    dd_real   det () const;

    Matrix invert ()       const; // inverse matrix
    Vector invert (Vector) const; // solve linear equations
    Matrix invert (Matrix) const; // solve linear equations
  };

  /****************************************************************
   ******* LU Factorization for symmetric packed matrices *********
   ****************************************************************/
  
  class SymLU : private SymmetricMatrix {
    //
    mutable ::RefArr<int_t> _ipiv;

  public:
    //
    explicit SymLU (SymmetricMatrix);
    
    int_t   size () const { return SymmetricMatrix::size(); }
    dd_real  det () const;

    SymmetricMatrix invert ()       const; // inverse matrix
    Vector          invert (Vector) const; // solve linear equations
    Matrix          invert (Matrix) const; // solve linear equations
  };

  /****************************************************************
   ******************* Cholesky Factorization *********************
   ****************************************************************/

  // for positively defined matrices
  //
  class Cholesky : private SymmetricMatrix {

  public:
    //
    explicit Cholesky (SymmetricMatrix);
    
    int_t size () const { return SymmetricMatrix::size(); }

    SymmetricMatrix invert ()       const;  // inverse matrix
    Vector          invert (Vector) const;  // solve linear equations
    Matrix          invert (Matrix) const;  // solve linear equations
    
    dd_real det_sqrt () const;
  };

}// namespace Mpack_dd

#endif
#endif
