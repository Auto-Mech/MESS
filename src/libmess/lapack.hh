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

#ifdef USE_INT64
  typedef int64_t int_t;
#else
  typedef int32_t int_t;
#endif

  typedef std::complex<double> complex;

  class Matrix; // general matrix
  class SymmetricMatrix;// packed symmetric matrix

  /****************************************************************
   ************************** Vector ******************************
   ****************************************************************/

  class Vector : private RefArr<double> {
    //
   explicit Vector(const RefArr<double>& a) : RefArr<double>(a) {}

  public:
    //
    void _assert ()    const;

    typedef       double*       iterator;
    typedef const double* const_iterator;
    typedef       double      value_type;

    Vector () {}
    explicit Vector (int_t s)            : RefArr<double>(s)    {}
    Vector          (int_t s, double p)  : RefArr<double>(s, p) {}

    template<typename V> explicit Vector (const V& v) : RefArr<double>(v) {}
    
    bool isinit   () const { return RefArr<double>::isinit(); }
    int_t  size   () const { if(isinit()) return RefArr<double>::size(); return 0; }

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
    double operator* (Vector)             const;
    double operator* (const double*)      const;
    double operator* (ConstSlice<double>) const;
    double vdot      ()                   const;
    
    // block operations
    //
    Vector operator* (Matrix) const;
    Vector operator* (SymmetricMatrix) const;

    Vector operator+ (Vector v) const { return Vector(RefArr<double>::operator+(v)); }
    Vector operator- (Vector v) const { return Vector(RefArr<double>::operator-(v)); }

    Vector operator=  (const double* p) { RefArr<double>::operator= (p); return *this; }
    Vector operator+= (Vector m)        { RefArr<double>::operator+=(m); return *this; }
    Vector operator-= (Vector m)        { RefArr<double>::operator-=(m); return *this; }

    Vector operator-  ()         { RefArr<double>::operator-  (); return *this; }
    Vector operator=  (double d) { RefArr<double>::operator= (d); return *this; }
    Vector operator*= (double d) { RefArr<double>::operator*=(d); return *this; }
    Vector operator/= (double d) { RefArr<double>::operator/=(d); return *this; }
  };

  inline void Vector::_assert () const
  {
    const char funame [] = "Lapack::Vector::_assert: ";

    if(isinit())
      //
      return;

    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }
  
  inline double operator* (const double* p, Vector v) { return v * p; }
  
  inline double operator* (ConstSlice<double> p, Vector v) { return v * p; }
  
  /****************************************************************
   ************************** Matrix ******************************
   ****************************************************************/

  // Fortran style indexing
  //
  class Matrix : private RefArr<double> {
    //
    SharedPointer<int_t> _size1;
    
    SharedPointer<int_t> _size2;

  protected:
    //
    Matrix (Matrix, int_t);// copy constructor by value

  public:
    
    void _assert () const;
    
    void _assert (int_t, int_t)  const;

    void _assert   (Matrix) const;

    void _assert   (Vector) const;

    void resize (int_t, int_t);
    
    void resize (int_t s) { resize(s, s); }

    bool isinit () const { return _size1; }

    Matrix () {}
    
    explicit Matrix  (int_t s1)  { resize(s1); }
    
    Matrix (int_t s1, int_t s2)  { resize(s1, s2); }

    operator       double* ()       { return RefArr<double>::operator       double*(); }
    
    operator const double* () const { return RefArr<double>::operator const double*(); }

    Matrix copy () const { return Matrix(*this, 0); }

    Matrix operator= (SymmetricMatrix);
    
    Matrix (SymmetricMatrix m);

    template<typename V> Matrix (const ::Matrix<V>&);
    
    template<typename V> Matrix operator= (const ::Matrix<V>&);
    
    const double& operator() (int_t i1, int_t i2) const { _assert(i1, i2); return RefArr<double>::operator[](i1 + size1() * i2); }

    double&       operator() (int_t i1, int_t i2)       { _assert(i1, i2); return RefArr<double>::operator[](i1 + size1() * i2); }

    inline int_t size1 () const { if(_size1) return *_size1; return 0; }
    
    inline int_t size2 () const { if(_size2) return *_size2; return 0; }
    
    int_t size  () const; // square matrix size

    // matrix-matrix multiplication
    //
    Matrix operator* (Matrix)          const;
    Matrix operator* (SymmetricMatrix) const;

    // matrix-vector multiplication
    //
    Vector operator* (const double*)   const;
    Vector operator* (Vector v)        const { _assert(v); return *this * (const double*)v; }

    // arithmetic operations with matrix
    //
    Matrix operator+= (Matrix m) { _assert(m); RefArr<double>::operator+=(m); return *this; }
    Matrix operator-= (Matrix m) { _assert(m); RefArr<double>::operator-=(m); return *this; }

    Matrix operator+  (Matrix m) const { Matrix res(*this, 0); res += m; return res;  }
    Matrix operator-  (Matrix m) const { Matrix res(*this, 0); res -= m; return res;  }

    // arithmetic operations with double
    //
    Matrix operator-  ()         { _assert(); RefArr<double>::operator-  ();  return *this; }
    Matrix operator=  (double d) { _assert(); RefArr<double>::operator=  (d); return *this; }
    Matrix operator*= (double d) { _assert(); RefArr<double>::operator*= (d); return *this; }
    Matrix operator/= (double d) { _assert(); RefArr<double>::operator/= (d); return *this; }

    Slice<double> row      (int_t)    ;
    Slice<double> column   (int_t)    ;
    Slice<double> diagonal (int_t =0) ;

    ConstSlice<double> row      (int_t)    const ;
    ConstSlice<double> column   (int_t)    const ;
    ConstSlice<double> diagonal (int_t =0) const ;

    Matrix transpose () const;

    Matrix  invert ()       const;
    Vector  invert (Vector) const;
    Matrix  invert (Matrix) const;

    // orthogonalize column-wise
    //
    double orthogonalize (double = -1.);

    // null space of the matrix
    //
    Matrix kernel () const;

    // eigenvalues
    //
    Vector    eigenvalues (Matrix* =0) const;
  };
  
  Vector operator* (const double*, Matrix);

  inline void Matrix::_assert () const {
    //
    const char funame [] = "Lapack::Matrix::_assert: ";
    
    if(isinit())
      //
      return;
    
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }
  
  inline void Matrix::_assert (int_t i1, int_t i2) const 
  {
    const char funame [] = "Lapack::Matrix::_assert: ";

    _assert();

    if(i1 < 0 || i1 >= size1() || i2 < 0 || i2 >= size2()) {
      //
      std::cerr << funame << "out of range: i1 = " << i1 << ", i2 = " << i2 
		//
		<< ", size1 = " << size1() << ", size2 = " << size2() << "\n";
      
      throw Error::Range();
    }
  }

  inline void Matrix::_assert (Matrix m) const 
  {
    const char funame [] = "Lapack::Matrix::_assert: ";

    _assert();

    m._assert();
    
    if(m.size1() != size1() || m.size2() != size2()) {
      //
      std::cerr << funame << "dimensions mismatch: " << size1() << ", " << size2() << ", " << m.size1() << ", " << m.size2() << "\n";
      
      throw Error::Range();
    }
  }

  // matrix-vector product
  //
  inline void Matrix::_assert (Vector v) const 
  {
    const char funame [] = "Lapack::Matrix::_assert: ";

    _assert();

    v._assert();

    if(size2() != v.size()) {
      //
      std::cerr << funame << "dimensions mismatch: " << size2() << ", " << v.size() << "\n";
      
      throw Error::Range();
    }
  }

  template<typename V>
  inline Matrix::Matrix (const ::Matrix<V>& m)
  {
    const char funame [] = "Lapack::Matrix::Matrix: ";
    
    if(!m.isinit()) {
      //
      std::cerr << funame << "not initialized\n";

      throw Error::Init();
    }

    resize(m.size1(), m.size2());

    RefArr<double>::operator=((const ArrayWithDefault<V>&)m);
  }
  
  template<typename V>
  inline Matrix Matrix::operator= (const ::Matrix<V>& m)
  {
    const char funame [] = "Lapack::Matrix::operator=: ";
    
    if(!m.isinit()) {
      //
      std::cerr << funame << "not initialized\n";

      throw Error::Init();
    }

    if(!isinit())
      //
      resize(m.size1(), m.size2());

    if(size1() != m.size1() || size2() != m.size2()) {
      //
      std::cerr << funame << "dimensions mismatch: " << size1() << ", " << m.size1() << ", " << size2() << ", " << m.size2() << "\n";

      throw Error::Range();
    }

    RefArr<double>::operator=((const ArrayWithDefault<V>&)m);
  }
  
  // square matrix size
  //
  inline int_t Matrix::size () const  
  {
    const char funame [] = "Lapack::Matrix::size: ";

    if(size1() != size2()) {
      //
      std::cerr << funame << "not square\n";
      
      throw Error::Logic();
    }

    return size1();
  }

  inline void Matrix::resize(int_t s1, int_t s2) 
  {
    const char funame [] = "Lapack::Matrix::resize: ";

    if(s1 <= 0 || s2 <= 0) {
      //
      std::cerr << funame << "dimensions out of range: " << s1 << ", " << s2 << "\n";
      
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

    RefArr<double>::resize(s1 * s2);
  }

  // copy constructor by value
  //
  inline Matrix::Matrix (Matrix m, int_t)
  {
    if(m.isinit()) {
      //
      _size1.init(new int_t(m.size1()));
      
      _size2.init(new int_t(m.size2()));
      
      RefArr<double>::operator=(m.RefArr<double>::copy());
    }
  }

  /****************************************************************
   ********************* Symmetric Matrix *************************
   ****************************************************************/

  // packed symmetric matrix with upper triangle reference
  //
  class SymmetricMatrix : private RefArr<double> {
    //
    SharedPointer<int_t> _size;
    
    explicit SymmetricMatrix (SymmetricMatrix, int_t); // copy constructor by value

  public:
    //
    void _assert ()   const;

    void _assert (int_t) const;
    
    void resize (int_t) ;

    bool isinit () const { return _size; }

    SymmetricMatrix () {}
    
    explicit SymmetricMatrix (int_t s)  { resize(s); }
    
    explicit SymmetricMatrix (Matrix, char ='U');

    operator       double* ()       { return RefArr<double>::operator       double*(); }
    operator const double* () const { return RefArr<double>::operator const double*(); }

    SymmetricMatrix copy () const { return SymmetricMatrix(*this, 0); }

    int_t size () const { if(_size) return *_size; return 0; }

    // referencing as an upper triangle column-wise
    //
    double  operator() (int_t, int_t) const;
    double& operator() (int_t, int_t);

    // matrix-vector product
    //
    Vector operator* (Vector)      const;
    Vector operator* (const double*)      const;
    Vector operator* (ConstSlice<double>) const;

    Matrix operator* (Matrix)          const;
    Matrix operator* (SymmetricMatrix) const;
    
    SymmetricMatrix operator+= (SymmetricMatrix m) { _assert(); m._assert(); RefArr<double>::operator+=(m); return *this; }
    SymmetricMatrix operator-= (SymmetricMatrix m) { _assert(); m._assert(); RefArr<double>::operator-=(m); return *this; }

    SymmetricMatrix operator+ (SymmetricMatrix m) const { _assert(); m._assert(); SymmetricMatrix res(*this, 0); res += m; return res; }
    SymmetricMatrix operator- (SymmetricMatrix m) const { _assert(); m._assert(); SymmetricMatrix res(*this, 0); res -= m; return res; }

    SymmetricMatrix operator-  ()  { _assert(); RefArr<double>::operator-(); return *this; }

    SymmetricMatrix operator*= (double d) { _assert();    RefArr<double>::operator*=(d);  return *this; }
    SymmetricMatrix operator/= (double d) { _assert();    RefArr<double>::operator/=(d);  return *this; }

    SymmetricMatrix operator=  (double);
    SymmetricMatrix operator+= (double);
    SymmetricMatrix operator-= (double);

    Vector eigenvalues (Matrix* =0) const;

    SymmetricMatrix invert ()             const;
    SymmetricMatrix positive_invert ()    const;
  };

  inline Vector operator* (ConstSlice<double> v, SymmetricMatrix m) { return m * v; }
  
  inline Vector operator* (const double* v, SymmetricMatrix m) { return m * v; }

  inline void SymmetricMatrix::_assert () const
  {
    const char funame [] = "Lapack::SymmetricMatrix::_assert: ";

    if(isinit())
      //
      return;

    std::cerr << funame << "not iniitalized\n";

    throw Error::Init();
  }
  
  inline void SymmetricMatrix::_assert (int_t i) const
  {
    const char funame [] = "Lapack::SymmetricMatrix::_assert: ";

    _assert();

    if(i < 0 || i >= size()) {
      //
      std::cerr << funame << "out of range: " << i << ", " << size() << "\n";

      throw Error::Init();
    }
  }
  
  inline void SymmetricMatrix::resize(int_t s) 
  {
    const char funame [] = "Lapack::SymmetricMatrix::resize: ";

    if(s <= 0) {
      //
      std::cerr << funame << "non-positive size: " << s << "\n";
      
      throw Error::Range();
    }

    if(!_size) {
      //
      _size.init(new int_t(s));
    }
    else
      //
      *_size = s;
    
    RefArr<double>::resize(s * (s + 1) / 2);
  }

  // copy constructor by value
  //
  inline SymmetricMatrix::SymmetricMatrix (SymmetricMatrix m, int_t)
  {
    if(m.isinit()) {
      //
      _size.init(new int_t(m.size()));
      
      RefArr<double>::operator=(m.RefArr<double>::copy());
    }
  }

  inline double SymmetricMatrix::operator() (int_t i1, int_t i2) const
  {
    _assert(); _assert(i1); _assert(i2);

    if(i1 > i2)
      //
      std::swap(i1, i2);

    return RefArr<double>::operator[](i1 + i2 * (i2 + 1) / 2);
  }

  inline double& SymmetricMatrix::operator() (int_t i1, int_t i2)
  {
    _assert(); _assert(i1); _assert(i2);

    if(i1 > i2)
      //
      std::swap(i1, i2);

    return RefArr<double>::operator[](i1 + i2 * (i2 + 1) / 2);
  }

  inline SymmetricMatrix SymmetricMatrix::operator= (double d)
  {
    _assert();

    RefArr<double>::operator=(0.);
    
    double* p = *this;
    
    int_t n = size() + 2;
    
    for(int_t i = 2; i < n; ++i) {
      //
      *p = d;
      
      p += i;
    }
    
    return *this;
  }

  inline SymmetricMatrix SymmetricMatrix::operator+= (double d)
  {
    _assert();

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
    _assert();

    double* p = *this;
    int_t n = size() + 2;
    for(int_t i = 2; i < n; ++i) {
      *p -= d;
      p += i;
    }
    return *this;
  }

  // Generalized eigenvalue problem
  //
  Vector diagonalize(SymmetricMatrix, SymmetricMatrix, Matrix* = 0) ;

  /****************************************************************
   ******************* Band Symmetric Matrix **********************
   ****************************************************************/

  class BandMatrix : private Matrix {
    //
    void _assert () const ;

    BandMatrix (BandMatrix m, int_t) : Matrix(m, 0) { } // copy construction by value

  public:
    //
    BandMatrix () {}
    
    BandMatrix   (int_t s, int_t b)  : Matrix(b, s) { _assert(); }
    
    void  resize (int_t s, int_t b) { Matrix::resize(b, s); _assert(); }

    bool isinit () const { return Matrix::isinit(); }

    BandMatrix copy () const { return BandMatrix(*this, 0); }

    int_t size      () const  { return size2(); }
    int_t band_size () const  { return size1(); }

    double  operator() (int_t, int_t) const;
    double& operator() (int_t, int_t) ;

    BandMatrix& operator= (double d) { Matrix::operator=(d); return *this; }

    Vector eigenvalues (Matrix* =0) const ;
  };

  inline void BandMatrix::_assert () const 
  {
    const char funame [] = "Lapack::BandMatrix::_assert: ";
    
    if(band_size() <= size())
      //
      return;

    std::cerr << funame << "band size bigger than the matrix size\n";
    
    throw Error::Range();
  }

  /****************************************************************
   ********************** LU Factorization ************************
   ****************************************************************/
 
  int_t parity (RefArr<int_t>);

  class LU : private Matrix {
    Array<int_t> _ipiv;

  public:
    //
    explicit LU (Matrix) ;
    
    int_t size ()            const { return Matrix::size1(); }
    Array<int_t> ipiv ()     const { return _ipiv; }
    double det ()          const;

    Matrix invert ()       const ; // inverse matrix
    Vector invert (Vector) const ; // solve linear equations
    Matrix invert (Matrix) const ; // solve linear equations
  };

  /****************************************************************
   ******* LU Factorization for symmetric packed matrices *********
   ****************************************************************/
  class SymLU : private SymmetricMatrix {
    Array<int_t> _ipiv;

  public:
    //
    explicit SymLU (SymmetricMatrix) ;
    int_t size () const { return SymmetricMatrix::size(); }
    double det () const;

    SymmetricMatrix invert ()              const ; // inverse matrix
    Vector          invert (Vector) const ; // solve linear equations
    Matrix          invert (Matrix) const ; // solve linear equations
  };

  /****************************************************************
   ******************* Cholesky Factorization *********************
   ****************************************************************/

  // for positively defined matrices
  class Cholesky : private SymmetricMatrix {

  public:
    explicit Cholesky (SymmetricMatrix);
    int_t size () const { return SymmetricMatrix::size(); }

    SymmetricMatrix invert ()       const; // inverse matrix
    Vector          invert (Vector) const; // solve linear equations
    Matrix          invert (Matrix) const; // solve linear equations
    
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
      //
      std::cerr << funame << "out of range: i1 = " << i1 << ", i2 = " << i2 
	//
		<< ", size1 = " << size1() << ", size2 = " << size2() << "\n";

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
      //
      std::cerr << funame << "out of range: i1 = " << i1 
		//
		<< ", i2 = " << i2 << ", size = " << size() << "\n";

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
		//
		<< ", i2 = " << i2 << ", size =" << size() << "\n";

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

  ComplexVector fourier_transform (Vector,                                const MultiIndexConvert&);
  ComplexVector fourier_transform (ComplexVector,                         const MultiIndexConvert&);
  ComplexVector fourier_transform (const std::map<int, Lapack::complex>&, const MultiIndexConvert&);

  // find a complimentary orthogonal basis set to the non-orthogonal vector set
  //
  double orthogonalize (Matrix basis, int_t vsize);

  // solve linear equations by svd-decomposition
  //
  Vector svd_solve (Matrix a, Vector b, double pres = -1., double (*weight)(double) = 0);
  //
  Matrix svd_solve (Matrix a, Matrix b, int_t& rank, double pres = -1.);

  // singular value decomposition, A = U * S * V**T
  //
  void svd (Matrix a, Vector s, Matrix u, Matrix v);
  //
}// namespace Lapack

#endif
