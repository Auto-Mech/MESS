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

#ifndef D3_HH
#define D3_HH

#include "linpack.hh"
#include "array.hh"

/******************************************
 ************ 3-D Space  *****************
 ******************************************/

class Quaternion;

namespace D3 {
  class Matrix;

  class Vector
  {
    double _begin [3];
    double* _end;

  public:
    typedef       double*       iterator;
    typedef const double* const_iterator;
    typedef       double      value_type;
    
    Vector ();
    Vector (const Vector&);
    Vector& operator= (const Vector&);

    explicit Vector (double);
    explicit Vector (const double*);
    template <typename V>
    explicit Vector (const V&) ;

    double& operator [] (int i)       { return _begin[i]; }
    double  operator [] (int i) const { return _begin[i]; }

    operator       double* ()       { return _begin; }
    operator const double* () const { return _begin; }

    const double* begin () const { return _begin; }
    double*       begin ()       { return _begin; }
    const double*   end () const { return _end; }
    double*         end ()       { return _end; }
    
    int size () const { return 3; }

    // block operations
    template <typename V> Vector& operator=  (const V&) ;
    template <typename V> Vector& operator+= (const V&) ;
    template <typename V> Vector& operator-= (const V&) ;

    Vector& operator=  (const double*);
    Vector& operator+= (const double*);
    Vector& operator-= (const double*);

    Vector operator+ (const double*)   const;
    Vector operator+ (const Vector& v) const { return *this + (const double*)v; }
    Vector operator- (const double*)   const;
    Vector operator- (const Vector& v) const { return *this - (const double*)v; }
    Vector operator* (double)          const;
    Vector operator* (const Matrix&)   const;

    Vector& operator=  (double);
    Vector& operator*= (double);
    Vector& operator/= (double);

    Vector& operator*= (const Matrix&); // orthogonal transformation M*v

    double vlength () const;

    double normalize ();
    double orthogonalize (const double*) ;
  };

  template <typename V>
  Vector::Vector (const V& v) 
  {
    const char funame [] = "D3::Vector::Vector: ";

    if(v.size() != 3) {
      std::cerr << funame << "dimensions mismatch\n";
      throw Error::Range();
    }

    _end = _begin + 3;

    typename V::const_iterator vit = v.begin();
    for(iterator it = begin(); it != end(); ++it, ++vit)
      *it = *vit;
  }

  template <typename V>
  Vector& Vector::operator= (const V& v) 
  {
    const char funame [] = "D3::Vector::operator=: ";

    if(v.size() != 3) {
      std::cerr << funame << "dimensions mismatch\n";
      throw Error::Range();
    }

    typename V::const_iterator vit = v.begin();
    for(iterator it = begin(); it != end(); ++it, ++vit)
      *it = *vit;

    return *this;
  }

  template <typename V>
  Vector& Vector::operator+= (const V& v) 
  {
    const char funame [] = "D3::Vector::operator+=: ";

    if(v.size() != 3) {
      std::cerr << funame << "dimensions mismatch\n";
      throw Error::Range();
    }

    typename V::const_iterator vit = v.begin();
    for(iterator it = begin(); it != end(); ++it, ++vit)
      *it += *vit;

    return *this;
  }

  template <typename V>
  Vector& Vector::operator-= (const V& v) 
  {
    const char funame [] = "D3::Vector::operator-=: ";

    if(v.size() != 3) {
      std::cerr << funame << "dimensions mismatch\n";
      throw Error::Range();
    }

    typename V::const_iterator vit = v.begin();
    for(iterator it = begin(); it != end(); ++it, ++vit)
      *it -= *vit;

    return *this;
  }

  inline Vector Vector::operator+ (const double* p) const
  {
    Vector res = *this;
    res += p;
    return res;
  }

  inline Vector Vector::operator- (const double* p) const
  {
    Vector res = *this;
    res -= p;
    return res;
  }

  inline Vector Vector::operator* (double d) const
  {
    Vector res = *this;
    res *= d;
    return res;
  }

  inline double Vector::vlength () const 
  {
    return ::vlength(_begin, 3);
  }

  inline double Vector::normalize () 
  {
    return ::normalize(_begin, 3);
  }

  inline double Vector::orthogonalize (const double* n) 
  {
    return ::orthogonalize(_begin, n, 3);
  }

  Vector operator+ (const double*, const Vector&);
  Vector operator- (const double*, const Vector&);
  Vector operator* (double,        const Vector&);

  std::istream& operator>> (std::istream&,       Vector&);
  std::ostream& operator<< (std::ostream&, const Vector&);

  Vector   vprod (const double*, const double*);          // vector-vector product
  void     vprod (const double*, const double*, double*); // vector-vector product
  void     vprod (const Matrix&, const double*, double*); // matrix-vector product
  void     vprod (const double*, const Matrix&, double*); // vector-matrix product

  double volume (const double*, const double*, const double*); // signed volume

  void find_orth(const double*, Vector*); // get two vectors orthogonal to the given one

  /******************************************* 3-D matrix *****************************************/

  // C style matrix
  class Matrix 
  {
    double _begin [9];
    void   _init ();

  public:
    Matrix () { _init(); }
    explicit Matrix (double a) { _init(); diagonal() = a; }
    Matrix (Vector, Vector) ;               // standard orientation
    explicit Matrix (const Quaternion&) ;

    double&        operator() (int, int)       ;   // C style indexing
    const double&  operator() (int, int) const ;   

    operator Quaternion () const ;

    Matrix& operator= (double a) { _init(); diagonal() = a; return *this; }

    Slice<double>      column (int i)       ;
    ConstSlice<double> column (int i) const ;
    Slice<double>      row    (int i)       ;
    ConstSlice<double> row    (int i) const ;

    Slice<double>          diagonal ();
    ConstSlice<double>     diagonal () const;

    Matrix transpose () const;
    Matrix operator- () const;

    Matrix& operator *= (double);
    Matrix& operator /= (double a) { return operator*=(1. / a); }
    
    
    Matrix operator* (const Matrix&) const; // matrix multiplication
    Vector operator* (const double*) const; // matrix vector product M*v
    Vector operator* (const Vector& v) const { return *this * (const double*)v; }

    void orthogonalize ();
    void orthogonality_check (double = -1.) const ;

    int inversion () const;
    void     chirality_check () const ;
  };

  inline Slice<double> Matrix::diagonal () 
  {
    return Slice<double>(_begin, 3, 4);
  }

  inline ConstSlice<double> Matrix::diagonal () const
  {
    return ConstSlice<double>(_begin, 3, 4);
  }

  class Plane
  {
    Vector _normal;
    double _dist;
    Vector _orth [2];
  public:
	
    void set (const Vector& n, double d);

    Plane ();
    Plane(const Vector& n, double d) { set(n, d); }

    const Vector& normal ()   const { return _normal; }
    double dist ()            const { return _dist; }
    const Vector& orth(int i) const { return _orth[i]; }
  };

  inline std::ostream& operator<< (std::ostream& to, const Plane& p)
  {
    to << "{normal = " << p.normal() << ", offset = " << p.dist() << "}";
    return to;
  }
   

  // reference frame
  struct Frame
  {
    Vector origin; // displacement
    Matrix orient; // orientation matrix

    void set(const Vector& v1, const Vector& v2, const Vector& v3) 
    { origin = v1; orient = Matrix(v2 - v1, v3 - v1); }

    Frame () : origin(0.), orient(1.) { }
    Frame (const Vector& v1, const Vector& v2, const Vector& v3) { set(v1, v2, v3); }

    void fv2lv (const double*, double*) const; // frame vector to lab vector
    void lv2fv (const double*, double*) const; // lab vector to frame vector
    void fp2lp (const double*, double*) const; // frame position to lab position
    void lp2fp (const double*, double*) const; // lab position to frame position
  };
} // D3 namespace

/*************************************************************************
 *************************** QUATERNION **********************************
 *************************************************************************/

class Quaternion : private Array<double> {

public:

  enum { NOCHECK = 1 };
  static double tolerance;

  typedef const double* const_iterator;
  typedef       double*       iterator;
  typedef       double      value_type;

  Quaternion () : Array<double>(4, 0.) {}

  explicit Quaternion (double d)        : Array<double>(4, 0.) { *begin() = d; }
  explicit Quaternion (const double* p) : Array<double>(4) { Array<double>::operator=(p); }
  explicit Quaternion (const D3::Matrix&, int =0) ;

  template<class V>
  explicit Quaternion (const V&) ;

  int size () const { return Array<double>::size(); }

  operator       double* ()       { return Array<double>::begin(); }
  operator const double* () const { return Array<double>::begin(); }

  operator D3::Matrix () const ;

  double*       begin ()       { return Array<double>::begin(); }
  const double* begin () const { return Array<double>::begin(); }
  double*         end ()       { return Array<double>::end(); }
  const double*   end () const { return Array<double>::end(); }

  double&       operator[] (int i)       { return Array<double>::operator[](i); }
  const double& operator[] (int i) const { return Array<double>::operator[](i); }

  //vector operations
  template<class V> Quaternion& operator=  (const V&) ;
  template<class V> Quaternion& operator+= (const V& v) { Array<double>::operator+=(v); return *this; }
  template<class V> Quaternion& operator-= (const V& v) { Array<double>::operator-=(v); return *this; }
  
  Quaternion& operator=   (const double* p) { Array<double>::operator= (p); return *this; } 
  Quaternion& operator+=  (const double* p) { Array<double>::operator+=(p); return *this; } 
  Quaternion& operator-=  (const double* p) { Array<double>::operator-=(p); return *this; } 
  
  Quaternion& operator*=  (double d) { Array<double>::operator*=(d); return *this; } 
  Quaternion& operator/=  (double d) { Array<double>::operator/=(d); return *this; } 

  //quaternion operations
  Quaternion operator* (const double*) const;
  Quaternion operator/ (const double*) const;

  Quaternion operator+ (const Quaternion& q) const { Quaternion res(*this); res += q; return res; }
  Quaternion operator- (const Quaternion& q) const { Quaternion res(*this); res -= q; return res; }

  void normalize ();

  static void qprod (const double*, const double*, double*);
};

// quaternion-to-matrix transformation
void quat2mat (const double* q, D3::Matrix& m) ;

// (orthogonal) matrix-to-quaternion transformation
void mat2quat (const D3::Matrix& mat, double* q, int =0) ;

inline D3::Matrix::operator Quaternion () const 
{
  Quaternion res;
  mat2quat(*this, res);

  return res;
}

inline D3::Matrix::Matrix (const Quaternion& q) 
{
  quat2mat(q, *this);
}

inline Quaternion::operator D3::Matrix () const  
{
  D3::Matrix res;
  quat2mat(*this, res);
  return res;
}

inline Quaternion::Quaternion (const D3::Matrix& m, int flags) 
  : Array<double>(4)
{ 
  mat2quat(m, *this, flags);
}

inline Quaternion Quaternion::operator* (const double* q) const
{
  Quaternion res;

  qprod(*this, q, res);

  return res;
}

inline Quaternion Quaternion::operator/ (const double* q) const
{
  Quaternion qtemp(q);
  for(iterator it = qtemp.begin() + 1; it != qtemp.end(); ++it)
    *it = - *it;

  Quaternion res;
  qprod(*this, qtemp, res);
  res /= vdot(qtemp);

  return res;
}

template<class V>
Quaternion::Quaternion (const V& v)  : Array<double>(v)
{
  const char funame [] = "Quaternion::Quaternion: ";

  if(v.size() != 4) {
    std::cerr << funame << "wrong initializer dimension (" << v.size() << ")/n";
    throw Error::Range();
  }  
}

template<class V>
Quaternion& Quaternion::operator= (const V& v) 
{
  const char funame [] = "Quaternion::operator=: ";

  if(v.size() != 4) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Range();
  } 

  Array<double>::operator=(v);

  return *this;
}

void axis_quaternion (int, double, double*);

void polar2cart (const double*, double*);

#endif
