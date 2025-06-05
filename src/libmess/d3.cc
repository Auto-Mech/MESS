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

#include "d3.hh"
#include <iostream>
#include <cmath>

double Quaternion::tolerance = 1.e-12;

// transformation from euler angles in appropriate convention to quaternion
//
void euler2quat (const double* euler, double* quat, char conv)
{
  const char funame [] = "euler2quat: ";

  double theta = euler[0], phi, psi;

  switch(conv) {
    //
  case 'X':
    //
    phi   = euler [1];
    psi   = euler [2];

    break;
    //
  case 'Y':
    //
    phi   = euler [1] + M_PI_2;
    psi   = euler [2] - M_PI_2;

    break;
    //
  default:
    //
    std::cerr << funame << "unknown convention: " << conv << "\n";

    throw Error::Range();
  }

   double s = (phi + psi) / 2.0, d = (phi - psi) / 2.0, t = theta / 2.;

   quat[0] = std::cos(t) * std::cos(s);
   quat[1] = std::sin(t) * std::cos(d);
   quat[2] = std::sin(t) * std::sin(d);
   quat[3] = std::cos(t) * std::sin(s);
}

/*********************** 3-D reference frame ***********************/

void D3::Frame::fv2lv (const double* fv, double* lv) const // frame vector to lab vector
{
  vprod(fv, orient, lv);
}

void D3::Frame::fp2lp (const double* fp, double* lp) const // frame to lab
{
  vprod(fp, orient, lp);
  for(int i = 0; i < 3; ++i)
    lp[i] += origin[i];
}

void D3::Frame::lv2fv (const double* lv, double* fv) const // lab to frame
{
  vprod(orient, lv, fv);
}

void D3::Frame::lp2fp (const double* lp, double* fp) const // lab to frame
{
  vprod(orient, lp - origin, fp);
}


/*********************** 3-D plane ***********************/

void D3::Plane::set (const Vector& n, double d)
{
  _normal = n;
  _dist = d / _normal.normalize();
  find_orth(_normal, _orth);
}

D3::Plane::Plane () 
{
  _normal = 0.;
  _normal[0] = 1.;
  _dist = 1.;
  find_orth(_normal, _orth);
}

/************************************* 3-D rotational matrix *******************************/

D3::Matrix D3::Matrix::transpose () const
{
  Matrix res;
  for(int i = 0; i < 3; ++i)
    res.row(i) = column(i);

  return res;
}

D3::Matrix::Matrix (Vector n1, Vector n2)  
{
  n1.normalize();
  row(0) = n1;

  n2.orthogonalize(n1);
  n2.normalize();  
  row(1) = n2;

  row(2) = vprod(n1, n2);  
}

void D3::Matrix::_init ()
{
  double* end = _begin + 9;
  
  for(double* it = _begin; it != end; ++it)
    *it = 0.;
}

D3::Matrix& D3::Matrix::operator *= (double a)
{
  double* end = _begin + 9;
  
  for(double* it = _begin; it != end; ++it)
    *it *= a;
  
  return *this;
}

D3::Matrix D3::Matrix::operator- () const
{
  D3::Matrix res;

  const double* end = _begin + 9;
  
  double* rit = res._begin;
  for(const double* it = _begin; it != end; ++it)
    *rit++ = - *it;
  
  return res;
}

double& D3::Matrix::operator() (int i, int j) 
{
  const char funame [] = "D3::Matrix::operator() (int, int):";

#ifdef DEBUG

  if(i < 0 || i > 2 || j < 0 || j > 2) {
    std::cerr << funame << "indices out of range\n";
    throw Error::Range();
  }

#endif

  return *(_begin + (i * 3 +  j));
}

const double& D3::Matrix::operator() (int i, int j) const 
{
  const char funame [] = "D3::Matrix::operator() (int, int) const:";

#ifdef DEBUG

  if(i < 0 || i > 2 || j < 0 || j > 2) {
    std::cerr << funame << "indices out of range\n";
    throw Error::Range();
  }

#endif

  return *(_begin + (i * 3 +  j));
}

Slice<double> D3::Matrix::column (int i)  
{
  const char funame [] = "D3::Matrix::column (int):";

#ifdef DEBUG

  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

#endif

  return Slice<double>(_begin + i, 3, 3);
}

ConstSlice<double> D3::Matrix::column (int i) const 
{
  const char funame [] = "D3::Matrix::column (int) const:";

#ifdef DEBUG

  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

#endif

  return ConstSlice<double>(_begin + i, 3, 3);
}

Slice<double> D3::Matrix::row (int i) 
{
  const char funame [] = "D3::Matrix::row (int):";

#ifdef DEBUG

  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

#endif

  return Slice<double>(_begin + i * 3, 3);

}

ConstSlice<double> D3::Matrix::row (int i) const 
{
  const char funame [] = "D3::Matrix::row (int) const:";

#ifdef DEBUG

  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

#endif

  return ConstSlice<double>(_begin + i * 3, 3);
}

D3::Matrix D3::Matrix::operator* (const Matrix& r) const
{
  Matrix res;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      res(i, j) = vdot(row(i), r.column(j));

  return res;
}

D3::Vector D3::Matrix::operator* (const double* v) const
{
  Vector res;
  for(int i = 0; i < 3; ++i)
    res[i] = vdot(row(i), v);

  return res;
}

void D3::Matrix::orthogonality_check (double toll) const 
{
  const char funame [] = "D3::Matrix::orthogonality_check: ";

  double dtemp;

  if(toll <= 0.)
    toll = Quaternion::tolerance;

  for(int i = 0; i < 3; ++i)
    for(int j = i; j < 3; ++j) {
      dtemp = vdot(row(i), row(j));
      if(i == j)
	dtemp -= 1.;
      if(dtemp < -toll || dtemp > toll) {
	std::cerr << funame << "failed\n";
	throw Error::Range();
      }
    } 
}

void D3::Matrix::orthogonalize ()
{
  const char funame [] = "D3::Matrix::orthogonalize: ";

  double dtemp;

  for(int i = 0; i < 3; ++i) {
    Slice<double> x = row(i);
    for(int j = 0; j < i; ++j) {
      dtemp = -vdot(x, row(j));
      x.add(dtemp, row(j));
    }
    normalize(x);
  } 
}

int D3::Matrix::inversion () const
{
  const char funame [] = "D3::Matrix::inversion: ";

  if(volume(_begin, _begin + 3, _begin + 6) < 0.)
    return 1;
  
  return 0;
}

void D3::Matrix::chirality_check () const 
{
  const char funame [] = "D3::Matrix::chirality_check: ";

  if(inversion()) {
    std::cerr << funame << "failed\n";
    throw Error::Range();
  }
}

/******************************************
 ************ 3-D Real Vector *************
 ******************************************/

D3::Vector::Vector () 
{
  _end = _begin + 3;

  for(iterator it = begin(); it != end(); ++it)
    *it = 0.;
}

D3::Vector::Vector(const Vector& v)
{
  _end = _begin + 3;

  const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
}

D3::Vector& D3::Vector::operator= (const Vector& v)
{
  const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;

  return *this;
}

D3::Vector::Vector (double a)
{
  _end = _begin + 3;

  for(iterator it = begin(); it != end(); ++it)
    *it = a;
}

D3::Vector::Vector (const double* p)
{
  _end = _begin + 3;

  for(iterator it = begin(); it != end(); ++it)
    *it = *p++;
}

D3::Vector& D3::Vector::operator= (const double* p)
{
  for(iterator it = begin(); it != end(); ++it)
    *it = *p++;

  return *this;
}

D3::Vector& D3::Vector::operator+= (const double* p)
{
  for(iterator it = begin(); it != end(); ++it)
    *it += *p++;

  return *this;
}

D3::Vector& D3::Vector::operator-= (const double* p)
{
  for(iterator it = begin(); it != end(); ++it)
    *it -= *p++;

  return *this;
}

D3::Vector& D3::Vector::operator= (double val)
{
  for(iterator it = begin(); it != end(); ++it)
    *it = val;

  return *this;
}

D3::Vector& D3::Vector::operator*= (double val)
{
  for(iterator it = begin(); it != end(); ++it)
    *it *= val;

  return *this;
}

D3::Vector& D3::Vector::operator/= (double val)
{
  for(iterator it = begin(); it != end(); ++it)
    *it /= val;

  return *this;
}

D3::Vector& D3::Vector::operator*= (const Matrix& m) //rotation
{
  *this = m * (*this);
  return *this;
}

D3::Vector D3::Vector::operator* (const Matrix& m) const
{
  Vector res;
  vprod(*this, m, res);
  return res;
}

D3::Vector D3::operator* (double d, const Vector& a)
{
  return a * d;
}

D3::Vector D3::operator+ (const double* b, const D3::Vector& a)
{
  return a + b;
}

D3::Vector D3::operator- (const double* p, const D3::Vector& v)
{
  Vector res;
  Vector::const_iterator vit = v.begin();
  for(Vector::iterator rit = res.begin(); rit != res.end(); ++rit, ++vit, ++p)
    *rit = *p - *vit;

  return res;
}

std::ostream& D3::operator<< (std::ostream& to, const Vector& v)
{
  int old_precision = to.precision(3);

  to << "{";
  for(int i = 0; i < 3; ++i) {
    if(i)
      to << ", ";
    to << v[i];
  }
  to << "}";

  to.precision(old_precision);

  return to;
}

std::istream& D3::operator>> (std::istream& from, Vector& v)
{
  for(int i = 0; i < 3; ++i)
    from >> v[i];
  return from;
}


D3::Vector D3::vprod (const double* a1, const double* a2)
{
  D3::Vector res;
  vprod(a1, a2, (double*)res);
  return res;
}

void D3::vprod (const double* a1, const double* a2, double* res)
{
  int i1, i2;
  for(int i = 0; i < 3; ++i) {
    i1 = (i + 1) % 3;
    i2 = (i + 2) % 3;
    res[i] = a1[i1] * a2[i2] - a1[i2] * a2[i1];
  }
}

void D3::vprod (const D3::Matrix& r, const double* v, double* res)
{
  for(int i = 0; i < 3; ++i)
    res[i] = vdot(r.row(i), v);
}

void D3::vprod (const double* v, const D3::Matrix& r, double* res)
{
  for(int i = 0; i < 3; ++i)
    res[i] = vdot(r.column(i), v);
}

void D3::find_orth (const double* n0, D3::Vector* ort) // get two vectors orthogonal to the given one
{
  int imin;
  double dtemp, dmin;
  for(int i = 0; i < 3; ++i) {
    dtemp = n0[i] > 0. ? n0[i] : -n0[i];
    if(!i || dtemp < dmin) {
      imin = i;
      dmin = dtemp;
    }
  }
  
  ort[0] = 0.;
  ort[0][imin] = 1.;

  ort[0].orthogonalize(n0);
  ort[0].normalize();

  vprod(n0, ort[0], ort[1]);
  ort[1].normalize();
}

double D3::volume (const double* a, const double* a1, const double* a2)
{
  double res = 0.;

  for(int i = 0; i < 3; ++i) {
    int i1 = (i + 1) % 3;
    int i2 = (i + 2) % 3;

    res += a[i] * (a1[i1] * a2[i2] - a1[i2] * a2[i1]);
  }

  return res;
}

/*************************************************************************
 *************************** Quaternions *********************************
 *************************************************************************/


/********************************************************************
*
template <typename T>
void quat2mat (const double* quat, T& mat)
{// normalized quaternion-to-matrix transformation
 
  register     double    q0, q1, q2, q3, qq0, qq1, qq2, qq3;
  register     double    a01, a02, a03, a12, a13, a23;
  
  q0 = quat[0];  q1 = quat[1]; q2 = quat[2];  q3 = quat[3];

  qq0 = q0 * q0;  qq1 = q1 * q1;  qq2 = q2 * q2;  qq3 = q3 * q3;
 
  a01 = 2.0 * q0 * q1;  a02 = 2.0 * q0 * q2;  a03 = 2.0 * q0 * q3;
  a12 = 2.0 * q1 * q2;  a13 = 2.0 * q1 * q3;
  a23 = 2.0 * q2 * q3;

  mat (0, 1) = a12 + a03;
  mat (1, 0) = a12 - a03;
  mat (0, 2) = a13 - a02;
  mat (2, 0) = a13 + a02;
  mat (1, 2) = a23 + a01;
  mat (2, 1) = a23 - a01;
  
  mat (0, 0) = (qq0 + qq1 - qq2 - qq3);
  mat (1, 1) = (qq0 - qq1 + qq2 - qq3);
  mat (2, 2) = (qq0 - qq1 - qq2 + qq3);
}
*****************************************************************/

// quaternion-to-matrix transformation
//
void quat2mat (const double* q, D3::Matrix& m) 
{
  const char funame [] = "quat2mat: ";

  double qq = vdot(q, 4);

  if(qq == 0.) {
    //
    std::cerr << funame << "quaternion length is 0\n";

    throw Error::Range();
  }

  m = 0.;

  double qterm;

  for(int i = 0; i < 4; ++i)
    //
    for(int j = i; j < 4; ++j) {
      //
      qterm = q[i] * q[j] / qq;

      // q_0 * q_0 terms
      //
      if(!j) {
	//
	m.diagonal() += qterm;
      }
      // q_i  * q_i terms
      //
      else if(i == j) {
	//
	for(int k = 0; k < 3; ++k)
	  //
	  if(k + 1 == i) {
	    //
	    m(k, k) += qterm;
	  }
	  else
	    //
	    m(k, k) -= qterm;
      }
      // q_0 * q_i terms
      //
      else if(!i) {
	//
	int j1 = j % 3;
	int j2 = (j + 1) % 3;

	qterm *= 2.;
	
	m(j1, j2) += qterm;
	m(j2, j1) -= qterm;
      }
      // q_i * q_j terms
      //
      else {
	//
	int i1 = i - 1;
	int j1 = j - 1;

	qterm *= 2.;

	m(i1, j1) += qterm;
	m(j1, i1) += qterm;
      }
    }
}

// (orthogonal) matrix-to-quaternion transformation
//
void mat2quat (const D3::Matrix& mat, double* quat, int flags) 
{
  const char funame [] = "mat2quat: ";

  double dtemp;
  int    itemp;
  
  if(!(flags & Quaternion::NOCHECK)) {
    mat.orthogonality_check();
    mat.chirality_check();
  }


  const double qq = (sum(mat.diagonal()) + 1.) / 4.;

  int i1, i2;
  if(qq > 0.25) {
    *quat = std::sqrt(qq);

    for(int i = 1; i <= 3; ++i) {
      i1 = i % 3;
      i2 = (i + 1) % 3;
      quat[i] = (mat(i1, i2) - mat(i2, i1)) / 4. / *quat;
    }

    return;
  }

  double dmax;
  int    imax;
  for(int i = 1; i <= 3; ++i) {
    i1 = i % 3;
    i2 = (i + 1) % 3;
    dtemp = qq - (mat(i1, i1) + mat(i2, i2)) / 2.;

    if(i == 1 || dtemp > dmax) {
      dmax = dtemp;
      imax = i;
    }
  }

  dmax = std::sqrt(dmax);
 
  for(int i = 1; i <= 3; ++i) 
    if(i != imax)
      quat[i] = (mat(i - 1, imax - 1) + mat(imax - 1, i - 1)) / 4. / dmax;
    else
      quat[i] = dmax;

  i1 = imax % 3;
  i2 = (imax + 1) % 3;
    
  *quat = (mat(i1, i2) - mat(i2, i1)) / 4. / dmax;
}

// quaternion-quaternion product, q1 * q2
//
void Quaternion::qprod (const double* q1, const double* q2, double* q) 
{

  const double* v1 = q1 + 1;
  const double* v2 = q2 + 1;
  double*       v  = q  + 1;

  D3::vprod(v1, v2, v);

  for(int i = 0; i < 3; ++i)
    //
    v[i] += *q1 * v2[i] + *q2 * v1[i];

  *q = *q1 * *q2 - vdot(v1, v2, 3);
}

void axis_quaternion (int axis, double angle, double* q)  // qvq* rotates vector v about axis x on angle a
{
  q[0] = std::cos(angle / 2.);

  for(int i = 1; i < 4; ++i)
    if(i == axis + 1)
      q[i] = std::sin(angle / 2.);
    else
      q[i] = 0.;
}

void polar2cart (const double* polar_ang, double* n)
{// polar_angles-to-directional_cosines transformation

   const double& theta = polar_ang[0];
   const double& phi   = polar_ang[1];

   const double sin_theta = std::sin(theta);

   n[0] = sin_theta * std::cos(phi);
   n[1] = sin_theta * std::sin(phi);
   n[2] = std::cos(theta);
}
