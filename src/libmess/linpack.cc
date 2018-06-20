/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include "linpack.hh"
#include <iostream>
#include <cmath>

/*********************** Vector operations *******************************/

void multiply (double* v, double val, int size, int step)
{
  for (double* end = v + step * size; v != end; v += step)
    *v *= val;
}

double normalize (double* v, int size, int step) throw(Error::General)
{
  const char funame [] = "normalize: ";

  double norm = vlength(v, size, step);

  if (norm == 0.) {
#ifdef DEBUG
    std::cerr << funame << "WARNING: the vector length is zero\n";
#endif
    return 0.;
  }

  for (double* end = v + step * size; v != end; v += step)
    *v /= norm;

  return norm;
}

double orthogonalize (double* vstart, const double* nstart, int size, int vstep, int nstep) throw(Error::General)
{
  const char funame [] = "orthogonalize: ";
  
  static const double tol = 1.e-12;
  
  if(vstep <= 0 || nstep <= 0) {
    std::cerr << funame << "step out of range\n";
    throw Error::Range();
  }

  if(size < 0) {
    std::cerr << funame << "size out of range\n";
    throw Error::Range();
  }

  double*       vtemp;
  const double* ntemp;
  double*       vend = vstart + size * vstep;

  double nlen = 0.0, proj = 0.0, vlen = 0.0;
  for (vtemp = vstart, ntemp = nstart; vtemp != vend; vtemp += vstep, ntemp += nstep) {
    nlen += *ntemp * *ntemp;
    vlen += *vtemp * *vtemp;
    proj += *ntemp * *vtemp;
  }

  if(nlen == 0.) {
    std::cerr << funame << "the ortogonalizing vector length is zero\n";
    throw Error::Range();
  }

  if(vlen == 0.) {
#ifdef DEBUG
    std::cerr << funame << "WARNING: the ortogonalized vector length is zero\n";
#endif
    return 0.;
  }

  proj /= nlen;
  for (vtemp = vstart, ntemp = nstart; vtemp != vend; vtemp += vstep, ntemp += nstep)
    *vtemp -= proj * *ntemp;

  double cos2 = proj * proj * nlen / vlen;


  return cos2;
}

double parallel_orthogonalize (double* v, const double* n, int size, int vstep, int nstep)
  throw(Error::General)
{
  const char funame [] = "parallel_orthogonalize: ";
  
  static const double tol = 1.e-12;
  
  double vval, nval;

  if(vstep <= 0 || nstep <= 0) {
    std::cerr << funame << "step out of range\n";
    throw Error::Range();
  }

  if(size < 0) {
    std::cerr << funame << "size out of range\n";
    throw Error::Range();
  }

  double nlen = 0.0, proj = 0.0, vlen = 0.0;

#pragma omp parallel for default(shared) reduction(+: nlen, proj, vlen) schedule(static)
	
  for(int i = 0; i < size; ++i) {
    double nval = n[i * nstep];
    double vval = v[i * vstep];

    nlen += nval * nval;
    vlen += vval * vval;
    proj += nval * vval;
  }

  if(nlen == 0.) {
    std::cerr << funame << "the ortogonalizing vector length is zero\n";
    throw Error::Range();
  }

  if(vlen == 0.) {
#ifdef DEBUG
    std::cerr << funame << "WARNING: the ortogonalized vector length is zero\n";
#endif
    return 0.;
  }

  proj /= nlen;

#pragma omp parallel for default(shared) schedule(static)
	
  for(int i = 0; i < size; ++i) {
    v[i * vstep] -= proj * n[i * nstep];
  }

  double cos2 = proj * proj * nlen / vlen;

#ifdef DEBUG

  if(cos2 > 1. - tol) {
    std::cerr << funame << "vectors are colinear, sin**2 = " << 1. - cos2 << "\n";
    throw Error::Range();
  }

#endif

  return cos2;
}

double vdistance (const double* v1, const double* v2, int size, int step1, int step2) throw(Error::General)
{
  const char funame [] = "vdistance: ";

  if(step1 <= 0 || step2 <= 0) {
    std::cerr << funame << "step out of range\n";
    throw Error::Range();
  }

  if(size < 0) {
    std::cerr << funame << "size out of range\n";
    throw Error::Range();
  }

  double dtemp;

  double res = 0.;

  for (const double* end = v1 + step1 * size; v1 != end; v1 += step1, v2 += step2) {
    dtemp = *v1 - *v2;
    res += dtemp * dtemp; 
  }

  return std::sqrt(res);
}

double vdot (const double* v1, const double* v2, int size, int step1, int step2) throw(Error::General)
{
  const char funame [] = "vdot: ";

  if(step1 <= 0 || step2 <= 0) {
    std::cerr << funame << "step out of range\n";
    throw Error::Range();
  }

  if(size < 0) {
    std::cerr << funame << "size out of range\n";
    throw Error::Range();
  }

  double res = 0.;

  for (const double* end = v1 + step1 * size; v1 != end; v1 += step1, v2 += step2)
    res += *v1 * *v2;
  

  return res;
}

double parallel_vdot (const double* v1, const double* v2, int size, int step1, int step2) throw(Error::General)
{
  const char funame [] = "vdot: ";

  if(step1 <= 0 || step2 <= 0) {
    std::cerr << funame << "step out of range\n";
    throw Error::Range();
  }

  if(size < 0) {
    std::cerr << funame << "size out of range\n";
    throw Error::Range();
  }

  double res = 0.;

#pragma omp parallel for default(shared) reduction(+: res) schedule(static)
	
  for (int i = 0; i < size; ++i) {
    res += v1[i * step1] * v2[i * step2];
  }

  return res;
}

double vdot (const double* v, int size, int step) throw(Error::General)
{
  const char funame [] = "vdot: ";

  if(step <= 0) {
    std::cerr << funame << "step out of range\n";
    throw Error::Range();
  }

  if(size < 0) {
    std::cerr << funame << "size out of range\n";
    throw Error::Range();
  }

  double res = 0.;

  for (const double* end = v + step * size; v != end; v += step)
    res += *v * *v;

  return res;
}

double vlength (const double* v, int size, int step) throw(Error::General)
{
  return std::sqrt(vdot(v, size, step));
}
