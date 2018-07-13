

#ifndef LINPACK_HH
#define LINPACK_HH

#include "error.hh"
#include <iostream>
#include <cmath>

/********************** General Vector Operations ***************************/

double normalize (double*, int, int =1) throw(Error::General);
double orthogonalize (double*, const double*, int, int =1, int =1) throw(Error::General);
double parallel_orthogonalize (double*, const double*, int, int =1, int =1) throw(Error::General);
double vdot      (const double*, int, int =1) throw(Error::General);
double vlength   (const double*, int, int =1) throw(Error::General);
double vdot      (const double*, const double*, int, int =1, int =1) throw(Error::General);
double parallel_vdot  (const double*, const double*, int, int =1, int =1) throw(Error::General);
double vdistance (const double*, const double*, int, int =1, int =1) throw(Error::General);

void multiply (double*, double, int, int= 1);

template <typename V>
typename V::value_type vdot (const V& v)
{
  typename V::value_type res = 0;
  for(typename V::const_iterator it = v.begin(); it != v.end(); ++it)
    res += *it * *it;

  return res;
}

template <typename V>
typename V::value_type vdot (const V& v1, const V& v2) throw(Error::General)
{
  const char funame [] = "vdot: ";
  
  if(v1.size() != v2.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }
  
  typename V::value_type res = 0;
  typename V::const_iterator i1, i2;
  for(i1 = v1.begin(), i2 = v2.begin(); i1 != v1.end(); ++i1, ++i2)
    res += *i1 * *i2;

  return res;
}

template <typename V>
typename V::value_type vdot (const V& v, const typename V::value_type* vit, int stride = 1);

template <typename V>
typename V::value_type vdot (const V& v, const typename V::value_type* vit, int stride)
{
  typename V::value_type res = 0;
  
  for(typename V::const_iterator it = v.begin(); it != v.end(); ++it, vit += stride)
    res += *it * *vit;

  return res;
}

template <typename V>
typename V::value_type vlength (const V& v)
{
  return std::sqrt(vdot(v));
}

template <typename V>
typename V::value_type normalize (V& v)
{
  const char funame [] = "normalize: ";

  typename V::value_type norm = vdot(v);

  if (norm == 0) {
#ifdef DEBUG
    std::cerr << funame << "WARNING: the vector length is zero\n";
#endif
    return norm;
  }

  norm = std::sqrt(norm);

  for (typename V::iterator it = v.begin(); it != v.end(); ++it)
    *it /= norm;

  return norm;
}

template <typename V>
typename V::value_type vdistance (const V& v1, const V& v2) throw(Error::General)
{
  const char funame [] = "vdistance: ";
  
  if(v1.size() != v2.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }
  
  typename V::value_type res = 0;
  typename V::value_type   vtemp;
  
  typename V::const_iterator i2 = v2.begin();
  for(typename V::const_iterator i1 = v1.begin(); i1 != v1.end(); ++i1, ++i2) {
    vtemp = *i1 - *i2;
    res += vtemp * vtemp;
  }

  return std::sqrt(res);
}

template <typename V> typename V::value_type max (const V&, int* p = 0);
template <typename V> typename V::value_type min (const V&, int* p = 0);

template <typename V>
typename V::value_type max (const V& v, int* p)
{
  typename V::value_type res;
  int count = 0;
  for(typename V::const_iterator it = v.begin(); it != v.end(); ++it, ++count)
    if(it == v.begin() || *it > res) {
      res = *it;
      if(p)
	*p = count;
    }
  return res;
}

template <typename V>
typename V::value_type min (const V& v, int* p)
{
  typename V::value_type res;
  int count = 0;
  for(typename V::const_iterator it = v.begin(); it != v.end(); ++it, ++count)
    if(it == v.begin() || *it < res) {
      res = *it;
      if(p)
	*p = count;
    }
  return res;
}

template <typename V>
typename V::value_type sum (const V& v)
{
  typename V::value_type res = 0;
  for(typename V::const_iterator it = v.begin(); it != v.end(); ++it)
    res += *it;
  return res;
}

template <typename V>
typename V::value_type product (const V& v)
{
  typename V::value_type res = 1;

  for(typename V::const_iterator it = v.begin(); it != v.end(); ++it)
    res *= *it;

  return res;
}

#endif
