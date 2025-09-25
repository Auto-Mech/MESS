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

#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include "error.hh"

/**************************************************************************
 *********************** ARRAY WITH DEFAULT VALUE *************************
 **************************************************************************/

template<typename T>
class ArrayWithDefault {
  class _D {
    static const T* _default;
    T               _data;

  public:

    static void set_default (const T& t) { delete _default; _default = new T(t); }
    static const T& default_value () { return *_default; }

    _D () : _data(*_default) {}

    operator       T& ()       { return _data; }
    operator const T& () const { return _data; }
  
    T*       operator& ()       { return &_data; }
    const T* operator& () const { return &_data; }
  };

  _D*       _begin;
  _D*       _end;
  int64_t       _size;
  const T*  _default;

  void _delete ();
  void _create (int64_t);
  void _isinit () const;

public:

  class const_iterator;
  class iterator {
    _D* _value;

  public:

    iterator (_D* v) : _value(v) {}

    T* operator-> () { return &(*_value); }
    T& operator*  () { return *_value; }

    iterator operator++ ()    { return ++_value; }
    iterator operator++ (int) { return _value++; }

    bool operator != (iterator       it) { return _value != it._value; }
    bool operator != (const_iterator it) { return _value != it._value; }

    iterator operator+  (int64_t i) { return _value +  i; }
    iterator operator-  (int64_t i) { return _value -  i; }

    iterator& operator+= (int64_t i) { _value += i; return *this; }
    iterator& operator-= (int64_t i) { _value -= i; return *this; }

    friend class const_iterator;
  };

  class const_iterator {
    const _D* _value;

  public:

    const_iterator (const _D* v) : _value(v) {}
    const_iterator (iterator v) : _value(v._value) {}
    const_iterator (const const_iterator& v) : _value(v._value) {}

    const T* operator-> () { return &(*_value); }
    const T& operator*  () { return *_value; }

    const_iterator operator++ ()    { return ++_value; }
    const_iterator operator++ (int) { return _value++; }

    bool operator != (const_iterator it) { return _value != it._value; }
    bool operator != (iterator       it) { return _value != it._value; }

    const_iterator operator+  (int64_t i) { return _value +  i; }
    const_iterator operator-  (int64_t i) { return _value -  i; }

    const_iterator& operator+= (int64_t i) { _value += i; return *this; }
    const_iterator& operator-= (int64_t i) { _value -= i; return *this; }

    friend class iterator;
  };

  typedef T value_type;

  ArrayWithDefault (int64_t =0);

  ArrayWithDefault (int64_t s, const T& t) : _begin(0), _end(0), _size(0), _default(new T(t)) { _create(s); }

  ~ArrayWithDefault () { _delete(); delete _default; }

  ArrayWithDefault& operator= (const ArrayWithDefault& a);

  ArrayWithDefault (const ArrayWithDefault& a) : _begin(0), _end(0), _size(0), _default(new T(*a._default)) { *this = a; }

  int64_t  size () const { return _size; }

  void resize (int64_t =0);
  void resize (int64_t s, const T& t) { delete _default; _default = new T(t); resize(s); }

  const_iterator begin () const { return _begin; }
  iterator       begin ()       { return _begin; }

  const_iterator   end () const { return _end; }
  iterator         end ()       { return _end; }

  const T& operator[] (int64_t i) const { return *(_begin + i); }
  T&       operator[] (int64_t i)       { return *(_begin + i); }

  operator const T* () const { _isinit(); return &(*_begin); }
  operator       T* ()       { _isinit(); return &(*_begin); }

  const T& front () const { _isinit(); return *_begin; }
  T&       front ()       { _isinit(); return *_begin; }

  const T& back () const { _isinit(); return *(_end - 1); }
  T&       back ()       { _isinit(); return *(_end - 1); }

  // block opperations
  template <typename V> void operator=  (const V&);
  template <typename V> void operator+= (const V&);
  template <typename V> void operator-= (const V&);
  template <typename V> void operator*= (const V&);
  template <typename V> void operator/= (const V&);

  template <typename V> void add (const T&, const V&);
  template <typename F, typename V> void apply (const F&, const V&);

  void operator=  (const T*);
  void operator+= (const T*);
  void operator-= (const T*);
  void operator*= (const T*);
  void operator/= (const T*);

  void operator=  (const T&);
  void operator+= (const T&); 
  void operator-= (const T&); 
  void operator*= (const T&); 
  void operator/= (const T&);

  // apply function-like operation to array
  template<typename F> void apply (const F&);
}; 

template<typename T>
const T* ArrayWithDefault<T>::_D::_default = new T;

template<typename T>
void ArrayWithDefault<T>::_isinit () const
{
  const char funame [] = "ArrayWithDefault<T>::_isinit: ";
 
  if(!size()) {
    std::cerr << funame << "not initialized\n";
    throw std::exception();
  }
}

template<typename T>
void ArrayWithDefault<T>::_delete () {
  if(!size())
    return;

  delete[] _begin;
  _size  = 0;
  _begin = 0;
  _end   = 0;
}

template<typename T>
void ArrayWithDefault<T>::_create (int64_t s) {
  const char funame [] = "ArrayWithDefault<T>::_create: ";

  if(size()) {
    std::cerr << funame << "already exists\n";
    throw std::exception();
  }

  if(s < 0) {
    std::cerr << funame << "negative dimension\n";
    throw std::exception();
  }

  if(!s)
    return;

  _D::set_default(*_default);

  _size = s;
  _begin = new _D[s];
  _end = _begin + s;
}


template<typename T>
ArrayWithDefault<T>::ArrayWithDefault (int64_t s) : _begin(0), _end(0), _size(0), _default(new T(_D::default_value())) {
  _create(s);
}

template<typename T>
ArrayWithDefault<T>& ArrayWithDefault<T>::operator= (const ArrayWithDefault& a) {
  const char funame [] = "ArrayWithDefault<T>::operator=: ";

  if(!a.size()) {
    _delete();
    return *this;
  }

  if(size() != a.size()) {
    _delete();
    _create(a.size());
  }

  _D* rit = _begin;
  for(const _D* it = a._begin; it != a._end; ++it, ++rit)
    *rit = *it;

  return  *this;
}

template<typename T>
void ArrayWithDefault<T>::resize (int64_t s) {
  const char funame [] = "ArrayWithDefault<T>::resize: ";

  if(!s) {
    _delete();
    return;
  }

  if(!size()) {
    _create(s);
    return;
  }

  if(s < 0) {
    std::cerr << funame << "negative dimension\n";
    throw std::exception();
  }

  _D::set_default(*_default);

  _D* new_begin = new _D[s];
  _D* new_end   = new_begin + s;

  if(s < size()) {
    const _D* it = _begin;
    for(_D* nit = new_begin; nit != new_end; ++nit, ++it)
      *nit = *it;
  }
  else {
    _D* nit = new_begin;
    for(const _D* it = _begin; it != _end; ++nit, ++it)
      *nit = *it;
  }

  delete[] _begin;

  _size  = s;
  _begin = new_begin;
  _end   = new_end;  
}

template <typename T>
template <typename V>
void ArrayWithDefault<T>::operator= (const V& v)
{
  const char funame [] = "ArrayWithDefault<T>::operator=: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
}

template <typename T>
template <typename V>
void ArrayWithDefault<T>::operator+= (const V& v)
{
  const char funame [] = "ArrayWithDefault<T>::operator+=: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it += *vit;
}

template <typename T>
template <typename V>
void ArrayWithDefault<T>::operator-= (const V& v)
{
  const char funame [] = "ArrayWithDefault<T>::operator-=: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it -= *vit;
}

template <typename T>
template <typename V>
void ArrayWithDefault<T>::operator*= (const V& v)
{
  const char funame [] = "ArrayWithDefault<T>::operator*=: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it *= *vit;
}

template <typename T>
template <typename V>
void ArrayWithDefault<T>::operator/= (const V& v)
{
  const char funame [] = "ArrayWithDefault<T>::operator/=: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it /= *vit;
}

template <typename T>
template <typename V>
void ArrayWithDefault<T>::add (const T& f, const V& v)
{
  const char funame [] = "ArrayWithDefault<T>::add: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it += *vit * f;
}

template <typename T>
template <typename F, typename V>
void ArrayWithDefault<T>::apply (const F& f, const V& v) 
{
  const char funame [] = "ArrayWithDefault<T>::apply: ";

  if (size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw std::exception();
  }

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = f(*vit);
}

template <typename T>
void ArrayWithDefault<T>::operator= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
}

template <typename T>
void ArrayWithDefault<T>::operator+= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it += *vit;
}

template <typename T>
void ArrayWithDefault<T>::operator-= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it -= *vit;
}

template <typename T>
void ArrayWithDefault<T>::operator*= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it *= *vit;
}

template <typename T>
void ArrayWithDefault<T>::operator/= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it /= *vit;
}

template <typename T>
void ArrayWithDefault<T>::operator= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it = t;
}

template <typename T>
void ArrayWithDefault<T>::operator+= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it += t;
}

template <typename T>
void ArrayWithDefault<T>::operator-= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it -= t;
}

template <typename T>
void ArrayWithDefault<T>::operator*= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it *= t;
}

template <typename T>
void ArrayWithDefault<T>::operator/= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it /= t;
}

template <typename T>
template <typename F>
void ArrayWithDefault<T>::apply (const F& f) 
{
  for(iterator it = begin(); it != end(); ++it)
    *it = f(*it);
}

/**************************************************************************
 ******************************* ARRAY ************************************
 **************************************************************************/

template <typename T>
//
class Array {
  //
  int64_t _capacity;
  int64_t _size;
  T*  _begin;
  T*  _end;

  void _init () {  _capacity = 0; _size = 0; _begin = 0; _end = 0; }

  int64_t _compare (const Array&) const;
  
public:
  //
  typedef       T*       iterator;
  typedef const T* const_iterator;
  typedef       T      value_type;

  Array& operator=  (const T&);

  void resize  (int64_t);
  void reserve (int64_t);
  void compact ();

  Array (const Array&);
  Array& operator= (const Array&);

  Array		 ()		   { _init();		 }
  explicit Array (int64_t	s) { _init(); resize(s); }
  explicit Array (int		s) { _init(); resize(s); }
  explicit Array (unsigned	s) { _init(); resize(s); }
  explicit Array (long unsigned s) { _init(); resize(s); }

  Array (int64_t  s, const T& t) { _init(); resize(s); (*this) = t; }

  template <typename V>
  explicit Array (const V&);

  operator std::vector<T> () const;

  ~Array () { if(_begin) delete[] _begin; }

  T*       begin ()       { return _begin; }
  T*       end   ()       { return _end; }
  const T* begin () const { return _begin; }
  const T* end   () const { return _end; }

  operator       T* ()       { return _begin; }
  operator const T* () const { return _begin; }

  bool operator> (const Array& v) const { if(_compare(v) > 0) return true; return false; }
  bool operator< (const Array& v) const { if(_compare(v) < 0) return true; return false; }

  bool operator>= (const Array& v) const { if(_compare(v) >= 0) return true; return false; }
  bool operator<= (const Array& v) const { if(_compare(v) <= 0) return true; return false; }

  bool operator== (const Array& v) const { if(!_compare(v)) return true; return false; }
  bool operator!= (const Array& v) const { if(_compare(v))  return true; return false; }
  
  T&       front ();
  T&       back  ();
  const T& front () const;
  const T& back  () const;

  int64_t size     () const { return _size; }
  int64_t capacity () const { return _capacity; }

  T&       operator[] (int64_t);
  const T& operator[] (int64_t) const;

  // block operations
  //
  Array& operator- ();

  template <typename V> Array& operator=  (const V&);
  template <typename V> Array& operator+= (const V&) ;
  template <typename V> Array& operator-= (const V&) ;

  template <typename V> void add (const T&, const V&) ;

  Array& operator=  (const T*);
  Array& operator+= (const T*);
  Array& operator-= (const T*);

  Array& operator+= (const T&); 
  Array& operator-= (const T&); 
  Array& operator*= (const T&); 
  Array& operator/= (const T&);

  // apply function-like operation to array
  template <typename F>
  void apply (const F&);

};// class Array

template <typename T>
//
Array<T>::operator std::vector<T>() const
{
  std::vector<T> res(size());

  for(int64_t i = 0; i < size(); ++i)
    //
    res[i] = (*this)[i];

  return res;
}

template <typename T>
//
int64_t Array<T>::_compare (const Array& v) const
{
  if(size() > v.size())
    return 1;
  else if(size() < v.size())
    return -1;

  const_iterator vit = v.begin();
  for(const_iterator it = begin(); it != end(); ++it, ++vit)
    if(*it > *vit)
      return 1;
    else if(*it < *vit)
      return -1;

  return 0;
}
  
template <typename T>
Array<T>::Array (const Array& v)
  : _size(v.size()), _capacity(v.size())
{
  if(!_size) {
    _end = _begin = 0;
    return;
  }
  
  _begin = new T[_size];
  _end = _begin + _size;

  const_iterator vit = v.begin();
  for(T* it = begin(); it != end(); ++it, ++vit)
    *it = *vit; 
}

template <typename T>
Array<T>& Array<T>::operator= (const Array& v)
{
  const char funame [] = "Array<T>::operator=: ";

  _size = v.size();

  if(!v.size()) {
    _end = _begin;
    return *this;
  }


  if(v.size() > _capacity) {
    _capacity = v.size();
    if(_begin)
      delete[] _begin;
    _begin = new T[_size];
  }

  _end  = _begin + _size;

  const T* vit = v.begin();
  for(T* it = begin(); it != end(); ++it, ++vit)
    *it = *vit;

  return *this;
}

template <typename T>
template <typename V>
Array<T>::Array (const V& v)
  : _size(v.size()), _capacity(v.size())
{
  if(!_size) {
    //
    _end = _begin = 0;

    return;
  }
  
  _begin = new T[_size];
  _end = _begin + _size;

  typename V::const_iterator vit = v.begin();
  for(T* it = begin(); it != end(); ++it, ++vit)
    *it = *vit; 
}

template <typename T>
void Array<T>::compact ()
{
    if(_capacity == _size)
	return;

    if(!_size) {
	delete[] _begin;
	_end = _begin = 0;
	_capacity = 0;
	return;
    }

    T* old_begin = _begin;

    _begin = new T[_size];
    _end = _begin + _size;
    _capacity = _size;

    const T* vit = old_begin;
    for(T* it = begin(); it != end(); ++it, ++vit)
	*it = *vit;

    delete[] old_begin;
}

template <typename T>
void Array<T>::resize (int64_t s) 
{
  const char funame [] = "Array<T>::resize: ";

  if(!_begin && !s)
    //
    return;

  if(s <= 0) {
    //
    std::cerr << funame << "out of range: " << s << "\n";

    throw Error::Range();
  }

  if(s <= _capacity) {
    _size = s;
    _end = _begin + _size;
    return;
  }

  T* new_begin = new T[s];

  if(_begin) {
    T* vit = new_begin;
    for(const T* it = begin(); it != end(); ++it, ++vit)
      *vit = *it;
    delete[] _begin;
  } else {
    // Initialize new memory to zero for primitive types
    std::fill(new_begin, new_begin + s, T(0));
  }

  _capacity = _size = s;
  _begin = new_begin;
  _end = _begin + _size;
}

template <typename T>
void Array<T>::reserve (int64_t s) 
{
  const char funame [] = "Array<T>::reserve: ";

  if(s <= _capacity)
      return;

  _capacity  = s;

  T* old_begin = _begin;
  _begin = new T[s];
  _end = _begin + _size;

  if(old_begin) {
    const T* vit = old_begin;
    for(T* it = begin(); it != end(); ++it, ++vit)
      *it = *vit;
    delete[] old_begin;
  }
}

template <typename T>
T& Array<T>::operator[] (int64_t i) 
{
  const char funame [] = "Array<T>::operator[]: ";

#ifdef DEBUG
  if (i >= _size || i < 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
#endif

  return _begin[i];
}

template <typename T>
const T& Array<T>::operator[] (int64_t i) const 
{
  const char funame [] = "Array<T>::operator[]: ";

#ifdef DEBUG
  if (i >= _size || i < 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
#endif

  return _begin[i];
}

template <typename T>
T& Array<T>::front () 
{
  if(!_size) {
    std::cerr << "Array<T>::front(): array is empty\n";
    throw Error::Range();
  }
  return *_begin;
}

template <typename T>
const T& Array<T>::front () const 
{
  if(!_size) {
    std::cerr << "Array<T>::front(): array is empty\n";
    throw Error::Range();
  }
  return *_begin;
}

template <typename T>
T& Array<T>::back () 
{
  if(!_size) {
    std::cerr << "Array<T>::back(): array is empty\n";
    throw Error::Range();
  }
  return _begin[_size-1];
}

template <typename T>
const T& Array<T>::back () const 
{
  if(!_size) {
    std::cerr << "Array<T>::back(): array is empty\n";
    throw Error::Range();
  }
  return _begin[_size-1];
}

template <typename T>
Array<T>& Array<T>::operator- ()
{
  for (T* it = begin(); it != end(); ++it)
    *it *= -*it;
  return *this;
}

template <typename T>
//
template <typename V>
//
Array<T>& Array<T>::operator= (const V& v)
{
  const char funame [] = "Array<T>::operator=: ";

  if (!v.size()) {
    //
    std::cerr << funame << "empty source\n";
    
    throw Error::Init();
  }

  resize(v.size());

  typename V::const_iterator vit = v.begin();
  
  for(T* it = begin(); it != end(); ++it, ++vit)
    //
    *it = *vit;

  return *this;
}

template <typename T>
//
template <typename V>
//
Array<T>& Array<T>::operator+= (const V& v) 
{
  const char funame [] = "Array<T>::operator+=: ";

  if (size() != v.size()) {
    //
    std::cerr << funame << "sizes mismatch: " << size() << " vs. " << v.size() << "\n";
    
    throw Error::Range();
  }

  typename V::const_iterator vit = v.begin();
  
  for(T* it = begin(); it != end(); ++it, ++vit)
    //
    *it += *vit;
  
  return *this;
}

template <typename T>
//
template <typename V>
//
void Array<T>::add (const T& f, const V& v) 
{
  const char funame [] = "Array<T>::add: ";

  if (size() != v.size()) {
    //
    std::cerr << funame << "sizes mismatch: " << size() << " vs. " << v.size() << "\n";
    
    throw Error::Range();
  }

  typename V::const_iterator vit = v.begin();
  
  for(T* it = begin(); it != end(); ++it, ++vit)
    //
    *it += *vit * f;
}

template <typename T>
//
template <typename V>
//
Array<T>& Array<T>::operator-= (const V& v) 
{
  const char funame [] = "Array<T>::operator-=: ";

  if (size() != v.size()) {
    //
    std::cerr << funame << "sizes mismatch: " << size() << " vs. " << v.size() << "\n";
    
    throw Error::Range();
  }

  typename V::const_iterator vit = v.begin();
  
  for(T* it = begin(); it != end(); ++it, ++vit)
    //
    *it -= *vit;

  return *this;
}


template <typename T>
Array<T>& Array<T>::operator= (const T* vit)
{
  for (T* it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const T* vit)
{
  for (T* it = begin(); it != end(); ++it, ++vit)
      *it += *vit;
  
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const T* vit)
{
  for (T* it = begin(); it != end(); ++it, ++vit)
      *it -= *vit;
  
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it = val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it += val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it -= val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator*= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it *= val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator/= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it /= val;
  return *this;
}

template <typename T>
template <typename F>
void Array<T>::apply (const F& f) 
{
  for(iterator it = begin(); it != end(); ++it)
    *it = f(*it);
}

/*****************************************************************************************
 ***************************************** STRIDE POINTER ********************************
 *****************************************************************************************/
template<typename T> class ConstStridePointer;

template<typename T>
class StridePointer {
  T* _pnt;
  int64_t _stride;

public:
  StridePointer () : _pnt(0), _stride(0) {}
  StridePointer (T*, int64_t =1) ;
  int64_t stride () const { return _stride; }

  T& operator*  () const { return *_pnt; }
  T* operator-> () const { return  _pnt; }
  operator T*   () const { return  _pnt; }
  
  operator ConstStridePointer<T> () const { return ConstStridePointer<T>(_pnt, _stride); }

  StridePointer operator++ ()      { _pnt += _stride; return *this; }
  StridePointer operator-- ()      { _pnt -= _stride; return *this; }
  StridePointer operator++ (int)   { StridePointer p = *this; _pnt += _stride; return p; }
  StridePointer operator-- (int)   { StridePointer p = *this; _pnt -= _stride; return p; }
  void          operator+= (int64_t d) { _pnt += d * _stride; }
  void          operator-= (int64_t d) { _pnt -= d * _stride; }

  StridePointer operator+  (int64_t d)         const { return StridePointer(_pnt + d * _stride, _stride); }
  StridePointer operator-  (int64_t d)         const { return StridePointer(_pnt - d * _stride, _stride); }

  int64_t           operator-  (StridePointer) const ;
  bool          operator!= (StridePointer) const ;
};

template <typename T> 
StridePointer<T> operator+ (int64_t d, StridePointer<T> p) 
{ 
  return p + d;
}

template <typename T> 
int64_t distance (StridePointer<T> first, StridePointer<T> last)
{
  return last - first;
}

template <typename T> 
bool StridePointer<T>::operator!= (StridePointer p) const 
{
  const char funame [] = "StridePointer<T>::operator!=: ";

#ifdef DEBUG

  if(_stride != p._stride) {
    std::cerr << funame << "strides mismatch\n";
    throw Error::Range();
  }

  if(((_pnt - p._pnt) % _stride)) {
    std::cerr << funame << "slices mismatch\n";
    throw Error::Range();
  }
 
#endif

  if(_pnt != p._pnt)
    return true;

  return false;
}

template<typename T> 
int64_t StridePointer<T>::operator- (StridePointer p) const 
{
  const char funame [] = "StridePointer<T>::operator-: ";

#ifdef DEBUG

  if(_stride != p._stride) {
    std::cerr << funame << "strides mismatch\n";
    throw Error::Range();
  }

  if(((_pnt - p._pnt) % _stride)) {
    std::cerr << funame << "slices mismatch\n";
    throw Error::Range();
  }
 
#endif

  return (_pnt - p._pnt) / _stride;
}

template<typename T>
StridePointer<T>::StridePointer (T* p, int64_t s)  : _pnt(p), _stride(s) 
{
  const char funame [] = "StridePointer<T>::StridePointer: ";

#ifdef DEBUG

  if(!p) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(s <= 0) {
    std::cerr << funame << "non-positive stride = " << s << std::endl;
    throw Error::Range();
  }

#endif

}

/*************************************************************************************
 ****************************** CONSTANT STRIDE POINTER ******************************
 *************************************************************************************/

template<typename T> 
class ConstStridePointer {
  const T*  _pnt;
  int64_t    _stride;

public:
  ConstStridePointer () : _pnt(0), _stride(0) {}
  ConstStridePointer (const T*, int64_t =1) ;
  ConstStridePointer (StridePointer<T> p) : _pnt(p), _stride(p.stride()) {}
  int64_t stride () const { return _stride; }

  const T& operator*  () const { return *_pnt; }
  const T* operator-> () const { return  _pnt; }
  operator const T*   () const { return  _pnt; }

  ConstStridePointer operator++ ()      { _pnt += _stride; return *this; }
  ConstStridePointer operator-- ()      { _pnt -= _stride; return *this; }
  ConstStridePointer operator++ (int)   { ConstStridePointer p = *this; _pnt += _stride; return p; }
  ConstStridePointer operator-- (int)   { ConstStridePointer p = *this; _pnt -= _stride; return p; }
  void               operator+= (int64_t d) { _pnt += d * _stride; }
  void               operator-= (int64_t d) { _pnt -= d * _stride; }

  ConstStridePointer operator+ (int64_t d) const { return ConstStridePointer(_pnt + d * _stride, _stride); }
  ConstStridePointer operator- (int64_t d) const { return ConstStridePointer(_pnt - d * _stride, _stride); }

  int64_t  operator-  (ConstStridePointer) const ;
  bool operator!= (ConstStridePointer) const ;
};

template <typename T> 
int64_t distance (ConstStridePointer<T> first, ConstStridePointer<T> last)
{
  return last - first;
}

template <typename T> 
bool ConstStridePointer<T>::operator!= (ConstStridePointer p) const 
{
  const char funame [] = "ConstStridePointer<T>::operator!=: ";

#ifdef DEBUG

  if(_stride != p._stride) {
    std::cerr << funame << "strides mismatch\n";
    throw Error::Range();
  }

  if(((_pnt - p._pnt) % _stride)) {
    std::cerr << funame << "slices mismatch\n";
    throw Error::Range();
  }
 
#endif

  if(_pnt != p._pnt)
    return true;

  return false;
}

template <typename T> 
ConstStridePointer<T> operator+ (int64_t d, ConstStridePointer<T> p) 
{ 
  return p + d;
}

template <typename T> 
int64_t ConstStridePointer<T>::operator- (ConstStridePointer p) const 
{
  const char funame [] = "ConstStridePointer::operator-: ";

#ifdef DEBUG

  if(_stride != p._stride) {
    std::cerr << funame << "strides are different\n";
    throw Error::Range();
  }

  if(((_pnt - p._pnt) % _stride)) {
    std::cerr << funame << "different slices\n";
    throw Error::Range();
  }

#endif
 
  return (_pnt - p._pnt) / _stride;

}

template <typename T>
ConstStridePointer<T>::ConstStridePointer (const T* p, int64_t s)  
    : _pnt(p), _stride(s) 
{
  const char funame [] = "ConstStridePointer<T>::ConstStridePointer: ";

#ifdef DEBUG

  if(!p) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(s <= 0) {
    std::cerr << funame << "non-positive stride = " << s << std::endl;
    throw Error::Range();
  }

#endif

}

/**************************************************************************
 **************************** CONSTANT SLICE ******************************
 **************************************************************************/

template <typename T>
class ConstSlice {
  const T* _begin;
  const T*   _end;
  int64_t       _size;
  int64_t     _stride;

  ConstSlice& operator= (const ConstSlice&);

public:
  typedef ConstStridePointer<T> const_iterator;
  typedef                    T      value_type;

  ConstSlice (const T* st, int64_t sz, int64_t sd =1) ;

  const_iterator begin () const { return ConstStridePointer<T>(_begin, _stride); }
  const_iterator end   () const { return ConstStridePointer<T>(_end,   _stride); }

  int64_t   size () const { return   _size; }
  int64_t stride () const { return _stride; }

  const T&  operator[] (int64_t) const ;

};

template <typename T>
ConstSlice<T>::ConstSlice(const T* st, int64_t sz, int64_t sd) 
  : _begin(st), _size(sz), _stride(sd), _end(st + sz * sd)
{
  const char funame [] = "ConstSlice<T>::ConstSlice: ";

#ifdef DEBUG

  if(!st) {
    std::cerr << funame << "no data\n";
    throw Error::Range();
  }

  if(sz < 0) {
    std::cerr << funame << "negative size = " << sz << std::endl;
    throw Error::Range();
  }

  if(sd <= 0) {
    std::cerr << funame << "non-positive stride = " << sd << std::endl;
    throw Error::Range();
  }

#endif

}

template <typename T>
const T& ConstSlice<T>::operator[] (int64_t i) const 
{
  const char funame [] = "ConstSlice<T>::operator[]: ";

#ifdef DEBUG

  if(i < 0 || i >= _size) {
    std::cerr << funame << "out of range: index = " << i << std::endl;
    throw Error::Range();
  }

#endif

  return _begin[i * _stride];
}

/**************************************************************************
 ********************************** SLICE *********************************
 **************************************************************************/

template <typename T>
class Slice 
{
  T*  _begin;
  T*  _end;
  int64_t _size;
  int64_t _stride;

public:
  typedef      StridePointer<T>       iterator;
  typedef ConstStridePointer<T> const_iterator;
  typedef                    T      value_type;

  void operator= (const Slice&) ;
  
  Slice(T* st, int64_t sz, int64_t sd =1) ;

  iterator       begin ()       { return      StridePointer<T>(_begin, _stride); }
  iterator         end ()       { return      StridePointer<T>(_end,   _stride); }
  const_iterator begin () const { return ConstStridePointer<T>(_begin, _stride); }
  const_iterator   end () const { return ConstStridePointer<T>(_end,   _stride); }

  operator ConstSlice<T> () { return ConstSlice<T>(_begin, _size, _stride); }
  
  int64_t size   () const { return _size; }
  int64_t stride () const { return _stride; }

  const T&  operator[] (int64_t) const ;
  T&        operator[] (int64_t)       ;

  // block operations
  template <typename V>  void operator=  (const V&) ;
  template <typename V>  void operator+= (const V&) ;
  template <typename V>  void operator-= (const V&) ;
  template <typename V>  void operator*= (const V&) ;
  template <typename V>  void operator/= (const V&) ;

  template <typename V>  void add (const T&, const V&) ;

  void operator=  (const T*);
  void operator+= (const T*);
  void operator-= (const T*);
  void operator*= (const T*);
  void operator/= (const T*);

  void operator-  ();
  void operator=  (const T&);
  void operator+= (const T&);
  void operator-= (const T&);
  void operator*= (const T&);
  void operator/= (const T&);

  // apply function-like operation to array
  template <typename F>
  void apply (const F&);
};

template <typename T>
Slice<T>::Slice(T* st, int64_t sz, int64_t sd) 
  : _begin(st), _size(sz), _stride(sd), _end(st + sz * sd)
{
  const char funame [] = "Slice<T>::Slice: ";

#ifdef DEBUG
  if(!st) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(sz < 0) {
    std::cerr << funame << "negative size = " << sz << std::endl;
    throw Error::Range();
  }

  if(sd <= 0) {
    std::cerr << funame << "negative stride = " << sd << std::endl;
    throw Error::Range();
  }
#endif

}

template <typename T>
void Slice<T>::operator= (const Slice& v) 
{
  const char funame [] = "Slice<T>::operator=: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
}

template <typename T>
const T& Slice<T>::operator[] (int64_t i) const 
{
  const char funame [] = "Slice<T>::operator[]: ";

#ifdef DEBUG

  if(i < 0 || i >= size()) {
    std::cerr << funame << "out of range: index = " << i << std::endl;
    throw Error::Range();
  }

#endif

  return _begin[i * _stride];
}

template <typename T>
T& Slice<T>::operator[] (int64_t i) 
{
  const char funame [] = "Slice<T>::operator[]: ";

#ifdef DEBUG

  if(i < 0 || i >= size()) {
    std::cerr << funame << "out of range: index = " << i << std::endl;
    throw Error::Range();
  }

#endif

  return _begin[i * _stride];
}

template <typename T>
template <typename V>
void Slice<T>::operator= (const V& v) 
{
  const char funame [] = "Slice<T>::operator=: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
}

template <typename T>
template <typename V>
void Slice<T>::operator+= (const V& v) 
{
  const char funame [] = "Slice<T>::operator+=: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it += *vit;
}

template <typename T>
template <typename V>
void Slice<T>::add (const T& f, const V& v) 
{
  const char funame [] = "Slice<T>::add: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it +=  *vit * f;
}

template <typename T>
template <typename V>
void Slice<T>::operator-= (const V& v) 
{
  const char funame [] = "Slice<T>::operator-=: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it -= *vit;
}

template <typename T>
template <typename V>
void Slice<T>::operator*= (const V& v) 
{
  const char funame [] = "Slice<T>::operator*=: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it *= *vit;
}

template <typename T>
template <typename V>
void Slice<T>::operator/= (const V& v) 
{
  const char funame [] = "Slice<T>::operator/=: ";
  
#ifdef DEBUG

  if(size() != v.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

#endif  

  typename V::const_iterator vit = v.begin();
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it /= *vit;
}

template <typename T>
void Slice<T>::operator= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it = *vit;
}

template <typename T>
void Slice<T>::operator+= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it += *vit;
}

template <typename T>
void Slice<T>::operator-= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it -= *vit;
}

template <typename T>
void Slice<T>::operator*= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it *= *vit;
}

template <typename T>
void Slice<T>::operator/= (const T* vit)
{
  for(iterator it = begin(); it != end(); ++it, ++vit)
    *it /= *vit;
}

template <typename T>
void Slice<T>::operator- ()
{
  for(iterator it = begin(); it != end(); ++it)
    *it = -*it;
}

template <typename T>
void Slice<T>::operator= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it = t;
}

template <typename T>
void Slice<T>::operator+= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it += t;
}

template <typename T>
void Slice<T>::operator-= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it -= t;
}

template <typename T>
void Slice<T>::operator*= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it *= t;
}

template <typename T>
void Slice<T>::operator/= (const T& t)
{
  for(iterator it = begin(); it != end(); ++it)
    *it /= t;
}

template <typename T>
template <typename F>
void Slice<T>::apply (const F& f) 
{
  for(iterator it = begin(); it != end(); ++it)
    *it = f(*it);
}
 
/**************************************************************************
 *************************** ARRAY BY REFERENCE  **************************
 **************************************************************************/

template <typename T>
class RefArr {

  Array<T>* _data;
  int*      _count;  // number of references

  void delete_ref ();
  void create_ref (const RefArr&);

  void _assert () const;
  
public:

  bool isinit () const { return _data; }

  RefArr ()				  : _data(0), _count(0) {}
  explicit RefArr (int64_t	 s)	  : _data(0), _count(0) { if(!s) return; _data = new Array<T>(s); _count = new int(1); }
  explicit RefArr (int		 s)	  : _data(0), _count(0) { if(!s) return; _data = new Array<T>(s); _count = new int(1); }
  explicit RefArr (unsigned	 s)	  : _data(0), _count(0) { if(!s) return; _data = new Array<T>(s); _count = new int(1); }
  explicit RefArr (long unsigned s)	  : _data(0), _count(0) { if(!s) return; _data = new Array<T>(s); _count = new int(1); }

  RefArr (int64_t s, const T& t)	  : _data(0), _count(0) { if(!s)	return; _data = new Array<T>(s, t); _count = new int(1); }

  template<typename V> RefArr(const V& v) : _data(0), _count(0) { if(!v.size()) return; _data = new Array<T>(v);    _count = new int(1); }

  RefArr (const RefArr& a) { create_ref(a); }
  ~RefArr () { delete_ref(); }

  RefArr& operator= (const RefArr& a) { if(_data != a._data) { delete_ref(); create_ref(a); } return *this; }
  RefArr copy() const;// make a new copy
    
  void resize  (int64_t s);
  void reserve (int64_t s);
  void compact ();

  int64_t size      () const;
  int64_t capacity  () const;
  int64_t ref_count () const;

  T*       begin ();
  const T* begin () const;
  T*       end   ();
  const T* end   () const;

  T&       front ();
  const T& front () const;
  T&       back  ();
  const T& back  () const;

  // conversions
  operator       T* ();
  operator const T* () const;

  // indexing
  const T& operator[] (int64_t i) const ;
  T&       operator[] (int64_t i)       ;

  // comparisons
  bool operator== (const RefArr& a) const { return _data == a._data; }
  bool operator!= (const RefArr& a) const { return _data != a._data; }

  // arithmetic operations
  RefArr operator+ (const RefArr&) const ;
  RefArr operator- (const RefArr&) const ;

  RefArr operator+= (const RefArr& a) ;
  RefArr operator-= (const RefArr& a) ;

  template<typename V> RefArr operator= (const V&);

  RefArr operator= (const T*);

  RefArr operator-  ()           ;
  RefArr operator=  (const T& t) ;
  RefArr operator*= (const T& t) ;
  RefArr operator/= (const T& t) ;

};// class RefArr

template<typename T>
inline void RefArr<T>::_assert () const
{
  const char funame [] = "RefArr<T>::_assert: ";

  if(_data)
    //
    return;

  std::cerr << funame << "not initialized\n";

  throw Error::Init();
}
template <typename T>
inline void RefArr<T>::create_ref (const RefArr& a)
{
    _data = a._data;
    _count =  a._count;
    if(_count)
      ++(*_count);
}

template <typename T>
inline void RefArr<T>::delete_ref ()
{
  if(!_count)
    return;

  if(*_count > 1)
      --(*_count);
  else {
      delete _data;
      delete _count;
  }
}

template <typename T>
RefArr<T> RefArr<T>::copy() const
{
  const char funame [] = "RefArr<T>::copy: ";

  if(!isinit())
    //
    return RefArr<T>();

  RefArr<T> res(size());
  const T* p1 = *this;
  for(T* p = res.begin(); p != res.end(); ++p, ++p1)
    *p = *p1;

  return res;
}

template <typename T>
inline int64_t  RefArr<T>::size ()  const 
{ 
  const char funame [] = "RefArr<T>::size: ";

  if(_data)
    return _data->size();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline void RefArr<T>::resize  (int64_t s) 
{ 
  if(!_data) {
    _data = new Array<T>(s);
    _count= new int(1);
  }
  else
    _data->resize(s);  
}

template <typename T>
inline void RefArr<T>::reserve (int64_t s) 
{ 
  if(!_data) {
    _data = new Array<T>(0);
    _count= new int(1);
  }
  _data->reserve(s); 
}

template <typename T>
inline void  RefArr<T>::compact () 
{ 
  const char funame [] = "RefArr<T>::compact: ";

  if(_data)
    _data->compact();
  
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline int64_t  RefArr<T>::capacity () const 
{ 
  const char funame [] = "RefArr<T>::capacity: ";

  if(_data)
    return _data->capacity(); 

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline int64_t  RefArr<T>::ref_count () const 
{ 
  const char funame [] = "RefArr<T>::ref_count: ";

  if(_count) 
    return *_count;

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline T&  RefArr<T>::back () 
{
  const char funame [] = "RefArr<T>::back: ";

  if(_data)
    return _data->back();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline const T&  RefArr<T>::back () const 
{
  const char funame [] = "RefArr<T>::back: ";

  if(_data)
    return _data->back();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline T&  RefArr<T>::front () 
{
  const char funame [] = "RefArr<T>::front: ";

  if(_data)
    return _data->front();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline const T&  RefArr<T>::front () const 
{
  const char funame [] = "RefArr<T>::front: ";

  if(_data)
    return _data->front();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline T*  RefArr<T>::begin () 
{
  const char funame [] = "RefArr<T>::begin: ";

  if(_data)
    return _data->begin();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline T*  RefArr<T>::end () 
{
  const char funame [] = "RefArr<T>::end: ";

  if(_data)
    return _data->end();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline const T*  RefArr<T>::begin () const 
{
  const char funame [] = "RefArr<T>::begin: ";

  if(_data)
    return _data->begin();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline const T*  RefArr<T>::end () const 
{ 
  const char funame [] = "RefArr<T>::end: ";

  if(_data)
    return _data->end();

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

// conversions
template <typename T>
inline RefArr<T>::operator T* () 
{ 
  if(_data)
    return *_data; 
  return 0;
}

template <typename T>
inline RefArr<T>::operator const T* () const 
{ 
  if(_data) 
    return *_data;
  return 0;
}

  // indexing
template <typename T>
const T&  RefArr<T>::operator[] (int64_t i) const 
{
  const char funame [] = "RefArr<T>::operator[]: ";

  if(_data)
    return (*_data)[i];

  std::cerr << funame << "not initialized\n";
  throw Error::Init();  
}

template <typename T>
T& RefArr<T>::operator[] (int64_t i) 
{
  const char funame [] = "RefArr<T>::operator[]: ";

  if(_data)
    return (*_data)[i];

  std::cerr << funame << "not initialized\n";
  throw Error::Init();  
}

template <typename T>
RefArr<T>  RefArr<T>::operator+= (const RefArr& a) 
{
  const char funame [] = "RefArr<T>::operator+=: ";

  if(_data && a._data) {
    *_data += *a._data;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
RefArr<T>  RefArr<T>::operator-= (const RefArr& a) 
{
  const char funame [] = "RefArr<T>::operator-=: ";

  if(_data && a._data) {
    *_data -= *a._data;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline RefArr<T>  RefArr<T>::operator- () 
{
  const char funame [] = "RefArr<T>::operator-: ";

  if(_data) {
    -(*_data);
    return *this;
  }

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline RefArr<T>  RefArr<T>::operator=  (const T& t)  
{
  const char funame [] = "RefArr<T>::operator=: ";

  if(_data) {
    *_data =  t;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline RefArr<T>  RefArr<T>::operator=  (const T* t)  
{
  const char funame [] = "RefArr<T>::operator=: ";

  if(_data) {
    *_data =  t;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
template <typename V>
inline RefArr<T>  RefArr<T>::operator=  (const V& v)  
{
  const char funame [] = "RefArr<T>::operator=: ";

  if(_data) {
    *_data =  v;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline RefArr<T>  RefArr<T>::operator*= (const T& t)  
{
  const char funame [] = "RefArr<T>::operator*=: ";

  if(_data) {
    *_data *= t;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
inline RefArr<T>  RefArr<T>::operator/= (const T& t) 
{ 
  const char funame [] = "RefArr<T>::operator/=: ";

  if(_data) {
    *_data /= t;
    return *this;
  }
    
  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

template <typename T>
RefArr<T> RefArr<T>::operator+ (const RefArr& a) const 
{
  const char funame [] = "RefArr<T>::operator+: ";
  
  if(!isinit() || !a.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(a.size() != size()) {
    std::cerr << funame << "size mismatch\n";
    throw Error::Range();
  }

  RefArr<T> res(size());
  const T* p1 = *this;
  const T* p2 = a;
  for(T* p = res.begin(); p != res.end(); ++p, ++p1, ++p2)
    *p = *p1 + *p2;

  return res;
}

template <typename T>
RefArr<T> RefArr<T>::operator- (const RefArr& a) const 
{
  const char funame [] = "RefArr<T>::operator-: ";

  if(!isinit() || !a.isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(a.size() != size()) {
    std::cerr << funame << "size mismatch\n";
    throw Error::Range();
  }

  RefArr<T> res(size());
  const T* p1 = *this;
  const T* p2 = a;
  for(T* p = res.begin(); p != res.end(); ++p, ++p1, ++p2)
    *p = *p1 - *p2;

  return res;
}

template<class C>
class Array_2 : public Array<C>
{// class Array_2

  int64_t dd_1;
  int64_t dd_2;
  
public:

  Array_2 (int64_t d1, int64_t d2) : Array<C>(d1 * d2), dd_1(d1), dd_2(d2) {}
  
  void resize (int64_t d1, int64_t d2) { Array<C>::resize(d1 * d2); dd_1 = d1; dd_2 = d2; }

  const C& operator () (int64_t i1, int64_t i2) const
  {
    return (*this) [i1 + dd_1 * i2];
  }
  
  C& operator () (int64_t i1, int64_t i2)
  {
    return (*this) [i1 + dd_1 * i2];
  }

  int64_t dim_1 () const { return dd_1; }
  int64_t dim_2 () const { return dd_2; }
};

template <class C>
class Array_3 : public Array<C>
{// class Array_3

  int64_t dd_1;
  int64_t dd_2;
  int64_t dd_3;

public:

  Array_3 (int64_t d1, int64_t d2, int64_t d3)
    : Array<C>(d1 * d2 * d3), dd_1(d1), dd_2(d2), dd_3(d3)
  {}

  void resize (int64_t d1, int64_t d2, int64_t d3) { Array<C>::resize(d1 * d2 * d3); dd_1 = d1; dd_2 = d2; dd_3 = d3; }

  const C& operator () (int64_t i1, int64_t i2, int64_t i3) const
  {
    return (*this)[i1 + dd_1 * i2 + dd_1 * dd_2 * i3];
  }

  C& operator () (int64_t i1, int64_t i2, int64_t i3)
  {
    return (*this)[i1 + dd_1 * i2 + dd_1 * dd_2 * i3];
  }

  int64_t dim_1 () const { return dd_1; }
  int64_t dim_2 () const { return dd_2; }
  int64_t dim_3 () const { return dd_3; }

};

/**************************************************************************
 ******************************* Matrix ************************************
 **************************************************************************/

template <typename T>
class Vector : public ArrayWithDefault<T> {
  //
public:
  //
  Vector () {}

  explicit Vector (int64_t i) : ArrayWithDefault<T>(i) {}

  template<typename V> explicit Vector (const V& v) : ArrayWithDefault<T>(v) {}
};

template <typename T>
class Matrix : public ArrayWithDefault<T> {
  //
  int64_t _size1, _size2;

  void _resize (int64_t s1, int64_t s2);

  void _assert () const;

  void _assert (int64_t i1, int64_t i2) const;

  void _assert (const Matrix&) const;
  
public:

  bool isinit () const { return _size1; }
  
  void resize (int64_t s1, int64_t s2)             { _resize(s1, s2); ArrayWithDefault<T>::resize(s1 * s2);    }
  void resize (int64_t s1, int64_t s2, const T& t) { _resize(s1, s2); ArrayWithDefault<T>::resize(s1 * s2, t); }

  void resize (int64_t s) { resize(s, s); }
  
  Matrix () : _size1(0), _size2(0) {} 

  Matrix (int64_t s1, int64_t s2)             { resize(s1, s2);    }
  Matrix (int64_t s1, int64_t s2, const T& t) { resize(s1, s2, t); }

  Matrix (int64_t s)             { resize(s, s);    }
  Matrix (int64_t s, const T& t) { resize(s, s, t); }

  int64_t size1 () const { return _size1; }
  int64_t size2 () const { return _size2; }

  int64_t size () const;

  operator       T* ()       { return ArrayWithDefault<T>::operator       T*(); }
  operator const T* () const { return ArrayWithDefault<T>::operator const T*(); }

  ConstSlice<T> column (int64_t i) const { _assert(0, i); return ConstSlice<T>((const T*)*this + i * size1(), size1()); }
  Slice<T>      column (int64_t i)       { _assert(0, i); return      Slice<T>((T*)*this + i * size1(), size1()); }

  ConstSlice<T> row (int64_t i) const { _assert(i, 0); return ConstSlice<T>((const T*)*this + i, size2(), size1()); }
  Slice<T>      row (int64_t i)       { _assert(i, 0); return      Slice<T>((T*)*this + i, size2(), size1()); }

  ConstSlice<T> diagonal (int64_t =0) const;
  Slice<T>      diagonal (int64_t =0);

  const T& operator () (int64_t i1, int64_t i2) const { _assert(i1, i2); return (*this)[i1 + _size1 * i2]; }
  T&       operator () (int64_t i1, int64_t i2)       { _assert(i1, i2); return (*this)[i1 + _size1 * i2]; }

  // arithmetic operations with matrix
  //
  Matrix& operator+= (const Matrix& m) { _assert(m); ArrayWithDefault<T>::operator+=(m); return *this; }
  Matrix& operator-= (const Matrix& m) { _assert(m); ArrayWithDefault<T>::operator-=(m); return *this; }

  Matrix operator+  (const Matrix& m) const { Matrix res(*this); return res += m; }
  Matrix operator-  (const Matrix& m) const { Matrix res(*this); return res -= m; }

  // arithmetic operations with double
  //
  Matrix& operator-  ()           { _assert();  ArrayWithDefault<T>::operator-  ();  return *this; }
  Matrix& operator=  (const T& d) { _assert();  ArrayWithDefault<T>::operator=  (d); return *this; }
  Matrix& operator*= (const T& d) { _assert();  ArrayWithDefault<T>::operator*= (d); return *this; }
  Matrix& operator/= (const T& d) { _assert();  ArrayWithDefault<T>::operator/= (d); return *this; }
};

template <typename T>
void Matrix<T>::_assert () const
{
  const char funame [] = " Matrix<T>: ";

  if(isinit())
    //
    return;
  
  std::cerr << funame << "not initialized\n";
  
  throw Error::Init();
}

template <typename T>
void Matrix<T>::_assert (int64_t i1, int64_t i2) const
{
  const char funame [] = " Matrix<T>::_assert: ";

  _assert();

  if(i1 >= 0 && i1 < _size1 && i2 >= 0 && i2 < _size2)
    //
    return;
  
  std::cerr << funame << "indices out of range: " << i1 << ", " << i2 << "\n";
    
  throw Error::Range();
}

template <typename T>
inline void Matrix<T>::_assert (const Matrix& m) const 
{
  const char funame [] = "Matrix<T>::_assert: ";

  _assert();

  m._assert();
    
  if(m.size1() != size1() || m.size2() != size2()) {
    //
    std::cerr << funame << "dimensions mismatch: " << size1() << ", " << size2() << ", " << m.size1() << ", " << m.size2() << "\n";
      
    throw Error::Range();
  }
}

template <typename T>
void Matrix<T>::_resize (int64_t s1, int64_t s2)
{
  const char funame [] = " Matrix<T>::_resize: ";

  if(s1 <= 0 || s2 <= 0) {
    //
    std::cerr << funame << "out of range: " << s1 << ", " << s2 << "\n";
    
    throw Error::Range();
  }
  _size1 = s1;
  _size2 = s2;
}

template<typename T>
int64_t Matrix<T>::size () const
{
  const char funame [] = " Matrix<T>::size: ";

  if(size1() != size2()) {
    //
    std::cerr << funame << "sizes mismatch: " << size1() << ", " << size2() << "\n";
    
    throw Error::Logic();
  }

  return size1();
}  

template <typename T>
Slice<T> Matrix<T>::diagonal (int64_t i)
{
  const char funame [] = "Matrix<T>::diagonal: ";
  
  int64_t sz, itemp;

  if(i >= size2() || i <= -size1()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";

    throw Error::Range();
  }

  // upper diagonal
  //
  if(i > 0) {
    //
    itemp = size2() - i;
    
    sz = itemp < size1() ? itemp : size1();
    
    return Slice<double>((T*)*this + i * size1(), sz, size1() + 1);
  }
  // lower diagonal
  //
  else {
    //
    itemp = size1() + i;
    
    sz = itemp < size2() ? itemp : size2();
    
    return Slice<double>((T*)*this - i, sz, size1() + 1);
  }
}

template <typename T>
ConstSlice<T> Matrix<T>::diagonal (int64_t i) const
{
  const char funame [] = "Matrix<T>::diagonal: ";
  
  int64_t sz, itemp;

  if(i >= size2() || i <= -size1()) {
    //
    std::cerr << funame << "out of range: " << i << "\n";

    throw Error::Range();
  }

  // upper diagonal
  //
  if(i > 0) {
    //
    itemp = size2() - i;
    
    sz = itemp < size1() ? itemp : size1();
    
    return ConstSlice<double>((const T*)*this + i * size1(), sz, size1() + 1);
  }
  // lower diagonal
  //
  else {
    //
    itemp = size1() + i;
    
    sz = itemp < size2() ? itemp : size2();
    
    return ConstSlice<double>((const T*)*this - i, sz, size1() + 1);
  }
}

#endif
