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

#ifndef SHARED_HH
#define SHARED_HH

#include "error.hh"
#include <iostream>

/************************************************************************************
 ********************************* Shared Pointer ***********************************
 ****************** creates a reference for a new dynamical object  *****************
 ***********************************************************************************/
template<typename T>
class ConstSharedPointer;

template <typename T>
class SharedPointer 
{
private:
  T*   _pnt;
  int* _count;

  void _delete_ref ();
  void _create_ref (const SharedPointer&);

public:
  explicit SharedPointer (T* =0);
  SharedPointer (const SharedPointer& s) { _create_ref(s); }
  SharedPointer& operator= (const SharedPointer& s) { _delete_ref(); _create_ref(s); return *this; }
  ~SharedPointer () { _delete_ref(); }

  void init (T*) ;

  T* operator-> () const ;
  T& operator*  () const ;
  operator T*   () const { return _pnt; }

  operator bool  () const { if(_pnt) return  true; return false; }
  bool operator! () const { if(_pnt) return false; return true;  }

  bool operator== (SharedPointer<T> p) const { return _pnt == p._pnt; }
  bool operator!= (SharedPointer<T> p) const { return _pnt != p._pnt; }
  bool operator<  (SharedPointer<T> p) const { return _pnt <  p._pnt; }
  bool operator>  (SharedPointer<T> p) const { return _pnt >  p._pnt; }
  bool operator<= (SharedPointer<T> p) const { return _pnt <= p._pnt; }
  bool operator>= (SharedPointer<T> p) const { return _pnt >= p._pnt; }
  
  int count () const { if(_count) return *_count; return 0; }

  friend class ConstSharedPointer<T>;
};

template <typename T>
void SharedPointer<T>::_delete_ref()
{
  if(!_count)
    return;
  if(*_count > 1)
    --(*_count);
  else {
    delete _pnt;
    delete _count;
  }
}

template <typename T>
void SharedPointer<T>::_create_ref(const SharedPointer& s)
{
  _pnt = s._pnt;
  _count = s._count;
  if(_count)
    ++(*_count);
}

template <typename T>
void SharedPointer<T>::init (T* p) 
{
  const char funame [] = " SharedPointer<T>::init: ";

  if(!p) {
    std::cerr << funame << "zero initialization pointer\n";
    throw Error::Range();
  }
  
  if(_count) {
    std::cerr << funame << "already initialized\n";
    throw Error::Init();
  }

  _pnt = p;
  _count = new int(1);
}

template <typename T>
SharedPointer<T>::SharedPointer (T* p) : _pnt(p)
{
  if(p)
    _count = new int(1);
  else
    _count = 0;
}

template <typename T>
T* SharedPointer<T>::operator-> () const  
{ 
  const char funame [] = "SharedPointer::operator->: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

  return _pnt; 
}

template <typename T>
T& SharedPointer<T>::operator*  () const 
{
  const char funame [] = "SharedPointer::operator*: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

 return *_pnt; 
}

template <typename T>
class ConstSharedPointer 
{
private:
  const T*   _pnt;
  int* _count;

  void _delete_ref ();
  void _create_ref (const ConstSharedPointer&);
  void _create_ref (const SharedPointer<T>&);

public:
  explicit ConstSharedPointer (const T* =0);
  ConstSharedPointer (const ConstSharedPointer& s) { _create_ref(s); }
  ConstSharedPointer& operator= (const ConstSharedPointer& s) 
  { _delete_ref(); _create_ref(s); return *this; }
  ConstSharedPointer (const SharedPointer<T>& s) { _create_ref(s); }
  ConstSharedPointer& operator= (const SharedPointer<T>& s) 
  { _delete_ref(); _create_ref(s); return *this; }

  ~ConstSharedPointer () { _delete_ref(); }

  void init (const T*) ;

  operator  bool () const { if(_pnt) return  true; return false; }
  bool operator! () const { if(_pnt) return false; return true;  }

  bool operator== (ConstSharedPointer<T> p) const { return _pnt == p._pnt; }
  bool operator!= (ConstSharedPointer<T> p) const { return _pnt != p._pnt; }
  bool operator<  (ConstSharedPointer<T> p) const { return _pnt <  p._pnt; }
  bool operator>  (ConstSharedPointer<T> p) const { return _pnt >  p._pnt; }
  bool operator<= (ConstSharedPointer<T> p) const { return _pnt <= p._pnt; }
  bool operator>= (ConstSharedPointer<T> p) const { return _pnt >= p._pnt; }
  
  int count () const { if(_count) return *_count; return 0; }

  const T* operator-> () const ;
  const T& operator*  () const ;
  operator const T* () const { return _pnt; }

};

template <typename T>
void ConstSharedPointer<T>::_delete_ref()
{
  if(!_count)
    return;
  if(*_count > 1)
    --(*_count);
  else {
    delete _pnt;
    delete _count;
  }
}

template <typename T>
void ConstSharedPointer<T>::_create_ref(const ConstSharedPointer& s)
{
  _pnt = s._pnt;
  _count = s._count;
  if(_count)
    ++(*_count);
}

template <typename T>
void ConstSharedPointer<T>::_create_ref(const SharedPointer<T>& s)
{
  _pnt = s._pnt;
  _count = s._count;
  if(_count)
    ++(*_count);
}

template <typename T>
void ConstSharedPointer<T>::init (const T* p) 
{
  const char funame [] = " ConstSharedPointer<T>::init: ";

  if(!p) {
    std::cerr << funame << "zero initialization pointer\n";
    throw Error::Range();
  }
  
  if(_count) {
    std::cerr << funame << "already initialized\n";
    throw Error::Init();
  }

  _pnt = p;
  _count = new int(1);
}

template <typename T>
ConstSharedPointer<T>::ConstSharedPointer (const T* p) : _pnt(p)
{
  if(p)
    _count = new int(1);
  else
    _count = 0;
}

template <typename T>
const T* ConstSharedPointer<T>::operator-> () const  
{ 
  const char funame [] = "ConstSharedPointer::operator->: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

  return _pnt; 
}

template <typename T>
const T& ConstSharedPointer<T>::operator*  () const 
{
  const char funame [] = "ConstSharedPointer::operator*: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

 return *_pnt; 
}

#endif
