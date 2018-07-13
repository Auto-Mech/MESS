

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

  void init (T*) throw(Error::General);

  T* operator-> () const throw(Error::General);
  T& operator*  () const throw(Error::General);
  operator T*   () const { return _pnt; }

  operator bool  () const { if(_pnt) return  true; return false; }
  bool operator! () const { if(_pnt) return false; return true;  }

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
void SharedPointer<T>::init (T* p) throw(Error::General)
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
T* SharedPointer<T>::operator-> () const  throw(Error::General)
{ 
  const char funame [] = "SharedPointer::operator->: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

  return _pnt; 
}

template <typename T>
T& SharedPointer<T>::operator*  () const throw(Error::General)
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

  void init (const T*) throw(Error::General);

  operator  bool () const { if(_pnt) return  true; return false; }
  bool operator! () const { if(_pnt) return false; return true;  }

  int count () const { if(_count) return *_count; return 0; }

  const T* operator-> () const throw(Error::General);
  const T& operator*  () const throw(Error::General);
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
void ConstSharedPointer<T>::init (const T* p) throw(Error::General)
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
const T* ConstSharedPointer<T>::operator-> () const  throw(Error::General)
{ 
  const char funame [] = "ConstSharedPointer::operator->: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

  return _pnt; 
}

template <typename T>
const T& ConstSharedPointer<T>::operator*  () const throw(Error::General)
{
  const char funame [] = "ConstSharedPointer::operator*: ";

  if(!_pnt) {
    std::cerr << funame << "has not been initialized\n";
    throw Error::Init();
  }

 return *_pnt; 
}

#endif
