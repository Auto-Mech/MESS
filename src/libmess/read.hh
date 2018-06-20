/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef READ_HH
#define READ_HH

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "error.hh"
#include "shared.hh"
#include "io.hh"

/*************************************************************************
 *      Base abstract class to read objects from the input stream        *
 *************************************************************************/

class ReadBase
{
public:
  enum State {NOVAL, DEFAULT, READ};// reading status

protected:
  State    _state;

public:
  ReadBase ()      : _state(NOVAL)   {}
  ReadBase (int) : _state(DEFAULT) {}

  State state () const { return _state; }
  virtual void operator() (std::istream&)       throw(Error::General) =0;
  virtual void operator() (std::ostream&) const throw(Error::General) =0;
  virtual ~ReadBase () {}
};

/********************************************************************************
 *                            Class to read an integer                          *
 ********************************************************************************/

class ReadInt : public ReadBase
{
  int& _data;

public:
  explicit ReadInt (int& v) : _data(v) {}
  ReadInt (int& v, int s) : ReadBase(s), _data(v) {}
  void operator() (std::istream&) throw(Error::General);
  void operator() (std::ostream& to) const throw(Error::General);
  ~ReadInt () {}
};

inline void ReadInt::operator() (std::ostream& to) const throw(Error::General)
{
  const char funame [] =  "ReadInt::operator() (std::ostream&): ";

  if(state() != NOVAL)
    to << _data;
  else
    to << funame << "has not been initialized\n";
}

/********************************************************************************
 *                            Class to read long integer                          *
 ********************************************************************************/

class ReadLong : public ReadBase
{
  long& _data;

public:
  explicit ReadLong (long& v) : _data(v) {}
  ReadLong (long& v, long s) : ReadBase(s), _data(v) {}
  void operator() (std::istream&) throw(Error::General);
  void operator() (std::ostream& to) const throw(Error::General);
  ~ReadLong () {}
};

inline void ReadLong::operator() (std::ostream& to) const throw(Error::General)
{
  const char funame [] =  "ReadInt::operator() (std::ostream&): ";

  if(state() != NOVAL)
    to << _data;
  else
    to << funame << "has not been initialized\n";
}

/******************************************************************************
 *                            Class to read a double                          *
 ******************************************************************************/

class ReadDouble : public ReadBase
{
  double& _data;

public:
  explicit ReadDouble (double& v) : _data(v) {}
  ReadDouble (double& v, int s) : ReadBase(s), _data(v) {}
  void operator() (std::istream&) throw(Error::General);    
  void operator() (std::ostream&) const throw(Error::General);
  ~ReadDouble () {}
};

inline void ReadDouble::operator() (std::ostream& to) const throw(Error::General)
{
  const char funame [] =  "ReadDouble::operator() (std::ostream&): ";

  if(state() != NOVAL)
    to << _data;
  else
    to << funame << "has not been initialized\n";
}

/******************************************************************************
 *                            Class to read a string                          *
 ******************************************************************************/

class ReadString : public ReadBase 
{
  std::string& _data;

public:
  explicit ReadString (std::string& v) : _data(v) {}
  ReadString (std::string& v, int s) : ReadBase(s), _data(v) {}
  void operator() (std::istream&) throw(Error::General);    
  void operator() (std::ostream&) const throw(Error::General);
  ~ReadString () {}
};

inline void ReadString::operator() (std::ostream& to) const throw(Error::General)
{
  const char funame [] =  "ReadString::operator() (std::ostream&): ";

  if(state() != NOVAL)
    to << _data;
  else
    to << funame << "has not been initialized\n";
}

/******************************************************************************
 *                       Class to read array of integers                      *
 ******************************************************************************/

class ReadIarr : public ReadBase 
{
  std::vector<int>& _data;

public:
  explicit ReadIarr (std::vector<int>& v) : _data(v) {}
  ReadIarr (std::vector<int>& v, int s) : ReadBase(s), _data(v) {}
  void operator() (std::istream&) throw(Error::General);    
  void operator() (std::ostream&) const throw(Error::General);
  ~ReadIarr () {}
};

inline void ReadIarr::operator () (std::ostream& to) const throw(Error::General)
{
  const char funame [] =  "ReadIarr::operator() (std::ostream&): ";

  if(state() != NOVAL)
    for(int i = 0; i < _data.size(); ++i)
      to << _data[i] << " ";
  else
    to << funame << "has not been initialized\n";

}

/******************************************************************************
 *                        Class to read array of doubles                      *
 ******************************************************************************/

class ReadDarr : public ReadBase 
{
  std::vector<double>& _data;

public:
  explicit ReadDarr (std::vector<double>& v) : _data(v) {}
  ReadDarr (std::vector<double>& v, int s) : ReadBase(s), _data(v) {}
  void operator() (std::istream&) throw(Error::General);    
  void operator() (std::ostream&) const throw(Error::General);
  ~ReadDarr () {}
};

inline void ReadDarr::operator () (std::ostream& to) const throw(Error::General)
{
  const char funame [] =  "ReadDarr::operator() (std::ostream&): ";

  if(state() != NOVAL)
    for(int i = 0; i < _data.size(); ++i)
      to << _data[i] << " ";
  else
    to << funame << "has not been initialized\n";

}

/******************************************************************************
 *                    Class to read with reading functions                    *
 ******************************************************************************/

class ReadFun : public ReadBase 
{
public:
  typedef void (*fun_t) (std::istream&);

private:
    fun_t _data;

public:
  explicit ReadFun (fun_t v) : _data(v) {}
  ReadFun (fun_t v, int s) : ReadBase(s), _data(v) {}
  void operator () (std::istream& from) throw(Error::General)
  { _state = READ; _data(from); }    
  void operator () (std::ostream&) const throw(Error::General);
  ~ReadFun () {}
};

inline void ReadFun::operator () (std::ostream& to) const throw(Error::General) 
{
    if(_state == READ)
	to << "has been processed";
    else
	to << "has not been processed";
}

/******************************************************************************
 *  Class to read to objects which can read themselves from the input stream  *
 ******************************************************************************/

class ReadRead : public ReadBase 
{
  IO::Read* _data;

public:
  explicit ReadRead (IO::Read& v) : _data(&v) {}
  ReadRead (IO::Read& v, int s) : ReadBase(s), _data(&v) {}
  void operator () (std::istream& from) throw(Error::General)
  { _state = READ; _data->read(from); }    
  void operator () (std::ostream&) const throw(Error::General);
  ~ReadRead () {}
};

inline void ReadRead::operator () (std::ostream& to) const throw(Error::General) 
{
    if(_state == READ)
	to << "has been read";
    else
	to << "has not been read";
}

/*****************************************************************************************
 *                                  Wrapper reading class                                *
 *****************************************************************************************/

class Read
{
  SharedPointer<ReadBase> _read;

  // check if the pointer has been initialized
  void _check () const throw(Error::General);

public:
  Read () {}

  explicit Read(int& var) 
    : _read(new ReadInt(var)) {} 
  Read (int& var, int val) 
    : _read(new ReadInt(var, 0)) { var = val; }
 
  explicit Read(long& var) 
    : _read(new ReadLong(var)) {} 
  Read (long& var, long val) 
    : _read(new ReadLong(var, 0)) { var = val; }
 
  explicit Read(double& var) 
    : _read(new ReadDouble(var)) {} 
  Read(double& var, double val) 
    : _read(new ReadDouble(var, 0)) { var = val; }
  
  explicit Read(std::string& var) 
    : _read(new ReadString(var)) {} 
  Read (std::string& var, const std::string& val) 
    : _read(new ReadString(var, 0)) { var = val; }

  explicit Read(std::vector<double>& var) 
    : _read(new ReadDarr(var)) {} 
  Read (std::vector<double>& var, const std::vector<double>& val) 
    : _read(new ReadDarr(var, 0)) { var = val; }
  
  explicit Read(std::vector<int>& var) 
    : _read(new ReadIarr(var)) {} 
  Read (std::vector<int>& var, const std::vector<int>& val) 
    : _read(new ReadIarr(var, 0)) { var = val; }
  
  explicit Read(ReadFun::fun_t var) 
    : _read(new ReadFun(var)) {} 
  Read  (ReadFun::fun_t var, int) 
    : _read(new ReadFun(var, 0)) {} 

  explicit Read(IO::Read& var) 
    : _read(new ReadRead(var)) {} 
  Read  (IO::Read& var, int) 
    : _read(new ReadRead(var, 0)) {} 

  void operator () (std::ostream& o) const
    throw(Error::General) { _check(); (*_read)(o); }
  void operator () (std::istream& i)
    throw(Error::General) { _check(); (*_read)(i); }

  bool is_read    () const throw(Error::General);
  bool is_init    () const throw(Error::General);
  bool is_default () const throw(Error::General);

};

inline void Read::_check () const throw(Error::General)
{
  const char funame []  = "Read::_check: ";

  if(!_read) {
    std::cerr << funame << "Read object should be initialized with the non-default constructor\n";
    throw Error::Init();
  }

}

inline std::ostream& operator<< (std::ostream& s, const Read& r) throw(Error::General) { r(s); return s; }
inline std::istream& operator>> (std::istream& s,       Read& r) throw(Error::General) { r(s); return s; }

inline bool Read::is_read () const throw(Error::General)
{
  _check();

  if(_read->state() == ReadBase::READ)
    return true;
  else
    return false;
}

inline bool Read::is_init () const throw(Error::General)
{
  _check();

  if(_read->state() == ReadBase::NOVAL)
    return false;
  else
    return true;
}

inline bool Read::is_default () const throw(Error::General)
{
  _check();

  if(_read->state() == ReadBase::DEFAULT)
    return true;
  else
    return false;
}

#endif
