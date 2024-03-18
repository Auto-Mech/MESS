/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2017, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef IO_HH
#define IO_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <iomanip>
#include <vector>
#include <ctime>

#include "error.hh"

namespace IO {
  //
  extern int mpi_rank;

  const std::string& comment_symbol ();

  std::string white_space (int);

  void toupper (std::string&);
  
  void tolower (std::string&);

  inline std::string end_key () { return "End"; }

  bool is_break (const std::string&);

  char skip_space(std::istream&);
  
  int skip_comment (const std::string&, std::istream&);

  // verbosity level
  //
  enum log_t {ERROR, WARNING, NOTICE, INFO, DEVEL};

  log_t loglevel     ();
  
  void  set_loglevel (const std::string&);

  /************************************************************************
   ****************************** OUTPUT **********************************
   ************************************************************************/

  class LogOut : public std::ofstream {};

  template <typename T>
  //
  LogOut& operator<< (LogOut& to, const T& t)
  {
    // only the master node can print
    //
    if(mpi_rank)
      //
      return to;

    if(to.is_open()) {
      //
      (std::ostream&)to << t;
    }
    else
      //
      std::cout << t;

    return to;
  }

  extern LogOut log;
  
  extern LogOut out;

  extern LogOut aux;
  
  extern std::ofstream null;
  
  /***********************************************************************************
   ****************************** LINE INPUT STREAM **********************************
   ***********************************************************************************/

  class LineInput : public std::istringstream {
    //
    LineInput();
    
    LineInput(const LineInput&);
    
    LineInput& operator=(const LineInput&);

  public:
    //
    LineInput(std::istream& from);

    void read_line(std::istream& from);
  };

  /*************************************************************************************
   ****************************** ABSTRACT READ CLASS **********************************
   *************************************************************************************/

  // abstract class which objects read from the input stream
  //
  class Read {
    //
  public:
    //
    // read object from input stream
    //
    virtual void read (std::istream&) =0;
    
    virtual ~Read () {}
  };

  inline std::istream& operator>> (std::istream& from, Read& r) { r.read(from); return from; }

  /*************************************************************************
   ****************************** OFFSETS **********************************
   *************************************************************************/

  extern std::string    first_offset;
  
  extern std::string   second_offset;

  // numerical offset
  //
  class Offset {
    //
    int _value;
    
    int _step;

  public:
    //
    explicit Offset (int s) : _value(0), _step(s) {}

    void increase () { _value += _step; }
    
    void decrease () { _value -= _step; }

    operator int () const { return _value; }
  };

  inline std::ostream& operator<< (std::ostream& to, const Offset& off)
  {
    return to << std::setw(off) << "";
  }

  extern Offset log_offset;

  /***********************************************************************************
   ****************************** KEY BUFFER STREAM **********************************
   ***********************************************************************************/

  // facility to put keyword back to the input stream
  //
  class KeyBufferStream : public std::ifstream {
    //
    std::vector<std::string> _buffer;

  public:
    //
    KeyBufferStream (const char* f) : std::ifstream(f) {}
    
    KeyBufferStream () {}

    void put_back (const std::string&) ;

    template <typename T>
    //
    friend KeyBufferStream& operator>> (KeyBufferStream& from, T& t);
  };

  template <typename T>
  //
  KeyBufferStream& operator>> (KeyBufferStream& from, T&  t) { (std::ifstream&)from >> t; return from; }

  template <>
  //
  KeyBufferStream& operator>> (KeyBufferStream& from, std::string& t);
  

  /************************************************************************************
   ************************************ INPUT MARKER **********************************
   ************************************************************************************/

  class Marker {
    //
    std::string    _header;

    std::clock_t  _start_cpu;

    std::time_t   _start_time;

    int           _flags;

    std::ostream* _to;

    bool          _isinit;

  public:
    //
    Marker () : _isinit(false) {}
    
    void init (const char*, int =0, std::ostream* =0);
    
    Marker    (const char*, int =0, std::ostream* =0);
    
    ~Marker ();
    
    enum {
      NOTIME   = 1,
      ONE_LINE = 2,
      NOPRINT  = 4
    };
  };

  inline Marker::Marker (const char* h, int f, std::ostream* o) : _isinit(false) { init(h, f, o); }
  
  /****************************************************************************************
   ************************************ STRING CONVERTER **********************************
   ****************************************************************************************/

  class String : public std::string {

  public:
    //
    String () {}
    
    template<typename T>
    //
    String(const T& t) : std::string(t) {}

    explicit String (int i);
    
    operator double() const;
    
    operator    int() const;
  };

  inline String::String (int i)
  {
    std::ostringstream oss;

    oss << i;

    *this = oss.str();
  }
  //
  //
} // IO

/************************************************************************************
 ************************************ ERROR OUTPUT **********************************
 ************************************************************************************/

class ErrOut : public std::ostringstream {
  //
  ErrOut         (const ErrOut&);

  void operator= (const ErrOut&);

public:
  //
  ErrOut ();

  ~ErrOut ();
};

template<typename T>
//
ErrOut& operator<< (ErrOut& err_out,const T& t) 
{ 
  (std::ostringstream&)err_out << t; 

  return err_out;
}
  
inline ErrOut::ErrOut () 
{
  *this << "Error: ";
}

inline ErrOut::~ErrOut () 
{ 
  IO::log << str() << std::endl; 

  std::cerr << str() << std::endl; 

  //throw Error::General();
  
  std::terminate();
}

#endif
