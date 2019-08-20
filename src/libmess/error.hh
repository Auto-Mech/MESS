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

#ifndef ERROR_HH
#define ERROR_HH

#include <string>
#include <iostream>
#include <sstream>

namespace Exception {

  class Base {

    std::string _mesg;

  public:
    Base () {}
    Base (const Base& b)        : _mesg(b._mesg)  {}
    Base (const std::string& s) : _mesg(s)        {}
    Base (const char*        s) : _mesg(s)        {}
    ~Base ()  {}

    template<typename T>
    Base& operator<< (const T&);

    operator std::string () const  { return _mesg;         }
    operator const char* () const  { return _mesg.c_str(); }

    int size () const { return _mesg.size() + 1; }

    void print (std::ostream& =std::cerr) const;
    void add_front(const std::string& s) { _mesg = s + _mesg; }

    friend std::ostream& operator<< (std::ostream&, const Base&);
  };

  inline std::ostream& operator<< (std::ostream& out, const Base& e) 
  {
    return out << e._mesg;
  }

  inline void Base::print (std::ostream& out) const
  {
    out << _mesg << "\n";
  }

  template<typename T>
  Base& Base::operator<< (const T& t)
  {
    std::ostringstream out;
    out << _mesg << t;
    _mesg = out.str();
    return *this;
  }

  // end of file
  class Eof : public Base {
  public:
    Eof () {}
    Eof (const Eof&  e) : Base(e) {}
    Eof (const Base& b) : Base(b) {}

    template<typename T>
    Eof& operator<< (const T& t) { (Base&)(*this) << t; return *this; }
  };

  // range
  class Range : public Base {
  public:
    Range () {}
    Range (const Range&  r) : Base(r) {}
    Range (const Base&   b) : Base(b) {}

    template<typename T>
    Range& operator<< (const T& t) { (Base&)(*this) << t; return *this; }
  };

  // open file 
  class Open : public Base {
  public:
    Open () {}
    Open (const Open& o) : Base(o) {}
    Open (const Base& b) : Base(b) {}

    template<typename T>
    Open& operator<< (const T& t) { (Base&)(*this) << t; return *this; }
  };

}// exception namespace

// Error classes
namespace Error {
  class General {};

  // File handling errors
  class  File: public General {};
  class  Open: public File{}; // File opening errors
  class  Form: public File{}; // Wrong format
  class Input: public File{}; // Wrong input
  class   Eof: public File{}; // End-of-file encountered

  // Math errors
  class Math: public General {};
  // Indexing errors
  class Range: public General {};
  // Search errors
  class Find: public General {};
  // Initialization errors
  class Init: public General {};
  // Run errors
  class Run: public General {};
  // Logic errors
  class Logic: public General {};

  class Molpro: public General, private std::string {

  public:
    Molpro () {}
    explicit Molpro (const std::string& s) : std::string(s) {}

    friend std::ostream& operator<< (std::ostream&, const Molpro&);
  };

  inline std::ostream& operator<< (std::ostream& to, const Molpro& m)
  {
    return to << (const std::string&)m;
  }
}

#endif
