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

#include "io.hh"

#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>

namespace IO {
  //
  int mpi_rank = 0;

  std::string   _comment_symbol = "#!";
  std::string separator_symbol = ",;";

  const std::string& comment_symbol () { return _comment_symbol; }

  LogOut log;
  LogOut out;
  LogOut aux;

  std::ofstream null("/dev/null");

  Offset log_offset(3);

  log_t _loglevel = NOTICE;
  //
  log_t loglevel () { return _loglevel; }

  std::map<std::string, log_t> _loglevel_map;

  std::string  first_offset = "     ";
  std::string second_offset = "         ";
}

std::string IO::white_space (int n) {
  std::ostringstream res;
  res << std::setw(n) << "";
  return res.str();
}

void IO::toupper (std::string& s)
{
  for(int i = 0; i < s.size(); ++i)
    s[i] = std::toupper(s[i]);
}

void IO::tolower (std::string& s)
{
  for(int i = 0; i < s.size(); ++i)
    s[i] = std::tolower(s[i]);
}

void IO::set_loglevel (const std::string& s)
{
  Exception::Base funame = "IO::set_loglevel: ";
  
  // initialization
  if(!_loglevel_map.size()) {
    _loglevel_map["error"  ] = ERROR; 
    _loglevel_map["warning"] = WARNING; 
    _loglevel_map["notice" ] = NOTICE; 
    _loglevel_map["info"   ] = INFO; 
    _loglevel_map["debug"  ] = DEVEL; 
  }

  if(_loglevel_map.find(s) != _loglevel_map.end()) {
    _loglevel = _loglevel_map[s];
  }
  else {
    funame << "wrong loglevel: " << s << ": available loglevels: ";
    for(std::map<std::string, log_t>::const_iterator it = _loglevel_map.begin(); it != _loglevel_map.end(); ++it) {
      if(it != _loglevel_map.begin())
	funame << ", ";
      funame << it->first;
    }
    throw funame;
  }
}

bool IO::is_break (const std::string& token) {
  //
  static const std::string break_symbol = "@+";

  static const int length_min = 3;

  if(token.size() < length_min)
    //
    return false;

  for(std::string::const_iterator cit = break_symbol.begin(); cit != break_symbol.end(); ++cit) {
    //
    bool btemp = true;
  
    for(int i = 0; i < length_min; ++i)
      //
      if(token[i] != *cit) {
	//
	btemp = false;
	
	break;
      }

    if(btemp)
      //
      return true;
  }

  return false;
}

int IO::skip_comment (const std::string& token, std::istream& from) 
{
  std::string comment;
  
  // check if the first first token character is a comment symbol
  //
  for(int i = 0; i < _comment_symbol.size(); ++i)
    //
    if(token[0] == _comment_symbol[i]) {
      //
      std::getline(from, comment);
      
      return 0;
    }

  return 1;
}

char IO::skip_space(std::istream& from)
{
  Exception::Base funame = "IO::skip_space: ";

  char        next;
  std::string comment;

  while(from.get(next)) {
    bool btemp = false;
    for(int i = 0; i < separator_symbol.size(); ++i)
      if(next == separator_symbol[i]) {
	btemp = true;
	break;
      }

    if(btemp)
      continue;

    if(isspace(next))
      continue;

    btemp = false;
    for(int i = 0; i < _comment_symbol.size(); ++i)
      if(next == _comment_symbol[i]) {
	btemp = true;
	break;
      }

    if(btemp)
      getline(from, comment);
    else
      break;
  }

  if(from.eof())
    return EOF;  

  if(from.fail())
    throw funame << "input stream is corrupted";

  from.putback(next);
  return next;
}

/***********************************************************************************
 ****************************** KEY BUFFER STREAM **********************************
 ***********************************************************************************/

void IO::KeyBufferStream::put_back (const std::string& s) 
{
  const char funame [] = "IO::KeyBufferStream::put_back: ";

  //IO::Marker funame_marker(funame);

  //IO::log << IO::log_offset << s << std::endl;
  
  if(_buffer.size()) {
    std::cerr << funame << "buffer already full\n";
    throw Error::Logic();
  }

  _buffer.push_back(s);
}

template <>
IO::KeyBufferStream& IO::operator>> (KeyBufferStream& from, std::string& s)
{
  const char funame [] = "IO::operator>>(KeyBufferStream&, std::string&): ";

  //IO::Marker funame_marker(funame);
    
  if(from._buffer.size()) {
    //
    s = from._buffer.back();

    //IO::log << IO::log_offset << "pop back: " << s << std::endl;
      
    from._buffer.pop_back();
  }
  else
    //
    (std::ifstream&)from >> s;
  
  return from;
}

/***********************************************************************************
 ****************************** LINE INPUT STREAM **********************************
 ***********************************************************************************/

IO::LineInput::LineInput(std::istream& from)
{
  Exception::Base funame = "IO::LineInput::LineInput: ";

  std::string line;
  if(!std::getline(from, line))
    throw Exception::Eof(funame << "end of file");

  std::istringstream::str(line);
}

void IO::LineInput::read_line(std::istream& from)
{
  Exception::Base funame = "IO::LineInput::read_line: ";

  std::string line;
  if(!std::getline(from, line))
    throw Exception::Eof(funame << "end of file");

  std::istringstream::clear();
  std::istringstream::str(line);
}

/************************************************************************************
 ************************************ INPUT MARKER **********************************
 ************************************************************************************/

void IO::Marker::init(const char* h, int f, std::ostream* out)
{
  const char funame [] = "IO::Marker::init: ";
  
  if(_isinit) {
    //
    std::cerr << funame << "already initialized\n";

    throw Error::Init();
  }

  _isinit = true;
  
  _header = h;

  _start_time = std::time(0);

  _start_cpu = std::clock();

  _flags = f;

  // only master node can print
  //
  if(mpi_rank)
    //
    return;

  if(out) {
    //
    _to = out;
  }
  else if(_flags & NOPRINT) {
    //
    _to = &null;
  }
  else if(log.is_open()) {
    //
    _to = &log;
  }
  else
    //
    _to = &std::cout;

  *_to << log_offset << _header;

  if(_flags & ONE_LINE) {
    //
    *_to << " ..." << std::flush;
  }
  else {
    //
    *_to << " starts" << std::endl;

    log_offset.increase();
  }
}

IO::Marker::~Marker ()
{
  // only master node can print
  //
  if(!_isinit || mpi_rank)
    //
    return;
  
  if(!(_flags & ONE_LINE)) {
    //
    log_offset.decrease();

    *_to << log_offset << _header;
  }

  *_to << " done";

  if(_flags & NOTIME) {
    //
    *_to << "\n";
  }
  else {
    //
    *_to << ", cpu time[sec] = " << double(std::clock() - _start_cpu) / CLOCKS_PER_SEC 
	 << ", elapsed time[sec] = "<< std::time(0) - _start_time
	 << std::endl;
  }
}

/**************************************************************************************************
 ************************************ STRING-TO-NUMBER CONVERTER **********************************
 **************************************************************************************************/

IO::String::operator double() const
{
  const char funame [] = "IO::String::double: ";

  char* endptr;

  double res = std::strtod(c_str(), &endptr);
  
  if(std::strlen(endptr)) {
    //
    ErrOut err_out;
    
    err_out << funame << "conversion failure: residual string: " << endptr;
  }
  
  return res;
}

IO::String::operator int() const
{
  const char funame [] = "IO::String::int: ";

  char* endptr;

  int res = std::strtol(c_str(), &endptr, 10);
  
  if(std::strlen(endptr)) {
    //
    ErrOut err_out;
    
    err_out << funame << "conversion failure: residual string: " << endptr;
  }

  return res;
}

