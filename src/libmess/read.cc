/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include <sstream>
#include <cmath>
#include <cstdio>

#include "read.hh"
#include "units.hh"
#include "io.hh"

void ReadInt::operator () (std::istream& from) throw(Error::General)
{
  const char funame [] = "ReadInt::operator () (std::istream&): ";

  _state = READ;

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cerr << funame << "cannot read integer\n";
    throw Error::Form();
  }

  std::string comment;
  iss >> comment;
  if(iss && comment[0] != '#') {
    std::cerr << funame << "comment does not start with #\n";
    throw Error::Form();
  }
}

void ReadLong::operator () (std::istream& from) throw(Error::General)
{
  const char funame [] = "ReadLong::operator () (std::istream&): ";

  _state = READ;

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cerr << funame << "cannot read integer\n";
    throw Error::Form();
  }

  std::string comment;
  iss >> comment;
  if(iss && comment[0] != '#') {
    std::cerr << funame << "comment does not start with #\n";
    throw Error::Form();
  }
}

void ReadIarr::operator () (std::istream& from) throw(Error::General)
{
  const char funame [] = "ReadIarr::operator () (std::istream&): ";
  
  _state = READ;

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  int itemp;
  while(iss >> itemp)
    _data.push_back(itemp);

  if(!_data.size()) {
    std::cerr << funame << "cannot read integer array\n";
    throw Error::Form();
  }
}

void ReadDouble::operator () (std::istream& from) throw(Error::General)
{
  const char funame [] = "ReadDouble::operator () (std::istream&): ";

  _state = READ;

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cerr << funame << "cannot read double\n";
    throw Error::Form();
  }

  std::string unit;
  iss >> unit;
  if(iss && unit[0] !='#')
    _data *= Phys_const::str2fac(unit);
}

void ReadString::operator () (std::istream& from) throw(Error::General)
{
  const char funame [] = "ReadString::operator () (std::istream&): ";

  _state = READ;

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  iss >> _data;
  if(!iss) {
    std::cerr << funame << "cannot read string\n";
    throw Error::Form();
  }
}

void ReadDarr::operator () (std::istream& from) throw(Error::General)
{
  const char funame [] = "ReadDarr::operator () (std::istream&): ";

  _state = READ;

  std::string line;
  std::getline(from, line);
  std::istringstream iss(line);  

  char next = IO::skip_space(iss);

  if(next == EOF) {
    std::cerr << funame << "no array specification is found\n";
    throw Error::Form();
  }

  // progression with variable step (old fashion input)
  if(next == '.' || next == '-' || next == '+' || isdigit(next)) { 
    double d, dd, incr; 
    int n;

    iss >> d >> dd >> incr >> n;

    if(!iss) {
      std::cerr << funame << "old style array specification is corrupted\n";
      throw Error::Form();
    }

    std::string unit;
    iss >> unit;
    if(iss && unit[0] != '#') {
      double fac = Phys_const::str2fac(unit);
      d *= fac;
      dd *= fac;
    }

    double val = d;
    double step = dd;
    for(int i = 0; i < n; ++i) {
      _data.push_back(val);
      val += step;
      step *= incr;
    }
    return;
  }
  
  // new style input
  while(next == '@' || next  == '*') {
    bool is_geom = true;;
    if(next == '@')// arithmetic progression
      is_geom = false;

    double term;
    iss >> term;
    if(!iss) {
      std::cerr << funame << "cannot read the progression first term\n";
      throw Error::Form();
    }
    
    bool is_step = false;
    bool is_size = false;
    bool is_end  = false;
    double step = 0.;
    double  end = 0.;
    int    size = 0;

    for(int i = 0; i < 2; ++i) {
      next = IO::skip_space(iss);
      if(next == EOF) {
	std::cerr << funame << "unexpected end of stream\n";
	throw Error::Form();
      }
      iss.get();
      switch(toupper(next)) {
      case 'N':
	if(is_size) {
	  std::cerr << "the progression size has been allready defined\n";
	  throw Error::Form();
	}
	is_size = true;
	iss >> size;
	if(!iss) {
	  std::cerr << funame << "cannot read the progression size\n";
	  throw Error::Form();
	}
	break;
      case 'E':
	if(is_end) {
	  std::cerr << "the progression last term has been allready defined\n";
	  throw Error::Form();
	}
	is_end = true;
	iss >> end;
	if(!iss) {
	  std::cerr << funame << "cannot read the progression last term\n";
	  throw Error::Form();
	}	
	break;
      case 'S':
	if(is_step) {
	  std::cerr << "the progression step has been allready defined\n";
	  throw Error::Form();
	}
	is_step = true;
	iss >> step;
	if(!iss) {
	  std::cerr << funame << "cannot read the progression step\n";
	  throw Error::Form();
	}
	break;
      default:
	std::cerr << funame << "unknown control character: " << next << "\n";
	throw Error::Form();
      }
    }
    
    if(!is_step && size > 1) {
      if(is_geom)
	step = exp(log(end/term)/(size-1));
      else
	step = (end-term)/(size-1);
    }
    if(!is_size) {
      if(is_geom) {
	if(step == 1.) {
	  std::cerr << funame << "infinite geometrical progression\n";
	  throw Error::Form();
	}
	size = int(log(end/term) / log(step)) + 1;
      }
      else {
	if(step == 0.) {
	  std::cerr << funame << "infinite arithmetical progression\n";
	  throw Error::Form();
	}
	size = int((end-term)/step) + 1;
      }
    }

    for(int i = 0; i < size; ++i) {
      _data.push_back(term);
      if(is_geom)
	term *= step;
      else
	term += step;
    }
  }

  if(next == EOF || next == '#')
    return;
  
  std::string unit;
  iss >> unit;
  double fac = Phys_const::str2fac(unit);
  for(int i = 0; i < _data.size(); ++i)
    _data[i] *= fac;
}

