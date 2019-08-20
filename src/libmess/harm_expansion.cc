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

#include "structure.hh"
#include "harmonic.hh"
#include "key.hh"

#include <fstream>
#include <sstream>

int main (int argc, char* argv [])
{
  const char funame [] = "harm_expansion: ";

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  if (argc < 2) {
    std::cout << funame  << "usage: harm_expansion input_file\n";
    return 1;
  }

  // base name
  std::string base_name = argv[1];
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    base_name.resize(base_name.size() - 4);

  // log file
  IO::log.open((base_name + ".log").c_str());

  // layout
  std::vector<int> vtemp(2);
  vtemp[0] = 3;
  vtemp[1] = 4;
  Configuration::State::layout.set(vtemp);

  // input
  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    throw Error::Input();
  }

  // input parameters
  int rmin = -1, rmax = -1, qmin = -1, qmax = -1;

  KeyGroup MainGroup;
 
  Key struc_key("Structure"  );
  Key  rmin_key("OrbSizeMin" );
  Key  rmax_key("OrbSizeMax" );
  Key  qmin_key("FragSizeMin");
  Key  qmax_key("FragSizeMax");

  std::string token, comment;
  while(from >> token) {
    // molecular structure
    if(struc_key == token) {
      std::getline(from, comment);
      Structure::init(from);
    }
    // orbital expansion dimension minimum
    else if(rmin_key == token) {
      if(rmin > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> rmin)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(rmin < 2) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // orbital expansion dimension maximum
    else if(rmax_key == token) {
      if(rmax > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> rmax)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(rmax < 2) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // fragment expansion dimension minimum
    else if(qmin_key == token) {
      if(qmin > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> qmin)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(qmin < 2) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // fragment expansion dimension maximum
    else if(qmax_key == token) {
      if(qmax > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> qmax)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(qmax < 2) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
  from.close();
  from.clear();

  if(qmin < 2 || qmax < 2 || rmin < 2 || rmax < 2 || qmin > qmax || rmin > rmax) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  // symmetry group
  std::vector<Symmetry::SpaceGroup> sg;
  for(int frag = 0; frag < 2; ++frag)
    sg.push_back(*Structure::fragment(frag).symmetry_group);

  Configuration::DoubleSpaceGroup symm_group(sg);

  // harmonic expansion
  for(int r = rmin; r <= rmax; ++r)
    for(int q = qmin; q <= qmax; ++q) {

      HarmonicExpansion hex(symm_group, r, q);
      std::ostringstream hex_name;
      hex_name << base_name << "_" << r << "_" << q << ".vec";
      
      std::ofstream to(hex_name.str().c_str());
      to << std::setprecision(17) << std::scientific <<  hex;
    }

  return 0.;
}

