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

#include <fstream>

#include "atom.hh"
#include "key.hh"
#include "io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "main: ";

  std::string stemp;
  double      dtemp;
  bool        btemp;
  int         itemp;

  if (argc < 2) {
    std::cerr << "usage: symmetry_number input_file\n";
    return -1;
  }
  
  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    throw Error::Input();
  }

  KeyGroup SymmetryNumberGroup;

  Key geom_key("Geometry");             // molecular geometry
  Key isot_key("IgnoreIsotope");          // 180-degree rotation around flip axis in the standard orientation
  Key unit_key("Angstrom");             // use angstrom as distance units
  Key  tol_key("Tolerance");             // use angstrom as distance units

  std::string token, comment;
  int flags = 0;
  bool angstrom = false;
  std::vector<Atom> molecule;
  double tolerance = 1.e-3;

  while(from >> token) {
    // geometry section
    if(token == geom_key) {
      IO::LineInput lin(from);
      if(!(lin >> itemp)) {
	std::cerr << funame << "could not read number of atoms\n";
	throw Error::Form();
      }

      molecule.resize(itemp);
      for(int a = 0; a < molecule.size(); ++a)
	from >> molecule[a];
    }
    // angstrom
    else if(token == unit_key) {
      angstrom = true;
      std::getline(from, comment);
    }
    //  geometry tolerance
    else if(token == tol_key) {
      if(!(from >> tolerance)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(tolerance > 1. || tolerance <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // ignore isotope
    else if(token == isot_key) {
      flags |= IGNORE_ISOTOPE;
      std::getline(from, comment);
    }
    // unknown key
    else {
      std::cerr << funame << "unknown key: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }  
  }

  if(!molecule.size()) {
    std::cerr << funame << "molecule not defined\n";
    throw Error::Init();
  }

  if(angstrom)
    for(int a = 0; a < molecule.size(); ++a)
      molecule[a] *= Phys_const::angstrom;

  std::pair<int, int> res = symmetry_number(molecule, tolerance, flags);

  std::cout << res.first << "\t" << res.second << "\n";

  return 0;
}
