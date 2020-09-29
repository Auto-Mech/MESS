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

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<set>

#include "libmess/model.hh"
#include "libmess/key.hh"
#include "libmess/units.hh"
#include "libmess/io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "gumbo: ";

  if (argc < 2) {
    std::cout << "usage: gumbo input_file\n";
    return 0;
  }

  int                 itemp;
  double              dtemp;
  std::string         stemp;

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key     spec_key("Species"                    );
  Key     emax_key("ModelEnergyLimit[kcal/mol]" );
  Key  adm_bor_key("AtomDistanceMin[bohr]"      );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );

  std::vector<SharedPointer<Model::Species> > species;

  // base name
  //
  std::string base_name = argv[1];
  
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    //
    base_name.resize(base_name.size() - 4);

  IO::KeyBufferStream from(argv[1]);
  //
  if(!from && base_name == argv[1]) {
    //
    // try inp extension
    //
    stemp = base_name + ".inp";
    
    from.open(stemp.c_str());
  }

  if(!from) {
    //
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    
    return 1;
  }
  
  //IO::out.open((base_name + ".out").c_str());
  
  IO::log.open((base_name + ".log").c_str());
  
  IO::aux.open((base_name + ".aux").c_str());

  std::string token, comment, line, name;

  std::set<std::string> name_pool;
  
  while(from >> token) {
    //
    // model energy limit
    //
    if(emax_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      Model::set_energy_limit(dtemp * Phys_const::kcal);
    }
    // minimal interatomic distance
    //
    else if(adm_bor_key == token || adm_ang_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
        std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
        throw Error::Range();
      }
      
      if(adm_ang_key == token)
	//
	dtemp *= Phys_const::angstrom;

      Model::atom_dist_min = dtemp;
    }
    // species
    //
    else if(spec_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> name)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      if(!name_pool.insert(name).second) {
	//
	std::cerr << funame << token << ": " << name << ": name already in use\n";

	throw Error::Input();
      }
      
      species.push_back(Model::new_species(from, name, Model::NOSTATES));
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      std::cerr << funame << "unknown keyword: " << token << "\n";
      //
      Key::show_all(std::cerr);
      
      std::cerr << "\n";
      
      throw Error::Init();
    }
  }

  /*************************** CHECKING ******************************************/

  if(!species.size()) {
    std::cerr << funame << "no species\n";
    throw Error::Init();
  }

  return 0;
}
