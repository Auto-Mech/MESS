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

#include "model.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "mc_test: ";

  if (argc < 2) {
    std::cout << "usage: mc_test input_file\n";
    return 0;
  }

  int                 itemp;
  double              dtemp;
  std::string         stemp;

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key spec_key("Species"                  );
  Key temp_key("TemperatureList[K]"       );
  Key freq_key("Frequencies[1/cm]"        );
  Key   rc_key("RotationalConstants[1/cm]");
  
  Key emax_key("ModelEnergyLimit[kcal/mol]"   );
  Key  adm_bor_key("AtomDistanceMin[bohr]"    );
  Key  adm_ang_key("AtomDistanceMin[angstrom]");


  std::vector<double> temperature;
  SharedPointer<Model::MonteCarlo> species;
  std::vector<double> freq, rc;

  double temp_rel_incr = 0.001;

  // base name
  std::string base_name = argv[1];
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    base_name.resize(base_name.size() - 4);

  IO::KeyBufferStream from(argv[1]);
  if(!from && base_name == argv[1]) {
    // try inp extension
    stemp = base_name + ".inp";
    from.open(stemp.c_str());
  }

  IO::log.open((base_name + ".log").c_str());
  IO::out.open((base_name + ".out").c_str());

  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string token, comment, line, name;

  while(from >> token) {
    //
    // model energy limit
    //
    if(emax_key == token) {
      //
      if(!(from >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";

        throw Error::Input();
      }

      std::getline(from, comment);

      Model::set_energy_limit(dtemp * Phys_const::kcal);
    }
    // minimal interatomic distance
    //
    else if(adm_bor_key == token || adm_ang_key == token) {
      //
      if(!(from >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";

        throw Error::Input();
      }

      std::getline(from, comment);

      if(dtemp <= 0) {
	//
        std::cerr << funame << token << ": out of range\n";

        throw Error::Range();
      }
      
      if(adm_ang_key == token)
	//
	dtemp *= Phys_const::angstrom;

      Model::atom_dist_min = dtemp;
    }
    // monte-carlo species
    //
    else if(spec_key == token) {
      //
      if(species) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      IO::LineInput lin(from);

      if(!(lin >> name)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      species.init(new Model::MonteCarlo(from, name, Model::NOSTATES));
    }
    // temperatures list
    //
    else if(temp_key == token) {
      //
      if(temperature.size()) {
	//
	std::cerr << funame << "temperature list has been already defined\n";
	
	throw Error::Input();
      }
      
      IO::LineInput line_input(from);
      
      while(line_input >> dtemp)
	//
	temperature.push_back( dtemp * Phys_const::kelv);
      
      if(!temperature.size()) {
	//
        std::cerr << funame << token << ": corrupted\n";

        throw Error::Input();
      }
    }
    // frequencies at minimum
    //
    else if(freq_key == token) {
      //
      if(freq.size()) {
	//
	std::cerr << funame << token << ": already defined\n";
	
	throw Error::Input();
      }
      
      IO::LineInput line_input(from);
      
      while(line_input >> dtemp)
	//
	freq.push_back(dtemp * Phys_const::incm);
      
      if(!freq.size()) {
	//
        std::cerr << funame << token << ": corrupted\n";

        throw Error::Input();
      }
    }
    // rotational constants at minimum
    //
    else if(rc_key == token) {
      //
      if(rc.size()) {
	//
	std::cerr << funame << token << ": already defined\n";
	
	throw Error::Input();
      }
      
      IO::LineInput line_input(from);
      
      while(line_input >> dtemp) {
	//
	if(dtemp <= 0.) {
	  //
	  std::cerr << funame << token << ": " << rc.size() + 1 << "-th rotational constant out of range: " << dtemp << "\n";

	  throw Error::Range();
	}
	
	rc.push_back(dtemp * Phys_const::incm);
      }
      
      if(rc.size() != 3) {
	//
        std::cerr << funame << token << ": # of rotational constants out of range: " << rc.size() << "\n";

        throw Error::Input();
      }
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      std::cerr << funame << "unknown keyword: " << token << "\n";

      Key::show_all(std::cerr);

      std::cerr << "\n";
      
      throw Error::Init();
    }
  }

  /*************************** CHECKING ******************************************/

  if(!temperature.size()) {
    //
    std::cerr << funame << "temperature list has not been initialized\n";

    throw Error::Init();
  }

  if(!species) {
    //
    std::cerr << funame << "no species\n";

    throw Error::Init();
  }

  dtemp = std::sqrt(M_PI);

  for(int i = 0; i < rc.size(); ++i)
    //
    dtemp /= std::sqrt(rc[i]);

  for(int i = 0; i < freq.size(); ++i)
    //
    dtemp /= freq[i];

  const double ho_fac = dtemp;
  
  /***************** PARTITION FUNCTION CALCULATION AND OUTPUT *******************/

  IO::out << std::setw(13) << std::left << "T, K" << std::right
	  << std::setw(13) << species->name()
	  << std::setw(13) << "Error, %";

  if(rc.size() && freq.size())
    //
    IO::out << std::setw(13) << "Z/Z_HO";
  
  IO::out << "\n";
  
  for(int t = 0; t < temperature.size(); ++t) {
    //
    double err;
    double z = species->weight_with_error(temperature[t], err);
    
    IO::out << std::setw(13) << std::left << temperature[t] / Phys_const::kelv << std::right
	    << std::setw(13) << z
	    << std::setw(13) << err;
    
    if(rc.size() && freq.size())
      //
      IO::out << std::setw(13) << z / ho_fac / std::pow(temperature[t], (double)freq.size() + (double)rc.size() / 2.);
    
    IO::out << "\n";
  }
  
  return 0;
}
