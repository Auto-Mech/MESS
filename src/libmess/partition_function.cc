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
  const char funame [] = "partition_function: ";

  if (argc < 2) {
    std::cout << "usage: partition_function input_file\n";
    return 0;
  }

  int                 itemp;
  double              dtemp;
  std::string         stemp;

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key spec_key("Species"                   );
  Key list_key("TemperatureList[K]"        );
  Key grid_key("Temperature(step[K],size)" );
  Key emax_key("ModelEnergyLimit[kcal/mol]");
  Key tincr_key("RelativeTemperatureIncrement");
  Key  adm_bor_key("AtomDistanceMin[bohr]"      );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );


  std::vector<double> temperature;
  std::vector<SharedPointer<Model::Species> > species;

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
  IO::out.open((base_name + ".dat").c_str());

  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string token, comment, line, name;
  while(from >> token) {
    // model energy limit
    if(emax_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      Model::set_energy_limit(dtemp * Phys_const::kcal);
    }
    // relative temperature increment
    else if(tincr_key == token) {
      if(!(from >> temp_rel_incr)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      if(temp_rel_incr <= 0. || temp_rel_incr >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // minimal interatomic distance
    else if(adm_bor_key == token || adm_ang_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0) {
        std::cerr << funame << token << ": out of range\n";
        throw Error::Range();
      }
      
      if(adm_ang_key == token)
	dtemp *= Phys_const::angstrom;

      Model::atom_dist_min = dtemp;
    }
    // species
    else if(spec_key == token) {
      if(!(from >> name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      species.push_back(Model::new_species(from, name, Model::NOSTATES));
    }
    // temperature list
    else if(list_key == token) {
      if(temperature.size()) {
	std::cerr << funame << "temperature list has been already defined\n";
	throw Error::Input();
      }
      IO::LineInput line_input(from);
      while(line_input >> dtemp)
	temperature.push_back( dtemp * Phys_const::kelv);
      
      if(!temperature.size()) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
    }
    // temperature grid
    else if(grid_key == token) {
      if(temperature.size()) {
	std::cerr << funame << "temperature list has been already defined\n";
	throw Error::Input();
      }

      IO::LineInput line_input(from);
      if(!(line_input >> dtemp >> itemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }

      if(dtemp <= 0. || itemp <= 0) {
        std::cerr << funame << token << ": out of range\n";
        throw Error::Range();
      }

      dtemp *=  Phys_const::kelv;
      for(int i = 0; i < itemp; ++i)
	temperature.push_back(double(i + 1) * dtemp);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  /*************************** CHECKING ******************************************/

  if(!temperature.size()) {
    std::cerr << funame << "temperature list has not been initialized\n";
    throw Error::Init();
  }

  if(!species.size()) {
    std::cerr << funame << "no species\n";
    throw Error::Init();
  }

  // add room temperature
  temperature.push_back(298.2 * Phys_const::kelv);

  /***************** PARTITION FUNCTION CALCULATION AND OUTPUT *******************/

  const double volume_unit = Phys_const::cm * Phys_const::cm * Phys_const::cm;

  //  IO::out << "Partition function (relative to the ground,1/cm^3):\n"
  IO::out << "Partition function (log) and its derivatives:\n"
	  << std::left << std::setw(5) << "T, K" << std::right;
  for(int s = 0; s < species.size(); ++s)
    IO::out << std::setw(13) << species[s]->name() << std::setw(26);
  IO::out << "\n";
  
  for(int t = 0; t < temperature.size(); ++t) {

    const double tval = temperature[t];

    double temp_incr = tval * temp_rel_incr;

    double tt[3];
    tt[0] = tval;
    tt[1] = tval - temp_incr;
    tt[2] = tval + temp_incr;

    temp_incr /= Phys_const::kelv;

    double zz[3];

    IO::out << std::left << std::setw(5) << temperature[t] / Phys_const::kelv << std::right; 
    for(int s = 0; s < species.size(); ++s) {
      
      for(int i = 0; i < 3; ++i) {
	//
	dtemp = species[s]->weight(tt[i]) * std::pow(species[s]->mass() * tt[i] / 2. / M_PI, 1.5) * volume_unit;

	if(dtemp <= 0.) {
	  //
	  ErrOut err_out;

	  err_out << funame << "negative weight: " << species[s]->weight(tt[i]);
	}
	
	zz[i] = std::log(dtemp);
      }
      
      IO::out << std::setw(13) << zz[0]
	      << std::setw(13) << (zz[2] - zz[1]) / 2. / temp_incr
	      << std::setw(13) << (zz[2] + zz[1] - 2. * zz[0]) / temp_incr / temp_incr;

      //IO::out << std::setw(13) << species[s]->weight(temperature[t]) 
      //* std::exp((species[s]->real_ground() - species[s]->ground())/ temperature[t])
      //* std::pow(species[s]->mass() * temperature[t] / 2. / M_PI, 1.5) * volume_unit;
    }
    IO::out << "\n";
  }

  return 0;
}
