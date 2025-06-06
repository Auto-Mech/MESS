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

#include "libmess/model.hh"
#include "libmess/key.hh"
#include "libmess/units.hh"
#include "libmess/io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "abstraction: ";

  if (argc < 2) {
    std::cout << "usage: abstraction input_file\n";
    return 0;
  }

  int                 itemp;
  double              dtemp;
  std::string         stemp;

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key temp_key("TemperatureList[K]"        );
  Key  bar_key("Barrier"                   );
  Key  bim_key("Bimolecular"               );
  Key rate_key("RateOutput"                );
  Key  log_key("LogOutput"                 );
  Key emax_key("ModelEnergyLimit[kcal/mol]");

  Key  adm_bor_key("AtomDistanceMin[bohr]"      );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );

  std::vector<double> temperature;
  SharedPointer<Model::Species> barrier;
  std::vector<SharedPointer<Model::Bimolecular> > product;

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

  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string token, comment, line, name;
  while(from >> token) {
    // barrier
    if(bar_key == token) {
      if(!(from >> name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      barrier = Model::new_species(from, name, Model::NOSTATES);
    }
    // bimolecular
    else if(bim_key == token) {
      if(!(from >> name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      product.push_back(Model::new_bimolecular(from, name));
    }
    // rate output    
    else if(rate_key == token) {
      if(IO::out.is_open()) {
        std::cerr << funame << token << ": output stream allready open\n";
        throw Error::Input();
      }      
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      IO::out.open(stemp.c_str());
      if(!IO::out) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // log output
    else if(log_key == token) {
      if(IO::log.is_open()) {
        std::cerr << funame << token << ": log stream allready open\n";
        throw Error::Input();
      }      

      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      IO::log.open(stemp.c_str());
      if(!IO::log) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // temperature
    else if(temp_key == token) {
      if(temperature.size()) {
	std::cerr << funame << "temperature list has been already defined\n";
	throw Error::Input();
      }
      IO::LineInput data_stream(from);
      while(data_stream >> dtemp)
	temperature.push_back( dtemp * Phys_const::kelv);
      
      if(!temperature.size()) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
    }
    // model energy limit
    else if(emax_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      Model::set_energy_limit(dtemp * Phys_const::kcal);
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
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
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

  if(!barrier) {
    std::cerr << funame << "no barrier\n";
    throw Error::Init();
  }

  if(product.size() != 2) {
    std::cerr << funame << "no products\n";
    throw Error::Init();
  }

  // default log output	
  if(!IO::log.is_open()) {
    stemp = base_name + ".log";
    IO::log.open(stemp.c_str());
    if(!IO::log) {
      std::cerr << funame << token << ": cannot open " << stemp << " file\n";
      throw Error::Input();
    }
  }

  // default rate output
  if(!IO::out.is_open()) {
    stemp = base_name + ".out";
    IO::out.open(stemp.c_str());
    if(!IO::out) {
      std::cerr << funame << token << ": cannot open " << stemp << " file\n";
      throw Error::Input();
    }
  }

  // rate calculation and output

  typedef std::vector<double>::const_iterator It;
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;
  IO::out << "Rate Units: cm^3/sec\n\n";
  IO::out << std::setw(5) << "T, K";
  for(int p = 0; p < 2; ++p) {
    stemp = product[p]->name() + "->" + product[1-p]->name();
    IO::out << std::setw(13) << stemp;
  }
  IO::out << "\n";
  
  for(It t = temperature.begin(); t != temperature.end(); ++t) {
    IO::out << std::setw(5) << *t / Phys_const::kelv;;
    for(int p = 0; p < 2; ++p)
      if(!product[p]->dummy()) {
	dtemp = *t / 2. / M_PI * barrier->weight(*t) / product[p]->weight(*t)
	  * std::exp((product[p]->ground() - barrier->ground()) / *t) / bru;
	IO::out << std::setw(13) << dtemp;
      }
      else 
	IO::out << std::setw(13) << "***";
    IO::out << "\n";
  }
  return 0;
}
