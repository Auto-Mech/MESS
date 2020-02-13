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

#include "model.hh"
#include "key.hh"
#include "io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "multirotor_test: ";

  KeyGroup Main;

  Key  temp_key("TemperatureList[K]"  );
  Key  ener_key("EnergyList[kcal/mol]");
  Key  geom_key("Geometry[angstrom]"  );
  Key   log_key("LogOutput"           );
  Key multi_key("MultiRotor"          );
  Key  hind_key("HinderedRotor"       );
  Key  umbr_key("Umbrella"            );

  std::vector<double> temperature;
  std::map<std::string, Model::MultiRotor*>        multi;
  std::map<std::string, Model::HinderedRotor*>  hindered;
  std::map<std::string, Model::Umbrella>        umbrella;
  
  if (argc < 2) {
    std::cout << "usage: multirotor_test input_file\n";
    return 0;
  }

  int    itemp;
  double dtemp;
  bool   btemp;

  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string token, comment, line, stemp;
  
  while(from >> token) {
    // temperature
    if(temp_key == token) {
      if(temperature.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      IO::LineInput temperature_input(from);
      while(temperature_input >> dtemp)
	temperature.push_back( dtemp * Phys_const::kelv);
    }
    // multi-rotor
    else if(mrot_key == token) {
      if(multi_rotor) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      multi_rotor = new Model::MultiRotor(from);
    }
    // log output
    else if(log_key == token) {
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": bad input\n";
        throw Error::Input();
      }
      if(IO::log.is_open()) {
        std::cerr << funame << token << ": log stream allready open\n";
        throw Error::Input();
      }      
      IO::log.open(stemp.c_str());
      if(!IO::log) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // unknown keyword
    else {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  /*************************** CHECKING ******************************************/

  if(!multi_rotor) {
    std::cerr << funame << "multi-rotor has not been defined\n";
    throw Error::Init();
  }

  std::cout << "\nstates density:\n"
	    << std::setw(15) << "E, kcal/mol"
	    << std::setw(15) << "D, mol/kcal"
	    << "\n";

  for(double e = 0.5; e < 5;  e += 0.5) {
    std::cout << std::setw(15) << e
	      << std::setw(15) << multi_rotor->states(e * Phys_const::kcal) * Phys_const::kcal
	      << "\n";
  }

  return 0;
}
