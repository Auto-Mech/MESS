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

  Key  temp_key("TemperatureList[K]");
  Key   ang_key("Geometry[angstrom]");
  Key  bohr_key("Geometry[bohr]"    );
  Key  mrot_key("MultiRotor"        );
  Key   log_key("LogOutput"         );
  Key estep_key("EnergyGridStep[kcal/mol]");
  Key estart_key("EnergyGridStart[kcal/mol]");
  Key esize_key("EnergyGridSize");
  Key prop_key("Property");

  std::vector<double> temperature;
  std::vector<Atom> atom;
  SharedPointer<Model::MultiRotor>  rotor;

  if (argc < 2) {
    std::cout << "usage: multirotor_test input_file\n";
    return 0;
  }

  int    itemp;
  double dtemp;
  bool   btemp;

  double estep  =  0.;
  double estart =  0.;
  int    esize   = -1;

  IO::KeyBufferStream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string token, comment, line, stemp;
  int prop = Model::DENSITY;
  while(from >> token) {
    // geometry
    if(ang_key == token || bohr_key == token) {
      if(atom.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      if(ang_key == token)
	Model::read_geometry(from, atom, Model::ANGSTROM);
      else
	Model::read_geometry(from, atom, Model::BOHR);
    }
    // property
    else if(prop_key == token) {
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      if(stemp == "density")
	prop = Model::DENSITY;
      else if(stemp == "number")
	prop = Model::NUMBER;
      else {
	std::cerr << funame << token << ": unknown property: " << stemp << "\n";
	throw Error::Range();
      }
    }
    // multi-rotor
    else if(mrot_key == token) {
      if(rotor) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      if(!atom.size()) {
	std::cerr << funame << token << ": geometry not defined\n";
	throw Error::Init();
      }
      rotor = SharedPointer<Model::MultiRotor>(new Model::MultiRotor(from, atom, prop));
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
      IO::log.open(stemp.c_str());
      if(!IO::log) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // temperature
    else if(temp_key == token) {
      if(temperature.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      IO::LineInput temperature_input(from);
      while(temperature_input >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}      
	temperature.push_back( dtemp * Phys_const::kelv);
      }
    }
    // energy grid step
    else if(estep_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      if(dtemp == 0.) {
        std::cerr << funame << token << ": should not be zero\n";
        throw Error::Range();
      }      
      estep = dtemp;
    }
    // energy grid start
    else if(estart_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      if(dtemp <= 0.) {
        std::cerr << funame << token << ": should be positive\n";
        throw Error::Range();
      }      
      estart = dtemp;
    }
    // energy grid size
    else if(esize_key == token) {
      if(!(from >> esize)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      if(esize <= 0) {
        std::cerr << funame << token << ": should be positive\n";
        throw Error::Range();
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

  if(!rotor) {
    std::cerr << funame << "multirotor has not been defined\n";
    throw Error::Init();
  }


  if(esize > 0 && estep != 0.) {

    double stat_unit = 1.;
    if(rotor->mode() == Model::DENSITY) {
      IO::log << "density of states:\n"
	      << std::setw(15) << "E, kcal/mol"
	      << std::setw(15) << "D, cm"
	      << "\n";
      stat_unit = Phys_const::incm;
    }
    else if(rotor->mode() == Model::NUMBER)
      IO::log << "number of states:\n"
	      << std::setw(15) << "E, kcal/mol"
	      << std::setw(15) << "N"
	      << "\n";

    double ener = estart;
    for(int e = 0; e < esize;  ++e, ener += estep) {
      IO::log << std::setw(15) << ener
	      << std::setw(15) << rotor->states(ener * Phys_const::kcal) * stat_unit
	      << "\n";
    }
  }
  return 0;
}
