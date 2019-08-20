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

#include "random.hh"
#include "read.hh"
#include "crossrate.hh"
#include "potential.hh"
#include "key.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "crossrate: ";

  std::string stemp;
  double      dtemp;
  bool        btemp;
  int         itemp;

  if (argc < 2) {
    std::cerr << funame  << "usage: crossrate input_file\n";
    throw Error::Input();
  }

  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    throw Error::Input();
  }

  std::string base_name = argv[1];
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    base_name.resize(base_name.size() - 4);

  KeyGroup CrossRateDriver;

  Key   struc_key("Structure"             );
  Key    main_key("Main"                  );
  Key     pot_key("Potential"             );
  Key     var_key("Variables"             );
  Key    tran_key("TransitionState"       );
  Key exclude_key("ExcludeRegion"         );
  Key    reac_key("ReactiveRegion"        );
  Key     log_key("LogOutput"             );
  
  // input parameters
  std::vector<DivSur::MultiSur> tran_surf;
  Potential::Wrap pot;
  std::vector<CrossRate::ReactiveRegion> reactive_region;
  std::map<std::string, std::vector<double> > var_array;

  std::string token, comment;
  while(from >> token) {
    // log output
    if(log_key == token) {
      if(IO::log.is_open()) {
        std::cerr << funame << token << ": allready opened\n";
        throw Error::Init();
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
    // molecular structure
    else if(struc_key == token) {
      // default log output	
      if(!IO::log.is_open()) {
	stemp = base_name + ".log";
	IO::log.open(stemp.c_str());
	if(!IO::log) {
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	  throw Error::Input();
	}
      }

      std::getline(from, comment);
      Structure::init(from);
    }
    // crossrate initialization
    else if(main_key == token) {
      std::getline(from, comment);
      CrossRate::init(from);
    }
    // potential energy surface
    else if(pot_key == token) {
      from >> pot;
    }
    // variables for transition state surface initialization
    else if(var_key == token) {
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      for(int i = 0; i < itemp; ++i) {
	from >> stemp;
	IO::LineInput lin(from);
	while(lin >> dtemp)
	  var_array[stemp].push_back(dtemp);
      }
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      for(std::map<std::string, std::vector<double> >::const_iterator it = var_array.begin(); it != var_array.end(); ++it) 
	if(it != var_array.begin() && it-second.size() != var_array.begin()->second.size()) {
	  std::cerr << funame << token << ": wrong size\n";
	  throw Error::Range();
	}
    }
    // transition state (multiple species) surface
    else if(tran_key == token) {
      if(var_array.size() && var_array.begin()-second.size())
	tran_surf.resize(var_array.begin()-second.size());
      for(int s = 0; s < tran_surf.size(); ++s) {

	from >> tran_surf;
      }
    }
    // exclude region
    else if(exclude_key == token) {
      if(Dynamic::exclude_region) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      std::getline(from, comment);

      Dynamic::exclude_region = Dynamic::CCP(new DivSur::BaseSur(from));
    }
    // reactive region
    else if(reac_key == token) {
      std::getline(from, comment);
      reactive_region.push_back(CrossRate::ReactiveRegion(from));
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

  // initialization check
  pot.isinit();
  if(!tran_surf.isinit()) {
    std::cerr << funame << "transition state surface not initialized\n";
    throw Error::Init();
  }
  
  if(!reactive_region.size() && CrossRate::reactive_energy() < -10.) {
    std::cerr << funame << "either reactive region or reactive energy should be active\n";
    throw Error::Init();
  }

  if(!Dynamic::exclude_region)
    std::cerr << funame << "WARNING: exclude region not initialized\n";
    
  // Random numbers initialization
  if(CrossRate::seed())
    Random::init(CrossRate::seed());
  else
    Random::init();


  IO::log << "\n" << IO::log_offset;
  switch(CrossRate::job()) {
  case CrossRate::DYN_JOB:
    IO::log << "CLASSICAL TRAJECTORIES ";
    break;
  case CrossRate::STAT_JOB:
    IO::log << "TRANSITION STATE THEORY ";
    break;
  case CrossRate::TEST_JOB:
    IO::log << "TEST ";
    break;
  }

  switch(CrossRate::mode()) {
  case CrossRate::T_MODE:
    IO::log << "CANONICAL CALCULATION\n";
    break;
  case CrossRate::E_MODE:
    IO::log << "MICROCANONICAL CALCULATION\n";
    break;
  case CrossRate::J_MODE:
    IO::log << "E,J-RESOLVED CALCULATION\n";
    break;
  }
  IO::log << "\n";

  //IO::log << IO::log_offset << "Transition state:\n";
  //tran_surf.print(IO::log, std::string(IO::log_offset + 3,' '));
  //IO::log << "\n";

  if(CrossRate::job() == CrossRate::TEST_JOB)
    return 0;

  // stop condition
  Dynamic::CCP stop(new Dynamic::DistanceCondition(CrossRate::diss_dist()));

  if(CrossRate::reactive_energy() > -10.)
    stop |= Dynamic::CCP(new Potential::Condition(pot, CrossRate::reactive_energy()));

  for(int r = 0; r < reactive_region.size(); ++r)
    stop |= reactive_region[r].condition(pot);
  
  //CrossRate::MultiArray cross_rate(tran_surf, pot, stop);

  return 0;
}
