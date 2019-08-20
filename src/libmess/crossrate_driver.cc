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
#include <cstring>

#include "random.hh"
#include "read.hh"
#include "crossrate.hh"
#include "potential.hh"
#include "key.hh"

Potential::Wrap pot;

Dynamic::CCP reactive_region (std::istream& from)
{
  const char funame [] = "reactive_region: ";

  IO::Marker funame_marker(funame);
  
  std::string stemp;
  double      dtemp;
  bool        btemp;
  int         itemp;

  KeyGroup ReactiveRegionGroup;

  Key  prim_key("ZmatPrimitive"           );
  Key space_key("SpaceRegion"             );
  Key  ener_key("ReactiveEnergy[kcal/mol]");  


  Dynamic::CCP res;

  bool isener   = false;
  bool isregion = false;

  std::string token, comment;
  while(from >> token) {
    if(IO::end_key() == token) {
      break;
    }
    // primitive (distance[angstrom], angle[degrees], or dihedral[degrees])
    else if(prim_key == token) {
      IO::LineInput lin(from);
      // limits
      std::string limit_str[2];
      for(int i = 0; i < 2; ++i)
	lin >> limit_str[i];

      char* endptr;
      
      std::map<int, double> limit;
      for(int i = 0; i < 2; ++i) {
	dtemp = std::strtod(limit_str[i].c_str(), &endptr);
	if(!std::strlen(endptr))
	  limit[i] = dtemp;
      }
      
      if(!limit.size()) {
	std::cerr << funame << token << ": at least one limit should be defined\n";
	throw Error::Range();
      }

      if(limit.size() == 2 && limit[0] >= limit[1]) {
	std::cerr << funame << token << ": lower limit should be smaller than the upper one\n";
	throw Error::Range();
      }
      
      std::vector<int> atom;
      while(lin >> itemp)
	atom.push_back(itemp);

      switch(atom.size()) {
      case 2: // distance
	for(std::map<int, double>::const_iterator it = limit.begin(); it != limit.end(); ++it)
	  if(it->second <= 0.) {
	    std::cerr << funame << token << ": distance limit out of range: " << it->second << "\n";
	    throw Error::Range();
	  }

	IO::log << IO::log_offset << "distance";
	for(int a = 0; a < atom.size(); ++a)
	  IO::log << " " <<  atom[a] / 2 + 1 << "(" << atom[a] % 2 + 1 << ")";
	IO::log << " [angstrom] between ";
	break;
	
      case 3: // angle
	for(std::map<int, double>::const_iterator it = limit.begin(); it != limit.end(); ++it)
	  if(it->second <= 0. || it->second >= 180.) {
	    std::cerr << funame << token << ": angle limit out of range: " << it->second << "\n";
	    throw Error::Range();
	  }

	IO::log << IO::log_offset << "angle";
	for(int a = 0; a < atom.size(); ++a)
	  IO::log << " " <<  atom[a] / 2 + 1 << "(" << atom[a] % 2 + 1 << ")";
	IO::log << " [degree] between ";
	break;
	
      case 4: // dihedral
	IO::log << IO::log_offset << "dihedral angle";
	for(int a = 0; a < atom.size(); ++a)
	  IO::log << " " <<  atom[a] / 2 + 1 << "(" << atom[a] % 2 + 1 << ")";
	IO::log << " [degree] between ";
	break;
	
      default:
	  std::cerr << funame << token << ": number of atoms out of range: " << atom.size() << "\n";
	throw Error::Range();
      }

      for(int l = 0; l < 2; ++l) {
	if(limit.find(l) != limit.end())
	  IO::log << limit[l];
	else
	  IO::log << "***";
	
	if(!l)
	  IO::log << " and ";
      }
      IO::log << "\n";

      // scale limits
      switch(atom.size()) {
      case 2: // distance
	
	for(std::map<int, double>::iterator it = limit.begin(); it != limit.end(); ++it)
	  it->second *= Phys_const::angstrom;
	break;
	
      case 3: // angle

	for(std::map<int, double>::iterator it = limit.begin(); it != limit.end(); ++it)
	  it->second *= M_PI / 180.;
	break;
	
      case 4: // dihedral
	
	for(std::map<int, double>::iterator it = limit.begin(); it != limit.end(); ++it)
	  it->second *= M_PI / 180.;
	break;
      }
      
      ConstSharedPointer<Dynamic::ObservBase> observ(new Dynamic::ZmatObserv(atom));
      
      Dynamic::CCP prim;
      if(limit.size() == 2)
	prim.init(new Dynamic::ObservCondition(observ, limit[0], limit[1]));
      else if(limit.begin()->first == 0)
	prim.init(new Dynamic::ObservCondition(observ, limit[0], Dynamic::ObservCondition::LOWER_BOUND));
      else
	prim.init(new Dynamic::ObservCondition(observ, limit[1], Dynamic::ObservCondition::UPPER_BOUND));

      if(res)
	res &= prim;
      else
	res  = prim;
    }
    // space region
    else if(space_key == token) {
      if(isregion) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      isregion = true;

      const DivSur::OneSur* reg = new DivSur::OneSur(from);

      IO::log << IO::log_offset << "SpaceRegion:\n";
      reg->print(IO::log, std::string(IO::log_offset + 3, ' '));

      if(res)
	res &= Dynamic::CCP(reg);
      else
	res  = Dynamic::CCP(reg);
    }
    // reactive energy
    else if(ener_key == token) {
      if(isener) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      isener = true;

      IO::LineInput lin(from);
      if(!(lin >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      
      IO::log << IO::log_offset << "reactive energy [kcal/mol] = " << dtemp << "\n";
      
      dtemp *= Phys_const::kcal;

      if(res)
	res &= Dynamic::CCP(new Potential::Condition(pot, dtemp));
      else
	res  = Dynamic::CCP(new Potential::Condition(pot, dtemp));
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Input();
  }

  if(!res) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
  
  return res;
}

Dynamic::CCP stop_condition (std::istream& from)
{
  const char funame [] = "stop_condition: ";

  IO::Marker funame_marker(funame);
  
  std::string stemp;
  double      dtemp;
  bool        btemp;
  int         itemp;

  KeyGroup StopConditionGroup;

  Key dist_key("DissociationDistance[angstrom]");
  Key ener_key("ReactiveEnergy[kcal/mol]"      );
  Key reac_key("ReactiveRegion"                );

  Dynamic::CCP res;
  
  double diss_dist = -1.;
  
  std::string token, comment;
  while(from >> token) {
    if(IO::end_key() == token) {
      break;
    }
    // separation distance
    else if(dist_key == token) {
      if(diss_dist > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      if(!(from >> diss_dist)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(diss_dist <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      IO::log << IO::log_offset << "dissociation distance [angstrom] = " << diss_dist << "\n";

      diss_dist *= Phys_const::angstrom;
      
      if(res)
	res |= Dynamic::CCP(new Dynamic::DistanceCondition(diss_dist));
      else
	res  = Dynamic::CCP(new Dynamic::DistanceCondition(diss_dist));
    }
    // reactive energy
    else if(ener_key == token) {
      if(CrossRate::reactive_energy() < 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(dtemp >= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      IO::log << IO::log_offset << "reactive energy [kcal/mol] = " << dtemp << "\n";
      
      CrossRate::set_reactive_energy(dtemp * Phys_const::kcal);
      
      if(res)
	res |= Dynamic::CCP(new Potential::Condition(pot, CrossRate::reactive_energy()));
      else
	res  = Dynamic::CCP(new Potential::Condition(pot, CrossRate::reactive_energy()));
	
    }
    // reactive region
    else if(reac_key == token) {
      std::getline(from, comment);
      
      if(res)
	res |= reactive_region(from);
      else
	res  = reactive_region(from);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Input();
  }
  
  if(CrossRate::reactive_energy() > 0.) {
    std::cerr << funame << "reactive energy not initialized\n";
    throw Error::Init();
  }

  if(diss_dist < 0.) {
    std::cerr << funame << "dissociation distance not initialized\n";
    throw Error::Init();
  }

  return res;
}

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

  Key    main_key("Main"                  );
  Key   struc_key("Structure"             );
  Key     pot_key("Potential"             );
  Key    tran_key("TransitionState"       );
  Key exclude_key("ExcludeRegion"         );
  Key    stop_key("StopCondition"         );
  Key     log_key("LogOutput"             );
  
  DivSur::MultiSur tran;
  Dynamic::CCP     stop;
  
  // input
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
    // transition state (multiple species) surface
    else if(tran_key == token) {
      from >> tran;
    }
    // exclude region
    else if(exclude_key == token) {
      if(Dynamic::exclude_region) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      std::getline(from, comment);

      Dynamic::exclude_region = Dynamic::CCP(new DivSur::OneSur(from));
    }
    // stop condition
    else if(stop_key == token) {
      if(!pot) {
	std::cerr << funame << token << ": potential should be initialized first\n";
	throw Error::Init();
      }

      if(stop) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      std::getline(from, comment);
      stop = stop_condition(from);
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
  if(!tran.isinit()) {
    std::cerr << funame << "transition state surface not initialized\n";
    throw Error::Init();
  }
  
  if(!stop) {
    std::cerr << funame << "stop condition not initialized\n";
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

  IO::log << IO::log_offset << "Transition state:\n";
  tran.print(IO::log, std::string(IO::log_offset + 3,' '));
  IO::log << "\n";

  if(CrossRate::job() == CrossRate::TEST_JOB)
    return 0;

  CrossRate::MultiArray cross_rate(tran, pot, stop);

  return 0;
}
