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

#include "molpro.hh"
#include "system.hh"
#include "error.hh"
#include "units.hh"
#include "key.hh"

#include <dlfcn.h>
#include <cstdlib>
#include <cctype>
#include <unistd.h>
#include <cstdio>
#include <cerrno>
#include <cstring>

#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>

// available potential corrections
namespace Molpro {

  bool _isinit = false;

  std::string input_template;
  std::string executable;
  std::string base_name;
  std::string scratch_dir;

  int geom_start = -1;
  bool wfu = false;

  std::string ener_pattern;
  std::string fail_pattern;

  int fail_counter = 0;

  void backup (const std::string&) ;
}

bool Molpro::isinit () 
{ 
  return _isinit; 
}

void Molpro::backup (const std::string& mess) 
{
  const char funame [] = "Molpro::backup: ";

  IO::log << "Molpro::pot: " << mess << "\n";

  std::ostringstream backup_file;
  backup_file << base_name << ".back." << fail_counter++;

  if(std::rename((base_name + ".out").c_str(), backup_file.str().c_str()))
    IO::log << funame << "failed" << std::endl;
  else 
    IO::log << funame <<  "check " << backup_file.str() << " file" << std::endl;

  throw Error::Molpro(mess);
}

void Molpro::set_scratch_dir (const std::string& s) 
{
  scratch_dir = s;
}

void Molpro::remove_wfu ()
{
  const char funame [] = "Molpro::remove_wfu: ";

  if(!scratch_dir.size()) {
    std::cerr << funame << "scratch directory is not initialized\n";
    throw Error::Init();
  }
  
  if(!wfu) {
    std::cerr << funame << "wfu not in template\n";
    throw Error::Init();
  }

  std::string stemp = scratch_dir + "/" + base_name + ".wfu";
  std::remove(stemp.c_str());
}

void Molpro::init(std::istream& from) 
{
  const char funame [] = "Molpro::init: ";

  if(_isinit) {
    std::cout << funame << "molpro environment is allready initialized\n";
    throw Error::Init();
  }
  _isinit = true;

  std::string stemp;
  int         itemp;
  double      dtemp;

  KeyGroup MolproInit;

  Key temp_key("TemplateFile"   );
  Key base_key("BaseFileName"   );
  Key exec_key("MolproExec"     );
  Key geom_key("GeometryPattern");
  Key ener_key("EnergyPattern"  );
  Key fail_key("FailurePattern" );

  std::string token, comment;
  while(from >> token) {
    // end input 
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // input template
    else if(temp_key == token) {
      if(input_template.size()) {
	std::cerr << funame <<  token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(!base_name.size())
	base_name = stemp;

      std::getline(from, comment);

      std::ifstream temp_in(stemp.c_str());
      if(!temp_in) {
	std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	throw Error::Input();
      }

      while(std::getline(temp_in, stemp))
	input_template += stemp + "\n";
    }
    // molpro executable
    else if(exec_key == token) {
      if(executable.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> executable)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      if(access(executable.c_str(), X_OK)) {
	std::cerr << funame << token << ": no access to " << executable << "\n";
	throw Error::Init();
      }
    }
    // base file name
    else if(base_key == token) {
      if(!(from >> base_name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);
    }
    // geometry pattern
    else if(geom_key == token) {
      if(geom_start >= 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      if(!input_template.size()) {
	std::cerr << funame << token << ": input template not initialized yet\n";
	throw Error::Init();
      }

      geom_start = input_template.find(stemp);
      if(geom_start < 0 || geom_start >= input_template.size()) {
	std::cerr << funame << token << ": cannot find " << stemp << " in the input file template\n";
	throw Error::Input();
      }
      input_template.erase(geom_start, stemp.size());
    }
    // energy pattern
    else if(ener_key == token) {
      if(ener_pattern.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> ener_pattern)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      if(input_template.find(ener_pattern) >= input_template.size()) {
	std::cerr << funame << token << " : cannot find " << ener_pattern
		  << " pattern in the input file template\n";
	throw Error::Input();
      }

      // convert to upper case
      for(int i = 0; i < ener_pattern.size(); ++i)
	ener_pattern[i] = std::toupper(ener_pattern[i]);

    }
    // fail pattern
    else if(fail_key == token) {
      if(fail_pattern.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> fail_pattern)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      if(input_template.find(fail_pattern) >= input_template.size()) {
	std::cerr << funame << token << " : cannot find " << fail_pattern
		  << " pattern in the input template\n";
	throw Error::Input();
      }

      // convert to upper case
      for(int i = 0; i < fail_pattern.size(); ++i)
	fail_pattern[i] = std::toupper(fail_pattern[i]);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  if(!from) {
    std::cerr << funame << " input corrupted\n";
    throw Error::Input();
  }

  if(!input_template.size()) {
    std::cerr << funame << "input template not initialized\n";
    throw Error::Init();
  }

  if(!ener_pattern.size()) {
    std::cerr << funame << "energy pattern not initialized\n";
    throw Error::Init();
  }

  if(geom_start < 0) {
    std::cerr << funame << "geometry pattern not initialized\n";
    throw Error::Init();
  }

  // check if wfu is defined
  stemp = base_name + ".wfu";
  if(input_template.find(stemp) < input_template.size())
    wfu = true;

} // Molpro::init

void Molpro::pot (const std::vector<Atom>& molec, Array<double>& ener_data, int flags) 
{
  const char funame [] = "Molpro::pot: ";

  double      dtemp;
  int         itemp;
  std::string stemp;

  // create a geometry string
  std::ostringstream geom;
  for(std::vector<Atom>::const_iterator at = molec.begin(); at != molec.end(); ++at) {
    if(at != molec.begin())
      geom << "\n";
    geom << *at;
  }

  // create input string
  std::string input = input_template;
  input.insert(geom_start, geom.str());

  // create an input file for molpro
  const std::string inp_name = base_name + ".inp";
  std::ofstream to(inp_name.c_str());
  if(!to) {
    IO::log << funame << "cannot open " << inp_name << " file\n";
    throw Error::Open();
  }
  to << input;
  to.close();

  // run molpro executable
  const std::string out_name = base_name + ".out";
  std::remove(out_name.c_str());

  if(scratch_dir.size())
    itemp = System::call_exe(executable.c_str(), "--scratch", scratch_dir.c_str(), inp_name.c_str(), (char*) 0);
  else
    itemp = System::call_exe(executable.c_str(), inp_name.c_str(), (char*) 0);

  if(itemp)
    backup("run failure");

  // read molpro results
  std::ifstream from(out_name.c_str());
  if(!from) {
    IO::log << funame << "cannot open " << out_name << "\n";
    throw Error::Molpro("cannot open output file");
  }

  // scan molpro output
  int ener_count = 0;
  while(std::getline(from, stemp)) {
    std::istringstream lin(stemp);
    
    // variable assignment
    if(lin >> stemp && stemp == "SETTING") {

      if(!(lin >> stemp))
	  backup("cannot read output");

      // read energy
      if(ener_pattern == stemp) {

	if(!(lin >> stemp) || stemp != "=" || !(lin >> dtemp))
	  backup("cannot read energy");
	
	IO::log << funame << "molpro energy [kcal/mol] = " << dtemp / Phys_const::kcal << std::endl;

	if(ener_count == ener_data.size())
	  backup("too many energies");

	ener_data[ener_count++] = dtemp;
      }// read energy

      // convergence failure
      if(fail_pattern == stemp)
	backup("convergence failure");

    }// variables assignment
  }// scan molpro output

  if(ener_count != ener_data.size())
    backup("too few energies");
}
