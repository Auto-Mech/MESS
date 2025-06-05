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

namespace Molpro {
  //
  int mute = 0;
  
  bool _isinit = false;
  bool isinit () { return _isinit; }

  std::string scratch_dir;
  void set_scratch_dir (const std::string& s) { scratch_dir = s; }
  
  std::string executable;
  std::string base_name;
  
  std::string main_template;
  int main_pos = -1;

  std::string guess_template;
  int guess_pos = -1;

  std::string relax_template;
  int relax_pos = -1;
  
  std::string ener_pattern;
  std::string fail_pattern;

  bool is_guess () { return guess_template.size(); }

  // backup facility and fail counter
  //
  int fail_counter = 0;
  
  void backup (const std::string&, int =0);

  void read (Array<double>&, int =0);
  
  std::string make_geom (const std::vector<Atom>&, int =0);
}

// backup
//
void Molpro::backup (const std::string& mess, int flags) 
{
  const char funame [] = "Molpro::backup: ";

  IO::log << mess;

  std::ostringstream backup_file;

  backup_file << base_name << ".back." << fail_counter++;

  std::string out_file = base_name + ".out";
  
  if(std::rename(out_file.c_str(), backup_file.str().c_str())) {
    //
    IO::log << funame << "failed" << std::endl;
  }
  else
    //
    IO::log << funame <<  "check " << backup_file.str() << " file" << std::endl;

  throw Error::Molpro(mess);
}

void Molpro::remove_wfu ()
{
  const char funame [] = "Molpro::remove_wfu: ";

  if(!scratch_dir.size()) {
    //
    IO::log << funame << "scratch directory is not initialized\n";
    throw Error::Init();
  }
  
  std::string stemp = scratch_dir + "/" + base_name + ".wfu";
  std::remove(stemp.c_str());
}

void Molpro::init(std::istream& from) 
{
  const char funame [] = "Molpro::init: ";

  if(_isinit) {
    //
    if(!mute)
      //
      std::cerr << funame << "molpro environment is allready initialized\n";
    
    throw Error::Init();
  }
  _isinit = true;

  std::string stemp;
  int         itemp;
  double      dtemp;

  KeyGroup MolproInit;

  Key  main_key("MainTemplate"   );
  Key guess_key("GuessTemplate"  );
  Key relax_key("RelaxTemplate"  );
  Key  exec_key("MolproExec"     );
  Key  geom_key("GeometryPattern");
  Key  ener_key("EnergyPattern"  );
  Key  fail_key("FailurePattern" );

  std::string token, comment, line, geom_pattern;
  
  while(from >> token) {
    //
    // end input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // main input template
    //
    else if(main_key == token) {
      //
      if(main_template.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame <<  token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> base_name)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      std::ifstream fin(base_name.c_str());
      
      if(!fin) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": cannot open " << base_name << " file\n";
	
	throw Error::Input();
      }

      while(std::getline(fin, line))
	//
	main_template += line + "\n";
    }
    // geometry relaxation input template
    //
    else if(relax_key == token) {
      //
      if(relax_template.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame <<  token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      std::ifstream fin(stemp.c_str());
      
      if(!fin) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	
	throw Error::Input();
      }

      while(std::getline(fin, line))
	//
	relax_template += line + "\n";
    }
    // initial guess input template
    //
    else if(guess_key == token) {
      //
      if(guess_template.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame <<  token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      std::ifstream fin(stemp.c_str());
      
      if(!fin) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	
	throw Error::Input();
      }

      while(std::getline(fin, line))
	//
	guess_template += line + "\n";
    }
    // molpro executable
    //
    else if(exec_key == token) {
      //
      if(executable.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> executable)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      if(access(executable.c_str(), X_OK)) {
	//
	if(!mute)
	  //
	std::cerr << funame << token << ": no access to " << executable << "\n";
	
	throw Error::Init();
      }
    }
    // geometry pattern
    //
    else if(geom_key == token) {
      //
      if(geom_pattern.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> geom_pattern)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
    }
    // energy pattern
    //
    else if(ener_key == token) {
      //
      if(ener_pattern.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> ener_pattern)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      if(main_template.find(ener_pattern) >= main_template.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << " : cannot find " << ener_pattern << " pattern in the input file template\n";
	
	throw Error::Input();
      }

      // convert to upper case
      //
      for(int i = 0; i < ener_pattern.size(); ++i)
	//
	ener_pattern[i] = std::toupper(ener_pattern[i]);
    }
    // fail pattern
    //
    else if(fail_key == token) {
      //
      if(fail_pattern.size()) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> fail_pattern)) {
	//
	if(!mute)
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      // convert to upper case
      //
      for(int i = 0; i < fail_pattern.size(); ++i)
	//
	fail_pattern[i] = std::toupper(fail_pattern[i]);
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      if(!mute) {
	//
	std::cerr << funame << "unknown keyword: " << token << "\n";
      
	Key::show_all(std::cerr);
	std::cerr << "\n";
      }
      throw Error::Init();
    }
  }

  if(!from) {
    //
    if(!mute)
      //
      std::cerr << funame << " input corrupted\n";
    
    throw Error::Input();
  }
  if(!main_template.size()) {
    //
    if(!mute)
      //
      std::cerr << funame << "main input template not initialized\n";
    
    throw Error::Init();
  }
  if(!ener_pattern.size()) {
    //
    if(!mute)
      //
      std::cerr << funame << "energy pattern not initialized\n";
    
    throw Error::Init();
  }

  // geometry input point in the main input template
  //
  main_pos = main_template.find(geom_pattern);
  
  if(main_pos < 0 || main_pos >= main_template.size()) {
    //
    if(!mute)
      //
      std::cerr << funame << token << ": cannot find " << geom_pattern << " in the main input template\n";
    
    throw Error::Input();
  }
  main_template.erase(main_pos, geom_pattern.size());
  
  // geometry input point in the initial guess input template
  //
  if(guess_template.size()) {
    //
    guess_pos = guess_template.find(geom_pattern);
  
    if(guess_pos < 0 || guess_pos >= guess_template.size()) {
      //
      if(!mute)
	//
	std::cerr << funame << token << ": cannot find " << geom_pattern << " in the initial guess input template\n";
      
      throw Error::Input();
    }
    guess_template.erase(guess_pos, geom_pattern.size());
  }
  //
  // geometry input point in the geometry relaxation input template
  //
  if(relax_template.size()) {
    //
    relax_pos = relax_template.find(geom_pattern);
  
    if(relax_pos < 0 || relax_pos >= relax_template.size()) {
      //
      if(!mute)
	//
	std::cerr << funame << token << ": cannot find " << geom_pattern << " in the geometry relaxation input template\n";
      
      throw Error::Input();
    }
    relax_template.erase(guess_pos, geom_pattern.size());
  }
}// Molpro::init

// make geometry section in molpro input file
//
std::string Molpro::make_geom (const std::vector<Atom>& molec, int flags)
{
  const char funame [] = "Molpro::make_geom: ";

  if(!molec.size()) {
    //
    IO::log << funame << "no atoms\n";

    throw Error::Init();
  }
  
  std::ostringstream geom;

  if(flags & WFU)
    //
    geom << "file, 2, " << base_name << ".wfu\n\n";
  
  geom << "symmetry, nosym\norient, noorient\ngeomtyp=xyz\n\ngeometry={\n"
    //
       << molec.size() << "\ncrossrate\n";
  
  for(std::vector<Atom>::const_iterator at = molec.begin(); at != molec.end(); ++at) {
    //
    if(at != molec.begin())
      //
      geom << "\n";
    
    geom << *at;
  }
  geom << "\n}";

  return geom.str();
}

// potential calculation
//
void Molpro::pot (const std::vector<Atom>& molec, Array<double>& ener, int flags) 
{
  const char funame [] = "Molpro::pot: ";

  double      dtemp;
  int         itemp;
  std::string stemp;

  std::ostringstream oss;
  
  // create a geometry string
  //
  std::string geom = make_geom(molec, flags);
  
  // create input string
  //
  std::string input;
  
  if(flags & GUESS) {
    //
    if(!guess_template.size()) {
      //
      oss.str("");
      
      oss << funame << "guess template not initialized\n";

      IO::log << oss.str();
      
      throw Error::Molpro(oss.str());
    }
    
    input = guess_template;
    
    input.insert(guess_pos, geom);
  }
  else if(flags & RELAX) {
    //
    if(!relax_template.size()) {
      //
      oss.str("");
      
      oss << funame << "geometry relaxation template not initialized\n";

      IO::log << oss.str();
      
      throw Error::Molpro(oss.str());
    }
    
    input = relax_template;
    
    input.insert(relax_pos, geom);
  }
  else {
    //
    input = main_template;
    
    input.insert(main_pos, geom);
  }

  // create an input file for molpro
  //
  const std::string inp_name = base_name + ".inp";
  
  std::ofstream to(inp_name.c_str());
  
  if(!to) {
    //
    oss.str("");
    
    oss << funame << "cannot open " << inp_name << " file for input\n";

    IO::log << oss.str();
    
    throw Error::Molpro(oss.str());
  }
  to << input;
  to.close();

  // run molpro executable
  //
  std::remove((base_name + ".out").c_str());

  std::remove((base_name + ".log").c_str());

  if(scratch_dir.size()) {
    //
    itemp = System::call_exe(executable.c_str(), "--scratch", scratch_dir.c_str(), inp_name.c_str(), (char*) 0);
  }
  else
    //
    itemp = System::call_exe(executable.c_str(), inp_name.c_str(), (char*) 0);

  if(itemp) {
    //
    oss.str("");

    oss << funame << "run failure\n";

    backup(oss.str());
  }

  if(!(flags & GUESS))
    //
    read(ener, flags);
}

// read molpro energies
//
void Molpro::read(Array<double>& ener, int flags)
{
  const char funame [] = "Molpro::read: ";

  double dtemp;
  int    itemp;

  std::ostringstream oss;
  
  std::string restart = " Reading variables from file 2";
  
  std::string stemp, token, line, name;
  
  const std::string out_name = base_name + ".out";
  
  std::ifstream from(out_name.c_str());
  
  if(!from) {
    //
    oss.str("");
    
    oss << funame << "cannot open " << out_name << " file\n";

    IO::log << oss.str();
    
    throw Error::Molpro(oss.str());
  }

  // scan molpro output
  //
  int count = 0;
  
  while(std::getline(from, line)) {
    //
    // removing variables read from saved wave function file
    //
    if(line == restart) {
      //
      std::getline(from, line); // empty line

      while(1) {
	//
	std::getline(from, line);

	if(line == "")
	  //
	  break;
      }
      continue;
    }
    
    std::istringstream lin(line);
    
    // variable assignment
    //
    if(!(lin >> token))
      //
      continue;
    
    if(token == ener_pattern) {

      if(!(lin >> stemp)) {
	//
	oss.str("");

	oss << funame << "cannot read the rest of the line: " << line << "\n";
	
	backup(oss.str());
      }
      if(stemp != "=") {
	//
	oss.str("");

	oss << funame << "wrong token: "<< stemp << ": shoud be =: full line: " << line << "\n";
	
	backup(oss.str());
      }
      if(!(lin >> dtemp)) {
	//
	oss.str("");

	oss << funame << "cannot read energy value: full line: " << line << "\n";
	
	backup(oss.str());
      }
      //IO::log << funame << "molpro energy [au] = " << dtemp << std::endl;
      
      if(count == ener.size()) {
	//
	oss.str("");

	oss << funame << "# of reported energies in molpro output file is more than requested: " << ener.size() << "\n";
	
	backup(oss.str());
      }

      ener[count++] = dtemp;
      //
    }
    // molpro failure
    //
    else if(token == fail_pattern) {
      //
      oss.str("");
      
      oss << funame << "failure pattern encountered: " << fail_pattern << "\n";
      
      backup(oss.str());
      //
    }
    //
  }// scan molpro output

  if(count != ener.size()) {
    //
    oss.str("");

    oss << funame << "# of reported energies in molpro output file, " << count << ", is less than requested, " << ener.size() << "\n";
    
    backup(oss.str());
  }
}

Lapack::Vector Molpro::gradients (int atom_size)
{
  const char funame [] = "Molpro::gradients: ";

  double      dtemp;
  int         itemp;
  bool        btemp;
  std::string stemp;

  const int size = atom_size * 3;
  
  std::ostringstream oss;
  
  std::string token, comment, name, line;

  // hassian string from molpro file
  //
  std::string grad_pattern = " Atom          dE/dx               dE/dy               dE/dz";

  // read molpro results
  //
  std::string out_name = base_name + ".out";
  
  std::ifstream from(out_name.c_str());
  
  if(!from) {
    //
    oss.str("");
    
    oss << funame << "cannot open " << out_name << " file\n";

    IO::log << oss.str();
    
    throw Error::Molpro(oss.str());
  }

  while(from) {
    //
    std::getline(from, line);

    if(line == grad_pattern)
      //
      break;
  }

  if(!from) {
    //
    oss.str("");
    
    oss << funame << "cannot find gradients seach pattern: " << grad_pattern << "\n";

    backup(oss.str());
  }

  Lapack::Vector res(size);

  std::getline(from, comment);
    
  for(int a = 0; a < atom_size; ++a) {
    //
    IO::LineInput lin(from);

    lin >> itemp;

    for(int i = 0; i < 3; ++i)
      //
      lin >> res[i + 3 * a];

    if(!lin){
      //
      oss.str("");
	
      oss << funame << a << "-th atom: reading failure\n";
      
      backup(oss.str());
    }
  }

  return res;
}

// hessian
//
Lapack::SymmetricMatrix Molpro::hessian (int atom_size)
{
  const char funame [] = "Molpro::hessian: ";

  double      dtemp;
  int         itemp;
  bool        btemp;
  std::string stemp;

  const int size = 3 * atom_size;

  std::ostringstream oss;
  
  std::string token, comment, name, line;

  // hassian string from molpro file
  //
  std::string hessian_pattern = " Mass weighted Second Derivative Matrix";

  // read molpro results
  //
  std::string out_name = base_name + ".out";
  
  std::ifstream from(out_name.c_str());
  
  if(!from) {
    //
    oss.str("");
    
    oss << funame << "cannot open " << out_name << " file\n";

    IO::log << oss.str();
    
    throw Error::Molpro(oss.str());
  }

  while(from) {
    //
    std::getline(from, line);

    if(line == hessian_pattern)
      //
      break;
  }

  if(!from) {
    //
    oss.str("");
    
    oss << funame << "cannot find hessian seach pattern: " << hessian_pattern << "\n";

    backup(oss.str());
  }

  Lapack::SymmetricMatrix res(size), hit(size);

  res = 0.;

  hit = 0.;
  
  int column = 0;

  while(column < size) {
    //
    std::getline(from, comment);

    for(int row = column; row < size; ++row) {
      //
      IO::LineInput lin(from);

      lin >> stemp;

      itemp = column;

      while(lin >> dtemp) {
	//
	if(itemp > row) {
	  //
	  oss.str("");
	  
	  oss << funame << "column index, " << itemp << ", went beyond row one, " << row << ": full line: " << lin.str() << "\n";
	  
	  backup(oss.str());
	}
	
	res(row, itemp) = dtemp;
	hit(row, itemp) += 1.;

	++itemp;
      }
    }
    column += 5;
  }

  // testing
  //
  itemp = 0;
  
  for(int i = 0; i < size; ++i)
    //
    for(int j = 0; j <=i; ++j)
      //
      if((int)std::floor(hit(i, j) + 0.5) != 1)
	//
	itemp = 1;

  if(itemp) {
    //
    oss.str("");
    
    oss << funame << "reading failed: reading hits matrix:\n";
    
    for(int i = 0; i < size; ++i) {
      //
      for(int j = 0; j <=i; ++j)
	//
	oss << std::setw(2) << (int)std::floor(hit(i, j) + 0.5);

      oss << "\n";
    }

    backup(oss.str());
  }
  
  res /= Phys_const::amu;

  return res;
}
