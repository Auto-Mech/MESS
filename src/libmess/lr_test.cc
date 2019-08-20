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

#include "lr.hh"
#include "read.hh"

#include <cmath>

int main (int argc, char* argv [])
{
  const char funame [] = "lr: ";

  int itemp;
  bool btemp;
  double dtemp;

  if (argc < 2) {
    std::cout << funame  << "usage: lr input_file\n";
    return 1;
  }

  std::map<std::string, Read> input;
  std::map<std::string, Read>::iterator idit;
 
  double energy, distance;

  input ["Structure"             ] = Read(Structure::init);
  input ["Main"                  ] = Read(LongRange::init);
  input ["Potential"             ] = Read(LongRange::pot);
  input ["Energy[kcal/mol]"      ] = Read(energy);
  input ["Distance[au]"          ] = Read(distance);

  // actual input
  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string key, comment;
  while(from >> key) {
    idit = input.find(key);
    if (idit == input.end()) {
      std::cerr << funame << "WARNING: did not find the key: " << key
		<< ", taking it as a comment\n";
      getline(from, comment);
    }
    else
      from >> idit->second;
  }
  from.close();
  from.clear();

  // check if all parameters were initialized
  btemp = true;
  for(idit = input.begin(); idit != input.end(); ++idit)
    if(!idit->second.is_init()) {
      std::cerr << funame << idit->first << " is not initialized\n";
      btemp = false;
    }
  if(!btemp)
    return 1;

  // checking for default values
  std::cout << funame << "Default parameters:\n" << std::left;
  for(idit = input.begin(); idit != input.end(); ++idit)
    if(idit->second.is_default())
      std::cout << "   " << std::setw(20) << idit->first << " = " << idit->second << "\n";
  std::cout << "\n" << std::right;

  energy *= Phys_const::kcal;

  LongRange::ENumberFix N_e(distance, energy);
  std::cout << "e-resolved number of states = "  << N_e() << "\n" << std::endl; 

  LongRange::ENumberOpt NO_e(distance, energy);
  std::cout << "e-optimized number of states = "  << NO_e() << "\n" << std::endl; 

  LongRange::JNumberOpt N_j(distance);
  LongRange::JIntegral  N_ej(&N_j, energy);
  std::cout << "\ne,j-optimized number of states = "  << N_ej() << "\n" << std::endl; 
  
  LongRange::MNumberOpt N_m(distance);
  LongRange::MIntegral  N_jm(&N_m);
  LongRange::JIntegral  N_ejm(&N_jm, energy);
  std::cout << "\ne,j,m-optimized number of states = "  << N_ejm() << "\n" << std::endl; 
  
  LongRange::KNumberOpt N_k(distance);
  LongRange::KIntegral  N_mk(&N_k);
  LongRange::MIntegral  N_jmk(&N_mk);
  LongRange::JIntegral  N_ejmk(&N_jmk, energy);
  std::cout << "\ne,j,m,k-optimized number of states = "  << N_ejmk() << std::endl; 

  return 0;
}
