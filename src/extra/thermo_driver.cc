/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2017, Yuri Georgievski <ygeorgi@anl.gov>

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

#include "thermo.hh"
#include "graph_omp.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "thermo: ";

  if (argc < 2) {
    std::cout << "usage: thermo input_file\n";
    return 0;
  }

  int                 itemp;
  double              dtemp;
  std::string         stemp;

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup ThermoGroup;

  Key      list_key("TemperatureList[K]"            );
  Key      grid_key("TemperatureGrid(step[K],size)" );
  Key      spec_key("Species[angstrom]"             );
  Key       pot_key("Potential"                     );

  Key   smp_max_key("SamplingNumberMax"             );
  Key   smp_min_key("SamplingNumberMin"             );
  Key     potex_key("PotentialExpansionOrder"       );
  Key  der_step_key("DerivativeStep[bohr]"          );
  Key  smp_step_key("SamplingStep[bohr]"            );
  Key      qlow_key("QFactorFrequencyLow"           );
  Key     qhigh_key("QFactorFrequencyHigh"          );

  Key      bond_key("BondNumberMax"                 );
  Key      ftol_key("FrequencyTolerance"            );
  Key      scut_key("FourierSumCutoff"              );
  Key      redt_key("GraphExFrequencyHigh"          );
  Key      lowf_key("GraphExFrequencyLow"           );
  Key      keep_key("KeepPermutedGraphs"            );

  // temperature list
  //
  std::vector<double> temperature;

  // potential
  //
  Thermo::PotWrap  pot;

  // molecular specification
  //
  std::vector<Atom>  molecule;
  
  // base name
  //
  std::string base_name = argv[1];
  //
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    //
    base_name.resize(base_name.size() - 4);

  // input stream
  //
  std::ifstream from(argv[1]);
  //
  if(!from && base_name == argv[1]) {
    // try inp extension
    //
    stemp = base_name + ".inp";
    //
    from.open(stemp.c_str());
  }

  if(!from) {
    //
    ErrOut err_out;

    err_out << funame << "input file " << argv[1] << " is not found";
  }
  
  IO::log.open((base_name + ".log").c_str());
  //
  IO::out.open((base_name + ".out").c_str());

  std::string token, comment, line, name;
  //
  while(from >> token) {
    //
    // temperature list
    //
    if(list_key == token) {
      //
      if(temperature.size()) {
	//
	ErrOut err_out;

	err_out << funame << token << ": already initialized";
      }
      
      IO::LineInput line_input(from);
      //
      while(line_input >> dtemp)
	//
	temperature.push_back( dtemp * Phys_const::kelv);
      
      if(!temperature.size()) {
	//
	ErrOut err_out;

        err_out << funame << token << ": corrupted\n";
      }
    }
    // temperature grid
    //
    else if(grid_key == token) {
      //
      if(temperature.size()) {
	//
	ErrOut err_out;

	err_out << funame << token <<  ": already initialized";
      }
      
      IO::LineInput lin(from);
      //
      if(!(lin >> dtemp >> itemp))  {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      if(dtemp <= 0. || itemp <= 0) {
	//
	ErrOut err_out;

        err_out << funame << token << ": out of range";
      }
      
      dtemp *=  Phys_const::kelv;
      //
      for(int i = 0; i < itemp; ++i)
	//
	temperature.push_back(double(i + 1) * dtemp);
    }
    // molecular specification 
    //
    else if(spec_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      std::getline(from, comment);
      
      if(itemp < 2) {
	//
	ErrOut err_out;

	err_out << funame << token << ": wrong number of atoms: " << itemp;
      }

      molecule.resize(itemp);
      //
      for(int i = 0; i < molecule.size(); ++i) {
	//
	if(!(from >> molecule[i])) {
	  //
	  ErrOut err_out;

	  err_out << funame << token << ": cannot read " << i + 1 << "-th atom";
	}

	molecule[i] *= Phys_const::angstrom;
      }
    }
    // potential
    //
    else if(pot_key == token) {
      //
      pot = Thermo::new_pot(from);
    }
    // maximal number of importance samplings
    //
    else if(smp_max_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(itemp < 1) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << itemp;
      }

      Thermo::Species::count_max = itemp;
    }
    // importance sampling number to thermolize
    //
    else if(smp_min_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(itemp < 1) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << itemp;
      }

      Thermo::Species::count_min = itemp;
    }
    // maximal potential expansion order
    //
    else if(potex_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(itemp < 3) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << itemp;
      }

      Graph::potex_max = itemp;
    }
    // numerical derivative step (bohr)
    //
    else if(der_step_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(dtemp > 1. || dtemp <= 0.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }

      Thermo::Species::numd_step = dtemp;
    }
    // importance sampling step (bohr)
    //
    else if(smp_step_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(dtemp > 1. || dtemp <= 0.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }

      Thermo::Species::samp_step = dtemp;
    }
    // qfactor low frequency threshold
    //
    else if(qlow_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(dtemp > 1. || dtemp <= 0.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }

      Thermo::Species::low_freq_thresh = dtemp;
    }
    // qfactor high frequency threshold
    //
    else if(qhigh_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      std::getline(from, comment);
      
      if(dtemp < 1.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }

      Thermo::Species::high_freq_thresh = dtemp;
    }
    // frequency tolerance
    //
    else if(ftol_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      if(dtemp > 1. || dtemp < 0.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }
      
      Graph::Expansion::freq_tol = dtemp;

      std::getline(from, comment);
    }	
    // maximal number of bonds in graphs perturbation theory
    //
    else if(bond_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      if(itemp <= 0) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< itemp;
      }
      
      Graph::bond_max = itemp;
      
      std::getline(from, comment);
    }
    // keep permuted graphs in the databases
    //
    else if(keep_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << "corrupted";
      }
      
      if(itemp) {
	//
	Graph::Expansion::mod_flag |= Graph::Expansion::KEEP_PERM;
      }
      else if(Graph::Expansion::mod_flag & Graph::Expansion::KEEP_PERM) {
	//
	Graph::Expansion::mod_flag ^= Graph::Expansion::KEEP_PERM;
      }

      std::getline(from, comment);
    }
    // low temperature / high frequency graph reduction threshold
    //
    else if(redt_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      if(dtemp <= 1.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }
      
      Graph::FreqGraph::red_thresh = dtemp;

      std::getline(from, comment);
    }	
    // fourier sum graph evaluation cutoff
    //
    else if(scut_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      if(itemp < 2) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< itemp;
      }
      
      Graph::FreqGraph::four_cut = itemp;

      std::getline(from, comment);
    }	
    // low frequency threshold
    //
    else if(lowf_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }
      
      if(dtemp  < 0. || dtemp > 1.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< dtemp;
      }
      
      Graph::Expansion::low_freq_thresh = dtemp;

      std::getline(from, comment);
    }	
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      ErrOut err_out;

      err_out << funame << "unknown keyword: " << token << "\n" << Key::show_all();
    }
  }

  if(!temperature.size()) {
    //
    ErrOut err_out;

    err_out << funame << "temperature list not initialized";
  }

  if(!molecule.size()) {
    //
    ErrOut err_out;

    err_out << funame << "molecule not specified";
  }

  if(!pot) {
    //
    ErrOut err_out;

    err_out << funame << "potential not initialized";
  }

  // initialize random number generator
  //
  srand48(std::time(0));

  // initialize graph set
  //
  Graph::init();

  // initialize thermochemical species
  //
  Thermo::CSpec species(molecule, pot);

  for(int i = 0; i < temperature.size(); ++i) {
    //
    double qfac;

    std::map<int, double> tot_corr = species.weight(temperature[i], qfac);

    if(!i) {
      //
      IO::out << std::setw(7)  << "T, K"
	      << std::setw(12) << "QFactor";

      for(std::map<int, double>::const_iterator cit = tot_corr.begin(); cit != tot_corr.end(); ++cit)
	//
	IO::out << std::setw(10)  << "BO = " << std::setw(2) << cit->first;

      IO::out << "\n";
    }
    
    IO::out << std::setw(7)  << temperature[i] / Phys_const::kelv
	    << std::setw(12) << qfac;

    for(std::map<int, double>::const_iterator cit = tot_corr.begin(); cit != tot_corr.end(); ++cit)
      //
      IO::out << std::setw(12)  << cit->second;

    IO::out << std::endl;
  }

  return 0;
}
