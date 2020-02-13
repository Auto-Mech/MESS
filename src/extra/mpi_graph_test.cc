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

#include<mpi.h>

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>

#include "graph_mpi.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "mpi_graph_test: ";

  int         itemp;
  double      dtemp;
  std::string stemp;
  
  MPI::Init(argc, argv);

  const int mpi_size = MPI::COMM_WORLD.Get_size();
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();

  IO::mpi_rank = mpi_rank;

  if (argc < 2) {
    if(!mpi_rank)
      std::cout << "usage: mpi_graph_test input_file" << std::endl;
    
    MPI::Finalize();
    return 0;
  }

  if(mpi_size < 6) {
    if(!mpi_rank)
      std::cout << funame << "can run only with at least six nodes" << std::endl;
    
    MPI::Finalize();
    return 0;
  }

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup GraphGroup;

  Key      list_key("TemperatureList[K]"              );
  Key      grid_key("TemperatureGrid(step[K],size)"   );
  Key      freq_key("Frequencies[1/cm]"               );
  Key     potex_key("PotentialExpansion[1/cm]"        );
  Key      bond_key("BondNumberMax"                   );
  Key      scut_key("FourierSumCutoff"                );
  Key      redt_key("ReductionThreshold"              );
  Key      ftol_key("FrequencyTolerance"              );
  Key      lowf_key("LowFrequencyThreshold"           );
  
  //Key       drv_key("DriversNumber"                   );

  // temperature list
  //
  std::vector<double>  temperature;

  // frequencies
  //
  std::vector<double> frequency;

  // potential expansion
  //
  Graph::potex_t potex;
  
  // base name
  std::string base_name = argv[1];
  //
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    //
    base_name.resize(base_name.size() - 4);

  std::ifstream from(argv[1]);
  //
  if(!from && base_name == argv[1]) {
    //
    // try inp extension
    //
    stemp = base_name + ".inp";

    from.open(stemp.c_str());
  }

  if(!mpi_rank)
    //
    IO::log.open((base_name + ".log").c_str());

  if(!from) {
    //
    ErrOut err_out;

    err_out << funame << "input file " << argv[1] << " is not found";
  }
  
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

	err_out << funame << "temperature list already initialized";
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

	err_out << funame << "temperature list already initialized";
      }
      
      IO::LineInput line_input(from);
      //
      if(!(line_input >> dtemp >> itemp))  {
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
    // frequencies
    //
    else if(freq_key == token) {
      //
      if(frequency.size()) {
	//
	ErrOut err_out;

	err_out << funame << token << ": already initialized";
      }
      
      IO::LineInput line_input(from);
      //
      while(line_input >> dtemp)
	//
	frequency.push_back( dtemp * Phys_const::incm);
      
      if(!frequency.size()) {
	//
	ErrOut err_out;

        err_out << funame << token << ": corrupted";
      }
    }
    // potential expansion
    //
    else if(potex_key == token) {
      //
      if(!frequency.size()) {
	//
	ErrOut err_out;

	err_out << funame << token << ": frequencies should be initialized first";
      }
      
      std::getline(from, comment);
      
      Graph::read_potex(frequency, from, potex);
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
      
      if(itemp < 2) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< itemp;
      }
      
      Graph::bond_max = itemp;
      
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
    // number of drivers
    /*
    else if(drv_key == token) {
      if(!(from >> itemp)) {
	ErrOut err_out;
	err_out << funame << token << ": corrupted";
      }
      
      if(itemp < 1) {
	ErrOut err_out;
	err_out << funame << token << ": out of range: "<< itemp;
      }
      
      Graph::Expansion::WORK_NODE = Graph::Expansion::SUM_SERV + itemp + 1;

      if( Graph::Expansion::WORK_NODE >= mpi_size) {
	ErrOut err_out;
	err_out << funame << token << ": there are no working nodes: increase mpi size or decrease number of drivers";
      }

      std::getline(from, comment);
    }
    */
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      ErrOut err_out;

      err_out << funame << "unknown keyword: " << token << "\n";

      Key::show_all(err_out);
    }
  }

  if(!mpi_rank) {
    //
    IO::log << IO::log_offset << "Number of working nodes = " << mpi_size - 1 << std::endl;

    if(temperature.size()) {
      //
      IO::log << IO::log_offset << "Temperatures[K]:";

      for(int i = 0; i < temperature.size(); ++i)
	//
	IO::log << "   " << temperature[i] / Phys_const::kelv;

      IO::log << "\n";
    }
  
    if(!frequency.size()) {
      //
      ErrOut err_out;

      err_out << funame << "frequencies have not been initialized";
    }
    else {
      //
      IO::log << IO::log_offset << "Frequencies[1/cm]:";
      
      for(int i = 0; i < frequency.size(); ++i)
	//
	IO::log << "   " << frequency[i] / Phys_const::incm;

      IO::log << "\n\n";
    }
  }

  itemp = 0;
  //
  for(Graph::potex_t::const_iterator pit = potex.begin(); pit != potex.end(); ++pit)
    //
    if(pit->first.size() > itemp)
      //
      itemp = pit->first.size();

  Graph::potex_max = itemp;

  Graph::init();

  Graph::Expansion graphex(frequency, potex);

  graphex.correction();

  graphex.correction(Graph::Expansion::CENTROID);

  for(int i = 0; i < temperature.size(); ++i)
    //
    graphex.correction(Graph::Expansion::CENTROID, temperature[i]);

  MPI::Finalize();
  return 0;
}
