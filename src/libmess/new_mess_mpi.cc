/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2019, Yuri Georgievski <ygeorgi@anl.gov>

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
#include <sys/resource.h>

#include "key.hh"
#include "units.hh"
#include "io.hh"

#include "new_mess.hh"
#include "new_comm.hh"
#include "mpack.hh"

#include <mpi.h>

int main (int argc, char* argv [])
{
  const char funame [] = "master_equation: ";

  using namespace MasterEquation;
  
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &IO::mpi_rank);

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;

  // pressure list
  //
  std::vector<double> prs_list;

  // temperature list
  //
  std::vector<double> tmr_list;

  try {
    //
    if (argc < 2) {
      //
      std::cout << "usage: mess input_file\n";
    
      throw Error::Init();
    }

    // base name
    //
    std::string base_name = argv[1];
  
    if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
      //
      base_name.resize(base_name.size() - 4);

    IO::KeyBufferStream from(argv[1]);
  
    if(!from && base_name == argv[1]) {
      //
      // try inp extension
      //
      stemp = base_name + ".inp";
    
      from.open(stemp.c_str());
    }

    if(!from) {
      //
      if(!IO::mpi_rank)
	//
	std::cout << funame << "input file " << argv[1] << " is not found\n";
    
      throw Error::Init();
    }

    // auxiliary output
    //
    IO::aux.open((base_name + ".aux").c_str());
  
    IO::log.open((base_name + ".log").c_str());
  
    IO::out.open((base_name + ".out").c_str());
  
    IO::log << std::setprecision(Model::log_precision);

    IO::out << std::setprecision(Model::out_precision);

    IO::Marker init_marker("mess initialization");
    
    // no core dumps
    //
    rlimit no_core;

    getrlimit(RLIMIT_CORE, &no_core);

    if(no_core.rlim_max > 0) {
      //
      no_core.rlim_cur = 0;
      no_core.rlim_max = 0;

      if(setrlimit(RLIMIT_CORE, &no_core)) {
	//
	IO::log << IO::log_offset << funame << "setrlimit failed\n";

	throw Error::Run();
      }
    }

    /********************************************************************************************
     ************************************* INPUT ************************************************
     ********************************************************************************************/

    KeyGroup Main;

    Key prec_num_key("DecimalPrecision"           );
    Key    float_key("FloatType"                  );
    Key   use_mp_key("UseMultiPrecision"          );
    Key bar_pres_key("PressureList[bar]"          );
    Key tor_pres_key("PressureList[torr]"         );
    Key atm_pres_key("PressureList[atm]"          );
    Key tmp_list_key("TemperatureList[K]"         );
    Key tmp_prog_key("TemperatureProgression[K]"  );
    Key bar_prog_key("PressureProgression[bar]"   );
    Key     etot_key("EnergyStepOverTemperature"  );
    Key dist_wid_key("ThermalWidthFactor"         );
    Key well_ext_key("WellExtension"              );
    Key well_cut_key("WellCutoff"                 );
    Key ext_corr_key("ExtensionCorrection"        );
    Key eval_max_key("ChemicalEigenvalueMax"      );
    Key chem_trs_key("ChemicalThreshold"          );
    Key eval_min_key("ChemicalEigenvalueMin"      );
    Key chem_tol_key("ChemicalTolerance"          );
    Key proj_thd_key("WellProjectionThreshold"    );
    Key proj_tol_key("WellProjectionTolerance"    );
    Key  red_max_key("WellReductionThreshold"     );
    Key  red_min_key("WellReductionTolerance"     );
    Key red_incr_key("WellReductionIncrement"     );
    Key rate_max_key("MicroRateMax[1/sec]"        );
    Key time_max_key("TimePropagationLimit"       );
    Key time_out_key("TimeOutputInterval"         );
    Key time_stp_key("TimePropagationStep"        );
    Key    model_key("Model"                      );
    Key prec_out_key("OutPrecision"               );
    Key prec_log_key("LogPrecision"               );
    Key     emax_key("ModelEnergyLimit[kcal/mol]" );
    Key     xtun_key("TunnelActionMax"            );
    Key  adm_bor_key("AtomDistanceMin[bohr]"      );
    Key  adm_ang_key("AtomDistanceMin[angstrom]"  );

    // input cycle
    //
    std::string token, comment, line;

    int prec_num = 0;

    while(from >> token) {
      //
      // main input group
      //
      if(model_key == token) {
	//
	Model::init(from);
      
	break;
      }
      // out stream precision
      //
      else if(prec_out_key == token) {
	//
	if(Model::isinit()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 2 || itemp > 10) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range[2-10]: " << itemp << "\n";

	  throw Error::Range();
	}

	Model::out_precision = itemp;
      
	IO::out << std::setprecision(Model::out_precision);
      }
      // log stream precision
      //
      else if(prec_log_key == token) {
	//
	if(Model::isinit()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 2 || itemp > 6) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range[2-6]: " << itemp << "\n";

	  throw Error::Range();
	}

	Model::log_precision = itemp;
      
	IO::log << std::setprecision(Model::log_precision);
      }
      // decimal multiple precision
      //
      else if(prec_num_key == token) {
	//
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 16) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range (>15): " << itemp << "\n";

	  throw Error::Range();
	}

	prec_num = itemp;
      }
      // temperature
      //
      else if(tmp_list_key == token) {
	//
	if(tmr_list.size()) {
	  //
	  IO::log << IO::log_offset << funame << token <<  ": already initialized\n";
	
	  throw Error::Init();
	}
      
	IO::LineInput data_input(from);
      
	std::set<double> data;
      
	while(data_input >> dtemp) {
	  //
	  if(dtemp <= 0.) {
	    //
	    IO::log << IO::log_offset << funame << token << ": should be positive\n";
	  
	    throw Error::Range();
	  }
	
	  data.insert( dtemp * Phys_const::kelv);
	}
      
	if(!data.size()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": no data\n";
	
	  throw Error::Init();
	}

	for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it)
	  //
	  tmr_list.push_back(*it);
      }
      // temperature geometric progression
      //
      else if(tmp_prog_key == token) {
	//
	std::string format = ": format: <initial_value> <final_value> <n>";
      
	IO::LineInput lin(from);

	double start, stop;

	int n;

	if(!(lin >> start >> stop >> n)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted" << format << "\n";

	  throw Error::Input();
	}

	if(start <= 0. || stop <= start || n < 2) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << start << " " << stop << " " << n << format << "\n";

	  throw Error::Range();
	}

	double step = std::exp(std::log(stop / start) / double(n - 1));

	double t = start * Phys_const::kelv;
      
	for(int i = 0; i < n; ++i, t *= step) 
	  //
	  tmr_list.push_back(t);
      
      }
      // pressure[bar] geometric progression
      //
      else if(bar_prog_key == token) {
	//
	pressure_unit = BAR;
      
	std::string format = ": format: <initial_value> <final_value> <n>";
      
	IO::LineInput lin(from);

	double start, stop;

	int n;

	if(!(lin >> start >> stop >> n)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted" << format << "\n";

	  throw Error::Input();
	}

	if(start <= 0. || stop <= start || n < 2) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << start << " " << stop << " " << n << format << "\n";

	  throw Error::Range();
	}

	double step = std::exp(std::log(stop / start) / double(n - 1));

	double p = start * Phys_const::bar;
      
	for(int i = 0; i < n; ++i, p *= step) 
	  //
	  prs_list.push_back(p);
      }
      // pressure
      //
      else if(bar_pres_key == token || tor_pres_key == token || atm_pres_key == token) {
	//
	if(prs_list.size()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": already initialized\n";
	
	  throw Error::Init();
	}

	IO::LineInput data_input(from);
      
	std::set<double> data;

	while(data_input >> dtemp) {
	  //
	  if(dtemp <= 0.) {
	    //
	    IO::log << IO::log_offset << funame << token << ": should be positive\n";
	  
	    throw Error::Range();
	  }

	  if(bar_pres_key == token) {
	    //
	    dtemp *= Phys_const::bar;
	  
	    pressure_unit = BAR;
	  }
	
	  if(atm_pres_key == token) {
	    //
	    dtemp *= Phys_const::atm;
	  
	    pressure_unit = ATM;
	  }
	
	  if(tor_pres_key == token) {
	    //
	    dtemp *= Phys_const::tor;
	  
	    pressure_unit = TORR;
	  }

	  data.insert(dtemp);
	}

	if(!data.size()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": no data\n";
	
	  throw Error::Init();
	}

	for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it)
	  //
	  prs_list.push_back(*it);
      }
      // energy step over temperature
      //
      else if(etot_key == token) {
	//
	if(energy_step_over_temperature > 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": already defined\n";
	
	  throw Error::Input();
	}
      
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(dtemp <= 0. || dtemp >= 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range [0, 1]: " << dtemp << "\n";
	
	  throw Error::Range();
	}

	energy_step_over_temperature = dtemp;
      }
      // thermal distribution width factor (for automatic excess energy calculation)
      //
      else if(dist_wid_key == token) {
	//
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(dtemp <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be positive: " << dtemp << "\n";
	
	  throw Error::Range();
	}

	dist_width_factor = dtemp;
      }
      // model energy limit
      //
      else if(emax_key == token) {
	//
	if(Model::isinit()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	Model::set_energy_limit(dtemp * Phys_const::kcal);
      }
      // minimal interatomic distance
      //
      else if(adm_bor_key == token || adm_ang_key == token) {
	//
	if(Model::isinit()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range\n";
	
	  throw Error::Range();
	}
      
	if(adm_ang_key == token)
	  //
	  dtemp *= Phys_const::angstrom;

	Model::atom_dist_min = dtemp;
      }
      // maximum tunneling exponent
      //
      else if(xtun_key == token) {
	//
	if(Model::isinit()) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}

	Model::Tunnel::set_action_max(dtemp);
      }
      // microcanonical rate limit
      //
      else if(rate_max_key == token) {
	//
	if(!(from >> micro_rate_max)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(micro_rate_max <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}

	micro_rate_max *= Phys_const::herz;
      }
      // time propagation limit in collision frequency units
      //
      else if(time_max_key == token) {
	//
	if(!(from >> time_limit)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(time_limit <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}
      }
      // time output interval in collision frequency units
      //
      else if(time_out_key == token) {
	//
	if(!(from >> time_output)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(time_output <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}
      }
      // time propagation step in fastest relaxation time-scale = collision frequency * reduction threshold
      //
      else if(time_stp_key == token) {
	//
	if(!(from >> time_step)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(time_step <= 0.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}
      }
      // well cutoff
      //
      else if(well_cut_key == token) {
	//
	IO::LineInput lin(from);

	if(!(lin >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";

	  throw Error::Input();
	}
	
	if(dtemp < 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << dtemp << "\n";
	
	  throw Error::Range();
	}

	well_cutoff = dtemp;
      }
      // well extention
      //
      else if(well_ext_key == token) {
	//
	IO::LineInput lin(from);

	if(lin >> dtemp) {
	  //
	  if(dtemp < 0. || dtemp >= 1.) {
	    //
	    IO::log << IO::log_offset << funame << token << ": out of range: " << dtemp << "\n";
	
	    throw Error::Range();
	  }

	  well_extension = dtemp;
	}
	else
	  //
	  well_extension = 0.;
      }
      // well extention correction factor
      //
      else if(ext_corr_key == token) {
	//
	IO::LineInput lin(from);

	dtemp = 0.;

	if(!(lin >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";

	  throw Error::Input();
	}
      
	if(dtemp <= 0. || dtemp >= 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << dtemp << "\n";
	
	  throw Error::Range();
	}

	well_ext_corr = dtemp;
      }
      // fast isomerization threshold
      //
      else if(red_max_key == token) {
	//
	if(!(from >> reduction_threshold)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(reduction_threshold < 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << reduction_threshold << "\n";

	  throw Error::Range();
	}
      }
      // slow isomerization limit
      //
      else if(red_min_key == token) {
	//
	if(!(from >> reduction_tolerance)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(reduction_tolerance > 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << reduction_tolerance << "\n";

	  throw Error::Range();
	}
      }
      // gap between highest slow isomerization eigenvalue and the lowest active eigenvalue
      //
      else if(red_incr_key == token) {
	//
	if(!(from >> reduction_increment)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(reduction_increment < 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range: " << reduction_increment << "\n";

	  throw Error::Range();
	}
      }
      // float type
      //
      else if(float_key == token) {
	//
#if not defined(WITH_MPACK) && not defined(WITH_MPLAPACK)

	IO::log << IO::log_offset << funame << token << ": to use multiple precision the code should be compiled with WITH_MPACK or WITH_MPLAPACK macro defined\n";

	throw Error::Logic();
      
#endif
      
	IO::LineInput lin(from);

	if(!(lin >> stemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(stemp == "double") {
	  //
	  Mpack::mp_type = Mpack::DOUBLE;
	}
	else if(stemp == "dd") {
	  //
	  Mpack::mp_type = Mpack::DD;
	}
	else if(stemp == "qd") {
	  //
	  Mpack::mp_type = Mpack::QD;
	}
	else if(stemp == "mpfr") {
	  //
	  Mpack::mp_type = Mpack::MPFR;
	}
	else if(stemp == "gmp") {
	  //
	  Mpack::mp_type = Mpack::GMP;
	}
	else if(stemp == "float128") {
	  //
	  Mpack::mp_type = Mpack::FLOAT128;
	}
	else if(stemp == "float64x") {
	  //
	  Mpack::mp_type = Mpack::FLOAT64X;
	}
	else {
	  //
	  IO::log << IO::log_offset << funame << token << ": unknown type: " << stemp << ": available types: double (default), dd, qd, mpfr, gmp, float128, float64x\n";

	  throw Error::Init();
	}
      }
      // use multi-precision for global kinetic matrix diagonalization
      //
      else if(use_mp_key == token) {
	//
	use_mp = 1;
      }
      // maximal chemical eigenvalue a.k.a. chemical threshold
      //
      else if(eval_max_key == token || chem_trs_key == token) {
	//
	IO::LineInput lin(from);
      
	if(!(lin >> chemical_threshold)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      }
      // minumal chemical eigenvalue a.k.a. chemical tolerance
      //
      else if(eval_min_key == token || chem_tol_key == token) {
	//
	if(!(from >> chemical_tolerance)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);
      }
      // well projection threshold
      //
      else if(proj_thd_key == token) {
	//
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0. || dtemp >= 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range\n";
	
	  throw Error::Range();
	}
	
	well_projection_threshold = dtemp;
      }
      // well projection tolerance
      //
      else if(proj_tol_key == token) {
	//
	if(!(from >> dtemp)) {
	  //
	  IO::log << IO::log_offset << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0. || dtemp >= 1.) {
	  //
	  IO::log << IO::log_offset << funame << token << ": out of range\n";
	
	  throw Error::Range();
	}
	
	well_projection_tolerance = dtemp;
      }
      // unknown keyword
      //
      else if(IO::skip_comment(token, from)) {
	//
	IO::log << IO::log_offset << funame << "unknown keyword " << token << "\n";
      
	Key::show_all(IO::log, IO::log_offset);
      
	throw Error::Init();
      }
    }

    if(prec_num)
      //
      set_precision(prec_num);
  
    IO::log << IO::log_offset << "numeric precision: digits10 = " << get_precision() << "\n";

    /*************************** CHECKING ******************************************/

    if(!from) {
      //
      IO::log << IO::log_offset << funame << "corrupted\n";
    
      throw Error::Input();
    }

    if(!Model::isinit()) {
      //
      IO::log << IO::log_offset << funame << "model is not initialized\n";
    
      throw Error::Init();
    }

    if(!tmr_list.size()) {
      //
      IO::log << IO::log_offset << funame << "temperature list has not been initialized\n";
    
      throw Error::Init();
    }

    if(!prs_list.size()) {
      //
      IO::log << IO::log_offset << funame << "pressure list has not been initialized\n";
    
      throw Error::Init();
    }

    if(energy_step_over_temperature < 0.) {
      //
      IO::log << IO::log_offset << funame << "energy step over temperature has not been initialized\n";
    
      throw Error::Init();
    }

    if(reduction_threshold < 0.) {
      //
      IO::log << IO::log_offset << funame << "fast isomerization reduction threshold has not been initialized\n";

      throw Error::Init();
    }

    time_step /= reduction_threshold;
    
    if(!IO::mpi_rank) {
      //
      Model::pes_print();

      Model::pf_print();
    
      Model::names_translation(IO::out);
  
      Model::names_translation(IO::log);
    }
  }
  // initialization failed
  //
  catch(Error::General) {
    //
    if(!IO::mpi_rank)
      //
      std::cout << funame << "mess initialization failed, see the log file\n";

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
  }
  
  if(Model::no_run()) {
    //
    if(!IO::mpi_rank)
      //
      std::cout << funame << "no rate calculation required, exiting\n";

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
  }

  int job;

  MPI_Status stat;

  IO::log.clear();

  IO::aux.clear();

  /**********************************************************************************
   ************************************* SLAVE **************************************
   **********************************************************************************/
  
  if(IO::mpi_rank) {
    //
    IO::log_offset.increase();
    
    while(1) {
      //
      MPI_Recv(&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

      if(stat.MPI_TAG == Comm::STOP_TAG) {
	//
	break;
      }

      try {
	//
	temperature = tmr_list[job / prs_list.size()];

	pressure    = prs_list[job % prs_list.size()];

	std::list<RateData> rate_data;

	get_rate_data(rate_data);

	MPI_Send(&job, 1, MPI_INT, 0, Comm::DATA_TAG, MPI_COMM_WORLD);

	Comm::send_rate_data(rate_data);
      }
      // non-critical error
      //
      catch(Error::Run) {
	//
	MPI_Send(&job, 1, MPI_INT, 0, Comm::FAIL_TAG, MPI_COMM_WORLD);
      }
      // unrecoverable error
      //
      catch(Error::General) {
	//
	MPI_Send(&job, 1, MPI_INT, 0, Comm::STOP_TAG, MPI_COMM_WORLD);
      }
      
      Comm::send_log();
    }
  }
  /***********************************************************************************
   ************************************* MASTER **************************************
   ***********************************************************************************/
  else {
    //
    IO::Marker rate_marker("rate calculation");

    const int job_size = tmr_list.size() * prs_list.size();

    std::vector<RateData> rate_data(job_size);

    MPI_Comm_size(MPI_COMM_WORLD, &itemp);

    const int mpi_size = itemp;
    
    IO::log << IO::log_offset << "nodes # = " << mpi_size << "\n";
    
    IO::log << IO::log_offset << "jobs  # = " << job_size << "\n";
    
    int running_nodes = mpi_size - 1;

    int stop_work = 0;
    
    if(job_size < running_nodes) {
      //
      IO::log << IO::log_offset << "WARNING: number of jobs (" << job_size
	//
	      << ") is less than the number of working nodes (" << running_nodes << ")\n";
      
      for(int i = job_size; i < running_nodes; ++i)
	//
	MPI_Send(&i, 1, MPI_INT, i + 1, Comm::STOP_TAG, MPI_COMM_WORLD);
    
      running_nodes = job_size;
    }

    for(int i = 0; i < running_nodes; ++i)
      //
      MPI_Send(&i, 1, MPI_INT, i + 1, Comm::RUN_TAG, MPI_COMM_WORLD);
    
    int new_job = running_nodes;

    // new job cycle
    //
    while(running_nodes) {
      //
      MPI_Recv(&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
      
      const int node = stat.MPI_SOURCE;

      const int tag = stat.MPI_TAG;

      if(tag == Comm::DATA_TAG) {
	//
	std::list<RateData> rate_data_list;
	
	Comm::recv_rate_data(node, rate_data_list);

	Comm::recv_log(node);

	try {
	  //
	  IO::log << "\n";

	  aggregate_rate_data(rate_data_list, rate_data[job]);

	  IO::log << IO::log_offset << job << "-th job is done (node = " << node << ")\n\n";
	}
	catch(Error::General) {
	  //
	  IO::log << "\n" << IO::log_offset << "rate data aggregation for " << job << "-th job failed: stop work request initiated\n\n";
	  
	  stop_work = 1;
	}
      }
      // job failed 
      //
      else if(tag == Comm::FAIL_TAG) {
	//
	Comm::recv_log(node);
      
	IO::log << "\n" << IO::log_offset << job << "-th job failed (node = " << node << ")\n\n";
      }
      // unrecoverable error
      //
      else if(tag == Comm::STOP_TAG) {
	//
	Comm::recv_log(node);

	IO::log << "\n" << IO::log_offset << job << "-th job (" << node << "-th node): unrecoverable error\n\n";
	
	stop_work = 1;
      }
      
      // remove working node
      //
      if(new_job == job_size || stop_work) {
	//
	--running_nodes;

	MPI_Send(&new_job, 1, MPI_INT, node, Comm::STOP_TAG, MPI_COMM_WORLD);
      }
      // start new job
      //
      else {
	//
	MPI_Send(&new_job, 1, MPI_INT, node, Comm::RUN_TAG, MPI_COMM_WORLD);

	++new_job;
      }//
      //
    }// new job cycle

    if(stop_work) {
      //
      IO::log << IO::log_offset << "calculation failed => exitting\n";

      std::cerr << funame << "calculation failed, see the log file\n";
    }    
    /***********************************************************
     ************************* OUTPUT **************************
     ***********************************************************/
    else {
      //
      IO::log << IO::log_offset << "all jobs are done => exitting\n\n";

      const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

      IO::out << "Species-Species Tables:\n\n";
      
      for(job = 0; job < job_size; ++job) {
	//
	pressure    = prs_list[job % prs_list.size()];

	temperature = tmr_list[job / prs_list.size()];
	
	IO::out << "pressure = ";
  
	switch(pressure_unit) {
	  //
	case BAR:
	  //
	  IO::out << pressure / Phys_const::bar << " bar";
    
	  break;
	  //
	case TORR:
	  //
	  IO::out << pressure / Phys_const::tor << " torr";
    
	  break;
	  //
	case ATM:
	  //
	  IO::out << pressure / Phys_const::atm << " atm";
    
	  break;
	}

	IO::out << "\t temperature = " <<  temperature / Phys_const::kelv << " K\n\n";

	std::map<std::set<int>, int>::const_iterator fit, git;

	const std::map<std::set<int>, int>& map = rate_data[job].group_index_map;

	IO::out << std::setw(Model::out_precision + 7) << "S\\S";
	
	for(git = map.begin(); git != map.end(); ++git)
	  //
	  IO::out << std::setw(Model::out_precision + 7) << Model::group_name(git->first);

	for(int b = 0; b < Model::bimolecular_size(); ++b)
	  //
	  IO::out << std::setw(Model::out_precision + 7) << Model::bimolecular(b).name();

	IO::out << "\n";
	
	for(fit = map.begin(); fit != map.end(); ++fit) {
	  //
	  IO::out << std::setw(Model::out_precision + 7) << Model::group_name(fit->first);
	
	  for(git = map.begin(); git != map.end(); ++git)
	    //
	    IO::out << std::setw(Model::out_precision + 7) << rate_data[job].ww_rate(fit->second, git->second)
	      //
	      / Model::weight(fit->first, temperature) / Phys_const::herz;

	  for(int b = 0; b < Model::bimolecular_size(); ++b)
	    //
	    IO::out << std::setw(Model::out_precision + 7) << rate_data[job].wb_rate(fit->second, b)
	      //
	      / Model::weight(fit->first, temperature) / Phys_const::herz;

	  IO::out << "\n";
	}

	for(int c = 0; c < Model::bimolecular_size(); ++c) {
	  //
	  IO::out << std::setw(Model::out_precision + 7) << Model::bimolecular(c).name();
	
	  for(git = map.begin(); git != map.end(); ++git)
	    //
	    if(Model::bimolecular(c).dummy()) {
	      //
	      IO::out << std::setw(Model::out_precision + 7) << "***";
	    }
	    else
	      //
	      IO::out << std::setw(Model::out_precision + 7) << rate_data[job].bw_rate(c, git->second) / bru
		//
		/ Model::bimolecular(c).weight(temperature) * std::exp(Model::bimolecular(c).ground() / temperature);

	  for(int b = 0; b < Model::bimolecular_size(); ++b)
	    //
	    if(Model::bimolecular(c).dummy()) {
	      //
	      IO::out << std::setw(Model::out_precision + 7) << "***";
	    }
	    else
	      //
	      IO::out << std::setw(Model::out_precision + 7) << rate_data[job].bb_rate(c, b) / bru
		//
		/ Model::bimolecular(c).weight(temperature) * std::exp(Model::bimolecular(c).ground() / temperature);
	      
	  IO::out << "\n";
	}
	IO::out << "\n";
      }//
      //
    }//
    //
  }//master node

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  
  return 0;
}
