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
#include <cstdlib>

//#undef INT
#include <mpi.h>

#include "libmess/mess.hh"
#include "libmess/key.hh"
#include "libmess/units.hh"
#include "libmess/io.hh"
#include "libmess/offload.hh"
#include "libmess/comm.hh"
#include "libmess/mpack.hh"
#include "libmess/limits.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "master_equation: ";

  if (argc < 2) {
    std::cout << "usage: mess input_file\n";
    return 0;
  }

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
      std::cerr << funame << "setrlimit failed\n";

      return 1;
    }
  }

  MPI_Init(&argc, &argv);

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  std::pair<int, int> ptemp;

  MPI_Comm_rank(MPI_COMM_WORLD, &IO::mpi_rank);

  MPI_Comm_size(MPI_COMM_WORLD, &itemp);

  const int mpi_size = itemp;
    
  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key    model_key("Model"                      );
  Key prec_out_key("OutPrecision"               );
  Key prec_ped_key("PedPrecision"               );
  Key prec_log_key("LogPrecision"               );
  Key prec_num_key("NumericPrecision"           );
  Key bar_pres_key("PressureList[bar]"          );
  Key tor_pres_key("PressureList[torr]"         );
  Key atm_pres_key("PressureList[atm]"          );
  Key tmp_list_key("TemperatureList[K]"         );
  Key tmp_prog_key("TemperatureProgression[K]"  );
  Key bar_prog_key("PressureProgression[bar]"   );
  Key     etot_key("EnergyStepOverTemperature"  );
  Key     xtot_key("ExcessEnergyOverTemperature");
  Key ener_cut_key("EnergyCutoffOverTemperature");
  Key cut_kcal_key("GlobalCutoff[kcal/mol]"     );
  Key cut_incm_key("GlobalCutoff[1/cm]"         );
  Key   cut_kj_key("GlobalCutoff[kJ/mol]"       );
  Key lim_kcal_key("ModelEnergyLimit[kcal/mol]" );
  Key lim_incm_key("ModelEnergyLimit[1/cm]"     );
  Key   lim_kj_key("ModelEnergyLimit[kJ/mol]"   );
  Key  adm_bor_key("AtomDistanceMin[bohr]"      );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );
  Key     xtun_key("TunnelActionMax"            );
  Key well_cut_key("WellCutoff"                 );
  Key well_ext_key("WellExtension"              );
  Key ext_corr_key("ExtensionCorrection"        );
  Key eval_max_key("ChemicalEigenvalueMax"      );
  Key chem_trs_key("ChemicalThreshold"          );
  Key eval_min_key("ChemicalEigenvalueMin"      );
  Key chem_tol_key("ChemicalTolerance"          );
  Key well_red_key("WellReductionThreshold"     );
  Key    float_key("FloatType"                  );
  Key   use_mp_key("UseMultiPrecision"          );
  Key     calc_key("CalculationMethod"          );
  Key  rat_red_key("ReductionMethod"            );
  Key      wpm_key("WellPartitionMethod"        );
  Key      wpt_key("WellProjectionThreshold"    );
  Key eval_out_key("EigenvalueOutput"           );
  Key evec_num_key("EigenvectorNumber"          );
  Key evec_out_key("EigenvectorOutput"          );
  Key  red_out_key("ReductionNumber"            );
  Key ped_spec_key("PEDSpecies"                 );
  Key  ped_out_key("PEDOutput"                  );
  Key    react_key("Reactant"                   );
  Key  def_red_key("DefaultReductionScheme"     );
  Key def_chem_key("DefaultChemicalSize"        );
  Key hot_kcal_key("HotEnergies[kcal/mol]"      );
  Key hot_incm_key("HotEnergies[1/cm]"          );
  Key   hot_kj_key("HotEnergies[kJ/mol]"        );
  Key rate_max_key("MicroRateMax[1/sec]"        );
  Key mic_rout_key("MicroRateOutput"            );
  Key mic_emax_key("MicroEnerMax[kcal/mol]"     );
  Key mic_emin_key("MicroEnerMin[kcal/mol]"     );
  Key mic_step_key("MicroEnerStep[kcal/mol]"    );
  Key tim_evol_key("TimeEvolution"              );
  Key       sl_key("StateLandscape"             );
  Key save_mat_key("SaveKineticMatrix"          );
  Key gpu_size_key("GpuSize"                    );
  Key  exp_max_key("ExpArgMax"                  );

  int gpu_size = 0;
  std::vector<std::string> ped_spec;// product energy distribution pairs verbal
  std::vector<std::string> reduction_scheme;
  std::vector<double> pressure;
  std::vector<double> temperature;
  double xtot   = -1.; // exsess energy over temperature
  MasterEquation::Method method = MasterEquation::direct_diagonalization_method;
  std::string micro_rate_file;
  double micro_ener_max  = 0.;
  double micro_ener_min  = 0.;
  double micro_ener_step = -1.;
  std::string state_landscape;

  // input cycle
  //
  std::string token, comment, line;

  int prec_num = 0;

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

  try {
    //
    if(!from) {
      //
      std::cerr << funame << "input file " << argv[1] << " is not found\n";
    
      throw Error::Init();
    }

    // auxiliary output
    //
    IO::aux.open((base_name + ".aux").c_str());
  
    IO::log.open((base_name + ".log").c_str());
  
    IO::log << std::setprecision(Model::log_precision);

    IO::out.open((base_name + ".out").c_str());
  
    IO::out << std::setprecision(Model::out_precision);

    while(from >> token) {
      //
      // main input group
      //
      if(model_key == token) {
	//
	Model::init(from);
      
	break;
      }
      // exponent argument max
      //
      else if(exp_max_key == token) {
	//
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted" << std::endl;

	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0.) {
	  //
	  std::cerr << funame << token << ": out of range: " << dtemp << std::endl;

	  throw Error::Range();
	}
        
	Limits::set_exp_pow_max(dtemp);
      }
      // gpu size
      //
      else if(gpu_size_key == token) {
	//
	if(!(from >> gpu_size)) {
	  //
	  std::cerr << funame << token << ": corrupted" << std::endl;

	  throw Error::Input();
	}

	if(gpu_size < 1) {
	  //
	  std::cerr << funame << token << ": out of range: " << gpu_size << std::endl;

	  throw Error::Range();
	}
	
	std::getline(from, comment);
      }
      // out stream precision
      //
      else if(prec_out_key == token) {
	//
	if(Model::isinit()) {
	  //
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 2 || itemp > 10) {
	  //
	  std::cerr << funame << token << ": out of range[2-10]: " << itemp << "\n";

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
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 2 || itemp > 6) {
	  //
	  std::cerr << funame << token << ": out of range[2-6]: " << itemp << "\n";

	  throw Error::Range();
	}

	Model::log_precision = itemp;
      
	IO::log << std::setprecision(Model::log_precision);
      }
      // calculation precision
      //
      else if(prec_num_key == token) {
	//
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 16) {
	  //
	  std::cerr << funame << token << ": out of range (>15): " << itemp << "\n";

	  throw Error::Range();
	}

	prec_num = itemp;
      }
      // ped stream precision
      //
      else if(prec_ped_key == token) {
	//
	if(Model::isinit()) {
	  //
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	IO::LineInput lin(from);

	if(!(lin >> itemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(itemp < 2 || itemp > 6) {
	  //
	  std::cerr << funame << token << ": out of range[2-6]: " << itemp << "\n";

	  throw Error::Range();
	}

	Model::ped_precision = itemp;
      }
      // eigenvalue output
      //
      else if(eval_out_key == token) {
	//
	if(MasterEquation::eval_out.is_open()) {
	  //
	  std::cerr << funame << token << ": allready opened\n";
	
	  throw Error::Init();
	}      
	if(!(from >> stemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	MasterEquation::eval_out.open(stemp.c_str());
      
	if(!MasterEquation::eval_out) {
	  //
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	
	  throw Error::Input();
	}
      }
      // eigenvector output
      //
      else if(evec_out_key == token) {
	//
	if(MasterEquation::evec_out.is_open()) {
	  //
	  std::cerr << funame << token << ": allready opened\n";
	
	  throw Error::Init();
	}
      
	if(!(from >> stemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	MasterEquation::evec_out.open(stemp.c_str());
      
	if(!MasterEquation::evec_out) {
	  //
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	
	  throw Error::Input();
	}
      }
      // product energy distribution output
      //
      else if(ped_out_key == token) {
	//
	if(MasterEquation::ped_out.is_open()) {
	  //
	  std::cerr << funame << token << ": already opened\n";
	
	  throw Error::Open();
	}
      
	if(!(from >> stemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	MasterEquation::ped_out.open(stemp.c_str());
      
	if(!MasterEquation::ped_out.is_open()) {
	  //
	  std::cerr << funame << token << ": cannot open the " << stemp << " file\n";
	
	  throw Error::Open();
	}
      }
      // reactants and products for product energy distribution output
      //
      else if(ped_spec_key == token) {
	//
	IO::LineInput ped_input(from);
      
	while(ped_input >> stemp)
	  //
	  ped_spec.push_back(stemp);
      }
      // bimolecular reactant to be used as a reference
      //
      else if(react_key == token) {
	//
	if(Model::isinit()) {
	  //
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> Model::reactant)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);      
      }
      // eigenvector number to output
      //
      else if(evec_num_key == token) {
	//
	if(!(from >> MasterEquation::evec_out_num)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(MasterEquation::evec_out_num < 0) {
	  //
	  std::cerr << funame << token << ": out of range\n";
	
	  throw Error::Range();
	}
      }
      // temperature
      //
      else if(tmp_list_key == token) {
	//
	if(temperature.size()) {
	  //
	  std::cerr << funame << token <<  ": already initialized\n";
	
	  throw Error::Init();
	}
      
	IO::LineInput data_input(from);
      
	std::set<double> data;
      
	while(data_input >> dtemp) {
	  //
	  if(dtemp <= 0.) {
	    //
	    std::cerr << funame << token << ": should be positive\n";
	  
	    throw Error::Range();
	  }
	
	  data.insert( dtemp * Phys_const::kelv);
	}
      
	if(!data.size()) {
	  //
	  std::cerr << funame << token << ": no data\n";
	
	  throw Error::Init();
	}

	for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it)
	  //
	  temperature.push_back(*it);
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
	  std::cerr << funame << token << ": corrupted" << format << "\n";

	  throw Error::Input();
	}

	if(start <= 0. || stop <= start || n < 2) {
	  //
	  std::cerr << funame << token << ": out of range: " << start << " " << stop << " " << n << format << "\n";

	  throw Error::Range();
	}

	double step = std::exp(std::log(stop / start) / double(n - 1));

	double t = start * Phys_const::kelv;
      
	for(int i = 0; i < n; ++i, t *= step) 
	  //
	  temperature.push_back(t);
      
      }
      // pressure[bar] geometric progression
      //
      else if(bar_prog_key == token) {
	//
	std::string format = ": format: <initial_value> <final_value> <n>";
      
	IO::LineInput lin(from);

	double start, stop;

	int n;

	if(!(lin >> start >> stop >> n)) {
	  //
	  std::cerr << funame << token << ": corrupted" << format << "\n";

	  throw Error::Input();
	}

	if(start <= 0. || stop <= start || n < 2) {
	  //
	  std::cerr << funame << token << ": out of range: " << start << " " << stop << " " << n << format << "\n";

	  throw Error::Range();
	}

	double step = std::exp(std::log(stop / start) / double(n - 1));

	double p = start * Phys_const::bar;
      
	for(int i = 0; i < n; ++i, p *= step) 
	  //
	  pressure.push_back(p);
      }
      // pressure
      //
      else if(bar_pres_key == token || tor_pres_key == token || atm_pres_key == token) {
	//
	if(pressure.size()) {
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	  throw Error::Init();
	}

	IO::LineInput data_input(from);
      
	std::set<double> data;

	while(data_input >> dtemp) {
	  //
	  if(dtemp <= 0.) {
	    //
	    std::cerr << funame << token << ": should be positive\n";
	  
	    throw Error::Range();
	  }

	  if(bar_pres_key == token) {
	    //
	    dtemp *= Phys_const::bar;
	  
	    MasterEquation::pressure_unit = MasterEquation::BAR;
	  }
	
	  if(atm_pres_key == token) {
	    //
	    dtemp *= Phys_const::atm;
	  
	    MasterEquation::pressure_unit = MasterEquation::ATM;
	  }
	
	  if(tor_pres_key == token) {
	    //
	    dtemp *= Phys_const::tor;
	  
	    MasterEquation::pressure_unit = MasterEquation::TORR;
	  }

	  data.insert(dtemp);
	}

	if(!data.size()) {
	  //
	  std::cerr << funame << token << ": no data\n";
	
	  throw Error::Init();
	}

	for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it)
	  //
	  pressure.push_back(*it);
      }
      // energy step over temperature
      //
      else if(etot_key == token) {
	//
	if(MasterEquation::energy_step_over_temperature > 0.) {
	  //
	  std::cerr << funame << token << ": already defined\n";
	
	  throw Error::Input();
	}
      
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(dtemp <= 0. || dtemp >= 1.) {
	  //
	  std::cerr << funame << token << ": out of range: [0, 1]\n";
	
	  throw Error::Range();
	}

	MasterEquation::energy_step_over_temperature = dtemp;
      }
      // excess energy over temperature
      //
      else if(xtot_key == token) {
	//
	if(xtot > 0.) {
	  //
	  std::cerr << funame << "reference energy has been already defined\n";
	
	  throw Error::Input();
	}
      
	if(!(from >> xtot)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);

	if(xtot <= 0.) {
	  //
	  std::cerr << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}
      }
      // energy cutoff over temperature (inactive)
      //
      else if(ener_cut_key == token) {
	//
	if(MasterEquation::energy_cutoff_over_temperature > 0.) {
	  //
	  std::cerr << funame << token << ": already initialized\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);
      
	if(dtemp <= 0.) {
	  //
	  std::cerr << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}
	
	MasterEquation::energy_cutoff_over_temperature = dtemp;
      }
      // model energy limit
      //
      else if(lim_kcal_key == token || lim_incm_key == token || lim_kj_key == token) {
	//
	if(Model::isinit()) {
	  //
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(lim_kcal_key == token)
	  //
	  dtemp *= Phys_const::kcal;

	if(lim_kj_key == token)
	  //
	  dtemp *= Phys_const::kjoul;

	if(lim_incm_key == token)
	  //
	  dtemp *= Phys_const::incm;

	Model::set_energy_limit(dtemp);
      }
      // global energy cutoff
      //
      else if(cut_incm_key == token || cut_kcal_key == token || cut_kj_key == token) {
	//
	if(MasterEquation::is_global_cutoff) {
	  //
	  std::cerr << funame << token << ": already initialized\n";

	  throw Error::Init();
	}

	MasterEquation::is_global_cutoff = true;

	IO::LineInput lin(from);

	if(!(lin >> dtemp)) {
	  //
	  std::cerr << funame << token << ": lower cutoff corrupted\n";
	
	  throw Error::Input();
	}
	    
	if(cut_incm_key == token)
	  //
	  dtemp *= Phys_const::incm;

	if(cut_kcal_key == token)
	  //
	  dtemp *= Phys_const::kcal;

	if(cut_kj_key == token)
	  //
	  dtemp *= Phys_const::kjoul;

	MasterEquation::lower_global_cutoff = dtemp;
      }
      // minimal interatomic distance
      //
      else if(adm_bor_key == token || adm_ang_key == token) {
	//
	if(Model::isinit()) {
	  //
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0) {
	  //
	  std::cerr << funame << token << ": out of range\n";
	
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
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0.) {
	  //
	  std::cerr << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}

	Model::Tunnel::set_action_max(dtemp);
      }
      // calculation method
      //
      else if(calc_key == token) {
	//
	if(!(from >> stemp)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(stemp == "direct")
	  method = MasterEquation::direct_diagonalization_method;
	else if(stemp == "low-eigenvalue")
	  method = MasterEquation::low_eigenvalue_method;
	else if(stemp == "well-reduction")
	  method = MasterEquation::well_reduction_method;
	else {
	  std::cerr << funame << token << ": unknown method: " << stemp << ": available methods: direct, low-eigenvalue, well-reduction\n";
	  throw Error::Input();
	}
      }
      // microcanonical rate limit
      //
      else if(rate_max_key == token) {
	if(!(from >> MasterEquation::rate_max)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(MasterEquation::rate_max <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}

	MasterEquation::rate_max *= Phys_const::herz;
      }
      // well cutoff
      //
      else if(well_cut_key == token) {
	//
	if(!(from >> MasterEquation::well_cutoff)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(MasterEquation::well_cutoff <= 0.) {
	  //
	  std::cerr << funame << token << ": should be positive\n";
	
	  throw Error::Range();
	}
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
	    std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	    throw Error::Range();
	  }

	  MasterEquation::well_extension = dtemp;
	}
	else
	  //
	  MasterEquation::well_extension = 0.;
      }
      // well extention correction factor
      //
      else if(ext_corr_key == token) {
	//
	IO::LineInput lin(from);

	dtemp = 0.;

	if(!(lin >> dtemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

	  throw Error::Input();
	}
      
	if(dtemp <= 0. || dtemp >= 1.) {
	  //
	  std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	  throw Error::Range();
	}

	MasterEquation::well_ext_corr = dtemp;
      }
      // chemical relaxation to collision frequency ratio
      //
      else if(well_red_key == token) {
	//
	if(!(from >> MasterEquation::reduction_threshold)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
	std::getline(from, comment);
      }
      // use multi-precision for global kinetic matrix diagonalization
      //
      else if(use_mp_key == token) {
	//
	MasterEquation::use_mp = 1;

	std::getline(from, comment);
      }
      // float type
      //
      else if(float_key == token) {
	//
#if not defined(WITH_MPACK) && not defined(WITH_MPLAPACK)

	std::cerr << funame << token << ": to use multiple precision the code should be compiled with WITH_MPACK or WITH_MPLAPACK macro defined\n";

	throw Error::Logic();
      
#endif
      
	IO::LineInput lin(from);

	if(!(lin >> stemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

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
	  std::cerr << funame << token << ": unknown type: " << stemp << ": available types: double (default), dd, qd, mpfr, gmp, float128, float64x\n";

	  throw Error::Init();
	}
      }
      // save kinetic matrix
      //
      else if(save_mat_key == token) {
	//
	IO::LineInput lin(from);
      
	if(!(lin >> MasterEquation::save_kinetic_matrix)) {
	  //
	  MasterEquation::save_kinetic_matrix = base_name + ".skm";
	}

	std::remove(MasterEquation::save_kinetic_matrix.c_str());
      }
      // maximal chemical eigenvalue a.k.a. chemical threshold
      //
      else if(eval_max_key == token || chem_trs_key == token) {
	//
	IO::LineInput lin(from);
      
	if(!(lin >> MasterEquation::chemical_threshold)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      }
      // minumal chemical eigenvalue
      //
      else if(eval_min_key == token || chem_tol_key == token) {
	//
	if(!(from >> MasterEquation::chemical_tolerance)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	std::getline(from, comment);
      }
      // number of closest reductions to print
      //
      else if(red_out_key == token) {
	if(!(from >> MasterEquation::red_out_num)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(MasterEquation::red_out_num < 1) {
	  std::cerr << funame << token << ": out of range\n";
	  throw Error::Range();
	}
      }
      // species reduction algorithm for low-eigenvalue method
      else if(rat_red_key == token) {
	if(!(from >> stemp)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(stemp == "diagonalization")
	  MasterEquation::reduction_method = MasterEquation::DIAGONALIZATION;
	else if(stemp == "projection")
	  MasterEquation::reduction_method = MasterEquation::PROJECTION;
	else {
	  std::cerr << funame << token << ": unknown reduction method: " << stemp 
		    << "; available methods: diagonalization, projection\n";
	  throw Error::Range();
	}
      }
      // well partition method
      else if(wpm_key == token) {
	if(!(from >> stemp)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
	std::getline(from, comment);

	MasterEquation::set_well_partition_method(stemp);
      }
      // well partition threshold
      else if(wpt_key == token) {
	if(!(from >> dtemp)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
	std::getline(from, comment);

	if(dtemp <= 0. || dtemp >= 1.) {
	  std::cerr << funame << token << ": out of range\n";
	  throw Error::Range();
	}
	
	MasterEquation::well_projection_threshold = dtemp;
      }
      // default reduction scheme
      else if(def_red_key == token) {
	IO::LineInput scheme_input(from);
	while(scheme_input >> stemp)
	  reduction_scheme.push_back(stemp);
      }
      // default chemical subspace size
      //
      else if(def_chem_key == token) {
	//
	IO::LineInput lin(from);
      
	if(!(lin >> itemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      
	MasterEquation::set_default_chem_size(itemp);
      }
      // hot energies
      //
      else if(hot_kcal_key == token || hot_incm_key == token || hot_kj_key == token) {
	//
	IO::LineInput lin(from);

	int count = -1;

	if(lin >> count && count <= 0) {
	  //
	  std::cerr << funame << token << ": count out of range: " << count << "\n";

	  throw Error::Range();
	}

	std::string name;
      
	while(count) {
	  //
	  //std::cerr << "count = " << count << "\n";
	
	  lin.read_line(from);
	  
	  if(!(lin >> name)) // empty line
	    //
	    continue;

	  if(IO::end_key() == name)
	    //
	    break;
	
	  while(lin >> stemp) {
	    //
	    // use a separator in the case of range
	    //
	    const char sep [] = ":;,&%$@#";

	    std::string::size_type pos = stemp.find_first_of(sep);

	    if(pos != std::string::npos) {
	      //
	      double start = (double)IO::String(stemp.substr(0, pos));

	      stemp.erase(0, pos + 1);

	      pos = stemp.find_first_of(sep);

	      if(pos == std::string::npos) {
		//
		std::cerr << funame << token << ": did not find second separator: " << sep <<
		  //
		  ": available separators: <" << sep << ">: the grammar is start<sep>step<sep>finish\n";

		throw Error::Input();
	      }
	    
	      double stride = (double)IO::String(stemp.substr(0, pos));

	      if(stride <= 0.) {
		//
		std::cerr << funame << token << ": stride should be positive: " << stride << "\n";

		throw Error::Range();
	      }

	      stemp.erase(0, pos + 1);

	      double finish = (double)IO::String(stemp);

	      //std::cerr << "start, stride, finish = " << start << ", " << stride << ", " << finish << "\n";

	      for(double i = start; i <= finish; i += stride) {
		//
		dtemp = i;
		//std::cerr << "hot energy = " << dtemp << "\n";
	      
		if(hot_kcal_key == token) {
		  //
		  dtemp *= Phys_const::kcal;
		}
		else if(hot_incm_key == token) {
		  //
		  dtemp *= Phys_const::incm;
		}
		else if(hot_kj_key == token) {
		  //
		  dtemp *= Phys_const::kjoul;
		}
	  
		MasterEquation::hot_energy[name].insert(dtemp);
	      }
	    }
	    else {
	      //
	      dtemp = (double)IO::String(stemp);
	    
	      if(hot_kcal_key == token) {
		//
		dtemp *= Phys_const::kcal;
	      }
	      else if(hot_incm_key == token) {
		//
		dtemp *= Phys_const::incm;
	      }
	      else if(hot_kj_key == token) {
		//
		dtemp *= Phys_const::kjoul;
	      }
	  
	      MasterEquation::hot_energy[name].insert(dtemp);
	    }
	  }
	
	  if(!MasterEquation::hot_energy[name].size()) {
	    //
	    std::cerr << funame << token << ": " <<  name << ": no energies provided\n";
	  
	    throw Error::Init();
	  }

	  count--;
	}
      }
      // microscopic rates file
      //
      else if(mic_rout_key == token) {
	//
	if(!(from >> micro_rate_file)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      }
      // microscopic rate energy maximum
      //
      else if(mic_emax_key == token) {
	//
	if(!(from >> micro_ener_max)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}

	micro_ener_max *= Phys_const::kcal;
      }
      // microscopic rate energy minimum
      //
      else if(mic_emin_key == token) {
	//
	if(!(from >> micro_ener_min)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}

	micro_ener_min *= Phys_const::kcal;
      }
      // microscopic rate energy step
      else if(mic_step_key == token) {
	//
	if(!(from >> micro_ener_step)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}

	micro_ener_step *= Phys_const::kcal;
      }
      // time evolution
      //
      else if(tim_evol_key == token) {
	//
	if(Model::isinit()) {
	  //
	  std::cerr << funame << token << ": should be initialized before model initialization\n";

	  throw Error::Init();
	}
      
	if(Model::time_evolution) {
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	  throw Error::Input();
	}
      
	Model::time_evolution.init(new Model::TimeEvolution(from));
      }
      // state landscape output
      //
      else if(sl_key == token) {
	//
	if(state_landscape.size()) {
	  //
	  std::cerr << funame << token << ": already initialized\n";
	
	  throw Error::Init();
	}
      
	IO::LineInput lin(from);
      
	if(!(lin >> state_landscape)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";
	
	  throw Error::Input();
	}
      }
      // unknown keyword
      //
      else if(IO::skip_comment(token, from)) {
	//
	std::cerr << funame << "unknown keyword " << token << "\n";
      
	Key::show_all(std::cerr);
      
	std::cerr << "\n";
      
	throw Error::Init();
      }
    }

    if(prec_num)
      //
      MasterEquation::set_precision(prec_num);
  
    IO::log << IO::log_offset << "numerical precision: digits10 = " << MasterEquation::get_precision() << "\n";

    /*************************** CHECKING ******************************************/

    if(!from) {
      std::cerr << funame << "corrupted\n";
      throw Error::Input();
    }

    if(!Model::isinit()) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    if(!temperature.size()) {
      std::cerr << funame << "temperature list has not been initialized\n";
      throw Error::Init();
    }

    if(!pressure.size()) {
      std::cerr << funame << "pressure list has not been initialized\n";
      throw Error::Init();
    }

    if(MasterEquation::energy_step_over_temperature < 0.) {
      //
      std::cerr << funame << "energy step over temperature has not been initialized\n";
    
      throw Error::Init();
    }

    if(xtot < 0.)
      //
      IO::log << IO::log_offset << "WARNING: reference energy has not been initialized: using the default\n";

    if(reduction_scheme.size())
      //
      MasterEquation::set_default_partition(reduction_scheme);

    if(ped_spec.size())
      //
      MasterEquation::set_ped_pair(ped_spec);

    if(MasterEquation::ped_out.is_open()) {
      //
      MasterEquation::ped_out << std::setprecision(Model::ped_precision);
    
      Model::names_translation(MasterEquation::ped_out);
    }
    
    /************************** MICROSCOPIC RATE COEFFICIENTS **********************************/

    if(!IO::mpi_rank) {
      //
      if(micro_rate_file.size()) {
	//
	if(micro_ener_max <= micro_ener_min || micro_ener_step <= 0.) {
	  //
	  std::cerr << funame << "microscopic rate output: out of range\n";
      
	  throw Error::Range();
	}

	std::ofstream micro_out(micro_rate_file.c_str());
    
	if(!micro_out.is_open()) {
	  //
	  std::cerr << funame << "microscopic rate output: cannot open " << micro_rate_file << " file\n";
      
	  throw Error::Open();
	}

	Model::names_translation(micro_out);
    
	// well cycle
	//
	for(int w = 0; w < Model::well_size(); ++w) {
	  //
	  micro_out << std::setw(15) << "E, kcal/mol"
		    << std::setw(15) << "D, mol/kcal";

	  for(int b = 0; b < Model::inner_barrier_size(); ++b) {
	    //
	    const int w1 = Model::inner_connect(b).first;
	
	    const int w2 = Model::inner_connect(b).second;
	
	    if(w1 == w) {
	      //
	      stemp = Model::well(w).short_name() + "->" + Model::well(w2).short_name();
	  
	      micro_out << std::setw(15) << stemp;
	    }
	    else if(w2 == w) {
	      //
	      stemp = Model::well(w).short_name() + "->" + Model::well(w1).short_name();
	  
	      micro_out << std::setw(15) << stemp;
	    }
	  }

	  for(int b = 0; b < Model::outer_barrier_size(); ++b)
	    //
	    if(Model::outer_connect(b).first == w) {
	      //
	      const int p = Model::outer_connect(b).second;
	  
	      stemp = Model::well(w).short_name() + "->" + Model::bimolecular(p).short_name();
	  
	      micro_out << std::setw(15) << stemp;
	    }
      
	  micro_out << "\n";

	  // energy cycle
	  //
	  for(double ener = micro_ener_min; ener <= micro_ener_max; ener += micro_ener_step) {
	    //
	    double states = Model::well(w).states(ener);
	
	    micro_out << std::setw(15) << ener   / Phys_const::kcal
		      << std::setw(15) << states * Phys_const::kcal;

	    for(int b = 0; b < Model::inner_barrier_size(); ++b)
	      //
	      if(Model::inner_connect(b).first == w || Model::inner_connect(b).second == w) {
		//
		if(states != 0.) {
		  //
		  micro_out << std::setw(15) << Model::inner_barrier(b).states(ener) / states / 2. / M_PI / Phys_const::herz;
		}
		else {
		  //
		  micro_out << std::setw(15) << "***";
		}
	      }

	    for(int b = 0; b < Model::outer_barrier_size(); ++b)
	      //
	      if(Model::outer_connect(b).first == w) {
		//
		if(states != 0.) {
		  //
		  micro_out << std::setw(15) << Model::outer_barrier(b).states(ener) / states / 2. / M_PI / Phys_const::herz;
		}
		else {
		  //
		  micro_out << std::setw(15) << "***";
		}
	      }
	
	    micro_out << "\n";
	    //
	  }// energy cycle
	  //
	}// well cycle

	// bimolecular density of states
	//
	itemp = 0;
    
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  if(!Model::bimolecular(p).dummy())
	    //
	    for(int f = 0; f < 2; ++f)
	      //
	      if(Model::bimolecular(p).mode(f) == Model::DENSITY)
		//
		itemp = 1;

	if(itemp) {
	  //
	  micro_out << "Bimolecular fragments density of states, mol/kcal:\n";

	  micro_out << std::setw(15) << "E, kcal/mol";

	  std::vector<int> name_size(Model::bimolecular_size());
      
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    if(!Model::bimolecular(p).dummy())
	      //
	      for(int f = 0; f < 2; ++f)
		//
		if(Model::bimolecular(p).mode(f) == Model::DENSITY) {
		  //
		  itemp = Model::bimolecular(p).short_name().size() + 3;

		  name_size[p] = itemp > 15 ? itemp : 15;
	      
		  micro_out << std::setw(name_size[p]) << Model::bimolecular(p).short_name() + "_" + IO::String(f);
		}
      
	  micro_out << "\n";

	  for(double ener = micro_ener_min; ener <= micro_ener_max; ener += micro_ener_step) {
	    //
	    micro_out << std::setw(15) << ener / Phys_const::kcal;
      
	    for(int p = 0; p < Model::bimolecular_size(); ++p)
	      //
	      if(!Model::bimolecular(p).dummy())
		//
		for(int f = 0; f < 2; ++f)
		  //
		  if(Model::bimolecular(p).mode(f) == Model::DENSITY)
		    //
		    micro_out << std::setw(name_size[p]) << Model::bimolecular(p).states(f, ener) * Phys_const::kcal;

	    micro_out << "\n";
	  }//
	  //
	}//
	//
      }// micro output

      Model::pes_print();

      Model::pf_print();
    
    }// master node
  }
  // handling error
  //
  catch(Error::General) {
    //
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Finalize();

    return 0;
  }
  
  if(Model::no_run()) {
    //
    MPI_Finalize();
    
    return 0;
  }

#ifdef OFFLOAD_CUDA

  Offload::Cuda::Init cuda_init(IO::mpi_rank);

#else
    
  if(gpu_size && IO::mpi_rank) {
    //
    std::ostringstream gpu;

    gpu << IO::mpi_rank % gpu_size;
    
    setenv("CUDA_VISIBLE_DEVICES", gpu.str().c_str(), 1);
  }
  
#endif    
  
  std::map<std::pair<int, int>, double> rate_data, hp_rate_data;
  std::map<int, double>                            capture_data;
  MasterEquation::Partition                      partition_data;

  int job;

  MPI_Status stat;

  IO::log.clear();

  IO::aux.clear();
  
  // working node
  //
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

      dtemp = temperature[job / pressure.size()];
      
      MasterEquation::set_temperature(dtemp);
	
      // reference energy
      //
      itemp = 0;
      
      if(xtot > 0.) {
	//
	dtemp = std::floor((MasterEquation::temperature() * xtot + Model::maximum_barrier_height())
			   //
			   / Phys_const::incm + 0.5) * Phys_const::incm;

	MasterEquation::set_energy_reference(dtemp);
      }
      else
	//
	itemp = MasterEquation::DEFAULT_EREF;

      // set barriers, wells, and bimolecular species
      //
      hp_rate_data.clear();

      capture_data.clear();
      
      MasterEquation::set(hp_rate_data, capture_data, itemp);

      // pressure dependent rate coefficients
      //
      MasterEquation::set_pressure(pressure[job % pressure.size()]);
	
      // rate calculation
      //
      rate_data.clear();
      
      partition_data.clear();
      
      method(rate_data, partition_data, 0);

      // send rate, capture, and partition data to master process
      //
      MPI_Send(&job, 1, MPI_INT, 0, Comm::DATA_TAG, MPI_COMM_WORLD);
      
      Comm::send_rate_data(hp_rate_data);
      
      Comm::send_rate_data(rate_data);
      
      Comm::send_capture_data(capture_data);
      
      Comm::send_partition_data(partition_data);

      Comm::send_log();
    }
  }
  // master node
  //
  else {
    //
    Model::names_translation(IO::out);
  
    Model::names_translation(IO::log);
  
  
    /*********** PRESSURE AND TEMPERATURE DEPENDENT RATE COEFFICIENTS CALCULATION*************/

    if(MasterEquation::eval_out.is_open()) {
      //
      MasterEquation::eval_out << "*F - collisional frequency\n"
			       << "*Q - minimal relaxational eigenvalue\n"
			       << "*E - eigenvalue\n"
			       << "*P - eigenvector projection squared on the relaxational subspace\n"
			       << std::setw(13) << "Temperature,"
			       << std::setw(13) << "Pressure,"
			       << std::setw(13) << "*F,"
			       << std::setw(13) << "*Q/F";

      int eval_max = Model::well_size() + MasterEquation::evec_out_num;
    
      for(int l = 0; l < eval_max; ++l)
	//
	MasterEquation::eval_out << std::setw(13) << "*E/F" << std::setw(13) << "*P";

    
      // units
      //
      MasterEquation::eval_out << "\n"
			       << std::setw(13) << "K"
			       << std::setw(13);

      switch(MasterEquation::pressure_unit) {
      case MasterEquation::BAR:
	MasterEquation::eval_out << "bar";
	break;
      case MasterEquation::TORR:
	MasterEquation::eval_out << "torr";
	break;
      case MasterEquation::ATM:
	MasterEquation::eval_out << "atm";
	break;
      }

      MasterEquation::eval_out << std::setw(13) << "1/sec"
			       << std::setw(13) << "  ";
      for(int l = 0; l < eval_max; ++l)
	MasterEquation::eval_out << std::setw(13) << l
				 << std::setw(13) << l;

      MasterEquation::eval_out << "\n";

    }

    IO::Marker rate_marker("rate calculation");

    const int job_size = temperature.size() * pressure.size();

    IO::log << IO::log_offset << "nodes # = " << mpi_size << "\n";
    
    IO::log << IO::log_offset << "jobs  # = " << job_size << "\n";
    
    std::vector<MasterEquation::Partition>            well_partition(job_size);
    std::vector<std::map<std::pair<int, int>, double> >    rate_coef(job_size);
    std::vector<std::map<std::pair<int, int>, double> > hp_rate_coef(temperature.size());
    std::vector<std::map<int, double> >                      capture(temperature.size());

    int running_nodes = mpi_size - 1;

    if(job_size < running_nodes) {
      //
      IO::log << IO::log_offset << "WARNING: number of jobs (" << job_size
	//
	      << ") is less than the number of working nodes (" << running_nodes << ")\n";
      
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
      MPI_Recv(&job, 1, MPI_INT, MPI_ANY_SOURCE, Comm::DATA_TAG, MPI_COMM_WORLD, &stat);
      
      const int node = stat.MPI_SOURCE;

      // high pressure rate coefficients
      //
      Comm::recv_rate_data(node, hp_rate_data);

      hp_rate_coef[job / pressure.size()] = hp_rate_data;

      // phenomenological rate coefficients
      //
      Comm::recv_rate_data(node, rate_data);

      rate_coef[job] = rate_data;

      // capture rate constants
      //
      Comm::recv_capture_data(node, capture_data);

      capture[job / pressure.size()] = capture_data;

      // well partitioning
      //
      Comm::recv_partition_data(node, partition_data);

      well_partition[job] = partition_data;

      // receve log
      //
      Comm::recv_log(node);
      
      IO::log << "\n" << IO::log_offset << job << "-th job is done (node = " << node << ")\n\n";
      
      // remove working node
      //
      if(new_job == job_size) {
	//
	--running_nodes;
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

    /***********************************************************
     ************************* OUTPUT **************************
     ***********************************************************/
    //
    std::vector<std::string> spec_name;
  
    for(int w = 0; w < Model::well_size(); ++w)
      //
      spec_name.push_back(Model::well(w).short_name());
  
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      spec_name.push_back(Model::bimolecular(p).short_name());

    for(int i = 0; i < spec_name.size(); ++i)
      //
      if(!i || spec_name[i].size() > itemp)
	//
	itemp = spec_name[i].size();

    const int max_spec_size = itemp > 8 ? itemp : 8;

    std::vector<int> spec_field(spec_name.size());

    for(int i = 0; i < spec_name.size(); ++i)
      //
      spec_field[i] = Model::out_precision + 7 > spec_name[i].size() + 1 ? Model::out_precision + 7 : spec_name[i].size() + 1;

    std::vector<int> escape_field(Model::escape_size());
  
    for(int e = 0; e < Model::escape_size(); ++e)
      //
      escape_field[e] = Model::out_precision + 7 > Model::escape_name(e).size() + 1 ? Model::out_precision + 7 : Model::escape_name(e).size() + 1;
	
    IO::out << "Unimolecular Rate Units: 1/sec;  Bimolecular Rate Units: cm^3/sec\n\n"
	    << "______________________________________________________________________________________\n\n"
	    << "Species-Species Rate Tables:\n\n";


    // Output
    //
    IO::out << "______________________________________________________________________________________\n\n"
	    << "High Pressure Rate Coefficients (Temperature-Species Rate Tables):\n\n";

    for(int i = 0; i < spec_name.size(); ++i) {
      //
      IO::out << "Reactant = " << spec_name[i] << "\n";
    
      IO::out << std::setw(7) << "T(K)";
    
      for(int j = 0; j < spec_name.size(); ++j)
	//
	IO::out << std::setw(spec_field[j]) << spec_name[j];
    
      IO::out << "\n";
    
      for(int t = 0; t < temperature.size(); ++t) {
	//
	IO::out << std::setw(7) << temperature[t] / Phys_const::kelv;
      
	for(int j = 0; j < spec_name.size(); ++j)
	  //
	  if(hp_rate_coef[t].find(std::make_pair(i, j)) != hp_rate_coef[t].end())
	    //
	    IO::out << std::setw(spec_field[j]) << hp_rate_coef[t][std::make_pair(i, j)];
	  else
	    //
	    IO::out << std::setw(spec_field[j]) << "***";
      
	IO::out << "\n";
      }
      IO::out << "\n";
    }
  
    IO::out << "Capture/Escape Rate Coefficients:\n\n";
  
    IO::out << std::setw(7) << "T(K)";
  
    for(int i = 0; i < spec_name.size(); ++i)
      //
      IO::out << std::setw(spec_field[i]) << spec_name[i];
  
    IO::out << "\n";
  
    for(int t = 0; t < temperature.size(); ++t) {
    
      IO::out << std::setw(7) << temperature[t] / Phys_const::kelv;
    
      for(int i = 0; i < spec_name.size(); ++i)
	//
	if(capture[t].find(i) != capture[t].end())
	  //
	  IO::out << std::setw(spec_field[i]) << capture[t][i];
	else
	  //
	  IO::out << std::setw(spec_field[i]) << "***";
    
      IO::out << "\n";
    }
    IO::out << "\n";

    IO::out << "______________________________________________________________________________________\n\n"
	    << "Pressure-Species Rate Tables:\n\n";

    // pressure dependence
    //
    for(int i = 0; i < spec_name.size(); ++i) {
      //
      for(int t = 0; t < temperature.size(); ++t) {
	// title
	IO::out << "Reactant = " << spec_name[i] << "   Temperature = " << temperature[t] / Phys_const::kelv << " K\n\n";
      
	IO::out << std::setw(9);
	switch(MasterEquation::pressure_unit) {
	case MasterEquation::BAR:
	  IO::out << "P(bar)";
	  break;
	case MasterEquation::TORR:
	  IO::out << "P(torr)";
	  break;
	case MasterEquation::ATM:
	  IO::out << "P(atm)";
	  break;
	}
      
	for(int j = 0; j < spec_name.size(); ++j)
	  //
	  if(j != i)
	    //
	    IO::out << std::setw(spec_field[j]) << spec_name[j];

	IO::out << std::setw(Model::out_precision + 7) << "Loss";

	// escape
	//
	for(int e = 0; e < Model::escape_size(); ++e)
	  //
	  IO::out << std::setw(escape_field[e]) << Model::escape_name(e);
      
	IO::out << "\n";

	// rates
	//
	for(int p = 0; p < pressure.size(); ++p) {
	  IO::out << std::setw(9);
	  switch(MasterEquation::pressure_unit) {
	  case MasterEquation::BAR:
	    IO::out << pressure[p]    / Phys_const::bar;
	    break;
	  case MasterEquation::TORR:
	    IO::out << pressure[p]    / Phys_const::tor;
	    break;
	  case MasterEquation::ATM:
	    IO::out << pressure[p]    / Phys_const::atm;
	    break;
	  }
	  itemp = p + t * pressure.size();
	
	  // to products
	  //
	  for(int j = 0; j < spec_name.size(); ++j)
	    if(j != i) {
	      //
	      if(rate_coef[itemp].find(std::make_pair(i, j)) != rate_coef[itemp].end())
		//
		IO::out << std::setw(spec_field[j]) << rate_coef[itemp][std::make_pair(i, j)];
	      else
		//
		IO::out << std::setw(spec_field[j]) << "***";
	    }
	  // reactant loss
	  //
	  if(rate_coef[itemp].find(std::make_pair(i, i)) != rate_coef[itemp].end())
	    //
	    IO::out << std::setw(Model::out_precision + 7) << rate_coef[itemp][std::make_pair(i, i)];
	  else
	    //
	    IO::out << std::setw(Model::out_precision + 7) << "***";
	
	  // escape rates
	  //
	  for(int e = 0; e < Model::escape_size(); ++e) {
	    //
	    ptemp = std::make_pair(i, Model::well_size() + Model::bimolecular_size() + e);
	  
	    if(rate_coef[itemp].find(ptemp) != rate_coef[itemp].end()) {
	      //
	      IO::out << std::setw(escape_field[e]) << rate_coef[itemp][ptemp];
	    }
	    else
	      //
	      IO::out << std::setw(escape_field[e]) << "***";
	  }
	  IO::out << "\n";
	}
      
	IO::out << std::setw(9) << "O-O";
      
	for(int j = 0; j < spec_name.size(); ++j)
	  //
	  if(j != i) {
	    if(hp_rate_coef[t].find(std::make_pair(i, j)) != hp_rate_coef[t].end())
	      IO::out << std::setw(spec_field[j]) << hp_rate_coef[t][std::make_pair(i, j)];
	    else
	      IO::out << std::setw(spec_field[j]) << "***";
	  }
      
	if(capture[t].find(i) != capture[t].end())
	  //
	  IO::out << std::setw(Model::out_precision + 7) << capture[t][i];
	else
	  //
	  IO::out << std::setw(Model::out_precision + 7) << "***";
	IO::out << "\n\n";
      }
    }

    IO::out << "______________________________________________________________________________________\n\n"
	    << "Temperature-Species Rate Tables:\n\n";

    // temperature dependence
    //
    for(int i = 0; i < spec_name.size(); ++i) {
      //
      for(int p = 0; p < pressure.size(); ++p) {
	//
	IO::out << "Reactant = " << spec_name[i] << "   Pressure = ";
      
	switch(MasterEquation::pressure_unit) {
	  //
	case MasterEquation::BAR:
	  //
	  IO::out << pressure[p] / Phys_const::bar << " bar\n\n";
	
	  break;
	case MasterEquation::TORR:
	  //
	  IO::out << pressure[p] / Phys_const::tor << " torr\n\n";
	
	  break;
	case MasterEquation::ATM:
	  //
	  IO::out << pressure[p] / Phys_const::atm << " atm\n\n";
	
	  break;
	}
	IO::out << std::setw(7) << "T(K)";
	//
	for(int j = 0; j < spec_name.size(); ++j)
	  //
	  if(j != i)
	    //
	    IO::out << std::setw(spec_field[j]) << spec_name[j];
      
	IO::out << std::setw(Model::out_precision + 7) << "Loss"
	  //
		<< std::setw(Model::out_precision + 7) << "Capture";

	// escape
	//
	for(int e = 0; e < Model::escape_size(); ++e)
	  //
	  IO::out << std::setw(escape_field[e]) << Model::escape_name(e);
      
	IO::out << "\n";

	// rates
	//
	for(int t = 0; t < temperature.size(); ++t) {
	  //
	  IO::out << std::setw(7) << temperature[t] / Phys_const::kelv;
	
	  itemp = p + t * pressure.size();
	
	  // to products
	  //
	  for(int j = 0; j < spec_name.size(); ++j)
	    //
	    if(j != i) {
	      //
	      if(rate_coef[itemp].find(std::make_pair(i, j)) != rate_coef[itemp].end())
		//
		IO::out << std::setw(spec_field[j]) << rate_coef[itemp][std::make_pair(i, j)];
	      else
		//
		IO::out << std::setw(spec_field[j]) << "***";
	    }
	
	  // reactant loss
	  //
	  if(rate_coef[itemp].find(std::make_pair(i, i)) != rate_coef[itemp].end())
	    //
	    IO::out << std::setw(Model::out_precision + 7) << rate_coef[itemp][std::make_pair(i, i)];
	  else
	    //
	    IO::out << std::setw(Model::out_precision + 7) << "***";
	
	  // capture
	  //
	  if(capture[t].find(i) != capture[t].end())
	    //
	    IO::out << std::setw(Model::out_precision + 7) << capture[t][i];
	  else
	    //
	    IO::out << std::setw(Model::out_precision + 7) << "***";

	  // escape rates
	  //
	  for(int e = 0; e < Model::escape_size(); ++e) {
	    //
	    ptemp = std::make_pair(i, Model::well_size() + Model::bimolecular_size() + e);
	  
	    if(rate_coef[itemp].find(ptemp) != rate_coef[itemp].end()) {
	      //
	      IO::out << std::setw(escape_field[e]) << rate_coef[itemp][ptemp];
	    }
	    else
	      //
	      IO::out << std::setw(escape_field[e]) << "***";
	  }
	
	  IO::out << "\n";
	}
	IO::out << "\n";
      }
    }

    IO::out << "______________________________________________________________________________________\n\n"
	    << "Temperature-Pressure Rate Tables:\n\n";

    // T-P table
    for(int i = 0; i < spec_name.size(); ++i)
      for(int j = 0; j < spec_name.size(); ++j)
	if(j != i) {
	  // title
	  std::pair<int, int> proc(i, j);
	  stemp = spec_name[i] + "->" + spec_name[j];
	  IO::out << stemp << "\n\n";

	  IO::out << std::left << std::setw(7) << "P\\T" << std::right;
	  for(int t = 0; t < temperature.size(); ++t)
	    IO::out << std::setw(Model::out_precision + 7) << temperature[t] / Phys_const::kelv; 
	  IO::out << "\n";

	  for(int p = 0; p < pressure.size(); ++p) {
	    IO::out << std::left << std::setw(7);
	    switch(MasterEquation::pressure_unit) {
	    case MasterEquation::BAR:
	      IO::out << pressure[p] / Phys_const::bar;
	      break;
	    case MasterEquation::TORR:
	      IO::out << pressure[p] / Phys_const::tor;
	      break;
	    case MasterEquation::ATM:
	      IO::out << pressure[p] / Phys_const::atm;
	      break;
	    }
	    IO::out << std::right;
	    for(int t = 0; t < temperature.size(); ++t) {
	      itemp = p + t * pressure.size();
	      if(rate_coef[itemp].find(proc) != rate_coef[itemp].end())
		IO::out << std::setw(Model::out_precision + 7) << rate_coef[itemp][proc];
	      else
		IO::out << std::setw(Model::out_precision + 7) << "***";
	    }
	    IO::out << "\n";
	  }
	  IO::out << std::left << std::setw(7) << "O-O" << std::right;
	  for(int t = 0; t < temperature.size(); ++t)
	    if(hp_rate_coef[t].find(proc) != hp_rate_coef[t].end())
	      IO::out << std::setw(Model::out_precision + 7) << hp_rate_coef[t][proc];
	    else
	      IO::out << std::setw(Model::out_precision + 7) << "***";
	  IO::out << "\n\n";
	}

    // loss
    //
    for(int i = 0; i < spec_name.size(); ++i) {
      //
      IO::out << spec_name[i] + "-> " << "\n\n";

      IO::out << std::left << std::setw(7) << "P\\T" << std::right;
      for(int t = 0; t < temperature.size(); ++t)
	IO::out << std::setw(Model::out_precision + 7) << temperature[t] / Phys_const::kelv; 
      IO::out << "\n";

      for(int p = 0; p < pressure.size(); ++p) {
	IO::out << std::left << std::setw(7);
	switch(MasterEquation::pressure_unit) {
	case MasterEquation::BAR:
	  IO::out << pressure[p] / Phys_const::bar;
	  break;
	case MasterEquation::TORR:
	  IO::out << pressure[p] / Phys_const::tor;
	  break;
	case MasterEquation::ATM:
	  IO::out << pressure[p] / Phys_const::atm;
	  break;
	}
	IO::out << std::right;
	for(int t = 0; t < temperature.size(); ++t) {
	  itemp = p + t * pressure.size();
	  if(rate_coef[itemp].find(std::make_pair(i, i)) != rate_coef[itemp].end())
	    IO::out << std::setw(Model::out_precision + 7) << rate_coef[itemp][std::make_pair(i, i)];
	  else
	    IO::out << std::setw(Model::out_precision + 7) << "***";
	}
	IO::out << "\n";
      }
      IO::out << std::left << std::setw(7) << "O-O" << std::right;
      
      for(int t = 0; t < temperature.size(); ++t)
	//
	if(capture[t].find(i) != capture[t].end()) {
	  //
	  IO::out << std::setw(Model::out_precision + 7) << capture[t][i];
	}
	else
	  //
	  IO::out << std::setw(Model::out_precision + 7) << "***";
      
      IO::out << "\n\n";
      //
    }// species cycle
    
    IO::log << IO::log_offset << "all jobs are done => exitting\n\n";

    for(int i = 1; i < mpi_size; ++i)
      //
      MPI_Send(&new_job, 1, MPI_INT, i, Comm::STOP_TAG, MPI_COMM_WORLD);

  }//master node

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  
  return 0;
}
