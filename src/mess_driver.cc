/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>

#include "mess.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "master_equation: ";

  if (argc < 2) {
    std::cout << "usage: mess input_file\n";
    return 0;
  }

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  std::pair<int, int> ptemp;

  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key    model_key("Model"                      );
  Key bar_pres_key("PressureList[bar]"          );
  Key tor_pres_key("PressureList[torr]"         );
  Key atm_pres_key("PressureList[atm]"          );
  Key     temp_key("TemperatureList[K]"         );
  Key    estep_key("EnergyStep[1/cm]"           );
  Key     etot_key("EnergyStepOverTemperature"  );
  Key     xtot_key("ExcessEnergyOverTemperature");
  Key ref_incm_key("ReferenceEnergy[1/cm]"      );
  Key ref_kcal_key("ReferenceEnergy[kcal/mol]"  );
  Key   ref_kj_key("ReferenceEnergy[kJ/mol]"    );
  Key cut_kcal_key("GlobalCutoff[kcal/mol]"     );
  Key cut_incm_key("GlobalCutoff[1/cm]"         );
  Key   cut_kj_key("GlobalCutoff[kJ/mol]"       );
  Key     emax_key("ModelEnergyLimit[kcal/mol]" );
  Key  adm_bor_key("AtomDistanceMin[bohr]"      );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );
  Key     xtun_key("TunnelingActionCutoff"      );
  Key well_cut_key("WellCutoff"                 );
  Key eval_max_key("ChemicalEigenvalueMax"      );
  Key eval_min_key("ChemicalEigenvalueMin"      );
  Key well_red_key("WellReductionThreshold"     );
  Key     calc_key("CalculationMethod"          );
  Key  rat_red_key("ReductionMethod"            );
  Key      wpm_key("WellPartitionMethod"        );
  Key      wpt_key("WellProjectionThreshold"    );
  Key rate_out_key("RateOutput"                 );
  Key  log_out_key("LogOutput"                  );
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

  std::vector<std::string> ped_spec;// product energy distribution pairs verbal
  std::vector<std::string> reduction_scheme;
  std::vector<double> pressure;
  std::vector<double> temperature;
  double estep  = -1.; // energy step
  double etot   = -1.; // energy step over temperature
  double xtot   = -1.; // exsess energy over temperature
  double eref;
  bool iseref = false;
  MasterEquation::Method method = MasterEquation::direct_diagonalization_method;
  std::string micro_rate_file;
  double micro_ener_max  = 0.;
  double micro_ener_min  = 0.;
  double micro_ener_step = -1.;
  std::string state_landscape;

  // base name
  std::string base_name = argv[1];
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    base_name.resize(base_name.size() - 4);

  IO::KeyBufferStream from(argv[1]);
  if(!from && base_name == argv[1]) {
    // try inp extension
    stemp = base_name + ".inp";
    from.open(stemp.c_str());
  }

  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string token, comment, line;
  while(from >> token) {
    // main input group 
    if(model_key == token) {
      // default log output	
      if(!IO::log.is_open()) {
	stemp = base_name + ".log";
	IO::log.open(stemp.c_str());
	if(!IO::log) {
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	  throw Error::Input();
	}
      }

      // default rate output
      if(!IO::out.is_open()) {
	stemp = base_name + ".out";
	IO::out.open(stemp.c_str());
	if(!IO::out) {
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	  throw Error::Input();
	}
      }

      // global energy limit check
      if(!Model::is_energy_limit()) {
	std::cerr << funame << "model energy limit should be set BEFORE model input section\n";
	throw Error::Init();
      }

      // main initialization
      try {
	Model::init(from);
      }  
      catch(Error::General) {
	IO::log << std::flush;
	throw;
      }
      break;
    }
    // rate output
    else if(rate_out_key == token) {
      if(IO::out.is_open()) {
        std::cerr << funame << token << ": allready opened\n";
        throw Error::Input();
      }      
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": bad input\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      IO::out.open(stemp.c_str());
      if(!IO::out) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // log output
    else if(log_out_key == token) {
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
    // eigenvalue output
    else if(eval_out_key == token) {
      if(MasterEquation::eval_out.is_open()) {
        std::cerr << funame << token << ": allready opened\n";
        throw Error::Init();
      }      
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      MasterEquation::eval_out.open(stemp.c_str());
      if(!MasterEquation::eval_out) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // eigenvector output
    else if(evec_out_key == token) {
      if(MasterEquation::evec_out.is_open()) {
        std::cerr << funame << token << ": allready opened\n";
        throw Error::Init();
      }      
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      MasterEquation::evec_out.open(stemp.c_str());
      if(!MasterEquation::evec_out) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // product energy distribution output
    else if(ped_out_key == token) {
      if(MasterEquation::ped_out.is_open()) {
	std::cerr << funame << token << ": already opened\n";
	throw Error::Open();
      }
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      MasterEquation::ped_out.open(stemp.c_str());
      if(!MasterEquation::ped_out.is_open()) {
	std::cerr << funame << token << ": cannot open the " << stemp << " file\n";
	throw Error::Open();
      }
    }
    // reactants and products for product energy distribution output
    else if(ped_spec_key == token) {
      IO::LineInput ped_input(from);
      while(ped_input >> stemp)
	ped_spec.push_back(stemp);
    }
    // bimolecular reactant to be used as a reference
    else if(react_key == token) {
      if(!(from >> Model::reactant)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);      
    }
    // eigenvector number to output
    else if(evec_num_key == token) {
      if(!(from >> MasterEquation::evec_out_num)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(MasterEquation::evec_out_num < 0) {
        std::cerr << funame << token << ": out of range\n";
        throw Error::Range();
      }
    }
    // temperature
    else if(temp_key == token) {
      if(temperature.size()) {
	std::cerr << funame << token <<  ": already initialized\n";
	throw Error::Init();
      }
      IO::LineInput data_input(from);
      std::set<double> data;
      while(data_input >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}
	data.insert( dtemp * Phys_const::kelv);
      }
      
      if(!data.size()) {
        std::cerr << funame << token << ": no data\n";
        throw Error::Init();
      }

      for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it)
	temperature.push_back(*it);
    }
    // pressure
    else if(bar_pres_key == token || tor_pres_key == token || atm_pres_key == token) {
      if(pressure.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      IO::LineInput data_input(from);
      std::set<double> data;

      while(data_input >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}

	if(bar_pres_key == token) {
	  dtemp *= Phys_const::bar;
	  MasterEquation::pressure_unit = MasterEquation::BAR;
	}
	if(atm_pres_key == token) {
	  dtemp *= Phys_const::atm;
	  MasterEquation::pressure_unit = MasterEquation::ATM;
	}
	if(tor_pres_key == token) {
	  dtemp *= Phys_const::tor;
	  MasterEquation::pressure_unit = MasterEquation::TORR;
	}

	data.insert(dtemp);
      }

      if(!data.size()) {
        std::cerr << funame << token << ": no data\n";
        throw Error::Init();
      }

      for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it)
	pressure.push_back(*it);
    }
    // energy step
    else if(estep_key == token) {
      if(estep > 0. || etot > 0.) {
	std::cerr << funame << "energy step has been already defined\n";
	throw Error::Input();
      }
      if(!(from >> estep)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(estep <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      estep *= Phys_const::incm;
    }
    // energy step over temperature
    else if(etot_key == token) {
      if(estep > 0. || etot > 0.) {
	std::cerr << funame << "energy step has been already defined\n";
	throw Error::Input();
      }
      if(!(from >> etot)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(etot <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
    }
    // excess energy over temperature
    else if(xtot_key == token) {
      if(xtot > 0. || iseref) {
	std::cerr << funame << "reference energy has been already defined\n";
	throw Error::Input();
      }
      if(!(from >> xtot)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(xtot <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
    }
    // reference energy
    else if(ref_incm_key == token || ref_kcal_key == token || ref_kj_key == token) {
      if(xtot > 0. || iseref) {
	std::cerr << funame << "reference energy has been already defined\n";
	throw Error::Input();
      }
      iseref = true;
      if(!(from >> eref)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);
      
      if(ref_incm_key == token)
	eref *= Phys_const::incm;
      if(ref_kcal_key == token)
	eref *= Phys_const::kcal;
      if(ref_kj_key == token)
	eref *= Phys_const::kjoul;
    }
    // global energy cutoff
    else if(cut_incm_key == token || cut_kcal_key == token || cut_kj_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);
      
      if(cut_incm_key == token)
	dtemp *= Phys_const::incm;
      if(cut_kcal_key == token)
	dtemp *= Phys_const::kcal;
      if(cut_kj_key == token)
	dtemp *= Phys_const::kjoul;

      MasterEquation::set_global_cutoff(dtemp);
    }
    // model energy limit
    else if(emax_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      Model::set_energy_limit(dtemp * Phys_const::kcal);
    }
    // minimal interatomic distance
    else if(adm_bor_key == token || adm_ang_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0) {
        std::cerr << funame << token << ": out of range\n";
        throw Error::Range();
      }
      
      if(adm_ang_key == token)
	dtemp *= Phys_const::angstrom;

      Model::atom_dist_min = dtemp;
    }
    // maximum tunneling exponent
    else if(xtun_key == token) {
      if(!(from >> dtemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0.) {
        std::cerr << funame << token << ": should be positive\n";
        throw Error::Range();
      }

      Model::Tunnel::set_action_max(dtemp);
    }
    // calculation method
    else if(calc_key == token) {
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
    else if(well_cut_key == token) {
      if(!(from >> MasterEquation::well_cutoff)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(MasterEquation::well_cutoff <= 0.) {
        std::cerr << funame << token << ": should be positive\n";
        throw Error::Range();
      }
    }
    // chemical relaxation to collision frequency ratio
    else if(well_red_key == token) {
      if(!(from >> MasterEquation::reduction_threshold)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);
    }
    // maximal chemical eigenvalue
    else if(eval_max_key == token) {
      if(!(from >> MasterEquation::chemical_threshold)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);
    }
    // minumal chemical eigenvalue
    else if(eval_min_key == token) {
      if(!(from >> MasterEquation::min_chem_eval)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      if(MasterEquation::min_chem_eval <= 0. || MasterEquation::min_chem_eval >= 1.) {
        std::cerr << funame << token << ": out of range\n";
        throw Error::Range();
      }
    }
    // number of closest reductions to print
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
    else if(def_chem_key == token) {
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      MasterEquation::set_default_chem_size(itemp);
    }
    // hot energies 
    else if(hot_kcal_key == token || hot_incm_key == token || hot_kj_key == token) {
      int line_size;
      if(!(from >> line_size)) {
	std::cerr << funame << token <<  ": cannot read number of input lines\n";
	throw Error::Input();
      }

      if(line_size <= 0) {
	std::cerr << funame << token << ": negative number of input lines\n";
	throw Error::Range();
      }
      std::getline(from, comment);

      for(int l = 0; l < line_size; ++l) {
	IO::LineInput line_input(from);
	std::string well_name;
	if(!(line_input >> well_name)) {
	  std::cerr << funame << token << ": " <<  l << "-th input line is corrupted\n";
	  throw Error::Input();
	}

	while(line_input >> dtemp) {

	  if(hot_kcal_key == token)
	    dtemp *= Phys_const::kcal;
	  if(hot_incm_key == token)
	    dtemp *= Phys_const::incm;
	  if(hot_kj_key == token)
	    dtemp *= Phys_const::kjoul;

	  MasterEquation::hot_energy[well_name].push_back(dtemp);
	}

	if(!MasterEquation::hot_energy[well_name].size()) {
	  std::cerr << funame << token << ": " <<  well_name << ": no energies\n";
	  throw Error::Init();
	}
      }
    }
    // microscopic rates file
    else if(mic_rout_key == token) {
      if(!(from >> micro_rate_file)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
    }
    // microscopic rate energy maximum
    else if(mic_emax_key == token) {
      if(!(from >> micro_ener_max)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      micro_ener_max *= Phys_const::kcal;
    }
    // microscopic rate energy minimum
    else if(mic_emin_key == token) {
      if(!(from >> micro_ener_min)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      micro_ener_min *= Phys_const::kcal;
    }
    // microscopic rate energy step
    else if(mic_step_key == token) {
      if(!(from >> micro_ener_step)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      micro_ener_step *= Phys_const::kcal;
    }
    // time evolution
    else if(tim_evol_key == token) {
      if(Model::time_evolution) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Input();
      }
      Model::time_evolution.init(new Model::TimeEvolution(from));
    }
    // state landscape output
    else if(sl_key == token) {
      if(state_landscape.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      IO::LineInput lin(from);
      if(!(lin >> state_landscape)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  /*************************** CHECKING ******************************************/

  if(!from) {
    std::cerr << funame << "corrupted\n";
    throw Error::Input();
  }

  if(!Model::isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Input();
  }

  if(!temperature.size()) {
    std::cerr << funame << "temperature list has not been initialized\n";
    throw Error::Input();
  }

  if(!pressure.size()) {
    std::cerr << funame << "pressure list has not been initialized\n";
    throw Error::Input();
  }

  if(estep <= 0. && etot <= 0.) {
    std::cerr << funame << "energy step has not been initialized\n";
    throw Error::Input();
  }

  if(xtot <= 0. && !iseref) {
    std::cerr << funame << "reference energy has not been initialized\n";
    throw Error::Input();
  }

  if(reduction_scheme.size())
    MasterEquation::set_default_partition(reduction_scheme);

  if(ped_spec.size())
    MasterEquation::set_ped_pair(ped_spec);

  /************************** MICROSCOPIC RATE COEFFICIENTS **********************************/

  if(micro_rate_file.size()) {

    if(micro_ener_max <= micro_ener_min || micro_ener_step <= 0.) {
      std::cerr << funame << "microscopic rate output: out of range\n";
      throw Error::Range();
    }

    std::ofstream micro_out(micro_rate_file.c_str());
    if(!micro_out.is_open()) {
      std::cerr << funame << "microscopic rate output: cannot open " << micro_rate_file << " file\n";
      throw Error::Open();
    }

    // well cycle
    for(int w = 0; w < Model::well_size(); ++w) {
      micro_out << std::setw(15) << "E, kcal/mol"
		<< std::setw(15) << "D, mol/kcal";

      for(int b = 0; b < Model::inner_barrier_size(); ++b) {
	const int w1 = Model::inner_connect(b).first;
	const int w2 = Model::inner_connect(b).second;
	if(w1 == w) {
	  stemp = Model::well(w).name() + "->" + Model::well(w2).name();
	  micro_out << std::setw(15) << stemp;
	}
	else if(w2 == w) {
	  stemp = Model::well(w).name() + "->" + Model::well(w1).name();
	  micro_out << std::setw(15) << stemp;
	}
      }

      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	if(Model::outer_connect(b).first == w) {
	  const int p = Model::outer_connect(b).second;
	  stemp = Model::well(w).name() + "->" + Model::bimolecular(p).name();
	  micro_out << std::setw(15) << stemp;
	}
      micro_out << "\n";

      // energy cycle
      for(double ener = micro_ener_min; ener <= micro_ener_max; ener += micro_ener_step) {
	double states = Model::well(w).states(ener);
	micro_out << std::setw(15) << ener   / Phys_const::kcal
		  << std::setw(15) << states * Phys_const::kcal;

	for(int b = 0; b < Model::inner_barrier_size(); ++b)
	  if(Model::inner_connect(b).first == w || Model::inner_connect(b).second == w)
	    if(states != 0.)
	      micro_out << std::setw(15) << Model::inner_barrier(b).states(ener) / states / 2. / M_PI / Phys_const::herz;
	    else
	      micro_out << std::setw(15) << "***";

	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  if(Model::outer_connect(b).first == w)
	    if(states != 0.)
	      micro_out << std::setw(15) << Model::outer_barrier(b).states(ener) / states / 2. / M_PI / Phys_const::herz;
	    else
	      micro_out << std::setw(15) << "***";
	micro_out << "\n";
      } // energy cycle
    } // well cycle
  } // micro output

  if(Model::no_run())
    return 0;
  
  /*********** PRESSURE AND TEMPERATURE DEPENDENT RATE COEFFICIENTS CALCULATION*************/

  if(estep > 0.)
    etot = estep / *temperature.begin();

  if(iseref) {
    eref += Model::energy_shift();
    xtot = (eref - Model::maximum_barrier_height()) / *temperature.begin();
  }

  if(MasterEquation::eval_out.is_open()) {
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
      MasterEquation::eval_out << std::setw(13) << "*E/F"
			       << std::setw(13) << "*P";

    
    // units
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

  typedef std::vector<double>::const_iterator It;

  std::map<std::pair<int, int>, double> rate_data;
  std::map<int, double> capture_data;
  std::vector<MasterEquation::Partition> well_partition;
  std::vector<std::map<std::pair<int, int>, double> > rate_coef;
  std::vector<std::map<std::pair<int, int>, double> > hp_rate_coef;
  std::vector<std::map<int, double> > capture;

  std::vector<std::string> spec_name;
  for(int w = 0; w < Model::well_size(); ++w)
    spec_name.push_back(Model::well(w).name());
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    spec_name.push_back(Model::bimolecular(p).name());

  IO::out << "Unimolecular Rate Units: 1/sec;  Bimolecular Rate Units: cm^3/sec\n\n"
	  << "______________________________________________________________________________________\n\n"
	  << "Species-Species Rate Tables:\n\n";


  //  try {
  {

    IO::Marker rate_marker("rate calculation");

    for(It t = temperature.begin(); t != temperature.end(); ++t) {// temperature cycle
      MasterEquation::set_temperature(*t);

      // energy step
      if(estep > 0.)
	MasterEquation::set_energy_step(estep);
      else
	MasterEquation::set_energy_step(nearbyint(*t * etot / Phys_const::incm) * Phys_const::incm);
      
      // reference energy
      if(iseref)
	MasterEquation::set_energy_reference(eref);
      else
	MasterEquation::set_energy_reference(nearbyint((*t * xtot + Model::maximum_barrier_height())
						       / Phys_const::incm) * Phys_const::incm);

      // set barriers, wells, and bimolecular species
      MasterEquation::set(rate_data, capture_data);
      hp_rate_coef.push_back(rate_data);
      capture.push_back(capture_data);

      // output
      IO::out << "Temperature = " << *t / Phys_const::kelv  << " K\n\n";
      IO::out << "High Pressure Rate Coefficients:\n\n"
	      << std::left << std::setw(8) << "From\\To" << std::right;
      for(int j = 0; j < spec_name.size(); ++j)
	IO::out << std::setw(13) << spec_name[j];
      IO::out << "\n";
      for(int i = 0; i < spec_name.size(); ++i) {
	IO::out << std::left << std::setw(8) << spec_name[i] << std::right;
	for(int j = 0; j < spec_name.size(); ++j) {
	  ptemp = std::make_pair(i, j);
	  if(rate_data.find(ptemp) != rate_data.end())
	    IO::out << std::setw(13) << rate_data[ptemp];
	  else
	    IO::out << std::setw(13) << "***";
	}
	IO::out << "\n";
      }
      IO::out << "\n";

      // pressure dependent rate coefficients
      for(It p = pressure.begin(); p != pressure.end(); ++p) {// pressure cycle
	MasterEquation::set_pressure(*p);
	// rate calculation
	if(method) {
	  well_partition.push_back(MasterEquation::Partition());
	  method(rate_data, well_partition.back(), 0);
	  rate_coef.push_back(rate_data);
	}

	// output
	IO::out << "Temperature = " << *t / Phys_const::kelv <<  " K    Pressure = ";
	switch(MasterEquation::pressure_unit) {
	case MasterEquation::BAR:
	  IO::out << *p / Phys_const::bar << " bar";
	  break;
	case MasterEquation::TORR:
	  IO::out << *p / Phys_const::tor << " torr";
	  break;
	case MasterEquation::ATM:
	  IO::out << *p / Phys_const::atm << " atm";
	  break;
	}
	IO::out << "\n\n";
	IO::out << std::left << std::setw(8) << "From\\To" << std::right;
	for(int j = 0; j < spec_name.size(); ++j)
	  IO::out << std::setw(13) << spec_name[j];
	
	for(int w = 0; w < Model::well_size(); ++w)
	  if(Model::well(w).escape())
	    IO::out << std::setw(13) << Model::well(w).name(); 
	IO::out << "\n";

	for(int i = 0; i < spec_name.size(); ++i) {
	  IO::out << std::left << std::setw(8) << spec_name[i] << std::right;
	  for(int j = 0; j < spec_name.size(); ++j) {
	    ptemp = std::make_pair(i, j);
	    if(rate_data.find(ptemp) != rate_data.end())
	      IO::out << std::setw(13) << rate_data[ptemp];
	    else
	      IO::out << std::setw(13) << "***";
	  }
	  // escape rates
	  for(int w = 0; w < Model::well_size(); ++w)
	    if(Model::well(w).escape()) {
	      ptemp = std::make_pair(i, Model::well_size() + Model::bimolecular_size() + w);
	      if(rate_data.find(ptemp) != rate_data.end())
		IO::out << std::setw(13) << rate_data[ptemp];
	      else
		IO::out << std::setw(13) << "***";
	    }
	  IO::out << "\n";
	}
	IO::out << "\n";
      }// pressure cycle
    }// temperature cycle
  }
  //catch(Error::General) {
  // IO::log << std::flush;
  //throw;
  // }

  // Output
  IO::out << "______________________________________________________________________________________\n\n"
	  << "High Pressure Rate Coefficients (Temperature-Species Rate Tables):\n\n";

  for(int i = 0; i < spec_name.size(); ++i) {
    IO::out << std::setw(7) << "T(K)";
    for(int j = 0; j < spec_name.size(); ++j) {
      stemp = spec_name[i] + "->" + spec_name[j];
      IO::out << std::setw(13) << stemp;
    }
    IO::out << "\n";
    for(int t = 0; t < temperature.size(); ++t) {
      IO::out << std::setw(7) << temperature[t] / Phys_const::kelv; 
      for(int j = 0; j < spec_name.size(); ++j)
	if(hp_rate_coef[t].find(std::make_pair(i, j)) != hp_rate_coef[t].end())
	  IO::out << std::setw(13) << hp_rate_coef[t][std::make_pair(i, j)];
	else
	  IO::out << std::setw(13) << "***";
      IO::out << "\n";
    }
    IO::out << "\n";
  }
  
  IO::out << "Capture/Escape Rate Coefficients:\n\n";
  IO::out << std::setw(7) << "T(K)";
  for(int i = 0; i < spec_name.size(); ++i)
    IO::out << std::setw(13) << spec_name[i];
  IO::out << "\n";
  for(int t = 0; t < temperature.size(); ++t) {
    IO::out << std::setw(7) << temperature[t] / Phys_const::kelv; 
    for(int i = 0; i < spec_name.size(); ++i)
      if(capture[t].find(i) != capture[t].end())
	IO::out << std::setw(13) << capture[t][i];
      else
	IO::out << std::setw(13) << "***";
    IO::out << "\n";
  }
  IO::out << "\n";

  IO::out << "______________________________________________________________________________________\n\n"
	  << "Pressure-Species Rate Tables:\n\n";

  // pressure dependence
  for(int i = 0; i < spec_name.size(); ++i) {
    for(int t = 0; t < temperature.size(); ++t) {
      // title
      IO::out << "   Temperature = " << temperature[t] / Phys_const::kelv << " K\n\n";
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
	if(j != i) {
	  stemp = spec_name[i] + "->" + spec_name[j];
	  IO::out << std::setw(13) << stemp;
	}
      IO::out << std::setw(13) << spec_name[i] + "->";

      // escape
      for(int w = 0; w < Model::well_size(); ++w)
	if(Model::well(w).escape())
	  IO::out << std::setw(13) << Model::well(w).name(); 
      IO::out << "\n";

      // rates
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
	for(int j = 0; j < spec_name.size(); ++j)
	  if(j != i) {
	    if(rate_coef[itemp].find(std::make_pair(i, j)) != rate_coef[itemp].end())
	      IO::out << std::setw(13) << rate_coef[itemp][std::make_pair(i, j)];
	    else
	      IO::out << std::setw(13) << "***";
	  }
	// reactant loss
	if(rate_coef[itemp].find(std::make_pair(i, i)) != rate_coef[itemp].end())
	  IO::out << std::setw(13) << rate_coef[itemp][std::make_pair(i, i)];
	else
	  IO::out << std::setw(13) << "***";
	// escape rates
	for(int w = 0; w < Model::well_size(); ++w)
	  if(Model::well(w).escape()) {
	    ptemp = std::make_pair(i, Model::well_size() + Model::bimolecular_size() + w);
	    if(rate_coef[itemp].find(ptemp) != rate_coef[itemp].end())
	      IO::out << std::setw(13) << rate_coef[itemp][ptemp];
	    else
	      IO::out << std::setw(13) << "***";
	  }
	IO::out << "\n";
      }
      IO::out << std::setw(7) << "O-O";
      for(int j = 0; j < spec_name.size(); ++j)
	if(j != i) {
	  if(hp_rate_coef[t].find(std::make_pair(i, j)) != hp_rate_coef[t].end())
	    IO::out << std::setw(13) << hp_rate_coef[t][std::make_pair(i, j)];
	  else
	    IO::out << std::setw(13) << "***";
	}
      if(capture[t].find(i) != capture[t].end())
	IO::out << std::setw(13) << capture[t][i];
      else
	IO::out << std::setw(13) << "***";
      IO::out << "\n\n";
    }
  }

  IO::out << "______________________________________________________________________________________\n\n"
	  << "Temperature-Species Rate Tables:\n\n";

  // temperature dependence
  for(int i = 0; i < spec_name.size(); ++i) {
    for(int p = 0; p < pressure.size(); ++p) {
      // title
      IO::out << "   Pressure = ";
      switch(MasterEquation::pressure_unit) {
      case MasterEquation::BAR:
	IO::out << pressure[p] / Phys_const::bar << " bar\n\n";
	break;
      case MasterEquation::TORR:
	IO::out << pressure[p] / Phys_const::tor << " torr\n\n";
	break;
      case MasterEquation::ATM:
	IO::out << pressure[p] / Phys_const::atm << " atm\n\n";
	break;
      }
      IO::out << std::setw(7) << "T(K)";
      for(int j = 0; j < spec_name.size(); ++j)
	if(j != i) {
	  stemp = spec_name[i] + "->" + spec_name[j];
	  IO::out << std::setw(13) << stemp;
	}
      IO::out << std::setw(13) << spec_name[i] + "->"
	      << std::setw(13) << "Capture";

      // escape
      for(int w = 0; w < Model::well_size(); ++w)
	if(Model::well(w).escape())
	  IO::out << std::setw(13) << Model::well(w).name(); 
      IO::out << "\n";

      // rates
      for(int t = 0; t < temperature.size(); ++t) {
	IO::out << std::setw(7) << temperature[t] / Phys_const::kelv; 
	itemp = p + t * pressure.size();
	// to products
	for(int j = 0; j < spec_name.size(); ++j)
	  if(j != i) {
	    if(rate_coef[itemp].find(std::make_pair(i, j)) != rate_coef[itemp].end())
	      IO::out << std::setw(13) << rate_coef[itemp][std::make_pair(i, j)];
	    else
	      IO::out << std::setw(13) << "***";
	  }
	// reactant loss
	if(rate_coef[itemp].find(std::make_pair(i, i)) != rate_coef[itemp].end())
	  IO::out << std::setw(13) << rate_coef[itemp][std::make_pair(i, i)];
	else
	  IO::out << std::setw(13) << "***";
	// capture
	if(capture[t].find(i) != capture[t].end())
	  IO::out << std::setw(13) << capture[t][i];
	else
	  IO::out << std::setw(13) << "***";

	// escape rates
	for(int w = 0; w < Model::well_size(); ++w)
	  if(Model::well(w).escape()) {
	    ptemp = std::make_pair(i, Model::well_size() + Model::bimolecular_size() + w);
	    if(rate_coef[itemp].find(ptemp) != rate_coef[itemp].end())
	      IO::out << std::setw(13) << rate_coef[itemp][ptemp];
	    else
	      IO::out << std::setw(13) << "***";
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
	  IO::out << std::setw(13) << temperature[t] / Phys_const::kelv; 
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
	      IO::out << std::setw(13) << rate_coef[itemp][proc];
	    else
	      IO::out << std::setw(13) << "***";
	  }
	  IO::out << "\n";
	}
	IO::out << std::left << std::setw(7) << "O-O" << std::right;
	for(int t = 0; t < temperature.size(); ++t)
	  if(hp_rate_coef[t].find(proc) != hp_rate_coef[t].end())
	    IO::out << std::setw(13) << hp_rate_coef[t][proc];
	  else
	    IO::out << std::setw(13) << "***";
	IO::out << "\n\n";
      }

  // escape rates
  if(Model::escape_size()) {
    for(int spec = 0; spec < spec_name.size(); ++spec)
      for(int w = 0; w < Model::well_size(); ++w) {

      }

  }
  
  // loss
  for(int i = 0; i < spec_name.size(); ++i) {
    // title
    stemp = spec_name[i] + "-> ";
    IO::out << stemp << "\n\n";

    IO::out << std::left << std::setw(7) << "P\\T" << std::right;
    for(int t = 0; t < temperature.size(); ++t)
      IO::out << std::setw(13) << temperature[t] / Phys_const::kelv; 
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
	  IO::out << std::setw(13) << rate_coef[itemp][std::make_pair(i, i)];
	else
	  IO::out << std::setw(13) << "***";
      }
      IO::out << "\n";
    }
    IO::out << std::left << std::setw(7) << "O-O" << std::right;
    for(int t = 0; t < temperature.size(); ++t)
      if(capture[t].find(i) != capture[t].end())
	IO::out << std::setw(13) << capture[t][i];
      else
	IO::out << std::setw(13) << "***";
    IO::out << "\n\n";
  }// species cycle


  // state landscape output
  if(state_landscape.size()) {
    std::ofstream slout((base_name + ".gpi").c_str());

    for(int p = 0; p < pressure.size(); ++p) {// pressure cycle

      int tmax;
      for(tmax = temperature.size() - 1; tmax >= 0; --tmax)
	if(well_partition[p + tmax * pressure.size()].size())
	  break;

      std::map<int, MasterEquation::Group> well_map;
      for(int t = 0; t <= tmax; ++t) {
	//MasterEquation::partition_function;
      }

      for(int t = 1; t <= tmax; ++t) {
	const int t0 = p + t       * pressure.size();
	const int t1 = p + (t - 1) * pressure.size();

	std::set<MasterEquation::Group> gpool;
	for(int g1 = 0; g1 < well_partition[t1].size(); ++g1)
	  gpool.insert(well_partition[t1][g1]);

	MasterEquation::Partition new_part;
	    
	for(int g0 = 0; g0 < well_partition[t0].size(); ++g0) {
	  MasterEquation::Group new_group;
	  for(MasterEquation::Group::const_iterator it = well_partition[t0][g0].begin(); it != well_partition[t0][g0].end(); ++it)
	    for(std::set<MasterEquation::Group>::iterator g1 = gpool.begin(); g1 != gpool.end(); ++g1)
	      if(g1->find(*it) != g1->end()) {
		new_group.insert(*g1);
		gpool.erase(g1);
		break;
	      }
	  if(new_group.size())
	    new_part.push_back(new_group);
	}
	well_partition[t0] = new_part;
      }

      slout << "set term postscript portrait enhanced color\n";
      slout << "set output \"" << state_landscape << "_" << p << ".ps\"\n";
      slout << "set title \"Pressure=";

      switch(MasterEquation::pressure_unit) {
      case MasterEquation::BAR:
	slout << pressure[p] / Phys_const::bar << " bar";
	break;
      case MasterEquation::TORR:
	slout << pressure[p] / Phys_const::tor << " torr";
	break;
      case MasterEquation::ATM:
	slout << pressure[p] / Phys_const::atm << " atm";
	break;
      }
      slout << "\"\n";

      slout << "set ylabel \"T, K\"\n";

      // set ranges
      slout << "set xrange [-1:" << Model::well_size() << "]\n"; 
      slout << "set yrange [" << 100 * (int)std::floor(temperature.front() / Phys_const::kelv / 100.) << ":" 
	    << 100 * (int)std::ceil(temperature[tmax] / Phys_const::kelv / 100.) << "]\n";
    
      itemp = p + tmax * pressure.size();
      MasterEquation::Partition prev_part = well_partition[itemp];

      for(int t = tmax - 1; t >= 0; --t) {
	int tind = p + t * pressure.size();
      
	std::set<MasterEquation::Group> gpool;
	for(int i = 0; i < well_partition[tind].size(); ++i)
	  gpool.insert(well_partition[tind][i]);

	MasterEquation::Partition curr_part;

	int prev_vertex = 0;
	int curr_vertex = 0;
	for(int pi = 0; pi < prev_part.size(); prev_vertex += prev_part[pi++].size()) {
	  std::set<MasterEquation::Group>::iterator cit = gpool.begin();
	  while(cit != gpool.end())
	    if(prev_part[pi].contain(*cit)) {
	      curr_part.push_back(*cit);
	      // ... draw arrow
	      slout << "set arrow from " << prev_vertex << "," << temperature[t + 1] / Phys_const::kelv
		    << " to " << curr_vertex << "," << temperature[t] / Phys_const::kelv << " nohead lw 5\n";


	      curr_vertex += cit->size();
	      gpool.erase(cit);
	      cit = gpool.begin();
	    }
	    else
	      ++cit;
	}

	for(std::set<MasterEquation::Group>::const_iterator cit = gpool.begin(); cit != gpool.end(); ++cit)
	  curr_part.push_back(*cit);
      
	prev_part = curr_part;
      }

      // set xtics
      slout << "set xtics(";
      itemp = 0;
      for(int g = 0; g < prev_part.size(); ++g) 
	for(MasterEquation::Group::const_iterator it = prev_part[g].begin(); it != prev_part[g].end(); ++it, ++itemp) {
	  if(itemp)
	    slout << ", ";
	  slout << "\"" << Model::well(*it).name() << "\" " << itemp;
	}
      slout << ")\n";


      slout << "plot -1000 notitle\n\n";

    }// pressure cycle
  }

  return 0;
}
