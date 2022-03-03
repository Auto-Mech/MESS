/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2021, Yuri Georgievski <ygeorgi@anl.gov>

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

#include "mess.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"


//auxiliary structures for rate constants storage and output
//
struct Reactant : public std::pair<int, int> {
  //
  Reactant (int t, int i) : std::pair<int, int>(t, i) {}

  int type  () const { return first;  }
  int index () const { return second; }

  bool dummy () const;
  
  const std::string& name () const;

  double weight (double) const;
};

std::ostream& operator<< (std::ostream& to, const Reactant& r) { return to << r.name(); }

std::ostream& operator<< (std::ostream& to, const std::pair<Reactant, Reactant>& r)
{
  std::string s;

  s = r.first.name() + "->" + r.second.name();

  return to << s;
}

const std::string& Reactant::name () const
{
  const char funame [] = "Reactant::name: ";
  
  switch(first) {
    //
  case Model::WELL:
    //
    return Model::well(second).name();

  case Model::BIMOLECULAR:
    //
    return Model::bimolecular(second).name();

  default:
    //
    std::cerr << funame << "wrong reactant type: " << first << "\n";

    throw Error::Logic();
  }
}
    
bool Reactant::dummy () const
{
  const char funame [] = "Reactant::name: ";
  
  switch(first) {
    //
  case Model::WELL:
    //
    return false;

  case Model::BIMOLECULAR:
    //
    return Model::bimolecular(second).dummy();

  default:
    //
    std::cerr << funame << "wrong reactant type: " << first << "\n";

    throw Error::Logic();
  }
}
    
double Reactant::weight (double t) const
{
  const char funame [] = "Reactant::weight: ";
  
  switch(first) {
    //
  case Model::WELL:
    //
    return Model::well(second).weight(t) / std::exp(Model::well(second).ground() / t) * Phys_const::herz;

  case Model::BIMOLECULAR:
    //
    if(Model::bimolecular(second).dummy())
      //
      return -1.;
    
    return Model::bimolecular(second).weight(t) / std::exp(Model::bimolecular(second).ground() / t)
      //
      * Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  default:
    //
    std::cerr << funame << "wrong reactant type: " << first << "\n";

    throw Error::Logic();
  }
}
    
struct RateData: public std::map<std::pair<int, int>, std::map<std::set<std::pair<int, int> >, double> > {
  //
  double rate (int t, int p, Reactant r1, Reactant r2) const;
};

double RateData::rate (int t, int p, Reactant r1, Reactant r2) const
{
  const char funame [] = "RateData:rate: ";
  
  const_iterator i = find(std::make_pair(t, p));

  if(i == end()) {
    //
    std::cerr << funame << "temperature and pressure indices out of range\n";

    throw Error::Logic();
  }

  std::set<std::pair<int, int> > rr;

  rr.insert(r1);

  rr.insert(r2);
  
   std::map<std::set<std::pair<int, int> >, double>::const_iterator j = i->second.find(rr);

  if(j != i->second.end())
    //
    return j->second;

  return 0.;
}

int main (int argc, char* argv [])
{
  const char funame [] = "cool: ";

  if (argc < 2) {
    //
    std::cout << "usage: mess_test input_file\n";
    
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

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  std::pair<int, int> ptemp;

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
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    
    return 1;
  }

  // output streams
  //
  IO::out.open((base_name + ".out").c_str());
  
  if(!IO::out.is_open()) {
    //
    std::cerr << funame << "cannot open rate output file: " << base_name + ".out" << "\n";

    throw Error::Init();
  }
  
  IO::log.open((base_name + ".log").c_str());

  if(!IO::log.is_open()) {
    //
    std::cerr << funame << "cannot open log file: " << base_name + ".log" << "\n";

    throw Error::Init();
  }
  
  IO::aux.open((base_name + ".aux").c_str());
  
  if(!IO::aux.is_open()) {
    //
    std::cerr << funame << "cannot open rate auxiliary output file: " << base_name + ".aux" << "\n";

    throw Error::Init();
  }
  
  IO::Marker funame_marker(funame);
  
  /********************************************************************************************
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key    model_key("Model"                      );
  Key bar_pres_key("PressureList[bar]"          );
  Key tor_pres_key("PressureList[torr]"         );
  Key atm_pres_key("PressureList[atm]"          );
  Key     temp_key("TemperatureList[K]"         );
  Key     tema_key("TemperatureArray[K]"        );
  Key     rtem_key("ReferenceTemperature[K]"    );
  Key    rpres_key("ReferencePressure[atm]"     );
  Key     etot_key("EnergyStepOverTemperature"  );
  Key     xtot_key("ExcessEnergyOverTemperature");
  Key ener_cut_key("EnergyCutoffOverTemperature");
  Key ener_win_key("EnergyWindowOverTemperature");
  Key     emax_key("ModelEnergyLimit[kcal/mol]" );
  Key  adm_bor_key("AtomDistanceMin[bohr]"      );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );
  Key     xtun_key("TunnelingActionCutoff"      );
  Key well_cut_key("WellCutoff"                 );
  Key well_ext_key("WellExtension"              );
  Key ext_corr_key("ExtensionCorrection"        );
  Key eval_max_key("ChemicalEigenvalueMax"      );
  Key eval_min_key("ChemicalEigenvalueMin"      );
  Key well_red_key("WellReductionThreshold"     );
  Key time_ste_key("TimePropagationStep"        );
  Key time_lim_key("TimePropagationLimit"       );
  Key    float_key("FloatType"                  );
  Key     calc_key("CalculationMethod"          );
  Key  rat_red_key("ReductionMethod"            );
  Key      wpm_key("WellPartitionMethod"        );
  Key      wpt_key("WellProjectionThreshold"    );
  Key     wtol_key("WellProjectionTolerance"    );
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

  std::vector<std::string> ped_spec;// product energy distribution pairs, verbal specification
  
  std::vector<std::string> reduction_scheme;
  
  std::vector<double>      pressure;
  
  std::vector<double>      temperature;
  
  std::string              micro_rate_file;
  
  double                   micro_ener_max  = 0.;
  
  double                   micro_ener_min  = 0.;
  
  double                   micro_ener_step = -1.;
  
  std::string              state_landscape;

  std::string token, comment, line;
  
  while(from >> token) {
    //
    // main input group
    //
    if(model_key == token) {
      //
      // global energy limit check
      //
      if(!Model::is_energy_limit()) {
	//
	std::cerr << funame << "model energy limit should be set BEFORE model input section\n";
	
	throw Error::Init();
      }

      // main initialization
      //
      try {
	//
	Model::init(from);
      }  
      catch(Error::General) {
	//
	IO::log << std::flush;
	
	throw;
      }
      break;
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

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

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

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
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
      
      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
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
      IO::LineInput lin(from);
      
      if(!(lin >> Model::reactant)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
    }
    // eigenvector number to output
    //
    else if(evec_num_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> itemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(itemp < 0) {
	//
        std::cerr << funame << token << ": out of range: " << itemp << "\n";
	
        throw Error::Range();
      }

      MasterEquation::evec_out_num = itemp;
    }
    // temperature
    //
    else if(temp_key == token) {
      //
      
      if(temperature.size()) {
	//
	std::cerr << funame << token <<  ": already initialized\n";
	
	throw Error::Init();
      }
      
      IO::LineInput lin(from);
      
      std::set<double> data;
      
      while(lin >> dtemp) {
	//
	if(dtemp <= 0.) {
	  //
	  std::cerr << funame << token << ": should be positive\n";
	  
	  throw Error::Range();
	}
	data.insert(dtemp * Phys_const::kelv);
      }
      
      if(!data.size()) {
	//
        std::cerr << funame << token << ": no data\n";

        throw Error::Init();
      }

      temperature.resize(data.size());

      itemp = 0;
      
      for(std::set<double>::const_iterator it = data.begin(); it != data.end(); ++it, ++itemp)
	//
	temperature[itemp] = *it;
    }
    // temperature
    //
    else if(tema_key == token) {
      //
      if(temperature.size()) {
	//
	std::cerr << funame << token <<  ": already initialized\n";
	
	throw Error::Init();
      }
      
      std::string format = ": format: <initial_value>(int) <final_value>(int) <step>(int)";
      
      IO::LineInput lin(from);

      int start, stop, step;

      if(!(lin >> start >> stop >> step)) {
	//
	std::cerr << funame << token << ": corrupted" << format << "\n";

	throw Error::Input();
      }

      if(start <= 0 || stop < start || step <= 0) {
	//
	std::cerr << funame << token << ": out of range: " << start << " " << stop << " " << step << format << "\n";

	throw Error::Range();
      }

      for(int t = start; t <= stop; t += step) 
	//
	temperature.push_back(t * Phys_const::kelv);
      
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

      IO::LineInput lin(from);
      
      std::set<double> data;

      while(lin >> dtemp) {
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
	else if(atm_pres_key == token) {
	  //
	  dtemp *= Phys_const::atm;
	  
	  MasterEquation::pressure_unit = MasterEquation::ATM;
	}
	else if(tor_pres_key == token) {
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

      pressure.resize(data.size());

      itemp = 0;
      
      for(std::set<double>::const_iterator p = data.begin(); p != data.end(); ++p, ++itemp)
	//
	pressure[itemp] = *p;
    }
    // reference temperature
    //
    else if(rtem_key == token) {
      //
      if(MasterEquation::reference_temperature > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Input();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp <<"\n";
	
	throw Error::Range();
      }

      MasterEquation::reference_temperature = dtemp * Phys_const::kelv;
    }
    // reference pressure
    //
    else if(rpres_key == token) {
      //
      if(MasterEquation::reference_pressure > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Input();
      }
      
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp <<"\n";
	
	throw Error::Range();
      }

      MasterEquation::reference_pressure = dtemp * Phys_const::atm;
    }
    // energy step over temperature
    //
    else if(etot_key == token) {
      //
      if(MasterEquation::energy_step_over_temperature > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Input();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
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
      if(MasterEquation::excess_energy_over_temperature > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Input();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": should be positive\n";
	
	throw Error::Range();
      }

      MasterEquation::excess_energy_over_temperature = dtemp;
    }
    // excess energy over temperature
    //
    else if(ener_win_key == token) {
      //
      if(MasterEquation::energy_window_over_temperature > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Input();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": should be positive\n";
	
	throw Error::Range();
      }

      MasterEquation::energy_window_over_temperature = dtemp;
    }
    // energy cutoff over temperature
    //
    else if(ener_cut_key == token) {
      //
      if(MasterEquation::energy_cutoff_over_temperature > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
        std::cerr << funame << token << ": should be positive\n";
	
        throw Error::Range();
      }
      
      MasterEquation::energy_cutoff_over_temperature = dtemp;
    }
    // model energy limit
    //
    else if(emax_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      Model::set_energy_limit(dtemp * Phys_const::kcal);
    }
    // minimal interatomic distance
    //
    else if(adm_bor_key == token || adm_ang_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

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
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
        std::cerr << funame << token << ": should be positive\n";
	
        throw Error::Range();
      }

      Model::Tunnel::set_action_max(dtemp);
    }
    // calculation method (inactive)
    //
    else if(calc_key == token) {
      //
      std::getline(from, comment);
    }
    // microcanonical rate limit
    //
    else if(rate_max_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> MasterEquation::rate_max)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(MasterEquation::rate_max <= 0.) {
	//
        std::cerr << funame << token << ": should be positive\n";
	
        throw Error::Range();
      }

      MasterEquation::rate_max *= Phys_const::herz;
    }
    // well cutoff (inactive)
    //
    else if(well_cut_key == token) {
      //
      std::getline(from, comment);
    }
    // well extention
    //
    else if(well_ext_key == token) {
      //
      IO::LineInput lin(from);

      if(lin >> dtemp) {
	//
	if(dtemp < 0. || dtemp > 1.) {
	  //
	  std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	  throw Error::Range();
	}

	MasterEquation::well_extension = dtemp;
      }
      else {
	//
	std::cerr << funame << token << ": well extension parameter has not been provided: suggested 0.2\n";

	throw Error::Init();
      }
    }
    // well extention correction factor (inactive)
    //
    else if(ext_corr_key == token) {
      //
      std::getline(from, comment);
    }
    // chemical relaxation to collision frequency ratio
    //
    else if(well_red_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	throw Error::Range();
      }

      MasterEquation::reduction_threshold = dtemp;
    }
    // time-propagation step over fastest relaxation time-scale
    //
    else if(time_ste_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	throw Error::Range();
      }

      MasterEquation::time_propagation_step = dtemp;
    }
    // time-propagation limit over collisional frequency time-scale
    //
    else if(time_lim_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	throw Error::Range();
      }

      MasterEquation::time_propagation_limit = dtemp;
    }
    // float type (inactive)
    //
    else if(float_key == token) {
      //
      std::getline(from, comment);
    }
    // maximal chemical eigenvalue
    //
    else if(eval_max_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }
      
      if(dtemp < -1.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	throw Error::Range();
      }

      MasterEquation::chemical_threshold = dtemp;
    }
    // minumal chemical eigenvalue (inactive)
    //
    else if(eval_min_key == token) {
      //
      std::getline(from, comment);
    }
    // number of closest reductions to print (inactive)
    //
    else if(red_out_key == token) {
      //
      std::getline(from, comment);
    }
    // species reduction algorithm for low-eigenvalue method (inactive)
    //
    else if(rat_red_key == token) {
      //
      std::getline(from, comment);
    }
    // well partition method (inactive)
    //
    else if(wpm_key == token) {
      //
      std::getline(from, comment);
    }
    // well projection threshold
    //
    else if(wpt_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp <= 0. || dtemp >= 1.) {
	//
        std::cerr << funame << token << ": out of range: " << dtemp << ": suggested value: 0.2\n";
	
        throw Error::Range();
      }
	
      MasterEquation::well_projection_threshold = dtemp;
    }
    // well projection tolerance
    //
    else if(wtol_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp < 0. || dtemp >= 1.) {
	//
        std::cerr << funame << token << ": out of range: " << dtemp << ": suggested value: 0.0001\n";
	
        throw Error::Range();
      }
	
      MasterEquation::well_projection_tolerance = dtemp;
    }
    // default reduction scheme (inactive)
    //
    else if(def_red_key == token) {
      //
      std::getline(from, comment);
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

      lin >> count;

      while(count) {
	//
	lin.read_line(from);
	  
	if(!(lin >> stemp)) // empty line
	  //
	  continue;

	if(IO::end_key() == stemp)
	  //
	  break;
	
	while(lin >> dtemp) {
	  //
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
	  
	  MasterEquation::hot_energy[stemp].insert(dtemp);
	}

	if(!MasterEquation::hot_energy[stemp].size()) {
	  //
	  std::cerr << funame << token << ": " <<  stemp << ": no energies provided\n";
	  
	  throw Error::Init();
	}

	--count;
      }
    }
    // microscopic rates file
    //
    else if(mic_rout_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> micro_rate_file)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
    }
    // microscopic rate energy maximum
    //
    else if(mic_emax_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> micro_ener_max)) {
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
      IO::LineInput lin(from);
      
      if(!(lin >> micro_ener_min)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      micro_ener_min *= Phys_const::kcal;
    }
    // microscopic rate energy step
    //
    else if(mic_step_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> micro_ener_step)) {
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

  /*************************** CHECKING ******************************************/

  if(!from) {
    //
    std::cerr << funame << "corrupted\n";
    
    throw Error::Input();
  }

  if(!Model::isinit()) {
    //
    std::cerr << funame << "model not initialized\n";
    
    throw Error::Init();
  }

  if(!temperature.size()) {
    //
    std::cerr << funame << "temperature list has not been initialized\n";
    
    throw Error::Init();
  }

  if(!pressure.size()) {
    //
    std::cerr << funame << "pressure list has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::reference_temperature < 0.) {
    //
    std::cerr << funame << rtem_key << ": reference temperature has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::reference_pressure < 0.) {
    //
    std::cerr << funame << rpres_key << ": reference pressure has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::energy_step_over_temperature < 0.) {
    //
    std::cerr << funame << etot_key << ": energy step over temperature has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::energy_window_over_temperature < 0.) {
    //
    std::cerr << funame << ener_win_key << ": energy window over temperature has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::excess_energy_over_temperature < 0.) {
    //
    std::cerr << funame << xtot_key << ": excess energy over temperature (upper part of the energy window)"
      //
	      << " has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::energy_cutoff_over_temperature < 0.) {
    //
    std::cerr << funame << ener_cut_key << ": energy cutoff over temperature (lower part of the energy window)"
      //
	      << " has not been initialized\n";
    
    throw Error::Init();
  }

  if(MasterEquation::reduction_threshold < 0.) {
    //
    std::cerr << funame << well_red_key << ": well reduction threshold has not been initialized\n";

    throw Error::Init();
  }
  
  if(MasterEquation::well_projection_threshold < 0.) {
    //
    std::cerr << funame << wpt_key << ": well projection threshold has not been initialized\n";

    throw Error::Init();
  }
  
  if(MasterEquation::well_projection_tolerance < 0.) {
    //
    std::cerr << funame << wtol_key << ": well projection tolerance has not been initialized\n";

    throw Error::Init();
  }
  
  if(MasterEquation::chemical_threshold < -1.) {
    //
    std::cerr << funame << eval_max_key << ": chemical threshold has not been initialized\n";

    throw Error::Init();
  }
  
  if(reduction_scheme.size())
    //
    MasterEquation::set_default_partition(reduction_scheme);

  if(ped_spec.size())
    //
    MasterEquation::set_ped_pair(ped_spec);

  /************************** MICROSCOPIC RATE COEFFICIENTS **********************************/

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

    // well cycle
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      micro_out << std::setw(15) << "E, kcal/mol"
	//
		<< std::setw(15) << "D, mol/kcal";

      for(int b = 0; b < Model::inner_barrier_size(); ++b) {
	//
	const int w1 = Model::inner_connect(b).first;
	
	const int w2 = Model::inner_connect(b).second;
	
	if(w1 == w) {
	  //
	  stemp = Model::well(w).name() + "->" + Model::well(w2).name();
	  
	  micro_out << std::setw(15) << stemp;
	}
	else if(w2 == w) {
	  //
	  stemp = Model::well(w).name() + "->" + Model::well(w1).name();
	  
	  micro_out << std::setw(15) << stemp;
	}
      }

      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(Model::outer_connect(b).first == w) {
	  //
	  const int p = Model::outer_connect(b).second;
	  
	  stemp = Model::well(w).name() + "->" + Model::bimolecular(p).name();
	  
	  micro_out << std::setw(15) << stemp;
	}
      
      micro_out << "\n";

      // energy cycle
      //
      for(double ener = micro_ener_min; ener <= micro_ener_max; ener += micro_ener_step) {
	//
	double states = Model::well(w).states(ener);
	
	micro_out << std::setw(15) << ener   / Phys_const::kcal
	  //
		  << std::setw(15) << states * Phys_const::kcal;

	for(int b = 0; b < Model::inner_barrier_size(); ++b)
	  //
	  if(Model::inner_connect(b).first == w || Model::inner_connect(b).second == w)
	    //
	    if(states != 0.) {
	      //
	      micro_out << std::setw(15) << Model::inner_barrier(b).states(ener) / states / 2. / M_PI / Phys_const::herz;
	    }
	    else
	      //
	      micro_out << std::setw(15) << "***";

	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(Model::outer_connect(b).first == w)
	    //
	    if(states != 0.) {
	      //
	      micro_out << std::setw(15) << Model::outer_barrier(b).states(ener) / states / 2. / M_PI / Phys_const::herz;
	    }
	    else
	      //
	      micro_out << std::setw(15) << "***";
	
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
    
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(!Model::bimolecular(p).dummy())
	  //
	  for(int f = 0; f < 2; ++f)
	    //
	    if(Model::bimolecular(p).mode(f) == Model::DENSITY)
	      //
	      micro_out << std::setw(15) << Model::bimolecular(p).fragment_name(f);
      
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
		micro_out << std::setw(15) << Model::bimolecular(p).states(f, ener) * Phys_const::kcal;

	micro_out << "\n";
      }//
      //
    }//
    //
  }// micro output

  if(Model::no_run())
    //
    return 0;
  
  /*********** PRESSURE AND TEMPERATURE DEPENDENT RATE COEFFICIENTS CALCULATION*************/

  // units
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;
  
  // abbreviations
  //
  typedef std::set<std::pair<int, int> >                              reac_t;

  typedef std::map<reac_t, double>                                    rate_t;

  IO::out << std::setprecision(3) << "\n";

  RateData storage;
  
  for(std::vector<double>::const_iterator ti = temperature.begin(); ti != temperature.end(); ++ti) {
    //
    // thermal landscape
    //
    Model::landscape_t landscape = Model::kinetic_landscape(Model::thermal_energy_map(*ti));
    
    for(std::vector<double>::const_iterator pi = pressure.begin(); pi != pressure.end(); ++pi) {
      //
      IO::log << IO::log_offset << "Temperature = " << *ti / Phys_const::kelv << " K"  << "   Pressure = ";
  
      switch(MasterEquation::pressure_unit) {
	//
      case MasterEquation::BAR:
	//
	IO::log << *pi / Phys_const::bar << " bar";
    
	break;
	//
      case MasterEquation::TORR:
	//
	IO::log << *pi / Phys_const::tor << " torr";
    
	break;
	//
      case MasterEquation::ATM:
	//
	IO::log << *pi / Phys_const::atm << " atm";
    
	break;
      }
      
      IO::log << "\n\n";

      // barrier-specific rate constant calculation
      //
      rate_t reaction_map;
      
      for(Model::landscape_t::const_iterator li = landscape.begin(); li != landscape.end(); ++li) {
	//
	const int& type = li->first.first;
	
	const int& b    = li->first.second;

	IO::log << IO::log_offset;

	switch(type) {
	  //
	case Model::INNER:
	  //
	  if(li->second.size() != 2) {
	    //
	    std::cerr << "there should be exactly two chemical graphs for the "
	      //
		      << Model::inner_barrier(b).name()
	      //
		      << " inner barrier\n";

	    throw Error::Logic();
	  }
	    
	  IO::log << Model::inner_barrier(b).name();
	    
	  IO::log  << ": rate-limiting barrier for chemical groups: ";

	  for(std::set<int>::const_iterator w = li->second.front().well_set.begin(); w != li->second.front().well_set.end(); ++w) {
	    //
	    if(w != li->second.front().well_set.begin())
	      //
	      IO::log << ", ";
	      
	    IO::log << "  " << Model::well(*w).name();
	  }
	    
	  IO::log << " and ";
	  
	  for(std::set<int>::const_iterator w = li->second.back().well_set.begin(); w != li->second.back().well_set.end(); ++w) {
	    //
	    if(w != li->second.back().well_set.begin())
	      //
	      IO::log << ", ";
	      
	    IO::log << "  " << Model::well(*w).name();
	  }
	    
	  break;
	  //
	case Model::OUTER:
	  //
	  if(li->second.size() != 1) {
	    //
	    std::cerr << "there should be exactly one chemical graph for the "
	      //
		      << Model::outer_barrier(b).name()
	      //
		      << " outer barrier\n";

	    throw Error::Logic();
	  }
	    
	  IO::log << Model::outer_barrier(b).name();
	  
	  IO::log  << ": rate-limiting barrier for the chemical group: ";

	  for(std::set<int>::const_iterator w = li->second.front().well_set.begin(); w != li->second.front().well_set.end(); ++w) {
	    //
	    if(w != li->second.front().well_set.begin())
	      //
	      IO::log << ", ";
	      
	    IO::log << "  " << Model::well(*w).name();
	  }
	    
	  IO::log << " and bimolecular " << Model::bimolecular(Model::outer_connect(b).second).name();
	  
	  break;
	  //
	default:
	  //
	  std::cerr << funame << "wrong barrier type: " << type << "\n";

	  throw Error::Logic();
	}

	IO::log << "\n\n";

	// reactive complex
	//
	std::list<std::map<std::pair<int, int>, double> > bf;

	double er, ec;
	
	switch(type) {
	  //
	case Model::INNER:
	  //
	  dtemp = Model::inner_barrier(b).thermal_energy(*ti);

	  er = dtemp + MasterEquation::excess_energy_over_temperature * *ti;

	  ec = dtemp - MasterEquation::energy_cutoff_over_temperature * *ti;

	  for(std::list<Model::ChemGraph>::const_iterator g = li->second.begin(); g != li->second.end(); ++g) {
	    //
	    MasterEquation::ReactiveComplex rc(*ti, *pi, er, ec, *g);

	    bf.push_back(rc.branching_fraction(li->first));
	  }

	  dtemp = *ti / 2. / M_PI * Model::inner_barrier(b).weight(*ti) / std::exp(Model::inner_barrier(b).ground() / *ti);

	  for(std::map<std::pair<int, int>, double>::const_iterator r = bf.front().begin(); r != bf.front().end(); ++r) {
	    //
	    for(std::map<std::pair<int, int>, double>::const_iterator p = bf.back().begin(); p != bf.back().end(); ++p) {
	      //
	      std::set<std::pair<int, int> > rp;

	      rp.insert(r->first);

	      rp.insert(p->first);
		
	      reaction_map[rp] += dtemp * r->second * p->second;
	    }
	  }
	    
	  break;
	  //
	case Model::OUTER:
	  //
	  dtemp = Model::outer_barrier(b).thermal_energy(*ti);

	  er = dtemp + MasterEquation::excess_energy_over_temperature * *ti;

	  ec = dtemp - MasterEquation::energy_cutoff_over_temperature * *ti;

	  for(std::list<Model::ChemGraph>::const_iterator g = li->second.begin(); g != li->second.end(); ++g) {
	    //
	    MasterEquation::ReactiveComplex rc(*ti, *pi, er, ec, *g);

	    bf.push_back(rc.branching_fraction(li->first));
	  }

	  dtemp = *ti / 2. / M_PI * Model::outer_barrier(b).weight(*ti) / std::exp(Model::outer_barrier(b).ground() / *ti);

	  for(std::map<std::pair<int, int>, double>::const_iterator r = bf.front().begin(); r != bf.front().end(); ++r) {
	    //
	    std::set<std::pair<int, int> > rp;

	    rp.insert(std::make_pair((int)Model::BIMOLECULAR, Model::outer_connect(b).second));

	    rp.insert(r->first);
		
	    reaction_map[rp] += dtemp * r->second;
	  }
	    
	  break;
	  //
	default:
	  //
	  std::cerr << funame << "wrong barrier type: " << type << "\n";

	  throw Error::Logic();
	}
      }
      
      storage[std::make_pair(int(ti - temperature.begin()), int(pi - pressure.begin()))] = reaction_map;
      //
    }// pressure cycle
    //
  }// temperature cycle

  std::vector<Reactant> reactant;

  for(int w = 0; w < Model::well_size(); ++w)
    //
    if(Model::well_exclude_group.find(Model::well(w).name()) == Model::well_exclude_group.end())
      //
      reactant.push_back(Reactant(Model::WELL, w));

  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    reactant.push_back(Reactant(Model::BIMOLECULAR, p));
  
  /*********************************************************************************************
   *************************************** RATE OUTPUT *****************************************
   *********************************************************************************************/

  // maximal print length
  //
  const int count_max = 80;

  // field width
  //
  const int width = 10;

  // initial width
  //
  const int first = 7;

  int start, count;
  
  IO::out << "Unimolecular Rate Units: 1/sec;  Bimolecular Rate Units: cm^3/sec\n\n"
    //
	  << "______________________________________________________________________________________\n\n";

  IO::out<< "Species-Species Rate Tables:\n\n";

  for(int p = 0; p < pressure.size(); ++p) {
    //
    for(int t = 0; t < temperature.size(); ++t) {
      //
      IO::out << "Temperature = " << temperature[t] / Phys_const::kelv <<  " K    Pressure = ";
      
      switch(MasterEquation::pressure_unit) {
	//
      case MasterEquation::BAR:
	//
	IO::out << pressure[p] / Phys_const::bar << " bar";
	  
	break;
	//
      case MasterEquation::TORR:
	//
	IO::out << pressure[p] / Phys_const::tor << " torr";
	  
	break;
	//
      case MasterEquation::ATM:
	//
	IO::out << pressure[p] / Phys_const::atm << " atm";
	  
	break;
	//
      default:
	//
	std::cerr << funame << "wrong pressure unit index: " << MasterEquation::pressure_unit << "\n";

	throw Error::Logic();
      }
      
      IO::out << "\n\n";

      start = 0;

      while(start < reactant.size()) {
	//
	int i;
	
	IO::out << std::left << std::setw(first) << "From\\to" << std::right;

	count = first;

	for(i = start; count < count_max && i < reactant.size(); ++i, count += width)
	  //
	  IO::out << std::setw(width) << reactant[i];

	IO::out << "\n";
	    
	for(int r = 0; r < reactant.size(); ++r) {
	  //
	  const double weight = reactant[r].weight(temperature[t]);

	  if(reactant[r].dummy())
	    //
	    continue;
	
	  IO::out << std::left << std::setw(first) << reactant[r] << std::right;
	
	  count = first;
	    
	  for(i = start; count < count_max && i < reactant.size() ; ++i, count += width)
	    //
	    IO::out << std::setw(width) << storage.rate(t, p, reactant[r], reactant[i]) / weight;

	  IO::out << "\n";
	}

	start = i;
	
	IO::out << "\n";
      }
    }
  }

  IO::out<< "Temperature-Species Rate Tables:\n\n";

  for(int p = 0; p < pressure.size(); ++p) {

    for(int r = 0; r < reactant.size(); ++r) {
      //
      if(reactant[r].dummy())
	//
	continue;
	
      IO::out << "Pressure = ";
      
      switch(MasterEquation::pressure_unit) {
	//
      case MasterEquation::BAR:
	//
	IO::out << pressure[p] / Phys_const::bar << " bar";
	  
	break;
      case MasterEquation::TORR:
	//
	IO::out << pressure[p] / Phys_const::tor << " torr";
	  
	break;
      case MasterEquation::ATM:
	//
	IO::out << pressure[p] / Phys_const::atm << " atm";
	  
	break;
      default:
	//
	std::cerr << funame << "wrong pressure unit index: " << MasterEquation::pressure_unit << "\n";
	
	throw Error::Logic();
      }
      
      IO::out << "   Reactant = " << reactant[r] << "\n\n";

      start = 0;

      while(start < reactant.size()) {
	//
	int i;
	
	IO::out << std::left << std::setw(first) << "T, K" << std::right;

	count = first;
	
	for(i = start; count < count_max && i < reactant.size(); ++i, count += width)
	  //
	  IO::out << std::setw(width) << reactant[i];

	IO::out << "\n";
	    
	for(int t = 0; t < temperature.size(); ++t) {
	  //
	  //
	  const double weight = reactant[r].weight(temperature[t]);

	  IO::out << std::left << std::setw(first) << temperature[t] / Phys_const::kelv << std::right;
	
	  count = first;
	    
	  for(i = start; count < count_max && i < reactant.size() ; ++i, count += width)
	    //
	    IO::out << std::setw(width) << storage.rate(t, p, reactant[r], reactant[i]) / weight;

	  IO::out << "\n";
	}

	start = i;
	
	IO::out << "\n";
      }

      IO::out << "\n";
    }
  }
  
  IO::out<< "Pressure-Species Rate Tables:\n\n";

  for(int t = 0; t < temperature.size(); ++t) {
    //
    for(int r = 0; r < reactant.size(); ++r) {
      //
      if(reactant[r].dummy())
	//
	continue;
	
      const double weight = reactant[r].weight(temperature[t]);

      IO::out << "Temperature = " << temperature[t] / Phys_const::kelv << " K   Reactant = " << reactant[r] << "\n\n";

      start = 0;

      while(start < reactant.size()) {
	//
	int i;
	
	stemp = "P, ";

	switch(MasterEquation::pressure_unit) {
	  //
	case MasterEquation::BAR:
	  //
	  stemp += "bar";
	  
	  break;
	  //
	case MasterEquation::TORR:
	  //
	  stemp += "torr";
	  
	  break;
	case MasterEquation::ATM:
	  //
	  stemp += " atm";
	  
	  break;
	  //
	default:
	  //
	  std::cerr << funame << "wrong pressure unit index: " << MasterEquation::pressure_unit << "\n";

	  throw Error::Logic();
	}

	IO::out << std::left << std::setw(first) << stemp << std::right;
	
	count = first;
	
	for(i = start; count < count_max && i < reactant.size(); ++i, count += width)
	  //
	  IO::out << std::setw(width) << reactant[i];

	IO::out << "\n";
	    
	for(int p = 0; p < pressure.size(); ++p) {
	  //
	  IO::out << std::left << std::setw(first);

	  switch(MasterEquation::pressure_unit) {
	    //
	  case MasterEquation::BAR:
	    //
	    IO::out << pressure[p] / Phys_const::bar;
	  
	    break;
	    //
	  case MasterEquation::TORR:
	    //
	    IO::out << pressure[p] / Phys_const::tor;
	  
	    break;
	    //
	  case MasterEquation::ATM:
	    //
	    IO::out << pressure[p] / Phys_const::atm;
	    
	    break;
	    //
	  default:
	    //
	    std::cerr << funame << "wrong pressure unit index: " << MasterEquation::pressure_unit << "\n";

	    throw Error::Logic();
	  }

	  IO::out << std::right;
	
	  count = first;
	    
	  for(i = start; count < count_max && i < reactant.size() ; ++i, count += width)
	    //
	    IO::out << std::setw(width) << storage.rate(t, p, reactant[r], reactant[i]) / weight;

	  IO::out << "\n";
	}

	start = i;
	
	IO::out << "\n";
      }

      IO::out << "\n";
    }
  }
  
  return 0;
}
