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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/resource.h>

#include "libmess/mess.hh"
#include "libmess/key.hh"
#include "libmess/units.hh"
#include "libmess/io.hh"


//auxiliary structures for rate constants storage and output
//
struct Reactant : public std::pair<int, int> {
  //
  Reactant (int t, int i) : std::pair<int, int>(t, i) {}

  int type  () const { return first;  }
  int index () const { return second; }

  bool dummy () const;
  
  std::string name () const;

  double weight (double) const;
};

std::ostream& operator<< (std::ostream& to, const Reactant& r) { return to << r.name(); }

std::ostream& operator<< (std::ostream& to, const std::pair<Reactant, Reactant>& r)
{
  std::string s;

  s = r.first.name() + "->" + r.second.name();

  return to << s;
}

std::string Reactant::name () const
{
  const char funame [] = "Reactant::name: ";
  
  switch(first) {
    //
  case Model::WELL:
    //
    return Model::well(second).short_name();

  case Model::BIMOLECULAR:
    //
    return Model::bimolecular(second).short_name();

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

class DoublePositive {
  //
  double _val;
  
public:
  //
  DoublePositive (double d) : _val(d) {};

  DoublePositive operator/ (double d) const { return _val / d; }

  friend std::ostream& operator<< (std::ostream&, const DoublePositive&);
};

std::ostream& operator<< (std::ostream& to, const DoublePositive& d)
{
  if(d._val > 0.)
    //
    return to << d._val;
 
  return to << "***";
}

struct RateData: public std::map<std::pair<int, int>, std::map<std::set<std::pair<int, int> >, double> > {
  //
  DoublePositive rate (int t, int p, Reactant r1, Reactant r2) const;
};

DoublePositive RateData::rate (int t, int p, Reactant r1, Reactant r2) const
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

  return -1.;
}

int main (int argc, char* argv [])
{
  const char funame [] = "mess_test: ";

  using Model::operator<<;
  
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
  IO::log.open((base_name + ".log").c_str());

  if(!IO::log.is_open()) {
    //
    std::cerr << funame << "cannot open log file: " << base_name + ".log" << "\n";

    throw Error::Init();
  }

  IO::log << std::setprecision(Model::log_precision);
  
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
  Key prec_out_key("OutPrecision"               );
  Key prec_log_key("LogPrecision"               );
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
  Key     rout_key("RateOutput"                 );
  Key  adm_ang_key("AtomDistanceMin[angstrom]"  );
  Key     xtun_key("TunnelingActionCutoff"      );
  Key well_cut_key("WellCutoff"                 );
  Key well_ext_key("WellExtension"              );
  Key chem_max_key("ChemicalThreshold"          );
  Key chem_min_key("ChemicalTolerance"          );
  Key well_red_key("WellReductionThreshold"     );
  Key prop_stp_key("TimePropagationStep"        );
  Key prop_int_key("TimePropagationInterval"    );
  Key prop_lim_key("TimePropagationLimit"       );
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
      Model::init(from);
      
      break;
    }
    // out stream precision
    //
    else if(prec_out_key == token) {
      //
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
    }
    // log stream precision
    //
    else if(prec_log_key == token) {
      //
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

      IO::log << std::setprecision(itemp);
    }
    // custom rate output
    //
    else if(rout_key == token) {
      //
      if(IO::out.is_open()) {
	//
        std::cerr << funame << token << ": allready opened\n";
	
        throw Error::Input();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
        std::cerr << funame << token << ": bad input\n";
	
        throw Error::Input();
      }
      
      IO::out.open(stemp.c_str());
      
      if(!IO::out) {
	//
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	
        throw Error::Input();
      }
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
    // well cutoff
    //
    else if(well_cut_key == token) {
      //
      IO::LineInput lin(from);

      if(!(lin >> dtemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << ": suggested value: 10\n";
	
	throw Error::Range();
      }

      MasterEquation::well_cutoff = dtemp;
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
	std::cerr << funame << token << ": has not been provided: suggested value: 0.5\n";

	throw Error::Init();
      }
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
    // time-propagation step over collisional frequency time-scale
    //
    else if(prop_stp_key == token) {
      //
      IO::LineInput lin(from);
      
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

      MasterEquation::time_propagation_step = dtemp;
    }
    // time-propagation interval over collisional frequency time-scale
    //
    else if(prop_int_key == token) {
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

      MasterEquation::time_propagation_interval = dtemp;
    }
    // time-propagation limit over collisional frequency time-scale
    //
    else if(prop_lim_key == token) {
      //
      IO::LineInput lin(from);
      
      if(!(lin >> dtemp)) {
	//
        std::cerr << funame << token << ": corrupted\n";
	
        throw Error::Input();
      }

      if(dtemp < 1.) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << "\n";
	
	throw Error::Range();
      }

      MasterEquation::time_propagation_limit = dtemp;
    }
    // maximal chemical eigenvalue
    //
    else if(chem_max_key == token) {
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
    // minumal chemical eigenvalue 
    //
    else if(chem_min_key == token) {
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

      MasterEquation::chemical_tolerance = dtemp;
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

  if(MasterEquation::energy_step_over_temperature < 0.) {
    //
    std::cerr << funame << etot_key << ": energy step over temperature has not been initialized: suggested value: 0.2\n";
    
    throw Error::Init();
  }

  if(MasterEquation::excess_energy_over_temperature < 0.) {
    //
    std::cerr << funame << xtot_key << ": excess energy over temperature (upper part of the energy window)"
      //
	      << " has not been initialized: suggested value: 10\n";
    
    throw Error::Init();
  }

  if(MasterEquation::energy_cutoff_over_temperature < 0.) {
    //
    std::cerr << funame << ener_cut_key << ": energy cutoff over temperature (lower part of the energy window)"
      //
	      << " has not been initialized: suggested value: 20\n";
    
    throw Error::Init();
  }

  /*
    if(MasterEquation::well_cutoff < 0.) {
    //
    std::cerr << funame << well_cut_key << ": well cutoff has not been initialized: suggested value: 10\n";

    //MasterEquation::well_cutoff = MasterEquation::energy_cutoff_over_temperature;
    
    throw Error::Init();
    }
  */
    
  if(MasterEquation::reduction_threshold < 0.) {
    //
    std::cerr << funame << well_red_key
      //
	      << ": well reduction threshold has not been initialized: suggested value: 10\n";

    throw Error::Init();
  }
  
  if(MasterEquation::time_propagation_step < 0.) {
    //
    std::cerr << funame << prop_stp_key
      //
	      << ": propagation time-step over fastest relaxation time-scale has not been initialized: suggested value 0.1\n";

    throw Error::Init();
  }
  
  if(MasterEquation::time_propagation_limit < 0.) {
    //
    std::cerr << funame << prop_lim_key
      //
	      << ": propagation time over collisional frequency time-scale has not been initialized: suggested value: 50\n";

    throw Error::Init();
  }
  
  if(MasterEquation::well_projection_threshold < 0.) {
    //
    std::cerr << funame << wpt_key << ": well projection threshold has not been initialized: sugested value: 0.1\n";

    throw Error::Init();
  }
  
  if(MasterEquation::well_projection_tolerance < 0.) {
    //
    std::cerr << funame << wtol_key << ": well projection tolerance has not been initialized: suggested value: 1.e-4\n";

    throw Error::Init();
  }
  
  if(MasterEquation::chemical_threshold < -1.) {
    //
    std::cerr << funame << chem_max_key << ": chemical threshold has not been initialized: sugested value: 0.1\n";

    throw Error::Init();
  }
  
  if(MasterEquation::chemical_tolerance < 0.) {
    //
    std::cerr << funame << chem_min_key << ": chemical tolerance has not been initialized: suggested value: 1.e-12\n";

    throw Error::Init();
  }
  
  if(reduction_scheme.size())
    //
    MasterEquation::set_default_partition(reduction_scheme);

  if(ped_spec.size())
    //
    MasterEquation::set_ped_pair(ped_spec);


  if(MasterEquation::reference_temperature < 0.)
    //
    IO::log << IO::log_offset << "WARNING: reference temperature has not been initialized\n";

  if(MasterEquation::reference_pressure < 0.)
    //
    IO::log << IO::log_offset << "WARNING: reference pressure has not been initialized\n";

  // setting reference lumping scheme
  //
  if(!Model::lump_scheme.size() && !Model::well_exclude_group.size() &&
     //
     MasterEquation::reference_temperature > 0. && MasterEquation::reference_pressure > 0.) {
    //
    IO::log << IO::log_offset << "setting reference lumping scheme:\n";
    
    std::pair<std::list<std::set<int> >, std::set<int> > lumping_scheme =
      //
      MasterEquation::ReactiveComplex::lumping_scheme(MasterEquation::reference_temperature, MasterEquation::reference_pressure);

    IO::log << "\n";
      
    IO::log << IO::log_offset << "reference lumping scheme:";

    for(std::list<std::set<int> >::const_iterator g = lumping_scheme.first.begin(); g != lumping_scheme.first.end(); ++g)
      //
      IO::log << "  " << *g;

    IO::log << "\n";
      
    IO::log << IO::log_offset << "reference exclude group: " << lumping_scheme.second << "\n\n";

    for(std::set<int>::const_iterator i = lumping_scheme.second.begin(); i != lumping_scheme.second.end(); ++i)
      //
      Model::well_exclude_group.insert(Model::well(*i).name());

    std::vector<std::set<int> > wp;

    for(std::list<std::set<int> >::const_iterator i = lumping_scheme.first.begin(); i != lumping_scheme.first.end(); ++i)
      //
      wp.push_back(*i);

    Model::reset(wp);
  }

  if(Model::use_short_names) {
    //
    IO::log << IO::log_offset << "Translation Tables:\n";

    if(Model::well_size()) {
      //
      IO::log << IO::log_offset << "Wells:\n";

      for(int i = 0; i < Model::well_size(); ++i)
	//
	IO::log << IO::log_offset << std::setw(5) << Model::well(i).short_name() << "  " << Model::well(i).name() << "\n";
    }
  
    if(Model::bimolecular_size()) {
      //
      IO::log << IO::log_offset << "Bimolecular:\n";

      for(int i = 0; i < Model::bimolecular_size(); ++i)
	//
	IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(i).short_name() << "  " << Model::bimolecular(i).name() << "\n";
    }
  
    if(Model::inner_barrier_size()) {
      //
      IO::log << IO::log_offset << "Inner Barriers:\n";

      for(int i = 0; i < Model::inner_barrier_size(); ++i)
	//
	IO::log << IO::log_offset << std::setw(5) << Model::inner_barrier(i).short_name() << "  " << Model::inner_barrier(i).name() << "\n";
    }
  
    if(Model::outer_barrier_size()) {
      //
      IO::log << IO::log_offset << "Outer Barriers:\n";

      for(int i = 0; i < Model::outer_barrier_size(); ++i)
	//
	IO::log << IO::log_offset << std::setw(5) << Model::outer_barrier(i).short_name() << "  " << Model::outer_barrier(i).name() << "\n";
    }
  }
  
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

  // default rate output
  //
  if(!IO::out.is_open()) {
    //
    IO::out.open((base_name + ".out").c_str());
  
    if(!IO::out.is_open()) {
      //
      std::cerr << funame << "cannot open rate output file: " << base_name + ".out" << "\n";

      throw Error::Init();
    }
  }

  IO::out << std::setprecision(Model::out_precision);

  Model::pes_print();

  Model::pf_print();

  if(Model::no_run())
    //
    return 0;

  /*********** PRESSURE AND TEMPERATURE DEPENDENT RATE COEFFICIENTS CALCULATION*************/

  // units
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;
  
  // abbreviations
  //
  typedef Model::reac_t                                                                       reac_t;

  typedef std::list<std::pair<std::list<reac_t>, MasterEquation::ReactiveComplex::RateData> > data_t;

  typedef std::map<reac_t, double>                                                            rate_t;

  IO::out << "\n";

  RateData storage;

  // temperature and pressure cycles
  //
  for(std::vector<double>::const_iterator pi = pressure.begin(); pi != pressure.end(); ++pi) {
    //
    std::map<int, double> well_extension_cap;
    
    for(std::vector<double>::const_iterator ti = temperature.begin(); ti != temperature.end(); ++ti) {
      //
      IO::log << IO::log_offset << "Temperature = " << int(*ti / Phys_const::kelv) << " K"  << "   Pressure = ";
  
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

      // barrier - reaction map
      //
      Model::rmap_t barrier_reaction_map = Model::barrier_reaction_map(*ti);

      for(Model::rmap_t::const_iterator ri = barrier_reaction_map.begin(); ri != barrier_reaction_map.end(); ++ri)
	//
	IO::log << IO::log_offset << "reactions controlled by " << ri->first << " barrier:\n" << ri->second;
	
      IO::log << "\n";

      // well lumping
      //
      std::pair<std::list<std::set<int> >, std::set<int> > lumping_scheme =
	//
	MasterEquation::ReactiveComplex::lumping_scheme(*ti, *pi, MasterEquation::ReactiveComplex::NOPRINT);

      //IO::log << "\n";
      
      /*
	IO::log << IO::log_offset << "lumping scheme:";

	for(std::list<std::set<int> >::const_iterator g = lumping_scheme.first.begin(); g != lumping_scheme.first.end(); ++g)
	//
	IO::log << "  " << *g;

	IO::log << "\n";
      
	IO::log << IO::log_offset << "exclude group: " << lumping_scheme.second << "\n\n";
      */
      
      // well extension cap for qualified wells from well lumping scheme
      //
      if(MasterEquation::well_extension >= 0.) {
	//
	std::set<int> mg;
	
	std::vector<std::set<int> > merged(lumping_scheme.first.size());

	for(int w = 0; w < Model::well_size();++w) {
	  //
	  if(Model::well_exclude_group.find(Model::well(w).name()) != Model::well_exclude_group.end())
	    //
	    continue;

	  if(lumping_scheme.second.find(w) != lumping_scheme.second.end() &&  well_extension_cap.find(w) == well_extension_cap.end()) {
	    //
	    well_extension_cap[w] = Model::well(w).thermal_energy(*ti);

	    mg.insert(w);
	    
	    continue;
	  }

	  int count = 0;

	  for(std::list<std::set<int> >::const_iterator i = lumping_scheme.first.begin(); i != lumping_scheme.first.end(); ++i, ++count)
	    //
	    if(i->find(w) != i->end())
	      //
	      merged[count].insert(w);
	}

	for(int g = 0; g < merged.size(); ++g)
	  //
	  if(merged[g].size() > 1)
	    //
	    for(std::set<int>::const_iterator w = merged[g].begin(); w != merged[g].end(); ++w)
	      //
	      if(well_extension_cap.find(*w) == well_extension_cap.end()) {
		//
		well_extension_cap[*w] = Model::well(*w).thermal_energy(*ti);

		mg.insert(*w);
	      }

	if(mg.size()) {
	  //
	  IO::log << IO::log_offset << "adding well extension ";

	  if(mg.size() != 1) {
	    //
	    IO::log << "caps for the ";
	  }
	  else
	    //
	    IO::log << "cap for the ";
	
	  for(std::set<int>::const_iterator w = mg.begin(); w != mg.end(); ++w) {
	    //
	    if(w != mg.begin())
	      //
	      IO::log << ", ";
	    
	    IO::log << Model::well(*w).short_name();
	  }

	  if(mg.size() != 1) {
	    //
	    IO::log << " wells";
	  }
	  else
	    //
	    IO::log << " well";
	  
	  IO::log << "\n\n";
	}
      }

      /*
	IO::log << IO::log_offset << "Temperature = " << int(*ti / Phys_const::kelv) << " K"  << "   Pressure = ";
  
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
      */

      data_t rate_data;

      // barrier-specific rate constant calculation
      //
      for(Model::rmap_t::const_iterator ri = barrier_reaction_map.begin(); ri != barrier_reaction_map.end(); ++ri) {
	//
	const int& type = ri->first.first;
	
	const int& b    = ri->first.second;

	IO::log << IO::log_offset << ri->first << ": rate-limiting barrier for reactions:\n" << ri->second << "\n";

	// dissociation energy minimum
	//
	bool   diss_ener_def = false;
	
	double diss_ener_min;

	if(Model::bimolecular_size()) {
	  //
	  std::vector<double> diss_ener = Model::dissociation_energy_map(*ti);
      
	  btemp = true;
	
	  for(std::list<reac_t>::const_iterator r = ri->second.begin(); r != ri->second.end(); ++r)
	    //
	    for(reac_t::const_iterator s = r->begin(); s != r->end(); ++s)
	      //
	      if(s->first == Model::WELL && (btemp || diss_ener[s->second] < dtemp)) {
		//
		btemp = false;
	      
		dtemp = diss_ener[s->second];
	      }

	  if(btemp)
	    //
	    IO::log << IO::log_offset << "there are no bound reactants for " << ri->first << " barrier\n";
	
	  diss_ener_def = !btemp;
	
	  diss_ener_min = dtemp;
	}

	// control barrier thermal energy
	//
	const double be = Model::thermal_energy(ri->first, *ti);
	
	IO::log << IO::log_offset << ri->first << ": thermal energy = "
	  //
		<< std::ceil(be / Phys_const::kcal * 10.) / 10. << " kcal/mol\n\n";
	
	// smallest reactant dissociation energy
	//
	IO::log << IO::log_offset << "smallest reactant dissociation energy = "
	  //
		<< std::ceil(diss_ener_min / Phys_const::kcal * 10.) / 10. << " kcal/mol\n\n";

	// energy reference
	//
	double er = be + *ti * MasterEquation::excess_energy_over_temperature;

	// energy cutoff
	//
	double ec = be - *ti * MasterEquation::energy_cutoff_over_temperature;
 
	if(MasterEquation::well_cutoff > 0. && diss_ener_def) {
	  //
	  IO::log << IO::log_offset << "hard core energy cutoff = "
	    //
		  << std::ceil(ec / Phys_const::kcal * 10.) / 10. << " kcal/mol\n\n";

	  dtemp = std::min(diss_ener_min, be) - MasterEquation::well_cutoff * *ti;

	  ec = std::max(dtemp, ec);
	}

	Model::ChemGraph cg(er, *ti);

	std::list<Model::ChemGraph> gl = cg.factorize();

	std::list<Model::ChemGraph>::const_iterator li;
	  
	for(li = gl.begin(); li != gl.end(); ++li) {
	  //
	  if(type == Model::INNER && li->inner_set.find(b) != li->inner_set.end())
	    //
	    break;
	  
	  if(type == Model::OUTER && li->outer_set.find(b) != li->outer_set.end())
	    //
	    break;
	}

	if(li == gl.end()) {
	  //
	  std::cerr << funame << "barrier is not in the list\n";

	  throw Error::Logic();
	}

	// actual rate constant calculation
	//
	MasterEquation::ReactiveComplex rc(*ti, *pi, er, ec, *li, well_extension_cap);

	rc.control_barrier = ri->first;
	
	IO::log << "\n";
	
	rate_data.push_back(std::make_pair(ri->second, rc.well_reduction_method()));

	IO::log << "\n";
	//
      }// rate-limiting barrier reaction cycle

      // branching fraction database
      //
      typedef std::map<std::set<int>, std::map<int, double> > bfdb_t;

      bfdb_t bfdb;
      
      for(data_t::const_iterator ri = rate_data.begin(); ri != rate_data.end(); ++ri) {
	//
	const std::set<int> missing = ri->second.missing_bf(ri->first);

	for(std::set<int>::const_iterator g = missing.begin(); g != missing.end(); ++g)
	  //
	  MasterEquation::ReactiveComplex::branching_fraction(*ti, *pi, ri->second.graph, ri->second.well_partition[*g], bfdb);
      }

      IO::log << "\n";

      if(bfdb.size()) {
	//
	IO::log << IO::log_offset << "branching fractions summary:\n";// << std::setprecision(2);

	itemp = 0;
      
	for(bfdb_t::const_iterator dbi = bfdb.begin(); dbi != bfdb.end(); ++dbi) {
	  //
	  stemp.clear();
	
	  for(std::set<int>::const_iterator w = dbi->first.begin(); w != dbi->first.end(); ++w) {
	    //
	    if(w != dbi->first.begin())
	      //
	      stemp += "+";

	    stemp += Model::well(*w).short_name();
	  }

	  if(stemp.size() > itemp)
	    //
	    itemp = stemp.size();
	}

	for(bfdb_t::const_iterator dbi = bfdb.begin(); dbi != bfdb.end(); ++dbi) {
	  //
	  stemp.clear();
	
	  for(std::set<int>::const_iterator w = dbi->first.begin(); w != dbi->first.end(); ++w) {
	    //
	    if(w != dbi->first.begin())
	      //
	      stemp += "+";

	    stemp += Model::well(*w).short_name();
	  }

	  IO::log << IO::log_offset << std::setw(itemp) << stemp << ":";

	  for(std::map<int, double>::const_iterator bfi = dbi->second.begin(); bfi != dbi->second.end(); ++bfi)
	    //
	    IO::log << "  " << bfi->second << "/" << Model::well(bfi->first).short_name();

	  IO::log << "\n";
	}

	IO::log << "\n";// << std::setprecision(6);
      }
      
      // reaction-rate map
      //
      rate_t reaction_map;

      std::set<int> missing;
      
      for(data_t::const_iterator r = rate_data.begin(); r != rate_data.end(); ++r) {
	//
	std::set<int> m = r->second.add_rate_data(reaction_map, r->first, bfdb);

	for(std::set<int>::const_iterator i = m.begin(); i != m.end(); ++i)
	  //
	  missing.insert(*i);
      }
	
      if(missing.size()) {
	//
	IO::log << IO::log_offset << "some rate constants for ";

	for(std::set<int>::const_iterator w = missing.begin(); w != missing.end(); ++w) {
	  //
	  if(w != missing.begin())
	    //
	    IO::log << ", ";
	    
	  IO::log << Model::well(*w).short_name();
	}

	IO::log << " well(s) are missing\n\n";
      }
	  
      storage[std::make_pair(int(ti - temperature.begin()), int(pi - pressure.begin()))] = reaction_map;
      //
    }// temperature cycle
    //
  }// pressure cycle
  
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
  const int count_max = 120;

  // field width
  //
  const int width = Model::out_precision + 7;

  // initial width
  //
  const int first = 7;

  int start, count;

  if(Model::use_short_names) {
    //
    IO::out << "Well Names Translation:\n";

    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::out << std::setw(3) << Model::well(w).short_name() << "  " << Model::well(w).name() << "\n";

    IO::out << "End\n";
  
    IO::out << "Bimolecular Names Translation:\n";

    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::out << std::setw(3) << Model::bimolecular(p).short_name() << "  " << Model::bimolecular(p).name() << "\n";

    IO::out << "End\n";
  }
  
  IO::out << "Unimolecular Rate Units: 1/sec;  Bimolecular Rate Units: cm^3/sec\n\n"
    //
	  << "______________________________________________________________________________________\n\n";

  IO::out << "Capture rates [cm^3/sec]:\n\n";

  IO::out << std::setw(first) << "T, K";

  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
    //
    const int& p = Model::outer_connect(b).second;
      
    if(!Model::bimolecular(p).dummy())
      //
      IO::out << std::setw(width) << Model::bimolecular(p).short_name();
  }
  
  IO::out << "\n";
      
  for(int t = 0; t < temperature.size(); ++t) {
    //
    IO::out << std::setw(first) << temperature[t] / Phys_const::kelv;

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int& p = Model::outer_connect(b).second;
      
      if(!Model::bimolecular(p).dummy())
	//
	IO::out << std::setw(width) << temperature[t] / 2. / M_PI / bru * Model::outer_barrier(b).weight(temperature[t])
	  //
	  / Model::bimolecular(p).weight(temperature[t]) * std::exp((Model::bimolecular(p).ground() - Model::outer_barrier(b).ground()) / temperature[t]);
    }

    IO::out  << "\n";
  }

  IO::out << "\n";
    
  IO::out << "High Pressure rates [1/sec]:\n\n";

  IO::out << std::setw(first) << "T, K";

  for(int b = 0; b < Model::inner_barrier_size(); ++b) {
    //
    const int& w1 = Model::inner_connect(b).first;
      
    const int& w2 = Model::inner_connect(b).second;
      
    IO::out << std::setw(width) << std::make_pair(Reactant(Model::WELL, w1), Reactant(Model::WELL, w2));
    
    IO::out << std::setw(width) << std::make_pair(Reactant(Model::WELL, w2), Reactant(Model::WELL, w1));
  }
  
  IO::out << "\n";
      
  for(int t = 0; t < temperature.size(); ++t) {
    //
    IO::out << std::setw(first) << temperature[t] / Phys_const::kelv;

    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      const int& w1 = Model::inner_connect(b).first;
      
      const int& w2 = Model::inner_connect(b).second;
      
      IO::out << std::setw(width) << temperature[t] / 2. / M_PI / Phys_const::herz * Model::inner_barrier(b).weight(temperature[t])
	//
	/ Model::well(w1).weight(temperature[t]) * std::exp((Model::well(w1).ground() - Model::inner_barrier(b).ground()) / temperature[t]);

      IO::out << std::setw(width) << temperature[t] / 2. / M_PI / Phys_const::herz * Model::inner_barrier(b).weight(temperature[t])
	//
	/ Model::well(w2).weight(temperature[t]) * std::exp((Model::well(w2).ground() - Model::inner_barrier(b).ground()) / temperature[t]);
    }

    IO::out  << "\n";
  }

  IO::out << "\n";
    

  IO::out<< "Species-Species Rate Tables:\n\n";

  for(int p = 0; p < pressure.size(); ++p) {
    //
    for(int t = 0; t < temperature.size(); ++t) {
      //
      IO::out << "Temperature = " << int(temperature[t] / Phys_const::kelv) <<  " K    Pressure = ";
      
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
	
      IO::out << "Reactant = " << reactant[r] << "   Pressure = ";
      
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
      
      IO::out << "\n\n";

      start = 0;

      while(start < reactant.size()) {
	//
	int i;
	
	IO::out << std::left << std::setw(first) << "T" << std::right;

	for(i = start, count = first; count < count_max && i < reactant.size(); ++i, count += width)
	  //
	  if(i != r)
	    //
	    IO::out << std::setw(width) << reactant[i];

	IO::out << "\n";
	    
	for(int t = 0; t < temperature.size(); ++t) {
	  //
	  //
	  const double weight = reactant[r].weight(temperature[t]);

	  IO::out << std::left << std::setw(first) << int(temperature[t] / Phys_const::kelv) << std::right;
	
	  for(i = start, count = first; count < count_max && i < reactant.size() ; ++i, count += width)
	    //
	    if(i != r)
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

      IO::out << "Reactant = " << reactant[r] << "   Temperature = " << int(temperature[t] / Phys_const::kelv) << " K\n\n";
      
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
	
	for(i = start, count = first; count < count_max && i < reactant.size(); ++i, count += width)
	  //
	  if(i != r)
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
	
	  for(i = start, count = first; count < count_max && i < reactant.size() ; ++i, count += width)
	    //
	    if(i != r)
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
