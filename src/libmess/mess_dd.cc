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

#include <iostream>
#include<ios>
#include <iomanip>
#include<cmath>
#include <list>
#include <map>
#include <ctime>

#include "mess.hh"
#include "units.hh"
#include "key.hh"
#include "io.hh"
#include "shared.hh"

#ifdef WITH_MPACK

#include "mpack.hh"

#endif

namespace MasterEquation {
  /*********************************** AUXILIARY OUTPUT *************************************/

  std::ofstream eval_out;// eigenvalues output
  int  red_out_num = 5; // number of reduction schemes to print 

  std::ofstream evec_out;
  int evec_out_num = 0;

  std::ofstream arr_out; // arrhenius 

  /********************************* USER DEFINED PARAMETERS ********************************/


#ifdef WITH_MPACK

  typedef dd_real dd_t;

#else

  typedef double dd_t;

#endif
  
  bool well_reduction_correction = true;

  int float_type = DOUBLE;
  
  // temperature
  //
  double                                                     _temperature     = -1.;

  // pressure
  //
  double                                                     _pressure        = -1.;

  // energy grid step
  //
  double                                                     _energy_step     = -1.;
  
  // zero energy (also upper limit)
  //
  double                                                     _energy_reference;
  
  // well depth cutoff parameter
  //
  double                                                     well_cutoff       = -1.;
  
  // highest chemecial eigenvalue
  //
  double                                                     chemical_threshold     = 0.;
  
  // smallest chemical eigenvalue
  //
  double                                                     min_chem_eval     = -1.;
  
  // maximal chemical relaxation eigenvalue to collision frequency ratio
  //
  double                                                     reduction_threshold    = 10.;
  
  // species reduction algorithm for low-eigenvalue method
  //
  int                                                        reduction_method  = DIAGONALIZATION;
  
  // default reduction scheme
  //
  Partition                                                  default_partition;
  
  // default chemical subspace size
  //
  int                                                        _default_chem_size = -1;

  // microcanonical rate maximum
  //
  double                                                     rate_max          = -1.;
  
  // global energy cutoff
  //
  double                                                     global_cutoff;
  
  bool                                                       is_global_cutoff = false;
  
  void set_global_cutoff (double e) { global_cutoff = e; is_global_cutoff = true; }

  // pressure units
  //
  int                                                        pressure_unit = BAR;

  // well partition method
  //
  double (*well_partition_method) (const Lapack::Matrix&, Partition&, Group&, const std::vector<double>&) = threshold_well_partition;

  // well partition threshold
  //
  double                                                     well_projection_threshold = 0.2;

  // well extention (in temperature units)
  //
  double                                                     well_extension = -1.;

  // well extension correction
  //
  double                                                     well_ext_corr = 0.;
  
  /********************************* INTERNAL PARAMETERS ************************************/

  // inner setting check
  bool  _isset = false;

  // collisional frequency
  //double _collision_frequency_factor;
  //double _collision_frequency;
  //std::vector<double> _kernel_fraction;

  // numerical accuracy
  const double epsilon = 1.e-10;

  /*******************************************************************************************/

  std::vector<double> _thermal_factor;

  std::vector<SharedPointer<Well> >  _well;
  const Well& well (int w) { return *_well[w]; } 

  std::vector<SharedPointer<Bimolecular> >  _bimolecular;
  const Bimolecular& bimolecular (int p) { return *_bimolecular[p]; }

  std::vector<SharedPointer<Barrier> > _inner_barrier;
  const Barrier& inner_barrier (int p) { return *_inner_barrier[p]; }

  std::vector<SharedPointer<Barrier> > _outer_barrier;
  const Barrier& outer_barrier (int p) { return *_outer_barrier[p]; }
  
  // Boltzmann factor: exp(-E/T)
  double        thermal_factor (int); 
  void   resize_thermal_factor (int);
  void    reset_thermal_factor (int = 0);

  bool isset () { return _isset; }
  
  double energy_reference () { return _energy_reference; }

  double energy_step ()
  {
    const char funame [] = "MasterEquation::energy_step: ";
    
    if(_energy_step <= 0.) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }
    
    return _energy_step;
  }

  double energy_bin (int i) { return energy_reference() - i * energy_step(); }
  
  double temperature ()
  {
    const char funame [] = "MasterEquation::temperature: ";
    
    if(_temperature <= 0.) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }
    
    return _temperature;
  }
  
  double pressure                ()
  {
    const char funame [] = "MasterEquation::pressure: ";
    
    if(_pressure <= 0.) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }
    
    return _pressure;
  }
  //double collision_frequency     () { return _collision_frequency;    }
  //double kernel_fraction    (int i) { return _kernel_fraction[i];     }

  // cumulative number of states for each well
  //
  std::vector<Array<double> >  cum_stat_num;

  // capture probabilities
  //
  std::map<std::string, std::set<double> > hot_energy;
  
  std::set<std::pair<int, int> > hot_index;

  // product energy distributions
  //
  std::ofstream  ped_out;// product energy distribution output stream
  
  std::vector<std::pair<int, int> > ped_pair; // product energy distribution reactants and products indices
}

void MasterEquation::set_temperature (double t) 
{ 
  double dtemp;

  _isset  = false;
  _temperature = t;

  reset_thermal_factor();
  /*
  _collision_frequency_factor = 0.;
  _kernel_fraction.resize(Model::buffer_size());
  for(int b = 0; b < Model::buffer_size(); ++b) {
    dtemp = (*Model::collision(b))(t) * Model::buffer_fraction(b);
    _kernel_fraction[b] = dtemp;
    _collision_frequency_factor += dtemp;
  }

  for(int b = 0; b < Model::buffer_size(); ++b)
    _kernel_fraction[b] /= _collision_frequency_factor;

  */
  
}

void MasterEquation::set_energy_step (double e) 
{ 
  _isset = false; 
  _energy_step = e; 
  reset_thermal_factor(); 
}

void MasterEquation::set_energy_reference (double e) 
{
  _isset = false;
  _energy_reference = e;
}

void MasterEquation::set_pressure (double p) 
{
  _pressure = p;
  //_collision_frequency = p * _collision_frequency_factor; 
}

/********************************************************************************************
 ************************************* THERMAL FACTOR ***************************************
 ********************************************************************************************/

double MasterEquation::thermal_factor(int e)
{
  if(e >= _thermal_factor.size())
    resize_thermal_factor(e + 1);
    return _thermal_factor[e];
}

void  MasterEquation::resize_thermal_factor(int s)
{
  if(s <= _thermal_factor.size())
    return;

  int emin = _thermal_factor.size();
  _thermal_factor.resize(s);

  double x = energy_step() / temperature();
  for(int e = emin; e < _thermal_factor.size(); ++e)
    _thermal_factor[e] = std::exp(double(e) * x);
}

void  MasterEquation::reset_thermal_factor(int s)
{
  _thermal_factor.resize(s);
  if(!s)
    return;

  double x = energy_step() / temperature();
  for(int e = 0; e < _thermal_factor.size(); ++e)
    _thermal_factor[e] = std::exp(double(e) * x);
}

/************************************* PRODUCT ENERGY DISTRIBUTION PAIRS ********************************************/

void MasterEquation::set_ped_pair(const std::vector<std::string>& ped_spec) 
{
  const char funame [] = "MasterEquation::set_ped_pair: ";

  IO::Marker funame_marker(funame, IO::Marker::NOTIME);

  std::string::size_type itemp, start, next;
  int spec;

  for(std::vector<std::string>::const_iterator ped = ped_spec.begin(); ped != ped_spec.end(); ++ped) {
    std::vector<int> pair;
    start = 0;
    while(start < ped->size()) {
      next = ped->size();
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	itemp = ped->find(Model::bimolecular(p).name(), start);
	if(itemp < next || (itemp == next && itemp < ped->size() &&
			    Model::bimolecular(p).name().size() > Model::bimolecular(spec).name().size())) {
	  next = itemp;
	  spec = p;
	}
      }
      if(next >= ped->size())
	break;
      pair.push_back(spec);
      start = next + Model::bimolecular(spec).name().size();
    }
    if(pair.size() != 2) {
      std::cerr << funame << "wrong number of bimolecular species (" << pair.size() << ") in " << *ped << " token\n";
      throw Error::Init();
    }
    ped_pair.push_back(std::make_pair(pair.front(), pair.back()));
  }

  if(!ped_pair.size()) {
    std::cerr << funame << "no pairs\n";
    throw Error::Init();
  }
  
  IO::log << IO::log_offset << "ped pairs:\n";
  for(int i = 0; i < ped_pair.size(); ++i)
    IO::log << IO::log_offset << "   " 
	    << std::setw(5) << Model::bimolecular(ped_pair[i].first).name()
	    << " -> " 
	    << Model::bimolecular(ped_pair[i].second).name()
	    << "\n";
}

/***************************************** DEFAULT REDUCTION SCHEME *************************************************/

void MasterEquation::set_default_partition(const std::vector<std::string>& scheme) 
{
  const char funame [] = "MasterEquation::set_default_partition: ";

  IO::Marker funame_marker(funame, IO::Marker::NOTIME);

  std::string::size_type itemp, start, next;
  int spec;

  std::set<int> pool;
  for(std::vector<std::string>::const_iterator g = scheme.begin(); g != scheme.end(); ++g) {
    Group group;
    start = 0;
    while(start < g->size()) {
      next = g->size();
      for(int w = 0; w < Model::well_size(); ++w) {
	itemp = g->find(Model::well(w).name(), start);
	if(itemp < next || (itemp == next && itemp < g->size() &&
			    Model::well(w).name().size() > Model::well(spec).name().size())) {
	  next = itemp;
	  spec = w;
	}
      }
      if(next >= g->size())
	break;
      
      if(!pool.insert(spec).second) {
	std::cerr << funame << "duplicated well name: " << Model::well(spec).name() << "\n";
	throw Error::Init();
      }
      
      group.insert(spec);      
      start = next + Model::well(spec).name().size();
    }

    if(!group.size()) {
      std::cerr << funame << "empty " << *g << " group\n";
      throw Error::Init();
    }
    default_partition.push_back(group);
  }
  if(default_partition.size() >= Model::well_size()) {
    std::cerr << funame << "default reduction scheme: no reduction\n";
    throw Error::Init();   
  }

  // output
  IO::log << IO::log_offset << "default reduction scheme: ";
  for(Pit g = default_partition.begin(); g != default_partition.end(); ++g) {
    if(g != default_partition.begin())
      IO::log << " ";
    for(Git w = g->begin(); w != g->end(); ++w) {
      if(w != g->begin())
	IO::log << "+";
      IO::log << Model::well(*w).name();
    }
  }
  IO::log << "\n";
}

void MasterEquation::set_default_chem_size (int s)
{
  if(s < 0) {
    std::cerr << "MasterEquation::set_default_chem_size: out of range\n";
    throw Error::Range();
  }

  _default_chem_size = s;
}

/********************************************************************************************
 ********************* SETTING WELLS, BARRIERS, AND BIMOLECULAR *****************************
 ********************************************************************************************/

void MasterEquation::set (std::map<std::pair<int, int>, double>& rate_data, std::map<int, double>& capture, int flags)
  
{
  const char funame [] = "MasterEquation::set: ";

  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  double dtemp;
  int    itemp;
  bool   btemp;
  
  IO::Marker funame_marker(funame);

  // set default energy reference
  //
  if(flags & DEFAULT_EREF) {

    double e, s, c;
    
    btemp = true;

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      Model::outer_barrier(b).esc_parameters(temperature(), e, s, c);

      IO::log << IO::log_offset << Model::outer_barrier(b).name() << ": e = " << e / Phys_const::kcal << " kcal/mol, c = " << c << std::endl;
      
      e += Model::outer_barrier(b).ground() + 5. * temperature() * std::sqrt(c);
      
      if(btemp || _energy_reference < e) {
	//
	btemp = false;
	
	_energy_reference = e;
      }
    }
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      Model::inner_barrier(b).esc_parameters(temperature(), e, s, c);
      
      IO::log << IO::log_offset << Model::inner_barrier(b).name() << ": e = " << e / Phys_const::kcal << " kcal/mol, c = " << c << std::endl;
      
      e += Model::inner_barrier(b).ground() + 5. * temperature() * std::sqrt(c);
      
      if(btemp || _energy_reference < e) {
	//
	btemp = false;
	
	_energy_reference = e;
      }
    }

    //_energy_reference = nearbyint(_energy_reference / Phys_const::incm) * Phys_const::incm;
  }

  if(energy_reference() > Model::energy_limit()) {
    //
    std::cerr << funame << "model energy limit (" << Model::energy_limit() / Phys_const::kcal
      //
	      << " kcal/mol) is lower than the energy reference (" << energy_reference() / Phys_const::kcal
      //
	      << " kcal/mol)\n";
    
    throw Error::Range();
  }

  rate_data.clear();


  _isset = true;

  IO::log << IO::log_offset << "Temperature      = " << temperature()      / Phys_const::kelv << " K\n";
  //IO::log << IO::log_offset << "Collision rate   = " << _collision_frequency_factor * temperature() / bru << "cm^3/sec\n";
  IO::log << IO::log_offset << "Energy reference = " << std::ceil(energy_reference() / Phys_const::incm) << " 1/cm\n";
  IO::log << IO::log_offset << "Energy step      = " << energy_step()      / Phys_const::incm << " 1/cm\n";

  {
    IO::Marker set_marker("setting wells, barriers, and bimolecular");

    _well.resize(Model::well_size());
    // wells
    for(int w = 0; w < _well.size(); ++w)
      _well[w] = SharedPointer<Well>(new Well(Model::well(w)));

    // well-to-well barriers
    _inner_barrier.resize(Model::inner_barrier_size());
    for(int b = 0; b < _inner_barrier.size(); ++b)
      _inner_barrier[b] = SharedPointer<Barrier>(new Barrier(Model::inner_barrier(b)));

    // well-to-bimolecular barriers
    _outer_barrier.resize(Model::outer_barrier_size());
    for(int b = 0; b < _outer_barrier.size(); ++b)
      _outer_barrier[b] = SharedPointer<Barrier>(new Barrier(Model::outer_barrier(b)));

    // bimolecular products
    _bimolecular.resize(Model::bimolecular_size());
    for(int p = 0; p < _bimolecular.size(); ++p)
      _bimolecular[p] = SharedPointer<Bimolecular>(new Bimolecular(Model::bimolecular(p)));

  }

  // checking if wells are deeper than the barriers between them
  
  // inner barriers
  //
  for(int b = 0; b < Model::inner_barrier_size(); ++b) {
    itemp = well(Model::inner_connect(b).first).size() <  well(Model::inner_connect(b).second).size() ?
      well(Model::inner_connect(b).first).size() : well(Model::inner_connect(b).second).size();
    if(itemp  < inner_barrier(b).size()) {
      IO::log << IO::log_offset << Model::inner_barrier(b).name() 
	      << " barrier top is lower than the bottom of one of the wells it connects => truncating\n"; 
      _inner_barrier[b]->truncate(itemp);
    }
  }
  
  // outer barriers
  //
  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
    itemp = well(Model::outer_connect(b).first).size(); 
    if(itemp  < outer_barrier(b).size()) {
      IO::log << IO::log_offset << Model::outer_barrier(b).name()
	      << " barrier top is lower than the bottom of the well it connects => truncating\n"; 
      _outer_barrier[b]->truncate(itemp);
    }
  }

  // clipping the number of states with the maximum rate constant
  
  // inner barriers
  //
  if(rate_max > 0.) {
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      int w1 = Model::inner_connect(b).first;
      int w2 = Model::inner_connect(b).second;
      for(int i = 0; i < _inner_barrier[b]->size(); ++i) {
	dtemp = well(w1).state_density(i) < well(w2).state_density(i) ? 
	  well(w1).state_density(i) : well(w2).state_density(i);
	dtemp *= rate_max * 2. * M_PI;
	
	if(_inner_barrier[b]->state_number(i) > dtemp)
	  _inner_barrier[b]->state_number(i) = dtemp;
      }
    }
    
    // outer barriers
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      int w = Model::outer_connect(b).first;
      for(int i = 0; i < _outer_barrier[b]->size(); ++i) {
	dtemp = rate_max * 2. * M_PI * well(w).state_density(i);

	if(_outer_barrier[b]->state_number(i) > dtemp)
	  _outer_barrier[b]->state_number(i) = dtemp;
      }
    }
  }

  // cumulative number of states for each well
  //
  cum_stat_num.resize(Model::well_size());
  
  for(int w = 0; w < Model::well_size(); ++w) {// well cycle
    //
    itemp = 0;
    
    for(int b = 0; b < Model::inner_barrier_size(); ++b)
      //
      if(Model::inner_connect(b).first == w || Model::inner_connect(b).second == w)
	//
	itemp = inner_barrier(b).size() > itemp ? inner_barrier(b).size() : itemp;
    
    for(int b = 0; b < Model::outer_barrier_size(); ++b)
      //
      if(Model::outer_connect(b).first == w)
	//
	itemp = outer_barrier(b).size() > itemp ? outer_barrier(b).size() : itemp;

    cum_stat_num[w].resize(itemp);
    
    cum_stat_num[w] = 0.;
    
    for(int b = 0; b < Model::inner_barrier_size(); ++b)
      //
      if(Model::inner_connect(b).first == w || Model::inner_connect(b).second == w)
	//
	for(int i = 0; i < inner_barrier(b).size(); ++i)
	  //
	  cum_stat_num[w][i] += inner_barrier(b).state_number(i);
    
    for(int b = 0; b < Model::outer_barrier_size(); ++b)
      //
      if(Model::outer_connect(b).first == w)
	//
	for(int i = 0; i < outer_barrier(b).size(); ++i)
	  //
	  cum_stat_num[w][i] += outer_barrier(b).state_number(i);
    
  }// well cycle


  /************************************** OUTPUT ***************************************/

  if(evec_out.is_open()) {
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      if(!w || well(w).size() > itemp)
	//
	itemp = well(w).size();
    
    const int well_size_max = itemp;
    
    evec_out << "Temperature = " << temperature() / Phys_const::kelv << "K\n"
      //
	     <<"THERMAL DISTRIBUTIONS:\n"
      //
	     << std::setw(13) << "E, kcal/mol";
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      evec_out << std::setw(13) << Model::well(w).name();
    
    evec_out << "\n";
    
    for(int i = 0; i < well_size_max; ++i) {
      //
      evec_out << std::setw(13) << energy_bin(i) / Phys_const::kcal;
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	if(i < well(w).size()) {
	  //
	  evec_out << std::setw(13) << well(w).state_density(i) * thermal_factor(i) / well(w).weight_sqrt();
	}
	else
	  //
	  evec_out << std::setw(13) << 0;
      
      evec_out << "\n";
    }
    evec_out << "\n";    
  }

  IO::log << std::setprecision(3);

  if(Model::well_size()) {
    //
    IO::log << IO::log_offset << "Bound Species (D1/D2 - density of states at dissociation energy (DE)/reference energy):\n"
      //
	    << IO::log_offset << std::setw(5)  << "Name" << std::setw(7)  << "Depth" << std::setw(7)  << "DE"
      //
	    << std::setw(9) << "*D1" << std::setw(9) << "*D2"  << "\n"
      //
	    << IO::log_offset  << std::setw(5) << "" << std::setw(7)  << "1/cm" << std::setw(7)  << "1/cm"
      //
	    << std::setw(9) << "cm" << std::setw(9) << "cm"  << "\n";

    for(int w = 0; w < Model::well_size(); ++w) {
      //
      IO::log << IO::log_offset
	//
	      << std::setw(5) << Model::well(w).name();
      
      dtemp = energy_bin(well(w).size());
      
      IO::log << std::setw(7) << (int)std::ceil(dtemp / Phys_const::incm);
      
      dtemp = energy_bin(cum_stat_num[w].size());
      
      IO::log << std::setw(7) << (int)std::ceil(dtemp  / Phys_const::incm);

      if(cum_stat_num[w].size() != well(w).size()) {
	//
	itemp = cum_stat_num[w].size();
      }
      else
	//
	itemp = well(w).size() - 1;

      IO::log << std::setw(9) << well(w).state_density(itemp) * Phys_const::incm
	//
	      << std::setw(9) << well(w).state_density(0) * Phys_const::incm
	      << "\n";
    }
  }

  if(Model::inner_barrier_size()) {
    //
    IO::log << IO::log_offset << "Well-to-Well Barriers (N - number of states at the reference energy):\n"
      //
	    << IO::log_offset << std::setw(5)  << "Name" << std::setw(7)  << "Height" << std::setw(9) << "*N" << "\n"
      //
	    << IO::log_offset << std::setw(5)  << ""     << std::setw(7)  << "1/cm"    << std::setw(9) << ""  << "\n";

    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::inner_barrier(b).name();
      
      dtemp = energy_bin(inner_barrier(b).size());
      
      IO::log << std::setw(7) << (int)std::ceil(dtemp / Phys_const::incm)
	//
	      << std::setw(9) << inner_barrier(b).state_number(0)
	//
	      << "\n";
    }
  }

  if(Model::outer_barrier_size()) {
    //
    IO::log << IO::log_offset
      //
	    << "Well-to-Bimolecular Barriers (N - number of states at the reference energy):\n"
      //
	    << IO::log_offset
      //
	    << std::setw(5)  << "Name"
      //
	    << std::setw(7)  << "Height"
      //
	    << std::setw(9) << "*N"
      //
	    << "\n"
      //
	    << IO::log_offset
      //
	    << std::setw(5)  << ""
      //
	    << std::setw(7)  << "1/cm"
      //
	    << std::setw(9) << ""
      //
	    << "\n";

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::outer_barrier(b).name();
      //
      dtemp = energy_bin(outer_barrier(b).size());
      
      IO::log << std::setw(7) <<  (int)std::ceil(dtemp / Phys_const::incm)
	//
	      << std::setw(9) << outer_barrier(b).state_number(0)
	//
	      << "\n";
    }
  }

  // bimolecular partition function units
  //
  const double bpu = Phys_const::cm * Phys_const::cm * Phys_const::cm;

  IO::log << IO::log_offset << "Effective equilibrium constants(bimolecular units - cm^3):\n"
    //
	  << IO::log_offset << std::setw(5) << "Q/Q";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    if(!Model::bimolecular(p).dummy())
      //
      IO::log << std::setw(10) << Model::bimolecular(p).name();
  
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(5) << Model::well(w).name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1)
      //
      IO::log << std::setw(10) << well(w).weight() / well(w1).weight();
    
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy())
	//
	IO::log << std::setw(10) << well(w).weight() / bimolecular(p1).weight() / bpu;
    
    IO::log << "\n";
  }
  
  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    if(Model::bimolecular(p).dummy())
      //
      continue;
    
    IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1)
      //
      IO::log << std::setw(10) << bimolecular(p).weight() * bpu / well(w1).weight();
    
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy())
	//
	IO::log << std::setw(10) << bimolecular(p).weight() / bimolecular(p1).weight();
    
    IO::log << "\n";
  }
  
  IO::log << IO::log_offset << "Real equilibrium constants(bimolecular units - cm^3):\n"
    //
	  << IO::log_offset << std::setw(5) << "Q/Q";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    if(!Model::bimolecular(p).dummy())
      //
      IO::log << std::setw(10) << Model::bimolecular(p).name();
  
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(5) << Model::well(w).name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1)
      //
      IO::log << std::setw(10) << well(w).real_weight() / well(w1).real_weight();
    
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy())
	//
	IO::log << std::setw(10) << well(w).real_weight() / bimolecular(p1).weight() / bpu;
    
    IO::log << "\n";
  }
  
  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    if(Model::bimolecular(p).dummy())
      //
      continue;
    
    IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1)
      //
      IO::log << std::setw(10) << bimolecular(p).weight() * bpu / well(w1).real_weight();
    
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy())
	//
	IO::log << std::setw(10) << bimolecular(p).weight() / bimolecular(p1).weight();
    
    IO::log << "\n";
  }
  
  IO::log << std::setprecision(6);

  // High pressure well-to-well rate coefficients
  //
  if(Model::well_size()) {
    //
    Lapack::Matrix ww_rate(Model::well_size());
    
    ww_rate = 0.;
    
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      const int& w1 = Model::inner_connect(b).first;
      
      const int& w2 = Model::inner_connect(b).second;
      
      rate_data[std::make_pair(w1, w2)] = temperature() / 2. / M_PI / Phys_const::herz
	//
	* inner_barrier(b).real_weight() / well(w1).real_weight();

      rate_data[std::make_pair(w2, w1)] = temperature() / 2. / M_PI / Phys_const::herz
	//
	* inner_barrier(b).real_weight() / well(w2).real_weight();
    }
  }

  // High pressure well-to-bimolecular and bimolecular-to-well rate coefficients
  //
  if(Model::well_size() && Model::bimolecular_size()) {
    //
    // well-to-bimolecular
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int& w = Model::outer_connect(b).first;
      
      rate_data[std::make_pair(w, Model::outer_connect(b).second + Model::well_size())] =
	//
	temperature() / 2. / M_PI / Phys_const::herz * outer_barrier(b).real_weight() / well(w).real_weight();
    }

    // bimolecular-to-well
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int& p = Model::outer_connect(b).second;
      
      if(bimolecular(p).weight() > 0.)
	//
	rate_data[std::make_pair(p + Model::well_size(), Model::outer_connect(b).first)] =
	  //
	  temperature() / 2. / M_PI / bru * outer_barrier(b).real_weight() / bimolecular(p).weight();
    }
  }

  // capture
  //
  capture.clear();
  
  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
    //
    const int& w = Model::outer_connect(b).first;
    
    const int& p = Model::outer_connect(b).second;
    
    if(bimolecular(p).weight() > 0.)
      //
      capture[p + Model::well_size()] += temperature() / 2. / M_PI / bru
	//
	* outer_barrier(b).real_weight() / bimolecular(p).weight();
    
    capture[w] += temperature() / 2. / M_PI / Phys_const::herz
      //
      * outer_barrier(b).real_weight() / well(w).real_weight();
  }

  for(int b = 0; b < Model::inner_barrier_size(); ++b) {
    //
    const int& w1 = Model::inner_connect(b).first;
    
    const int& w2 = Model::inner_connect(b).second;
    
    capture[w1] += temperature() / 2. / M_PI / Phys_const::herz
      //
      * inner_barrier(b).real_weight() / well(w1).real_weight();
    
    capture[w2] += temperature() / 2. / M_PI / Phys_const::herz
      //
      * inner_barrier(b).real_weight() / well(w2).real_weight();
  }

  // hot energies
  //
  if(hot_energy.size()) {
    //
    hot_index.clear();
    
    for(std::map<std::string, std::set<double> >::const_iterator hit = hot_energy.begin(); hit != hot_energy.end(); ++hit) {
      //
      const int w = Model::well_by_name(hit->first);
      
      // check if hot energy wells exist
      //
      if(w < 0) {
	//
	std::cerr << funame << hit->first << " well does not exist\n";

	throw Error::Init();
      }
      
      for(std::set<double>::const_iterator it = hit->second.begin(); it != hit->second.end(); ++it) {
	//
	itemp = (energy_reference() - *it) / energy_step();
	  
	if(itemp >= 0 && itemp < well(w).size())
	  //
	  hot_index.insert(std::make_pair(w, itemp));
      }//
      //
    }//
    //
  }// hot energies
  //
}

/********************************************************************************************
 ************************************** SETTING WELL ****************************************
 ********************************************************************************************/

// set state density, relaxation mode basis, collisional energy transfer kernal,
// partition function, etc.

void  MasterEquation::Well::_set_state_density (const Model::Well& model)
{
  const char funame [] = "MasterEquation::Well::_set_state_density: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  // enthalpy (average energy), entropy, and thermal capacity
  //
  double hval, sval, cval;

  model.esc_parameters(temperature(), hval, sval, cval);

  hval += model.ground();
  
  IO::log << IO::log_offset << model.name() << " Well: average energy = "
    //
	  << std::ceil(hval / Phys_const::kcal * 10.) / 10.
    //
	  << " kcal/mol\n";
  
  /***************************** SETTING DENSITY OF STATES ************************************/

  double ener, base_ener, ext_ener = -1.;
  
  if(well_extension >= 0. && well_cutoff > 0.) {
    //
    //dtemp = hval > model.dissociation_limit ? hval : model.dissociation_limit;

    dtemp = hval;
    
    base_ener = dtemp * (1. - well_extension) + model.ground() * well_extension;

    // well extension cap
    //
    if(model.well_ext_cap > 0. && model.well_ext_cap < base_ener)
      //
      base_ener = model.well_ext_cap + well_ext_corr * (base_ener - model.well_ext_cap);
     
    ext_ener = well_cutoff * temperature();
	
    if(is_global_cutoff) {
      //
      if(global_cutoff >= base_ener) {
	//
	base_ener = global_cutoff;

	ext_ener = -1.;
      }
      else {
	//
	dtemp = base_ener - global_cutoff;

	ext_ener = dtemp < ext_ener ? dtemp : ext_ener;
      }
    }
  }
  //
  else {
    //
    if(well_cutoff > 0.)
      //
      base_ener = model.dissociation_limit - well_cutoff * temperature();

    if(well_cutoff <= 0. || model.ground() > base_ener)
      //
      base_ener = model.ground();

    if(is_global_cutoff && global_cutoff > base_ener)
      //
      base_ener = global_cutoff;
  }

  int base_size = std::ceil((energy_reference() - base_ener) / energy_step());

  int ext_size = ext_ener < 0. ? 0 : std::ceil(ext_ener / energy_step());
  
  if(ext_size)
    //
    IO::log << IO::log_offset << model.name() << " Well: extension size = " << ext_size << "\n";
											   
  itemp = base_size + ext_size;
  
  _state_density.resize(itemp);
  
  resize_thermal_factor(itemp);

  ener = energy_reference();
  
  for(int e = 0; e < base_size; ++e, ener -= energy_step()) {
    //
    dtemp = model.states(ener);
    
    if(dtemp <= 0.) {
      //
      IO::log << IO::log_offset << model.name()  << " Well: WARNING: nonpositive density at " 
	      << ener / Phys_const::incm << " 1/cm, truncating\n";

      _state_density.resize(e);      

      return;    
    }
    
    _state_density[e] = dtemp;
  }

  //double fac = std::exp(energy_step() / temperature() / 2.);
  
  for(int e = base_size; e < size(); ++e)
    //
    _state_density[e] = dtemp;
}

void MasterEquation::Well::_set_kernel (const Model::Well& model) 
{
  const char funame [] = "MasterEquation::Well::_set_kernel: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  double a, c;

  btemp = false;

  do {
    //
    _kernel.resize(size());
    
    _kernel = 0.;

    Lapack::Matrix tmp_kernel(size());

    for(int b = 0; b < Model::buffer_size(); ++b) {

      /********************* SETTING COLLISIONAL ENERGY TRANSFER KERNEL ****************************/

      tmp_kernel = 0.;

      // collisional energy transfer on the grid
      //
      itemp = (int)std::ceil(model.kernel(b)->cutoff_energy(temperature()) / energy_step());

      IO::log << IO::log_offset << model.name() << " Well: collisional kernel bandwidth = " << itemp << "\n";

      std::vector<double> energy_transfer(itemp);
	
      for(int i = 0; i < energy_transfer.size(); ++i)
	//
	energy_transfer[i] = (*model.kernel(b))((double)i * energy_step(), temperature());
      
      if(!b || kernel_bandwidth < energy_transfer.size())
	//
	kernel_bandwidth = energy_transfer.size();

      // energy transfer UP probability functional form predefined
      //
      if(Model::Kernel::flags() & Model::Kernel::UP) {
	//
	for(int i = size() - 1; i >= 0; --i) {// energy grid cycle
	  //
	  itemp = i + energy_transfer.size();
	  
	  const int jmax = itemp < size() ? itemp : size();
	  
	  itemp = i - energy_transfer.size() + 1;
	  
	  const int jmin = itemp  > 0 ? itemp : 0;

	  // normalization constant
	  //
	  c = 0.;
	  
	  for(int j = jmin; j <= i; ++j) {
	    //
	    dtemp = energy_transfer[i - j];
	    
	    if(Model::Kernel::flags() & Model::Kernel::DENSITY)
	      //
	      dtemp *= state_density(j);
	    
	    tmp_kernel(i, j) = -dtemp;
	    
	    c += dtemp;
	  }

	  // down-transitions contribution
	  //
	  a = kernel_fraction(b);
	  
	  for(int j = i + 1; j < jmax; ++j) {
	    //
	    a += tmp_kernel(i, j);
	  }
    
	  if(a <= 0.) {
	    //
	    std::cerr << model.name() << " Well: cannot satisfy the constant collision rate at energy = "
		      << energy_bin(i) / Phys_const::incm
		      << " 1/cm\n";
	    
	    throw Error::Logic();
	  }
	  else {
	    //
	    a /= c;

	    for(int j = jmin; j < i; ++j) {
	      //
	      tmp_kernel(i, j) *= a;
	      
	      tmp_kernel(j, i) = tmp_kernel(i, j) * state_density(i) / state_density(j) * thermal_factor(i - j);
	    }
      
	    tmp_kernel(i, i) = tmp_kernel(i, i) * a + kernel_fraction(b);
	  }
	}// energy grid cycle
	//
      }// recursion up
      //
      // energy transfer DOWN probability functional form predefined
      //
      else if(Model::Kernel::flags() & Model::Kernel::DOWN) {
	//
	for(int i = 0; i < size(); ++i) {// energy grid cycle
	  //
	  itemp = i + energy_transfer.size();
	  
	  const int jmax = itemp < size() ? itemp : size();
	  
	  itemp = i - energy_transfer.size() + 1;
	  
	  const int jmin = itemp  > 0 ? itemp : 0;

	  // normalization constant
	  //
	  c = 0.;
	  
	  for(int j = i; j < jmax; ++j) {
	    //
	    dtemp = energy_transfer[j - i];
	    
	    if(Model::Kernel::flags() & Model::Kernel::DENSITY)
	      //
	      dtemp *= state_density(j);
	    
	    tmp_kernel(i, j) = -dtemp;
	    
	    c += dtemp;
	  }

	  // up-transitions contribution
	  //
	  a = kernel_fraction(b);
	  
	  for(int j = jmin; j < i; ++j)
	    //
	    a += tmp_kernel(i, j);
	  
    
	  if(a < 0.) {
	    //
	    IO::log << IO::log_offset << model.name()
	      //
		    << " Well: cannot satisfy the constant collision frequency at energy = "
	      //
		    << energy_bin(i) / Phys_const::incm
	      //
		    << " 1/cm: ";
	    
	    if(Model::Kernel::flags() & Model::Kernel::NOTRUN) {
	      //
	      dtemp = kernel_fraction(b) - a;
	      
	      IO::log << "renormalizing collision frequency: " << dtemp << "\n";
	      
	      tmp_kernel(i, i) = dtemp;
	      
	      for(int j = i + 1; j < jmax; ++j) {
		//
		tmp_kernel(i, j) = 0.;
		
		tmp_kernel(j, i) = 0.;
	      }
	    }
	    else {
	      //
	      IO::log << "truncating the well\n";
	      
	      _state_density.resize(i);
	      
	      break;
	    }
	  }
	  else {
	    //
	    a /= c;
	  
	    for(int j = i + 1; j < jmax; ++j) {
	      //
	      tmp_kernel(i, j) *= a;
	      
	      tmp_kernel(j, i) = tmp_kernel(i, j) * state_density(i) / state_density(j) / thermal_factor(j - i);
	    }
	    
	    tmp_kernel(i, i) = tmp_kernel(i, i) * a + kernel_fraction(b);
	  }
	}// energy grid cycle
	//
      }// recursion down
      //
      // default
      //
      else {
	//
	for(int e = 0; e < size(); ++e) {
	  //
	  dtemp = energy_transfer[0];
	  
	  for(int i = 1; i < energy_transfer.size(); ++i) {
	    //
	    itemp = e + i;

	    if(itemp == size())
	      //
	      break;
	    
	    tmp_kernel(e, itemp) = energy_transfer[i];

	    tmp_kernel(itemp, e) = energy_transfer[i] * _state_density[e] / _state_density[itemp]
	      //
	      / std::exp(i * energy_step() / temperature());

	    dtemp += tmp_kernel(e, itemp) + tmp_kernel(itemp, e);
	  }

	  dtemp /= kernel_fraction(b);
	  
	  for(int i = 1; i < energy_transfer.size(); ++i) {
	    //
	    itemp = e + i;

	    if(itemp == size())
	      //
	      break;
	    
	    tmp_kernel(e, itemp) /= -dtemp;

	    tmp_kernel(itemp, e) /= -dtemp;
	  }
	}

	for(int e = 0; e < size(); ++e)
	  //
	  for(int i = 1; i < energy_transfer.size(); ++i) {
	    //
	    itemp = e - i;

	    if(itemp >= 0.)
	      //
	      tmp_kernel(e, e) -= tmp_kernel(e, itemp);

	    itemp = e + i;

	    if(itemp < size())
	      //
	      tmp_kernel(e, e) -= tmp_kernel(e, itemp);
	  }
      }
      
      if(_kernel.size() != size()) {
	//
	break;
      }
      else
	//
	_kernel += tmp_kernel;
      //
    }//
    //
  } while(_kernel.size() != size());

#ifdef DEBUG
  
  IO::log << IO::log_offset << model.name() << " well: kernel diagonal elements:\n";
  
  for(int i = 0; i < size(); ++i)
    //
    IO::log << IO::log_offset << std::setw(5) << i  << std::setw(15) << _kernel(i, i)  << "\n";
  
#endif

}

void MasterEquation::Well::_set_crm_basis ()
{
  const char funame [] = "MasterEquation::Well::_set_crm_basis: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  /**************** SETTING RELAXATION MODE BASIS & PARTITION FUNCTION ****************/

  _boltzman.resize(size());
  
  _boltzman_sqrt.resize(size());
  
  _crm_basis.resize(size(), size() - 1);
  
  _crm_basis = 0.;
  
  _weight = 0.;
  
  std::vector<double> nfac(crm_size());
  
  for(int i = 0; i < size(); ++i) {
    //
    dtemp = state_density(i) * thermal_factor(i);
    
    _boltzman[i] = dtemp;
    
    _boltzman_sqrt[i] = std::sqrt(dtemp);
    
    if(i) {
      //
      itemp = i - 1;
      
      _crm_basis(i, itemp) = - _weight;
      
      nfac[itemp] = std::sqrt((_weight / dtemp + 1.) * _weight);
    }
    
    for(int r = i; r < crm_size(); ++r)
      //
      _crm_basis(i, r) = dtemp;
    
    _weight += dtemp;    
  }

  for(int r = 0; r < crm_size(); ++r) {
    //
    itemp = r + 2;
    
    for(int i = 0; i < itemp; ++i)
      //
      _crm_basis(i, r) /= nfac[r];
  }

  _weight_sqrt = std::sqrt(_weight);
}

MasterEquation::Well::Well (const Model::Well& model)
{
  const char funame [] = "MasterEquation::Well: ";

  int                 itemp;
  double              dtemp;
  bool                btemp;
  Lapack::Vector      vtemp;
  std::string         stemp;
  
  std::clock_t start_time;

  _collision_factor = 0.;
  
  _kernel_fraction.resize(Model::buffer_size());
  
  for(int b = 0; b < Model::buffer_size(); ++b) {
    //
    dtemp = (*model.collision(b))(temperature()) * Model::buffer_fraction(b);
    
    _kernel_fraction[b] = dtemp;
    
    _collision_factor += dtemp;
  }

  for(int b = 0; b < Model::buffer_size(); ++b)
    //
    _kernel_fraction[b] /= _collision_factor;

  _real_weight = model.weight(temperature()) *
    //
    std::exp((energy_reference() - model.ground()) / temperature());

  // state density
  //
  start_time = std::clock();

  _set_state_density(model);

  IO::log << IO::log_offset << model.name() << " Well: density of states done, elapsed time[sec] = "
    //
	    << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;

  // collisiona relaxation kernel
  //
  start_time = std::clock();

  _set_kernel(model);

  IO::log << IO::log_offset << model.name() << " Well: collisional energy transfer kernel done, elapsed time[sec] = "
    //
	    << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;

  // CRM basis
  //
  start_time = std::clock();

  _set_crm_basis();

  IO::log << IO::log_offset << model.name() << " Well: relaxation modes basis done, elapsed time[sec] = "
    //
	    << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;

  // collisional relaxation kernel in CRM basis
  //
  start_time = std::clock();

  _crm_bra = _crm_basis.copy();
  
  for(int i = 0; i < size(); ++i)
    //
    _crm_bra.row(i) /= boltzman(i);

  _crm_kernel = Lapack::SymmetricMatrix(_crm_basis.transpose() * _kernel * _crm_bra);

  IO::log << IO::log_offset << model.name()
    //
	  << " Well: kernel in relaxation modes basis done, elapsed time[sec] = "
    //
	  << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;  

  /*
  Lapack::SymmetricMatrix ktest(size());

  for(int i = 0; i < size(); ++i) {
    //
    ktest(i, i) = _kernel(i, i);
    
    for(int j = i + 1; j < size(); ++j)
      //
      ktest(i, j) = _kernel(i, j) * boltzman_sqrt(i) / boltzman_sqrt(j);
  }

  Lapack::Matrix kvec(size());

  Lapack::Vector kval = ktest.eigenvalues(&kvec);

  IO::log << IO::log_offset << model.name() << " Well: energy transfer kernel lowest eigenvalues:\n";

  itemp = size() < 10 ? size() : 10;
  
  for(int i = 0; i < itemp; ++i)
    //
    IO::log << IO::log_offset << std::setw(15) << kval[i] << "\n";

  stemp = "kvec_" + model.name() + "_" + IO::String(int(temperature() / Phys_const::kelv)) + ".dat";

  std::ofstream kvec_out(stemp.c_str());
  
  for(int e = 0; e < size(); ++e) {
    //
    for(int i = 0; i < itemp; ++i)
      //
      kvec_out << std::setw(15) << kvec(e, i) / boltzman_sqrt(e) * boltzman_sqrt(size() - 1);

    kvec_out << "\n";
  }
  */
  
  // minimal collisional relaxation eigenvalue
  //
  start_time = std::clock();

  vtemp = _crm_kernel.eigenvalues();
  
  _min_relax_eval = vtemp.front();
  
  _max_relax_eval = vtemp.back();

  IO::log << IO::log_offset << model.name() << " Well: relaxation eigenvalues done, elapsed time[sec] = "
    //
	  << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;  

  IO::log << IO::log_offset << model.name() << " Well: minimal relaxation eigenvalue = "
    //
	  << _min_relax_eval << "\n";
  
  IO::log << IO::log_offset << model.name() << " Well: maximal relaxation eigenvalue = "
    //
	  << _max_relax_eval << "\n";

  /*
    #ifdef DEBUG
    
    start_time = std::clock();
    
    // symmetrized kernel eigenvalues
    Lapack::Vector f_0(size());
    for(int i = 0; i < size(); ++i)
    f_0[i] = std::sqrt(state_density(i) * thermal_factor(i));
    Lapack::SymmetricMatrix sym_kernel(size());
    for(int i = 0; i < size(); ++i)
    for(int j = i; j < size(); ++j)
    sym_kernel(i, j) = kernel(i, j) * f_0[i] / f_0[j];
    vtemp = sym_kernel.eigenvalues();
    
    IO::log << IO::log_offset << model.name() << " Well: symmetrized kernel eigenvalues done, elapsed time[sec] = "
    << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;  
    IO::log << IO::log_offset << model.name() << " Well: symmetrized kernel eigenvalues: "
    << std::setw(13) << vtemp[0] << std::setw(13) << vtemp[1] << std::setw(13) << vtemp.back() << "\n";
    
    #endif
  

  // escape rate
  //
  if(model.escape_size()) {
    //
    IO::log << IO::log_offset << model.name() << " Well: Escape rate:\n";
    
    _escape_rate.resize(size());
    
    double ener = energy_reference();

    for(int i = 0; i < size(); ++i, ener -= energy_step()) {
      //
      dtemp = model.escape_rate(ener);

      IO::log << IO::log_offset 
	      << "    E[kcal/mol] = " << std::setw(13) << ener / Phys_const::kcal
	      << "    rate[1/sec] = " << std::setw(13) << dtemp / Phys_const::herz << "\n";
      
      _escape_rate[i] = dtemp;
    }
  }
  */
    
  // radiational transitions
  //
  if(model.oscillator_size()) {
    //
    _crm_radiation_rate.resize(crm_size());
    
    _radiation_rate.resize(size());
    
    _radiation_rate = 0.;
    
    _crm_radiation_rate = 0.;

    for(int ue = 0; ue < size() - 1; ++ue) {
      //
      double ener = energy_bin(ue);
    
      for(int f = 0; f < model.oscillator_size(); ++f) {
	//
	itemp = (int)round(model.oscillator_frequency(f) / energy_step());
	
	if(!itemp)
	  //
	  continue;

	double rad_prob = model.transition_probability(ener, temperature(), f);
	if(rad_prob <= 0.)
	  continue;

	int le = itemp + ue;
	le = le < size() ? le : size() - 1; 

	dtemp = boltzman_sqrt(ue) / boltzman_sqrt(le);
	
	_radiation_rate(ue, le) -= rad_prob * dtemp;
	
	_radiation_rate(ue, ue) += rad_prob;
	
	_radiation_rate(le, le) += rad_prob * dtemp * dtemp;

#pragma omp parallel for default(shared) schedule(dynamic)
	
	for(int r1 = 0; r1 < crm_size(); ++r1) {
	  //
	  for(int r2 = r1; r2 < crm_size(); ++r2)
	    //
	    _crm_radiation_rate(r1, r2) += rad_prob * boltzman(ue)
	      //
	      * (crm_bra(ue, r1) - crm_bra(le, r1))
	      //
	      * (crm_bra(ue, r2) - crm_bra(le, r2));
	}
      }
    }
  }
  
  IO::log << IO::log_offset << model.name()
    //
	  << " Well:       grid size = " << size() << "\n"
    //
	  << IO::log_offset << model.name()
    //
	  << " Well:      real depth = " << int(model.ground() / Phys_const::incm) << " 1/cm\n"
    //
	  << IO::log_offset << model.name()
    //
	  << " Well: effective depth = "
    //
	  << int(energy_bin(size()) / Phys_const::incm) << " 1/cm\n";
}

/********************************************************************************************
 ************************************ SETTING BARRIER ***************************************
 ********************************************************************************************/

MasterEquation::Barrier::Barrier (const Model::Species& model)
{
  const char funame [] = "MasterEquation::Barrier::Barrier: ";

  int    itemp;
  double dtemp;

  int new_size = (int)std::ceil((energy_reference() - model.ground()) / energy_step());

  if(is_global_cutoff) {
    //
    itemp =(int)std::ceil((energy_reference() - global_cutoff) / energy_step());
    
    new_size = itemp < new_size ? itemp : new_size;      
  }

  _state_number.resize(new_size);
  
  resize_thermal_factor(new_size);

  double ener = energy_reference();
  
  _weight = 0.;
  
  for(int e = 0; e < size(); ++e, ener -= energy_step()) {
    //
    dtemp = model.states(ener);
    
    if(dtemp <= 0.) {
      //
      IO::log << IO::log_offset  << model.name() << " Barrier: nonpositive number of states at "
	//
	      << ener / Phys_const::incm  << " 1/cm => truncating\n";
      
      _state_number.resize(e);
      
      break;
    }
    
    _state_number[e] = dtemp;
    
    _weight += dtemp * thermal_factor(e);
  }
  
  _weight *= energy_step() / temperature();

  _real_weight = model.weight(temperature()) * std::exp((energy_reference() - model.ground()) / temperature());  

  IO::log << IO::log_offset << model.name()
    //
	  << " Barrier:        grid size = " << size() << "\n"
    //
	  << IO::log_offset << model.name()
    //
	  << " Barrier:      real height = " << int(model.ground() / Phys_const::incm) << " 1/cm\n"
    //
	  << IO::log_offset << model.name()
    //
	  << " Barrier: effective height = "
    //
	  << int(energy_bin(size()) / Phys_const::incm) << " 1/cm\n";
}

void MasterEquation::Barrier::truncate (int new_size) 
{
  const char funame [] = "MasterEquation::Barrier::truncate: ";

  double dtemp;
  
  if(new_size == size())
    //
    return;

  if(new_size <= 0 || new_size > size()) {
    //
    std::cerr << funame << "out of range: " << new_size << "\n";
    
    throw Error::Range();
  }

  dtemp = 0.;
  
  for(int i = new_size; i < size(); ++i)
    //
    dtemp += _state_number[i] * thermal_factor(i);
  
  dtemp *= energy_step() / temperature();

  _state_number.resize(new_size);
  
  _weight -= dtemp;
  
  IO::log << IO::log_offset << funame << "grid size = " << size() << "\n"
    //
	  << IO::log_offset << funame << "real weight = " << _real_weight << "\n"
    //
	  << IO::log_offset << funame << "effective weight = " << _weight * energy_step() << "\n";
}

/********************************************************************************************
 ********************************** SETTING BIMOLECULAR *************************************
 ********************************************************************************************/

MasterEquation::Bimolecular::Bimolecular (const Model::Bimolecular& model)
{ 
  _weight  = model.weight(temperature()) * std::exp((energy_reference() - model.ground()) / temperature());
}

/********************************************************************************************
 ********************************** CALCULATION METHODS *************************************
 ********************************************************************************************/

/********************************************************************************************
 ******************************* LOW CHEMICAL EIGENVALUE METHOD *****************************
 ********************************************************************************************/

// the method based on the assumption that the chemical eigenvalues are small in comparison
// to the relaxation eigenvalues

double triple_product ( const double* p1, const double* p2, const double* p3, int n)
{
  double res = 0.;

#pragma omp parallel for default(shared) reduction(+: res) schedule(static)
	
  for(int i = 0; i < n; ++i) {
    //
    res += p1[i] * p2[i] * p3[i];
  }

  return res;
}

void MasterEquation::low_eigenvalue_matrix (Lapack::SymmetricMatrix& k_11, Lapack::SymmetricMatrix& k_33, 
					    Lapack::Matrix& k_13, Lapack::Matrix& l_21) 
{
  const char funame [] = "MasterEquation::low_eigenvalue_matrix: ";

  IO::Marker funame_marker(funame);

  if(!isset()) {
    std::cerr << funame << "reactive complex is not set\n";
    throw Error::Init();
  }

  int            itemp;
  double         dtemp;
  bool           btemp;
  std::string    stemp;

  Lapack::Vector vtemp;
  Lapack::Matrix mtemp;  

  // total collisional relaxation modes dimension and index shifts for individual wells
  //
  std::vector<int> well_shift(Model::well_size());
  
  itemp = 0;
  
  for(int w = 0; w < Model::well_size(); itemp += well(w++).crm_size())
    //
    well_shift[w] = itemp;
  
  const int crm_size = itemp;
  
  //IO::log << IO::log_offset << "global relaxation matrix size =  " << crm_size << "\n";

  // kinetic relaxation matrices

  // chemical modes
  //
  k_11.resize(Model::well_size());
  k_11 = 0.;

  // chemical basis
  //
  k_13.resize(Model::well_size(),  Model::bimolecular_size());
  k_13 = 0.;

  // chemical-to-collision modes
  //
  Lapack::Matrix k_21(crm_size, Model::well_size());
  k_21 = 0.;

  // collision modes
  //
  Lapack::SymmetricMatrix k_22(crm_size);
  k_22 = 0.;

  // relaxational basis
  //
  Lapack::Matrix k_23(crm_size, Model::bimolecular_size());
  k_23 = 0.;
  
  /*********************************** GLOBAL MATRICES **************************************/

  {
    IO::Marker work_marker("setting up kinetic matrices", IO::Marker::ONE_LINE);

    // k_11 initialization
    //
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      const int& w1 = Model::inner_connect(b).first;
      
      const int& w2 = Model::inner_connect(b).second;
      
      dtemp = 0.;
      
      for(int i = 0; i < inner_barrier(b).size(); ++i)
	//
	dtemp += inner_barrier(b).state_number(i) * thermal_factor(i);
      
      k_11(w1, w2) = - dtemp / 2. / M_PI / well(w1).weight_sqrt() / well(w2).weight_sqrt();
    }
  
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = 0.;
      
      for(int i = 0; i < cum_stat_num[w].size(); ++i)
	//
	dtemp += cum_stat_num[w][i] * thermal_factor(i);
      
      dtemp /= 2. * M_PI * well(w).weight();
      
      k_11(w, w) = dtemp;
      
      if(Model::well(w).escape_size()) {
	//
	dtemp = 0.;
	
	for(int i = 0; i < well(w).size(); ++i)
	  //
	  dtemp += Model::well(w).escape_rate(energy_bin(i)) * well(w).boltzman(i);
	
	dtemp /= well(w).weight();
	
	k_11(w, w) += dtemp;
      }
    }

    // k_21 initialization
    //
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      int w1 = Model::inner_connect(b).first;
      
      int w2 = Model::inner_connect(b).second;
      
      vtemp.resize(inner_barrier(b).size());
      
      for(int dir = 0; dir < 2; ++dir) {
	//
	if(dir)
	  //
	  std::swap(w1, w2);

	for(int i = 0; i < inner_barrier(b).size(); ++i)
	  //
	  vtemp[i] = inner_barrier(b).state_number(i) / 2. / M_PI / well(w1).state_density(i)
	    //
	    / well(w2).weight_sqrt();

	for(int r = 0; r < well(w1).crm_size(); ++r)
	  //
	  k_21(r + well_shift[w1], w2) = -parallel_vdot(well(w1).crm_column(r), vtemp, vtemp.size()) ;
      }
    }
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      if(Model::well(w).escape_size()) {
	//
	vtemp.resize(well(w).size());
	
	for(int i = 0; i < well(w).size(); ++i)
	  //
	  vtemp[i] = Model::well(w).escape_rate(energy_bin(i)) / well(w).weight_sqrt();
      }
      else {
	//
	vtemp.resize(cum_stat_num[w].size());
	
	vtemp = 0.;
      }
      
      for(int i = 0; i < cum_stat_num[w].size(); ++i)
	//
	vtemp[i] += cum_stat_num[w][i] / 2. / M_PI / well(w).state_density(i) / well(w).weight_sqrt();
      
      for(int r = 0; r < well(w).crm_size(); ++r)
	//
	k_21(r + well_shift[w], w) =  parallel_vdot(well(w).crm_column(r), vtemp, vtemp.size());
    }

    // k_22 initialization
    // nondiagonal isomerization
    //
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      int w1 = Model::inner_connect(b).first;
      
      int w2 = Model::inner_connect(b).second;    

      vtemp.resize(inner_barrier(b).size());
      
      for(int i = 0; i < inner_barrier(b).size(); ++i)
	//
	vtemp[i] = inner_barrier(b).state_number(i) / 2. / M_PI / well(w1).state_density(i)
	  //
	  / well(w2).state_density(i) / thermal_factor(i);

      for(int r1 = 0; r1 < well(w1).crm_size(); ++r1)
	//
	for(int r2 = 0; r2 < well(w2).crm_size(); ++r2) {
	  //
	  k_22(r1 + well_shift[w1], r2 + well_shift[w2]) =
	    //
	    -triple_product(well(w1).crm_column(r1), well(w2).crm_column(r2), vtemp, vtemp.size());
	}
    }

    // diagonal isomerization
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      if(Model::well(w).escape_size()) {
	//
	vtemp.resize(well(w).size());
	
	for(int i = 0; i < well(w).size(); ++i)
	  //
	  vtemp[i] = Model::well(w).escape_rate(energy_bin(i)) / well(w).boltzman(i);
      }
      else {
	//
	vtemp.resize(cum_stat_num[w].size());
	
	vtemp = 0.;
      }
      for(int i = 0; i < cum_stat_num[w].size(); ++i) {
	//
	dtemp = well(w).state_density(i);
	
	vtemp[i] += cum_stat_num[w][i] / 2. / M_PI / dtemp / dtemp / thermal_factor(i);
      }

      for(int r1 = 0; r1 < well(w).crm_size(); ++r1)
	//
	for(int r2 = r1; r2 < well(w).crm_size(); ++r2) {
	  //
	  k_22(r1 + well_shift[w], r2 + well_shift[w]) =
	    //
	    triple_product(well(w).crm_column(r1), well(w).crm_column(r2), vtemp, vtemp.size());
	}
    }

    // collisional energy transfer
    //
    for(int w = 0; w < Model::well_size(); ++w) {

#pragma omp parallel for default(shared) schedule(dynamic)

      for(int r1 = 0; r1 < well(w).crm_size(); ++r1) {
	//
	for(int r2 = r1; r2 < well(w).crm_size(); ++r2)
	  //
	  k_22(r1 + well_shift[w], r2 + well_shift[w]) +=  well(w).collision_frequency() *
	    //
	    well(w).crm_kernel(r1, r2);
      }
    }

    // radiational transitions contribution
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      if(well(w).radiation()) {

#pragma omp parallel for default(shared) schedule(dynamic)
	
	for(int r1 = 0; r1 < well(w).crm_size(); ++r1) {
	  //
	  for(int r2 = r1; r2 < well(w).crm_size(); ++r2)
	    //
	    k_22(r1 + well_shift[w], r2 + well_shift[w]) +=
	      //
	      well(w).crm_radiation_rate(r1, r2);
	}
      }

    // bimolecular number of states
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int w = Model::outer_connect(b).first;
      
      const int p = Model::outer_connect(b).second;

      // chemical modes
      //
      dtemp = 0.;
      
      for(int i = 0; i < outer_barrier(b).size(); ++i)
	//
	dtemp += outer_barrier(b).state_number(i) * thermal_factor(i);
      
      k_13(w, p) = dtemp / 2. / M_PI / well(w).weight_sqrt();
    
      // relaxation modes
      //
      vtemp.resize(outer_barrier(b).size());
      
      for(int i = 0; i < outer_barrier(b).size(); ++i)
	//
	vtemp[i] = outer_barrier(b).state_number(i) / 2. / M_PI / well(w).state_density(i);
    
      for(int r = 0; r < well(w).crm_size(); ++r)
	//
	k_23(r + well_shift[w], p) = parallel_vdot(well(w).crm_column(r), vtemp, vtemp.size());
    }
  }

  /****************************** SOLVING MASTER EQUATION *************************************/

  {
    IO::Marker work_marker("inverting kinetic matrices", IO::Marker::ONE_LINE);
    
    Lapack::Cholesky l_22(k_22);
    
    l_21 = l_22.invert(k_21);

    // well-to-well rate coefficients
    //
    k_11 -= Lapack::SymmetricMatrix(k_21.transpose() * l_21);

    // bimolecular-to-bimolecular rate coefficients
    //
    if(Model::bimolecular_size()) {
      //
      k_33 = Lapack::SymmetricMatrix(k_23.transpose() * l_22.invert(k_23)); 

      // well-to-bimolecular rate coefficients
      //
      k_13 -= l_21.transpose() * k_23;
    }
  }
}

void MasterEquation::low_eigenvalue_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags) 
  
{
  const char funame [] = "MasterEquation::low_eigenvalue_method: ";

  IO::Marker funame_marker(funame);

  int            itemp;
  double         dtemp;
  bool           btemp;
  std::string    stemp;

  Lapack::Vector vtemp;
  Lapack::Matrix mtemp;  

  IO::log << IO::log_offset << "Pressure = ";
  
  switch(pressure_unit) {
    //
  case BAR:
    //
    IO::log << pressure() / Phys_const::bar << " bar";
    
    break;
    //
  case TORR:
    //
    IO::log << pressure() / Phys_const::tor << " torr";
    
    break;
    //
  case ATM:
    //
    IO::log << pressure() / Phys_const::atm << " atm";
    
    break;
  }

  IO::log << "\t Temperature = "
    //
	  << temperature() / Phys_const::kelv << " K\n";

  rate_data.clear();

  Lapack::SymmetricMatrix k_11;
  
  Lapack::Matrix          k_13;
  
  Lapack::SymmetricMatrix k_33;
  
  Lapack::Matrix          l_21;
  
  low_eigenvalue_matrix(k_11, k_33, k_13, l_21);

  std::vector<int> well_shift(Model::well_size());
  
  itemp = 0;
  
  for(int w = 0; w < Model::well_size(); itemp += well(w++).crm_size())
    //
    well_shift[w] = itemp;
  
  const int crm_size = itemp;

  // hot distribution
  //
  Lapack::Matrix hot_chem;
  
  if(hot_index.size()) {
    //
    hot_chem.resize(hot_index.size(), Model::well_size());
    
    hot_chem = 0.;
    
    vtemp.resize(crm_size);
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
      //
      hot_chem(count, hit->first) = 1. / well(hit->first).weight_sqrt();
	
      vtemp = 0.;
      
      itemp = 0;
      
      for(int r = 0; r < well(hit->first).crm_size(); ++r)
	//
	vtemp[r + well_shift[hit->first]] = well(hit->first).crm_bra(hit->second, r);
      
      hot_chem.row(count) -= vtemp * l_21;
      }  
  }
  
  // well escape
  //
  Lapack::Matrix escape_chem;
  
  if(Model::escape_size()) {
    //
    escape_chem.resize(Model::escape_size(), Model::well_size());
    
    escape_chem = 0.;
    
    vtemp.resize(crm_size);
    
    for(int count = 0; count < Model::escape_size(); ++count) {
      //
      const int ew = Model::escape_channel(count).first;

      const int ei = Model::escape_channel(count).second;

      dtemp = 0.;
      
      for(int i = 0; i < well(ew).size(); ++i)
	//
	dtemp += Model::well(ew).escape_rate(energy_bin(i), ei)  * well(ew).boltzman(i);
      
      dtemp /= well(ew).weight_sqrt();

      escape_chem(count, ew) = dtemp;

      vtemp = 0.;
      
      for(int r = 0; r < well(ew).crm_size(); ++r)
	//
	for(int i = 0; i < well(ew).size(); ++i)
	  //
	  vtemp[r + well_shift[ew]] += well(ew).crm_basis(i, r) * Model::well(ew).escape_rate(energy_bin(i), ei);
 
      escape_chem.row(count) -= vtemp * l_21;
    }
    //
  }// well escape

  // minimal energy relaxation eigenvalue
  //
  double min_relax_eval, max_relax_eval;
  
  int wmin, wmax;
  
  for(int w = 0; w < Model::well_size(); ++w) {
    //
    double rmin, rmax;

    // radiational transitions contribution
    //
    if(well(w).radiation()) {
      //
      Lapack::SymmetricMatrix erk(well(w).crm_size());
      
#pragma omp parallel for default(shared) schedule(dynamic)
	
      for(int r1 = 0; r1 < well(w).crm_size(); ++r1) {
	//
	  for(int r2 = r1; r2 < well(w).crm_size(); ++r2)
	    //
	    erk(r1, r2) =  well(w).crm_radiation_rate(r1, r2)
	      //
	      + well(w).collision_frequency() * well(w).crm_kernel(r1, r2);
	}
      
      vtemp = erk.eigenvalues();
      //
      rmin = vtemp.front();
      //
      rmax = vtemp.back();
    }
    else {
      //
      rmin = well(w).minimal_relaxation_eigenvalue() * well(w).collision_frequency();

      rmax = well(w).maximal_relaxation_eigenvalue() * well(w).collision_frequency();
    }

    if(!w || rmin < min_relax_eval) {
      //
      min_relax_eval = rmin;
      
      wmin = w;
    }
    if(!w || rmax > max_relax_eval) {
      //
      wmax = w;
      
      max_relax_eval = rmax;
    }
  }

  IO::log << IO::log_offset
    //
	  << "minimal relaxation eigenvalue / collisional frequency = "
    //
	  << min_relax_eval / well(wmin).collision_frequency() << "\n";

  // chemical eigenvalues and eigenvectors
  //
  Lapack::Matrix chem_evec(Model::well_size());

#ifdef WITH_MPACK
  
  Lapack::Vector chem_eval = Mpack::dd_eigenvalues(k_11, &chem_evec);

#else
  
  Lapack::Vector chem_eval = k_11.eigenvalues(&chem_evec);

#endif
  
  // relaxational projection of the chemical eigenvector
  //
  l_21 = l_21 * chem_evec;
  
  Lapack::Vector rel_proj(Model::well_size());
  
  for(int l = 0; l < Model::well_size(); ++l)
    //
    rel_proj[l] = vdot(l_21.column(l));

  // eigenvalues output
  //
  if(eval_out.is_open()) {
    //
    eval_out << std::setw(13) << temperature() / Phys_const::kelv;
    
    switch(pressure_unit) {
      //
    case BAR:
      //
      eval_out << pressure() / Phys_const::bar;
      
      break;
      //
    case TORR:
      //
      eval_out << pressure() / Phys_const::tor;
      
      break;
      //
    case ATM:
      //
      eval_out << pressure() / Phys_const::atm;
      
      break;
    }
    
    eval_out << std::setw(13) << well(wmin).collision_frequency() / Phys_const::herz
      //
	     << std::setw(13) << min_relax_eval / well(wmin).collision_frequency();
    
    for(int l = 0; l < Model::well_size(); ++l)
      //
      eval_out << std::setw(13) << chem_eval[l] / well(wmin).collision_frequency() << std::setw(13) << rel_proj[l];
    
    eval_out << "\n";
  }

  /************************************* EIGENVECTOR OUTPUT *****************************************/

  IO::log << std::setprecision(3)
    //
	  << IO::log_offset << "eigenvector populations normalized:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(10) << "*R"
    //
	  << std::setw(10) << "*P";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  IO::log << "\n";

  // maximal population
  //
  for(int l = 0; l < Model::well_size(); ++l) {
    //
    double pos_pop = 0.;
    
    double neg_pop = 0.;
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = chem_evec(w, l) * well(w).weight_sqrt();
      
      if(dtemp > 0.)
	//
	pos_pop += dtemp;
      
      if(dtemp < 0.)
	//
	neg_pop += dtemp;
    }
    
    double max_pop = pos_pop > -neg_pop ? pos_pop : neg_pop;
    
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(10) << chem_eval[l] / min_relax_eval
      //
	    << std::setw(10) << rel_proj[l];
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = chem_evec(w, l) * well(w).weight_sqrt() / max_pop;
      //
      IO::log << std::setw(10);
      
      if(dtemp > .01 || dtemp < -.01) {
	//
	IO::log << dtemp;
      }
      else
	//
	IO::log << 0;
    }
    
    IO::log << "\n";
  }

  IO::log << IO::log_offset << std::setw(5) << "*R"
    //
	  << " - eigenvalue over the relaxation limit\n"
    //
	  << IO::log_offset << std::setw(5) << "*P"
    //
	  << " - eigenvector projection squared on the relaxation subspace\n"
    //
	  << IO::log_offset << "eigenvector projections:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(10) << "*Q"
    //
	  << std::setw(10) << "*P";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  IO::log << "\n";

  for(int l = 0; l < Model::well_size(); ++l) {
    //
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(10) << chem_eval[l] / well(wmin).collision_frequency()
      //
	    << std::setw(10) << rel_proj[l];
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(10) << chem_evec(w, l);
    
    IO::log << "\n";
  }
  
  IO::log << IO::log_offset
    //
	  << std::setw(5) << "*Z"
    //
	  << std::setw(10) << "---"
    //
	  << std::setw(10) << "---";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << well(w).weight_sqrt();
  
  IO::log << "\n";

  IO::log << IO::log_offset << std::setw(5) << "*Q"
    //
	  << " - eigenvalue over the collision frequency\n"
    //
	  << IO::log_offset << std::setw(5) << "*P"
    //
	  << " - eigenvector projection squared on the relaxation subspace\n"
    //
	  << IO::log_offset << std::setw(5) << "*Z"
    //
	  << " - well partition function square root\n"
    //
	  << std::setprecision(6);
  
  // rate coefficients connecting eigenstates to bimolecular products
  //
  Lapack::Matrix chem_bim = chem_evec.transpose() * k_13;

#ifdef DEBUG

  IO::log << IO::log_offset << "eigenvector-to-bimolecular projections:\n";

  IO::log << IO::log_offset << std::setw(5) << "L\\P";
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    IO::log << std::setw(13) << Model::bimolecular(p).name();
  
  IO::log << "\n";
  
  for(int l = 0; l < Model::well_size(); ++l) {
    //
    IO::log << IO::log_offset << std::setw(5) << l;
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(13) << chem_bim(l, p);
    
    IO::log << "\n";
  }

#endif

  // bimolecular rate units
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  // bimolecular-to-bimolecular rate coefficients
  //
  for(int i = 0; i < Model::bimolecular_size(); ++i)
    //
    if(bimolecular(i).weight() > 0.)
      //
      for(int j = 0; j < Model::bimolecular_size(); ++j) {
	//
	dtemp = k_33(i, j) * energy_step();
	
	// bimolecular reactant loss
	//
	if(i == j) {
	  //
	  dtemp = -dtemp;
	  
	  for(int b = 0; b < Model::outer_barrier_size(); ++b)
	    //
	    if(Model::outer_connect(b).second == i)
	      //
	      dtemp += temperature() / 2. / M_PI * outer_barrier(b).weight();	  
	}
	
	dtemp /= bimolecular(i).weight() * bru;
	
	rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 	
      }
  
  /********************************** CHEMICAL SUBSPACE DIMENSION ****************************************/

  if(default_partition.size())
    //
    itemp = default_partition.size();
  
  else if(chemical_threshold > 0.) {
    //
    dtemp = 1. / chemical_threshold;
    
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(chem_eval[itemp] / min_relax_eval  > dtemp)
	//
	break;
  }
  else if(chemical_threshold < 0. ) {
    //
    dtemp = -1. / chemical_threshold;
    
    for(itemp = Model::well_size(); itemp > 0; --itemp)
      //
      if(chem_eval[itemp - 1] / min_relax_eval < dtemp && chem_eval[itemp - 1] / chem_eval[itemp] < dtemp)
	//
	break;
  }
  else
    //
    itemp = Model::well_size();
	
  const int chem_size = itemp;

  if(chem_size != Model::well_size())
    //
    IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n";

  if(!chem_size)
    //
    return;

  // reduction of species
  //
  if(chem_size < Model::well_size()) {
    //
    // chemical eigenvectors
    //
    Lapack::Matrix pop_chem(Model::well_size(), chem_size);
    //
    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = chem_evec.column(l);

    // partitioning wells
    //
    if(default_partition.size()) {
      //
      // default reduction scheme
      //
      well_partition = default_partition;

      IO::log << IO::log_offset << "using default reduction scheme, projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    }
    else {
      //
      Group bimolecular_group;
      //
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group, std::vector<double>());
    }

    // convert chemical eigenvectors in the new basis
    //
    pop_chem = well_partition.basis().transpose() * pop_chem;

    std::vector<int>    group_index = well_partition.group_index();
    
    std::vector<double>      weight = well_partition.weight();
    
    std::vector<double> real_weight = well_partition.real_weight();
    
    Lapack::Matrix            basis = well_partition.basis();

    Lapack::Matrix m_direct(pop_chem);
    
    Lapack::Matrix m_inverse = m_direct.invert();

#ifdef DEBUG
    
    Lapack::Matrix one(m_direct * m_inverse);
    //
    one.diagonal() -= 1.;
    
    double val_max = -1.;
    
    for(int i = 0; i < one.size1(); ++i)
      //
      for(int j = 0; j < one.size2(); ++j) {
	//
	dtemp = one(i, j);
	
	dtemp = dtemp < 0. ? -dtemp : dtemp;
	
	if(dtemp > epsilon && dtemp > val_max)
	  
	  val_max = dtemp;
      }
    
    if(val_max > 0.)
      //
      IO::log << IO::log_offset << funame
	//
	      << "WARNING: matrix inversion error = " << val_max
	//
	      << " exceeds numerical accuracy = " << epsilon
	//
	      << "\n";

#endif

    IO::log << IO::log_offset << "species:\n"
      //
	    << IO::log_offset << std::setw(2) << "#"  << std::setw(15)
      //
	    << "assigned name" << IO::first_offset
      //
	    << "group\n";
    
    for(Pit g = well_partition.begin(); g != well_partition.end(); ++g) {
      //
      itemp = g - well_partition.begin();
      
      IO::log << IO::log_offset << std::setw(2) << itemp << std::setw(15)
	//
	      << Model::well(group_index[itemp]).name() << IO::first_offset;
      
      for(Git w = g->begin(); w != g->end(); ++w) {
	//
	if(w != g->begin())
	  //
	  IO::log << "+";
	
	IO::log << Model::well(*w).name();
      }
      
      IO::log << "\n";
    }
  
    // rate constant reduction
    //
    switch(reduction_method) {
      //
    case PROJECTION:
      //
      // isomerization
      //
      k_11 = Lapack::SymmetricMatrix(basis.transpose() * k_11 * basis);
      
      // dissociation
      //
      k_13 = basis.transpose() * k_13;
      
      // new chemcal eigenvalues
      //
      IO::log << IO::log_offset << "new chemical eigenvalues:";
      
      vtemp = k_11.eigenvalues();
      
      for(int l = 0; l < chem_size; ++l)
	//
	IO::log << std::setw(13) << vtemp[l] / well(wmin).collision_frequency();
      
      IO::log << "\n";
      
      // well-to-well rate coefficients
      //
      for(int i = 0; i < chem_size; ++i)
	//
	for(int j = 0; j < chem_size; ++j) {
	  //
	  dtemp = k_11(i, j) * std::sqrt(weight[i] * weight[j])
	    //
	    * energy_step() / real_weight[i] / Phys_const::herz;
	  
	  if(i != j)
	    //
	    dtemp = -dtemp;
	  
	  rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp;
	}
      
      // well-to-bimolecular rate coefficients
      //
      for(int w = 0; w < chem_size; ++w)
	//
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = k_13(w, p) * std::sqrt(weight[w])
	    //
	    * energy_step() / real_weight[w] / Phys_const::herz;
      
      // bimolecular-to-well rate coefficients
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(bimolecular(p).weight() > 0.)
	  //
	  for(int w = 0; w < chem_size; ++w)
	    //
	    rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = k_13(w, p) * std::sqrt(weight[w])
	      //
	      * energy_step() / bimolecular(p).weight() / bru;
      
      break;
      //
    case DIAGONALIZATION:
      //
      // well-to-well rate coefficients
      //
      for(int i = 0; i < chem_size; ++i)
	//
	for(int j = 0; j < chem_size; ++j) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_direct(j, l) * m_inverse(l, i) * chem_eval[l];
	  
	  dtemp *= std::sqrt(weight[i] * weight[j]) * energy_step() / real_weight[i] / Phys_const::herz;
	  
	  if(i != j)
	    //
	    dtemp = -dtemp;
	  
	  rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp;
	}
      
      // well-to-bimolecular rate coefficients
      //
      for(int w = 0; w < chem_size; ++w)
	//
	for(int p = 0; p < Model::bimolecular_size(); ++p) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_inverse(l, w) * chem_bim(l, p);
	  
	  dtemp *= std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
	}
      // bimolecular-to-well rate coefficients
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(bimolecular(p).weight() > 0.)
	  //
	  for(int w = 0; w < chem_size; ++w) {
	    //
	    dtemp = 0.;
	    
	    for(int l = 0; l < chem_size; ++l)
	      //
	      dtemp += m_direct(w, l) * chem_bim(l, p);
	    
	    dtemp *= std::sqrt(weight[w]) * energy_step() / bimolecular(p).weight() / bru;
	    
	    rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	  }
      
      break;
      //
    default:
      //
      std::cerr << funame << "unknown reduction method\n";
      //
      throw Error::Logic();
    }//
    //
  }// reduction of species
  //
  // no reduction (chem_size == well_size)
  //
  else {
    //
    // well-to-well rate coefficients
    //
    for(int i = 0; i < Model::well_size(); ++i)
      //
      for(int j = 0; j < Model::well_size(); ++j) {
	//
	dtemp = k_11(i, j) * well(i).weight_sqrt() * well(j).weight_sqrt()
	  //
	  * energy_step() / well(i).real_weight() / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(i, j)] = dtemp;
      }
    
    // well-to-bimolecular rate coefficients
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	rate_data[std::make_pair(w, Model::well_size() + p)] = k_13(w, p) * well(w).weight_sqrt()
	  //
	  * energy_step() / well(w).real_weight() / Phys_const::herz;
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(bimolecular(p).weight() > 0.)
	//
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  rate_data[std::make_pair(Model::well_size() + p, w)] = k_13(w, p) * well(w).weight_sqrt()
	    //
	    * energy_step() / bimolecular(p).weight() / bru;
    //
  }// no reduction
  //
}// Low chemical eigenvalue method

/************************************************************************************************
 ********************************** DIVIDE AND CONQUER METHOD ***********************************
 ************************************************************************************************/

void MasterEquation::divide_and_conquer_method (const Model::ChemGraph& g)
{
  const char funame [] = "MasterEquation::devide_and_conquer_method: ";

  std::list<Model::ChemGraph> cl = g.factorize();

  if(cl.size() != 1) {
    //
    std::cerr << funame << "original chemical graph is not connected\n";

    throw Error::Logic();
  }
  
  std::multimap<double, std::pair<int, int> > energy_map;

  for(std::set<int>::const_iterator bit = g.inner_set.begin(); bit != g.inner_set.end(); ++bit) {
    //
    energy_map.insert(std::make_pair(Model::inner_barrier(*bit).real_ground(), std::make_pair(0, *bit)));
  }

  for(std::set<int>::const_iterator bit = g.outer_set.begin(); bit != g.outer_set.end(); ++bit) {
    //
    energy_map.insert(std::make_pair(Model::outer_barrier(*bit).real_ground(), std::make_pair(1, *bit)));
  }

}

/********************************************************************************************
 ********************************** WELL REDUCTION METHOD ***********************************
 ********************************************************************************************/

// fast horizontal relaxation eigenvectors/eigenvalues are removed from ME
//
// and used only for bimolecular-to-bimolecular contributions

struct KineticBasis {
  //
  int active_size;

#if defined WITH_MPACK && defined WITH_DD

  Mpack_dd::Vector eigenvalue;
  
  Mpack_dd::Matrix eigenvector;

#else
  
  Lapack::Vector eigenvalue;
  
  Lapack::Matrix eigenvector;

#endif
  
  std::map<int, int> well_index_map;
  
  std::vector<int>   index_well_map;

  int size () const { return eigenvalue.size(); }
};

void MasterEquation::well_reduction_method (std::map<std::pair<int, int>, double>& rate_data,
					    //
					    Partition& well_partition, int flags)
{
  const char funame [] = "MasterEquation::well_reduction_method: ";

  static const double machine_eps = 1.e-14;

  if(!isset()) {
    //
    std::cerr << funame << "reactive complex is not set\n";
    
    throw Error::Init();
  }

  rate_data.clear();

  IO::Marker funame_marker(funame);
  
  IO::log << IO::log_offset << "Pressure = ";
  
  switch(pressure_unit) {
    //
  case BAR:
    //
    IO::log << pressure() / Phys_const::bar << " bar";
    
    break;
    //
  case TORR:
    //
    IO::log << pressure() / Phys_const::tor << " torr";
    
    break;
    //
  case ATM:
    //
    IO::log << pressure() / Phys_const::atm << " atm";
    
    break;
  }
  
  IO::log << "\t  Temperature = " << temperature() / Phys_const::kelv << " K\n"
    //
	  << "\t  Collision Frequency = "
    //
	  << well(0).collision_frequency() / Phys_const::herz << " 1/sec\n";

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  Lapack::Vector      vtemp;
  Lapack::Matrix      mtemp;

  dd_t  dd_temp;
  
  for(int w = 0; w < Model::well_size(); ++w) {
    //
    if(!w || well(w).size() > itemp)
      //
      itemp = well(w).size();
  }
  
  const int ener_index_max = itemp;

  /*******************************************************************
   ********************* PARTITIONING THE WELLS **********************
   *******************************************************************/

  std::vector<KineticBasis> kinetic_basis(ener_index_max);

  {// kinetically active basis

    IO::Marker well_part_marker("kinetically active basis");

    for(int e = 0; e < ener_index_max; ++e) {// energy cycle
      //
      // available energy bins
      //
      std::vector<int> well_array;
     
      std::map<int, int> well_index;
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	if(e < well(w).size()) {
	  //
	  well_index[w] = well_array.size();
	  
	  well_array.push_back(w);
	}

      kinetic_basis[e].well_index_map = well_index;
      
      kinetic_basis[e].index_well_map = well_array;

      // microcanonical kinetic matrix
      //

      std::vector<dd_t> density_sqrt(well_array.size());

      for(int i = 0; i < well_array.size(); ++i)
	//
	density_sqrt[i] = (dd_t)std::sqrt(well(well_array(i)).state_density(e));

#ifdef WITH_MPACK
      
      Mpack_dd::SymmetricMatrix km(well_array.size());

#else

      Lapack::SymmetricMatrix km(well_array.size());
      
#endif
      
      km = 0.;

      // nondiagonal isomerization contribution
      //
      for(int b = 0; b < Model::inner_barrier_size(); ++b) {
	//
	if(e < inner_barrier(b).size()) {
	  //
	  const int& w1 = Model::inner_connect(b).first;
	  
	  const int& w2 = Model::inner_connect(b).second;
      
	  std::map<int, int>::const_iterator i1 = well_index.find(w1);
	  
	  std::map<int, int>::const_iterator i2 = well_index.find(w2);

	  dtemp = inner_barrier(b).state_number(e) / 2. / M_PI;
	   
	  dd_temp = (dd_t)dtemp;
	  
	  km(i1->second, i2->second) -= dd_temp / density_sqrt[i1->second] / density_sqrt[i2->second];

	  km(i1->second, i1->second) += dd_temp / density_sqrt[i1->second] / density_sqrt[i1->second];
	  
	  km(i2->second, i2->second) += dd_temp / density_sqrt[i2->second] / density_sqrt[i2->second];
	}
      }

      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(e < outer_barrier(b).size()) {
	  //
	  const int& w = Model::inner_connect(b).first;
	  
	  std::map<int, int>::const_iterator i = well_index.find(w);

	  dtemp = outer_barrier(b).state_number(e) / 2. / M_PI;
	  
	  dd_temp = (dd_t)dtemp;
	  
	  km(i->second, i->second) += dd_temp / density_sqrt[i->second] / density_sqrt[i->second];
	}

      
      for(int w = 0; w < Model::well_size(); ++w) {
	//
	std::map<int, int>::const_iterator i = well_index.find(w);

	if(Model::well(w).escape_size() && i != well_index.end())
	  //
	  km(i->second, i->second) += (dd_t)Model::well(w).escape_rate(energy_bin(e));
      }

      // relaxation eigenvalues
      //

#ifdef WITH_MPACK
      
      Mpack_dd::Matrix evec(well_array.size());

      Mpack_dd::Vector eval = km.eigenvalues(&evec);

#else

      Lapack::Matrix evec(well_array.size());

      Lapack::Vector eval = km.eigenvalues(&evec);

#endif
      
      kinetic_basis[e].eigenvalue  = eval;
      
      kinetic_basis[e].eigenvector = evec;
      
      // kinetically active subspace
      //
      for(itemp = 0; itemp < well_array.size(); ++itemp) {
	//
	const int& w = well_array[itemp];

	dtemp = well(w).collision_frequency() * reduction_threshold;
	
	if(eval[itemp] > (dd_t)dtemp)
	  //
	  break;
      }
      
      kinetic_basis[e].active_size = itemp;
      //
    }//energy cycle
    //
  }// kinetically active basis


  // global united species indexing
  //
  std::vector<int> well_shift(ener_index_max);
  
  itemp = 0;
  
  for(int e = 0; e < ener_index_max; itemp += kinetic_basis[e++].active_size)
    //
    well_shift[e] = itemp;
  
  const int global_size = itemp;

  // global kinetic matrix size without reduction
  //
  itemp = 0;
  //
  for(int w = 0; w < Model::well_size(); ++w)
    //
    itemp += well(w).size();

  IO::log << IO::log_offset << "original kinetic matrix size = " << itemp << "\n";
  
  IO::log << IO::log_offset << "reduced  kinetic matrix size = " << global_size << "\n";

  /************************************************************************************* 
   ********************************* KINETIC MATRIX ************************************
   *************************************************************************************/

  // kinetic relaxation matrix
  //
  Lapack::SymmetricMatrix kin_mat(global_size);
  
  kin_mat = 0.;
  
  // diagonal chemical relaxation
  //
  for(int e = 0; e < ener_index_max; ++e) {
    //
    for(int l = 0; l < kinetic_basis[e].active_size; ++l) {
      //
      itemp = well_shift[e] + l;
      
      kin_mat(itemp, itemp) = kinetic_basis[e].eigenvalue[l];
    }
  }
    
  // collisional energy relaxation
  //
  for(int e1 = 0; e1 < ener_index_max; ++e1) {
    //
    for(int e2 = e1; e2 < ener_index_max; ++e2) {
      //
      for(int l1 = 0; l1 < kinetic_basis[e1].active_size; ++l1) {
	//
	for(int l2 = 0; l2 < kinetic_basis[e2].active_size; ++l2) {
	  //
	  if(e1 == e2 && l2 < l1)
	    //
	    continue;

	  dtemp = 0.;
	  
	  for(int w = 0; w < Model::well_size(); ++w) {
	    //
	    std::map<int, int>::const_iterator i1 = kinetic_basis[e1].well_index_map.find(w);

	    std::map<int, int>::const_iterator i2 = kinetic_basis[e2].well_index_map.find(w);
	    
	    itemp = e2 - e1;
	    
	    if(i1 != kinetic_basis[e1].well_index_map.end() &&
	       //
	       i2 != kinetic_basis[e2].well_index_map.end() &&
	       //
	       itemp < well(w).kernel_bandwidth) {
	      //
	      dtemp +=  well(w).kernel(e1, e2) * well(w).collision_frequency()
		//
		* well(w).boltzman_sqrt(e1) / well(w).boltzman_sqrt(e2)
		//
		* kinetic_basis[e1].eigenvector(i1->second, l1)
		//
		* kinetic_basis[e2].eigenvector(i2->second, l2);
	    }
	  }
	  
	  kin_mat(well_shift[e1] + l1, well_shift[e2] + l2) += dtemp;
	  //
	}// l2 cycle
	//
      }// l1 cycle
      //
    }// e2 cycle
    //
  }// e1 cycle

  /***********************************************************************
   ************************* GLOBAL MATRICES *****************************
   ***********************************************************************/

  Lapack::Vector eigenval;
  
  Lapack::Matrix global_eigen(global_size);

  {
    IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);
    //
#ifdef WITH_MPACK
    
    switch(float_type) {
      //
    case DOUBLE:
      //
      eigenval = kin_mat.eigenvalues(&global_eigen);

      break;

    case DD:
      //
      eigenval = Mpack::dd_eigenvalues(kin_mat, &global_eigen);

      break;
	
    case QD:
      //
      eigenval = Mpack::qd_eigenvalues(kin_mat, &global_eigen);

      break;

    default:
      //
      std::cerr << funame << "unknown case\n";

      throw Error::Logic();
    }
    
#else
    
    eigenval = kin_mat.eigenvalues(&global_eigen);

#endif
	
  }

  const double& relax_eval_min = eigenval[Model::well_size()];
  
  const double& relax_eval_max = eigenval.back();

  IO::log << IO::log_offset << "minimal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_min / well(0).collision_frequency() << "\n";
  
  IO::log << IO::log_offset << "maximal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_max / well(0).collision_frequency() << "\n";

  // eigenvectors projection on thermal subspace;
  //
  Lapack::Matrix eigen_pop(global_size, Model::well_size());

  for(int l = 0; l < global_size; ++l) {
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      double prod = 0.;
      
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

	if(wit == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "well is not in the kinetic basis space\n";

	  throw Error::Logic();
	}

	dtemp = 0.;
	
	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(wit->second, k);

	prod += dtemp * well(w).boltzman_sqrt(e);
      }
      
      eigen_pop(l, w) = prod / well(w).weight_sqrt();
    }//
    //
  }// thermal
  
  // eigenvectors projection on the bimolecular subspace;
  //
  Lapack::Matrix eigen_bim;

  if(Model::bimolecular_size()) {
    //
    eigen_bim.resize(global_size, Model::bimolecular_size());

    eigen_bim = 0.;
  
    for(int l = 0; l < global_size; ++l)
      //
      for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	//
	const int& w = Model::outer_connect(b).first;
      
	const int& p = Model::outer_connect(b).second;

	for(int e = 0; e < outer_barrier(b).size(); ++e) {
	  //
	  std::map<int, int>::const_iterator cit = kinetic_basis[e].well_index_map.find(w);

	  if(cit == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "well is not in the kinetic basis space\n";

	    throw Error::Logic();
	  }

	  dtemp = 0.;
	
	  for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	    //
	    dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(cit->second, k);
	  }

	  eigen_bim(l, p) += dtemp * outer_barrier(b).state_number(e) / 2. / M_PI * thermal_factor(e) / well(w).boltzman_sqrt(e);
	}
      }
    //
  }// bimolecular
  
  // eigenvectors projection on the escape subspace;
  //
  Lapack::Matrix eigen_escape;

  if(Model::escape_size()) {
    //
    eigen_escape.resize(global_size, Model::escape_size());

    eigen_escape = 0.;
  
    for(int l = 0; l < global_size; ++l)
      //
      for(int s = 0; s < Model::escape_size(); ++s) {
	//
	const int& w = Model::escape_channel(s).first;

	const int& c = Model::escape_channel(s).second;
      
	for(int e = 0; e < well(w).size(); ++e) {
	  //
	  std::map<int, int>::const_iterator cit = kinetic_basis[e].well_index_map.find(w);

	  if(cit == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "well is not in the kinetic basis space\n";

	    throw Error::Logic();
	  }

	  dtemp = 0.;
	
	  for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	    //
	    dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(cit->second, k);
	  }

	  eigen_escape(l, s) += dtemp * Model::well(w).escape_rate(energy_bin(e), c)  * well(w).boltzman_sqrt(e);
	}
      }
    //
  }// escape

  // eigenvectors at hot energies
  //
  Lapack::Matrix eigen_hot;
  
  if(hot_index.size()) {
    //
    eigen_hot.resize(global_size, hot_index.size());
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
      //
      const int& w = hit->first;

      const int& e = hit->second;
	
      std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

      for(int l = 0; l < global_size; ++l) {
	//
	dtemp = 0.;
	
	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(wit->second, k);

	dtemp /= well(w).boltzman_sqrt(e);
	  
	eigen_hot(l, count) = dtemp;
      }//
      //
    }//
    //
  }// hot energies

  /************************************* EIGENVECTOR OUTPUT *****************************************/

  std::vector<double> relaxation_projection(Model::well_size());
  
  for(int l = 0; l < Model::well_size(); ++l)
    //
    relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));

 
  IO::log << std::setprecision(3)
    //
	  << IO::log_offset << "eigenstate populations:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(10) << "*R"
    //
	  << std::setw(10) << "*P";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  IO::log << "\n";

  for(int l = 0; l < Model::well_size(); ++l) {
    //
    double nfac;
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = std::fabs(eigen_pop(l, w) * well(w).weight_sqrt());

      if(!w || dtemp > nfac)
	//
	nfac = dtemp;
    }
    
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(10) << eigenval[l] / relax_eval_min
      //
	    << std::setw(10) << relaxation_projection[l];
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well(w).weight_sqrt() / nfac;
      
      IO::log << std::setw(10);
      
      if(dtemp > .01 || dtemp < -.01) {
	//
	IO::log << dtemp;
      }
      else
	//
	IO::log << 0;
    }
    
    IO::log << "\n";
  }

  IO::log << IO::log_offset
    //
	  << std::setw(5)  << "*Z"
    //
	  << std::setw(10) << "---"
    //
	  << std::setw(10) << "---";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    if(!w || dtemp < well(w).weight_sqrt())
      //
      dtemp = well(w).weight_sqrt();
  
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << well(w).weight_sqrt() / dtemp;
  
  IO::log << "\n";

  IO::log << IO::log_offset
    //
	  << "*R - eigenvalue over the relaxation limit\n"
    //
	  << IO::log_offset
    //
	  << "*P - eigenvector projection squared on the relaxation subspace (= 1 -F_ne)\n"
    //
	  << IO::log_offset
    //
	  << "*Z - well partition function square root (normalized)\n"
    //
	  << std::setprecision(6);
  
  /***************************************************************
   **************** CHEMICAL SUBSPACE DIMENSION ******************
   ***************************************************************/

  // predefined chemical subspace dimension
  //
  if(_default_chem_size >= 0) {
    //
    itemp = _default_chem_size;
  }
  //
  // default partitioning scheme
  //
  else if(default_partition.size()) {
    //
    itemp = default_partition.size();
  }
  //
  // absolute eigenvalue threshold
  //
  else if(chemical_threshold > 1.) {
    //
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(relax_eval_min / eigenval[itemp]  < chemical_threshold)
	//
	break;
  }
  //
  // relaxation projection threshold
  //
  else if(chemical_threshold < 1. && chemical_threshold > 0.) {
    //
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(relaxation_projection[itemp] > chemical_threshold)
	//
	break;
  }
  //
  // relative eigenvalue threshold
  //
  else if(chemical_threshold < -1. ) {
    //
    for(itemp = Model::well_size(); itemp > 0; --itemp)
      //
      if(eigenval[itemp] / eigenval[itemp - 1] > - chemical_threshold)
	//
	break;
  }
  //
  else
    //
    itemp = Model::well_size();
	
  const int chem_size = itemp;

  IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n";

  /***************************************************************
   ***************** HIGH EIGENSTATES CORRECTION *****************
   ***************************************************************/

  for(int l = 0; l < chem_size; ++l) {
    //
    for(int w1 = 0; w1 < Model::well_size(); ++w1) {
      //
      for(int e1 = 0; e1 < well(w1).size(); ++e1) {
	//
	std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);
	  
	dtemp = 0.;
	
	for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e1] + k, l) * kinetic_basis[e1].eigenvector(i11->second, k);
	 
	const double proj = dtemp;

	itemp = e1 + well(w1).kernel_bandwidth;

	const int e2_max = itemp < well(w1).size() ? itemp : well(w1).size();
	  
	itemp = e1 - well(w1).kernel_bandwidth + 1;

	const int e2_min = itemp > 0 ? itemp : 0;
	    
	for(int e2 = e2_min; e2 < e2_max; ++e2) {
	  //
	  std::map<int, int>::const_iterator i21 = kinetic_basis[e2].well_index_map.find(w1);

	  if(i21 == kinetic_basis[e2].well_index_map.end()) {
	    //
	    std::cerr << funame << "well_reduction_correction: logic error 2\n";

	    throw Error::Logic();
	  }
	  // population correction
	  //
	  for(int w2 = 0; w2 < Model::well_size(); ++w2) {
	    //
	    std::map<int, int>::const_iterator i22 = kinetic_basis[e2].well_index_map.find(w2);

	    if(i22 == kinetic_basis[e2].well_index_map.end())
	      //
	      continue;

	    dtemp = 0.;
	      
	    for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
	      //
	      dtemp += kinetic_basis[e2].eigenvector(i21->second, k)
		//
		* kinetic_basis[e2].eigenvector(i22->second, k)
		//
		/ kinetic_basis[e2].eigenvalue[k];

	    dtemp *= proj * well(w1).kernel(e1, e2) * well(w1).collision_frequency()
	      //
	      * well(w1).boltzman_sqrt(e1) / well(w1).boltzman_sqrt(e2);
	    
	    eigen_pop(l, w2) -= dtemp * well(w2).boltzman_sqrt(e2) / well(w2).weight_sqrt();
	    //
	  }// second well cycle
	  //
	  // bimolecular correction
	  //
	  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	    //
	    if(e2 >= outer_barrier(b).size())
	      //
	      continue;

	    const int& w2 = Model::outer_connect(b).first;
		
	    const int& p  = Model::outer_connect(b).second;
		
	    std::map<int, int>::const_iterator i22 = kinetic_basis[e2].well_index_map.find(w2);
		  
	    if(i22 == kinetic_basis[e2].well_index_map.end()) {
	      //
	      std::cerr << funame << "well_reduction_correction: logic error 3\n";

	      throw Error::Logic();
	    }

	    dtemp = 0.;
	      
	    for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
	      //
	      dtemp += kinetic_basis[e2].eigenvector(i21->second, k)
		//
		* kinetic_basis[e2].eigenvector(i22->second, k)
		//
		/ kinetic_basis[e2].eigenvalue[k];

	    dtemp *= proj * well(w1).kernel(e1, e2) * well(w1).collision_frequency()
	      //
	      * well(w1).boltzman_sqrt(e1) / well(w1).boltzman_sqrt(e2);
	    
	    eigen_bim(l, p) -= dtemp * outer_barrier(b).state_number(e2) / 2. / M_PI
	      //
	      * thermal_factor(e2) / well(w2).boltzman_sqrt(e2);
	    //
	  }// outer barrier cycle
	  //
	  // escape correction
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    const int& w2 = Model::escape_channel(s).first;

	    const int& c = Model::escape_channel(s).second;
      
	    std::map<int, int>::const_iterator i22 = kinetic_basis[e2].well_index_map.find(w2);

	    if(i22 == kinetic_basis[e2].well_index_map.end())
	      //
	      continue;

	    dtemp = 0.;
	
	    for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
	      //
	      dtemp += kinetic_basis[e2].eigenvector(i21->second, k)
		//
		* kinetic_basis[e2].eigenvector(i22->second, k)
		//
		/ kinetic_basis[e2].eigenvalue[k];

	    dtemp *= proj * well(w1).kernel(e1, e2) * well(w1).collision_frequency()
	      //
	      * well(w1).boltzman_sqrt(e1) / well(w1).boltzman_sqrt(e2);
	    
	    eigen_escape(l, s) -= dtemp * Model::well(w2).escape_rate(energy_bin(e2), c)  * well(w2).boltzman_sqrt(e2);
	    //
	  }// escape cycle
	  //
	}// second energy cycle
	//
	// hot energies correction
	//
	int count = 0;
	    
	for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
	  //
	  const int& w2 = hit->first;
	      
	  const int& e2 = hit->second;

	  if(e2 < e2_min || e2 >= e2_max)
	    //
	    continue;

	  std::map<int, int>::const_iterator i21 = kinetic_basis[e2].well_index_map.find(w1);

	  std::map<int, int>::const_iterator i22 = kinetic_basis[e2].well_index_map.find(w2);

	  dtemp = 0.;
	      
	  for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
	    //
	    dtemp += kinetic_basis[e2].eigenvector(i21->second, k)
	      //
	      * kinetic_basis[e2].eigenvector(i22->second, k)
	      //
	      / kinetic_basis[e2].eigenvalue[k];
	      
	    dtemp *= proj * well(w1).kernel(e1, e2) * well(w1).collision_frequency()
	      //
	      * well(w1).boltzman_sqrt(e1) / well(w1).boltzman_sqrt(e2);
	    
	  eigen_hot(l, count) -= dtemp / well(w2).boltzman_sqrt(e2);
	  //
	}// hot energies cycle
	//
      }// first energy cycle
      //
    }// first well cycle
    //
  }// eigenstate cycle
    
  /***************** BIMOLECULAR RATE CONSTANTS ******************/

  // collisional relaxation eigenvalues and eigenvectors
  //
  const int relax_size = global_size - chem_size;
  
  Lapack::Vector relax_lave(relax_size);
  
  for(int r = 0; r < relax_size; ++r) {
    //
    itemp = r + chem_size;

    relax_lave[r] = 1. / eigenval[itemp];
  }
      
  // bimolecular-to-bimolecular rate coefficients
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  Lapack::SymmetricMatrix bb_rate;
  
  if(Model::bimolecular_size()) {
    //
    bb_rate.resize(Model::bimolecular_size());
    
    for(int i = 0; i < Model::bimolecular_size(); ++i)
      //
      for(int j = i; j < Model::bimolecular_size(); ++j)
	//
    	bb_rate(i, j) = triple_product(&eigen_bim(chem_size, i), &eigen_bim(chem_size, j), relax_lave, relax_size);

    // high eigenvalue contribution
    //
    vtemp.resize(Model::bimolecular_size());
    
    for(int e = 0; e < ener_index_max; ++e)
      //
      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	//
	vtemp = 0.;
	
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	  //
	  const int p = Model::outer_connect(b).second;

	  if(e >= outer_barrier(b).size())
	    //
	    continue;
	  
	  std::map<int, int>::const_iterator bit = kinetic_basis[e].well_index_map.find(Model::outer_connect(b).first);

	  if(bit == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "bimolecular-to-bimolecular: logic error 1\n";
	    
	    throw Error::Logic();
	  }
	    
	  vtemp[p] += kinetic_basis[e].eigenvector(bit->second, k)
	    //
	    * outer_barrier(b).state_number(e) / 2. / M_PI
	    //
	    * thermal_factor(e) / well(bit->first).boltzman_sqrt(e);
	}
	
	for(int i = 0; i < Model::bimolecular_size(); ++i)
	  //
	  for(int j = i; j < Model::bimolecular_size(); ++j)
	    //
	    bb_rate(i, j) += vtemp[i] * vtemp[j] / kinetic_basis[e].eigenvalue[k];
      }

    
    for(int i = 0; i < Model::bimolecular_size(); ++i) {
      //
      if(Model::bimolecular(i).dummy())
	//
	continue;
      
      for(int j = 0; j < Model::bimolecular_size(); ++j) {
	//
	// bimolecular reactant loss
	//
	if(i == j) {
	  //
	  dtemp = 0.;
	    
	  for(int b = 0; b < Model::outer_barrier_size(); ++b)
	    //
	    if(Model::outer_connect(b).second == i)
	      //
	      for(int e = 0; e < outer_barrier(b).size(); ++e)
		//
		dtemp += outer_barrier(b).state_number(e) * thermal_factor(e);
	    
	  dtemp /= 2. * M_PI;
	    
	  dtemp -= bb_rate(i, i);
	}
	// crossrate
	//
	else
	  //
	  dtemp = bb_rate(i, j);

	dtemp *= energy_step() / bimolecular(i).weight() / bru;
	  
	rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 
      }
    }
    
    // bimolecular-to-escape rate coefficients
    //
    for(int s = 0; s < Model::escape_size(); ++s) {
      //
      const int& w = Model::escape_channel(s).first; // escape well

      const int& c = Model::escape_channel(s).second; // escape channel
	    
      vtemp = 0.;
	    
      // high eigenvalue contribution
      //
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator sit = kinetic_basis[e].well_index_map.find(w);

	if(sit == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "bimolecular-to-escape: logic error 1\n";
	    
	  throw Error::Logic();
	}

	for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	  //
	  const int& p = Model::outer_connect(b).second;

	  if(e >= outer_barrier(b).size() || Model::bimolecular(p).dummy())
	    //
	    continue;

	  std::map<int, int>::const_iterator bit = kinetic_basis[e].well_index_map.find(Model::outer_connect(b).first);
	    
	  if(bit == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "bimolecular-to-escape: logic error 2\n";
	      
	    throw Error::Logic();
	  }

	  dtemp = 0.;
	    
	  for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	    //
	    dtemp += kinetic_basis[e].eigenvector(sit->second, k)
	      //
	      * kinetic_basis[e].eigenvector(bit->second, k)
	      //
	      / kinetic_basis[e].eigenvalue[k];
	  }

	  dtemp *= Model::well(sit->first).escape_rate(energy_bin(e), c)
	    //
	    * well(sit->first).boltzman_sqrt(e)
	    //
	    * outer_barrier(b).state_number(e) / 2. / M_PI
	    //
	    * thermal_factor(e) / well(bit->first).boltzman_sqrt(e);

	  vtemp[p] += dtemp;
	}
      }
	
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	if(Model::bimolecular(p).dummy())
	  //
	  continue;
	  
	dtemp = vtemp[p] + triple_product(&eigen_bim(chem_size, p), &eigen_escape(chem_size, s), relax_lave, relax_size);

	dtemp *= energy_step() / bimolecular(p).weight() / bru;
	    
	rate_data[std::make_pair(Model::well_size() + p, Model::well_size() + Model::bimolecular_size() + s)] = dtemp; 
      }
    }
  }

  /***** PARTITIONING THE GLOBAL PHASE SPACE INTO THE CHEMICAL AND COLLISIONAL SUBSPACES *****/

  std::vector<int>  group_index;
  
  Lapack::Matrix m_direct;
  
  std::vector<double>      weight;
  
  std::vector<double> real_weight;

  if(chem_size) {
    //
    // projection of the chemical eigenvectors onto the thermal subspace
    //
    Lapack::Matrix pop_chem(Model::well_size(), chem_size);
    
    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

#ifdef DEBUG

    IO::log << IO::log_offset << "orthogonality check starts\n";
      
    IO::log_offset.increase();

    proj = 0.;
      
    for(int l = 0; l < chem_size; ++l) {
      //
      // normalize
      //
      //normalize(&pop_chem(0, l), Model::well_size());
      //
      for(int m = 0; m < l; ++m) {
	//
	dtemp = vdot(pop_chem.column(l), pop_chem.column(m));
	  
	dtemp = dtemp >= 0. ? dtemp : -dtemp;
	  
	proj = dtemp > proj ? dtemp : proj;
      }
    }

    IO::log << IO::log_offset << "maximal scalar product of different chemical eigenvectors = "
      //
	    << proj << "\n";

    for(int l = 0; l < chem_size; ++l) {
      //
      dtemp = vdot(pop_chem.column(l));
	
      if(!l || dtemp < proj)
	//
	proj = dtemp;
    }

    IO::log << IO::log_offset << "minimal chemical eigenvector square = "
      //
	    << proj << "\n";

    for(int l = 0; l < chem_size; ++l) {
      //
      dtemp = vdot(pop_chem.column(l));
	
      if(!l || dtemp > proj)
	//
	proj = dtemp;
    }

    IO::log << IO::log_offset << "maximal chemical eigenvector square = "
      //
	    << proj << "\n";

    IO::log_offset.decrease();
      
    IO::log << IO::log_offset << "orthogonality check done\n";
    

#endif

    // partitioning wells into equilibrated groups
    //
    if(default_partition.size()) {
      //
      // default reduction scheme
      //
      well_partition = default_partition;

      IO::log << IO::log_offset << "using default reduction scheme, projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    
      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }
    else if(chem_size == Model::well_size()) {
      //
      // no partitioning
      //
      well_partition.resize(Model::well_size());
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	well_partition[w].insert(w);

      IO::log << IO::log_offset << "projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    }
    else {
      //
      // well partitioning
      //
      Group bimolecular_group;
      
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group, std::vector<double>());

      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }

    group_index = well_partition.group_index();
    weight      = well_partition.weight();
    real_weight = well_partition.real_weight();

    // auxiliary output
    //
    IO::aux << "number of species = " << well_partition.size() << "\n";
    
    for(int p = 0; p < well_partition.size(); ++p) {
      //
      std::multimap<double, int> ww_map;
      
      for(Git git = well_partition[p].begin(); git != well_partition[p].end(); ++git)
	//
	ww_map.insert(std::make_pair(well(*git).real_weight(),*git));

      for(std::multimap<double, int>::const_reverse_iterator mit = ww_map.rbegin(); mit != ww_map.rend(); ++mit)
	//
	IO::aux << std::setw(15) << Model::well(mit->second).name();

      IO::aux << "\n";
      
      for(std::multimap<double, int>::const_reverse_iterator mit = ww_map.rbegin(); mit != ww_map.rend(); ++mit)
	//
	IO::aux << std::setw(15) << mit->first;
      
      IO::aux << "\n";
    }
    // output
    //
    if(chem_size != Model::well_size()) {
      //
      IO::log << IO::log_offset << "species:\n"
	//
	      << IO::log_offset << std::setw(2) << "#"  << std::setw(15)
	//
	      << "assigned name" << IO::first_offset
	//
	      << "group\n";
      
      for(Pit g = well_partition.begin(); g != well_partition.end(); ++g) {
	//
	itemp = g - well_partition.begin();
	
	IO::log << IO::log_offset << std::setw(2) << itemp << std::setw(15)
	  //
		<< Model::well(group_index[itemp]).name() << IO::first_offset;
	
	for(Git w = g->begin(); w != g->end(); ++w) {
	  //
	  if(w != g->begin())
	    //
	    IO::log << "+";
	  
	  IO::log << Model::well(*w).name();
	}
	
	IO::log << "\n";
      }
    }

    m_direct = pop_chem;
    
    Lapack::Matrix m_inverse = m_direct.invert();

#ifdef DEBUG
    
    Lapack::Matrix one(m_direct * m_inverse);
    
    one.diagonal() -= 1.;
    
    double val_max = -1.;
    
    for(int i = 0; i < one.size1(); ++i)
      //
      for(int j = 0; j < one.size2(); ++j) {
	//
	dtemp = one(i, j);
	
	dtemp = dtemp < 0. ? -dtemp : dtemp;
	
	if(dtemp > epsilon && dtemp > val_max)
	  //
	  val_max = dtemp;
      }
    
    if(val_max > 0.)
      //
      IO::log << IO::log_offset << funame
	//
	      << "WARNING: matrix inversion error = " << val_max
	//
	      << " exceeds numerical accuracy = " << epsilon
	//
	      << "\n";

#endif
    
    // well-to-well rate coefficients
    //
    Lapack::Matrix ww_rate(chem_size);
    
    ww_rate = 0.;
    
    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j)
	//
	for(int l = 0; l < chem_size; ++l)
	  //
	  ww_rate(i, j) += m_direct(j, l) * m_inverse(l, i) * eigenval[l];

    // well-to-bimolecular rate coefficients
    //
    Lapack::Matrix wb_rate, bw_rate;
    
    if(Model::bimolecular_size()) {
      //
      wb_rate.resize(chem_size, Model::bimolecular_size());
      
      for(int w = 0; w < chem_size; ++w)
	//
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  wb_rate(w, p) = vdot(m_inverse.column(w), &eigen_bim(0, p));

      // bimolecular-to-well rate coefficients
      //
      bw_rate.resize(Model::bimolecular_size(), chem_size);
      
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int w = 0; w < chem_size; ++w)
	  //
	  bw_rate(p, w) = vdot(m_direct.row(w), &eigen_bim(0, p));
    }
  
    // output
    //
    IO::log << std::setprecision(2)
      //
	    << IO::log_offset << "Wa->Wb/Wb->Wa rate constants ratios:\n"
      //
	    << IO::log_offset << std::setw(5) << "Wb\\Wa";
    
    for(int i = 0; i < chem_size; ++i)
      //
      IO::log << std::setw(10) << Model::well(group_index[i]).name();
    
    IO::log << "\n";
    
    for(int j = 0; j < chem_size; ++j) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(group_index[j]).name();
      
      for(int i = 0; i < chem_size; ++i)
	//
	if(i != j) {
	  //
	  if(ww_rate(j, i) != 0.) {
	    //
	    IO::log << std::setw(10) << ww_rate(i, j) / ww_rate(j, i);
	  }
	  else
	    //
	    IO::log << std::setw(10) << "***";
	}
	else
	  //
	  IO::log << std::setw(10) << "1";
      
      IO::log << "\n";
    }
    
    if(Model::bimolecular_size()) {
      //
      IO::log << IO::log_offset << "W->P/P->W rate constants ratios:\n"
	//
	      << IO::log_offset << std::setw(5) << "P\\W";
      
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(10) << Model::well(group_index[w]).name();
      
      IO::log << "\n";
    
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).name();
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  if(bw_rate(p, w) != 0.) {
	    //
	    IO::log << std::setw(10) << wb_rate(w, p) / bw_rate(p, w);
	  }
	  else
	    //
	    IO::log << std::setw(10) << "***";
	
	IO::log << "\n";
      }
    }

    IO::log << std::setprecision(6);

    // well-to-well rate coefficients
    //
    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j) {
	//
	dtemp = ww_rate(i, j) * std::sqrt(weight[i] * weight[j]) * energy_step() / real_weight[i] / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp; 	
      }

    // well-to-escape rate coefficients
    //
    if(Model::escape_size())
      //
      for(int w = 0; w < chem_size; ++w)
	//
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  dtemp = vdot(m_inverse.column(w), &eigen_escape(0, e))
	    //
	    * std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() +
				   //
				   Model::bimolecular_size() + e)] = dtemp;
	}

    // well-to-bimolecular rate coefficients
    //
    for(int w = 0; w < chem_size; ++w)
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = wb_rate(w, p) * std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	
	rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
      }
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(bimolecular(p).weight() > 0.)
	//
	for(int w = 0; w < chem_size; ++w) {
	  //
	  dtemp = bw_rate(p, w) * std::sqrt(weight[w]) * energy_step() / bimolecular(p).weight() / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	}
  }

  // hot energies branching ratios
  //
  if(hot_index.size()) {
    //
    IO::log << IO::log_offset << "hot energies branching fractions:\n"
      //
	    << IO::log_offset
      //
	    << std::setw(5)  << "Well"
      //
	    << std::setw(13) << "E, kcal/mol";
    
    for(int w = 0; w < chem_size; ++w)
      //
      IO::log << std::setw(13) << Model::well(group_index[w]).name();

    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(13) << Model::bimolecular(p).name();

    for(int s = 0; s < Model::escape_size(); ++s)
      //
      IO::log << std::setw(13) << Model::escape_name(s);
    
    IO::log << std::setw(13) << "total" << "\n";
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
      //
      const int w = hit->first;

      const int e = hit->second;

      std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

      if(wit == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "hot energies: logic error 1\n";

	throw Error::Logic();
      }
      
      IO::log << IO::log_offset
	//
	      << std::setw(5)  << Model::well(w).name()
	//
	      << std::setw(13) << energy_bin(e) / Phys_const::kcal;

      double tot = 0.;
      
      // bound species
      //
      for(int c = 0; c < chem_size; ++c) {
	//
	dtemp = vdot(m_direct.row(c), &eigen_hot(0, count)) * std::sqrt(weight[c]);
	
	IO::log << std::setw(13) << dtemp;

	tot += dtemp;
      }
      
      // bimolecular products
      //
      vtemp.resize(Model::bimolecular_size());

      vtemp = 0.;
      
      for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	//
	if(e >= outer_barrier(b).size())
	  //
	  continue;

	const int p = Model::outer_connect(b).second;
	
	std::map<int, int>::const_iterator bit = kinetic_basis[e].well_index_map.find(Model::outer_connect(b).first);

	if(bit == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "hot energies: logic error 2\n";

	  throw Error::Logic();
	}
      
	dtemp = 0.;
	
	for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
	  //
	  dtemp += kinetic_basis[e].eigenvector(wit->second, k)
	    //
	    * kinetic_basis[e].eigenvector(bit->second, k)
	    //
	    / kinetic_basis[e].eigenvalue[k];

	vtemp[p] += dtemp / well(wit->first).boltzman_sqrt(e)
	  //
	  * outer_barrier(b).state_number(e) / 2. / M_PI
	  //
	  * thermal_factor(e) / well(bit->first).boltzman_sqrt(e);
      }
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = vtemp[p]
	  //
	  + triple_product(&eigen_bim(chem_size, p), &eigen_hot(chem_size, count), relax_lave, relax_size);
	
	IO::log << std::setw(13) << dtemp;

	tot += dtemp;
      }

      // escape products
      //
      for(int s = 0; s < Model::escape_size(); ++s) {
	//
	std::map<int, int>::const_iterator sit = kinetic_basis[e].well_index_map.find(Model::escape_channel(s).first);

	dtemp = 0.;
	
	if(sit != kinetic_basis[e].well_index_map.end()) {
	  //
	  for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	    //
	    dtemp += kinetic_basis[e].eigenvector(wit->second, k)
	      //
	      * kinetic_basis[e].eigenvector(sit->second, k)
	      //
	      / kinetic_basis[e].eigenvalue[k];
	  }

	  dtemp *= Model::well(sit->first).escape_rate(energy_bin(e), Model::escape_channel(s).second)
	    //
	    * well(sit->first).boltzman_sqrt(e) / well(wit->first).boltzman_sqrt(e);
	}
	
	dtemp += triple_product(&eigen_escape(chem_size, s), &eigen_hot(chem_size, count), relax_lave, relax_size);
	
	IO::log << std::setw(13) << dtemp;

	tot += dtemp;
      }
      
      IO::log << std::setw(13) << tot << "\n";
    }
  }

  // prompt isomerization
  //
  IO::log << IO::log_offset << "prompt isomerization/dissociation:\n";

  IO::log << IO::log_offset << std::setw(7) << "W\\W,P,E";

  if(chem_size)
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(13) << Model::well(w).name();

  IO::log << std::setw(13) << "Total Diss.";

  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    IO::log << std::setw(13) << Model::bimolecular(p).name();

  for(int e = 0; e < Model::escape_size(); ++e)
    //
    IO::log << std::setw(13) << Model::escape_name(e);
    
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(7) << Model::well(w).name();

    double diss = 1.;

    if(chem_size)
      //
      for(int v = 0; v < Model::well_size(); ++v) {
	//
	dtemp = 0.;

	for(int l = 0; l < chem_size; ++l)
	  //
	  dtemp += eigen_pop(l, w) * eigen_pop(l, v);

	// renormalization
	//
	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	dtemp *= well(v).weight_sqrt() * well(w).weight_sqrt() * energy_step() / well(w).real_weight();

	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	diss -= dtemp;

	IO::log << std::setw(13) << dtemp;
      }

    IO::log << std::setw(13) << diss;

    // bimolecular products
    //
    vtemp.resize(Model::bimolecular_size());

    vtemp = 0.;

    for(int e = 0; e < well(w).size(); ++e) {
      //
      std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

      if(wit == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "prompt isomerization: logic error 1\n";

	throw Error::Logic();
      }
      
      for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	//
	if(e >= outer_barrier(b).size())
	  //
	  continue;

	const int p = Model::outer_connect(b).second;
	
	std::map<int, int>::const_iterator bit = kinetic_basis[e].well_index_map.find(Model::outer_connect(b).first);

	if(bit == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "prompt isomerizaion: logic error 2\n";

	  throw Error::Logic();
	}
      
	dtemp = 0.;
	
	for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
	  //
	  dtemp += kinetic_basis[e].eigenvector(wit->second, k)
	    //
	    * kinetic_basis[e].eigenvector(bit->second, k)
	    //
	    / kinetic_basis[e].eigenvalue[k];

	vtemp[p] += dtemp * well(wit->first).boltzman_sqrt(e)
	  //
	  * outer_barrier(b).state_number(e) / 2. / M_PI
	  //
	  * thermal_factor(e) / well(bit->first).boltzman_sqrt(e);
      }
    }
    
    for(int p = 0; p < Model::bimolecular_size(); ++p) {
      //
      dtemp = vtemp[p] + well(w).weight_sqrt()
	//
	* triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w), relax_lave, relax_size);

      IO::log << std::setw(13) << dtemp	* energy_step() / well(w).real_weight();
    }

    // escape products
    //
    for(int s = 0; s < Model::escape_size(); ++s) {
      //
      double corr = 0.;
    
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);
	
	if(wit == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "prompt isomerization: logic error 3\n";
	  
	  throw Error::Logic();
	}
      
	std::map<int, int>::const_iterator sit = kinetic_basis[e].well_index_map.find(Model::escape_channel(s).first);

	if(sit != kinetic_basis[e].well_index_map.end()) {
	  //
	  dtemp = 0.;
	
	  for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	    //
	    dtemp += kinetic_basis[e].eigenvector(wit->second, k)
	      //
	      * kinetic_basis[e].eigenvector(sit->second, k)
	      //
	      / kinetic_basis[e].eigenvalue[k];
	  }

	  dtemp *= Model::well(sit->first).escape_rate(energy_bin(e), Model::escape_channel(s).second)
	    //
	    * well(sit->first).boltzman_sqrt(e) * well(wit->first).boltzman_sqrt(e);

	  corr += dtemp;
	}
      }

      dtemp = corr + well(w).weight_sqrt()
	//
	* triple_product(&eigen_escape(chem_size, s), &eigen_pop(chem_size, w), relax_lave, relax_size);
      
      IO::log << std::setw(13) << dtemp * energy_step() / well(w).real_weight();
    }
    
    IO::log << "\n";
  }

  for(int w = 0; w < Model::well_size(); ++w)
    //
    if(!w || well(w).size() > itemp)
      //
      itemp = well(w).size();
  
  const int well_size_max = itemp;
    
  Lapack::Matrix eigen_well(global_size, Model::well_size());
  
  for(int l = 0; l < global_size; ++l)
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      double prod = 0.;
      
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

	if(wit == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "eigen_well: logic error\n";

	  throw Error::Logic();
	}

	dtemp = 0.;

	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(wit->second, k);

	prod += dtemp * dtemp;
      }
	
      eigen_well(l, w) = std::sqrt(prod);
    }
  
  // eigenvalues output
  //
  if(eval_out.is_open()) {
    //
    eval_out << std::setw(13) << temperature() / Phys_const::kelv
      //
	     << std::setw(13);
    
    switch(pressure_unit) {
      //
    case BAR:
      //
      eval_out << pressure() / Phys_const::bar;
      break;
      
    case TORR:
      //
      eval_out << pressure() / Phys_const::tor;
      break;
      
    case ATM:
      //
      eval_out << pressure() / Phys_const::atm;
      break;
    }
    
    eval_out << std::setw(13) << well(0).collision_frequency() / Phys_const::herz
      //
	     << std::setw(13) << relax_eval_min / well(0).collision_frequency();
    
    const int lmax = Model::well_size() + evec_out_num;
    
    for(int l = 0; l < lmax; ++l)
      //
      eval_out << std::setw(13) << eigenval[l] / well(0).collision_frequency()
	//
	       << std::setw(13) <<  1. - vdot(eigen_pop.row(l));
    
    eval_out << "\n";
  }
  
  // eigenvectors output
  //
  if(evec_out.is_open()) {
    //
    evec_out << "Pressure = ";
    
    switch(pressure_unit) {
    case BAR:
      evec_out << pressure() / Phys_const::bar << " bar";
      break;
    case TORR:
      evec_out << pressure() / Phys_const::tor << " torr";
      break;
    case ATM:
      evec_out << pressure() / Phys_const::atm << " atm";
      break;
    }
    
    evec_out << "\t Temperature = " << temperature() / Phys_const::kelv << " K\n";

    //
    evec_out << "EIGENVECTORS:\n";
    
    const int lmax = Model::well_size() + evec_out_num;
    
    for(int l = 0; l < lmax; ++l) {
      //
      evec_out << "l = " << l << "\n"
	//
	       << "eigenvalue / collision frequency = "
	//
	       << eigenval[l] / well(0).collision_frequency()
	//
	       << "\n";

      evec_out << std::setw(13) << "well length";
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << eigen_well(l, w);
      
      evec_out << "\n";

      evec_out << std::setw(13) << "E, kcal/mol";
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << Model::well(w).name() << std::setw(13) << Model::well(w).name();
      
      evec_out << "\n";

      for(int e = 0; e < well_size_max; ++e) {
	//
	evec_out << std::setw(13) << energy_bin(e) / Phys_const::kcal;
	
	for(int w = 0; w < Model::well_size(); ++w) {
	  //
	  std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

	  if(wit != kinetic_basis[e].well_index_map.end()) {
	    //
	    dtemp = 0.;

	    for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	      //
	      dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(wit->second, k);

	    evec_out << std::setw(13) << dtemp * well(w).boltzman_sqrt(e)
	      //
		     << std::setw(13) << dtemp / well(w).boltzman_sqrt(e) * well(w).weight_sqrt();
	  }
	  else
	    //
	    evec_out << std::setw(13) << "" <<  std::setw(13) << "";
	}
	
	evec_out << "\n";
      }
    }
    
    evec_out << "\n";
  }

  // products energy distributions
  //
  if(ped_out.is_open()) {
    //
    switch(pressure_unit) {
    case BAR:
      ped_out << "pressure[bar]        = " << pressure() / Phys_const::bar << "\n";
      break;
    case TORR:
      ped_out << "pressure[torr]       = " << pressure() / Phys_const::tor << "\n";
      break;
    case ATM:
      ped_out << "pressure[atm]        = " << pressure() / Phys_const::atm << "\n";
      break;
    }
    
    ped_out << "temperature[K]       = " << temperature() / Phys_const::kelv << "\n"
      //
	    << "energy step[1/cm]    = " << energy_step() / Phys_const::incm << "\n"
      //
	    << "max energy[kcal/mol] = " << energy_reference() / Phys_const::kcal << "\n\n";

    
    // bimolecular-to-bimolecular product energy distributions
    //
    if(Model::bimolecular_size()) {
      //
      ped_out << "Bimolecular-to-bimolecular product energy distributions:\n";

      // dimensions
      //
      itemp = 0;
      
      for(int count = 0; count < ped_pair.size(); ++count)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(Model::outer_connect(b).second == ped_pair[count].second && outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

      const int ener_index_max = itemp;
      
      Lapack::Matrix mtemp(itemp, ped_pair.size());

      mtemp = 0.;
      
      // distribution
      //
      for(int e = 0; e < ener_index_max; ++e) {
	//
	for(int count = 0; count < ped_pair.size(); ++count) {
	  //
	  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	    //
	    if(Model::outer_connect(b).second != ped_pair[count].second || e >= outer_barrier(b).size())
	      //
	      continue;
	    
	    std::map<int, int>::const_iterator bit = kinetic_basis[e].well_index_map.find(Model::outer_connect(b).first);

	    double ped = 0.;
	    
	    for(int l = chem_size; l < global_size; ++l) {
	      //
	      dtemp = 0.;

	      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
		//
		dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(bit->second, k);

	      dtemp *= eigen_bim(l, ped_pair[count].first) / eigenval[l];

	      ped += dtemp;
	    }

	    for(int c = 0; c < Model::outer_barrier_size(); ++c) {
	      //
	      if(Model::outer_connect(c).second != ped_pair[count].first || e >= outer_barrier(c).size())
		//
		continue;
		 
	      std::map<int, int>::const_iterator cit = kinetic_basis[e].well_index_map.find(Model::outer_connect(c).first);

	      dtemp = 0.;
	      
	      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
		//
		dtemp += kinetic_basis[e].eigenvector(bit->second, k)
		  //
		  * kinetic_basis[e].eigenvector(cit->second, k)
		  //
		  / kinetic_basis[e].eigenvalue[k];

	      dtemp *= outer_barrier(c).state_number(e) / 2. / M_PI
		//
		* thermal_factor(e) / well(cit->first).boltzman_sqrt(e);

	      ped += dtemp;
	    }

	    ped *= outer_barrier(b).state_number(e) / 2. / M_PI
	      //
	      * thermal_factor(e) / well(bit->first).boltzman_sqrt(e);
	      
	    mtemp(e, count) += ped;
	  }
	}
      }
      
      // normalization
      //
      for(int i = 0; i < mtemp.size2(); ++i)
	//
	mtemp.column(i) /= max(mtemp.column(i));

      // output
      //
      ped_out << std::setw(13) << "E, kcal/mol";
      
      for(int ped = 0; ped < ped_pair.size(); ++ ped)
	//
	ped_out << std::setw(13) << Model::bimolecular(ped_pair[ped].first).name() + "->"
	  //
	  + Model::bimolecular(ped_pair[ped].second).name();
      
      ped_out << "\n";

      for(int e = 0; e < ener_index_max; ++e) {
	//
	ped_out << std::setw(13) << energy_bin(e) / Phys_const::kcal;

	for(int i = 0; i < mtemp.size2(); ++i)
	  //
	  ped_out << std::setw(13) << mtemp(e, i);
	
	ped_out << "\n";
      }
      
      ped_out << "\n";
    }
  }

  IO::log << "\n";
}

/********************************************************************************************
 ************ THE DIRECT DIAGONALIZATION OF THE GLOBAL KINETIC RELAXATION MATRIX ************
 ********************************************************************************************/

void MasterEquation::direct_diagonalization_method (std::map<std::pair<int, int>, double>& rate_data,
						    //
						    Partition& well_partition, int flags)
  
{
  const char funame [] = "MasterEquation::direct_diagonalization_method: ";

  // bimolecular rate units
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  if(!isset()) {
    //
    std::cerr << funame << "reactive complex is not set\n";
    
    throw Error::Init();
  }

  IO::Marker funame_marker(funame);

  // pressure output
  //
  IO::log << IO::log_offset << "Pressure = ";
  
  switch(pressure_unit) {
  case BAR:
    IO::log << pressure() / Phys_const::bar << " bar";
    break;
  case TORR:
    IO::log << pressure() / Phys_const::tor << " torr";
    break;
  case ATM:
    IO::log << pressure() / Phys_const::atm << " atm";
    break;
  }

  // temperature output
  //
  IO::log << "\t Temperature = "
	  << temperature() / Phys_const::kelv << " K\n";

  IO::aux << "Pressure = " << pressure() / Phys_const::bar << " bar"
	  << "\t Temperature = " << temperature() / Phys_const::kelv << " K\n";
  
  // collisional frequency output
  //
  IO::log << IO::log_offset
	  << std::setw(10) << "Well"
	  << std::setw(20) << "Collision, 1/sec"
	  << "\n";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << IO::log_offset
	    << std::setw(10) << Model::well(w).name()
	    << std::setw(20) << well(w).collision_frequency() / Phys_const::herz
	    << "\n";

  if(evec_out.is_open()) {
    //
    evec_out << "Pressure = ";
    
    switch(pressure_unit) {
    case BAR:
      evec_out << pressure() / Phys_const::bar << " bar";
      break;
    case TORR:
      evec_out << pressure() / Phys_const::tor << " torr";
      break;
    case ATM:
      evec_out << pressure() / Phys_const::atm << " atm";
      break;
    }
    
    evec_out << "\t Temperature = " << temperature() / Phys_const::kelv << " K\n";
  }

  rate_data.clear();

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;

  Lapack::Vector      vtemp;
  Lapack::Matrix      mtemp;

  double proj;

  // total kinetic relaxation matrix dimension and index shifts for individual wells
  //
  std::vector<int> well_shift(Model::well_size());
  
  itemp = 0;
  
  for(int w = 0; w < Model::well_size(); itemp += well(w++).size())
    //
    well_shift[w] = itemp;

  const int global_size = itemp;

  IO::log << IO::log_offset << "global relaxation matrix dimension = " << global_size << "\n";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    if(!w || well(w).size() > itemp)
      //
      itemp = well(w).size();
  
  const int well_size_max = itemp;
    
  /********************************* SETTING GLOBAL MATRICES *********************************/

  // kinetic relaxation matrix
  //
  Lapack::SymmetricMatrix kin_mat(global_size); // kinetic relaxation matrix
  
  kin_mat = 0.;

  // bimolecular product vectors
  //
  Lapack::Matrix global_bim;
  
  if(Model::bimolecular_size()) {
    //
    global_bim.resize(global_size,  Model::bimolecular_size());
    
    global_bim = 0.;
  }

  // Boltzmann distributions
  //
  Lapack::Matrix global_pop(global_size, Model::well_size());
  
  global_pop = 0.;

  // well escape vectors
  //
  Lapack::Matrix global_escape;
  
  if(Model::escape_size()) {
    //
    global_escape.resize(global_size, Model::escape_size());
    
    global_escape = 0.;
  }

  {
    IO::Marker set_marker("setting global matrices", IO::Marker::ONE_LINE);

    // kin_mat initialization
    
    // nondiagonal isomerization contribution
    //
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      int w1 = Model::inner_connect(b).first;
      
      int w2 = Model::inner_connect(b).second;
      
      for(int i = 0; i < inner_barrier(b).size(); ++i)
	//
	kin_mat(i + well_shift[w1], i + well_shift[w2]) = - inner_barrier(b).state_number(i) / 2. / M_PI
	  //
	  / std::sqrt(well(w1).state_density(i) * well(w2).state_density(i));
    }

    // diagonal isomerization contribution
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      for(int i = 0; i < cum_stat_num[w].size(); ++i)
	//
	kin_mat(i + well_shift[w], i + well_shift[w]) = cum_stat_num[w][i] / 2. / M_PI
	  //
	  / well(w).state_density(i);
      
      if(Model::well(w).escape_size())
	//
	for(int i = 0; i < well(w).size(); ++i)
	  //
	  kin_mat(i + well_shift[w], i + well_shift[w]) += Model::well(w).escape_rate(energy_bin(i));
    }

    // collision relaxation contribution
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      for(int i = 0; i < well(w).size(); ++i) {
	//
	for(int j = i; j < well(w).size(); ++j)
	  //
	  kin_mat(i + well_shift[w], j + well_shift[w]) +=  well(w).collision_frequency() * well(w).kernel(i, j)
	    //
	    * well(w).boltzman_sqrt(i) / well(w).boltzman_sqrt(j);
      }
    }

    // radiational transitions contribution
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      if(well(w).radiation()) {

#pragma omp parallel for default(shared) schedule(dynamic)
	
	for(int i = 0; i < well(w).size(); ++i) {
	  //
	  for(int j = i; j < well(w).size(); ++j)
	    //
	    kin_mat(i + well_shift[w], j + well_shift[w]) +=  well(w).radiation_rate(i, j); 
	}
      }

    // bimolecular product vectors
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int w = Model::outer_connect(b).first;
      
      const int p = Model::outer_connect(b).second;

      for(int i = 0; i < outer_barrier(b).size(); ++i)
	//
	global_bim(i + well_shift[w], p) = outer_barrier(b).state_number(i) / 2. / M_PI
	  //
	  * thermal_factor(i) / well(w).boltzman_sqrt(i);
    }
  
    // thermal distributions
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      for(int i = 0; i < well(w).size(); ++i)
	//
	global_pop(i + well_shift[w], w) = well(w).boltzman_sqrt(i) / well(w).weight_sqrt();
    
    // well escape
    //
    for(int e = 0; e < Model::escape_size(); ++e) {
      //
      const int ew = Model::escape_channel(e).first;

      const int ei = Model::escape_channel(e).second;
      
      for(int i = 0; i < well(ew).size(); ++i)
	//
	global_escape(i + well_shift[ew], e) = Model::well(ew).escape_rate(energy_bin(i), ei) * well(ew).boltzman_sqrt(i);
      //
    }//
    //
  }// global matrices

  /******************** DIAGONALIZING THE GLOBAL KINETIC RELAXATION MATRIX ********************/

  Lapack::Vector eigenval;
  
  Lapack::Matrix eigen_global(global_size);

  {
    IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);

#ifdef WITH_MPACK
    
    switch(float_type) {
      //
    case DOUBLE:
	//
	eigenval = kin_mat.eigenvalues(&eigen_global);

	break;

      case DD:
	//
	eigenval = Mpack::dd_eigenvalues(kin_mat, &eigen_global);

	break;
	
      case QD:
	//
	eigenval = Mpack::qd_eigenvalues(kin_mat, &eigen_global);

	break;

      default:
	//
	std::cerr << funame << "unknown case\n";

	throw Error::Logic();
    }

#else

    eigenval = kin_mat.eigenvalues(&eigen_global);

#endif

    eigen_global = eigen_global.transpose();
  }

  const double min_relax_eval = eigenval[Model::well_size()];
  
  const double max_relax_eval = eigenval.back();

  //IO::log << IO::log_offset << "minimal relaxation eigenvalue / collision frequency = "
  //<< min_relax_eval / collision_frequency() << "\n";
  //IO::log << IO::log_offset << "maximal relaxation eigenvalue / collision frequency = "
  //<< max_relax_eval / collision_frequency() << "\n";

  // microcanonical rate coefficients at highest energy
  IO::log << IO::log_offset 
	  << "microscopic rate coefficients (at reference energy) over collision frequency:\n";
  
  // inner barriers
  //
  if(Model::inner_barrier_size()) {
    //
    IO::log << IO::log_offset;
    
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      int w1 = Model::inner_connect(b).first;
      
      int w2 = Model::inner_connect(b).second;
      
      IO::log << std::setw(13) << Model::well(w1).name() + "<->" + Model::well(w2).name();
    }
    
    IO::log << "\n" << IO::log_offset;

    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      int w1 = Model::inner_connect(b).first;
      
      int w2 = Model::inner_connect(b).second;
      
      dtemp = well(w1).state_density(0) < well(w2).state_density(0) ?
					  //
					  well(w1).state_density(0) : well(w2).state_density(0);
      
      IO::log << std::setw(13) << inner_barrier(b).state_number(0) / 2. / M_PI
	//
	/ dtemp / well(w1).collision_frequency();
    }
    
    IO::log << "\n";
    //
  }// inner barriers
  
  // outer barriers
  //
  if(Model::outer_barrier_size()) {
    //
    IO::log << IO::log_offset;
    
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      int w = Model::outer_connect(b).first;
      
      int p = Model::outer_connect(b).second;
      
      IO::log << std::setw(13) << Model::well(w).name() + "->" + Model::bimolecular(p).name();
    }
    
    IO::log << "\n" << IO::log_offset;

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      int w = Model::outer_connect(b).first;
      
      IO::log << std::setw(13) << outer_barrier(b).state_number(0) / 2. / M_PI
	//
	/ well(w).state_density(0) / well(w).collision_frequency();
    }
    
    IO::log << "\n";
    //
  }// outer barriers

  // low eigenvalue method
  //
  if(eigenval[0] / min_relax_eval < min_chem_eval) {
    //
    IO::log << IO::log_offset << "some eigenvalues are too small: using low eigenvalue method\n";
    
    IO::Marker low_eval_marker("low eigenvalue method");
    
    Lapack::SymmetricMatrix k_11;
    Lapack::SymmetricMatrix k_33;
    Lapack::Matrix k_13;
    Lapack::Matrix l_21;

    low_eigenvalue_matrix(k_11, k_33, k_13, l_21);

    Lapack::Matrix chem_evec(Model::well_size());
    
#ifdef WITH_MPACK
  
    Lapack::Vector chem_eval = Mpack::dd_eigenvalues(k_11, &chem_evec);

#else
    
    Lapack::Vector chem_eval = k_11.eigenvalues(&chem_evec);

#endif
    
    // low-eigenvalue chemical subspace
    //
    itemp = 1;
    
    while(eigenval[itemp] / min_relax_eval < min_chem_eval) { ++itemp; }
    
    const int chem_size = itemp < Model::well_size() ? itemp : Model::well_size();

    if(chem_size < Model::well_size()) {
      //
      Lapack::Matrix pop_chem(Model::well_size(), chem_size);
      
      for(int l = 0; l < chem_size; ++l)
	//
	pop_chem.column(l) = chem_evec.column(l);

      // partitioning wells
      //
      Partition low_well_partition;
      
      Group bimolecular_group;
      
      well_partition_method(pop_chem, low_well_partition, bimolecular_group, std::vector<double>());

      IO::log << IO::log_offset << "low eigenvalue bound species:";
      
      for(int g = 0; g < low_well_partition.size(); ++g) {
	//
	IO::log << " ";
	
	for(Git w = low_well_partition[g].begin(); w != low_well_partition[g].end(); ++w) {
	  //
	  if(w != low_well_partition[g].begin())
	    //
	    IO::log << "+";
	  
	  IO::log << Model::well(*w).name();
	}
      }
      
      IO::log << "\n" << IO::log_offset << "low eigenvalue bimolecular group:";
      
      for(Git w = bimolecular_group.begin(); w != bimolecular_group.end(); ++w)
	//
	IO::log << " " << Model::well(*w).name();
      
      IO::log << "\n";

      Lapack::SymmetricMatrix low_km(chem_size);
      
      low_km = 0.;
      
      for(int i = 0; i < chem_size; ++i)\
	//
	for(int j = i + 1; j < chem_size; ++j) {
	  //
	  dtemp = 0.;
	  
	  for(Git ig = low_well_partition[i].begin(); ig != low_well_partition[i].end(); ++ig)
	    //
	    for(Git jg = low_well_partition[j].begin(); jg != low_well_partition[j].end(); ++jg)
	      //
	      dtemp += k_11(*ig, *jg) * well(*ig).weight_sqrt() * well(*jg).weight_sqrt();
	  
	  low_km(i, j) = dtemp;
	}

      // diagonal terms
      //
      for(int i = 0; i < chem_size; ++i) {
	//
	// cross terms
	//
	dtemp = 0.;
	
	for(int j = 0; j < chem_size; ++j)
	  //
	  if(j != i)
	    //
	    dtemp -= low_km(i, j);
	
	low_km(i, i) = dtemp;

	for(Git ig = low_well_partition[i].begin(); ig != low_well_partition[i].end(); ++ig) {
	  //
	  dtemp = 0.;
	  
	  // bimolecular channel contribution
	  //
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    dtemp += k_13(*ig, p);
	  
	  // bimolecular group contribution
	  //
	  for(Git jg = bimolecular_group.begin(); jg != bimolecular_group.end(); ++jg)
	    //
	    dtemp -= k_11(*ig, *jg) * well(*jg).weight_sqrt();
	  
	  low_km(i, i) += dtemp * well(*ig).weight_sqrt();
	}
      }

      std::vector<double>      weight = low_well_partition.weight();
      
      Lapack::Matrix            basis = low_well_partition.basis();

      for(int i = 0; i < chem_size; ++i)
	//
	for(int j = i; j < chem_size; ++j)
	  //
	  low_km(i, j) /= std::sqrt(weight[i] * weight[j]);

      Lapack::Matrix low_evec(chem_size);
      
      Lapack::Vector low_eval = low_km.eigenvalues(&low_evec);
      
      IO::log << IO::log_offset << "low eigenvalues over minimal relaxation eigenvalue:\n"
	//
	      << IO::log_offset
	//
	      << std::setw(16) << "projection"
	//
	      << std::setw(16) << "diagonalization"
	      << "\n";

      for(int l = 0; l < chem_size; ++l) {
	//
	IO::log << IO::log_offset
	  //
		<< std::setw(16) << low_eval[l] / min_relax_eval
	  //
		<< std::setw(16) << chem_eval[l] / min_relax_eval
	  //
		<< "\n";
	/*
	  if(reduction_method == PROJECTION) {
	  chem_eval[l] = low_eval[l];
	  for(int w = 0; w < Model::well_size(); ++w)
	  chem_evec(w, l) = low_evec.column(l) * basis.row(w);
	  }
	*/
      }
    }

    // global eigenvector relaxational part
    //
    l_21 = l_21 * chem_evec;

    Lapack::Vector rel_proj(Model::well_size());
    
    for(int l = 0; l < Model::well_size(); ++l)
      //
      rel_proj[l] = vdot(l_21.column(l));

    // chemical eigenvector renormalization
    //
    for(int l = 0; l < Model::well_size(); ++l) {
      //
      dtemp = std::sqrt(1. + rel_proj[l]);
      
      chem_evec.column(l) /= dtemp;
      
      l_21.column(l) /= dtemp;
    }

    // direct-digonalization versus low-eigenvalue output
    //
    IO::log << IO::log_offset << "direct-diagonalization(DD)-versus-low-eigenvalue(LE) eigenvalues\n"
      //
	    << IO::log_offset << std::setw(5) << "L"
      //
	    << std::setw(13) << "DD eval"
      //
	    << std::setw(13) << "LE eval"
      //
	    << std::setw(13) << "LE proj\n";

    for(int l = 0; l < Model::well_size(); ++l) {
      //
      IO::log << IO::log_offset <<  std::setw(5) << l
	//
	      << std::setw(13) << eigenval[l] / min_relax_eval
	//
	      << std::setw(13) << chem_eval[l] / min_relax_eval
	//
	      << std::setw(13) << rel_proj[l];
      
      if(eigenval[l] / min_relax_eval >= min_chem_eval)
	//
	IO::log << std::setw(3) << "*";
      
      IO::log << "\n";
    }

    // eigenvalue and eigenvector substitution
    //
    for(int l = 0; l < chem_size; ++l) {
      //
      // eigenvalue
      //
      eigenval[l] = chem_eval[l];

      // global eigenvector
      //
      itemp = 0;
      
      for(int w = 0; w < Model::well_size(); itemp += well(w++).crm_size()) {
	//
	for(int i = 0; i < well(w).size(); ++i) {
	  //
	  dtemp = well(w).boltzman_sqrt(i);
	  //
	  eigen_global(l, i + well_shift[w]) = chem_evec(w, l) * dtemp / well(w).weight_sqrt()
	    //
	    - vdot(&l_21(itemp, l), well(w).crm_row(i), well(w).crm_size(), 1, well(w).size()) / dtemp;
	  //
	}//
	//
      }//
      //
    }//
    //
  }// low eigenvalue method

  Lapack::Matrix eigen_well(global_size, Model::well_size());
  
  for(int l = 0; l < global_size; ++l)
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      eigen_well(l, w) = vlength(&eigen_global(l, well_shift[w]), well(w).size(), global_size);
  
  // projection of the  eigenvectors onto the thermal subspace
  //
  Lapack::Matrix eigen_pop = eigen_global * global_pop;

  // eigenvector to bimolecular vector projection
  //
  Lapack::Matrix eigen_bim;
  
  if(Model::bimolecular_size())
    //
    eigen_bim = eigen_global * global_bim;

  // eigenvector to well escape projection
  //
  Lapack::Matrix eigen_escape;
  
  if(Model::escape_size())
    //
    eigen_escape = eigen_global * global_escape;

  /********************************* TIME EVOLUTION ********************************/

  if(Model::time_evolution) {
    //
    const int react = Model::time_evolution->reactant();

    // output header
    //
    Model::time_evolution->out << "Pressure = ";
    
    switch(pressure_unit) {
      //
    case BAR:
      //
      Model::time_evolution->out << pressure() / Phys_const::bar << " bar";
      
      break;
      
    case TORR:
      //
      Model::time_evolution->out << pressure() / Phys_const::tor << " torr";
      
      break;
      
    case ATM:
      //
      Model::time_evolution->out << pressure() / Phys_const::atm << " atm";
      
      break;
    }
    
    Model::time_evolution->out << "\t Temperature = " << temperature() / Phys_const::kelv << " K\n\n";

    Model::time_evolution->out << std::setw(13) << "time, sec";
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      Model::time_evolution->out << std::setw(13) << Model::well(w).name();
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      Model::time_evolution->out << std::setw(13) << Model::bimolecular(p).name();

    Model::time_evolution->out << "\n";
    
    // bimolecular reactants
    //
    if(Model::bimolecular_size() && Model::time_evolution->excess_reactant_concentration() > 0.) {
      //
      if(bimolecular(react).weight() < 0.) {
	//
	std::cerr << funame << "time evolution: reactants partition function not defined\n";
	
	throw Error::Init();
      }
      
      // normaziation factor
      //
      const double nfac = Model::time_evolution->excess_reactant_concentration() * energy_step()
	/ bimolecular(react).weight();

      double time_val = Model::time_evolution->start();
      for(int t = 0;  t < Model::time_evolution->size(); ++t, time_val *= Model::time_evolution->step()) {
	std::vector<double> well_pop(Model::well_size());
	std::vector<double> bim_pop(Model::bimolecular_size());

	for(int l = 0; l < global_size; ++l) {
	  dtemp = eigenval[l] * time_val;
	  if(dtemp > 50.)
	    dtemp = 1. / eigenval[l];
	  else
	    dtemp = (1. - std::exp(-dtemp)) /eigenval[l];
	
	  for(int w = 0; w < Model::well_size(); ++w)
	    well_pop[w] += eigen_bim(l, react) * eigen_pop(l, w) * dtemp;

	  dtemp = (time_val - dtemp) / eigenval[l]; 
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    bim_pop[p] += eigen_bim(l, react) * eigen_bim(l, p) * dtemp;
	}

	// normalization
	for(int w = 0; w < Model::well_size(); ++w)
	  well_pop[w] *= well(w).weight_sqrt() * nfac;
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  bim_pop[p] *= nfac;
	  
	// output
	Model::time_evolution->out << std::setw(13) << time_val * Phys_const::herz;
	for(int w = 0; w < Model::well_size(); ++w)
	  Model::time_evolution->out << std::setw(13) << well_pop[w];

	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  Model::time_evolution->out << std::setw(13) << bim_pop[p];
	
	Model::time_evolution->out << "\n";
	//
      }// time cycle
      //
    }// bimolecular reactants
    //
    // bound reactant
    //
    else {
      const double bfac = std::exp(energy_step() / Model::time_evolution->temperature() - energy_step() / temperature());
    
      // initial distribution
      Lapack::Vector init_dist(well(react).size());
      dtemp = 1.;
      double norm_fac = 0.;
      for(int e = 0; e < well(react).size(); ++e, dtemp *= bfac) {
	init_dist[e] = well(react).boltzman_sqrt(e) * dtemp;
	norm_fac    += well(react).boltzman(e)      * dtemp;
      }
      
      for(Lapack::Vector::iterator i = init_dist.begin(); i != init_dist.end(); ++i)
	*i /= norm_fac;

      std::vector<double> init_coef(global_size);
      for(int l = 0; l < global_size; ++l)
	init_coef[l] = parallel_vdot(init_dist, &eigen_global(l, well_shift[react]), well(react).size(), 1, global_size);

      double time_val = Model::time_evolution->start();
      for(int t = 0;  t < Model::time_evolution->size(); ++t, time_val *= Model::time_evolution->step()) {
	std::vector<double> well_pop(Model::well_size());
	std::vector<double> bim_pop(Model::bimolecular_size());

	for(int l = 0; l < global_size; ++l) {
	  dtemp = eigenval[l] * time_val;
	  if(dtemp > 100.)
	    dtemp = 0.;
	  else
	    dtemp = std::exp(-dtemp);
	
	  for(int w = 0; w < Model::well_size(); ++w)
	    well_pop[w] += init_coef[l] * eigen_pop(l, w) * dtemp;

	  dtemp = (1. - dtemp) / eigenval[l]; 
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    bim_pop[p] += init_coef[l] * eigen_bim(l, p) * dtemp;
	}

	// normalization
	//
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  well_pop[w] *= well(w).weight_sqrt();
	
	// output
	//
	Model::time_evolution->out << std::setw(13) << time_val * Phys_const::herz;
	
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  Model::time_evolution->out << std::setw(13) << well_pop[w];

	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  Model::time_evolution->out << std::setw(13) << bim_pop[p];
	
	Model::time_evolution->out << "\n";
	
      }// time cycle
      //
    }// bound reactant

    Model::time_evolution->out << "\n";

  }// time evolution
  
  // eigenvector distributions at hot energies
  //
  Lapack::Matrix eigen_hot;
  
  if(hot_index.size()) {
    //
    eigen_hot.resize(global_size, hot_index.size());
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
      //
      for(int l = 0; l < global_size; ++l)
	//
	eigen_hot(l, count) = eigen_global(l, well_shift[hit->first] + hit->second)
	  //
	  / well(hit->first).boltzman_sqrt(hit->second);
    }
  }

  std::vector<double> relaxation_projection(Model::well_size());
  
  for(int l = 0; l < Model::well_size(); ++l)
    //
    relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));

 
  /************************************* EIGENVECTOR OUTPUT *****************************************/

  IO::log << std::setprecision(3)
    //
	  << IO::log_offset << "eigenvector populations normalized:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(10) << "*R"
    //
	  << std::setw(10) << "*P";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  IO::log << "\n";

  // maximal population
  //
  for(int l = 0; l < Model::well_size(); ++l) {
    //
    double pos_pop = 0.;
    
    double neg_pop = 0.;
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well(w).weight_sqrt();
      
      if(dtemp > 0.) {
	//
	pos_pop += dtemp;
      }
      if(dtemp < 0.)
	//
	neg_pop += dtemp;
    }
    
    double max_pop = pos_pop > -neg_pop ? pos_pop : neg_pop;
    
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(10) << eigenval[l] / min_relax_eval
      //
	    << std::setw(10) << relaxation_projection[l];
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well(w).weight_sqrt() / max_pop;
      
      IO::log << std::setw(10);
      
      if(dtemp > .01 || dtemp < -.01) {
	//
	IO::log << dtemp;
      }
      else
	//
	IO::log << 0;
    }
    
    IO::log << "\n";
  }

  IO::log << IO::log_offset
    //
	  << "*R - eigenvalue over the relaxation limit\n"
    //
	  << IO::log_offset
    //
	  << "*P - eigenvector projection squared on the relaxation subspace (= 1 -F_ne)\n";
  
  IO::log << IO::log_offset << "eigenvector projections:\n"
    //
	  << IO::log_offset << std::setw(5)  << "L" << std::setw(10) << "*Q" << std::setw(10) << "*P";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << Model::well(w).name();
  
  IO::log << "\n";

  for(int l = 0; l < Model::well_size(); ++l) {
    //
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(10) << eigenval[l] / well(0).collision_frequency()
      //
	    << std::setw(10) << relaxation_projection[l];
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(10) << eigen_pop(l,w);
    
    IO::log << "\n";
  }
  
  IO::log << IO::log_offset
    //
	  << std::setw(5) << "*Z"
    //
	  << std::setw(10) << "---"
    //
	  << std::setw(10) << "---";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(10) << well(w).weight_sqrt();
  
  IO::log << "\n";

  IO::log << IO::log_offset
    //
	  << "*Q - eigenvalue over the collision frequency in first well\n"
    //
	  << IO::log_offset
    //
	  << "*P - eigenvector projection squared on the relaxation subspace (= 1 - F_ne)\n"
    //
	  << IO::log_offset
    //
	  << "*Z - well partition function square root\n"
    //
	  << std::setprecision(6);
  
  // eigenvalues output
  //
  if(eval_out.is_open()) {
    //
    eval_out << std::setw(13) << temperature() / Phys_const::kelv
      //
	     << std::setw(13);
    
    switch(pressure_unit) {
      //
    case BAR:
      //
      eval_out << pressure() / Phys_const::bar;
      break;
      
    case TORR:
      //
      eval_out << pressure() / Phys_const::tor;
      break;
      
    case ATM:
      //
      eval_out << pressure() / Phys_const::atm;
      break;
    }
    
    eval_out << std::setw(13) << well(0).collision_frequency() / Phys_const::herz
      //
	     << std::setw(13) << min_relax_eval / well(0).collision_frequency();
    
    int eval_max = Model::well_size() + evec_out_num;
    
    for(int l = 0; l < eval_max; ++l)
      //
      eval_out << std::setw(13) << eigenval[l] / well(0).collision_frequency()
	//
	       << std::setw(13) <<  1. - vdot(eigen_pop.row(l));
    
    eval_out << "\n";
  }

  // eigenvector output
  //
  if(evec_out.is_open()) {
    //
    evec_out << "EIGENVECTORS:\n";
    
    int evec_max = Model::well_size() + evec_out_num;
    
    for(int l = 0; l < evec_max; ++l) {
      //
      evec_out << "l = " << l << "\n"
	//
	       << "eigenvalue / collision frequency = "
	//
	       << eigenval[l] / well(0).collision_frequency()
	//
	       << "\n";

      evec_out << std::setw(13) << "well length";
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << eigen_well(l, w);
      
      evec_out << "\n";

      evec_out << std::setw(13) << "E, kcal/mol";
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << Model::well(w).name();
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << Model::well(w).name();
      
      evec_out << "\n";

      for(int i = 0; i < well_size_max; ++i) {
	//
	evec_out << std::setw(13) << energy_bin(i) / Phys_const::kcal;
	
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  if(i < well(w).size()) {
	    //
	    // micropopulational distribution
	    //
	    dtemp = eigen_global(l, well_shift[w] + i) * well(w).boltzman_sqrt(i);
	    
	    evec_out << std::setw(13) << dtemp;
	  }
	  else
	    //
	    evec_out << std::setw(13) << 0;
	
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  if(i < well(w).size()) {
	    //
	    // ratio to the thermal distribution
	    //
	    dtemp = eigen_global(l, well_shift[w] + i) / well(w).boltzman_sqrt(i)
	      //
	      * well(w).weight_sqrt();
	    
	    evec_out << std::setw(13) << dtemp;
	  }
	  else
	    //
	    evec_out << std::setw(13) << 0;
	
	evec_out << "\n";
      }
    }
    
    evec_out << "\n";
  }

  /********************************** CHEMICAL SUBSPACE DIMENSION ****************************************/

  //
  // predefined chemical subspace dimension
  //
  if(_default_chem_size >= 0) {
    //
    itemp = _default_chem_size;
  }
  //
  // default partitioning scheme
  //
  else if(default_partition.size()) {
    //
    itemp = default_partition.size();
  }
  //
  // absolute eigenvalue threshold
  //
  else if(chemical_threshold > 1.) {
    //
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(min_relax_eval / eigenval[itemp]  < chemical_threshold)
	//
	break;
  }
  //
  // relaxation projection threshold
  //
  else if(chemical_threshold < 1. && chemical_threshold > 0.) {
    //
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(relaxation_projection[itemp] > chemical_threshold)
	//
	break;
  }
  //
  // relative eigenvalue threshold
  //
  else if(chemical_threshold < -1. ) {
    //
    for(itemp = Model::well_size(); itemp > 0; --itemp)
      //
      if(eigenval[itemp] / eigenval[itemp - 1] > - chemical_threshold)
	//
	break;
  }
  //
  else
    //
    itemp = Model::well_size();
	
  const int chem_size = itemp;

  if(chem_size != Model::well_size())
    //
    IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n";

  /***** PARTITIONING THE GLOBAL PHASE SPACE INTO THE CHEMICAL AND COLLISIONAL SUBSPACES *****/

  // collisional relaxation eigenvalues and eigenvectors
  //
  const int relax_size = global_size - chem_size;
  
  Lapack::Vector relax_lave(relax_size);
  
  for(int r = 0; r < relax_size; ++r) {
    //
    itemp = r + chem_size;

    relax_lave[r] = 1. / eigenval[itemp];
  }

  // kinetic matrix modified

#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)
	
  for(int i = 0; i < global_size; ++i) {
    //
    for(int j = i; j < global_size; ++j) {
      //
      dtemp = 0.;
      
      for(int l = 0; l < chem_size; ++l)
	//
	dtemp += eigen_global(l, i) * eigen_global(l, j);
      
      kin_mat(i, j) += dtemp * well(0).collision_frequency();
    }
  }

  Lapack::Matrix proj_bim = global_bim.copy();
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    for(int l = 0; l < chem_size; ++l)
      //
      parallel_orthogonalize(&proj_bim(0, p), &eigen_global(l, 0), global_size, 1, global_size);

  Lapack::Matrix inv_proj_bim;
  
  if(Model::bimolecular_size())
    //
    inv_proj_bim = Lapack::Cholesky(kin_mat).invert(proj_bim);

  Lapack::Matrix proj_pop = global_pop.copy();
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    for(int l = 0; l < chem_size; ++l)
      //
      parallel_orthogonalize(&proj_pop(0, w), &eigen_global(l, 0), global_size, 1, global_size);
  
  // kappa matrix
  //
  if(Model::bimolecular_size()) {
    //
    Lapack::Matrix kappa = proj_pop.transpose() * inv_proj_bim;

    IO::log << std::setprecision(2)
      //
	    << IO::log_offset << "isomers-to-bimolecular equilibrium coefficients (kappa matrix):\n"
      //
	    << IO::log_offset << std::setw(5) << "W\\P";
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(9) << Model::bimolecular(p).name();
    
    IO::log << "\n";
   
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(w).name();
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	//dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w), relax_lave, relax_size) 
	// / well(w).weight_sqrt();
	//
	dtemp = kappa(w, p) / well(w).weight_sqrt();
	
	IO::log << std::setw(9);
	
	if(dtemp < 0.05 && dtemp > -0.05) {
	  //
	  IO::log << "0";
	}
	else
	  //
	  IO::log << dtemp;
      }
      
      IO::log << "\n";
    }

    IO::log << std::setprecision(6);
  }

  // bimolecular-to-bimolecular rate coefficients
  //
  if(Model::bimolecular_size()) {
    //
    Lapack::SymmetricMatrix bb_rate = Lapack::SymmetricMatrix(proj_bim.transpose() * inv_proj_bim);
    
    //  Lapack::SymmetricMatrix bb_rate(Model::bimolecular_size());
    //  for(int i = 0; i < Model::bimolecular_size(); ++i)
    //    for(int j = i; j < Model::bimolecular_size(); ++j)
    //	bb_rate(i, j) = triple_product(&eigen_bim(chem_size, i), &eigen_bim(chem_size, j), 
    //				       relax_lave, relax_size);

    for(int i = 0; i < Model::bimolecular_size(); ++i)
      //
      if(bimolecular(i).weight() > 0.)
	//
	for(int j = 0; j < Model::bimolecular_size(); ++j) {
	  //
	  // bimolecular reactant loss
	  //
	  if(i == j) {
	    //
	    dtemp = 0.;
	    
	    for(int b = 0; b < Model::outer_barrier_size(); ++b)
	      //
	      if(Model::outer_connect(b).second == i)
		//
		//dtemp += temperature() / 2. / M_PI * outer_barrier(b).weight();
		for(int e = 0; e < outer_barrier(b).size(); ++e)
		  //
		  dtemp += outer_barrier(b).state_number(e) * thermal_factor(e);
	    
	    dtemp /= 2. * M_PI;
	    
	    dtemp -= bb_rate(i, i);
	  }
	  // crossrate
	  //
	  else
	    //
	    dtemp = bb_rate(i, j);

	  dtemp *= energy_step() / bimolecular(i).weight() / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 
	  //arr_out;
	}
    
    // bimolecular-to-escape rate coefficients
    //
    if(Model::escape_size()) {
      /*
      Lapack::Matrix proj_escape = global_escape.copy();
      for(int count = 0; count < Model::escape_size(); ++count)
	for(int l = 0; l < chem_size; ++l)
	  parallel_orthogonalize(&proj_escape(0, count), &eigen_global(l, 0), global_size, 1, global_size);
    
      Lapack::Matrix escape_bim = proj_escape.transpose() * inv_proj_bim;
      */
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(!Model::bimolecular(p).dummy())
	  //
	  for(int e = 0; e < Model::escape_size(); ++e) {
	    //
	    dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_escape(chem_size, e),
				   //
				   relax_lave, relax_size) * energy_step() / bimolecular(p).weight() / bru;
	    
	    //dtemp = escape_bim(e, p) * energy_step() / bimolecular(p).weight() / bru;
	    
	    rate_data[std::make_pair(Model::well_size() + p, Model::well_size() +
				     //
				     Model::bimolecular_size() + e)] = dtemp; 
	  }
    }
  }

  // product energy distributions
  //
  if(ped_out.is_open()) {
    //
    switch(pressure_unit) {
    case BAR:
      ped_out << "pressure[bar]        = " << pressure() / Phys_const::bar << "\n";
      break;
    case TORR:
      ped_out << "pressure[torr]       = " << pressure() / Phys_const::tor << "\n";
      break;
    case ATM:
      ped_out << "pressure[atm]        = " << pressure() / Phys_const::atm << "\n";
      break;
    }
    
    ped_out << "temperature[K]       = " << temperature() / Phys_const::kelv << "\n"
	    << "energy step[1/cm]    = " << energy_step() / Phys_const::incm << "\n"
	    << "maximum energy[1/cm] = " << energy_reference() / Phys_const::incm << "\n\n";

    int ener_index_max;

    // bimolecular-to-bimolecular product energy distributions
    //
    if(Model::bimolecular_size()) {
      //
      ped_out << "Bimolecular-to-bimolecular product energy distributions:\n";

      // dimensions
      //
      itemp = 0;
      
      for(int ped = 0; ped < ped_pair.size(); ++ ped)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(Model::outer_connect(b).second == ped_pair[ped].second && outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

      ener_index_max = itemp;
      
      mtemp.resize(itemp, ped_pair.size());

      // distribution
      //
      for(int e = 0; e < ener_index_max; ++e)
	//
	for(int ped = 0; ped < ped_pair.size(); ++ ped) {
	  //
	  dtemp = 0.;
	  
	  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	    //
	    const int w = Model::outer_connect(b).first;
	    
	    const int p = Model::outer_connect(b).second;
	    
	    if(p == ped_pair[ped].second && e < outer_barrier(b).size())
	      //
	      dtemp += global_bim(e + well_shift[w], p) *
		//
		triple_product(&eigen_bim(chem_size, ped_pair[ped].first), relax_lave,
			       //
			       &eigen_global(chem_size, e + well_shift[w]), relax_size);
	  }
	  
	  mtemp(e, ped) = dtemp;
	}
    
      // normalization
      //
      for(int i = 0; i < mtemp.size2(); ++i)
	//
	mtemp.column(i) /= max(mtemp.column(i));

      // output
      //
      ped_out << std::setw(13) << "E, kcal/mol";
      
      for(int ped = 0; ped < ped_pair.size(); ++ ped)
	//
	ped_out << std::setw(13) << Model::bimolecular(ped_pair[ped].first).name() + "->"
	  //
	  + Model::bimolecular(ped_pair[ped].second).name();
      
      ped_out << "\n";

      for(int e = 0; e < ener_index_max; ++e) {
	//
	ped_out << std::setw(13) << energy_bin(e) / Phys_const::kcal;

	for(int i = 0; i < mtemp.size2(); ++i)
	  //
	  ped_out << std::setw(13) << mtemp(e, i);
	
	ped_out << "\n";
      }
      
      ped_out << "\n";

      // escape product energy distributions
      //
      if(Model::escape_size()) {
	//
	// Hot energies-to-escape channels distribution
	//
	if(hot_index.size()) {
	  //
	  ped_out << "Hot-to-escape product energy distributions:\n";

	  for(int e = 0; e < Model::escape_size(); ++e) {
	    //
	    const int ew = Model::escape_channel(e).first;

	    const int ei = Model::escape_channel(e).second;

	    vtemp.resize(well(ew).size());
	    
	    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit) {
	      //
	      const int hw = hit->first;
	      
	      const int he = hit->second;
		
	      ped_out << "Hot[" << Model::well(hw).name() << ", E = "
		//
		      << energy_bin(he) / Phys_const::kcal << " kcal/mol] ---> "
		//
		      << Model::escape_name(e) << "\n\n";

	      double nfac;
		
	      for(int i = 0; i < well(ew).size(); ++i) {
		//
		dtemp = triple_product(&eigen_global(chem_size, he + well_shift[hw]), relax_lave,
				       //
				       &eigen_global(chem_size, i  + well_shift[ew]), relax_size)
		  //
		  * Model::well(ew).escape_rate(energy_bin(i), ei) * well(ew).boltzman_sqrt(i);
		
		vtemp[i] = dtemp;
		  
		if(!i || std::fabs(dtemp) > std::fabs(nfac))
		  //
		  nfac = dtemp;
	      }
		
	      vtemp /= nfac;
		
	      // output
	      //
	      ped_out << std::setw(13) << "E, kcal/mol" << std::setw(13) << "PED, a.u." << "\n";

	      for(int i = 0; i < well(ew).size(); ++i) {
		//
		dtemp = energy_bin(i) / Phys_const::kcal;
		  
		ped_out << std::setw(13) << dtemp << std::setw(13) << vtemp[i] << "\n";
	      }

	      ped_out << "\n";
	    }
	  }
	}
	
	// bimolecular-to-escape product energy distributions
	//
	ped_out << "Bimolecular-to-escape product energy distributions:\n";

	// dimensions
	//
	itemp = 0;
	
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  const int w  = Model::escape_channel(e).first;
	  
	  if(well(w).size() > itemp)
	    //
	    itemp = well(w).size();
	}
	
	ener_index_max = itemp;

	itemp = 0;
	
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  if(!Model::bimolecular(p).dummy())
	    //
	    ++itemp;
	
	itemp *= Model::escape_size();

	mtemp.resize(ener_index_max, itemp);

	for(int e = 0; e < ener_index_max; ++e) {
	  //
	  itemp = 0;
	  
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    if(!Model::bimolecular(p).dummy())
	      //
	      for(int count = 0; count < Model::escape_size(); ++count, ++itemp) {
		//
		const int ew = Model::escape_channel(count).first;

		if(e < well(ew).size()) {
		  //
		  mtemp(e, itemp) =  global_escape(e + well_shift[ew], count) *
		    //
		    triple_product(&eigen_bim(chem_size, p), relax_lave,
				   //
				   &eigen_global(chem_size, e + well_shift[ew]), relax_size);
		}
		else
		  //
		  mtemp(e, itemp) = 0.;
	      }
	}

	// normalization
	//
	for(int i = 0; i < mtemp.size2(); ++i)
	  //
	  mtemp.column(i) /= max(mtemp.column(i));

	// output
	//
	ped_out << std::setw(13) << "E, kcal/mol";
	
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  if(!Model::bimolecular(p).dummy())
	    //
	    for(int e = 0; e < Model::escape_size(); ++e)
	      //
	      ped_out << std::setw(13) << Model::bimolecular(p).name() + "->" + Model::escape_name(e);
	    
	ped_out << "\n";

	for(int e = 0; e < ener_index_max; ++e) {
	  //
	  ped_out << std::setw(13) << energy_bin(e) / Phys_const::kcal;

	  for(int i = 0; i < mtemp.size2(); ++i)
	    //
	    ped_out << std::setw(13) << mtemp(e, i);
	  
	  ped_out << "\n";
	}
	ped_out << "\n";
	//
      }// escape output

      // hot product energy distributions
      //
      if(hot_index.size()) {
	//
	ped_out << "Hot product energy distributions:\n\n";

	// dimension
	//
	itemp = 0;
	
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();
      
	ener_index_max = itemp;
	//
	mtemp.resize(itemp, Model::bimolecular_size());

	// hot energy cycle
	//
	int count = 0;
	
	for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
	  //
	  ped_out << "Initial well: "<< Model::well(hit->first).name()
	      //
		  << "   Initial energy[kcal/mol] = " << energy_bin(hit->second) / Phys_const::kcal
	    //
		  << "\n";

	    // distribution
	    //
	    for(int e = 0; e < ener_index_max; ++e)
	      //
	      for(int p = 0; p < Model::bimolecular_size(); ++p) {
		//
		dtemp = 0.;
		
		for(int b = 0; b < Model::outer_barrier_size(); ++b) {
		  //
		  const int w = Model::outer_connect(b).first;
		  
		  if(p == Model::outer_connect(b).second && e < outer_barrier(b).size())
		    //
		    dtemp += global_bim(e + well_shift[w], p) *
		      //
		      triple_product(&eigen_global(chem_size, e + well_shift[w]), relax_lave,
				     //
				     &eigen_hot(chem_size, count), relax_size);
		}
		
		mtemp(e, p) = dtemp;
	      }
	  
	    // normalization
	    //
	    for(int i = 0; i < mtemp.size2(); ++i)
	      //
	      mtemp.column(i) /= max(mtemp.column(i));

	    //output
	    //
	    ped_out << std::setw(13) << "E, kcal/mol";
	    
	    for(int p = 0; p < Model::bimolecular_size(); ++p)
	      //
	      ped_out << std::setw(13) << Model::bimolecular(p).name();
	    
	    ped_out << "\n";
	
	    for(int e = 0; e < ener_index_max; ++e) {
	      //
	      ped_out << std::setw(13) << energy_bin(e) / Phys_const::kcal;

	      for(int i = 0; i < mtemp.size2(); ++i)
		//
		ped_out << std::setw(13) << mtemp(e, i);
	      
	      ped_out << "\n";
	    }
	    ped_out << "\n";
	    //
	  }// hot energy cycle
	//
      }// hot distributions
      //
    }//
    //
  }// ped output

  /*************************************** BOUND SPECIES *******************************************/

  std::vector<int>  group_index;
  Lapack::Matrix m_direct;
  std::vector<double>      weight;
  std::vector<double> real_weight;

  if(chem_size) {
    //
    // projection of the chemical eigenvectors onto the thermal subspace
    //
    Lapack::Matrix pop_chem(Model::well_size(), chem_size);
    
    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

#ifdef DEBUG

    if(chem_size) {
      //
      IO::log << IO::log_offset << "orthogonality check starts\n";
      
      IO::log_offset.increase();

      proj = 0.;
      
      for(int l = 0; l < chem_size; ++l) {
	//
	// normalize
	//
	//normalize(&pop_chem(0, l), Model::well_size());
	//
	for(int m = 0; m < l; ++m) {
	  //
	  dtemp = vdot(pop_chem.column(l), pop_chem.column(m));
	  
	  dtemp = dtemp >= 0. ? dtemp : -dtemp;
	  
	  proj = dtemp > proj ? dtemp : proj;
	}
      }

      IO::log << IO::log_offset << "maximal scalar product of different chemical eigenvectors = "
	//
	      << proj << "\n";

      for(int l = 0; l < chem_size; ++l) {
	//
	dtemp = vdot(pop_chem.column(l));
	
	if(!l || dtemp < proj)
	  //
	  proj = dtemp;
      }

      IO::log << IO::log_offset << "minimal chemical eigenvector square = "
	//
	      << proj << "\n";

      for(int l = 0; l < chem_size; ++l) {
	//
	dtemp = vdot(pop_chem.column(l));
	
	if(!l || dtemp > proj)
	  //
	  proj = dtemp;
      }

      IO::log << IO::log_offset << "maximal chemical eigenvector square = "
	//
	      << proj << "\n";

      IO::log_offset.decrease();
      
      IO::log << IO::log_offset << "orthogonality check done\n";
    }

#endif

    // partitioning wells into equilibrated groups
    //
    if(default_partition.size()) {
      //
      // default reduction scheme
      //
      well_partition = default_partition;

      IO::log << IO::log_offset << "using default reduction scheme, projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    
      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }
    else if(chem_size == Model::well_size()) {
      //
      // no partitioning
      //
      well_partition.resize(Model::well_size());
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	well_partition[w].insert(w);

      IO::log << IO::log_offset << "projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    }
    else {
      //
      // well partitioning
      //
      Group bimolecular_group;
      
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group, std::vector<double>());

      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }

    group_index = well_partition.group_index();
    weight      = well_partition.weight();
    real_weight = well_partition.real_weight();

    // auxiliary output
    //
    IO::aux << "number of species = " << well_partition.size() << "\n";
    
    for(int p = 0; p < well_partition.size(); ++p) {
      //
      std::multimap<double, int> ww_map;
      
      for(Git git = well_partition[p].begin(); git != well_partition[p].end(); ++git)
	//
	ww_map.insert(std::make_pair(well(*git).real_weight(),*git));

      for(std::multimap<double, int>::const_reverse_iterator mit = ww_map.rbegin(); mit != ww_map.rend(); ++mit)
	//
	IO::aux << std::setw(15) << Model::well(mit->second).name();

      IO::aux << "\n";
      
      for(std::multimap<double, int>::const_reverse_iterator mit = ww_map.rbegin(); mit != ww_map.rend(); ++mit)
	//
	IO::aux << std::setw(15) << mit->first;
      
      IO::aux << "\n";
    }
    // output
    //
    if(chem_size != Model::well_size()) {
      //
      IO::log << IO::log_offset << "species:\n"
	//
	      << IO::log_offset << std::setw(2) << "#"  << std::setw(15)
	//
	      << "assigned name" << IO::first_offset
	//
	      << "group\n";
      
      for(Pit g = well_partition.begin(); g != well_partition.end(); ++g) {
	//
	itemp = g - well_partition.begin();
	
	IO::log << IO::log_offset << std::setw(2) << itemp << std::setw(15)
	  //
		<< Model::well(group_index[itemp]).name() << IO::first_offset;
	
	for(Git w = g->begin(); w != g->end(); ++w) {
	  //
	  if(w != g->begin())
	    //
	    IO::log << "+";
	  
	  IO::log << Model::well(*w).name();
	}
	
	IO::log << "\n";
      }
    }

    m_direct = pop_chem;
    
    Lapack::Matrix m_inverse = m_direct.invert();

    Lapack::Matrix one(m_direct * m_inverse);
    
    one.diagonal() -= 1.;
    
    double val_max = -1.;
    
    for(int i = 0; i < one.size1(); ++i)
      //
      for(int j = 0; j < one.size2(); ++j) {
	//
	dtemp = one(i, j);
	
	dtemp = dtemp < 0. ? -dtemp : dtemp;
	
	if(dtemp > epsilon && dtemp > val_max)
	  //
	  val_max = dtemp;
      }
    
    if(val_max > 0.)
      //
      IO::log << IO::log_offset << funame
	//
	      << "WARNING: matrix inversion error = " << val_max
	//
	      << " exceeds numerical accuracy = " << epsilon
	//
	      << "\n";
  
    // well-to-well rate coefficients
    //
    Lapack::Matrix ww_rate(chem_size);
    
    ww_rate = 0.;
    
    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j)
	//
	for(int l = 0; l < chem_size; ++l)
	  //
	  ww_rate(i, j) += m_direct(j, l) * m_inverse(l, i) * eigenval[l];

    // well-to-bimolecular rate coefficients
    //
    Lapack::Matrix wb_rate, bw_rate;
    
    if(Model::bimolecular_size()) {
      //
      wb_rate.resize(chem_size, Model::bimolecular_size());
      
      for(int w = 0; w < chem_size; ++w)
	//
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  wb_rate(w, p) = vdot(m_inverse.column(w), &eigen_bim(0, p));

      // bimolecular-to-well rate coefficients
      //
      bw_rate.resize(Model::bimolecular_size(), chem_size);
      
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int w = 0; w < chem_size; ++w)
	  //
	  bw_rate(p, w) = vdot(m_direct.row(w), &eigen_bim(0, p));
    }

    // output
    //
    IO::log << std::setprecision(2)
      //
	    << IO::log_offset << "Wa->Wb/Wb->Wa rate constants ratios:\n"
      //
	    << IO::log_offset << std::setw(5) << "Wb\\Wa";
    
    for(int i = 0; i < chem_size; ++i)
      //
      IO::log << std::setw(10) << Model::well(group_index[i]).name();
    
    IO::log << "\n";
    
    for(int j = 0; j < chem_size; ++j) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(group_index[j]).name();
      
      for(int i = 0; i < chem_size; ++i)
	//
	if(i != j) {
	  //
	  if(ww_rate(j, i) != 0.) {
	    //
	    IO::log << std::setw(10) << ww_rate(i, j) / ww_rate(j, i);
	  }
	  else
	    //
	    IO::log << std::setw(10) << "***";
	}
	else
	  //
	  IO::log << std::setw(10) << "1";
      
      IO::log << "\n";
    }
    
    if(Model::bimolecular_size()) {
      //
      IO::log << IO::log_offset << "W->P/P->W rate constants ratios:\n"
	//
	      << IO::log_offset << std::setw(5) << "P\\W";
      
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(10) << Model::well(group_index[w]).name();
      
      IO::log << "\n";
    
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).name();
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  if(bw_rate(p, w) != 0.) {
	    //
	    IO::log << std::setw(10) << wb_rate(w, p) / bw_rate(p, w);
	  }
	  else
	    //
	    IO::log << std::setw(10) << "***";
	
	IO::log << "\n";
      }
    }

    IO::log << std::setprecision(6);

    // well-to-well rate coefficients
    //
    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j) {
	//
	dtemp = ww_rate(i, j) * std::sqrt(weight[i] * weight[j]) * energy_step() / real_weight[i] / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp; 	
      }

    // well-to-escape rate coefficients
    //
    if(Model::escape_size())
      //
      for(int w = 0; w < chem_size; ++w)
	//
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  dtemp = vdot(m_inverse.column(w), &eigen_escape(0, e))
	    //
	    * std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() +
				   //
				   Model::bimolecular_size() + e)] = dtemp;
	}

    // output of steady state distributions associated with individual wells
    //
    IO::aux << "steady state distributions:\n";
    
    for(int g = 0; g < chem_size; ++g) {
      //
      IO::aux << "Initial well: " << Model::well(group_index[g]).name() << "\n";

      IO::aux  << std::setw(13) << "E, kcal/mol";

      for(int w = 0; w < Model::well_size(); ++w)
	//
	IO::aux << std::setw(13) << Model::well(w).name();

      IO::aux << "\n";

      for(int i = 0; i < well_size_max; ++i) {
	//
	IO::aux << std::setw(13) << energy_bin(i) / Phys_const::kcal;

	for(int w = 0; w < Model::well_size(); ++w) {
	  //
	  if(i < well(w).size()) {
	    //
	    IO::aux << std::setw(13) << vdot(m_inverse.column(g), &eigen_global(0, i + well_shift[w]))
	      //
	      * well(w).boltzman_sqrt(i) / std::sqrt(weight[g]);
	  }
	  else
	    //
	    IO::aux<< std::setw(13) << "0";
	}

	IO::aux << "\n";
      }

      IO::aux << "\n";
    }
    
    // well-to-bimolecular rate coefficients
    //
    for(int w = 0; w < chem_size; ++w)
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = wb_rate(w, p) * std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	
	rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
      }
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(bimolecular(p).weight() > 0.)
	//
	for(int w = 0; w < chem_size; ++w) {
	  //
	  dtemp = bw_rate(p, w) * std::sqrt(weight[w]) * energy_step() / bimolecular(p).weight() / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	}

    // product energy distribution
    //
    if(ped_out.is_open()) {
      //
      // well-to-bimolecular distribution
      //
      itemp = 0;
      
      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(outer_barrier(b).size() > itemp)
	  //
	  itemp = outer_barrier(b).size();

      int ener_index_max = itemp;
      
      if(Model::bimolecular_size()) {
	//
	mtemp.resize(itemp, Model::bimolecular_size());

	ped_out << "Well-to-bimolecular product energy distributions:\n";
	
	for(int ww = 0; ww < chem_size; ++ww) {
	  //
	  for(int e = 0; e < ener_index_max; ++e)
	    //
	    for(int p = 0; p < Model::bimolecular_size(); ++p) {
	      //
	      dtemp = 0.;
	      
	      for(int b = 0; b < Model::outer_barrier_size(); ++b) {
		//
		const int w = Model::outer_connect(b).first;
		
		if(p == Model::outer_connect(b).second && e < outer_barrier(b).size())
		  //
		  dtemp += vdot(m_inverse.column(ww), &eigen_global(0, e + well_shift[w]))
		    //
		    * global_bim(e + well_shift[w], p);
	      }
	      
	      mtemp(e, p) = dtemp;
	    }

	  // normalization
	  //
	  for(int i = 0; i < mtemp.size2(); ++i)
	    //
	    mtemp.column(i) /= max(mtemp.column(i));

	  ped_out << std::setw(13) << "E, kcal/mol";
	  
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(13) << Model::well(group_index[ww]).name() + "->"
	      //
	      + Model::bimolecular(p).name();
	  
	  ped_out << "\n";

	  for(int e = 0; e < ener_index_max; ++e) {
	    //
	    ped_out << std::setw(13) << energy_bin(e) / Phys_const::kcal;

	    for(int i = 0; i < mtemp.size2(); ++i)
	      //
	      ped_out << std::setw(13) << mtemp(e, i);
	    
	    ped_out << "\n";
	  }
	  ped_out << "\n";
	  //
	}//well-to-bimolecular
      }

      // well-to-escape distributions
      //
      if(Model::escape_size()) {
	//
	itemp = 0;
	
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  const int w = Model::escape_channel(e).first;
	  
	  if(well(w).size() > itemp)
	    //
	    itemp = well(w).size();
	}
	
	ener_index_max = itemp;
	
	mtemp.resize(itemp, chem_size * Model::escape_size());

	// distribution
	//
	for(int e = 0; e < ener_index_max; ++e) {
	  //
	  itemp = 0;
	  
	  for(int w = 0; w < chem_size; ++w)
	    //
	    for(int count = 0; count < Model::escape_size(); ++count, ++itemp) {
	      //
	      const int ew = Model::escape_channel(count).first;
	      
	      if(e < well(ew).size()) {
		//
		mtemp(e, itemp) = vdot(m_inverse.column(w), &eigen_global(0, e + well_shift[ew]))
		  //
		  * global_escape(e + well_shift[ew], count);
	      }
	      else
		//
		mtemp(e, itemp) =  0.;
	    }
	}

	// normalization
	//
	for(int i = 0; i < mtemp.size2(); ++i)
	  //
	  mtemp.column(i) /= max(mtemp.column(i));

	// output
	//
	ped_out << "Well-to-escape product energy distributions:\n";
	
	ped_out << std::setw(13) << "E, kcal/mol";
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  for(int e = 0; e < Model::escape_size(); ++e)
	    //
	    ped_out << std::setw(13) << Model::well(group_index[w]).name() + "->" + Model::escape_name(e);
	
	ped_out << "\n";
      
	for(int e = 0; e < ener_index_max; ++e) {
	  //
	  dtemp = energy_bin(e) / Phys_const::kcal;
	  
	  ped_out << std::setw(13) << dtemp;

	  for(int i = 0; i < mtemp.size2(); ++i)
	    //
	    ped_out << std::setw(13) << mtemp(e, i);
	  
	  ped_out << "\n";
	}
	
	ped_out << "\n";
	//
      }// well-to-escape distributions
      //
    }// product energy distributions
    //
  }//

  // hot distribution branching ratios
  //
  if(hot_index.size()) {
    //
    IO::log << IO::log_offset << "Hot distribution branching ratios:\n"
      //
	    << IO::log_offset //<< std::setprecision(6)
      //
	    << std::setw(5)  << "Well"
      //
	    << std::setw(13) << "E, kcal/mol";
    
    for(int w = 0; w < chem_size; ++w)
      //
      IO::log << std::setw(13) << Model::well(group_index[w]).name();

    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(13) << Model::bimolecular(p).name();

    for(int e = 0; e < Model::escape_size(); ++e)
      //
      IO::log << std::setw(13) << Model::escape_name(e);
    
    IO::log << "\n";
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
	//
	IO::log << IO::log_offset
	  //
		<< std::setw(5)  << Model::well(hit->first).name()
	  //
		<< std::setw(13) << energy_bin(hit->second) / Phys_const::kcal;
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  IO::log << std::setw(13) << vdot(m_direct.row(w), &eigen_hot(0, count)) * std::sqrt(weight[w]);
	
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  IO::log << std::setw(13) << triple_product(&eigen_bim(chem_size, p), &eigen_hot(chem_size, count),
						     //
						     relax_lave, relax_size);
	
	for(int e = 0; e < Model::escape_size(); ++e)
	  //
	  IO::log << std::setw(13) << triple_product(&eigen_escape(chem_size, e), &eigen_hot(chem_size, count),
						     //
						     relax_lave, relax_size);
	
	IO::log << "\n";
      }
  }

  // prompt isomerization
  //
  IO::log << IO::log_offset << "prompt isomerization/dissociation:\n";

  IO::log << IO::log_offset << std::setw(7) << "W\\W,P,E";

  if(chem_size)
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(13) << Model::well(w).name();

  IO::log << std::setw(13) << "Total Diss.";

  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    IO::log << std::setw(13) << Model::bimolecular(p).name();

  for(int e = 0; e < Model::escape_size(); ++e)
    //
    IO::log << std::setw(13) << Model::escape_name(e);
    
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(7) << Model::well(w).name();

    double diss = 1.;

    if(chem_size)
      //
      for(int v = 0; v < Model::well_size(); ++v) {
	//
	dtemp = 0.;

	for(int l = 0; l < chem_size; ++l)
	  //
	  dtemp += eigen_pop(l, w) * eigen_pop(l, v);

	// renormalization
	//
	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	dtemp *= well(v).weight_sqrt() * well(w).weight_sqrt() * energy_step() / well(w).real_weight();

	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	diss -= dtemp;

	IO::log << std::setw(13) << dtemp;
      }

    IO::log << std::setw(13) << diss;
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(13) << triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w),
						 //
						 relax_lave, relax_size) * energy_step()
	//
	* well(w).weight_sqrt() / well(w).real_weight();

    for(int e = 0; e < Model::escape_size(); ++e)
      //
      IO::log << std::setw(13) << triple_product(&eigen_escape(chem_size, e), &eigen_pop(chem_size, w),
						 //
						 relax_lave, relax_size) * energy_step()
	//
	* well(w).weight_sqrt() / well(w).real_weight();
      
    IO::log << "\n";
  }

  IO::log << "\n";
}

/***********************************************************************************************
 ******************************* PARTITIONING THE WELLS ****************************************
 ***********************************************************************************************/

void MasterEquation::set_well_partition_method (const std::string& method)  
{
  const char funame [] = "MasterEquation::set_well_partition_method: ";

  // active well projection threshold method
  if(method == "threshold") {
    well_partition_method =    threshold_well_partition;
    return;
  }

  // direct sort out method
  if(method == "sort_out") {
    well_partition_method = sort_well_partition;
    return;
  }

  // incremental method
  if(method == "incremental") {
    well_partition_method = incremental_well_partition;
    return;
  }

  // sequential method
  if(method == "sequential") {
    well_partition_method =   sequential_well_partition;
    return;
  }

  std::cerr << funame << ": unknown well partition method: " << method
	    << "; available methods: threshold (default), sort_out, incremental, sequential\n";
  throw Error::Input();
}

double MasterEquation::threshold_well_partition (const Lapack::Matrix& pop_chem, Partition& well_partition,
						 Group& bimolecular_group, const std::vector<double>& weight) 
{
  const char funame [] = "MasterEquation::threshold_well_partition: ";
  
  std::clock_t start_cpu = std::clock();
  std::time_t  start_time = std::time(0);

  //IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  const int well_size = weight.size() ? weight.size() : Model::well_size();

  if(pop_chem.size1() != well_size) {
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1() 
	      << " differs from the number of wells = " << well_size << "\n";
    throw Error::Logic();
  }

  if(chem_size >= well_size) {
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";
    throw Error::Logic();
  }

  // projection-ordered wells
  //
  std::multimap<double, int> proj_well_map;
  
  bimolecular_group.clear();
  
  for(int w = 0; w < well_size; ++w) {
    //
    proj_well_map.insert(std::make_pair(Group().insert(w).projection(pop_chem, weight), w));
    
    bimolecular_group.insert(w);
  }

  // large projection wells map
  //
  std::vector<int> well_map;
  
  for(std::multimap<double, int>::const_reverse_iterator i = proj_well_map.rbegin(); i != proj_well_map.rend(); ++i) {
    //
    if(i->first < well_projection_threshold && well_map.size() >= chem_size)
      //
      break;

    dtemp = i->first;
    
    well_map.push_back(i->second);
    
    bimolecular_group.erase(i->second);
  }

  if(dtemp < well_projection_threshold)
    //
    IO::log << IO::log_offset << funame << "WARNING: large well projection is smaller than the well projection threshold\n";
  
  double pmax = -1.;
  
  for(PartitionGenerator pg(chem_size, well_map.size()); !pg.end(); ++pg) {
    //
    // initialize new partition
    //
    Partition p(pg, well_map);

    // partition projection
    //
    dtemp = p.projection(pop_chem, weight);

    if(dtemp > pmax) {
      //
      pmax = dtemp;
      
      well_partition = p;
    }
  }

  while(bimolecular_group.size()) {
    //
    int smax, wmax;
    
    pmax = -1.;
    
    for(Git w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
      //
      for(int s = 0; s < chem_size; ++s) {
	//
	dtemp = Group(well_partition[s]).insert(*w).projection(pop_chem, weight)
	  //
	  - well_partition[s].projection(pop_chem, weight);
	
	if(dtemp > pmax) {
	  //
	  pmax = dtemp;
	  
	  smax = s;
	  
	  wmax = *w;
	}
      }
    }
    
    if(pmax < 0.)
      break;

    well_partition[smax].insert(wmax);
    bimolecular_group.erase(wmax);
  }

  pmax = (double)chem_size - well_partition.projection(pop_chem, weight);

  // output
  if(!weight.size()) {

    IO::log << IO::log_offset << funame <<  "starts\n";
    IO::log_offset.increase();

    // single well projections
    IO::log << std::setprecision(3) << std::fixed;
    IO::log << IO::log_offset << "single well projections onto chemical subspace:\n";
    IO::log << IO::log_offset << std::left << std::setw(12) << "well:" << std::right;
    for(std::multimap<double, int>::const_reverse_iterator i = proj_well_map.rbegin(); i != proj_well_map.rend(); ++i)
      IO::log << std::setw(6) << Model::well(i->second).name();
    IO::log << "\n";
    IO::log << IO::log_offset << std::left << std::setw(12) << "projection:" << std::right;
    for(std::multimap<double, int>::const_reverse_iterator i = proj_well_map.rbegin(); i != proj_well_map.rend(); ++i)
      IO::log << std::setw(6) << i->first;
    IO::log << "\n";
    IO::log << std::setprecision(6);
    IO::log.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);

    IO::log << IO::log_offset << "well projection threshold  = " << well_projection_threshold << "\n";
    IO::log << IO::log_offset << "partition projection error = " << pmax << "\n";

    IO::log_offset.decrease();
    IO::log << IO::log_offset << funame
	    << "done, cpu time[sec] = " << double(std::clock() - start_cpu) / CLOCKS_PER_SEC 
	    << ", elapsed time[sec] = "<< std::time(0) - start_time
	    << std::endl;
  }

  return pmax;
}

double MasterEquation::sort_well_partition (const Lapack::Matrix& pop_chem, Partition& well_partition,
					    Group& bimolecular_group, const std::vector<double>& weight) 
{
  const char funame [] = "MasterEquation::sort_well_partition: ";
  
  //IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  const int well_size = weight.size() ? weight.size() : Model::well_size();

  if(pop_chem.size1() != well_size) {
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1() 
	      << " differs from the number of wells = " << well_size << "\n";
    throw Error::Logic();
  }

  if(chem_size >= well_size) {
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";
    throw Error::Logic();
  }

  std::multimap<double, std::pair<Partition, Group> > high_part;

  std::clock_t start_cpu = std::clock();
  std::time_t  start_time = std::time(0);

  // partitioning the wells into equilibrated groups
  int part_size = chem_size;
  for(PartitionGenerator pg(part_size, well_size); !pg.end(); ++pg) {
    // initialize new partition
    Partition part(pg);

    // partition projection
    dtemp = part.projection(pop_chem, weight);

    if(high_part.size() < red_out_num)
      high_part.insert(std::make_pair(dtemp, std::make_pair(part, Group())));
    else if(dtemp > high_part.begin()->first) {
      high_part.erase(high_part.begin());
      high_part.insert(std::make_pair(dtemp, std::make_pair(part, Group())));
    }
  }

  // the case when some of the wells equilibrate with bimolecular products
  // is solved by allowing the bimolecular group and then removing it
  part_size = chem_size + 1;
  for(PartitionGenerator pg(part_size, well_size); !pg.end(); ++pg) {
    // initialize new partition
    Partition part(pg);
    for(int g = 0; g < part_size; ++g) {
      Partition q = part;
      Group bg = q[g];
      // removing the bimolecular group
      q.erase(q.begin() + g);
      // partition projection
      dtemp = q.projection(pop_chem, weight);
      if(high_part.size() < red_out_num)
	high_part.insert(std::make_pair(dtemp, std::make_pair(q, bg)));
      else if(dtemp > high_part.begin()->first) {
	high_part.erase(high_part.begin());
	high_part.insert(std::make_pair(dtemp, std::make_pair(q, bg)));
      }
    }
  }

  well_partition    = high_part.rbegin()->second.first;
  bimolecular_group = high_part.rbegin()->second.second;

  // output
  if(!weight.size()) {

    IO::log << IO::log_offset << funame <<  "starts\n";
    IO::log_offset.increase();

    IO::log << std::setprecision(3) << std::fixed;
    IO::log << IO::log_offset << "single well projections onto chemical subspace:\n";
    IO::log << IO::log_offset << std::left << std::setw(6) << "well:" << std::right;
    for(int w = 0; w < Model::well_size(); ++w)
      IO::log << std::setw(6) << Model::well(w).name();
    IO::log << "\n";
    IO::log << IO::log_offset << std::left << std::setw(6) << "proj:" << std::right;
    for(int w = 0; w < Model::well_size(); ++w) {
      Group g;
      g.insert(w);
      IO::log << std::setw(6) << g.projection(pop_chem);
    }
    IO::log << "\n";
    IO::log << std::setprecision(6);
    IO::log.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);

    IO::log << IO::log_offset << "closest well partitions:\n";
    std::multimap<double, std::pair<Partition, Group> >::const_reverse_iterator rit;
    IO::log << IO::log_offset 
	    << std::setw(15) << "reduction error" 
	    << std::setw(50) << "equilibrated_groups"
	    << std::setw(50) << "wells_equilibrated_with_products"
	    << "\n";
    for(rit = high_part.rbegin(); rit != high_part.rend(); ++rit) {
      IO::log << IO::log_offset << std::setw(15) << (double)chem_size - rit->first;
    
      // equilibrated bound species
      stemp.clear();
      for(Pit p = rit->second.first.begin(); p != rit->second.first.end(); ++p) {   
	if(p != rit->second.first.begin())
	  stemp += " ";	
	for(Git g = p->begin(); g != p->end(); ++g) {
	  if(g != p->begin())
	    stemp += "+";
	  stemp += Model::well(*g).name();
	}
      }

      if(stemp.size())
	IO::log << std::setw(50) << stemp;
      else
	IO::log << std::setw(50) << "---";

      // wells equilibrated with the products
      if(rit->second.second.size()) {
	stemp.clear();
	for(Git g = rit->second.second.begin(); g != rit->second.second.end(); ++g) {
	  if(g != rit->second.second.begin())
	    stemp += " ";
	  stemp += Model::well(*g).name();
	}

	IO::log << std::setw(50) << stemp;    
      }
      else
	IO::log << std::setw(50) << "---";
      IO::log << "\n";
    }

    IO::log_offset.decrease();
    IO::log << IO::log_offset << funame
	    << "done, cpu time[sec] = " << double(std::clock() - start_cpu) / CLOCKS_PER_SEC 
	    << ", elapsed time[sec] = "<< std::time(0) - start_time
	    << std::endl;
  }

  return (double)chem_size - high_part.rbegin()->first;
}

double MasterEquation::incremental_well_partition (const Lapack::Matrix& pop_chem, Partition& well_partition,
						   Group& bimolecular_group, const std::vector<double>& weight) 
{
  const char funame [] = "MasterEquation::increment_well_partition: ";
  
  //IO::Marker funame_marker(funame, IO::Marker::ONE_LINE);

  std::clock_t start_cpu = std::clock();
  std::time_t  start_time = std::time(0);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  const int well_size = weight.size() ? weight.size() : Model::well_size();

  if(pop_chem.size1() != well_size) {
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1() 
	      << " differs from the number of wells = " << well_size << "\n";
    throw Error::Logic();
  }

  if(chem_size >= well_size) {
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";
    throw Error::Logic();
  }

  well_partition.clear();
  bimolecular_group.clear();

  // single well ordered projections
  std::multimap<double, int> proj_well_map;
  for(int w = 0; w < well_size; ++w) {
    proj_well_map.insert(std::make_pair(Group().insert(w).projection(pop_chem, weight), w));
    bimolecular_group.insert(w);
  }

  well_partition.resize(chem_size);

  itemp = 0;
  for(std::multimap<double, int>::const_reverse_iterator i = proj_well_map.rbegin(); itemp < chem_size; ++i, ++itemp) {
    well_partition[itemp].insert(i->second);
    bimolecular_group.erase(i->second);
  }
    
  while(bimolecular_group.size()) {
    int smax, wmax;
    double proj_max = -1.;
    for(int s = 0; s < chem_size; ++s) {
      double sp = well_partition[s].projection(pop_chem, weight);
      for(Git w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
	dtemp = Group(well_partition[s]).insert(*w).projection(pop_chem, weight) - sp;
	if(dtemp > proj_max) {
	  proj_max = dtemp;
	  smax = s;
	  wmax = *w;
	}
      }
    }
    if(proj_max < 0.)
      break;

    well_partition[smax].insert(wmax);
    bimolecular_group.erase(wmax);
  }

  if(!weight.size()) {
    IO::log << IO::log_offset << funame 
	    << "... done, cpu time[sec] = " << double(std::clock() - start_cpu) / CLOCKS_PER_SEC 
	    << ", elapsed time[sec] = "<< std::time(0) - start_time
	    << std::endl;
  }

  return (double)chem_size - well_partition.projection(pop_chem, weight);
}

double MasterEquation::sequential_well_partition (const Lapack::Matrix& pop_chem, Partition& well_partition,
						  Group& bimolecular_group, const std::vector<double>& weight) 
{
  const char funame [] = "MasterEquation::sequent_well_partition: ";
  
  //IO::Marker funame_marker(funame, IO::Marker::ONE_LINE);

  std::clock_t start_cpu = std::clock();
  std::time_t  start_time = std::time(0);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  const int well_size = weight.size() ? weight.size() : Model::well_size();

  if(pop_chem.size1() != well_size) {
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1() 
	      << " differs from the number of wells = " << well_size << "\n";
    throw Error::Logic();
  }

  if(chem_size >= well_size) {
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";
    throw Error::Logic();
  }

  well_partition.clear();
  bimolecular_group.clear();

  well_partition.resize(chem_size);

  for(int w = 0; w < well_size; ++w)
    bimolecular_group.insert(w);
  
  for(int spec = 0; spec < chem_size; ++spec) {
    double spec_proj = 0.;
    while(bimolecular_group.size()) {
      int wmax;
      double proj_diff = -1.;
      for(Git w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
	dtemp = Group(well_partition[spec]).insert(*w).projection(pop_chem, weight) - spec_proj;
	if(dtemp > proj_diff) {
	  proj_diff = dtemp;
	  wmax = *w;
	}
      }
      if(proj_diff < 0.)
	break;

      well_partition[spec].insert(wmax);
      bimolecular_group.erase(wmax);
      spec_proj += proj_diff;
    }

    if(!well_partition[spec].size()) {
      std::cerr << funame << "empty species\n";
      throw Error::Logic();
    }
  }

  if(!weight.size()) {
    IO::log << IO::log_offset << funame 
	    << "... done, cpu time[sec] = " << double(std::clock() - start_cpu) / CLOCKS_PER_SEC 
	    << ", elapsed time[sec] = "<< std::time(0) - start_time
	    << std::endl;
  }

  return (double)chem_size - well_partition.projection(pop_chem, weight);
}

/********************************************************************************************
 ************ HIGH PRESSURE ANALYSIS OF THE HIERARCHY OF THE RELAXATION TIME SCALES *********
 ********************************************************************************************/

void  MasterEquation::high_pressure_analysis () 
{
  const char funame [] = "MasterEquation::high_pressure_analysis: ";

  int    itemp;
  double dtemp;

  // groups of equilibrated wells at high pressure
  std::map<int, HPWell>       hp_well;
  // groups of the wells equilibrated with the corresponding bimolecular products at high pressure
  std::vector<Group> hp_prod(Model::bimolecular_size());
  // barriers separating well groups and bimolecular products at high pressure 
  std::list<HPBarrier>        hp_barr(Model::inner_barrier_size() + Model::outer_barrier_size());

  typedef std::list<HPBarrier>::iterator Bit;

  // initializing
  // wells
  for(int w = 0; w < Model::well_size(); ++w) {
    itemp = w + Model::bimolecular_size();
    hp_well[itemp].insert(w);
    hp_well[itemp].weight = well(w).real_weight(); 
  } 

  // inner and outer barriers
  for(Bit b = hp_barr.begin(); b != hp_barr.end(); ++b) {
    itemp = std::distance(hp_barr.begin(), b);
    if(itemp < Model::inner_barrier_size()) {
      b->first  = Model::inner_connect(itemp).first  + Model::bimolecular_size();
      b->second = Model::inner_connect(itemp).second + Model::bimolecular_size();
      b->weight = inner_barrier(itemp).real_weight();
    }
    else {
      itemp    -= Model::inner_barrier_size();
      b->first  = Model::outer_connect(itemp).first  + Model::bimolecular_size();
      b->second = Model::outer_connect(itemp).second;
      b->weight = outer_barrier(itemp).real_weight();
    }
  }

  // find equilibrated species
  double eq_rate;
  std::pair<int, int> eq_spec;
  double rate;
  for(Bit b = hp_barr.begin(); b != hp_barr.end(); ++b) {
    rate = b->weight / hp_well[b->first].weight;
    if(b->second >= Model::bimolecular_size()) //rate between two bound species
      rate += b->weight / hp_well[b->second].weight;
    
    if(hp_barr.begin() == b || rate > eq_rate) {
      eq_rate = rate;
      eq_spec = *b;
    }
  }

  if(eq_spec.second < Model::bimolecular_size()) {// bound and bimolecular species equilibrated
    hp_prod[eq_spec.second].insert(hp_well[eq_spec.first]);
    hp_well.erase(eq_spec.first);

 
    Bit b = hp_barr.begin(); 
    while(b != hp_barr.end()) {
      if(b->first == eq_spec.first) { 
	b->first = b->second;
	b->second = eq_spec.second;
      }
      if(b->second == eq_spec.first)
	b->second = eq_spec.second;
      if(b->first < Model::bimolecular_size())
	b = hp_barr.erase(b);
      else
	++b;
    }
  }
  else {// two bound species equilibrated
    itemp = hp_well.rbegin()->first + 1;
    hp_well[itemp].weight = hp_well[eq_spec.first].weight 
      + hp_well[eq_spec.second].weight;
    hp_well[itemp].insert(hp_well[eq_spec.first]);
    hp_well[itemp].insert(hp_well[eq_spec.second]);
    hp_well.erase(eq_spec.first);
    hp_well.erase(eq_spec.second);

    Bit b = hp_barr.begin();
    while(b != hp_barr.end()) {
      if(b->first == eq_spec.first || b->first == eq_spec.second) 
	b->first = hp_well.rbegin()->first;
      if(b->second == eq_spec.first || b->second == eq_spec.second) 
	b->second = hp_well.rbegin()->first;
      if(b->first == b->second)
	b = hp_barr.erase(b);
      else
	++b;
    }
  }

  // remove duplicate barriers 
  for(Bit b = hp_barr.begin(); b != hp_barr.end(); ++b) {
    Bit b1 = b;
    std::advance(b1, 1);
    while(b1 != hp_barr.end())
      if(b1->first == b->first  && b1->second == b->second ||
	 b1->first == b->second && b1->second == b->first) {
	b->weight += b1->weight;
	b1 = hp_barr.erase(b1);
      }
      else
	++b1;
  }
}

/***********************************************************************************************
 ************************************** HELPERS ************************************************
 ***********************************************************************************************/

/************** Generator of partitions of n objects into m groups ****************/

MasterEquation::PartitionGenerator::PartitionGenerator (int m, int n) 
  : _group_index(n), _frame(m + 1), _end(false), _partition_size(m)
{
  const char funame [] = "MasterEquation::PartitionGenerator::PartitionGenerator: ";
  
  if(m > n) {
    std::cerr << funame << "number of groups is bigger than the number of objects\n";
    throw Error::Range();
  }

  _frame[m] = n;  
  for(int i = 0; i < m; ++i)
    _group_index[i]= _frame[i] = i;
}

void MasterEquation::PartitionGenerator::operator++ ()
{
  if(_end)
    return;

  // iterate within a frame
  for(int i = _frame[1] + 1, j = 1; i < _group_index.size(); ++i)
    if(i == _frame[j + 1])
      ++j;
    else if(_group_index[i] == j)
      _group_index[i] = 0;
    else {
      ++_group_index[i];
      return;
    }
    
  // iterate over different frames
  for(int j = 1; j < _partition_size; ++j)
    if(_frame[j] + 1 == _frame[j + 1])
      _frame[j] = _frame[j - 1] + 1;
    else {
      ++_frame[j];
      // initialize the index array
      for(int i = 1, k = 1; i < _group_index.size(); ++i)
	if(i == _frame[k]) {
	  _group_index[i] = k;
	  ++k;
	}
	else
	  _group_index[i] = 0;
      return;
    }

  _end = true;
}

/*************** Generator of group of m elements from the pool of n elements *******************/

MasterEquation::GroupGenerator::GroupGenerator (int m, int n) 
  : _index(m + 1), _end(false), _size(m)
{
  const char funame [] = "MasterEquation::GroupGenerator::GroupGenerator: ";
  
  if(m > n) {
    std::cerr << funame << "group size is bigger than the pool size\n";
    throw Error::Range();
  }

  for(int i = 0; i < m; ++i)
    _index[i] = i;
  _index[m] = n;
}

void MasterEquation::GroupGenerator::operator++ ()
{
  bool btemp;

  if(_end)
    return;

  btemp = false;
  for(int i = 0; i < size(); ++i)
    if(_index[i] + 1 == _index[i + 1])
      if(i)
	_index[i] = _index[i - 1] + 1;
      else
	_index[i] = 0;
    else {
      ++_index[i];
      btemp = true;
      break;
    }

  if(btemp)
    return;
  _end = true;
}

/*********************************** Partitioning the wells *************************************/

MasterEquation::Partition::Partition (const PartitionGenerator& g, const std::vector<int>& map) : std::vector<Group>(g.partition_size())
{
  const char funame[] = "MasterEquation::Partition::Partition (const PartitionGenerator&): ";

  if(!map.size())
    for(int i = 0; i < g.size(); ++i)
      (*this)[g[i]].insert(i);
  else if(map.size() == g.size())
    for(int i = 0; i < g.size(); ++i)
      (*this)[g[i]].insert(map[i]);
  else {
    std::cerr << funame << "map dimension mismatch\n";
    throw Error::Range();
  }
}

std::vector<double> MasterEquation::Partition::weight () const
{
  std::vector<double> res;
  res.reserve(size());

  for(Pit g = begin(); g != end(); ++g)
    res.push_back(g->weight());
  return res;
}

std::vector<double> MasterEquation::Partition::real_weight () const
{
  std::vector<double> res;
  res.reserve(size());

  for(Pit g = begin(); g != end(); ++g)
    res.push_back(g->real_weight());
  return res;
}

std::vector<int> MasterEquation::Partition::group_index () const 
{
  std::vector<int>  res;
  res.reserve(size());

  for(Pit g = begin(); g != end(); ++g)
    res.push_back(g->group_index());
  
  return res;
}

double MasterEquation::Partition::projection (const Lapack::Matrix& pop_chem, const std::vector<double>& weight) const 
{
  double res = 0.;
  for(const_iterator g = begin(); g != end(); ++g)
    res += g->projection(pop_chem, weight);

  return res;
}
	
Lapack::Matrix MasterEquation::Partition::basis () const 
{
  double dtemp;

  Lapack::Matrix res(Model::well_size(), size());
  res = 0.;
  for(Pit g = begin(); g != end(); ++g)
    res.column(g - begin()) = g->basis();
  /*
    if(g->size() > 1) {
    dtemp = 0.;
    for(Git w = g->begin(); w != g->end(); ++w)
    dtemp += well(*w).weight();
    dtemp = std::sqrt(dtemp);
    for(Git w = g->begin(); w != g->end(); ++w)
    res(*w, g - begin()) = well(*w).weight_sqrt() / dtemp;
    }
    else
    res(*g->begin(), g - begin()) = 1.;
  */

  return res;
}

/**************************************** Well Group *******************************************/

MasterEquation::Group::Group (int w)
{
  insert(w);
}

double MasterEquation::Group::weight () const
{
  double res = 0.;
  for(const_iterator w = begin(); w != end(); ++w) 
    res += well(*w).weight();
  return res;
}

double MasterEquation::Group::real_weight () const
{
  double res = 0.;
  for(const_iterator w = begin(); w != end(); ++w) 
    res += well(*w).real_weight();
  return res;
}

int MasterEquation::Group::group_index () const
{
  double dtemp, wmax;

  int res = -1;
  for(const_iterator w = begin(); w != end(); ++w) {
    dtemp =  well(*w).real_weight();
    if(w == begin() || dtemp > wmax) {
      wmax = dtemp;
      res = *w;
    }
  }
  return res;
}

Lapack::Vector MasterEquation::Group::basis () const
{
  Lapack::Vector res(Model::well_size());
  res = 0.;

  if(size() == 1) {
    res[*begin()] = 1.;
    return res;
  }
    
  double dtemp = std::sqrt(weight());

  for(const_iterator w = begin(); w != end(); ++w)
    res[*w] = well(*w).weight_sqrt() / dtemp;
  return res;
}

MasterEquation::Group MasterEquation::Group::insert (const Group& group)
{
  for(const_iterator g = group.begin(); g != group.end(); ++g) 
    insert(*g);
  
  return *this;
}

MasterEquation::Group MasterEquation::Group::erase (const Group& group)
{
  for(const_iterator g = group.begin(); g != group.end(); ++g) 
    erase(*g);
  
  return *this;
}

double MasterEquation::Group::projection (const Lapack::Matrix& pop_chem, const std::vector<double>& weight) const
  
{
  
  const char funame [] = "MasterEquation::Group::projection: ";

  double res, dtemp;

  const int well_size = weight.size() ? weight.size() : Model::well_size();

  if(!size()) {
    std::cerr << funame << "group is empty\n";
    throw Error::Logic();
  }
  
  if(pop_chem.size1() != well_size) {
    std::cerr << funame << "pop_chem: well dimension = " << pop_chem.size1() << " mismatch\n";
    throw Error::Logic();
  }

  if(size() == 1)
    return vdot(pop_chem.row(*begin()));

  res = 0.;
  for(int l = 0; l < pop_chem.size2(); ++l) {
    dtemp = 0.;
    for(const_iterator w = begin(); w != end(); ++w)
      if(!weight.size())
	dtemp += well(*w).weight_sqrt() * pop_chem(*w, l);
      else
	dtemp += std::sqrt(weight[*w]) * pop_chem(*w, l);
    res += dtemp * dtemp;
  }

  // normalization
  dtemp = 0.;
  for(const_iterator w = begin(); w != end(); ++w)
    if(!weight.size())
      dtemp += well(*w).weight();
    else
      dtemp += weight[*w];

  res /= dtemp;
  return res;
}

bool MasterEquation::Group::contain (const Group& group) const
{
  for(Git w = group.begin(); w != group.end(); ++w)
    if(find(*w) == end())
      return false;
  return true;
}

double MasterEquation::Group::projection (const Lapack::Matrix& eigen_global, int chem_size, 
					  const std::vector<Group>& bound_group,
					  const std::vector<double>& group_weight_sqrt) const
  
{
  
  const char funame [] = "MasterEquation::Group::projection: ";

  double dtemp;

  if(!size()) {
    std::cerr << funame << "group is empty\n";
    throw Error::Logic();
  }
  
  const int global_size = eigen_global.size2();

  if(global_size != bound_group.size() || global_size != group_weight_sqrt.size()) {
    std::cerr << funame << "global dimensions mismatch\n";
    throw Error::Logic();
  }

  if(chem_size > eigen_global.size1()) {
    std::cerr << funame << "chemical subspace dimension out of range\n";
    throw Error::Range();      
  }

  double res = 0.;
  for(int l = 0; l < chem_size; ++l) {
    dtemp = 0.;
    for(int g = 0; g < global_size; ++g)
      if(contain(bound_group[g]))
	dtemp += group_weight_sqrt[g] * eigen_global(l, g);
    res += dtemp * dtemp;
  }

  // normalization
  dtemp = 0.;
  for(const_iterator w = begin(); w != end(); ++w)
    dtemp += well(*w).weight();

  res /= dtemp;
  return res;
}

int MasterEquation::Group::operator[] (int i) const
{
  const_iterator g = begin();
  std::advance(g, i);
  return *g;
}
