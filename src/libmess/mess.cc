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
#include <limits>

#include "mess.hh"
#include "units.hh"
#include "key.hh"
#include "io.hh"
#include "shared.hh"
#include "offload.hh"
#include "mpack.hh"
#include "limits.hh"

namespace MasterEquation {

#if defined(WITH_MPACK) || defined(WITH_MPLAPACK)
  
#pragma omp declare reduction (+: dd_real: omp_out += omp_in) initializer(omp_priv = 0.)
	
  std::ostream& operator<< (std::ostream& to, const dd_real& dd) { return to << dd._hi(); }

  IO::LogOut& operator<< (IO::LogOut& to, const dd_real& dd) { return to << dd._hi(); }

#define CONVERT_DD(a)  (a)._hi()

#else
  
#define CONVERT_DD(a)  a

#endif
  
  /*********************************** AUXILIARY OUTPUT *************************************/

  std::string save_kinetic_matrix;
  
  std::ofstream eval_out;// eigenvalues output
  
  int  red_out_num = 5; // number of reduction schemes to print 

  std::ofstream evec_out;
  
  int evec_out_num = 0;

  std::ofstream arr_out; // arrhenius 

  /********************************* USER DEFINED PARAMETERS ********************************/

  int                                                       use_mp = 0;

  double                                                    reference_temperature = -1.;

  double                                                    reference_pressure    = -1.;
  
  // temperature
  //
  double                                                     _temperature     = -1.;

  // pressure
  //
  double                                                     _pressure        = -1.;

  double                                                     energy_step_over_temperature = -1.;

  //  energy window
  //
  double                                                     energy_window_over_temperature = -1.;

  // upper part of the energy window
  //
  double                                                     excess_energy_over_temperature = -1.;

  // lower part of the energy window
  //
  double                                                     energy_cutoff_over_temperature = -1.;

  // global energy cutoff
  //
  bool                                                       is_global_cutoff = false;
  
  double                                                     lower_global_cutoff;

  // well cutoff parameter
  //
  double                                                     well_cutoff       = -1.;
  
  // energy reference and step
  //
  double                                                     _energy_reference;
  
  double                                                     _energy_step;
  
  // highest chemecial eigenvalue parameter
  //
  double                                                     chemical_threshold     = -2.;
  
  // smallest chemical eigenvalue
  //
  double                                                     chemical_tolerance     = -1.;
  
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
  
  // pressure units
  //
  int                                                        pressure_unit = BAR;

  // well partition method
  //
  double (*well_partition_method) (Lapack::Matrix, Partition&, Group&) = threshold_well_partition;

  // well partition thresholds
  //
  double                                                     well_projection_threshold = 0.2;

  double                                                     well_projection_tolerance = .0001;
  
  // well extention parameter
  //
  double                                                     well_extension = -1.;

  // well extension correction
  //
  double                                                     well_ext_corr = 0.;

  // minimal relaxational eigenvalue to collisional frequency ratio
  //
  double                                                     weak_collision_factor = 0.1;

  // high pressure regime threshold
  //
  double                                                     high_pressure_threshold = 100.;

  // time propagation step in maximal relaxation frequency units
  //
  double                                                     time_propagation_step  = 0.01;

  // time propagation limit in collision frequency units
  //
  double                                                     time_propagation_limit = 100.;

  // time-interval between printouts
  //
  double                                                     time_propagation_interval = -1.;

  // generate CRM basis
  //
  int                                                        with_crm_basis            = 0;
  
  /********************************* INTERNAL PARAMETERS ************************************/

  // inner setting check
  //
  bool  _isset = false;

  // numerical accuracy
  //
  const double epsilon = 1.e-10;

  /*******************************************************************************************/

  std::vector<double> _thermal_factor;

  std::vector<SharedPointer<Well> >  _well;
  const Well& well (int w) { return *_well[w]; } 

  std::vector<SharedPointer<Barrier> > _inner_barrier;
  const Barrier& inner_barrier (int p) { return *_inner_barrier[p]; }

  std::vector<SharedPointer<Barrier> > _outer_barrier;
  const Barrier& outer_barrier (int p) { return *_outer_barrier[p]; }
  
  // Boltzman exponent: exp(-E/T)
  //
  double        thermal_factor (int);
  void   resize_thermal_factor (int);
  void    reset_thermal_factor (int = 0);

  bool isset () { return _isset; }
  
  double energy_reference () { return _energy_reference; }

  double energy_bin (int i) { return energy_reference() - i * energy_step(); }

  // temperature
  //
  double temperature ()
  {
    const char funame [] = "MasterEquation::temperature: ";
    
    if(_temperature <= 0.) {
      //
      std::cerr << funame << "not initialized\n";
      
      throw Error::Init();
    }
    
    return _temperature;
  }

  // energy step
  //
  double energy_step ()
  {
    return _energy_step;
  }

  // pressure
  //
  double pressure ()
  {
    const char funame [] = "MasterEquation::pressure: ";
    
    if(_pressure <= 0.) {
      ///
      std::cerr << funame << "not initialized\n";
      
      throw Error::Init();
    }
    
    return _pressure;
  }
  
  void set_energy_step () 
  {
    _energy_step = energy_step_over_temperature * temperature();
  
    reset_thermal_factor(); 
  }

  void set_temperature (double t) 
  { 
    _isset  = false;
  
    _temperature = t;

    _energy_step = energy_step_over_temperature * t;
  }

  void set_energy_reference (double er) 
  {
    _isset = false;
  
    _energy_reference = er;
  }

  void set_pressure (double p) 
  {
    _pressure = p;
  }

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

int MasterEquation::get_precision ()
{
  const char funame [] = "MasterEquation::get_precision: ";
  
  static const double bits2digits = std::log10(2.);
  
  //IO::log << IO::log_offset << "numerical precision: digits10 = ";
  
#if defined(WITH_MPACK) || defined(WITH_MPLAPACK)

  switch(Mpack::mp_type) {
    //
  case Mpack::DOUBLE:
    //
    return std::numeric_limits<double>::digits10;

  case Mpack::DD:
    //
    return 31;

  case Mpack::QD:
    //
    return 62;
    
  case Mpack::MPFR:
      //
#ifdef WITH_MPACK
    
      return int(mpfr::mpreal::get_default_prec() * bits2digits);
      
#else

      return int(mpfr::mpreal::default_prec * bits2digits);
      
  case Mpack::FLOAT64X:
    //
    return 18;

#endif
      
  case Mpack::GMP:
    //
    return int(mpf_get_default_prec() * bits2digits);

  case Mpack::FLOAT128:
    //
    return 62;

  default:
    //
    std::cerr << funame << "should not be here\n";

    throw Error::Logic();
  }
  
#else

  return std::numeric_limits<double>::digits10;

#endif

}
  
void MasterEquation::set_precision (int prec)
{
  static const double bits2digits = std::log10(2.);
  
#if defined(WITH_MPACK) || defined(WITH_MPLAPACK)

  switch(Mpack::mp_type) {
    //
  case Mpack::MPFR:
      //
#ifdef WITH_MPACK
    
    mpfr::mpreal::set_default_prec(int(prec / bits2digits));
    
#else
    
    mpfr::mpreal::default_prec = int(prec / bits2digits);

#endif
    
    break;

  case Mpack::GMP:
    //
    mpf_set_default_prec(int(prec / bits2digits));

    break;
  }
  
#endif

}

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
    //
    return;

  int emin = _thermal_factor.size();
  //
  _thermal_factor.resize(s);

  for(int e = emin; e < _thermal_factor.size(); ++e)
    //
    _thermal_factor[e] = std::exp(double(e) * energy_step_over_temperature);
}

void  MasterEquation::reset_thermal_factor(int s)
{
  _thermal_factor.resize(s);
  
  if(!s)
    //
    return;

  for(int e = 0; e < _thermal_factor.size(); ++e)
    //
    _thermal_factor[e] = std::exp(double(e) * energy_step_over_temperature);
}

/************************************* PRODUCT ENERGY DISTRIBUTION PAIRS ********************************************/

void MasterEquation::set_ped_pair(const std::vector<std::string>& ped_spec) 
{
  const char funame [] = "MasterEquation::set_ped_pair: ";

  IO::Marker marker(funame, IO::Marker::NOTIME);

  std::string::size_type itemp, start, next;
  int spec;

  for(std::vector<std::string>::const_iterator ped = ped_spec.begin(); ped != ped_spec.end(); ++ped) {
    //
    std::vector<int> pair;
    
    start = 0;
    
    while(start < ped->size()) {
      //
      next = ped->size();
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	itemp = ped->find(Model::bimolecular(p).name(), start);
	
	if(itemp < next || itemp == next && itemp < ped->size() &&
	   //
	   Model::bimolecular(p).name().size() > Model::bimolecular(spec).name().size()) {
	  //
	  next = itemp;
	  
	  spec = p;
	}
      }
      
      if(next >= ped->size())
	//
	break;
      
      pair.push_back(spec);
      
      start = next + Model::bimolecular(spec).name().size();
    }
    
    if(pair.size() != 2) {
      //
      std::cerr << funame << "wrong number of bimolecular species (" << pair.size() << ") in " << *ped << " token\n";
      
      throw Error::Init();
    }
    
    ped_pair.push_back(std::make_pair(pair.front(), pair.back()));
  }

  if(!ped_pair.size()) {
    //
    std::cerr << funame << "no pairs\n";
    
    throw Error::Init();
  }
  
  IO::log << IO::log_offset << "ped pairs:\n";
  
  for(int i = 0; i < ped_pair.size(); ++i)
    //
    IO::log << IO::log_offset << "   "
      //
	    << Model::bimolecular(ped_pair[i].first).name()
      //
	    << "->"  << Model::bimolecular(ped_pair[i].second).name() << "\n";
}

/***************************************** DEFAULT REDUCTION SCHEME *************************************************/

void MasterEquation::set_default_partition(const std::vector<std::string>& scheme) 
{
  const char funame [] = "MasterEquation::set_default_partition: ";

  IO::Marker marker(funame, IO::Marker::NOTIME);

  std::string::size_type itemp, start, next;
  int spec;

  std::set<int> pool;
  
  for(std::vector<std::string>::const_iterator g = scheme.begin(); g != scheme.end(); ++g) {
    //
    Group group;
    
    start = 0;
    
    while(start < g->size()) {
      //
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
      //
      std::cerr << funame << "empty group: " << *g << "\n";
      
      throw Error::Init();
    }
    default_partition.push_back(group);
  }
  if(default_partition.size() >= Model::well_size()) {
    //
    std::cerr << funame << "default reduction scheme: no reduction\n";

    throw Error::Init();   
  }

  // output
  //
  IO::log << IO::log_offset << "default reduction scheme: ";
  
  for(int g = 0; g < default_partition.size(); ++g) {
    //
    if(!g)
      //
      IO::log << " ";
    
    for(Group::const_iterator w = default_partition[g].begin(); w != default_partition[g].end(); ++w) {
      //
      if(w != default_partition[g].begin())
	//
	IO::log << "+";
      
      IO::log << Model::well(*w).name();
    }
  }
  IO::log << "\n";
}

void MasterEquation::set_default_chem_size (int s)
{
  if(s < 0) {
    //
    std::cerr << "MasterEquation::set_default_chem_size: out of range\n";
    
    throw Error::Range();
  }

  _default_chem_size = s;
}

/********************************************************************************************
 ************************************** REACTIVE COMPLEX ************************************
 ********************************************************************************************/


bool MasterEquation::ReactiveComplex::there_are_bound_groups () const
{
  const char funame [] = "MasterEquation::ReactiveComplex::there_are_bound_groups: ";

#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#endif
  
  if(!isinit) {
    //
    std::cerr << funame << "reactive complex has not been initialized\n";
    
    throw Error::Init();
  }

  if(!kinetic_basis.size()) {
    //
    std::cerr << funame << "kinetic matrices have not been initialized\n";
    
    throw Error::Init();
  }

  if(chemical_threshold <= 0. || chemical_threshold >= 1.) {
    //
    std::cerr << funame << "chemical threshold out of range: " << chemical_threshold << "\n";

    throw Error::Logic();
  }
  
  if(global_size < well_size() + 1)
    //
    return false;
  
  LAPACK::Vector eigenval = kin_mat.eigenvalues();

  const float_t& relax_eval_min = eigenval[well_size()];
  
  if(eigenval[0] > chemical_threshold * relax_eval_min)
    //
    return false;

  return true;
}

MasterEquation::ReactiveComplex::RateData MasterEquation::ReactiveComplex::well_reduction_method (int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::well_reduction_method: ";

  using Model::operator<<;

  if(!isinit) {
    //
    std::cerr << funame << "reactive complex has not been initialized\n";
    
    throw Error::Init();
  }

  if(!kinetic_basis.size()) {
    //
    std::cerr << funame << "kinetic matrices have not been initialized\n";
    
    throw Error::Init();
  }
    
#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int                  itemp;
  float_t              dtemp;
  bool                 btemp;
  std::string          stemp;
  
  RateData rate_data;

  rate_data._bim_index_map = _bim_index_map;

  rate_data.graph = graph;

  // collision frequencies
  //
  if(!(flags & NOPRINT)) {
    //
    std::map<int, std::set<int> > m;

    for(int w = 0; w < well_size(); ++w) {
      //
      itemp = std::ceil(well(w).collision_frequency(pressure) / well(0).collision_frequency(pressure) * 100. - 0.5);

      m[itemp].insert(w);
    }

    if(m.size() > 1) {
      //
      IO::log << IO::log_offset << "collision frequencies[1/sec]:";

      for(std::map<int, std::set<int> >::const_iterator i = m.begin(); i != m.end(); ++i)
	//
	IO::log << "  " << well(0).collision_frequency(pressure) / Phys_const::herz * i->first / 100. << "/" << index_to_well(i->second);

      IO::log << "\n";
    }
    else
      //
      IO::log << IO::log_offset << "collision frequency = " << well(0).collision_frequency(pressure) / Phys_const::herz << " 1/sec\n";
  }
  
  // reactants groups
  //
  Model::ChemGraph rg1, rg2;

  Model::landscape_t landscape = Model::kinetic_landscape(temperature);
    
  Model::landscape_t::const_iterator cbi = landscape.find(control_barrier);

  switch(control_barrier.first) {
    //
  case Model::INNER:
    //
    if(cbi == landscape.end() || cbi->second.size() != 2) {
      //
      std::cerr << funame << "(inner) rate-limiting barrier index out of range: " << control_barrier.second << "\n";

      throw Error::Logic();
    }
    
    rg1 = cbi->second.front();
    
    for(std::set<int>::const_iterator b = graph.inner_set.begin(); b != graph.inner_set.end(); ++b)
      //
      if(rg1.well_set.find(Model::inner_connect(*b).first)  != rg1.well_set.end() &&
	 //
	 rg1.well_set.find(Model::inner_connect(*b).second) != rg1.well_set.end() &&
	 //
	 rg1.inner_set.find(*b) == rg1.inner_set.end())
	//
	rg1.inner_set.insert(*b);
	  
    for(std::set<int>::const_iterator b = graph.outer_set.begin(); b != graph.outer_set.end(); ++b)
      //
      if(rg1.well_set.find(Model::outer_connect(*b).first)  != rg1.well_set.end() && rg1.outer_set.find(*b) == rg1.outer_set.end())
	//
	rg1.outer_set.insert(*b);

    rg2 = cbi->second.back();

    for(std::set<int>::const_iterator b = graph.inner_set.begin(); b != graph.inner_set.end(); ++b)
      //
      if(rg2.well_set.find(Model::inner_connect(*b).first)  != rg2.well_set.end() &&
	 //
	 rg2.well_set.find(Model::inner_connect(*b).second) != rg2.well_set.end() &&
	 //
	 rg2.inner_set.find(*b) == rg2.inner_set.end())
	//
	rg2.inner_set.insert(*b);
	  
    for(std::set<int>::const_iterator b = graph.outer_set.begin(); b != graph.outer_set.end(); ++b)
      //
      if(rg2.well_set.find(Model::outer_connect(*b).first)  != rg2.well_set.end() && rg2.outer_set.find(*b) == rg2.outer_set.end())
	//
	rg2.outer_set.insert(*b);

    break;
    //
  case Model::OUTER:
    //
    if(cbi == landscape.end() || cbi->second.size() != 1) {
      //
      std::cerr << funame << "(outer) rate-limiting barrier index out of range: " << control_barrier.second << "\n";

      throw Error::Logic();
    }
    
    rg1 = cbi->second.front();

    for(std::set<int>::const_iterator b = graph.inner_set.begin(); b != graph.inner_set.end(); ++b)
      //
      if(rg1.well_set.find(Model::inner_connect(*b).first)  != rg1.well_set.end() &&
	 //
	 rg1.well_set.find(Model::inner_connect(*b).second) != rg1.well_set.end() &&
	 //
	 rg1.inner_set.find(*b) == rg1.inner_set.end())
	//
	rg1.inner_set.insert(*b);
	  
    for(std::set<int>::const_iterator b = graph.outer_set.begin(); b != graph.outer_set.end(); ++b)
      //
      if(rg1.well_set.find(Model::outer_connect(*b).first)  != rg1.well_set.end() && rg1.outer_set.find(*b) == rg1.outer_set.end())
	//
	rg1.outer_set.insert(*b);

  }
  
  /********************************************************************************************************
   ************************************ NO KINETICALLY ACTIVE STATES **************************************
   ********************************************************************************************************/
  
  if(global_size <= well_size()) {
    //
    group_t owg;

    IO::log << IO::log_offset << "WARNING: " << "global size: " << global_size << " <= # of wells: " << well_size() << "\n";
    
    for(int i = 0; i < well_size(); ++i)
      //
      owg.insert(index_to_well(i));

    rate_data.bimolecular_group = owg;

    if(!bimolecular_size()) {
      //
      std::cerr << "no kinetically active states & no bimolecular channels... Hm... bail out\n";

      throw Error::Logic();
    }

    // bimolecular-to-bimolecular rates
    //
    LAPACK::SymmetricMatrix bb_rate(bimolecular_size());
    
    bb_rate = 0.;
    
    for(int e = 0; e < ener_index_max; ++e) {
      //
      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	//
	LAPACK::Vector vtemp(bimolecular_size(), 0.);
    
	for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	  //
	  const int w = bit->first.first;
      
	  const int p = bit->first.second;

	  if(e >= bit->second.size())
	    //
	    continue;

	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "well size assiciated with the barrier smaller than the barrier size\n";

	    throw Error::Logic();
	  }
	  
	  vtemp[p] += kinetic_basis[e].eigenvector(i->second, k)
	    //
	    * bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e] / 2. / M_PI;
	}

	for(int p = 0; p < bimolecular_size(); ++p)
	  //
	  for(int q = p; q < bimolecular_size(); ++q)
	    //
	    bb_rate(p, q) += vtemp[p] * vtemp[q] / kinetic_basis[e].eigenvalue[k];
      }
    }

    // normalization
    //
    bb_rate *= energy_step / std::exp(energy_reference / temperature);

    rate_data.bb_rate.resize(bimolecular_size());

    for(int p = 0; p < bimolecular_size(); ++p)
      //
      for(int q = p; q < bimolecular_size(); ++q)
	//
	rate_data.bb_rate(p, q) = CONVERT_DD(bb_rate(p, q));
    
    return rate_data;
  }
  
  /*************************************************************************************************
   *************************************** GLOBAL MATRICES *****************************************
   *************************************************************************************************/

  LAPACK::Vector eigenval;
  
  LAPACK::Matrix global_eigen(global_size);

  {

#ifdef DEBUG
    //
    IO::Marker marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);
    //
#else
    //
    IO::Marker marker;

    if(!(flags & NOPRINT)) {
      //
      marker.init("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);
    }
    //
#endif
  
    eigenval = kin_mat.eigenvalues(&global_eigen);
  }

  const float_t& relax_eval_min = eigenval[well_size()];
  
  const float_t& relax_eval_max = eigenval.back();

  if(!(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "minimal relaxation eigenvalue (relaxation limit) over collision frequency = "
      //
	    << relax_eval_min / well(0).collision_frequency(pressure) << "\n";
  
    IO::log << IO::log_offset << "maximal relaxation eigenvalue over collision frequency = "
      //
	    << relax_eval_max / well(0).collision_frequency(pressure) << "\n";
  }
  
  /*************************************************************************************************
   *********************************** EIGENVECTOR PROJECTIONS *************************************
   *************************************************************************************************/

  // eigenvectors projection on thermal subspace;
  //
  LAPACK::Matrix eigen_pop(global_size, well_size());

  eigen_pop = 0.;

#pragma omp parallel for default(shared) private(dtemp, itemp) schedule(static)
  
  for(int l = 0; l < global_size; ++l) {
    //
    for(int w = 0; w < well_size(); ++w) {
      //
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	if(i == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "eigen_pop: well is not in the kinetic basis space\n";

	  throw Error::Logic();
	}

	dtemp = 0.;
	
	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	dtemp *= well(w).boltzman_sqrt[e] / well(w).weight_sqrt;

	eigen_pop(l, w) += dtemp;
      }
      //
    }//
    //
  }//

  std::vector<float_t> relaxation_projection(well_size());
  
  for(int l = 0; l < well_size(); ++l)
    //
    relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));

  if(!(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "eigenstate populations:\n"
      //
	    << IO::log_offset
      //
	    << std::setw(5)  << "L"
      //
	    << std::setw(Model::log_precision + 7) << "*R"
      //
	    << std::setw(Model::log_precision + 7) << "*P";
  
    for(int w = 0; w < well_size(); ++w)
      //
      IO::log << std::setw(7) << well_model(w).short_name();
  
    IO::log << "\n";

    for(int l = 0; l < well_size(); ++l) {
      //
      float_t nfac;
    
      for(int w = 0; w < well_size(); ++w) {
	//
	dtemp = eigen_pop(l, w) * well(w).weight_sqrt;

	dtemp = dtemp < 0. ? -dtemp : dtemp;

	if(!w || dtemp > nfac)
	  //
	  nfac = dtemp;
      }
    
      IO::log << IO::log_offset
	//
	      << std::setw(5)  << l
	//
	      << std::setw(Model::log_precision + 7) << CONVERT_DD(eigenval[l] / relax_eval_min)
	//
	      << std::setw(Model::log_precision + 7) << CONVERT_DD(relaxation_projection[l]);
    
      for(int w = 0; w < well_size(); ++w) {
	//
	dtemp = eigen_pop(l, w) * well(w).weight_sqrt / nfac;

	if(dtemp < -.01 || dtemp > .01) {
	  //
	  IO::log << std::setw(7) << std::ceil(CONVERT_DD(dtemp) * 100. - .5) / 100.;
	}
	else
	  //
	  IO::log << std::setw(7) << 0;
      }
      
      IO::log << "\n";
    }
    /*    
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << "*Z"
      //
	    << std::setw(Model::log_precision + 7) << "---"
      //
	    << std::setw(Model::log_precision + 7) << "---";

    for(int w = 0; w < well_size(); ++w)
      //
      if(!w || dtemp < well(w).weight_sqrt)
	//
	dtemp = well(w).weight_sqrt;
  
  
    for(int w = 0; w < well_size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << well(w).weight_sqrt / dtemp;
  
    IO::log << "\n";
    */
    IO::log << IO::log_offset
      //
	    << "*R - eigenvalue over the relaxation limit\n"
      //
	    << IO::log_offset
      //
	    << "*P - eigenvector projection squared on the relaxation subspace (=1-F_ne)\n\n";
      //
    //<< IO::log_offset
      //
    //<< "*Z - well partition function square root (normalized)\n\n";
  }
  
  /****************************************************************************************************
   ********************************* CHEMICAL SUBSPACE DIMENSION **************************************
   ****************************************************************************************************/
  
  // number of low eigenvalues
  //
  if(chemical_tolerance > 0.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(eigenval[itemp] > chemical_tolerance * relax_eval_min)
	//
	break;
    //
      else
	//
	eigenval[itemp] = 0.;
  }
  else
    //
    itemp = 0;
  
  if(!bimolecular_size() && itemp > 1 || bimolecular_size() && itemp) {
    //
    if(!(flags & NOPRINT)) {
      //
      if(itemp == 1) {
	//
	IO::log << IO::log_offset << "there is one low eigenvalue\n\n";
      }
      else
	//
	IO::log << IO::log_offset << "there are " << itemp << " low eigenvalues\n\n";
    }
  }
  else
    //
    itemp = 0;
  
  const int low_size = itemp;

  // number of chemical eigenvalues
  //
  if(chemical_threshold > 1.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(relaxation_projection[itemp] * chemical_threshold > 1.)
	//
	break;
  }
  // absolute eigenvalue threshold
  //
  else if(chemical_threshold < 1. && chemical_threshold > 0.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(eigenval[itemp] > chemical_threshold * relax_eval_min)
	//
	break;
  }
  // relative eigenvalue threshold
  //
  else if(chemical_threshold < 0. && chemical_threshold > -1.) {
    //
    for(itemp = well_size(); itemp > 0; --itemp)
      //
      if(eigenval[itemp] <= 0. || eigenval[itemp - 1] / eigenval[itemp] < -chemical_threshold)
	//
	break;
  }
  //
  else {
    //
    std::cerr << funame << "chemical threshold has not been initialized properly: " << chemical_threshold << "\n";
    //
    throw Error::Logic();
  }
  
  const int chem_size = itemp;

  if(chem_size < low_size) {
    //
    std::cerr << funame << "the number of chemical eigenvalues cannot be less than the number of low eigenvalues\n";

    throw Error::Logic();
  }

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n\n";

  // well partitioning infrastructure initialization
  //
  LAPACK::Matrix pop_chem;
    
  partition_t well_partition;
  
  group_t bimolecular_group;

  /*****************************************************************************************************
   *************************************** LOW EIGENVALUE REGIME ***************************************
   *****************************************************************************************************/
  
  if(!(flags & NOPRINT)) {
    //
    std::set<int> g1, g2;
    
    switch(control_barrier.first) {
      //
    case Model::INNER:
      //
      IO::log << IO::log_offset << "rate-limiting barrier " << Model::inner_barrier(control_barrier.second).short_name() << ": "
	//
	      <<  rg1  << "<--"  << Model::inner_barrier(control_barrier.second).short_name() << "-->" << rg2   << "\n";
    
      IO::log << IO::log_offset << "reactive groups: " << rg1.well_set << ", " << rg2.well_set << "\n";

      IO::log << IO::log_offset << "barriers, directly connecting reactive groups:";

      for(std::set<int>::const_iterator bi = graph.inner_set.begin(); bi != graph.inner_set.end(); ++bi)
	//
	if(rg1.well_set.find(Model::inner_connect(*bi).first)  != rg1.well_set.end() &&
	   //
	   rg2.well_set.find(Model::inner_connect(*bi).second) != rg2.well_set.end() ||
	   //
	   rg1.well_set.find(Model::inner_connect(*bi).second) != rg1.well_set.end() &&
	   //
	   rg2.well_set.find(Model::inner_connect(*bi).first) != rg2.well_set.end())
	  //
	  IO::log << "  " << Model::well(Model::inner_connect(*bi).first).short_name() << "<--"
	    //
		  << Model::inner_barrier(*bi).short_name() << "-->"
	    //
		  << Model::well(Model::inner_connect(*bi).second).short_name();
      
      IO::log << "\n";
      
      for(std::set<int>::const_iterator w = rg1.well_set.begin(); w != rg1.well_set.end(); ++w)
	//
	if(Model::well_exclude_group.find(Model::well(*w).name()) == Model::well_exclude_group.end())
	  //
	  g1.insert(*w);
      
      for(std::set<int>::const_iterator w = rg2.well_set.begin(); w != rg2.well_set.end(); ++w)
	//
	if(Model::well_exclude_group.find(Model::well(*w).name()) == Model::well_exclude_group.end())
	  //
	  g2.insert(*w);

      if(g1.size() && g2.size()) {
	//
	IO::log << IO::log_offset << "kinetically significant reactive groups: " << g1 << ", " << g2 << "\n";
      }
      else if(g1.size()) {
	//
	IO::log << IO::log_offset << "kinetically significant reactive group: " << g1 << "\n";
      }
      else if(g2.size()) {
	//
	IO::log << IO::log_offset << "kinetically significant reactive group: " << g2 << "\n";
      }
      
      IO::log << "\n";
      
      break;
      //
    case Model::OUTER:
      //
      IO::log << IO::log_offset << "rate-limiting barrier " <<  Model::outer_barrier(control_barrier.second).short_name() << ": "
	//
	      << rg1  << "<--" << Model::outer_barrier(control_barrier.second).short_name() << "-->"
	//
	      << Model::bimolecular(Model::outer_connect(control_barrier.second).second).short_name() << "\n";
      
      IO::log << IO::log_offset << "reactive groups: " << rg1.well_set << ", "
	//
	      << Model::bimolecular(Model::outer_connect(control_barrier.second).second).short_name()
	//
	      << "\n";
      
      IO::log << IO::log_offset << "barriers, directly connecting reactive groups:";

      for(std::set<int>::const_iterator bi = graph.outer_set.begin(); bi != graph.outer_set.end(); ++bi)
	//
	if(rg1.well_set.find(Model::outer_connect(*bi).first)  != rg1.well_set.end() &&
	   //
	   Model::outer_connect(*bi).second == Model::outer_connect(control_barrier.second).second)
	  //
	  IO::log << "  " << Model::well(Model::outer_connect(*bi).first).short_name() << "<--"
	    //
		  << Model::outer_barrier(*bi).short_name() << "-->"
	    //
		  << Model::bimolecular(Model::outer_connect(*bi).second).short_name();

      IO::log << "\n";
      
      for(std::set<int>::const_iterator w = rg1.well_set.begin(); w != rg1.well_set.end(); ++w)
	//
	if(Model::well_exclude_group.find(Model::well(*w).name()) == Model::well_exclude_group.end())
	  //
	  g1.insert(*w);
      
      if(g1.size())
	//
	IO::log << IO::log_offset << "kinetically significant reactive group: " << g1 << "\n";

      IO::log << "\n";
    }
  }

  //low eigenvalue regime
  //
  if(low_size) {
    //
#ifdef DEBUG
    //
    IO::Marker marker("low eigenvalue regime");
    //
#else
    //
    IO::Marker marker;

    if(!(flags & NOPRINT)) {
      //
      marker.init("low eigenvalue regime");
    }
    //
#endif
  
    pop_chem.resize(well_size(), low_size);
    
    for(int l = 0; l < low_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

    // low eigenvalue wells partition
    //
    if(low_size == well_size()) {
      //
      well_partition.resize(well_size());

      bimolecular_group.clear();
      
      for(int w = 0; w < well_size(); ++w)
	//
	well_partition[w].insert(w);

      dtemp = (float_t)chem_size - projection(well_partition, pop_chem);
    }
    else {
      //
      // partitioning wells into equilibrated groups
      //
      dtemp = threshold_well_partition(pop_chem, well_partition, bimolecular_group);

    }
    
    std::vector<std::set<int> > owp = index_to_well(well_partition);

    std::set<int>               owg = index_to_well(bimolecular_group);

    if(!(flags & NOPRINT)) {
      //
      IO::log << IO::log_offset << "partition projection error = " << dtemp << "\n";
    
      if(low_size != well_size())
	//
	IO::log << IO::log_offset << "well partition: "    << owp << "\n";

      if(owg.size())
	//
	IO::log << IO::log_offset << "bimolecular group: " << owg << "\n";
    }
    
    rate_data.bimolecular_group.clear();
    
    rate_data.well_partition.clear();

    rate_data.well_partition.resize(well_size());
    
    for(int i = 0; i < well_size(); ++i)
      //
      rate_data.well_partition[i].insert(index_to_well(i));
	
    rate_data._bim_index_map = _bim_index_map;

    rate_data.ww_rate.resize(well_size());

    rate_data.ww_rate = 0.;

    if(bimolecular_size()) {
      //
      rate_data.bb_rate.resize(bimolecular_size());
	
      rate_data.bb_rate = 0.;

      rate_data.wb_rate.resize(well_size(), bimolecular_size());

      rate_data.wb_rate = 0.;
    }

    // rate-limiting barrier 
    //
    graph.assert(control_barrier);
    
    const int& type = control_barrier.first;

    const int& b    = control_barrier.second;

    double flux = -1.;
    
    // inner barrier
    //
    if(type == Model::INNER) {
      //
      Model::ground(rg1.well_set, &itemp);
      
      const int w1 = itemp;
	
      Model::ground(rg2.well_set, &itemp);
      
      const int w2 = itemp;

      if(!(flags & NOPRINT))
	//
	IO::log << IO::log_offset << "rate-limiting barrier " << Model::inner_barrier(b).short_name()
	  //
		<< ": " <<  rg1 << "(" << Model::well(w1).short_name() << ")" << "<--"
	  //
		<< Model::inner_barrier(b).short_name()
	  //
		<< "-->" << rg2 << "(" << Model::well(w2).short_name() << ")" << "\n";
    
      // bound groups associated with the rate-limiting barrier
      //
      int g1 = -1, g2 = -1;
	
      for(int g = 0; g < owp.size(); ++g) {
	//
	if(owp[g].find(w1) != owp[g].end())
	  //
	  g1 = g;

	if(owp[g].find(w2) != owp[g].end())
	  //
	  g2 = g;
      }

      // controll barrier belongs to bimolecular group
      //
      if(g1 < 0 && g2 < 0) {
	//
	if(!(flags & NOPRINT))
	  //
	  IO::log << IO::log_offset << "rate-limiting barrier is an inner barrier of the bimolecular group: " << owg << "\n"
	    //
		  << IO::log_offset << "no low eigenvalue treatment is necessary\n";
      }
      // controll barrier is internal barrier of the group
      //
      else if(g1 == g2) {
	//
	if(!(flags & NOPRINT))
	  //
	  IO::log << IO::log_offset << "rate-limiting barrier is an inner barrier of the bound group: " << owp[g1] << "\n"
	    //
		  << IO::log_offset << "no low eigenvalue treatment is necessary\n";
      }
      // control barrier is actually rate-limiting barrier, which connects two different groups
      //
      else {
	//
	if(!(flags & NOPRINT)) {
	  //
	  
	  double ener = Model::inner_barrier(b).thermal_energy(temperature);

	  if(g1 >= 0 && g2 >= 0) {
	    //
	    IO::log << IO::log_offset << "rate-limiting barrier connects two bound groups: " << owp[g1]
	      //
		    << " and " << owp[g2] << "\n";

	    dtemp = std::min(Model::state_density(owp[g1], ener), Model::state_density(owp[g2], ener));

	    IO::log << IO::log_offset << "the microscopic rate constant at the distribution energy maximum over relaxation limit = "
	      //
		    << Model::inner_barrier(b).states(ener) / dtemp / 2. / M_PI / relax_eval_min
	      //
		    << "\n";
	  }
	  else if(g1 >= 0) {
	    //
	    IO::log << IO::log_offset << "rate-limiting barrier connects the bound group " << owp[g1]
	      
		    << " and the bimolecular group " << owg << "\n";

	    IO::log << IO::log_offset << "the microscopic rate constant at the distribution energy maximum over relaxation limit = "
	      //
		    << Model::inner_barrier(b).states(ener) / Model::state_density(owp[g1], ener) / 2. / M_PI / relax_eval_min
	      //
		    << "\n";
	  }
	  else {
	    //
	    IO::log << IO::log_offset << "rate-limiting barrier connects the bound group " << owp[g2]
	      
		    << " and the bimolecular group " << owg << "\n";

	    IO::log << IO::log_offset << "the microscopic rate at the distribution energy maximum over relaxation limit = "
	      //
		    << Model::inner_barrier(b).states(ener) / Model::state_density(owp[g2], ener) / 2. / M_PI / relax_eval_min
	      //
		    << "\n";
	  }
	}
	
	flux = Model::inner_barrier(b).weight(temperature) / std::exp(Model::inner_barrier(b).ground() / temperature)
	  //
	  * temperature / 2. / M_PI;
	  
	// branching fractions
	//
	std::map<std::pair<int, int>, double> bf = branching_fraction(control_barrier, flags);

	if(!(flags & NOPRINT)) {
	  //
	  IO::log << IO::log_offset << "branching fractions:";

	  for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
	    //
	    IO::log << "  " << i->second << "/" << i->first;

	  IO::log << "\n";
	}
	
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i) {
	  //
	  for(std::map<std::pair<int, int>, double>::const_iterator j = bf.begin(); j != bf.end(); ++j) {
	    //
	    if(i == j)
	      //
	      continue;
	    
	    if(i->first.first == Model::BIMOLECULAR && j->first.first == Model::BIMOLECULAR) {
	      //
	      rate_data.bb_rate(bim_to_index(i->first.second), bim_to_index(j->first.second)) = flux * i->second * j->second;
	    }
	    else if(i->first.first == Model::WELL && j->first.first == Model::BIMOLECULAR) {
	      //
	      rate_data.wb_rate(well_to_index(i->first.second), bim_to_index(j->first.second)) = flux * i->second * j->second;
	    }
	    else if(i->first.first == Model::WELL && j->first.first == Model::WELL) {
	      //
	      rate_data.ww_rate(well_to_index(i->first.second), well_to_index(j->first.second)) = flux * i->second * j->second;
	    }
	  }
	}//
	//
      }// flux through the rate-limiting barrier
      //
    }// inner barrier
    //
    // outer barrier
    //
    else if(type == Model::OUTER) {
      //
      Model::ground(rg1.well_set, &itemp);

      const int w1 = itemp;

      if(!(flags & NOPRINT))
	//
	IO::log << IO::log_offset << "rate-limiting barrier " <<  Model::outer_barrier(b).short_name() << ": "
	  //
		<< rg1 << "(" << Model::well(w1).short_name() << ")" << "<--" << Model::outer_barrier(b).short_name() << "-->"
	  //
		<< Model::bimolecular(Model::outer_connect(b).second).short_name() << "\n";
      
      int g;
	
      for(g = 0; g < owp.size(); ++g)
	//
	if(owp[g].find(w1) != owp[g].end())
	  //
	  break;

      // barrier connects to bimolecular group
      //
      if(g == owp.size()) {
	//
	if(!(flags & NOPRINT))
	  //
	  IO::log << IO::log_offset << "rate-limiting outer barrier is conected to the bimolecular group: " << owg << "\n"
	    //
		  << IO::log_offset << "no low igenvalue treatment is necessary\n";
      }
      // barrier connects to the bound group
      //
      else {
	//
	if(!(flags & NOPRINT)) {
	  //
	  double ener = Model::outer_barrier(b).thermal_energy(temperature);

	  IO::log << IO::log_offset << "rate-limiting barrier is connected to the bound group: " << owp[g] << "\n";

	  IO::log << IO::log_offset << "the microscopic rate constant at the distribution energy maximum over relaxation limit = "
	    //
		  << Model::outer_barrier(b).states(ener) / Model::state_density(owp[g], ener) / 2. / M_PI / relax_eval_min
	    //
		  << "\n";
	}
	
	flux = Model::outer_barrier(b).weight(temperature) / std::exp(Model::outer_barrier(b).ground() / temperature)
	  //
	  * temperature / 2. / M_PI;
	  
	// branching fractions
	//
	std::map<std::pair<int, int>, double> bf = branching_fraction(control_barrier, flags);

	if(!(flags & NOPRINT)) {
	  //
	  IO::log << IO::log_offset << "branching fractions:";

	  for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
	    //
	    IO::log << "  " << i->second << "/" << i->first;

	  IO::log << "\n";
	}
	
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i) {
	  //
	  if(i->first.first == Model::BIMOLECULAR) {
	    //
	    rate_data.bb_rate(bim_to_index(i->first.second),  bim_to_index(Model::outer_connect(b).second)) = flux * i->second;
	  }
	  else {
	    //
	    rate_data.wb_rate(well_to_index(i->first.second), bim_to_index(Model::outer_connect(b).second)) = flux * i->second;
	  }//
	  //
	}//
	//
      }// bound group
      //
    }// outer barrier

    if(flux > 0.) {
      //
      if(bimolecular_size())
	//
	rate_data.bw_rate = rate_data.wb_rate.transpose();

#ifdef DEBUG
      //
      if(!(flags & NOPRINT)) {
	//
	IO::log << IO::log_offset << "final well partition:    " << rate_data.well_partition << "\n";

	IO::log << IO::log_offset << "final bimolecular group: " << rate_data.bimolecular_group << "\n";

	IO::log << IO::log_offset << "bimolecular-to-index map:";

	for(std::map<int, int>::const_iterator p = rate_data._bim_index_map.begin(); p != rate_data._bim_index_map.end(); ++p)
	  //
	  IO::log << "  " << Model::bimolecular(p->first).short_name() << "/" << p->second;

	IO::log << "\n";
	
	IO::log << IO::log_offset << "ww_rate matrix:";
	
	for(int i = 0; i < rate_data.ww_rate.size(); ++i) {
	  //
	  IO::log << " {";

	  for(int j = 0; j < rate_data.ww_rate.size(); ++j) {
	    //
	    if(j)
	      //
	      IO::log << ", ";

	    IO::log << rate_data.ww_rate(i, j);
	  }

	  IO::log << "}";
	}

	IO::log << "\n";

	if(bimolecular_size()) {
	  //
	  IO::log << IO::log_offset << "wb_rate matrix:";
	
	  for(int i = 0; i < rate_data.wb_rate.size1(); ++i) {
	    //
	    IO::log << " {";

	    for(int j = 0; j < rate_data.wb_rate.size2(); ++j) {
	      //
	      if(j)
		//
		IO::log << ", ";

	      IO::log << rate_data.wb_rate(i, j);
	    }

	    IO::log << "}";
	  }

	  IO::log << "\n";
	
	  IO::log << IO::log_offset << "bb_rate matrix:";
	
	  for(int i = 0; i < rate_data.bb_rate.size(); ++i) {
	    //
	    IO::log << " {";

	    for(int j = 0; j < rate_data.bb_rate.size(); ++j) {
	      //
	      if(j)
		//
		IO::log << ", ";

	      IO::log << rate_data.bb_rate(i, j);
	    }

	    IO::log << "}";
	  }

	  IO::log << "\n";
	}
      }
      
#endif
      
      return rate_data;
    }
    //
  }// low eigenvalue regime
  
  well_partition.clear();

  bimolecular_group.clear();
  
  /*********************************************************************
   ***************** EIGENSTATES-BIMOLECULAR COUPLING  *****************
   *********************************************************************/

  // eigenvectors projection on the bimolecular subspace;
  //
  LAPACK::Matrix eigen_bim;

  if(bimolecular_size()) {
    //
    eigen_bim.resize(global_size, bimolecular_size());

    eigen_bim = 0.;

#pragma omp parallel for default(shared) private(dtemp, itemp) schedule(static)
  
    for(int l = 0; l < global_size; ++l) {
      //
      for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	//
	const int w = bit->first.first;
      
	const int p = bit->first.second;

	for(int e = 0; e < bit->second.size(); ++e) {
	  //
	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "eigen_bim: well is not in the kinetic basis space\n";

	    throw Error::Logic();
	  }

	  dtemp = 0.;
	
	  for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	    //
	    dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);
	  }

	  dtemp *= bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e];

	  eigen_bim(l, p) += dtemp / 2. / M_PI;
	}
      }
    }

    // high eigenstate correction
    //
    for(int l = 0; l < chem_size; ++l) {
      //
      for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	//
	const int w = bit->first.first;
	
	const int p  = bit->first.second;

	float_t dd_val = 0.;

#pragma omp parallel for default(shared) private(dtemp, itemp) reduction(+: dd_val) schedule(dynamic)
	
	for(int e = 0; e < bit->second.size(); ++e) {
	  //
	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);
		  
	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "high eigenstate eigen_bim correction: well w is not in the e map\n";

	    throw Error::Logic();
	  }

	  for(int w1 = 0; w1 < well_size(); ++w1) {
	    //
	    std::map<int, int>::const_iterator i1 = kinetic_basis[e].well_index_map.find(w1);
		  
	    if(i1 == kinetic_basis[e].well_index_map.end())
	      //
	      continue;
	    
	    float_t proj = 0.;
	      
	    for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
	      //
	      proj += kinetic_basis[e].eigenvector(i1->second, k)
		//
		* kinetic_basis[e].eigenvector(i->second, k)
		//
		/ kinetic_basis[e].eigenvalue[k];

	    proj *= bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e];

	    itemp = e + well(w1).kernel_bandwidth;

	    const int e1_max = itemp < well(w1).size() ? itemp : well(w1).size();
	  
	    itemp = e - well(w1).kernel_bandwidth + 1;

	    const int e1_min = itemp > 0 ? itemp : 0;
	    
	    for(int e1 = e1_min; e1 < e1_max; ++e1) {
	      
	      std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);
		  
	      if(i11 == kinetic_basis[e1].well_index_map.end()) {
		//
		std::cerr << funame << "high eigenstate eigen_bim correction: well w1 is not in the e1 map\n";

		throw Error::Logic();
	      }

	      dtemp = 0.;
	
	      for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
		//
		dtemp += global_eigen(well_shift[e1] + k, l) * kinetic_basis[e1].eigenvector(i11->second, k);
	 
	      dtemp *= proj * well(w1).kernel(e1, e) * well(w1).boltzman_sqrt[e1] / well(w1).boltzman_sqrt[e];
	    
	      dd_val -= dtemp;
	      //
	    }// e1 cycle 
	    //
	  } // w1 cycle
	    //
	}// e cycle

	dd_val *=  pressure / 2. / M_PI;
	
	eigen_bim(l, p) += dd_val;
	//
      }// outer barrier cycle
      //
    }// eigenstate cycle
    //
  }// bimolecular-eigenstate coupling
  
  /***********************************************************************************************
   ************************************* CHEMICAL SUBSPACE ***************************************
   ***********************************************************************************************/

  if(!chem_size) {
    //
    // no kinetically distinct species
    //
    group_t  owg;

    for(int i = 0; i < well_size(); ++i)
      //
      owg.insert(index_to_well(i));

    rate_data.bimolecular_group = owg;
    //
  }// no bound groups
  //
  // bound groups are available
  //
  else {
    //
    // projection of the chemical eigenvectors onto the thermal subspace
    //
    pop_chem.resize(well_size(), chem_size);
    
    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

    if(chem_size == well_size()) {
      //
      well_partition.resize(well_size());

      for(int w = 0; w < well_size(); ++w)
	//
	well_partition[w].insert(w);

      dtemp = (float_t)chem_size - projection(well_partition, pop_chem);

      if(!(flags & NOPRINT))
	//
	IO::log << IO::log_offset << "partition projection error = " << dtemp << "\n";
      
      // well partition in original well representation
      //
      partition_t owp = index_to_well(well_partition);

      rate_data.well_partition = owp;

      for(int g = 0; g < well_size(); ++g)
	//
	rate_data.bf[owp[g]][*owp[g].begin()] = 1.;
    }
    // partitioning wells into equilibrated groups
    //
    else  {
      //
      dtemp = threshold_well_partition(pop_chem, well_partition, bimolecular_group);
      
      // convert chemical eigenvectors in the new basis
      //
      pop_chem = basis(well_partition).transpose() * pop_chem;

      // well partition in original well representation
      //
      partition_t owp = index_to_well(well_partition);

      group_t     owg = index_to_well(bimolecular_group);

      if(!(flags & NOPRINT)) {
	//
	IO::log << IO::log_offset << "partition projection error = " << dtemp << "\n";

	if(chem_size != well_size())
	  //
	  IO::log << IO::log_offset << "well partition: " << owp << "\n";

	if(bimolecular_group.size())
	  //
	  IO::log << IO::log_offset << "bimolecular group: " << owg << "\n";
      }
    
      rate_data.well_partition = owp;

      rate_data.bimolecular_group = owg;

      for(int g = 0; g < chem_size; ++g) {
	//
	std::set<int> wg;
      
	for(std::set<int>::const_iterator w = owp[g].begin(); w != owp[g].end(); ++w)
	  //
	  if(well_extension >= 0. && Model::well_exclude_group.find(Model::well(*w).name()) == Model::well_exclude_group.end())
	    //
	    wg.insert(*w);

	if(wg.size() > 1 && !(flags & NOPRINT))
	  //
	  IO::log << IO::log_offset << "WARNING: " << wg <<" reactants belong to the same partition group\n";
      }
	
      std::map<std::pair<int, int>, double> bf1, bf2;

      double ener;

      float_t stat;
	
      switch(control_barrier.first) {
	//
      case Model::INNER:
	//
	// branching fractions
	//
	if(rg1.well_set.size() == 1 && !rg1.outer_set.size()) {
	  //
	  bf1[std::make_pair((int)Model::WELL, *rg1.well_set.begin())] = 1.;
	}
	else {
	  //
	  ReactiveComplex rc(temperature, pressure, energy_reference, energy_cutoff, rg1, well_extension_cap, flags | NOPRINT);

	  bf1 = rc.branching_fraction(control_barrier, flags);
	}

	if(rg2.well_set.size() == 1 && !rg2.outer_set.size()) {
	  //
	  bf2[std::make_pair((int)Model::WELL, *rg2.well_set.begin())] = 1.;
	}
	else {
	  //
	  ReactiveComplex rc(temperature, pressure, energy_reference, energy_cutoff, rg2, well_extension_cap, flags | NOPRINT);

	  bf2 = rc.branching_fraction(control_barrier, flags);
	}

	// branching fractions output
	//
	if(!(flags & NOPRINT)) {
	  //
	  IO::log << IO::log_offset << "branching fractions for the " << rg1 << " reactive group:\n" << IO::log_offset;
      
	  for(std::map<std::pair<int, int>, double>::const_iterator i = bf1.begin(); i != bf1.end(); ++i)
	    //
	    IO::log << i->second << "/" << i->first << "  ";

	  IO::log << "\n";

	  IO::log << IO::log_offset << "branching fractions for the " << rg2 << " reactive group:\n" << IO::log_offset;
      
	  for(std::map<std::pair<int, int>, double>::const_iterator i = bf2.begin(); i != bf2.end(); ++i)
	    //
	    IO::log << i->second << "/" << i->first << "  ";

	  IO::log << "\n";
	}

	ener = Model::inner_barrier(control_barrier.second).thermal_energy(temperature);

	stat = Model::inner_barrier(control_barrier.second).states(ener) / 2. / M_PI / relax_eval_min;

	for(int g = 0; g < chem_size; ++g) {
	  //
	  if(owp[g].size() == 1) {
	    //
	    rate_data.bf[owp[g]][*owp[g].begin()] = 1.;

	    continue;
	  }
	
	  std::set<int> wg1, wg2;
	
	  for(std::set<int>::const_iterator w = owp[g].begin(); w != owp[g].end(); ++w) {
	    //
	    if(Model::well_exclude_group.find(Model::well(*w).name()) == Model::well_exclude_group.end()) {
	      //
	      if(rg1.well_set.find(*w) != rg1.well_set.end())
		//
		wg1.insert(*w);
	    
	      if(rg2.well_set.find(*w) != rg2.well_set.end())
		//
		wg2.insert(*w);
	    }
	  }
	
	  // reactants and products are in the same partition group
	  //
	  if(wg1.size() && wg2.size()) {
	    //
	    if(!(flags & NOPRINT))
	      //
	      IO::log << IO::log_offset << "WARNING: well partition group " << owp[g]
		//
		      << " includes both " << wg1 << " reactants and " << wg2 << " products\n";
	  }
	  // partition group correlates with certain reactive group
	  //
	  else if(wg1.size() || wg2.size()) {
	    //
	    if(!(flags & NOPRINT))
	      //
	      IO::log << IO::log_offset << owp[g] << " group: "
		//
		      << "the microscopic rate constant at the distribution energy maximum over relaxation limit = "
		//
		      << stat / Model::state_density(owp[g], ener) << "\n";
	  
	    dtemp = 0.;
	  
	    for(std::set<int>::const_iterator w = owp[g].begin(); w != owp[g].end(); ++w) {
	      //
	      std::map<std::pair<int, int>, double>::const_iterator bfi;

	      if(wg1.size()) {
		//
		bfi = bf1.find(std::make_pair((int)Model::WELL, *w));

		if(bfi == bf1.end()) {
		  //
		  rate_data.bf[owp[g]][*w] = 0.;

		  continue;
		}
	      }
	      else {
		//
		bfi = bf2.find(std::make_pair((int)Model::WELL, *w));
	    
		if(bfi == bf2.end()) {
		  //
		  rate_data.bf[owp[g]][*w] = 0.;

		  continue;
		}
	      }

	      if(bfi->second <= 0.) {
		//
		IO::log << IO::log_offset << "WARNING: negative branching fraction: " << bfi->second << " in the " << owp[g] << " group\n";
	      
		rate_data.bf[owp[g]][*w] = 0.;

		continue;
	      }

	      dtemp += bfi->second;
	    
	      rate_data.bf[owp[g]][*w] = bfi->second;
	    }

	    // normalization
	    //
	    if(dtemp == 0.) {
	      //
	      if(!(flags & NOPRINT))
		//
		IO::log << IO::log_offset << "WARNING: zero branching share for the " << owp[g] << " group\n";
	  
	      rate_data.bf.erase(owp[g]);
	    }
	    else {
	      //
	      if(!(flags & NOPRINT))
		//
		IO::log << IO::log_offset << "total branching share for the " << owp[g] << " group = " << dtemp << "\n";
	  
	      for(std::map<int, double>::iterator w = rate_data.bf[owp[g]].begin(); w !=  rate_data.bf[owp[g]].end(); ++w)
		//
		w->second /= CONVERT_DD(dtemp);

	      if(!(flags & NOPRINT)) {
		//
		IO::log << IO::log_offset << "branching fractions for the " << owp[g] << " partition group:";

		for(std::map<int, double>::iterator w = rate_data.bf[owp[g]].begin(); w !=  rate_data.bf[owp[g]].end(); ++w)
		  //
		  IO::log << "  " << w->second << "/" << Model::well(w->first).short_name();

		IO::log << "\n";
	      }
	    }
	  }
	}
      
	break;
	//
      case Model::OUTER:
	//
	if(rg1.well_set.size() == 1 && !rg1.outer_set.size()) {
	  //
	  bf1[std::make_pair((int)Model::WELL, *rg1.well_set.begin())] = 1.;
	}
	else {
	  //
	  ReactiveComplex rc(temperature, pressure, energy_reference, energy_cutoff, rg1, well_extension_cap, flags | NOPRINT);

	  bf1 = rc.branching_fraction(control_barrier, flags);
	}

	// branching fractions output
	//
	if(!(flags & NOPRINT)) {
	  //
	  IO::log << IO::log_offset << "branching fractions for the " << rg1 << " reactive group:\n" << IO::log_offset;
      
	  for(std::map<std::pair<int, int>, double>::const_iterator i = bf1.begin(); i != bf1.end(); ++i)
	    //
	    IO::log << i->second << "/" << i->first << "  ";
      
	  IO::log << "\n";
	}
      
	ener = Model::outer_barrier(control_barrier.second).thermal_energy(temperature);

	stat = Model::outer_barrier(control_barrier.second).states(ener) / 2. / M_PI / relax_eval_min;

	for(int g = 0; g < chem_size; ++g) {
	  //
	  if(owp[g].size() == 1) {
	    //
	    rate_data.bf[owp[g]][*owp[g].begin()] = 1.;

	    continue;
	  }
	
	  std::set<int> wg1;
	
	  for(std::set<int>::const_iterator w = owp[g].begin(); w != owp[g].end(); ++w)
	    //
	    if(Model::well_exclude_group.find(Model::well(*w).name()) == Model::well_exclude_group.end())
	      //
	      if(rg1.well_set.find(*w) != rg1.well_set.end())
		//
		wg1.insert(*w);

	  if(wg1.size()) {
	    //
	    if(!(flags & NOPRINT))
	      //
	      IO::log << IO::log_offset << owp[g] << " group: "
		//
		      << "the microscopic rate constant at the distribution energy maximum over relaxation limit = "
		//
		      << stat / Model::state_density(owp[g], ener) << "\n";
	  
	    dtemp = 0.;
	  
	    for(std::set<int>::const_iterator w = owp[g].begin(); w != owp[g].end(); ++w) {
	      //
	      std::map<std::pair<int, int>, double>::const_iterator bfi = bf1.find(std::make_pair((int)Model::WELL, *w));;

	      if(bfi == bf1.end()) {
		//
		rate_data.bf[owp[g]][*w] = 0.;

		continue;
	      }
	    
	      if(bfi->second <= 0.) {
		//
		if(!(flags & NOPRINT))
		  //
		  IO::log << IO::log_offset << "WARNING: negative branching fraction: " << bfi->second << " for the "<< owp[g] << " group\n";
	      
		rate_data.bf[owp[g]][*w] = 0.;

		continue;
	      }
	    
	      dtemp += bfi->second;
	      
	      rate_data.bf[owp[g]][*w] = bfi->second;
	    }

	    // normalization
	    //
	    if(dtemp == 0.) {
	      //
	      if(!(flags & NOPRINT))
		//
		IO::log << IO::log_offset << "WARNING: zero branching share for the " << owp[g] << " group\n";
	  
	      rate_data.bf.erase(owp[g]);
	    }
	    else {
	      //
	      if(!(flags & NOPRINT))
		//
		IO::log << IO::log_offset << "total branching share for the " << owp[g] << " group = " << dtemp << "\n";
	  
	      for(std::map<int, double>::iterator w = rate_data.bf[owp[g]].begin(); w !=  rate_data.bf[owp[g]].end(); ++w)
		//
		w->second /= CONVERT_DD(dtemp);

	      if(!(flags & NOPRINT)) {
		//
		IO::log << IO::log_offset << "branching fractions for the " << owp[g] << " group:";

		for(std::map<int, double>::iterator w = rate_data.bf[owp[g]].begin(); w !=  rate_data.bf[owp[g]].end(); ++w)
		  //
		  IO::log << "  " << w->second << "/" << Model::well(w->first).short_name();

		IO::log << "\n";
	      }
	    }
	  }
	}
      }
    }
    
    /***************************************************************************
     ********************** WELL-TO-WELL RATE COEFFICIENTS *********************
     ***************************************************************************/

    LAPACK::Matrix  m_direct = pop_chem;
  
    LAPACK::Matrix m_inverse = m_direct.invert();

    std::vector<float_t>  weight = group_weight(well_partition);
  
    LAPACK::Matrix ww_rate(chem_size);

    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j) {
	//
	dtemp = 0.;
	
	for(int l = 0; l < chem_size; ++l)
	  //
	  dtemp += m_direct(j, l) * m_inverse(l, i) * eigenval[l];

	ww_rate(i, j) = -dtemp * SQRT(weight[i] * weight[j])
	  //
	  * energy_step / std::exp(energy_reference / temperature);
      }

    rate_data.ww_rate.resize(chem_size);

    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j)
	//
	rate_data.ww_rate(i, j) = CONVERT_DD(ww_rate(i, j));

#ifdef DEBUG
      
    IO::log << IO::log_offset << "isomerization rate constants ratios (G - well group index):\n"
      //
	    << IO::log_offset << std::setw(3) << "G\\G";
    
    for(int i = 0; i < chem_size; ++i)
      //
      IO::log << std::setw(Model::log_precision + 7) << i;
    
    IO::log << "\n";
    
    for(int j = 0; j < chem_size; ++j) {
      //
      IO::log << IO::log_offset << std::setw(3) << j;
      
      for(int i = 0; i < chem_size; ++i)
	//
	if(i != j) {
	  //
	  if(ww_rate(j, i) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << ww_rate(i, j) / ww_rate(j, i);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	}
	else
	  //
	  IO::log << std::setw(Model::log_precision + 7) << "1";
      
      IO::log << "\n";
    }
    
#endif
    
    /******************************************************************************************************************
     ********************** WELL-TO-BIMOLECULAR RATE COEFFICIENTS / WELL-TO-WELL BRANCHING RATIOS *********************
     ******************************************************************************************************************/

    if(bimolecular_size()) {
      //
      LAPACK::Matrix wb_rate(chem_size, bimolecular_size()), bw_rate(bimolecular_size(), chem_size);
    
      for(int g = 0; g < chem_size; ++g)
	//
	for(int p = 0; p < bimolecular_size(); ++p) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_inverse(l, g) * eigen_bim(l, p);
	  
	  wb_rate(g, p) = dtemp * SQRT(weight[g]) * energy_step / std::exp(energy_reference / temperature);
	}
      
      for(int p = 0; p < bimolecular_size(); ++p)
	//
	for(int g = 0; g < chem_size; ++g) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_direct(g, l) * eigen_bim(l, p);
	  
	  bw_rate(p, g) = dtemp * SQRT(weight[g]) * energy_step / std::exp(energy_reference / temperature);
	}

      rate_data.wb_rate.resize(chem_size, bimolecular_size());

      rate_data.bw_rate.resize(bimolecular_size(), chem_size);

      for(int g = 0; g < chem_size; ++g)
	//
	for(int p = 0; p < bimolecular_size(); ++p) {
	  //
	  rate_data.wb_rate(g, p) = CONVERT_DD(wb_rate(g, p));

	  rate_data.bw_rate(p, g) = CONVERT_DD(bw_rate(p, g));
	}


#ifdef DEBUG

      IO::log << IO::log_offset << "w->b/b->w ratios:\n";

      for(int p = 0; p < bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(6) << bimolecular(p).short_name();
	
	for(int g = 0; g < chem_size; ++g)
	  //
	  if(bw_rate(p, g) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << wb_rate(g, p) / bw_rate(p, g);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";

	IO::log << "\n";
      }
      
#endif
      
    }// well-to-bimolecular rates
    //
  }// bound species
  
  /*****************************************************************************************
   ************************** BIMOLECULAR-TO-BIMOLECULAR RATES *****************************
   *****************************************************************************************/

  if(bimolecular_size()) {
    //
    const int relax_size = global_size - chem_size;
        
    LAPACK::SymmetricMatrix bb_rate;
  
    bb_rate.resize(bimolecular_size());

    for(int p = 0; p < bimolecular_size(); ++p)
      //
      for(int q = p; q < bimolecular_size(); ++q) {
	//
	dtemp = 0.;
	
	for(int l = chem_size; l < global_size; ++l)
	  //
	  dtemp += eigen_bim(l, p) * eigen_bim(l, q) / eigenval[l];
	  
	bb_rate(p, q) = dtemp;
      }

    // high eigenvalue contribution
    //
    for(int e = 0; e < ener_index_max; ++e) {
      //
      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	//
	LAPACK::Vector vtemp(bimolecular_size(), 0.);
	
	for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	  //
	  const int w = bit->first.first;
	  
	  const int p = bit->first.second;

	  if(e >= bit->second.size())
	    //
	    continue;
	  
	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "size of the well associated with the barrier is smaller than the barrier size\n";

	    throw Error::Logic();
	  }
	  
	  vtemp[p] += kinetic_basis[e].eigenvector(i->second, k)
	    //
	    * bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e] / 2. / M_PI;
	}
	
	for(int p = 0; p < bimolecular_size(); ++p)
	  //
	  for(int q = p; q < bimolecular_size(); ++q)
	    //
	    bb_rate(p, q) += vtemp[p] * vtemp[q] / kinetic_basis[e].eigenvalue[k];
      }
    }
    
    bb_rate *= energy_step / std::exp(energy_reference / temperature);

    rate_data.bb_rate.resize(bimolecular_size());

    for(int p = 0; p < bimolecular_size(); ++p)
      //
      for(int q = p; q < bimolecular_size(); ++q)
	//
	rate_data.bb_rate(p, q) = CONVERT_DD(bb_rate(p, q));
    //
  }// bimolecular-to-bimolecular rates

  return rate_data;
}

void MasterEquation::ReactiveComplex::init (double _t, double _p, double _r, double _c, const Model::ChemGraph& _g, int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::init: ";

#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int     itemp;
  double  dtemp;
  
  if(isinit) {
    //
    std::cerr << funame << "already initialized\n";

    throw Error::Init();
  }

  isinit = true;

  graph            = _g;
  
  temperature      = _t;

  pressure         = _p;

  energy_reference = _r;

  energy_cutoff    = _c;

  energy_step      = energy_step_over_temperature * _t;

  control_barrier = std::make_pair((int)Model::UNKNOWN, -1);
  
  // checking graph connectivity
  //
  if(!_g.is_connected()) {
    //
    std::cerr << funame << "chemical graph not connected\n";
    
    throw Error::Logic();
  }

  std::vector<double> diss_ener;

  if(Model::bimolecular_size())
    //
    diss_ener = Model::dissociation_energy_map(_t);

  // output
  //
  if(!(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "<<<<<<<<<<<\n"
      //
	    << IO::log_offset << "temperature = " << int(_t / Phys_const::kelv)
      //
	    << " K  pressure = " << _p / Phys_const::atm
      //
	    << " atm  energy reference = " << std::ceil(_r / Phys_const::kcal * 10.) / 10.
      //
	    << " kcal/mol  energy cutoff = " << std::ceil(_c / Phys_const::kcal * 10.) / 10. << " kcal/mol\n\n";
  
    IO::log << IO::log_offset << "wells:";

    for(std::set<int>::const_iterator w = _g.well_set.begin(); w != _g.well_set.end(); ++w)
      //
      IO::log << " " << Model::well(*w).short_name();

    IO::log << "\n";

    if(_g.inner_set.size()) {
      //
      IO::log << IO::log_offset << "inner connections:";

      for(std::set<int>::const_iterator b = _g.inner_set.begin(); b != _g.inner_set.end(); ++b)
	//
	IO::log << " " << Model::well(Model::inner_connect(*b).first).short_name() << "<--"
	  //
		<< Model::inner_barrier(*b).short_name() << "-->" << Model::well(Model::inner_connect(*b).second).short_name();

      IO::log << "\n";
    }

    if(_g.outer_set.size()) {
      //
      IO::log << IO::log_offset << "outer connections:";

      for(std::set<int>::const_iterator b = _g.outer_set.begin(); b != _g.outer_set.end(); ++b)
	//
	IO::log << " " << Model::well(Model::outer_connect(*b).first).short_name() << "<--"
	  //
		<< Model::outer_barrier(*b).short_name() << "-->" << Model::bimolecular(Model::outer_connect(*b).second).short_name();

      IO::log << "\n";
    }

    IO::log << "\n";

    for(std::set<int>::const_iterator w = _g.well_set.begin(); w != _g.well_set.end(); ++w) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(*w).short_name()
	//
	      << " well:    thermal energy = "  << std::setw(7) 
	//
	      << std::ceil(Model::well(*w).thermal_energy(_t) / Phys_const::kcal * 10.) / 10. 
	//
	      << " kcal/mol  ground energy = " << std::setw(7)
	//
	      << std::ceil(Model::well(*w).ground() / Phys_const::kcal * 10.) / 10.
	//
	      << " kcal/mol";

      if(Model::bimolecular_size())
	//
	IO::log << "  dissociation energy = "
	  //
		<< std::setw(7) << std::ceil(diss_ener[*w] / Phys_const::kcal * 10.) / 10.
	  //
		<< " kcal/mol";

      IO::log << "\n";
    }
  
    for(std::set<int>::const_iterator b = _g.inner_set.begin(); b != _g.inner_set.end(); ++b) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::inner_barrier(*b).short_name()
	//
	      << " barrier: thermal energy = "  << std::setw(7) 
	//
	      << std::ceil(Model::inner_barrier(*b).thermal_energy(_t) / Phys_const::kcal * 10.) / 10. 
	//
	      << " kcal/mol  ground energy = " << std::setw(7)
	//
	      << std::ceil(Model::inner_barrier(*b).real_ground() / Phys_const::kcal * 10.) / 10.
	//
	      << " kcal/mol\n";
	
    }
  
    for(std::set<int>::const_iterator b = _g.outer_set.begin(); b != _g.outer_set.end(); ++b) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::outer_barrier(*b).short_name()
	//
	      << " barrier: thermal energy = " << std::setw(7) 
	//
	      << std::ceil(Model::outer_barrier(*b).thermal_energy(_t) / Phys_const::kcal * 10.) / 10. 
	//
	      << " kcal/mol  ground energy = " << std::setw(7)
	//
	      << std::ceil(Model::outer_barrier(*b).real_ground() / Phys_const::kcal * 10.) / 10.
	//
	      << " kcal/mol\n";
    }
  
    IO::log << "\n";
  }

  // well(s) initialization
  //
  _index_well_map.resize(_g.well_set.size());

  _well.resize(_g.well_set.size());
  
  int count = 0;
  
  for(std::set<int>::const_iterator w = _g.well_set.begin(); w != _g.well_set.end(); ++w, ++count) {
    //
    _well_index_map[*w] = count;

    _index_well_map[count] = *w;
    
    _well[count].set(_t, _r, _c, Model::well(*w), well_extension_cap, flags);
  }

  // inner barrier(s) initialization
  //
  for(std::set<int>::const_iterator bit = _g.inner_set.begin(); bit != _g.inner_set.end(); ++bit) {
    //
    const int w1 = Model::inner_connect(*bit).first;
    
    const int w2 = Model::inner_connect(*bit).second;
    
    const int i1 = well_to_index(w1);
    
    const int i2 = well_to_index(w2);
    
    std::set<int> ii;

    ii.insert(i1);
    
    ii.insert(i2);

    if(ii.size() != 2) {
      //
      std::cerr << funame << "identical indices: " << i1 << "\n";

      throw Error::Logic();
    }
    
    if(inner_barrier.find(ii) != inner_barrier.end()) {
      //
      std::cerr << funame << "well-well index pair is already in the map: (" << i1 << ", " << i2 << ")\n";
      
      throw Error::Logic();
    }
      
    itemp = well(i1).size() < well(i2).size() ? well(i1).size() : well(i2).size();
    
    inner_barrier[ii].set(_t, _r, itemp, Model::inner_barrier(*bit));
  }

  // outer barrier initialization
  //
  for(std::set<int>::const_iterator bit = _g.outer_set.begin(); bit != _g.outer_set.end(); ++bit) {
    //
    const int w = Model::outer_connect(*bit).first;
    
    const int p = Model::outer_connect(*bit).second;

    std::map<int, int>::const_iterator pit = _bim_index_map.find(p);

    int pi = _index_bim_map.size();
    
    if(pit != _bim_index_map.end()) {
      //
      pi = pit->second;
    }
    else {
      //
      _bim_index_map[p] = pi;

      _index_bim_map.push_back(p);
    }

    const int wi = well_to_index(w);

    if(outer_barrier.find(std::make_pair(wi, pi)) != outer_barrier.end()) {
      //
      std::cerr << funame << "well-bimolecular index pair is already in the map: (" << wi << ", " << pi << ")\n";

      throw Error::Logic();
    }
    
    outer_barrier[std::make_pair(wi, pi)].set(_t, _r, well(wi).size(), Model::outer_barrier(*bit));
  }
}

//  branching fractions for the thermal flux through the barrier
//
std::map<std::pair<int, int>, double> MasterEquation::ReactiveComplex::branching_fraction (std::pair<int, int> b, int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::branching_fraction: ";

  using Model::operator<<;
  
#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int itemp;

  const int& type = b.first;
  
  int w1 = -1, w2 = -1;

  Barrier eb;

  if(type == Model::INNER) {
    //
    std::map<int, int>::const_iterator i1 = _well_index_map.find(Model::inner_connect(b.second).first);

    std::map<int, int>::const_iterator i2 = _well_index_map.find(Model::inner_connect(b.second).second);

    
    if(i1 == _well_index_map.end() && i2 == _well_index_map.end()) {
       //
      std::cerr << funame << b << " barrier is not connected to the reactive complex " << graph << "\n";

      throw Error::Logic();
    }

    if(i1 != _well_index_map.end() && i2 != _well_index_map.end()) {
      //
      w1 = i1->second;

      w2 = i2->second;

      itemp = std::min(well(w1).size(), well(w2).size());
    }
    else if(i1 != _well_index_map.end()) {
      //
      w1 = i1->second;

      itemp = well(w1).size();
    }
    else {
      //
      w1 = i2->second;

      itemp = well(w1).size();
    }

    eb.set(temperature, energy_reference, itemp, Model::inner_barrier(b.second));
  }
  else if(type == Model::OUTER) {
    //
    std::map<int, int>::const_iterator i = _well_index_map.find(Model::outer_connect(b.second).first);

    if(i == _well_index_map.end()) {
       //
      std::cerr << funame << b << " barrier is not connected to the reactive complex " << graph << "\n";

      throw Error::Logic();
    }

    w1 = i->second;
    
    eb.set(temperature, energy_reference, well(w1).size(), Model::outer_barrier(b.second));
  }
  else {
    //
    std::cerr << funame << "unknown barrier type: " << type << "\n";

    throw Error::Logic();
  }

  LAPACK::Vector vtemp(eb.size());
  
  for(int e = 0; e < eb.size(); ++e)
    //
    vtemp[e] = eb.state_number[e] * thermal_factor(e) / eb.weight;
  
  std::map<int, LAPACK::Vector> pop_vector;

  pop_vector[index_to_well(w1)] = vtemp;

  if(w2 >= 0)
    //
    pop_vector[index_to_well(w2)] = vtemp;
  
  return propagate(pop_vector, flags);
}

// branching fraction for closed well system
//
std::map<int, double> MasterEquation::ReactiveComplex::branching_fraction (int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::branching_fraction: ";

  //IO::Marker funame_marker(funame);
  
  if(well_size() < 2) {
    //
    std::cerr << funame << "not enough wells: " << well_size() << "\n";

    throw Error::Logic();
  }

  if(bimolecular_size()) {
    //
    std::cerr << funame << "there are bimolecular channels: " << bimolecular_size() << "\n";

    throw Error::Logic();
  }

  LAPACK::Vector vtemp(1);

  vtemp[0] = 1.;

  std::map<int, LAPACK::Vector> pop_vector;

  pop_vector[index_to_well(0)] = vtemp;
  
  std::map<std::pair<int, int>, double> bf = propagate(pop_vector, flags);

  std::map<int, double> res;

  for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
    //
    res[i->first.second] = i->second;

  return res;
}
  
// time-propagation of the initial population vector 
//
std::map<std::pair<int, int>, double> MasterEquation::ReactiveComplex::propagate (const std::map<int, LAPACK::Vector>& pop_vector, int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::propagate: ";

  if(well_size() == 1 && !bimolecular_size()) {
    //
    std::map<std::pair<int, int>, double> res;

    res[std::make_pair((int)Model::WELL, index_to_well(0))] = 1.;

    return res;
  }
  
#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  float_t   dtemp;
  bool      btemp;
  int       itemp;

  // population normalization factor
  //
  btemp = true;
  
  float_t nfac = 0.;
  
  for(std::map<int, LAPACK::Vector>::const_iterator pit = pop_vector.begin(); pit != pop_vector.end(); ++pit) {
    //
    std::map<int, int>::const_iterator i = _well_index_map.find(pit->first);
    
    if(i == _well_index_map.end()) {
      //
      std::cerr << funame << "well is not part of the reactive complex\n";

      throw Error::Logic();
    }
      
    if(!pit->second.isinit() || !pit->second.size() || pit->second.size() > well(i->second).size()) {
      //
      std::cerr << funame << Model::well(pit->first).name() << " well: population vector size out of range\n";

      throw Error::Logic();
    }

    for(int e = 0; e < pit->second.size(); ++e) {
      //
      dtemp = pit->second[e];
      
      if(dtemp >= 0.) {
	//
	nfac += dtemp;
      }
      else if(btemp) {
	//
#ifdef DEBUG

	std::cerr << funame << "there are negative values in the population vector\n";

	throw Error::Range();

#endif
	if(!(flags & NOPRINT))
	  //
	  IO::log << IO::log_offset << "WARNING: there are negative values in the population vector: assumed to be zeros\n";

	btemp = false;
      }
    }
  }
  
  if(nfac <= 0.) {
    //
    std::cerr << funame << "zero population vector\n";
    
    throw Error::Logic();
  }

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "total population = " << nfac << "\n";
  
  // state vector initialization
  //
  LAPACK::Vector state(global_size, 0.);

  for(std::map<int, LAPACK::Vector>::const_iterator pit = pop_vector.begin(); pit != pop_vector.end(); ++pit) {
    //
#pragma omp parallel for default(shared) schedule(dynamic)
    
    for(int e = 0; e < pit->second.size(); ++e) {
      //
      if( pit->second[e] <= 0.)
	//
	continue;
	
      const int w = well_to_index(pit->first);
    
      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);
      
      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "well is not in the kinetic basis space\n";
	
	throw Error::Logic();
      }

      float_t wf = pit->second[e] / well(w).boltzman_sqrt[e];
	
      // state vector setting
      //
      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	//
	state[well_shift[e] + k] +=  wf * kinetic_basis[e].eigenvector(i->second, k);
    }
  }
  
  // bimolecular escape population initialization
  //
  LAPACK::Vector escape(bimolecular_size(), 0.);

  for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
    //
    const int p = bit->first.second;

    float_t dd_val = 0.;

#pragma omp parallel for default(shared) private(dtemp) reduction(+: dd_val) schedule(dynamic)
    
    for(int e = 0; e < bit->second.size(); ++e) {
      //
      std::map<int, int>::const_iterator i1 = kinetic_basis[e].well_index_map.find(bit->first.first);

      if(i1 == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "escape well is not in the kinetic basis space\n";
	
	throw Error::Logic();
      }

      float_t vfac = bit->second.state_number[e] * thermal_factor(e) / well(i1->first).boltzman_sqrt[e];
	
      for(std::map<int, LAPACK::Vector>::const_iterator pit = pop_vector.begin(); pit != pop_vector.end(); ++pit) {
	//
	if(e >= pit->second.size() || pit->second[e] <= 0.)
	  //
	  continue;
	
	const int w = well_to_index(pit->first);
    
	std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);
      
	if(i == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "well is not in the kinetic basis space\n";
	
	  throw Error::Logic();
	}

	// bimolecular escape initial population setting
	//
	dtemp = 0.;
      
	for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
	  //
	  dtemp += kinetic_basis[e].eigenvector(i->second, k)
	    //
	    * kinetic_basis[e].eigenvector(i1->second, k)
	    //
	    / kinetic_basis[e].eigenvalue[k];
      
	dtemp *=  vfac * pit->second[e] / well(w).boltzman_sqrt[e];
	
        dd_val += dtemp;
      }
    }

    escape[p] += dd_val / 2. / M_PI;
  }

  LAPACK::Vector vtemp;

  int      tout = time_propagation_interval                  / time_propagation_step;
  
  int      tmax = time_propagation_limit                     / time_propagation_step;

  float_t tstep = 1. / well(0).collision_frequency(pressure) * time_propagation_step;

  if(tout > 0) {
    //
    IO::log << IO::log_offset << "branching fractions:\n" << IO::log_offset;

    for(int w = 0; w < well_size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::well(index_to_well(w)).short_name();
  
    for(int p = 0; p < _bim_index_map.size(); ++p)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(index_to_bim(p)).short_name();

    IO::log << std::setw(Model::log_precision + 7) << "total" << "\n";
  }
  
  for(int t = 0; t < tmax; ++t) {
    //
    // log output
    //
    if(tout > 0 && !(t % tout)) {
      //
      IO::log << IO::log_offset;

      vtemp = population(state);

      for(int w = 0; w < well_size(); ++w)
	//
	IO::log << std::setw(Model::log_precision + 7) << std::ceil(CONVERT_DD(vtemp[w])  * 100. - 0.5) / 100.;
  
      for(int p = 0; p < _bim_index_map.size(); ++p)
	//
	IO::log << std::setw(Model::log_precision + 7) << std::ceil(CONVERT_DD(escape[p]) * 100. - 0.5) / 100.;

      dtemp = sum(vtemp) + sum(escape);

      IO::log << std::setw(Model::log_precision + 7) << std::ceil(CONVERT_DD(dtemp) * 100. - 0.5) / 100. << "\n";
    }
    
    // escape population contribution
    //
    vtemp = escape_flux(state);

    vtemp *= tstep;

    escape += vtemp;

    // second order time-propagation
    //
    vtemp = kin_mat * state;

    vtemp *= - tstep / 2.;

    vtemp += state;

    vtemp = kin_mat * vtemp;

    vtemp *= tstep;

    state -= vtemp;
  }
  
  std::map<std::pair<int, int>, double> res;

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "branching fractions:";

  // well capture branching fractions
  //
  vtemp = population(state);

  for(int w = 0; w < well_size(); ++w) {
    //
    itemp = index_to_well(w);
    
    if(!(flags & NOPRINT))
      //
      IO::log << "  " << vtemp[w] << "/" << Model::well(itemp).short_name();
    
    res[std::make_pair((int)Model::WELL, itemp)] = CONVERT_DD(vtemp[w]);
  }
  
  // bimolucular escape branching fractions
  //
  for(int p = 0; p < bimolecular_size(); ++p) {
    //
    itemp = index_to_bim(p);
    
    if(!(flags & NOPRINT))
      //
      IO::log << "  " << escape[p] << "/" << Model::bimolecular(itemp).short_name();
    
    res[std::make_pair((int)Model::BIMOLECULAR, itemp)] = CONVERT_DD(escape[p]);
  }

  if(!(flags & NOPRINT))
    //
    IO::log << "\n";
  
  return res;
}

// escape flux through the outer barrier
//
LAPACK::Vector MasterEquation::ReactiveComplex::escape_flux (LAPACK::Vector state) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::escape_flux: ";

  if(state.size() != global_size) {
    //
    std::cerr << funame << "wrong state vector dimension: " << state.size() << "/" << global_size << "\n";

    throw Error::Logic();
  }

  float_t dtemp;

  int itemp;
  
  LAPACK::Vector res(_bim_index_map.size(), 0.);

  for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
    //
    const int w = bit->first.first;
      
    const int p = bit->first.second;

    float_t dd_val = 0.;
    
#pragma omp parallel for default(shared) private(dtemp) reduction(+: dd_val) schedule(dynamic)

    for(int e = 0; e < bit->second.size(); ++e) {
      //
      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "well is not in the kinetic basis space\n";
	
	throw Error::Logic();
      }

      dtemp = 0.;
	
      for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	//
	dtemp += state[well_shift[e] + k] * kinetic_basis[e].eigenvector(i->second, k);
      }

      dtemp *= bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e];

      dd_val += dtemp;
    }

    dd_val /= 2. * M_PI;
    
    res[p] += dd_val;
  }

  for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
    //
    const int w = bit->first.first;
		
    const int p  = bit->first.second;

    float_t dd_val = 0.;
		
#pragma omp parallel for default(shared) private(dtemp, itemp) reduction(+: dd_val) schedule(dynamic)

    for(int e = 0; e < bit->second.size(); ++e) {

      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);
		  
      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "high eigenstate correction: well w is not in the e map\n";

	throw Error::Logic();
      }
	
      for(int w1 = 0; w1 < well_size(); ++w1) {
	//
	std::map<int, int>::const_iterator i1 = kinetic_basis[e].well_index_map.find(w1);

	if(i1 == kinetic_basis[e].well_index_map.end())
	  //
	  continue;

	float_t proj = 0.;

	for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
	  //
	  proj += kinetic_basis[e].eigenvector(i1->second, k)
	    //
	    * kinetic_basis[e].eigenvector(i->second, k)
	    //
	    / kinetic_basis[e].eigenvalue[k];

	proj *=  bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e];

	itemp = e + well(w1).kernel_bandwidth;

	const int e1_max = itemp < well(w1).size() ? itemp : well(w1).size();
	  
	itemp = e - well(w1).kernel_bandwidth + 1;

	const int e1_min = itemp > 0 ? itemp : 0;
	    
	for(int e1 = e1_min; e1 < e1_max; ++e1) {
	  //
	  std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);

	  if(i11 == kinetic_basis[e1].well_index_map.end()) {
	    //
	    std::cerr << funame << "high eigenstate correction: well w1 is not in the e1 map\n";

	    throw Error::Logic();
	  }
	      
	  dtemp = 0.;
      
	  for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
	    //
	    dtemp += state[well_shift[e1] + k] * kinetic_basis[e1].eigenvector(i11->second, k);
	      
	  dtemp *= proj * well(w1).kernel(e1, e) * well(w1).boltzman_sqrt[e1] / well(w1).boltzman_sqrt[e];
	    
	  dd_val -= dtemp;
	  //
	}// e1 cycle
	//
      }// w1 cycle
      //
    }// e2 cycle

    dd_val *= pressure / 2. / M_PI;
    
    res[p] += dd_val;
    //
  }// outer barrier cycle

  return res;
}

// population
//
LAPACK::Vector MasterEquation::ReactiveComplex::population (LAPACK::Vector state) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::population: ";

  if(state.size() != global_size) {
    //
    std::cerr << funame << "wrong state vector dimension: " << state.size() << "/" << global_size << "\n";

    throw Error::Logic();
  }

  float_t dtemp;
	
  LAPACK::Vector res(well_size(), 0.);

  for(int w = 0; w < well_size(); ++w) {
    
    for(int e = 0; e < well(w).size(); ++e) {
      //
      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	std::cerr << funame << "eigen_pop: well is not in the kinetic basis space\n";

	throw Error::Logic();
      }

      dtemp = 0.;
      
      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	//
	dtemp += state[well_shift[e] + k] * kinetic_basis[e].eigenvector(i->second, k);

      dtemp *= well(w).boltzman_sqrt[e];

      res[w] += dtemp;
    }
  }
  
  return res;
}

// dissociation channel
//
void   MasterEquation::ReactiveComplex::set_dissociation_channel (std::pair<int, int> _d, int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::set_dissociation_channel: ";
  
#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int     itemp;
  double  dtemp;
  bool    btemp;

  if(!isinit) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }
  
  switch(_d.first) {
    //
  case Model::INNER:
    //
    if(_d.second >= 0 || _d.second < Model::inner_barrier_size()) {
      //
      std::map<int, int>::const_iterator i1 = _well_index_map.find(Model::inner_connect(_d.second).first);
      
      std::map<int, int>::const_iterator i2 = _well_index_map.find(Model::inner_connect(_d.second).second);

      if(i1 == _well_index_map.end() && i2 == _well_index_map.end()) {
	//
	std::cerr << funame << "dissociation inner barrier is not connected to the reactive complex\n";

	throw Error::Logic();
      }

      if(i1 != _well_index_map.end() && i2 != _well_index_map.end()) {
	//
	std::cerr << funame << "dissociation inner barrier is an internal barrier of the reactive complex\n";

	throw Error::Logic();
      }

      if(i1 != _well_index_map.end()) {
	//
	dissociation_channel.first = i1->second;

	dissociation_channel.second.set(temperature, energy_reference, well(i1->second).size(), Model::inner_barrier(_d.second));
      }
      else {
	//
	dissociation_channel.first = i2->second;

	dissociation_channel.second.set(temperature, energy_reference, well(i2->second).size(), Model::inner_barrier(_d.second));
      }

      break;
    }
    
    std::cerr << funame << "inner barrier index out of range: " << _d.second << "\n";

    throw Error::Logic();


  case Model::OUTER:
    //
    if(_d.second >= 0 || _d.second < Model::outer_barrier_size()) {
      //
      std::map<int, int>::const_iterator wi = _well_index_map.find(Model::outer_connect(_d.second).first);

      std::map<int, int>::const_iterator pi = _bim_index_map.find(Model::outer_connect(_d.second).second);

      if(wi == _well_index_map.end()) {
	//
	std::cerr << funame << "dissociation outer barrier is not connected to the reactive complex\n";

	throw Error::Logic();
      }

      if(pi != _bim_index_map.end() && outer_barrier.find(std::make_pair(wi->second, pi->second)) != outer_barrier.end()) {
	//
	std::cerr << funame << "dissociation outer barrier is already part of the reactive complex\n";

	throw Error::Logic();
      }

      dissociation_channel.first = wi->second;

      dissociation_channel.second.set(temperature, energy_reference, well(wi->second).size(), Model::outer_barrier(_d.second));

      break;
    }
    
    std::cerr << funame << "outer barrier index out of range: " << _d.second << "\n";

    throw Error::Logic();

  default:
    //
    if(!(flags & NOPRINT))
      //
      IO::log << IO::log_offset  << "WARNING: unknown barrier type: " << _d.first << ": skipping initialization\n";
  }
}

void MasterEquation::ReactiveComplex::set_kinetic_matrix (int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::set_kinetic_matrix: ";
  
#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int     itemp;
  float_t dtemp;
  bool    btemp;
  
  /*********************************************************************
   ********************* KINETICALLY ACTIVE BASIS **********************
   *********************************************************************/

  // maximal energy index
  //
  for(int w = 0; w < well_size(); ++w) {
    //
    if(!w || well(w).size() > itemp)
      //
      itemp = well(w).size();
  }
  
  ener_index_max = itemp;

  kinetic_basis.resize(ener_index_max);

#pragma omp parallel for default(shared) private(itemp, dtemp, btemp) schedule(dynamic)

  for(int e = 0; e < ener_index_max; ++e) {// energy cycle
    //
    // available energy bins
    //
    std::vector<int>  index_well_map;
     
    std::map<int, int> well_index_map;
      
    for(int w = 0; w < well_size(); ++w)
      //
      if(e < well(w).size()) {
	//
	well_index_map[w] = index_well_map.size();
	  
	index_well_map.push_back(w);
      }

    if(!index_well_map.size()) {
      //
      std::cerr << funame << "empty basis\n";

      throw Error::Logic();
    }
    
    kinetic_basis[e].well_index_map = well_index_map;
      
    kinetic_basis[e].index_well_map = index_well_map;

    // microcanonical kinetic matrix
    //
    LAPACK::SymmetricMatrix km(index_well_map.size());
      
    km = 0.;

    // isomerization contribution
    //
    for(inner_t::const_iterator bit = inner_barrier.begin(); bit != inner_barrier.end(); ++bit) {
      //
      if(e >= bit->second.size()) 
	//
	continue;
	
      const int w1 = *bit->first.begin();
	  
      const int w2 = *bit->first.rbegin();
      
      std::map<int, int>::const_iterator i1 = well_index_map.find(w1);
	  
      std::map<int, int>::const_iterator i2 = well_index_map.find(w2);

      if(i1 == well_index_map.end() || i2 == well_index_map.end()) {
	//
	std::cerr << funame << "kinetic basis: well indices are not in the map\n";
	  
	throw Error::Logic();
      }
      
      km(i1->second, i2->second) -= bit->second.state_number[e] / 2. / M_PI
	//
	/ SQRT(well(w1).state_density[e] * well(w2).state_density[e]);
	  
      km(i1->second, i1->second) += bit->second.state_number[e] / 2. / M_PI
	//
	/ well(w1).state_density[e];

      km(i2->second, i2->second) += bit->second.state_number[e] / 2. / M_PI
	//
	/ well(w2).state_density[e];
    }

    // bimolecular contribution
    //
    for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
      //
      if(e >= bit->second.size())
	//
	continue;
	
      const int w = bit->first.first;
	  
      std::map<int, int>::const_iterator i = well_index_map.find(w);

      if(i == well_index_map.end()) {
	//
	std::cerr << funame << "kinetic basis: well index is not in the map\n";

	throw Error::Logic();
      }
	  
      km(i->second, i->second) += bit->second.state_number[e] / 2. / M_PI / well(w).state_density[e];
    }

    // dissociation channel
    //
    if(dissociation_channel.second.size() && e < dissociation_channel.second.state_number.size()) {
      //
      std::map<int, int>::const_iterator i = well_index_map.find(dissociation_channel.first);

      if(i == well_index_map.end()) {
	//
	std::cerr << funame << "dissociation channel: index out of range\n";
	
	throw Error::Logic();
      }
	
      km(i->second, i->second) += dissociation_channel.second.state_number[e] / 2. / M_PI / well(i->first).state_density[e];
    }
      
    LAPACK::Matrix evec(index_well_map.size());

    LAPACK::Vector eval;

    try {
      //
      eval = km.eigenvalues(&evec);
    }
    catch(Error::Lapack) {
      //
      std::cerr << funame << "kinetic basis initialization failed at "
	//
		<< std::ceil((energy_reference - e * energy_step) / Phys_const::kcal * 10.) / 10.
	//
		<< " kcal/mol failed\n";
      throw;
    }
      
    kinetic_basis[e].eigenvalue  = eval;
      
    kinetic_basis[e].eigenvector = evec;
      
    // kinetically active subspace
    //
    for(itemp = 0; itemp < index_well_map.size(); ++itemp) {
      //
      if(eval[itemp] > well(0).collision_frequency(pressure) * reduction_threshold)
	//
	break;
    }
      
    kinetic_basis[e].active_size = itemp;
    //
  }// kinetically active basis

  // global indexing
  //
  well_shift.resize(ener_index_max);
  
  itemp = 0;
  
  for(int e = 0; e < ener_index_max; itemp += kinetic_basis[e++].active_size)
    //
    well_shift[e] = itemp;
  
  global_size = itemp;

  // global kinetic matrix size without reduction
  //
  itemp = 0;
  //
  for(int w = 0; w < well_size(); ++w)
    //
    itemp += well(w).size();

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "original kinetic matrix size = " << itemp << "\n"
      
	    << IO::log_offset << "reduced  kinetic matrix size = " << global_size << "\n";

  /********************************************************************************************* 
   ********************************** GLOBAL KINETIC MATRIX ************************************
   *********************************************************************************************/

  kin_mat.resize(global_size);
  
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
  int emax = ener_index_max * ener_index_max;
  
#pragma omp parallel for default(shared) private(dtemp, itemp) schedule(dynamic)

  for(int e = 0; e < emax; ++e) {
    //
    const int e1 = e / ener_index_max;
    
    const int e2 = e % ener_index_max;
    
    for(int l1 = 0; l1 < kinetic_basis[e1].active_size; ++l1) {
      //
      for(int l2 = 0; l2 < kinetic_basis[e2].active_size; ++l2) {
	//
	if(e1 == e2 && l2 < l1)
	  //
	  continue;

	dtemp = 0.;
	  
	for(int w = 0; w < well_size(); ++w) {
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
	    dtemp +=  well(w).kernel(e1, e2) * pressure
	      //
	      * well(w).boltzman_sqrt[e1] / well(w).boltzman_sqrt[e2]
	      //
	      * kinetic_basis[e1].eigenvector(i1->second, l1)
	      //
	      * kinetic_basis[e2].eigenvector(i2->second, l2);
	  }
	}

	const int j1 = well_shift[e1] + l1;
	  
	const int j2 = well_shift[e2] + l2;
	  
	if(j1 != j2) {
	  //
	  kin_mat(j1, j2) = dtemp;

	  kin_mat(j2, j1) = dtemp;
	}
	else {
	  //
	  kin_mat(j1, j2) += dtemp;
	}
	//
      }// l2 cycle
      //
    }// l1 cycle
    //
  }// e cycle
}

const MasterEquation::ReactiveComplex::Well&  MasterEquation::ReactiveComplex::well (int w) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::well: ";

  if(w >= 0 && w < _well.size())
    //
    return _well[w];

  std::cerr << funame << "well index out of range: " << w << "\n";

  throw Error::Logic();
}

MasterEquation::float_t MasterEquation::ReactiveComplex::thermal_factor(int e)
{
  return EXP((float_t)e * (float_t)energy_step_over_temperature);
}

int MasterEquation::ReactiveComplex::RateData::bim_to_index (int p) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::RateData::bim_to_index: ";

  std::map<int, int>::const_iterator i = _bim_index_map.find(p);
  
  if(i == _bim_index_map.end()) {
    //
    std::cerr << funame << "wrong index\n";

    throw Error::Logic();
  }

  return i->second;
}

int MasterEquation::ReactiveComplex::well_to_index (int w) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::well_to_index: ";
  
  std::map<int, int>::const_iterator wit = _well_index_map.find(w);

  if(wit != _well_index_map.end())
    //
    return wit->second;

  std::cerr << funame << "well index not in the map: " << w << "\n";

  throw Error::Logic();
}

std::set<int> MasterEquation::ReactiveComplex::bim_to_index (const std::set<int>& g) const
{
  std::set<int> res;
  
  for(std::set<int>::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res.insert(bim_to_index(*w));

  return res;
}
  
std::set<int> MasterEquation::ReactiveComplex::well_to_index (const std::set<int>& g) const
{
  std::set<int> res;
  
  for(std::set<int>::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res.insert(well_to_index(*w));

  return res;
}
  
std::vector<std::set<int> > MasterEquation::ReactiveComplex::well_to_index (const std::vector<std::set<int> >& p) const
{
  std::vector<std::set<int> > res(p.size());

  for(int g = 0; g < p.size(); ++g)
    //
    res[g] = well_to_index(p[g]);
  
  return res;
}
  
std::vector<std::set<int> > MasterEquation::ReactiveComplex::bim_to_index (const std::vector<std::set<int> >& p) const
{
  std::vector<std::set<int> > res(p.size());

  for(int g = 0; g < p.size(); ++g)
    //
    res[g] = bim_to_index(p[g]);
  
  return res;
}
  
int MasterEquation::ReactiveComplex::index_to_well (int i) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::index_to_well: ";
  
  if(i >= 0 && i < _index_well_map.size())
    //
    return _index_well_map[i];

  std::cerr << funame << "index out of range: " << i << "\n";

  throw Error::Logic();
}

std::set<int> MasterEquation::ReactiveComplex::index_to_well (const std::set<int>& g) const
{
  std::set<int> res;
  
  for(std::set<int>::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res.insert(index_to_well(*w));

  return res;
}
  
std::set<int> MasterEquation::ReactiveComplex::index_to_bim (const std::set<int>& g) const
{
  std::set<int> res;
  
  for(std::set<int>::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res.insert(index_to_bim(*w));

  return res;
}
  
std::vector<std::set<int> > MasterEquation::ReactiveComplex::index_to_well (const std::vector<std::set<int> >& p) const
{
  std::vector<std::set<int> > res(p.size());

  for(int g = 0; g < p.size(); ++g)
    //
    res[g] = index_to_well(p[g]);
  
  return res;
}
  
std::vector<std::set<int> > MasterEquation::ReactiveComplex::index_to_bim (const std::vector<std::set<int> >& p) const
{
  std::vector<std::set<int> > res(p.size());

  for(int g = 0; g < p.size(); ++g)
    //
    res[g] = index_to_bim(p[g]);
  
  return res;
}
  
int MasterEquation::ReactiveComplex::bim_to_index (int p) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::bim_to_index: ";
  
  std::map<int, int>::const_iterator pit = _bim_index_map.find(p);

  if(pit != _bim_index_map.end())
    //
    return pit->second;

  std::cerr << funame << "bimolecular index not in the map: " << p << "\n";

  throw Error::Logic();
}

int MasterEquation::ReactiveComplex::index_to_bim (int i) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::index_to_bim: ";
  
  if(i >= 0 && i < _index_bim_map.size())
    //
    return _index_bim_map[i];

  std::cerr << funame << "index out of range: " << i << "\n";

  throw Error::Logic();
}

void MasterEquation::ReactiveComplex::Well::set (double temperature, double energy_reference, double energy_cutoff,
						 //
						 const Model::Well& model, const std::map<int, double>& well_extension_cap, int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::Well:set: ";

  double  dtemp;
  int     itemp;
  bool    btemp;

  float_t dd_temp;
  
  if(isinit) {
    //
    std::cerr << funame << "already initialized\n";

    throw Error::Init();
  }

  isinit = true;

  /******************************************************************************************
   *********************************** STATE DENSITY ****************************************
   ******************************************************************************************/

  double energy_step = energy_step_over_temperature * temperature;
  
  // enthalpy (average energy), entropy, and thermal capacity
  //
  double hval, sval, cval;

  model.esc_parameters(temperature, hval, sval, cval);

  double base_ener = model.ground();

  int ext_size = 0;

  std::map<int, double>::const_iterator wxi = well_extension_cap.find(Model::well_by_name(model.name()));

  // well extension not for excluded wells
  //
  if(well_extension >= 0. && !(flags & NO_WELL_EXTENSION) &&
     //
     Model::well_exclude_group.find(model.name()) == Model::well_exclude_group.end()) {
    //
    base_ener += hval;

    // well extension energy correction
    //
    if(wxi != well_extension_cap.end() && base_ener > wxi->second)
      //
      base_ener = wxi->second + well_extension * (base_ener - wxi->second);
    
    dtemp = base_ener - energy_cutoff;
    
    if(dtemp > 0.) {
      //
      ext_size = std::ceil(dtemp / energy_step);

      if(!(flags & NOPRINT))
	//
	IO::log << IO::log_offset << std::setw(4) << model.short_name() << " well: well extension: base energy = "
	  //
		<< std::ceil(base_ener / Phys_const::kcal * 10.) / 10.  << " kcal/mol"

		<< "   extension size = " << ext_size << "\n";
    }
  }
  
  if(base_ener < energy_cutoff) {
    //
    base_ener = energy_cutoff;
  }

  dtemp = energy_reference - base_ener;
  
  if(dtemp <= 0.) {
    //
    std::cerr << funame << model.name() << " well: base energy exceeds energy reference: energy reference = "
      //
	      << std::ceil(energy_reference / Phys_const::kcal * 10.) / 10. << " kcal/mol  base energy = "
      
	      << std::ceil(base_ener / Phys_const::kcal * 10.) / 10. << " kcal/mol\n";

    throw Error::Range();
  }

  const int base_size = std::ceil(dtemp / energy_step);

  state_density.resize(base_size + ext_size);

  double ener = energy_reference;
  
  for(int e = 0; e < base_size; ++e, ener -= energy_step) {
    //
    dtemp = model.states(ener);
    
    if(dtemp <= 0.) {
      //
      if(!e) {
	//
	std::cerr << funame << model.name() << " well: empty well\n";

	throw Error::Range();
      }

      if(!(flags & NOPRINT))
	//
	IO::log << IO::log_offset << std::setw(4) << model.short_name()  << " well: WARNING: nonpositive density at "
	  //
		<< std::ceil(ener / Phys_const::kcal * 10.) / 10. << " kcal/mol: truncating\n";

      state_density.resize(e);      

      break;    
    }
    
    state_density[e] = dtemp;
  }

  for(int e = base_size; e < state_density.size(); ++e)
    //
    state_density[e] = dtemp;

  // Boltzman distribuions
  //
  boltzman.resize(size());
  
  boltzman_sqrt.resize(size());
  
  weight = 0.;
  
  for(int e = 0; e < size(); ++e) {
    //
    dd_temp = state_density[e] * thermal_factor(e);

    boltzman[e] = dd_temp;
    
    boltzman_sqrt[e] = SQRT(dd_temp);

    weight += dd_temp;
  }
  
  weight_sqrt = SQRT(weight);

#ifdef DEBUG
  
  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << std::setw(4) << model.short_name() << " well: size = " << std::setw(4) << size() << "\n";

#endif
  
  /***********************************************************************************************
   **************************** COLLISIONAL ENERGY TRANSFER KERNEL *******************************
   ***********************************************************************************************/
  
  collision_factor = 0.;

  std::vector<double> kernel_fraction(Model::buffer_size());
  
  for(int b = 0; b < Model::buffer_size(); ++b) {
    //
    dtemp = model.collision(b)(temperature) * Model::buffer_fraction(b);
    
    kernel_fraction[b] = dtemp;
    
    collision_factor += dtemp;
  }

  //for(int b = 0; b < Model::buffer_size(); ++b)
    //
    //kernel_fraction[b] /= collision_factor;

  kernel.resize(size());
    
  kernel = 0.;

  kernel_t tmp_kernel(size());

  float_t nfac;

  for(int b = 0; b < Model::buffer_size(); ++b) {

    tmp_kernel = 0.;

    // collisional energy transfer on the grid
    //
    itemp = (int)std::ceil(model.kernel(b).cutoff_energy(temperature) / energy_step);

    if(itemp < 2) {
      //
      std::cerr << funame << "no collisional energy transfer\n";

      throw Error::Range();
    }
    
    std::vector<float_t> energy_transfer(itemp);
	
    for(int i = 0; i < energy_transfer.size(); ++i)
      //
      energy_transfer[i] = model.kernel(b)((double)i * energy_step, temperature);

    std::vector<float_t> aux_state_density(energy_transfer.size());

    for(int e = 1; e < energy_transfer.size(); ++e)
      //
      aux_state_density[e] = model.states(energy_reference + e * energy_step);
    
    if(!b || kernel_bandwidth < energy_transfer.size())
      //
      kernel_bandwidth = energy_transfer.size();

    // kernel escape
    //
    itemp = energy_transfer.size() - 1;

    const int esize = itemp < size() ? itemp : size();
    
    if(kernel_escape.size() < esize)
      //
      kernel_escape.resize(esize);

    for(int e = 1; e < energy_transfer.size(); ++e) {
      //
      nfac = energy_transfer[0];
      
      for(int i = 1; i < energy_transfer.size(); ++i) {
	//
	itemp = i - e;

	if(itemp == size())
	  //
	  break;
	    
	dd_temp = aux_state_density[e] / thermal_factor(i);

	if(itemp < 0) {
	  //
	  dd_temp /= aux_state_density[itemp];
	}
	else
	  //
	  dd_temp /= state_density[itemp];

	nfac += energy_transfer[i] * (1. + dd_temp);
      }
      
      nfac /= kernel_fraction[b];
	  
      for(int i = 0; i < esize; ++i) {
	//
	itemp = i + e;

	if(itemp >= energy_transfer.size())
	  //
	  break;
	
	kernel_escape[i] += energy_transfer[itemp] * aux_state_density[e]
	  //
	  / state_density[i] / thermal_factor(itemp) / nfac;
      }//
      //
    }// kernel escape
  
    for(int e = 0; e < size(); ++e) {
      //
      nfac = energy_transfer[0];
	  
      for(int i = 1; i < energy_transfer.size(); ++i) {
	//
	itemp = e + i;

	if(itemp == size())
	  //
	  break;
	    
	tmp_kernel(e, itemp) = energy_transfer[i];

	tmp_kernel(itemp, e) = energy_transfer[i] * state_density[e] / state_density[itemp] / thermal_factor(i);

	nfac += tmp_kernel(e, itemp) + tmp_kernel(itemp, e);
      }

      nfac /= kernel_fraction[b];
	  
      for(int i = 1; i < energy_transfer.size(); ++i) {
	//
	itemp = e + i;

	if(itemp == size())
	  //
	  break;
	    
	tmp_kernel(e, itemp) /= -nfac;

	tmp_kernel(itemp, e) /= -nfac;
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
      
    kernel += tmp_kernel;
  }
}
  
void MasterEquation::ReactiveComplex::Barrier::set (double temperature, double energy_reference, int size_max,
						    //
						    const Model::Species& model, int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::Barrier:set: ";

  int     itemp;
  float_t dtemp;

  if(isinit) {
    //
    std::cerr << funame << "already initialized\n";

    throw Error::Init();
  }

  isinit = true;

  weight = 0.;
  
  double energy_step = energy_step_over_temperature * temperature;

  dtemp = model.ground();

  if(is_global_cutoff && dtemp < lower_global_cutoff)
    //
    dtemp = lower_global_cutoff;
  
  itemp = std::ceil((energy_reference - (double)dtemp) / energy_step);

  itemp = itemp < size_max ? itemp : size_max;

  if(itemp <= 0)
    //
    return;
  
  state_number.resize(itemp);

  double ener = energy_reference;
  
  for(int e = 0; e < size(); ++e, ener -= energy_step) {
    //
    dtemp = model.states(ener);
    
    if(dtemp <= 0.) {
      //
      state_number.resize(e);
      
      break;
    }
    
    state_number[e] = dtemp;
    
    weight += state_number[e] * thermal_factor(e);
  }
  
  //weight *= energy_step_over_temperature;

#ifdef DEBUG
  
  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << std::setw(4) << model.short_name() << " barrier: size = " << std::setw(4) << size() << "\n";

#endif
}

MasterEquation::float_t MasterEquation::ReactiveComplex::threshold_well_partition (LAPACK::Matrix pop_chem,
								  //
								  partition_t&   well_partition,
								  //
								  group_t&       bimolecular_group, int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::threshold_well_partition: ";

#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int         itemp;
  float_t     dtemp;
  bool        btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  
  if(well_size() != pop_chem.size1()) {
    //
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1()
      //
	      << " differs from the number of wells = " << well_size() << "\n";
    
    throw Error::Logic();
  }

  if(chem_size >= well_size()) {
    //
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";
    
    throw Error::Logic();
  }

  bimolecular_group.clear();

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "well projection threshold  = " << well_projection_threshold << "\n";

  if(!(flags & NOPRINT) && flags & USE_PROJECTION_TOLERANCE)
    //
    IO::log << IO::log_offset << "well projection tolerance  = " << well_projection_tolerance << "\n";
  
  // projection-ordered wells
  //
  std::multimap<float_t, int> proj_well_map;

  std::set<int> exclude_group;
  
  for(int w = 0; w < well_size(); ++w) {
    //
    group_t g;

    g.insert(w);

    dtemp = projection(g, pop_chem);

    if(flags & USE_PROJECTION_TOLERANCE && dtemp < well_projection_tolerance)
      //
      exclude_group.insert(w);
    
    proj_well_map.insert(std::make_pair(dtemp, w));
    
    bimolecular_group.insert(w);
  }

  // large projection wells map
  //
  std::vector<int> prim_well_map;

  std::multimap<float_t, int>::const_reverse_iterator pit;

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "primary (large) well projection/well:";

  for(pit = proj_well_map.rbegin(); pit != proj_well_map.rend(); ++pit) {
    //
    if(pit->first < well_projection_threshold && prim_well_map.size() >= chem_size)
      //
      break;

    if(!(flags & NOPRINT))
      //
      IO::log << "   " << pit->first << "/" << pit->second;
    
    dtemp = pit->first;

    prim_well_map.push_back(pit->second);
    
    bimolecular_group.erase(pit->second);
  }
    
  if(!(flags & NOPRINT))
    //
    IO::log << "\n";

  if(pit != proj_well_map.rend() && !(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "secondary (small) well projection/well:";
    
    for(; pit != proj_well_map.rend(); ++pit)
      //
      IO::log << "   " << pit->first << "/" << pit->second;

    IO::log << "\n";
  }
    
  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "primary wells #/chemical groups #: "
      
	    << prim_well_map.size() << "/" << chem_size << "\n";
  
  if(dtemp < well_projection_threshold && !(flags & NOPRINT))
    //
    IO::log << IO::log_offset << "WARNING: primary well projection is smaller than the well projection threshold: "
      //
	    << dtemp << "\n"
      //
	    << IO::log_offset << "         temperature = " << (int)std::ceil(temperature / Phys_const::kelv)
      //
	    << " K  pressure = " << pressure / Phys_const::atm << " atm\n";
  
  float_t proj = -1.;
  
  for(PartitionGenerator pg(chem_size, prim_well_map.size()); !pg.end(); ++pg) {
    //
    // initialize new partition
    //
    partition_t p = pg.partition(prim_well_map);

    // partition projection
    //
    dtemp = projection(p, pop_chem);

    if(dtemp > proj) {
      //
      proj = dtemp;
      
      well_partition = p;
    }
  }

  while(bimolecular_group.size()) {
    //
    int gi, wi;
    
    proj = -1.;
    
    for(group_t::const_iterator w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
      //
      if(exclude_group.find(*w) != exclude_group.end())
	//
	continue;
      
      for(int s = 0; s < chem_size; ++s) {
	//
	group_t g = well_partition[s];

	g.insert(*w);
	
	dtemp = projection(g, pop_chem) - projection(well_partition[s], pop_chem);
	
	if(dtemp > proj) {
	  //
	  proj = dtemp;
	  
	  gi = s;
	  
	  wi = *w;
	}
      }
    }
    
    if(proj < 0.)
      //
      break;

    well_partition[gi].insert(wi);
    
    bimolecular_group.erase(wi);
  }

#ifdef DEBUG
  
  // partition output
  //
  if(!(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "well partition: ";
      
    for(int g = 0; g < well_partition.size(); ++g) {
      //
      IO::log << "   ";
    
      for(group_t::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w) {
	//
	if(w != well_partition[g].begin())
	  //
	  IO::log << "+";
      
	IO::log << *w;
      }
	
      IO::log << "/" << g;
    }

    IO::log << "\n";

    if(bimolecular_group.size()) {
      //
      IO::log << IO::log_offset << "bimolecular group:";
    
      for(group_t::const_iterator w = bimolecular_group.begin(); w != bimolecular_group.end(); ++w)
	//
	IO::log << "   " << *w;

      IO::log << "\n";
    }
  }

#endif
  
  dtemp = (float_t)chem_size - projection(well_partition, pop_chem);

  return dtemp;
}

MasterEquation::float_t MasterEquation::ReactiveComplex::projection (const partition_t& p, LAPACK::Matrix pop_chem) const 
{
  float_t res = 0.;
  
  for(int g = 0; g < p.size(); ++g)
    
    res += projection(p[g], pop_chem);

  return res;
}
	
MasterEquation::float_t MasterEquation::ReactiveComplex::projection (const group_t& g, LAPACK::Matrix pop_chem) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::projection: ";

  float_t dtemp;

  const int chem_size = pop_chem.size2();
  
  assert(g);
  
  if(well_size() != pop_chem.size1()) {
    //
    std::cerr << funame << "pop_chem: first dimension out of range: " << pop_chem.size1() << "\n";
    
    throw Error::Logic();
  }

  if(g.size() == 1)
    //
    return vdot(pop_chem.row(*g.begin()));

  float_t res = 0.;
  
  for(int l = 0; l < chem_size; ++l) {
    //
    dtemp = 0.;
    
    for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
      //
      dtemp += well(*w).weight_sqrt * pop_chem(*w, l);
    
    res += dtemp * dtemp;
  }

  // normalization
  //
  res /= group_weight(g);

  return res;
}

MasterEquation::float_t MasterEquation::ReactiveComplex::group_weight (const group_t& g) const
{
  assert(g);
  
  float_t res = 0.;
  
  for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res += well(*w).weight;

  return res;
}

std::vector<MasterEquation::float_t> MasterEquation::ReactiveComplex::group_weight (const partition_t& p) const
{
  std::vector<float_t> res(p.size());
  
  for(int g = 0; g < p.size(); ++g)
    //
    res[g] = group_weight(p[g]);

  return res;
}

LAPACK::Matrix MasterEquation::ReactiveComplex::basis (const partition_t& p) const 
{
  LAPACK::Matrix res(well_size(), p.size());
  
  for(int g = 0; g < p.size(); ++g)
    //
    res.column(g) = basis(p[g]);
  
  return res;
}

LAPACK::Vector MasterEquation::ReactiveComplex::basis (const group_t& g) const
{
  assert(g);

  LAPACK::Vector res(well_size(), 0.);
  
  if(g.size() == 1) {
    //
    res[*g.begin()] = 1.;
    
    return res;
  }
    
  float_t dtemp = SQRT(group_weight(g));

  for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res[*w] = well(*w).weight_sqrt / dtemp;
  
  return res;
}

void MasterEquation::ReactiveComplex::assert (const group_t& g) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::assert: ";

  if(!g.size() || *g.begin() < 0 || *g.rbegin() >= well_size()) {
    //
    std::cerr << funame << "group indices out of range\n";

    throw Error::Logic();
  }
}

std::set<int> MasterEquation::ReactiveComplex::RateData::missing_bf (const std::list<reac_t>& reactions) const
{
  const char funame [] = " MasterEquation::ReactiveComplex::RateData::missing_bf: ";

  std::map<int, int> wg_map;

  std::set<int> res;
  
  for(int g = 0; g < well_partition.size(); ++g)
    //
    for(group_t::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w)
      //
      wg_map[*w] = g;

  for(std::list<reac_t>::const_iterator r = reactions.begin(); r != reactions.end(); ++r) {
    //
    if(r->size() != 2) {
      //
      std::cerr << funame << "wrong number of reactants: " << r->size() << "\n";

      throw Error::Logic();
    }

    // well-to-bimolecular
    //
    if(r->begin()->first == Model::WELL && r->rbegin()->first == Model::BIMOLECULAR) {
      //
      const int& w = r->begin()->second;

      std::map<int, int>::const_iterator i = wg_map.find(w);

      if(i != wg_map.end() && bf.find(well_partition[i->second]) == bf.end())
	//
	res.insert(i->second);
    }
    // well-to-bimolecular
    //
    else if(r->begin()->first == Model::BIMOLECULAR && r->rbegin()->first == Model::WELL) {
      //
      const int& w = r->rbegin()->second;

      std::map<int, int>::const_iterator i = wg_map.find(w);

      if(i != wg_map.end()  && bf.find(well_partition[i->second]) == bf.end())
	//
	res.insert(i->second);
    }
    // well-to-well rate
    //
    else if(r->begin()->first == Model::WELL && r->rbegin()->first == Model::WELL) {
      //
      const int& w1 = r->begin()->second;

      const int& w2 = r->rbegin()->second;

      std::map<int, int>::const_iterator i1 = wg_map.find(w1);

      std::map<int, int>::const_iterator i2 = wg_map.find(w2);

      if(i1 != wg_map.end() && 	 i2 != wg_map.end() && i1->second != i2->second) {
	//
	if(bf.find(well_partition[i1->second]) == bf.end())
	  //
	  res.insert(i1->second);

	if(bf.find(well_partition[i2->second]) == bf.end())
	  //
	  res.insert(i2->second);
      }
    }
  }

  return res;
}

std::set<int> MasterEquation::ReactiveComplex::RateData::add_rate_data (std::map<std::set<std::pair<int, int> >, double>& rate_data,
									//
									const std::list<reac_t>& reactions,
									//
									const std::map<std::set<int>, std::map<int, double> >& bfdb) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::RateData::add_rate_data: ";

  using Model::operator<<;
  
  double dtemp;
  int    itemp;

  //IO::Marker funame_marker(funame);

  std::set<int> res;

  std::map<int, int> wg_map;
  
  bfdb_t::const_iterator bfi;
    
  for(int g = 0; g < well_partition.size(); ++g)
    //
    for(group_t::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w)
      //
      wg_map[*w] = g;

  for(std::list<reac_t>::const_iterator r = reactions.begin(); r != reactions.end(); ++r) {
    //
    if(rate_data.find(*r) != rate_data.end()) {
      //
      std::cerr << funame << "reaction already in the database\n";

      throw Error::Logic();
    }
    
    if(r->size() != 2) {
      //
      std::cerr << funame << "wrong number of reactants: " << r->size() << "\n";

      throw Error::Logic();
    }

    // bimolecular-to-bimolecular rate
    //
    if(r->begin()->first == Model::BIMOLECULAR && r->rbegin()->first == Model::BIMOLECULAR) {
      //
      rate_data[*r] = bb_rate(bim_to_index(r->begin()->second), bim_to_index(r->rbegin()->second));
      
    }
    // well-to-bimolecular
    //
    else if(r->begin()->first  == Model::BIMOLECULAR && r->rbegin()->first == Model::WELL ||
	    //
	    r->rbegin()->first == Model::BIMOLECULAR && r->begin()->first  == Model::WELL) {
      //
      int w, p;
      
      if(r->begin()->first  == Model::BIMOLECULAR) {
	//
	w = r->rbegin()->second;

	p = r->begin()->second;
      }
      else {
	//
	w = r->begin()->second;

	p = r->rbegin()->second;
      }

      if(bimolecular_group.find(w) != bimolecular_group.end()) {
	//
	res.insert(w);
	
	continue;
      }
      
      std::map<int, int>::const_iterator wgi = wg_map.find(w);

      if(wgi == wg_map.end()) {
	//
	std::cerr << funame << "well for " <<  Model::well(w).name() << "<--->" << Model::bimolecular(p).name()
	  //
		  << " reaction is not found in the "<< well_partition << " partition:\n"
	  //
		  << "bimolecular group: " << bimolecular_group << "\n";

	throw Error::Logic();
      }
      
      std::map<int, double>::const_iterator j;

      bfi = bf.find(well_partition[wgi->second]);

      if(bfi != bf.end()) {
	//
	j = bfi->second.find(w);

	if(j == bfi->second.end()) {
	  //
	  std::cerr << funame << Model::well(w).name() << " well is not in branching fraction database for the " << bfi->first << " group\n";
	  
	  throw Error::Logic();
	}
      }
      else {
	//
	bfi = bfdb.find(well_partition[wgi->second]);
	
	if(bfi == bfdb.end()) {
	  //
	  std::cerr << funame << well_partition[wgi->second] << " group is not in the branching fraction database\n";
	
	  throw Error::Logic();
	}

	j = bfi->second.find(w);

	if(j == bfi->second.end()) {
	  //
	  std::cerr << funame << Model::well(w).name() << " well is not in branching fraction database for the " << bfi->first << " group\n";
						  
	  throw Error::Logic();
	}
      }

      rate_data[*r] = wb_rate(wgi->second, bim_to_index(p)) * j->second;
    }
    // well-to-well rate
    //
    else if(r->begin()->first == Model::WELL && r->rbegin()->first == Model::WELL) {
      //
      const int& w1 = r->begin()->second;

      const int& w2 = r->rbegin()->second;

      if(bimolecular_group.find(w1) != bimolecular_group.end() || bimolecular_group.find(w2) != bimolecular_group.end()) {
	//
	if(bimolecular_group.find(w1) != bimolecular_group.end())
	  //
	  res.insert(w1);
	
	if(bimolecular_group.find(w2) != bimolecular_group.end())
	  //
	  res.insert(w2);
	
	continue;
      }
      
      std::map<int, int>::const_iterator wgi1 = wg_map.find(w1);

      std::map<int, int>::const_iterator wgi2 = wg_map.find(w2);

      if(wgi1 == wg_map.end() || wgi2 == wg_map.end()) {
	//
	std::cerr << funame << "one of the wells for " <<  Model::well(w1).name() << "<--->" << Model::well(w2).name()
	  //
		  << " reaction is not found in the " << well_partition << " partition:\n"
	  //
		  << funame << "bimolecular group: " << bimolecular_group << "\n";
	
	throw Error::Logic();
      }
      
      // wells belong to the same group
      //
      if(wgi1->second == wgi2->second) {
	//
	if(well_extension >= 0.) {
	  //
	  IO::log << IO::log_offset << funame << "WARNING: " << Model::well(w1).short_name() << " and " << Model::well(w2).short_name()
	    //
		  << " reactants belong to the same " << well_partition[wgi1->second]
	    //
		  << " group: check/increase the energy cutoff setting and/or well extension parameter\n";
	}

	res.insert(w1);

	res.insert(w2);
	
	continue;
      }
      
      std::map<int, double>::const_iterator j1, j2;
      
      bfi = bf.find(well_partition[wgi1->second]);

      if(bfi != bf.end()) {
	//
	j1 = bfi->second.find(w1);
      
	if(j1 == bfi->second.end()) {
	  //
	  std::cerr << funame << Model::well(w1).name() << " well is not in the branching fraction map for the " << bfi->first << " group\n";

	  throw Error::Logic();
	}
      }
      //
      // global branching fraction database
      //
      else {
	//
	bfi = bfdb.find(well_partition[wgi1->second]);
	
	if(bfi == bfdb.end()) {
	  //
	  std::cerr << funame << well_partition[wgi1->second] << " group is not in the branching fraction database\n";
	
	  throw Error::Logic();
	}
      
	j1 = bfi->second.find(w1);
      
	if(j1 == bfi->second.end()) {
	  //
	  std::cerr << funame << Model::well(w1).name() << " well is not in the branching fraction map for the " << bfi->first << " group\n";

	  throw Error::Logic();
	}
      }
      
      bfi = bf.find(well_partition[wgi2->second]);

      if(bfi != bf.end()) {
	//
	j2 = bfi->second.find(w2);
      
	if(j2 == bfi->second.end()) {
	  //
	  std::cerr << funame << Model::well(w2).name() << " well is not in the branching fraction map for the " << bfi->first << " group\n";

	  throw Error::Logic();
	}
      }
      //
      // global branching fraction database
      //
      else {
	//
	bfi = bfdb.find(well_partition[wgi2->second]);
	
	if(bfi == bfdb.end()) {
	  //
	  std::cerr << funame << well_partition[wgi2->second] << " group is not in the branching fraction databases\n";
	
	  throw Error::Logic();
	}
      
	j2 = bfi->second.find(w2);
      
	if(j2 == bfi->second.end()) {
	  //
	  std::cerr << funame << Model::well(w2).name() << " well is not in the branching fraction map for the " << bfi->first << " group\n";

	  throw Error::Logic();
	}
      }
      
      rate_data[*r] = ww_rate(wgi1->second, wgi2->second) * j1->second * j2->second;
    }
    else {
      //
      std::cerr << funame << "unknown reactant types: " << r->begin()->first << ", " << r->rbegin()->first << "\n";

      throw Error::Logic();
    }
  }

  return res;
}


// global lumping scheme
//
typedef MasterEquation::ReactiveComplex::lump_t lump_t;

lump_t MasterEquation::ReactiveComplex::lumping_scheme (double temperature, double pressure, int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::lumping_scheme: ";

  lump_t res;

  if(Model::bimolecular_size()) {
    //
    Model::bound_t bg = Model::bound_groups(temperature);
  
    for(Model::bound_t::const_iterator g = bg.begin(); g != bg.end(); ++g) {
      //
      lump_t ls = lumping_scheme(temperature, pressure, Model::ChemGraph(g->first), g->second, flags);

      for(std::list<std::set<int> >::const_iterator pi = ls.first.begin(); pi != ls.first.end(); ++pi)
	//
	res.first.push_back(*pi);

      for(std::set<int>::const_iterator gi = ls.second.begin(); gi != ls.second.end(); ++gi)
	//
	res.second.insert(*gi);
    }

    return res;
  }

  Model::ChemGraph full_graph;

  for(int w = 0; w < Model::well_size(); ++w)
    //
    full_graph.well_set.insert(w);

  for(int b = 0; b < Model::inner_barrier_size(); ++b)
    //
    full_graph.inner_set.insert(b);
  
  std::list<Model::ChemGraph> gl;

  int sb = full_graph.split(temperature, &gl);

  for(std::list<Model::ChemGraph>::const_iterator g = gl.begin(); g != gl.end(); ++g) {
      //
    lump_t ls = lumping_scheme(temperature, pressure, *g, std::make_pair((int)Model::INNER, sb), flags);

    for(std::list<std::set<int> >::const_iterator pi = ls.first.begin(); pi != ls.first.end(); ++pi)
      //
      res.first.push_back(*pi);

    for(std::set<int>::const_iterator gi = ls.second.begin(); gi != ls.second.end(); ++gi)
      //
      res.second.insert(*gi);
  }

  if(!res.first.size()) {
    //
    res.second.clear();

    res.first.push_back(full_graph.well_set);
  }
  
  return res;
}

lump_t MasterEquation::ReactiveComplex::lumping_scheme(double temperature, double pressure, const Model::ChemGraph& graph,
						       //
						       std::pair<int, int> diss, int flags)
{
  const char funame [] = "MasterEquation::ReaactiveComplex::lumping_scheme: ";

  using Model::operator<<;
  
#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  int    itemp;
  double dtemp;
  bool   btemp;

  lump_t res;

  if(!graph.is_connected()) {
    //
    std::cerr << funame << "graph is not connected: " << graph << "\n";

    throw Error::Logic();
  }
  
  const double be = Model::thermal_energy(diss, temperature);

  if(!(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "well_group: " << graph.well_set << "\n";

    IO::log << IO::log_offset << "dissociation barrier: " << diss << ": thermal energy = "
      //
	    << std::ceil(be / Phys_const::kcal * 10.) / 10. << " kcal/mol\n";
  }
  
  double er = be + excess_energy_over_temperature * temperature;
    
  double ec = be - energy_cutoff_over_temperature * temperature;

  ReactiveComplex rc(temperature, pressure, er, ec, graph, diss, flags | NO_WELL_EXTENSION | NOPRINT);

  if(!rc.there_are_bound_groups()) {
    //
    res.second = graph.well_set;

    return res;
  }

  if(graph.well_set.size() == 1) {
    //
    res.first.push_back(graph.well_set);

    return res;
  }
  
  std::list<Model::ChemGraph> gl;

  const int sb = graph.split(temperature, &gl);

  std::list<lump_t> ls;

  for(std::list<Model::ChemGraph>::const_iterator g = gl.begin(); g != gl.end(); ++g)
    //
    ls.push_back(lumping_scheme(temperature, pressure, *g, std::make_pair((int)Model::INNER, sb), flags));

  for(std::list<lump_t>::const_iterator l = ls.begin(); l != ls.end(); ++l) {
    //
    for(std::list<std::set<int> >::const_iterator g = l->first.begin(); g != l->first.end(); ++g)
      //
      res.first.push_back(*g);

    for(std::set<int>::const_iterator w = l->second.begin(); w != l->second.end(); ++w)
      //
      res.second.insert(*w);
  }

  if(!res.first.size()) {
    //
    res.second.clear();

    res.first.push_back(graph.well_set);
  }

  return res;
}
 
void MasterEquation::ReactiveComplex::branching_fraction(double  temperature, double pressure, const Model::ChemGraph& graph,
							 //
							 const std::set<int>& wg, bfdb_t& bfdb, int flags)
{
  const char funame [] = "MasterEquation::ReactiveComplex::branching_fraction: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  using Model::operator<<;
  
  if(bfdb.find(wg) != bfdb.end())
    //
    return;
  
  if(wg.size() == 1) {
    //
    std::map<int, double> res;

    res[*wg.begin()] = 1.;

    bfdb[wg] = res;

    return;
  }

  if(!graph.does_include(wg)) {
    //
    std::cerr << funame << wg << " partition group does not belong to the graph: " << graph << "\n";

    throw Error::Logic();
  }

#ifdef DEBUG
  //
  IO::Marker marker(funame);
  //
#else
  //
  IO::Marker marker;

  if(!(flags & NOPRINT)) {
    //
    marker.init(funame);
  }
  //
#endif
  
  Model::ChemGraph cg(wg);

  if(!cg.is_connected()) {

    cg = graph;

    btemp = true;

    // reduce the chemical graph to the smallest one that includes the considered partition group
    //
    while(btemp) {
      //
      std::list<Model::ChemGraph> gl;

      cg.split(temperature, &gl);

      btemp = false;

      // find if the splitted graphs include the considered partition group
      //
      for(std::list<Model::ChemGraph>::const_iterator g = gl.begin(); g != gl.end(); ++g)
	//
	if(g->does_include(wg)) {
	  //
	  btemp = true;

	  cg = *g;

	  break;
	}
    }

    if(!(flags & NOPRINT))
      //
      IO::log << IO::log_offset << "reactive complex:  " << graph << "\n"
	//
	      << IO::log_offset << "original group: " << wg << "\n"
	//
	      << IO::log_offset << "extended group: " << cg.well_set << "\n";
    
    cg.outer_set.clear();
  }
  
  // well-splitting barrier
  //
  const int sb = cg.split(temperature);
  
  // energy reference
  //
  const double er = std::max(Model::inner_barrier(sb).real_ground(), Model::inner_barrier(sb).thermal_energy(temperature))
    //
    + excess_energy_over_temperature * temperature;
 
  // energy cutoff
  //
  double ec = er;
  
  for(std::set<int>::const_iterator b = cg.inner_set.begin(); b != cg.inner_set.end(); ++b) {
    //
    dtemp = Model::inner_barrier(*b).thermal_energy(temperature);
    
    if(dtemp < ec)
      //
      ec = dtemp;
  }
  
  ec -= energy_cutoff_over_temperature * temperature;

  flags |= NO_WELL_EXTENSION;
  
  ReactiveComplex rc(temperature, pressure, er, ec, cg, flags);

  std::map<int, double> bf = rc.branching_fraction(flags);

  if(bf.size() == wg.size()) {
    //
    bfdb[wg] = bf;

    if(!(flags & NOPRINT)) {
      //
      IO::log << IO::log_offset << "branching fractions: ";
    
      for(std::map<int, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
	//
	IO::log << "  " << i->second << "/" << Model::well(i->first).short_name();
      
      IO::log << "\n";
    }
    
    return;
  }
  
  std::map<int, double> res;

  dtemp = 0.;
	
  for(std::map<int, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
    //
    if(wg.find(i->first) != wg.end())
      //
      dtemp += i->second;

  if(!(flags & NOPRINT))
    //
    IO::log << IO::log_offset << wg << " group total branching fraction = " << dtemp << "\n";

  if(dtemp < 1.e-12) {
    //
    if(!(flags & NOPRINT))
      //
      IO::log << IO::log_offset << "WARNING: total branching contribution of the  " << wg
	//
	      << " group is close to zero: assuming statistical branching\n";

    dtemp = Model::weight(wg, temperature);
      
    for(std::set<int>::const_iterator w = wg.begin(); w != wg.end(); ++w)
      //
      res[*w] = Model::well(*w).weight(temperature) / std::exp(Model::well(*w).ground() / temperature) / dtemp;

    bfdb[wg] = res;

    if(!(flags & NOPRINT)) {
      //
      IO::log << IO::log_offset << "branching fractions: ";

      for(std::map<int, double>::const_iterator i = res.begin(); i != res.end(); ++i)
	//
	IO::log << "  " << i->second << "/" << Model::well(i->first).short_name();
      
      IO::log << "\n";
    }
    
    return;
  }
	
  for(std::map<int, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
    //
    if(wg.find(i->first) != wg.end())
      //
      res[i->first] = i->second / dtemp;

  if(!(flags & NOPRINT)) {
    //
    IO::log << IO::log_offset << "branching fractions: ";

    for(std::map<int, double>::const_iterator i = res.begin(); i != res.end(); ++i)
      //
      IO::log << "  " << i->second << "/" << Model::well(i->first).short_name();

    IO::log << "\n";
  }
  
  bfdb[wg] = res;
  //
}// branching fraction

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

  set_energy_step();
  
  // set default energy reference
  //
  if(flags & DEFAULT_EREF) {

    double e, s, c;
    
    btemp = true;

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      Model::outer_barrier(b).esc_parameters(temperature(), e, s, c);

      IO::log << IO::log_offset << Model::outer_barrier(b).short_name() << ": e = " << e / Phys_const::kcal << " kcal/mol, c = " << c << std::endl;
      
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
      
      IO::log << IO::log_offset << Model::inner_barrier(b).short_name()
	//
	      << ": e = " << e / Phys_const::kcal << " kcal/mol, c = " << c << std::endl;
      
      e += Model::inner_barrier(b).ground() + 5. * temperature() * std::sqrt(c);
      
      if(btemp || _energy_reference < e) {
	//
	btemp = false;
	
	_energy_reference = e;
      }
    }
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

  IO::log << IO::log_offset << "Temperature      = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";
  
  IO::log << IO::log_offset << "Energy reference = " << std::ceil(energy_reference() / Phys_const::incm) << " 1/cm\n";
  
  IO::log << IO::log_offset << "Energy step      = " << energy_step()      / Phys_const::incm << " 1/cm\n";

  {
    IO::Marker set_marker("setting wells, barriers, and bimolecular");

    _well.resize(Model::well_size());
    
    // wells
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      _well[w] = SharedPointer<Well>(new Well(Model::well(w)));

    // well-to-well barriers
    //
    _inner_barrier.resize(Model::inner_barrier_size());
    
    for(int b = 0; b < Model::inner_barrier_size(); ++b)
      //
      _inner_barrier[b] = SharedPointer<Barrier>(new Barrier(Model::inner_barrier(b)));
    
    // well-to-bimolecular barriers
    //
    _outer_barrier.resize(Model::outer_barrier_size());

    for(int b = 0; b < Model::outer_barrier_size(); ++b)
      //
      _outer_barrier[b] = SharedPointer<Barrier>(new Barrier(Model::outer_barrier(b)));
  }
  // checking if wells are deeper than the barriers between them
  
  // inner barriers
  //
  for(int b = 0; b < Model::inner_barrier_size(); ++b) {
    //
    const int& w1 = Model::inner_connect(b).first;
    
    const int& w2 = Model::inner_connect(b).second;
    
    itemp = well(w1).size() <  well(w2).size() ? well(w1).size() : well(w2).size();
    
    if(itemp  < inner_barrier(b).size()) {
      //
      IO::log << IO::log_offset << Model::inner_barrier(b).short_name()
	//
	      << " barrier top is lower than the bottom of one of the wells it connects => truncating\n";
      
      _inner_barrier[b]->truncate(itemp);
    }
  }
  
  
  // outer barriers
  //
  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
    //
    itemp = well(Model::outer_connect(b).first).size();
    
    if(itemp  < outer_barrier(b).size()) {
      //
      IO::log << IO::log_offset << Model::outer_barrier(b).short_name()
	//
	      << " barrier top is lower than the bottom of the well it connects => truncating\n";
      
      _outer_barrier[b]->truncate(itemp);
    }
  }
  
  // clipping the number of states with the maximum rate constant
  
  // inner barriers
  //
  if(rate_max > 0.) {
    //
    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      const int& w1 = Model::inner_connect(b).first;
      
      const int& w2 = Model::inner_connect(b).second;
      
      for(int i = 0; i < _inner_barrier[b]->size(); ++i) {
	//
	dtemp = well(w1).state_density(i);

	dtemp = dtemp < well(w2).state_density(i) ? dtemp : well(w2).state_density(i);
	
	dtemp *= rate_max * 2. * M_PI;
	
	if(_inner_barrier[b]->state_number(i) > dtemp)
	  //
	  _inner_barrier[b]->state_number(i) = dtemp;
      }
    }
    
    // outer barriers
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int& w = Model::outer_connect(b).first;
      
      for(int i = 0; i < _outer_barrier[b]->size(); ++i) {
	//
	dtemp = rate_max * 2. * M_PI * well(w).state_density(i);

	if(_outer_barrier[b]->state_number(i) > dtemp)
	  //
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
    
    evec_out << "Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5)<< "K\n"
      //
	     <<"THERMAL DISTRIBUTIONS:\n"
      //
	     << std::setw(13) << "E, kcal/mol";
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      evec_out << std::setw(13) << Model::well(w).short_name();
    
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
  
  //IO::log << std::setprecision(3);

  if(Model::well_size()) {
    //
    IO::log << IO::log_offset << "Bound Species (D1/D2 - density of states at dissociation energy (DE)/reference energy):\n"
      //
	    << IO::log_offset << std::setw(5)  << "Name" << std::setw(7)  << "Depth" << std::setw(7)  << "DE"
      //
	    << std::setw(Model::log_precision + 7) << "*D1" << std::setw(Model::log_precision + 7) << "*D2"  << "\n"
      //
	    << IO::log_offset  << std::setw(5) << "" << std::setw(7)  << "1/cm" << std::setw(7)  << "1/cm"
      //
	    << std::setw(Model::log_precision + 7) << "cm" << std::setw(Model::log_precision + 7) << "cm"  << "\n";

    for(int w = 0; w < Model::well_size(); ++w) {
      //
      IO::log << IO::log_offset
	//
	      << std::setw(5) << Model::well(w).short_name();
      
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

      IO::log << std::setw(Model::log_precision + 7) << well(w).state_density(itemp) * Phys_const::incm
	//
	      << std::setw(Model::log_precision + 7) << well(w).state_density(0) * Phys_const::incm
	      << "\n";
    }
  }

  if(Model::inner_barrier_size()) {
    //
    IO::log << IO::log_offset << "Well-to-Well Barriers (N - number of states at the reference energy):\n"
      //
	    << IO::log_offset << std::setw(5)  << "Name" << std::setw(7)  << "Height" << std::setw(Model::log_precision + 7) << "*N" << "\n"
      //
	    << IO::log_offset << std::setw(5)  << ""     << std::setw(7)  << "1/cm"    << std::setw(Model::log_precision + 7) << ""  << "\n";

    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::inner_barrier(b).short_name();
      
      dtemp = energy_bin(inner_barrier(b).size());
      
      IO::log << std::setw(7) << (int)std::ceil(dtemp / Phys_const::incm)
	//
	      << std::setw(Model::log_precision + 7) << inner_barrier(b).state_number(0)
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
	    << std::setw(Model::log_precision + 7) << "*N"
      //
	    << "\n"
      //
	    << IO::log_offset
      //
	    << std::setw(5)  << ""
      //
	    << std::setw(7)  << "1/cm"
      //
	    << std::setw(Model::log_precision + 7) << ""
      //
	    << "\n";

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::outer_barrier(b).short_name();
      //
      dtemp = energy_bin(outer_barrier(b).size());
      
      IO::log << std::setw(7) <<  (int)std::ceil(dtemp / Phys_const::incm)
	//
	      << std::setw(Model::log_precision + 7) << outer_barrier(b).state_number(0)
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
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    if(!Model::bimolecular(p).dummy())
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(p).short_name();
  
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(5) << Model::well(w).short_name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1)
      //
      IO::log << std::setw(Model::log_precision + 7) << well(w).weight() / well(w1).weight();
    
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy()) {
	//
	dtemp = (Model::bimolecular(p1).ground() - energy_reference()) / temperature();

	IO::log << std::setw(Model::log_precision + 7);

	if(dtemp > -Limits::exp_pow_max()) {
	  //
	  IO::log << well(w).weight() / Model::bimolecular(p1).weight(temperature()) * std::exp(dtemp) / bpu;
	}
	else
	  //
	  IO::log << "***";
      }
    
    IO::log << "\n";
  }

  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    if(Model::bimolecular(p).dummy())
      //
      continue;
    
    IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).short_name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1) {
      //
      dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();
	
      IO::log << std::setw(Model::log_precision + 7);

      if(dtemp > -Limits::exp_pow_max()) {
	  //
	IO::log << Model::bimolecular(p).weight(temperature()) / std::exp(dtemp) * bpu / well(w1).weight();
      }
      else
	//
	IO::log << "***";
    }
  
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy()) {
	//
	dtemp = (Model::bimolecular(p).ground() - Model::bimolecular(p1).ground()) / temperature();

	IO::log << std::setw(Model::log_precision + 7);

	if(-Limits::exp_pow_max() < dtemp && dtemp < Limits::exp_pow_max()) {
	  //
	  IO::log << Model::bimolecular(p).weight(temperature()) / Model::bimolecular(p1).weight(temperature()) / std::exp(dtemp);
	}
	else
	  //
	  IO::log << "***";
      }
    
    IO::log << "\n";
  }

  IO::log << IO::log_offset << "Real equilibrium constants(bimolecular units - cm^3):\n"
    //
	  << IO::log_offset << std::setw(5) << "Q/Q";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    if(!Model::bimolecular(p).dummy())
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(p).short_name();
  
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(5) << Model::well(w).short_name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1) {
      //
      dtemp = (Model::well(w).ground() - Model::well(w1).ground()) / temperature();

      IO::log << std::setw(Model::log_precision + 7);
      
      if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	//
	IO::log << Model::well(w).weight(temperature()) / Model::well(w1).weight(temperature()) / std::exp(dtemp);
      }
      else
	//
	IO::log << "***";
    }
  
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy()) {
	//
	dtemp = (Model::well(w).ground() - Model::bimolecular(p1).ground()) / temperature();

	IO::log << std::setw(Model::log_precision + 7);
      
	if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	  //
	  IO::log << Model::well(w).weight(temperature()) / Model::bimolecular(p1).weight(temperature()) / std::exp(dtemp) / bpu;
	}
	else
	  //
	  IO::log << "***";
      }
    
    IO::log << "\n";
  }

  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    if(Model::bimolecular(p).dummy())
      //
      continue;
    
    IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).short_name();
    
    for(int w1 = 0; w1 < Model::well_size(); ++w1) {
      //
      dtemp = (Model::bimolecular(p).ground() - Model::well(w1).ground()) / temperature();

      IO::log << std::setw(Model::log_precision + 7);
      
      if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	//
	IO::log << Model::bimolecular(p).weight(temperature()) / Model::well(w1).weight(temperature()) / std::exp(dtemp) * bpu;
      }
      else
	//
	IO::log << "***";
    }
    
    for(int p1 = 0; p1 < Model::bimolecular_size(); ++p1)
      //
      if(!Model::bimolecular(p1).dummy()) {
	//
	dtemp = (Model::bimolecular(p).ground() - Model::bimolecular(p1).ground()) / temperature();

	IO::log << std::setw(Model::log_precision + 7);
      
	if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	  //
	  IO::log << Model::bimolecular(p).weight(temperature()) / Model::bimolecular(p1).weight(temperature()) / std::exp(dtemp);
	}
	else
	  //
	  IO::log << "***";
      }
    
    IO::log << "\n";
  }

  //IO::log << std::setprecision(6);

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
      
      dtemp = (Model::inner_barrier(b).ground() - Model::well(w1).ground()) / temperature();

      if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	//
	rate_data[std::make_pair(w1, w2)] = temperature() / 2. / M_PI / Phys_const::herz
	  //
	  * Model::inner_barrier(b).weight(temperature()) / Model::well(w1).weight(temperature()) / std::exp(dtemp);
      }
      
      dtemp = (Model::inner_barrier(b).ground() - Model::well(w2).ground()) / temperature();

      if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	//
	rate_data[std::make_pair(w2, w1)] = temperature() / 2. / M_PI / Phys_const::herz
	  //
	  * Model::inner_barrier(b).weight(temperature()) / Model::well(w2).weight(temperature()) / std::exp(dtemp);
      }
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
      
      dtemp = (Model::outer_barrier(b).ground() - Model::well(w).ground()) / temperature();

      if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	//
	rate_data[std::make_pair(w, Model::outer_connect(b).second + Model::well_size())] = temperature() / 2. / M_PI / Phys_const::herz
	  //
	  * Model::outer_barrier(b).weight(temperature()) / Model::well(w).weight(temperature()) / std::exp(dtemp);
      }
    }

    // bimolecular-to-well
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int& p = Model::outer_connect(b).second;
      
      if(!Model::bimolecular(p).dummy()) {
	//
	dtemp = (Model::outer_barrier(b).ground() - Model::bimolecular(p).ground()) / temperature();

	if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	  //
	  rate_data[std::make_pair(p + Model::well_size(), Model::outer_connect(b).first)] = temperature() / 2. / M_PI / bru
	    //
	    * Model::outer_barrier(b).weight(temperature()) / Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
	}
      }
    }
  }

  // capture
  //
  capture.clear();

  std::set<int> remove_capture;
  
  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
    //
    const int& w = Model::outer_connect(b).first;
    
    const int& p = Model::outer_connect(b).second;
    
    if(!Model::bimolecular(p).dummy()) {
      //
      dtemp = (Model::outer_barrier(b).ground() - Model::bimolecular(p).ground()) / temperature();

      if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
	  //
	capture[p + Model::well_size()] += temperature() / 2. / M_PI / bru
	  //
	  * Model::outer_barrier(b).weight(temperature()) / Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
      }
      else
	//
	remove_capture.insert(p + Model::well_size());
    }
    
    dtemp = (Model::outer_barrier(b).ground() - Model::well(w).ground()) / temperature();
    
    if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
      //
      capture[w] += temperature() / 2. / M_PI / Phys_const::herz
	//
	* Model::outer_barrier(b).weight(temperature()) / Model::well(w).weight(temperature()) / std::exp(dtemp);
    }
    else
      //
      remove_capture.insert(w);
  }

  for(int b = 0; b < Model::inner_barrier_size(); ++b) {
    //
    const int& w1 = Model::inner_connect(b).first;
    
    const int& w2 = Model::inner_connect(b).second;
    
    dtemp = (Model::inner_barrier(b).ground() - Model::well(w1).ground()) / temperature();
    
    if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
      //
      capture[w1] += temperature() / 2. / M_PI / Phys_const::herz
	//
	* Model::inner_barrier(b).weight(temperature()) / Model::well(w1).weight(temperature()) / std::exp(dtemp);
    }
    else
      //
      remove_capture.insert(w1);
    
    dtemp = (Model::inner_barrier(b).ground() - Model::well(w2).ground()) / temperature();
    
    if(dtemp < Limits::exp_pow_max() && dtemp > -Limits::exp_pow_max()) {
      //
      capture[w2] += temperature() / 2. / M_PI / Phys_const::herz
	//
	* Model::inner_barrier(b).weight(temperature()) / Model::well(w2).weight(temperature()) / std::exp(dtemp);
    }
    else
      //
      remove_capture.insert(w2);
  }

  for(std::set<int>::const_iterator cit = remove_capture.begin(); cit != remove_capture.end(); ++cit)
    //
    capture.erase(*cit);
  
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
}// set function

/********************************************************************************************
 ************************************** SETTING WELL ****************************************
 ********************************************************************************************/

// set state density, relaxation mode basis, collisional energy transfer kernel
//
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

  IO::log << IO::log_offset << model.short_name() << " Well: thermal energy = "
    //
	  << std::ceil(hval / Phys_const::kcal * 10.) / 10.
    //
	  << " kcal/mol\n";
  
  /***************************** SETTING DENSITY OF STATES ************************************/

  double ener, ext_ener = -1.;

  double base_ener = model.ground();

  // well extension
  //
  if(model.well_extension() && well_extension >= 0. && well_cutoff > 0.) {
    //
    base_ener += hval * (1. - well_extension);

    // well extension energy correction
    //
    dtemp = model.well_ext_cap + model.ground();
    
    if(model.well_ext_cap > 0. && dtemp < base_ener)
      //
      base_ener = dtemp + well_ext_corr * (base_ener - dtemp);

    ext_ener = well_cutoff * temperature();

    dtemp = base_ener - model.dissociation_limit;
    
    if(dtemp < 0.)
      //
      ext_ener += dtemp;
  }
  else if(well_cutoff > 0.) {
    //
    dtemp = model.dissociation_limit - well_cutoff * temperature();

    if(dtemp > base_ener)
      //
      base_ener = dtemp;
  }

  if(is_global_cutoff) {
    //
    dtemp = base_ener - lower_global_cutoff;

    if(dtemp < 0.) {
      //
      base_ener = lower_global_cutoff;

      ext_ener = -1.;
    }    
    else if(ext_ener > dtemp)
      //
      ext_ener = dtemp;          
  }
  
  int base_size = std::ceil((energy_reference() - base_ener) / energy_step());

  int ext_size = ext_ener < 0. ? 0 : std::ceil(ext_ener / energy_step());
  
  if(ext_size)
    //
    IO::log << IO::log_offset << model.short_name() << " Well: extension size = " << ext_size << "\n";
											   
  itemp = base_size + ext_size;
  
  _state_density.resize(itemp);
  
  resize_thermal_factor(itemp);

  ener = energy_reference();
  
  for(int e = 0; e < base_size; ++e, ener -= energy_step()) {
    //
    dtemp = model.states(ener);
    
    if(dtemp <= 0.) {
      //
      IO::log << IO::log_offset << model.short_name()  << " Well: WARNING: nonpositive density at " 
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
      itemp = (int)std::ceil(model.kernel(b).cutoff_energy(temperature()) / energy_step());

      IO::log << IO::log_offset << model.short_name() << " Well: collisional kernel bandwidth = " << itemp << "\n";

      std::vector<double> energy_transfer(itemp);
	
      for(int i = 0; i < energy_transfer.size(); ++i)
	//
	energy_transfer[i] = model.kernel(b)((double)i * energy_step(), temperature());
      
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
	    std::cerr << model.short_name() << " Well: cannot satisfy the constant collision rate at energy = "
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
	    IO::log << IO::log_offset << model.short_name()
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
  
  IO::log << IO::log_offset << model.short_name() << " well: kernel diagonal elements:\n";
  
  for(int i = 0; i < size(); ++i)
    //
    IO::log << IO::log_offset << std::setw(5) << i  << std::setw(Model::log_precision + 7) << _kernel(i, i)  << "\n";
  
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

  std::vector<double> nfac;
  
  // Always initialize CRM basis when this function is called
  _crm_basis.resize(size(), crm_size());
  
  _crm_basis = 0.;

  nfac.resize(crm_size());
  
  _weight = 0.;
  
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
  
  // Set the CRM basis flag to indicate initialization is complete
  with_crm_basis = 1;
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
    dtemp = model.collision(b)(temperature()) * Model::buffer_fraction(b);
    
    _kernel_fraction[b] = dtemp;
    
    _collision_factor += dtemp;
  }

  for(int b = 0; b < Model::buffer_size(); ++b)
    //
    _kernel_fraction[b] /= _collision_factor;

  // state density
  //
  start_time = std::clock();

  _set_state_density(model);

  IO::log << IO::log_offset << model.short_name() << " Well: density of states done, elapsed time[sec] = "
    //
	  << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;

  // collisiona relaxation kernel
  //
  start_time = std::clock();

  _set_kernel(model);

  IO::log << IO::log_offset << model.short_name() << " Well: collisional energy transfer kernel done, elapsed time[sec] = "
    //
	  << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;

  // CRM basis
  //
  _set_crm_basis();

  if(with_crm_basis) {
    //
    start_time = std::clock();

    IO::log << IO::log_offset << model.short_name() << " Well: relaxation modes basis done, elapsed time[sec] = "
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

    IO::log << IO::log_offset << model.short_name()
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

      IO::log << IO::log_offset << model.short_name() << " Well: energy transfer kernel lowest eigenvalues:\n";

      itemp = size() < 10 ? size() : 10;
  
      for(int i = 0; i < itemp; ++i)
      //
      IO::log << IO::log_offset << std::setw(Model::log_precision + 7) << kval[i] << "\n";

      stemp = "kvec_" + model.short_name() + "_" + IO::String(int(temperature() / Phys_const::kelv)) + ".dat";

      std::ofstream kvec_out(stemp.c_str());
  
      for(int e = 0; e < size(); ++e) {
      //
      for(int i = 0; i < itemp; ++i)
      //
      kvec_out << std::setw(Model::log_precision + 7) << kvec(e, i) / boltzman_sqrt(e) * boltzman_sqrt(size() - 1);

      kvec_out << "\n";
      }
    */
  
    // minimal collisional relaxation eigenvalue
    //
    start_time = std::clock();

    vtemp = _crm_kernel.eigenvalues();
  
    _min_relax_eval = vtemp.front();
  
    _max_relax_eval = vtemp.back();

    IO::log << IO::log_offset << model.short_name() << " Well: relaxation eigenvalues done, elapsed time[sec] = "
      //
	    << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;  

    IO::log << IO::log_offset << model.short_name() << " Well: minimal relaxation eigenvalue = "
      //
	    << _min_relax_eval << "\n";
  
    IO::log << IO::log_offset << model.short_name() << " Well: maximal relaxation eigenvalue = "
      //
	    << _max_relax_eval << "\n";
  }
  
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
    
    IO::log << IO::log_offset << model.short_name() << " Well: symmetrized kernel eigenvalues done, elapsed time[sec] = "
    << double(std::clock() - start_time) / CLOCKS_PER_SEC <<  std::endl;  
    IO::log << IO::log_offset << model.short_name() << " Well: symmetrized kernel eigenvalues: "
    << std::setw(Model::log_precision + 7) << vtemp[0] << std::setw(Model::log_precision + 7) << vtemp[1] << std::setw(Model::log_precision + 7) << vtemp.back() << "\n";
    
    #endif
  

    // escape rate
    //
    if(model.escape_size()) {
    //
    IO::log << IO::log_offset << model.short_name() << " Well: Escape rate:\n";
    
    _escape_rate.resize(size());
    
    double ener = energy_reference();

    for(int i = 0; i < size(); ++i, ener -= energy_step()) {
    //
    dtemp = model.escape_rate(ener);

    IO::log << IO::log_offset 
    << "    E[kcal/mol] = " << std::setw(Model::log_precision + 7) << ener / Phys_const::kcal
    << "    rate[1/sec] = " << std::setw(Model::log_precision + 7) << dtemp / Phys_const::herz << "\n";
      
    _escape_rate[i] = dtemp;
    }
    }
  */
    
  // radiational transitions
  //
  if(model.oscillator_size()) {
    //
    if(!with_crm_basis) {
      //
      std::cerr << funame << "CRM basis not initialized\n";

      throw Error::Init();
    }
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
  
  IO::log << IO::log_offset << model.short_name()
    //
	  << " Well:       grid size = " << size() << "\n"
    //
	  << IO::log_offset << model.short_name()
    //
	  << " Well:      real depth = " << int(model.ground() / Phys_const::incm) << " 1/cm\n"
    //
	  << IO::log_offset << model.short_name()
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

  dtemp = model.ground();

  if(is_global_cutoff && dtemp < lower_global_cutoff)
    //
    dtemp = lower_global_cutoff;
  
  int new_size = (int)std::ceil((energy_reference() - dtemp) / energy_step());

  if(new_size <= 0) {
    //
    std::cerr << funame << "non-positive size: " << new_size << ": check reference energy\n";

    throw Error::Range();
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
      IO::log << IO::log_offset  << model.short_name() << " Barrier: nonpositive number of states at "
	//
	      << ener / Phys_const::incm  << " 1/cm => truncating\n";
      
      _state_number.resize(e);
      
      break;
    }
    
    _state_number[e] = dtemp;
    
    _weight += dtemp * thermal_factor(e);
  }
  
  _weight *= energy_step() / temperature();

  IO::log << IO::log_offset << model.short_name()
    //
	  << " Barrier:        grid size = " << size() << "\n"
    //
	  << IO::log_offset << model.short_name()
    //
	  << " Barrier:      real height = " << int(model.ground() / Phys_const::incm) << " 1/cm\n"
    //
	  << IO::log_offset << model.short_name()
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

  _weight = 0.;
  
  for(int i = 0; i < new_size; ++i)
    //
    _weight += _state_number[i] * thermal_factor(i);
  
  _weight *= energy_step() / temperature();

  _state_number.resize(new_size);
  
  IO::log << IO::log_offset << funame << "grid size = " << size() << "\n"
    //
	  << IO::log_offset << funame << "effective weight = " << _weight << "\n";
}

/********************************************************************************************
 ******************************* LOW CHEMICAL EIGENVALUE METHOD *****************************
 ********************************************************************************************/

// the method based on the assumption that the chemical eigenvalues are small in comparison
//
// to the relaxation eigenvalues

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

  if(!with_crm_basis) {
    //
    std::cerr << funame << "crm basis is not initialized\n";

    throw Error::Init();
  }
  
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
	  << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";

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

  Lapack::Vector chem_eval;

  if(Mpack::mp_type == Mpack::DOUBLE) {
    //
    chem_eval = Offload::eigenvalues(k_11, &chem_evec);
  }
  else
    //
    chem_eval = Mpack::eigenvalues(k_11, &chem_evec);

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

  IO::log << IO::log_offset << "eigenvector populations normalized:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(Model::log_precision + 7) << "*R"
    //
	  << std::setw(Model::log_precision + 7) << "*P";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
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
	    << std::setw(Model::log_precision + 7) << chem_eval[l] / min_relax_eval
      //
	    << std::setw(Model::log_precision + 7) << rel_proj[l];
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = chem_evec(w, l) * well(w).weight_sqrt() / max_pop;
      //
      IO::log << std::setw(Model::log_precision + 7);
      
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
	  << std::setw(Model::log_precision + 7) << "*Q"
    //
	  << std::setw(Model::log_precision + 7) << "*P";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
  IO::log << "\n";

  for(int l = 0; l < Model::well_size(); ++l) {
    //
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(Model::log_precision + 7) << chem_eval[l] / well(wmin).collision_frequency()
      //
	    << std::setw(Model::log_precision + 7) << rel_proj[l];
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << chem_evec(w, l);
    
    IO::log << "\n";
  }
  
  IO::log << IO::log_offset
    //
	  << std::setw(5) << "*Z"
    //
	  << std::setw(Model::log_precision + 7) << "---"
    //
	  << std::setw(Model::log_precision + 7) << "---";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << well(w).weight_sqrt();
  
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
	  << " - well partition function square root\n";
  
  // rate coefficients connecting eigenstates to bimolecular products
  //
  Lapack::Matrix chem_bim = chem_evec.transpose() * k_13;

#ifdef DEBUG

  IO::log << IO::log_offset << "eigenvector-to-bimolecular projections:\n";

  IO::log << IO::log_offset << std::setw(5) << "L\\P";
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(p).short_name();
  
  IO::log << "\n";
  
  for(int l = 0; l < Model::well_size(); ++l) {
    //
    IO::log << IO::log_offset << std::setw(5) << l;
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(Model::log_precision + 7) << chem_bim(l, p);
    
    IO::log << "\n";
  }

#endif

  // bimolecular rate units
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  // bimolecular-to-bimolecular rate coefficients
  //
  for(int i = 0; i < Model::bimolecular_size(); ++i) {
    //
    if(Model::bimolecular(i).dummy())
      //
      continue;

    dtemp = (Model::bimolecular(i).ground() - energy_reference()) / temperature();

    if(dtemp < -Limits::exp_pow_max())
      //
      continue;

    const double fac = Model::bimolecular(i).weight(temperature()) / std::exp(dtemp);
    
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
	
      dtemp /= fac * bru;
	
      rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 	
    }
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
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group);
    }

    // convert chemical eigenvectors in the new basis
    //
    pop_chem = well_partition.basis().transpose() * pop_chem;

    std::vector<int>    group_index = well_partition.group_index();
    
    std::vector<double>      weight = well_partition.weight();
    
    Lapack::Matrix            basis = well_partition.basis();

    std::vector<double> real_weight = well_partition.real_weight(temperature());
    
    std::vector<double> real_ground = well_partition.real_ground();
    
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
	    << IO::log_offset << std::setw(2) << "#"  << std::setw(Model::log_precision + 7)
      //
	    << "new name" << IO::first_offset
      //
	    << "group\n";
    
    for(int g = 0; g < well_partition.size(); ++g) {
      //
      IO::log << IO::log_offset << std::setw(2) << g << std::setw(Model::log_precision + 7)
	//
	      << Model::well(group_index[g]).short_name() << IO::first_offset;
      
      for(Group::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w) {
	//
	if(w != well_partition[g].begin())
	  //
	  IO::log << "+";
	
	IO::log << Model::well(*w).short_name();
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
	IO::log << std::setw(Model::log_precision + 7) << vtemp[l] / well(wmin).collision_frequency();
      
      IO::log << "\n";
      
      // well-to-well rate coefficients
      //
      for(int i = 0; i < chem_size; ++i) {
	//
	dtemp = (real_ground[i] - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = real_weight[i] / std::exp(dtemp);
	
	for(int j = 0; j < chem_size; ++j) {
	  //
	  dtemp = k_11(i, j) * std::sqrt(weight[i] * weight[j]) * energy_step() / fac / Phys_const::herz;
	  
	  if(i != j)
	    //
	    dtemp = -dtemp;
	  
	  rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp;
	}
      }
      
      // well-to-bimolecular rate coefficients
      //
      for(int w = 0; w < chem_size; ++w) {
	//
	dtemp = (real_ground[w] - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = real_weight[w] / std::exp(dtemp);
	
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = k_13(w, p) * std::sqrt(weight[w]) * energy_step() / fac / Phys_const::herz;
      }
      
      // bimolecular-to-well rate coefficients
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(!Model::bimolecular(p).dummy()) {
	  //
	  dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	  if(dtemp < -Limits::exp_pow_max())
	    //
	    continue;
	  
	  const double fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
	  
	  for(int w = 0; w < chem_size; ++w)
	    //
	    rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = k_13(w, p) * std::sqrt(weight[w]) * energy_step() / fac / bru;
	}
      
      break;
      //
    case DIAGONALIZATION:
      //
      // well-to-well rate coefficients
      //
      for(int i = 0; i < chem_size; ++i) {
	//
	dtemp = (real_ground[i] - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = real_weight[i] / std::exp(dtemp);
	
	for(int j = 0; j < chem_size; ++j) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_direct(j, l) * m_inverse(l, i) * chem_eval[l];
	  
	  dtemp *= std::sqrt(weight[i] * weight[j]) * energy_step() / fac / Phys_const::herz;
	  
	  if(i != j)
	    //
	    dtemp = -dtemp;
	  
	  rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp;
	}
      }
      
      // well-to-bimolecular rate coefficients
      //
      for(int w = 0; w < chem_size; ++w) {
	//
	dtemp = (real_ground[w] - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = real_weight[w] / std::exp(dtemp);
	
	for(int p = 0; p < Model::bimolecular_size(); ++p) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_inverse(l, w) * chem_bim(l, p);
	  
	  dtemp *= std::sqrt(weight[w]) * energy_step() / fac / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
	}
      }
      
      // bimolecular-to-well rate coefficients
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(!Model::bimolecular(p).dummy()) {
	  //
	  dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	  if(dtemp < -Limits::exp_pow_max())
	    //
	    continue;
	  
	  const double fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
	  
	  for(int w = 0; w < chem_size; ++w) {
	    //
	    dtemp = 0.;
	    
	    for(int l = 0; l < chem_size; ++l)
	      //
	      dtemp += m_direct(w, l) * chem_bim(l, p);
	    
	    dtemp *= std::sqrt(weight[w]) * energy_step() / fac / bru;
	    
	    rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	  }
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
    for(int i = 0; i < Model::well_size(); ++i) {
      //
      dtemp = (Model::well(i).ground() - energy_reference()) / temperature();

      if(dtemp < -Limits::exp_pow_max())
	//
	continue;
      
      const double fac = Model::well(i).weight(temperature()) / std::exp(dtemp);
      
      for(int j = 0; j < Model::well_size(); ++j) {
	//
	dtemp = k_11(i, j) * well(i).weight_sqrt() * well(j).weight_sqrt() * energy_step() / fac / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(i, j)] = dtemp;
      }
    }
    
    // well-to-bimolecular rate coefficients
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = (Model::well(w).ground() - energy_reference()) / temperature();

      if(dtemp < -Limits::exp_pow_max())
	//
	continue;
      
      const double fac = Model::well(w).weight(temperature()) / std::exp(dtemp);
      
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	rate_data[std::make_pair(w, Model::well_size() + p)] = k_13(w, p) * well(w).weight_sqrt() * energy_step() / fac / Phys_const::herz;
    }
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(!Model::bimolecular(p).dummy()) {
	//
	dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
      
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  rate_data[std::make_pair(Model::well_size() + p, w)] = k_13(w, p) * well(w).weight_sqrt() * energy_step() / fac / bru;
      }
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
  
  Lapack::Vector eigenvalue;
  
  Lapack::Matrix eigenvector;
  
  std::map<int, int> well_index_map;
  
  std::vector<int>   index_well_map;

  int size () const { return eigenvalue.size(); }
};

void MasterEquation::well_reduction_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags)
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
  
  IO::log << "\t  Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n"
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

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    if(!w || well(w).size() > itemp)
      //
      itemp = well(w).size();
  }
  
  const int64_t ener_index_max = itemp;

  /*******************************************************************
   ********************* PARTITIONING THE WELLS **********************
   *******************************************************************/

  std::vector<KineticBasis> kinetic_basis(ener_index_max);

  {// kinetically active basis

    IO::Marker well_part_marker("kinetically active basis");

#pragma omp parallel for default(shared) private(itemp,dtemp) schedule(static)
  
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
      Lapack::SymmetricMatrix km(well_array.size());
      
      km = 0.;

      // nondiagonal isomerization contribution
      //
      for(int b = 0; b < Model::inner_barrier_size(); ++b) {
	//
	if(e < inner_barrier(b).size()) {
	  //
	  int w1 = Model::inner_connect(b).first;
	  
	  int w2 = Model::inner_connect(b).second;
      
	  std::map<int, int>::const_iterator p1 = well_index.find(w1);
	  
	  std::map<int, int>::const_iterator p2 = well_index.find(w2);
      
	  if(p1 == well_index.end() || p2 == well_index.end()) {
	    //
	    std::cerr << funame << "no density of states for wells connected with "
	      //
		      << Model::inner_barrier(b).name() << " barrier at "
	      //
		      <<  energy_bin(e) / Phys_const::kcal
	      //
		      << " kcal/mol\n";
	    
	    throw Error::Logic();
	  }

	  km(p1->second, p2->second) -= inner_barrier(b).state_number(e) / 2. / M_PI
	    //
	    / std::sqrt(well(w1).state_density(e) * well(w2).state_density(e));
	}
      }

      // diagonal isomerization contribution
      //
      for(int i = 0; i < well_array.size(); ++i) {
	//
	int w = well_array[i];
	
	if(e < cum_stat_num[w].size())
	  //
	  km(i, i) = cum_stat_num[w][e] / 2. / M_PI / well(w).state_density(e);
	
	if(Model::well(w).escape_size())
	  //
	  km(i, i) += Model::well(w).escape_rate(energy_bin(e));
      }

      // relaxation eigenvalues
      //
      Lapack::Matrix evec(well_array.size());

      Lapack::Vector eval = Mpack::eigenvalues(km, &evec);

      kinetic_basis[e].eigenvalue  = eval;
      
      kinetic_basis[e].eigenvector = evec;
      
      // kinetically active subspace
      //
      for(itemp = 0; itemp < well_array.size(); ++itemp) {
	//
	int w = well_array[itemp];
	
	if(eval[itemp] > well(w).collision_frequency() * reduction_threshold)
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

  IO::log << IO::log_offset << "original kinetic matrix size = " << itemp << std::endl;
  
  IO::log << IO::log_offset << "reduced  kinetic matrix size = " << global_size << std::endl;

  // bimolecular-to-bimolecular rate coefficients
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  Lapack::SymmetricMatrix bb_rate;
  
  if(Model::bimolecular_size()) {
    //
    bb_rate.resize(Model::bimolecular_size());
    
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

    // rate output
    //
    if(!global_size) {
      //
      for(int i = 0; i < Model::bimolecular_size(); ++i) {
	//
	if(Model::bimolecular(i).dummy())
	  //
	  continue;

	dtemp = (Model::bimolecular(i).ground() - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;

	const double fac = Model::bimolecular(i).weight(temperature()) / std::exp(dtemp);
	
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

	  dtemp *= energy_step() / fac / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 
	}
      }
    }
  }

  if(!global_size)
    //
    return;
  
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
  int64_t cycle_size = ener_index_max * ener_index_max;

#pragma omp parallel for default(shared) private(itemp,dtemp) schedule(static)
  
  for(int64_t cycle = 0; cycle < cycle_size; ++cycle) {
    //
    const int e1 = cycle / ener_index_max;

    const int e2 = cycle % ener_index_max;

    if(e2 < e1)
      //
      continue;
  
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
  }// cycle

  /***********************************************************************
   ************************* GLOBAL MATRICES *****************************
   ***********************************************************************/

  Lapack::Vector eigenval;
  
  Lapack::Matrix global_eigen(global_size);

  {
    IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);
    
    if(Mpack::mp_type == Mpack::DOUBLE) {
      //
      eigenval = Offload::eigenvalues(kin_mat, &global_eigen);
    }
    else if(use_mp) {
      //
      eigenval = Mpack::eigenvalues(kin_mat, &global_eigen);
    }
    else
      //
      eigenval = kin_mat.eigenvalues(&global_eigen);
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
	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "well is not in the kinetic basis space\n";

	    throw Error::Logic();
	  }

	  dtemp = 0.;
	
	  for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	    //
	    dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	  eigen_bim(l, p) += dtemp * outer_barrier(b).state_number(e) / 2. / M_PI
	    //
	    * thermal_factor(e) / well(w).boltzman_sqrt(e);
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

 
  IO::log << IO::log_offset << "eigenstate populations:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(Model::log_precision + 7) << "*R"
    //
	  << std::setw(Model::log_precision + 7) << "*P";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
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
	    << std::setw(Model::log_precision + 7) << eigenval[l] / relax_eval_min
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well(w).weight_sqrt() / nfac;
      
      IO::log << std::setw(Model::log_precision + 7);
      
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
	  << std::setw(Model::log_precision + 7) << "---"
    //
	  << std::setw(Model::log_precision + 7) << "---";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    if(!w || dtemp < well(w).weight_sqrt())
      //
      dtemp = well(w).weight_sqrt();
  
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << well(w).weight_sqrt() / dtemp;
  
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
	  << "*Z - well partition function square root (normalized)\n";
  
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

  /*******************************************************************
   ***************** CHEMICAL EIGENSTATES CORRECTION *****************
   *******************************************************************/

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
	      std::cerr << funame << Model::outer_barrier(b).name() << " barrier top below "
		//
			<< Model::well(w2).name() << " well bottom\n";

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
  }// chemical eigenstate cycle
    
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
  if(Model::bimolecular_size()) {
    //
    for(int i = 0; i < Model::bimolecular_size(); ++i)
      //
      for(int j = i; j < Model::bimolecular_size(); ++j)
	//
	bb_rate(i, j) += triple_product(&eigen_bim(chem_size, i), &eigen_bim(chem_size, j), relax_lave, relax_size);

    for(int i = 0; i < Model::bimolecular_size(); ++i) {
      //
      if(Model::bimolecular(i).dummy())
	//
	continue;

      dtemp = (Model::bimolecular(i).ground() - energy_reference()) / temperature();

      if(dtemp < - Limits::exp_pow_max())
	//
	continue;

      const double fac = Model::bimolecular(i).weight(temperature()) / std::exp(dtemp);
      
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

	dtemp *= energy_step() / fac / bru;
	  
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

	dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;

	const double fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
	
	dtemp = vtemp[p] + triple_product(&eigen_bim(chem_size, p), &eigen_escape(chem_size, s), relax_lave, relax_size);

	dtemp *= energy_step() / fac / bru;
	    
	rate_data[std::make_pair(Model::well_size() + p, Model::well_size() + Model::bimolecular_size() + s)] = dtemp; 
      }
    }
  }

  /***** PARTITIONING THE GLOBAL PHASE SPACE INTO THE CHEMICAL AND COLLISIONAL SUBSPACES *****/

  std::vector<int>  group_index;
  
  Lapack::Matrix m_direct;
  
  std::vector<double>      weight;
  
  std::vector<double> real_weight;

  std::vector<double> real_ground;

  Lapack::Matrix m_inverse;
  
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

    double proj = 0.;
      
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
      
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group);

      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }

    group_index = well_partition.group_index();
    weight      = well_partition.weight();
    real_weight = well_partition.real_weight(temperature());
    real_ground = well_partition.real_ground();

    // output
    //
    if(chem_size != Model::well_size()) {
      //
      IO::log << IO::log_offset << "species:\n"
	//
	      << IO::log_offset << std::setw(2) << "#"  << std::setw(Model::log_precision + 7)
	//
	      << "new name" << IO::first_offset
	//
	      << "group\n";
      
      for(int g = 0; g < well_partition.size(); ++g) {
	//
	IO::log << IO::log_offset << std::setw(2) << g << std::setw(Model::log_precision + 7)
	  //
		<< Model::well(group_index[g]).short_name() << IO::first_offset;
	
	for(Group::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w) {
	  //
	  if(w != well_partition[g].begin())
	    //
	    IO::log << "+";
	  
	  IO::log << Model::well(*w).short_name();
	}
	
	IO::log << "\n";
      }
    }

    m_direct = pop_chem;
    
    m_inverse = m_direct.invert();

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
    IO::log << IO::log_offset << "Wa->Wb/Wb->Wa rate constants ratios:\n"
      //
	    << IO::log_offset << std::setw(5) << "Wb\\Wa";
    
    for(int i = 0; i < chem_size; ++i)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::well(group_index[i]).short_name();
    
    IO::log << "\n";
    
    for(int j = 0; j < chem_size; ++j) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(group_index[j]).short_name();
      
      for(int i = 0; i < chem_size; ++i)
	//
	if(i != j) {
	  //
	  if(ww_rate(j, i) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << ww_rate(i, j) / ww_rate(j, i);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	}
	else
	  //
	  IO::log << std::setw(Model::log_precision + 7) << "1";
      
      IO::log << "\n";
    }
    
    if(Model::bimolecular_size()) {
      //
      IO::log << IO::log_offset << "W->P/P->W rate constants ratios:\n"
	//
	      << IO::log_offset << std::setw(5) << "P\\W";
      
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(Model::log_precision + 7) << Model::well(group_index[w]).short_name();
      
      IO::log << "\n";
    
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).short_name();
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  if(bw_rate(p, w) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << wb_rate(w, p) / bw_rate(p, w);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	
	IO::log << "\n";
      }
    }

    //IO::log << std::setprecision(6);

    // well-to-well rate coefficients
    //
    for(int i = 0; i < chem_size; ++i) {
      //
      dtemp = (real_ground[i] - energy_reference()) / temperature();

      if(dtemp < -Limits::exp_pow_max())
	//
	continue;
      
      const double fac = real_weight[i] / std::exp(dtemp);
      
      for(int j = 0; j < chem_size; ++j) {
	//
	dtemp = ww_rate(i, j) * std::sqrt(weight[i] * weight[j]) * energy_step() / fac / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp; 	
      }
    }
    
    // well-to-escape rate coefficients
    //
    if(Model::escape_size())
      //
      for(int w = 0; w < chem_size; ++w) {
	//
	dtemp = (real_ground[w] - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = real_weight[w] / std::exp(dtemp);
      
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  dtemp = vdot(m_inverse.column(w), &eigen_escape(0, e)) * std::sqrt(weight[w]) * energy_step() / fac / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() + Model::bimolecular_size() + e)] = dtemp;
	}
      }
    
    // well-to-bimolecular rate coefficients
    //
    for(int w = 0; w < chem_size; ++w) {
      //
      dtemp = (real_ground[w] - energy_reference()) / temperature();

      if(dtemp < -Limits::exp_pow_max())
	//
	continue;
      
      const double fac = real_weight[w] / std::exp(dtemp);
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = wb_rate(w, p) * std::sqrt(weight[w]) * energy_step() / fac / Phys_const::herz;
	
	rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
      }
    }
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(!Model::bimolecular(p).dummy()) {
	//
	dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	if(dtemp < -Limits::exp_pow_max())
	  //
	  continue;
	
	const double fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
      
	for(int w = 0; w < chem_size; ++w) {
	  //
	  dtemp = bw_rate(p, w) * std::sqrt(weight[w]) * energy_step() / fac / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	}
      }
  }

  // hot energies branching ratios
  //
  if(hot_index.size()) {
    //
    int well_name_size_max;

    std::vector<int> well_name_size(Model::well_size());
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      itemp = Model::well(w).short_name().size() + 1;

      if(!w || itemp > well_name_size_max)
	//
	well_name_size_max = itemp;
    
      if(itemp < Model::log_precision + 7)
	//
	itemp = Model::log_precision + 7;

      well_name_size[w] = itemp;
    }

    if(well_name_size_max < 5)
      //
      well_name_size_max = 5;

    std::vector<int> bim_name_size(Model::bimolecular_size());
    
    for(int p = 0; p < Model::bimolecular_size(); ++p) {
      //
      itemp = Model::bimolecular(p).short_name().size() + 1;

      if(itemp < Model::log_precision + 7)
	//
	itemp = Model::log_precision + 7;

      bim_name_size[p] = itemp;
    }

    IO::log << IO::log_offset << "hot energies branching fractions:\n"
      //
	    << IO::log_offset
      //
	    << std::setw(well_name_size_max)  << "Well"
      //
	    << std::setw(Model::log_precision + 7) << "E, kcal";
    
    for(int w = 0; w < chem_size; ++w)
      //
      IO::log << std::setw(well_name_size[group_index[w]]) << Model::well(group_index[w]).short_name();

    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();

    for(int s = 0; s < Model::escape_size(); ++s)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::escape_name(s);
    
    //IO::log << std::setw(Model::log_precision + 7) << "total" << "\n";

    IO::log << "\n";
      
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
	      << std::setw(well_name_size_max)  << Model::well(w).short_name()
	//
	      << std::setw(Model::log_precision + 7) << energy_bin(e) / Phys_const::kcal;

      double tot = 0.;
      
      // bound species
      //
      for(int c = 0; c < chem_size; ++c) {
	//
	dtemp = vdot(m_direct.row(c), &eigen_hot(0, count)) * std::sqrt(weight[c]);
	
	IO::log << std::setw(well_name_size[group_index[c]]) << dtemp;

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
	
	IO::log << std::setw(bim_name_size[p]) << dtemp;

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
	
	IO::log << std::setw(Model::log_precision + 7) << dtemp;

	tot += dtemp;
      }
      
      //IO::log << std::setw(Model::log_precision + 7) << tot << "\n";

      IO::log << "\n";
    }
  }

  // prompt isomerization
  //
  IO::log << IO::log_offset << "prompt isomerization/dissociation:\n";

  std::vector<int> well_name_size(Model::well_size());

  int well_name_size_max;

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    itemp = Model::well(w).short_name().size() + 1;

    if(!w || itemp > well_name_size_max)
      //
      well_name_size_max = itemp;
    
    well_name_size[w] = itemp > Model::log_precision + 7 ? itemp : Model::log_precision + 7;
  }

  if(well_name_size_max < 7)
    //
    well_name_size_max = 7;

  std::vector<int> bim_name_size(Model::bimolecular_size());
  
  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    itemp = Model::bimolecular(p).short_name().size() + 1;

    bim_name_size[p] = itemp > Model::log_precision + 7 ? itemp : Model::log_precision + 7;
  }

  // output
  //
  IO::log << IO::log_offset << std::setw(well_name_size_max) << "W\\W,P,E";

  if(chem_size)
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      IO::log << std::setw(well_name_size[w]) << Model::well(w).short_name();
    }

  IO::log << std::setw(Model::log_precision + 7) << "Total";

  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    IO::log << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
  
  for(int e = 0; e < Model::escape_size(); ++e)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::escape_name(e);
    
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    dtemp = (Model::well(w).ground() - energy_reference()) / temperature();

    if(dtemp < -Limits::exp_pow_max())
      //
      continue;
    
    const double fac = Model::well(w).weight(temperature()) / std::exp(dtemp);
      
    IO::log << IO::log_offset << std::setw(well_name_size_max) << Model::well(w).short_name();
      
    double diss = 1.;

    if(chem_size) {
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
      
	dtemp *= well(v).weight_sqrt() * well(w).weight_sqrt() * energy_step() / fac;

	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	diss -= dtemp;

	IO::log << std::setw(well_name_size[v]) << dtemp;
      }
    }
  
    IO::log << std::setw(Model::log_precision + 7) << diss;

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

      IO::log << std::setw(bim_name_size[p]) << dtemp	* energy_step() / fac;
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
      
      IO::log << std::setw(Model::log_precision + 7) << dtemp * energy_step() / fac;
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
    
    evec_out << "\t Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";

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
	evec_out << std::setw(13) << Model::well(w).short_name() << std::setw(13) << Model::well(w).short_name();
      
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
    if(!Model::bimolecular_size()) {
      //
      std::cerr << funame << "no bimolecular products\n";

      throw Error::Logic();
    }
    
    switch(pressure_unit) {
      //
    case BAR:
      //
      ped_out << "pressure[bar]        = " << pressure() / Phys_const::bar << "\n";
      break;
      
    case TORR:
      //
      ped_out << "pressure[torr]       = " << pressure() / Phys_const::tor << "\n";
      break;
      
    case ATM:
      //
      ped_out << "pressure[atm]        = " << pressure() / Phys_const::atm << "\n";
      break;
    }
    
    ped_out << "temperature[K]       = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << "\n"
      //
	    << "energy step[1/cm]    = " << energy_step() / Phys_const::incm << "\n"
      //
	    << "max energy[kcal/mol] = " << energy_reference() / Phys_const::kcal << "\n\n";


    std::vector<int> bim_name_size(Model::bimolecular_size());
    
    for(int p = 0; p < Model::bimolecular_size(); ++p) {
      //
      itemp = Model::bimolecular(p).short_name().size() + 1;

      bim_name_size[p] = itemp > Model::ped_precision + 7 ? itemp : Model::ped_precision + 7;
    }

    // bimolecular-to-bimolecular product energy distributions
    //
    if(ped_pair.size()) {
      //
      ped_out << "bimolecular PEDs:\n";

      // dimensions
      //
      itemp = 0;
      
      for(int p = 0; p < ped_pair.size(); ++p)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(Model::outer_connect(b).second == ped_pair[p].second && outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

      Lapack::Matrix ped(itemp, ped_pair.size());

      ped = 0.;
      
      for(std::vector<std::pair<int,int> >::const_iterator pi = ped_pair.begin(); pi != ped_pair.end(); ++pi) {
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	  //
	  if(Model::outer_connect(b).second != pi->second)
	    //
	    continue;
	    
	  for(int e = 0; e < outer_barrier(b).size(); ++e) {
	    //
	    std::map<int, int>::const_iterator bit = kinetic_basis[e].well_index_map.find(Model::outer_connect(b).first);

	    if(bit == kinetic_basis[e].well_index_map.end()) {
	      //
	      std::cerr << funame << "bimolecular PEDs: outer well is not in the kinetic basis space\n";

	      throw Error::Logic();
	    }
	    
	    double x = 0.;

	    // relaxational eigenstates contribution
	    //
	    for(int l = chem_size; l < global_size; ++l) {
	      //
	      dtemp = 0.;

	      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
		//
		dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(bit->second, k);

	      x += dtemp * eigen_bim(l, pi->first) / eigenval[l];
	    }

	    // high eigenstates contribution
	    //
	    for(int c = 0; c < Model::outer_barrier_size(); ++c) {
	      //
	      if(Model::outer_connect(c).second != pi->first || e >= outer_barrier(c).size())
		//
		continue;
		 
	      std::map<int, int>::const_iterator cit = kinetic_basis[e].well_index_map.find(Model::outer_connect(c).first);

	      if(cit == kinetic_basis[e].well_index_map.end()) {
		//
		std::cerr << funame << "bimolecular PEDs: inner well is not in the kinetic basis space\n";

		throw Error::Logic();
	      }
	    
	      dtemp = 0.;
	      
	      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
		//
		dtemp += kinetic_basis[e].eigenvector(bit->second, k)
		  //
		  * kinetic_basis[e].eigenvector(cit->second, k)
		  //
		  / kinetic_basis[e].eigenvalue[k];

	      x += dtemp * outer_barrier(c).state_number(e) / 2. / M_PI
		//
		* thermal_factor(e) / well(cit->first).boltzman_sqrt(e);
	    }

	    ped(e, pi - ped_pair.begin()) += x * outer_barrier(b).state_number(e) / 2. / M_PI
	      //
	      * thermal_factor(e) / well(bit->first).boltzman_sqrt(e);
	    //
	  }// final energy cycle
	  //
	}// final barrier cycle
	//
      }// PED pair cycle
      
      // output
      //
      std::vector<int> ped_name_size(ped_pair.size());
      
      ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";
      
      for(int p = 0; p < ped_pair.size(); ++p) {
	//
	stemp = " " + Model::bimolecular(ped_pair[p].first).short_name() + "->"
	  //
	  + Model::bimolecular(ped_pair[p].second).short_name();

	itemp = Model::ped_precision + 7;
	
	ped_name_size[p] = stemp.size() > itemp ? stemp.size() : itemp;
	
	ped_out << std::setw(ped_name_size[p]) << stemp;
      }
      
      ped_out << "\n";

      for(int e = 0; e < ped.size1(); ++e) {
	//
	ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	for(int p = 0; p < ped.size2(); ++p)
	  //
	  ped_out << std::setw(ped_name_size[p]) << ped(e, p);
	
	ped_out << "\n";
      }
      
      ped_out << "\n";
      //
    }// bimolecular PEDs
    
    // hot energies PEDs
    //
    if(hot_index.size()) {
      //
      ped_out << "hot energies PEDs:\n";
    
      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(!b || outer_barrier(b).size() > itemp)
	  //
	  itemp = outer_barrier(b).size();

      Lapack::Matrix ped(itemp, Model::bimolecular_size());

      // hot energy cycle
      //
      int count = 0;
	
      for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
	//
	const int& hw = hit->first;

	const int& he = hit->second;

	ped = 0.;
	  
	std::map<int, int>::const_iterator hi = kinetic_basis[he].well_index_map.find(hw);

	if(hi == kinetic_basis[he].well_index_map.end()) {
	  //
	  std::cerr << funame << "hot energies PEDs: well is not in the kinetic basis space\n";

	  throw Error::Logic();
	}
      
	// high eigenstates contribution
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {// outer barrier cycle
	  //
	  if(he >= outer_barrier(b).size())
	    //
	    continue;

	  const int& w = Model::outer_connect(b).first;
	    
	  const int& p = Model::outer_connect(b).second;
	
	  std::map<int, int>::const_iterator i = kinetic_basis[he].well_index_map.find(w);

	  if(i == kinetic_basis[he].well_index_map.end()) {
	    //
	    std::cerr << funame << "hot energies PEDs: high eigenstates: well is not in the kinetic basis space\n";

	    throw Error::Logic();
	  }
      
	  dtemp = 0.;
	
	  for(int k = kinetic_basis[he].active_size; k < kinetic_basis[he].size(); ++k)
	    //
	    dtemp += kinetic_basis[he].eigenvector(hi->second, k)
	      //
	      * kinetic_basis[he].eigenvector(i->second, k)
	      //
	      / kinetic_basis[he].eigenvalue[k];

	  ped(he, p) += dtemp / well(hw).boltzman_sqrt(he)
	    //
	    * outer_barrier(b).state_number(he) / 2. / M_PI
	    //
	    * thermal_factor(he) / well(w).boltzman_sqrt(he);
	  //
	}// outer barrier cycle

	// relaxational eigenstates contribution
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {// outer barrier cycle
	  //
	  const int& w = Model::outer_connect(b).first;
      
	  const int& p = Model::outer_connect(b).second;

	  for(int e = 0; e < outer_barrier(b).size(); ++e) {// energy cycle
	    //
	    std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	    if(i == kinetic_basis[e].well_index_map.end()) {
	      //
	      std::cerr << funame << "hot energies PEDs: relaxational eigenstates: well is not in the kinetic basis space\n";

	      throw Error::Logic();
	    }

	    double x = 0.;
	      
	    for(int l = chem_size; l < global_size; ++l) {
	      //
	      dtemp = 0.;
	
	      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
		//
		dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	      x += dtemp * eigen_hot(l, count) / eigenval[l];
	    }

	    ped(e, p) += x * outer_barrier(b).state_number(e) / 2. / M_PI
	      //
	      * thermal_factor(e) / well(w).boltzman_sqrt(e);
	    // 
	  }// energy cycle
	  //
	}// outer barrier cycle
	//
	// output
	//
	ped_out << "well: " << Model::well(hw).short_name()
	  //
		<< "   hot energy = " << energy_bin(he) / Phys_const::kcal
	  //
		<< " kcal/mol\n";

	ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";
	  
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  ped_out << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();

	ped_out << "\n";
      
	for(int e = 0; e < ped.size1(); ++e) {
	  //
	  ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(bim_name_size[p]) << ped(e, p);

	  ped_out << "\n";
	}

	ped_out << "\n";
	//
      }// hot energies cycle
      //
    }// hot energies PEDs
    
    // well-to-bimolecular PEDs
    //
    if(chem_size) {
      //
      ped_out << "well PEDs:\n";

      //dimensions
      //
      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(!b || outer_barrier(b).size() > itemp)
	  //
	  itemp = outer_barrier(b).size();

      Lapack::Matrix ped(itemp, Model::bimolecular_size());

      // chemical species cycle
      //
      for(int c = 0; c < chem_size; ++c) {
	//
	ped = 0.;

	for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	  //
	  const int& w = Model::outer_connect(b).first;

	  const int& p = Model::outer_connect(b).second;

	  for(int e = 0; e < outer_barrier(b).size(); ++e) {
	    //
	    std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	    if(i == kinetic_basis[e].well_index_map.end()) {
	      //
	      std::cerr << funame << "well PEDs: well is not in the kinetic basis space\n";

	      throw Error::Logic();
	    }

	    double x = 0.;
	    
	    for(int l = 0; l < chem_size; ++l) {
	    
	      dtemp = 0.;

	      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
		//
		dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	      x += m_inverse(l, c) * dtemp;
	    }

	    ped(e, p) += x * outer_barrier(b).state_number(e) / 2. / M_PI
	      //
	      * thermal_factor(e) / well(w).boltzman_sqrt(e);
	    //
	  }// energy cycle
	  //
	}// outer barrier cycle

	// high eigenstate correction
	//
	for(int l = 0; l < chem_size; ++l) {
	  //
	  for(int w1 = 0; w1 < Model::well_size(); ++w1) {
	    //
	    for(int e1 = 0; e1 < well(w1).size(); ++e1) {
	      //
	      std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);
	  
	      if(i11 == kinetic_basis[e1].well_index_map.end()) {
		//
		std::cerr << funame << "well PEDs: high eigenstate correction: i11 well is not in the kinetic basis space\n";

		throw Error::Logic();
	      }

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
		  std::cerr << funame << "well PEDs: high eigenstate correction: i21 well is not in the kinetic basis space\n";

		  throw Error::Logic();
		}
		
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
		    std::cerr << funame << "well PEDs: high eigenstate correction: i22 well is not in the kinetic basis space\n";

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
		    * well(w1).boltzman_sqrt(e1) / well(w1).boltzman_sqrt(e2)
		    //
		    * m_inverse(l, c);
	    
		  ped(e2, p) -= dtemp * outer_barrier(b).state_number(e2) / 2. / M_PI
		    //
		    * thermal_factor(e2) / well(w2).boltzman_sqrt(e2);
		  //
		}// outer barrier cycle
		//
	      }// second energy cycle
	      //
	    }// first energy cycle
	    //
	  }// first well cycle
	  //
	}// chemical eigenstate cycle

	// output
	//
	ped_out << "well: " << Model::well(group_index[c]).short_name() << "\n";
	
	ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  ped_out << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
	  
	ped_out << "\n";

	for(int e = 0; e < ped.size1(); ++e) {
	  //
	  ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(bim_name_size[p]) << ped(e, p);
	    
	  ped_out << "\n";
	}
	  
	ped_out << "\n";
	
      }// chemical species cycle
      //
    }// well PEDs
    //
  }// product energy distributions

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

  std::time_t start_time = std::time(0);

  IO::log << IO::log_offset << "time = 0\n";

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
    //
	  << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";

  IO::aux << "Pressure = " << pressure() / Phys_const::bar << " bar"
    //
	  << "\t Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";
  
  // collisional frequency output
  //
  IO::log << IO::log_offset
	  << std::setw(Model::log_precision + 7) << "Well"
	  << std::setw(20) << "Collision, 1/sec"
	  << "\n";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << IO::log_offset
	    << std::setw(Model::log_precision + 7) << Model::well(w).short_name()
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
    
    evec_out << "\t Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";
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
      
#pragma omp parallel for default(shared) schedule(static)

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
#pragma omp parallel for default(shared) schedule(static)

      for(int i = 0; i < cum_stat_num[w].size(); ++i)
	//
	kin_mat(i + well_shift[w], i + well_shift[w]) = cum_stat_num[w][i] / 2. / M_PI
	  //
	  / well(w).state_density(i);
      
      if(Model::well(w).escape_size())
	//
#pragma omp parallel for default(shared) schedule(static)

	for(int i = 0; i < well(w).size(); ++i)
	  //
	  kin_mat(i + well_shift[w], i + well_shift[w]) += Model::well(w).escape_rate(energy_bin(i));
    }

    // collision relaxation contribution
    //
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      const int64_t cycle_size = (int64_t)well(w).size() * (int64_t)well(w).size();
	
#pragma omp parallel for default(shared) schedule(static)

      for(int64_t cycle = 0; cycle < cycle_size; ++cycle) {
	//
	int i = cycle / well(w).size();

	int j = cycle % well(w).size();

	if(j < i)
	  //
	  continue;
	
	if(j - i >= well(w).kernel_bandwidth)
	  //
	  continue;
	  
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

	int64_t cycle_size = (int64_t)well(w).size() * (int64_t)well(w).size();
	
#pragma omp parallel for default(shared) schedule(static)

	for(int64_t cycle = 0; cycle < cycle_size; ++cycle) {
	  //
	  int i = cycle / well(w).size();

	  int j = cycle % well(w).size();

	  if(j < i)
	    //
	    continue;
	  
	  kin_mat(i + well_shift[w], j + well_shift[w]) +=  well(w).radiation_rate(i, j); 
	}
      }

    // bimolecular product vectors
    //
    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      const int w = Model::outer_connect(b).first;
      
      const int p = Model::outer_connect(b).second;

#pragma omp parallel for default(shared) schedule(static)

      for(int e = 0; e < outer_barrier(b).size(); ++e)
	//
	global_bim(e + well_shift[w], p) = outer_barrier(b).state_number(e) / 2. / M_PI
	  //
	  * thermal_factor(e) / well(w).boltzman_sqrt(e);
    }
  
    // thermal distributions
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
#pragma omp parallel for default(shared) schedule(static)

      for(int e = 0; e < well(w).size(); ++e)
	//
	global_pop(e + well_shift[w], w) = well(w).boltzman_sqrt(e) / well(w).weight_sqrt();
    
    // well escape
    //
    for(int e = 0; e < Model::escape_size(); ++e) {
      //
      const int ew = Model::escape_channel(e).first;

      const int ei = Model::escape_channel(e).second;
      
#pragma omp parallel for default(shared) schedule(static)

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
    IO::log << IO::log_offset << "time = " << std::time(0) - start_time << "\n";

    //IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);

    IO::Marker solve_marker("diagonalizing kinetic matrix");

    if(save_kinetic_matrix.size()) {
      //
      IO::log << IO::log_offset << "saving kinetic matrix ... ";
      
      std::ofstream to(save_kinetic_matrix.c_str(), std::ios_base::app);

      to.precision(3);
      
      to << "Pressure = ";
    
      switch(pressure_unit) {
	//
      case BAR:
	//
	to << pressure() / Phys_const::bar << " bar";
	
	break;
	
      case TORR:
	//
	to << pressure() / Phys_const::tor << " torr";
      
	break;
      
      case ATM:
	//
	to << pressure() / Phys_const::atm << " atm";
      
	break;
      }
    
      to << "   Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";

      to << kin_mat.size() << "\n";

      to.precision(14);
      
      for(int i = 0; i < kin_mat.size(); ++i)
	//
	for(int j = i; j < kin_mat.size(); ++j)
	  //
	  if(kin_mat(i, j) != 0.) {
	    //
	    to << std::setw(5) <<  i << " " << std::setw(5) << j << " " << kin_mat(i, j) << "\n";
	  }

      to << "\n";

      IO::log << "done\n";
    }
    
    if(Mpack::mp_type == Mpack::DOUBLE) {
      //
      eigenval = Offload::eigenvalues(kin_mat, &eigen_global);
    }
    else if(use_mp) {
      //
      eigenval = Mpack::eigenvalues(kin_mat, &eigen_global);
    }
    else
      //
      eigenval = kin_mat.eigenvalues(&eigen_global);

    eigen_global = eigen_global.transpose();
  }

  const double min_relax_eval = eigenval[Model::well_size()];
  
  const double max_relax_eval = eigenval.back();

  //IO::log << IO::log_offset << "minimal relaxation eigenvalue / collision frequency = "
  //<< min_relax_eval / collision_frequency() << "\n";
  //IO::log << IO::log_offset << "maximal relaxation eigenvalue / collision frequency = "
  //<< max_relax_eval / collision_frequency() << "\n";

  // microcanonical rate coefficients at highest energy
  
  IO::log << IO::log_offset << "time = " << std::time(0) - start_time << "\n";

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
      
      IO::log << std::setw(Model::log_precision + 7) << Model::well(w1).short_name() + "<->" + Model::well(w2).short_name();
    }
    
    IO::log << "\n" << IO::log_offset;

    for(int b = 0; b < Model::inner_barrier_size(); ++b) {
      //
      int w1 = Model::inner_connect(b).first;
      
      int w2 = Model::inner_connect(b).second;
      
      dtemp = well(w1).state_density(0) < well(w2).state_density(0) ?
					  //
					  well(w1).state_density(0) : well(w2).state_density(0);
      
      IO::log << std::setw(Model::log_precision + 7) << inner_barrier(b).state_number(0) / 2. / M_PI
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
      
      IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name() + "->" + Model::bimolecular(p).short_name();
    }
    
    IO::log << "\n" << IO::log_offset;

    for(int b = 0; b < Model::outer_barrier_size(); ++b) {
      //
      int w = Model::outer_connect(b).first;
      
      IO::log << std::setw(Model::log_precision + 7) << outer_barrier(b).state_number(0) / 2. / M_PI
	//
	/ well(w).state_density(0) / well(w).collision_frequency();
    }
    
    IO::log << "\n";
    //
  }// outer barriers

  /*
  // low eigenvalue method
  //
  if(eigenval[0] / min_relax_eval < chemical_tolerance) {
  //
  IO::log << IO::log_offset << "some eigenvalues are too small: using low eigenvalue method\n";
    
  IO::Marker low_eval_marker("low eigenvalue method");
    
  Lapack::SymmetricMatrix k_11;
  Lapack::SymmetricMatrix k_33;
  Lapack::Matrix k_13;
  Lapack::Matrix l_21;

  low_eigenvalue_matrix(k_11, k_33, k_13, l_21);

  Lapack::Matrix chem_evec(Model::well_size());
    
  Lapack::Vector chem_eval = Mpack::eigenvalues(k_11, &chem_evec);
    
  // low-eigenvalue chemical subspace
  //
  itemp = 1;
    
  while(eigenval[itemp] / min_relax_eval < chemical_tolerance) { ++itemp; }
    
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
      
  well_partition_method(pop_chem, low_well_partition, bimolecular_group);

  IO::log << IO::log_offset << "low eigenvalue bound species:";
      
  for(int g = 0; g < low_well_partition.size(); ++g) {
  //
  IO::log << " ";
	
  for(Group::const_iterator w = low_well_partition[g].begin(); w != low_well_partition[g].end(); ++w) {
  //
  if(w != low_well_partition[g].begin())
  //
  IO::log << "+";
	  
  IO::log << Model::well(*w).short_name();
  }
  }
      
  IO::log << "\n" << IO::log_offset << "low eigenvalue bimolecular group:";
      
  for(Group::const_iterator w = bimolecular_group.begin(); w != bimolecular_group.end(); ++w)
  //
  IO::log << " " << Model::well(*w).short_name();
      
  IO::log << "\n";

  Lapack::SymmetricMatrix low_km(chem_size);
      
  low_km = 0.;
      
  for(int i = 0; i < chem_size; ++i)\
  //
  for(int j = i + 1; j < chem_size; ++j) {
  //
  dtemp = 0.;
	  
  for(Group::const_iterator ig = low_well_partition[i].begin(); ig != low_well_partition[i].end(); ++ig)
  //
  for(Group::const_iterator jg = low_well_partition[j].begin(); jg != low_well_partition[j].end(); ++jg)
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

  for(Group::const_iterator ig = low_well_partition[i].begin(); ig != low_well_partition[i].end(); ++ig) {
  //
  dtemp = 0.;
	  
  // bimolecular channel contribution
  //
  for(int p = 0; p < Model::bimolecular_size(); ++p)
  //
  dtemp += k_13(*ig, p);
	  
  // bimolecular group contribution
  //
  for(Group::const_iterator jg = bimolecular_group.begin(); jg != bimolecular_group.end(); ++jg)
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
	
  if(reduction_method == PROJECTION) {
  chem_eval[l] = low_eval[l];
  for(int w = 0; w < Model::well_size(); ++w)
  chem_evec(w, l) = low_evec.column(l) * basis.row(w);
  }
	
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
  << std::setw(Model::log_precision + 7) << "DD eval"
  //
  << std::setw(Model::log_precision + 7) << "LE eval"
  //
  << std::setw(Model::log_precision + 7) << "LE proj\n";

  for(int l = 0; l < Model::well_size(); ++l) {
  //
  IO::log << IO::log_offset <<  std::setw(5) << l
  //
  << std::setw(Model::log_precision + 7) << eigenval[l] / min_relax_eval
  //
  << std::setw(Model::log_precision + 7) << chem_eval[l] / min_relax_eval
  //
  << std::setw(Model::log_precision + 7) << rel_proj[l];
      
  if(eigenval[l] / min_relax_eval >= chemical_tolerance)
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
  */
    
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
    
    Model::time_evolution->out << "\t Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n\n";

    Model::time_evolution->out << std::setw(13) << "time, sec";
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      Model::time_evolution->out << std::setw(13) << Model::well(w).short_name();
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      Model::time_evolution->out << std::setw(13) << Model::bimolecular(p).short_name();

    Model::time_evolution->out << "\n";
    
    // bimolecular reactants
    //
    if(Model::bimolecular_size() && Model::time_evolution->excess_reactant_concentration() > 0.) {
      //
      if(Model::bimolecular(react).dummy()) {
	//
	std::cerr << funame << "time evolution: reactants partition function not defined\n";
	
	throw Error::Init();
      }

      dtemp = (Model::bimolecular(react).ground() - energy_reference()) / temperature();

      if(dtemp > -Limits::exp_pow_max()) {
	//
	// normaziation factor
	//
	const double nfac = Model::time_evolution->excess_reactant_concentration() * energy_step()
	  //
	  / Model::bimolecular(react).weight(temperature()) * std::exp(dtemp);

	double time_val = Model::time_evolution->start();
      
	for(int t = 0;  t < Model::time_evolution->size(); ++t, time_val *= Model::time_evolution->step()) {
	  //
	  std::vector<double> well_pop(Model::well_size());
	
	  std::vector<double> bim_pop(Model::bimolecular_size());

	  for(int l = 0; l < global_size; ++l) {
	    //
	    dtemp = eigenval[l] * time_val;
	  
	    if(dtemp > 50.) {
	      //
	      dtemp = 1. / eigenval[l];
	    }
	    else
	      //
	      dtemp = (1. - std::exp(-dtemp)) /eigenval[l];
	
	    for(int w = 0; w < Model::well_size(); ++w)
	      //
	      well_pop[w] += eigen_bim(l, react) * eigen_pop(l, w) * dtemp;

	    dtemp = (time_val - dtemp) / eigenval[l];
	    //
	    for(int p = 0; p < Model::bimolecular_size(); ++p)
	      //
	      bim_pop[p] += eigen_bim(l, react) * eigen_bim(l, p) * dtemp;
	  }

	  // normalization
	  //
	  for(int w = 0; w < Model::well_size(); ++w)
	    //
	    well_pop[w] *= well(w).weight_sqrt() * nfac;
	
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    bim_pop[p] *= nfac;
	  
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
	  //
	}// time cycle
	//
      }//
      //
    }// bimolecular reactants
    //
    // bound reactant
    //
    else {
      //
      const double bfac = std::exp(energy_step() / Model::time_evolution->temperature() - energy_step() / temperature());
    
      // initial distribution
      //
      Lapack::Vector init_dist(well(react).size());
      
      dtemp = 1.;
      
      double norm_fac = 0.;
      
      for(int e = 0; e < well(react).size(); ++e, dtemp *= bfac) {
	//
	init_dist[e] = well(react).boltzman_sqrt(e) * dtemp;
	
	norm_fac    += well(react).boltzman(e) * dtemp;
      }

      init_dist /= norm_fac;
      
      std::vector<double> init_coef(global_size);
      
      for(int l = 0; l < global_size; ++l)
	//
	init_coef[l] = parallel_vdot(init_dist, &eigen_global(l, well_shift[react]), well(react).size(), 1, global_size);

      double time_val = Model::time_evolution->start();
      
      for(int t = 0;  t < Model::time_evolution->size(); ++t, time_val *= Model::time_evolution->step()) {
	//
	std::vector<double> well_pop(Model::well_size());
	
	std::vector<double> bim_pop(Model::bimolecular_size());

	for(int l = 0; l < global_size; ++l) {
	  //
	  dtemp = eigenval[l] * time_val;
	  
	  if(dtemp > 100.) {
	    //
	    dtemp = 0.;
	  }
	  else
	    //
	    dtemp = std::exp(-dtemp);
	
	  for(int w = 0; w < Model::well_size(); ++w)
	    //
	    well_pop[w] += init_coef[l] * eigen_pop(l, w) * dtemp;

	  dtemp = (1. - dtemp) / eigenval[l];
	  //
	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
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

  IO::log << IO::log_offset << "time = " << std::time(0) - start_time << "\n";

  IO::log  << IO::log_offset << "eigenvector populations normalized:\n"
    //
	   << IO::log_offset
    //
	   << std::setw(5)  << "L"
    //
	   << std::setw(Model::log_precision + 7) << "*R"
    //
	   << std::setw(Model::log_precision + 7) << "*P";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
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
	    << std::setw(Model::log_precision + 7) << eigenval[l] / min_relax_eval
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well(w).weight_sqrt() / max_pop;
      
      IO::log << std::setw(Model::log_precision + 7);
      
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
	  << IO::log_offset << std::setw(5)  << "L" << std::setw(Model::log_precision + 7) << "*Q" << std::setw(Model::log_precision + 7) << "*P";

  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
  IO::log << "\n";

  for(int l = 0; l < Model::well_size(); ++l) {
    //
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(Model::log_precision + 7) << eigenval[l] / well(0).collision_frequency()
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << eigen_pop(l,w);
    
    IO::log << "\n";
  }
  
  IO::log << IO::log_offset
    //
	  << std::setw(5) << "*Z"
    //
	  << std::setw(Model::log_precision + 7) << "---"
    //
	  << std::setw(Model::log_precision + 7) << "---";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << well(w).weight_sqrt();
  
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
	  << "*Z - well partition function square root\n";
  
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
	evec_out << std::setw(13) << Model::well(w).short_name();
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << Model::well(w).short_name();
      
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

  IO::log << IO::log_offset << "time = " << std::time(0) - start_time << "\n";

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

  /*
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
  */

  // kappa matrix
  //
  if(Model::bimolecular_size()) {
    //
    //Lapack::Matrix kappa = proj_pop.transpose() * inv_proj_bim;

    IO::log << IO::log_offset << "time = " << std::time(0) - start_time << "\n";

    IO::log << IO::log_offset << "isomers-to-bimolecular equilibrium coefficients (kappa matrix):\n"
      //
	    << IO::log_offset << std::setw(5) << "W\\P";
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(p).short_name();
    
    IO::log << "\n";
   
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(w).short_name();
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w), relax_lave, relax_size) 
	  / well(w).weight_sqrt();
	
	//dtemp = kappa(w, p) / well(w).weight_sqrt();
	
	IO::log << std::setw(Model::log_precision + 7);
	
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

    //    IO::log << std::setprecision(6);
  }

  // bimolecular-to-bimolecular rate coefficients
  //
  if(Model::bimolecular_size()) {
    //
    //Lapack::SymmetricMatrix bb_rate = Lapack::SymmetricMatrix(proj_bim.transpose() * inv_proj_bim);
    
    Lapack::SymmetricMatrix bb_rate(Model::bimolecular_size());

    for(int i = 0; i < Model::bimolecular_size(); ++i)
      //
      for(int j = i; j < Model::bimolecular_size(); ++j)
	//
      	bb_rate(i, j) = triple_product(&eigen_bim(chem_size, i), &eigen_bim(chem_size, j), relax_lave, relax_size);

    for(int i = 0; i < Model::bimolecular_size(); ++i) {
      //
      if(Model::bimolecular(i).dummy())
	//
	continue;

      dtemp = (Model::bimolecular(i).ground() - energy_reference()) / temperature();

      if(dtemp < -Limits::exp_pow_max())
	//
	continue;

      const double fac = Model::bimolecular(i).weight(temperature()) / std::exp(dtemp);
    
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

	dtemp *= energy_step() / fac / bru;
	  
	rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 
      }
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
	if(!Model::bimolecular(p).dummy()) {
	  //
	  dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	  if(dtemp > -Limits::exp_pow_max()) {
	    //
	    const double fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
	    
	    for(int e = 0; e < Model::escape_size(); ++e) {
	      //
	      dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_escape(chem_size, e),
				     //
				     relax_lave, relax_size) * energy_step() / fac / bru;
	      
	      //dtemp = escape_bim(e, p) * energy_step() / bimolecular(p).weight() / bru;
	    
	      rate_data[std::make_pair(Model::well_size() + p, Model::well_size() +
				       //
				       Model::bimolecular_size() + e)] = dtemp; 
	    }
	  }
	}
    }
  }

  // product energy distributions
  //
  if(ped_out.is_open()) {
    //
    switch(pressure_unit) {
      //
    case BAR:
      //
      ped_out << "pressure[bar]        = " << pressure() / Phys_const::bar << "\n";
      break;
      
    case TORR:
      //
      ped_out << "pressure[torr]       = " << pressure() / Phys_const::tor << "\n";
      break;
      
    case ATM:
      //
      ped_out << "pressure[atm]        = " << pressure() / Phys_const::atm << "\n";
      break;
    }
    
    ped_out << "temperature[K]           = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << "\n"
	    << "energy step[1/cm]        = " << energy_step() / Phys_const::incm << "\n"
	    << "maximum energy[kcal/mol] = " << energy_reference() / Phys_const::kcal << "\n\n";

    if(!Model::bimolecular_size()) {
      //
      std::cerr << funame << "no bimolecular products\n";

      throw Error::Logic();
    }
    
    int ener_index_max;

    // bimolecular-to-bimolecular PEDs
    //
    if(ped_pair.size()) {
      //
      ped_out << "bimolecular PEDs:\n";

      // dimensions
      //
      itemp = 0;
      
      for(int pi = 0; pi < ped_pair.size(); ++ pi)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(Model::outer_connect(b).second == ped_pair[pi].second && outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

      Lapack::Matrix ped(itemp, ped_pair.size());

      ped = 0.;
      
      // PED
      //
      for(std::vector<std::pair<int, int> >::const_iterator pi = ped_pair.begin(); pi != ped_pair.end(); ++pi)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {// outer barrier cycle
	  //
	  const int& w = Model::outer_connect(b).first;
	  
	  const int& p = Model::outer_connect(b).second;
	    
	  if(p != pi->second)
	    //
	    continue;
	    
	  for(int e = 0; e < outer_barrier(b).size(); ++e)
	    //
	    ped(e, pi - ped_pair.begin()) += global_bim(e + well_shift[w], p) *
	      //
	      triple_product(&eigen_bim(chem_size, pi->first),
			     //
			     &eigen_global(chem_size, e + well_shift[w]),
			     //
			     relax_lave, relax_size);
	
	}// outer barrier cycle
      
      // output
      //
      ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";
      
      std::vector<int> ped_name_size(ped_pair.size());
      
      for(int pi = 0; pi < ped_pair.size(); ++ pi) {
	//
	stemp = " " + Model::bimolecular(ped_pair[pi].first).short_name() + "->"
	  //
	  + Model::bimolecular(ped_pair[pi].second).short_name();

	itemp = Model::ped_precision + 7;
	
	ped_name_size[pi] = stemp.size() > itemp ? stemp.size() : itemp;
	
	ped_out << std::setw(ped_name_size[pi]) << stemp;
      }
      
      ped_out << "\n";

      for(int e = 0; e < ped.size1(); ++e) {
	//
	ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	for(int pi = 0; pi < ped_pair.size(); ++pi)
	  //
	  ped_out << std::setw(ped_name_size[pi]) << ped(e, pi);
	
	ped_out << "\n";
      }
      
      ped_out << "\n";
      //
    }// bimolecular PEDs
    
    // escape PEDs
    //
    if(Model::escape_size()) {
      //
      // hot energies-to-escape PEDs
      //
      if(hot_index.size()) {
	//
	ped_out << "hot-to-escape PEDs:\n";

	std::vector<int> escape_name_size(Model::escape_size());

	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  itemp = Model::ped_precision + 7;

	  escape_name_size[s] = itemp > Model::escape_name(s).size() + 1 ? itemp : Model::escape_name(s).size() + 1;
	}

	// dimensions
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  const int& w = Model::escape_channel(s).first;
	  
	  if(!s || itemp < well(w).size())
	    //
	    itemp = well(w).size();
	}
	
	Lapack::Matrix ped(itemp, Model::escape_size());

	// hot energies cycle
	//
	for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit) {
	  //
	  const int hw = hit->first;
	      
	  const int he = hit->second;

	  ped = 0.;

	  // PED
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    const int& w = Model::escape_channel(s).first;
	    
	    const int& c = Model::escape_channel(s).second;

	    for(int e = 0; e < well(w).size(); ++e)
	      //
	      ped(e, s) = triple_product(&eigen_global(chem_size, he + well_shift[hw]),
					 //
					 &eigen_global(chem_size, e + well_shift[w]),
					 //
					 relax_lave, relax_size)
		//
		* Model::well(w).escape_rate(energy_bin(e), c) * well(w).boltzman_sqrt(e);
	  }
	  
	  // output
	  //
	  ped_out << "well: " << Model::well(hw).short_name() << "  hot energy = "
	    //
		  << energy_bin(he) / Phys_const::kcal << " kcal/mol\n";

	  ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	  for(int s = 0; s < Model::escape_size(); ++s)
	    //
	    ped_out << std::setw(escape_name_size[s]) << Model::escape_name(s);

	  ped_out << "\n";
	
	  for(int e = 0; e < ped.size1(); ++e) {
	    //
	    ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	    for(int s = 0; s < Model::escape_size(); ++s)
	      //
	      ped_out << std::setw(escape_name_size[s]) << ped(e, s);

	    ped_out << "\n";
	  }

	  ped_out << "\n";
	  //
	}// hot energies cycle
	//
      }// hot energies-to-escape PEDs
	
      // bimolecular-to-escape PEDs
      //
      ped_out << "bimolecular-to-escape PEDs:\n";

      // dimensions
      //
      for(int s = 0; s < Model::escape_size(); ++s) {
	//
	const int& w  = Model::escape_channel(s).first;
	  
	if(!s || well(w).size() > itemp)
	  //
	  itemp = well(w).size();
      }
	
      Lapack::Matrix ped(itemp, Model::bimolecular_size() * Model::escape_size());

      ped = 0.;

      // PED
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  const int& w = Model::escape_channel(s).first;

	  const int& c = Model::escape_channel(s).second;

	  for(int e = 0; e < well(w).size(); ++e) {
	    //
	    itemp = s + Model::escape_size() * p;
	    
	    ped(e, itemp) = triple_product(&eigen_bim(chem_size, p),
					   //
					   &eigen_global(chem_size, e + well_shift[w]),
					   //
					   relax_lave, relax_size)
	      //
	      * Model::well(w).escape_rate(energy_bin(e), c) * well(w).boltzman_sqrt(e);
	  }
	}

      // output
      //
      ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

      std::vector<int> ped_name_size(Model::bimolecular_size() * Model::escape_size());
	
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  stemp = " " + Model::bimolecular(p).short_name() + "->" + Model::escape_name(s);

	  itemp = Model::ped_precision + 7;
	      
	  itemp = stemp.size() > itemp ? stemp.size() : itemp;
	      
	  ped_name_size[s + Model::escape_size() * p] = itemp;
	    
	  ped_out << std::setw(itemp) << stemp;
	}
	    
      ped_out << "\n";

      for(int e = 0; e < ped.size1(); ++e) {
	//
	ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	for(int i = 0; i < ped.size2(); ++i)
	  //
	  ped_out << std::setw(ped_name_size[i]) << ped(e, i);
	  
	ped_out << "\n";
      }
      
      ped_out << "\n";
      //
    }// escape PEDs
  
    // hot energies-to-bimolecular PEDs
    //
    if(hot_index.size()) {
      //
      ped_out << "hot energies PEDs:\n";

      std::vector<int> bim_name_size(Model::bimolecular_size());
	  
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	itemp = Model::bimolecular(p).short_name().size() + 1;

	bim_name_size[p] = itemp > Model::ped_precision + 7 ? itemp : Model::ped_precision + 7;
      }
	    
      // dimensions
      //
      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(!b || outer_barrier(b).size() > itemp)
	  //
	  itemp = outer_barrier(b).size();
      
      Lapack::Matrix ped(itemp, Model::bimolecular_size());

      // PED
      //
      int count = 0;
	
      for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
	//
	ped = 0.;
	
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	  //
	  const int& w = Model::outer_connect(b).first;

	  const int& p = Model::outer_connect(b).second;
	  
	  for(int e = 0; e < outer_barrier(b).size(); ++e)
	    //
	    ped(e, p) += triple_product(&eigen_global(chem_size, e + well_shift[w]),
					//
					&eigen_hot(chem_size, count),
					//
					relax_lave, relax_size)

	      * global_bim(e + well_shift[w], p);
	}
	  
	//output
	//
	ped_out << "well: "<< Model::well(hit->first).short_name()
	  //
		<< "   hot energy = " << energy_bin(hit->second) / Phys_const::kcal
	  //
		<< " kcal/mol\n";
	
	ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  ped_out << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
	    
	ped_out << "\n";
	
	for(int e = 0; e < ped.size1(); ++e) {
	  //
	  ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(bim_name_size[p]) << ped(e, p);
	      
	  ped_out << "\n";
	}

	ped_out << "\n";
	//
      }// hot energy cycle
      //
    }// hot energies PEDs
    //
  }// PEDs output

  /*************************************** BOUND SPECIES *******************************************/

  std::vector<int>  group_index;
  Lapack::Matrix m_direct;
  std::vector<double>      weight;
  std::vector<double> real_weight;
  std::vector<double> real_ground;

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
      
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group);

      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }

    group_index = well_partition.group_index();
    weight      = well_partition.weight();
    real_weight = well_partition.real_weight(temperature());
    real_ground = well_partition.real_ground();

    // output
    //
    if(chem_size != Model::well_size()) {
      //
      IO::log << IO::log_offset << "combined species:\n"
	//
	      << IO::log_offset << std::setw(2) << "#"  << std::setw(Model::log_precision + 7)
	//
	      << "new name" << IO::first_offset
	//
	      << "group\n";
      
      for(int g = 0; g < well_partition.size(); ++g) {
	//
	IO::log << IO::log_offset << std::setw(2) << g << std::setw(Model::log_precision + 7)
	  //
		<< Model::well(group_index[g]).short_name() << IO::first_offset;
	
	for(Group::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w) {
	  //
	  if(w != well_partition[g].begin())
	    //
	    IO::log << "+";
	  
	  IO::log << Model::well(*w).short_name();
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
    
    //std::cout << funame << "well-to-well rate contributions:\n";

    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j) {
	//
	//std::cout << i << " <---> " << j << ":";
	
	for(int l = 0; l < chem_size; ++l) {
	  //
	  dtemp = m_direct(j, l) * m_inverse(l, i) * eigenval[l];

	  //std::cout << " " << dtemp;
	  
	  ww_rate(i, j) += dtemp;
	}
	//std::cout << "\n";
      }
    
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
    IO::log << IO::log_offset << "Wa->Wb/Wb->Wa rate constants ratios:\n"
      //
	    << IO::log_offset << std::setw(5) << "Wb\\Wa";
    
    for(int i = 0; i < chem_size; ++i)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::well(group_index[i]).short_name();
    
    IO::log << "\n";
    
    for(int j = 0; j < chem_size; ++j) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(group_index[j]).short_name();
      
      for(int i = 0; i < chem_size; ++i)
	//
	if(i != j) {
	  //
	  if(ww_rate(j, i) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << ww_rate(i, j) / ww_rate(j, i);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	}
	else
	  //
	  IO::log << std::setw(Model::log_precision + 7) << "1";
      
      IO::log << "\n";
    }
    
    if(Model::bimolecular_size()) {
      //
      IO::log << IO::log_offset << "W->P/P->W rate constants ratios:\n"
	//
	      << IO::log_offset << std::setw(5) << "P\\W";
      
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(Model::log_precision + 7) << Model::well(group_index[w]).short_name();
      
      IO::log << "\n";
    
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).short_name();
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  if(bw_rate(p, w) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << wb_rate(w, p) / bw_rate(p, w);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	
	IO::log << "\n";
      }
    }

    //IO::log << std::setprecision(6);

    // well-to-well rate coefficients
    //
    for(int i = 0; i < chem_size; ++i) {
      //
      dtemp = (real_ground[i] - energy_reference()) / temperature();

      double fac;
      
      if(dtemp > -Limits::exp_pow_max()) {
	//
	fac = real_weight[i] / std::exp(dtemp);
      }
      else
	//
	continue;
      
      for(int j = 0; j < chem_size; ++j) {
	//
	dtemp = ww_rate(i, j) * std::sqrt(weight[i] * weight[j]) * energy_step() / fac / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp; 	
      }
    }
    
    // well-to-escape rate coefficients
    //
    if(Model::escape_size())
      //
      for(int w = 0; w < chem_size; ++w) {
	//
	dtemp = (real_ground[w] - energy_reference()) / temperature();

	double fac;
      
	if(dtemp > -Limits::exp_pow_max()) {
	  //
	  fac = real_weight[w] / std::exp(dtemp);
	}
	else
	  //
	  continue;
      
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  dtemp = vdot(m_inverse.column(w), &eigen_escape(0, e))
	    //
	    * std::sqrt(weight[w]) * energy_step() / fac / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() +
				   //
				   Model::bimolecular_size() + e)] = dtemp;
	}
      }
    
    // output of steady state distributions associated with individual wells
    //
    IO::aux << "steady state distributions:\n";
    
    for(int g = 0; g < chem_size; ++g) {
      //
      IO::aux << "Initial well: " << Model::well(group_index[g]).short_name() << "\n";

      IO::aux  << std::setw(13) << "E, kcal/mol";

      for(int w = 0; w < Model::well_size(); ++w)
	//
	IO::aux << std::setw(13) << Model::well(w).short_name();

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
    for(int w = 0; w < chem_size; ++w) {
      //
      dtemp = (real_ground[w] - energy_reference()) / temperature();

      double fac;
      
      if(dtemp > -Limits::exp_pow_max()) {
	//
	fac = real_weight[w] / std::exp(dtemp);
      }
      else
	//
	continue;
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = wb_rate(w, p) * std::sqrt(weight[w]) * energy_step() / fac / Phys_const::herz;
	
	rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
      }
    }
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(!Model::bimolecular(p).dummy()) {
	//
	dtemp = (Model::bimolecular(p).ground() - energy_reference()) / temperature();

	double fac;
      
	if(dtemp > -Limits::exp_pow_max()) {
	  //
	  fac = Model::bimolecular(p).weight(temperature()) / std::exp(dtemp);
	}
	else
	  //
	  continue;
      
	for(int w = 0; w < chem_size; ++w) {
	  //
	  dtemp = bw_rate(p, w) * std::sqrt(weight[w]) * energy_step() / fac / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	}
      }
    
    // product energy distributions
    //
    if(ped_out.is_open()) {
      //
      std::vector<int> bim_name_size(Model::bimolecular_size());
	  
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	itemp = Model::ped_precision + 7;

	bim_name_size[p] = Model::bimolecular(p).short_name().size() + 1 > itemp ?  Model::bimolecular(p).short_name().size() + 1 : itemp;
      }
      
      // well-to-bimolecular distributions
      //
      if(Model::bimolecular_size()) {
	//
	ped_out << "well PEDs:\n";

	// dimensions
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(!b || outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

	Lapack::Matrix ped(itemp, Model::bimolecular_size());

	// PED
	//
	for(int c = 0; c < chem_size; ++c) {
	  //
	  ped = 0.;
	  
	  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	    //
	    const int& w = Model::outer_connect(b).first;
	      
	    const int& p = Model::outer_connect(b).second;
	      
	    for(int e = 0; e < outer_barrier(b).size(); ++e) 
	      //
	      ped(e, p) += vdot(m_inverse.column(c), &eigen_global(0, e + well_shift[w]))
		//
		* global_bim(e + well_shift[w], p);
	  }

	  // output
	  //
	  ped_out << "well: " << Model::well(group_index[c]).short_name() << "\n";
	  
	  ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
	  
	  ped_out << "\n";

	  for(int e = 0; e < ped.size1(); ++e) {
	    //
	    ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	    for(int p = 0; p < Model::bimolecular_size(); ++p)
	      //
	      ped_out << std::setw(bim_name_size[p]) << ped(e, p);
	    
	    ped_out << "\n";
	  }
	  
	  ped_out << "\n";
	  //
	}//well cycle
	//
      }// well PEDs   

      // well-to-escape PEDs
      //
      if(Model::escape_size()) {
	//
	ped_out << "well-to-escape PEDs:\n";

	// dimensions
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  const int& w = Model::escape_channel(s).first;
	  
	  if(!s || well(w).size() > itemp)
	    //
	    itemp = well(w).size();
	}
	
	Lapack::Matrix ped(itemp, chem_size * Model::escape_size());

	ped = 0.;
	
	// PED
	//
	for(int c = 0; c < chem_size; ++c)
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    const int& w = Model::escape_channel(s).first;
	      
	    const int& d = Model::escape_channel(s).second;
	      
	    for(int e = 0; e < well(w).size(); ++e) {
	      //
	      itemp = s + Model::escape_size() * c;
	  
	      ped(e, itemp) = vdot(m_inverse.column(c), &eigen_global(0, e + well_shift[w]))
		//
		* Model::well(w).escape_rate(energy_bin(e), d) * well(w).boltzman_sqrt(e);
	    }
	  }

	// output
	//
	ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	std::vector<int> ped_name_size(chem_size * Model::escape_size());
	
	for(int c = 0; c < chem_size; ++c)
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    stemp = " " + Model::well(group_index[c]).short_name() + "->" + Model::escape_name(s);

	    itemp = Model::ped_precision + 7;

	    itemp = itemp > stemp.size() ? itemp : stemp.size();

	    ped_name_size[s + Model::escape_size() * c] = itemp;
	    
	    ped_out << std::setw(itemp) << stemp;
	  }
	
	ped_out << "\n";
      
	for(int e = 0; e < ped.size1(); ++e) {
	  //
	  ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	  for(int i = 0; i < ped.size2(); ++i)
	    //
	    ped_out << std::setw(ped_name_size[i]) << ped(e, i);
	  
	  ped_out << "\n";
	}
	
	ped_out << "\n";
	//
      }// well-to-escape PEDs
      //
    }// product energy distributions
    //
  }// bound species

  // hot distribution branching ratios
  //
  if(hot_index.size()) {
    //
    int well_name_size_max;

    std::vector<int> well_name_size(Model::well_size());
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      itemp = Model::well(w).short_name().size() + 1;

      if(!w || itemp > well_name_size_max)
	//
	well_name_size_max = itemp;
    
      if(itemp < Model::log_precision + 7)
	//
	itemp = Model::log_precision + 7;

      well_name_size[w] = itemp;
    }

    if(well_name_size_max < 5)
      //
      well_name_size_max = 5;

    IO::log << IO::log_offset << "hot energies branching fractions:\n"
      //
	    << IO::log_offset //<< std::setprecision(6)
      //
	    << std::setw(well_name_size_max)  << "Well"
      //
	    << std::setw(Model::log_precision + 7) << "E, kcal";
    
    for(int w = 0; w < chem_size; ++w)
      //
      IO::log << std::setw(well_name_size[group_index[w]]) << Model::well(group_index[w]).short_name();

    std::vector<int> bim_name_size(Model::bimolecular_size());
    
    for(int p = 0; p < Model::bimolecular_size(); ++p) {
      //
      itemp = Model::bimolecular(p).short_name().size() + 1;

      if(itemp < Model::log_precision + 7)
	//
	itemp = Model::log_precision + 7;

      bim_name_size[p] = itemp;
      
      IO::log << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
    }

    for(int e = 0; e < Model::escape_size(); ++e)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::escape_name(e);
    
    IO::log << "\n";
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
      //
      IO::log << IO::log_offset
	//
	      << std::setw(well_name_size_max)  << Model::well(hit->first).short_name()
	//
	      << std::setw(Model::log_precision + 7) << energy_bin(hit->second) / Phys_const::kcal;
	
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(well_name_size[group_index[w]]) << vdot(m_direct.row(w), &eigen_hot(0, count)) * std::sqrt(weight[w]);
	
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	IO::log << std::setw(bim_name_size[p]) << triple_product(&eigen_bim(chem_size, p), &eigen_hot(chem_size, count),
								 //
								 relax_lave, relax_size);
	
      for(int e = 0; e < Model::escape_size(); ++e)
	//
	IO::log << std::setw(Model::log_precision + 7) << triple_product(&eigen_escape(chem_size, e), &eigen_hot(chem_size, count),
									 //
									 relax_lave, relax_size);
	
      IO::log << "\n";
    }
  }

  // prompt isomerization
  //
  IO::log << IO::log_offset << "time = " << std::time(0) - start_time << "\n";

  IO::log << IO::log_offset << "prompt isomerization/dissociation:\n";

  std::vector<int> well_name_size(Model::well_size());

  int well_name_size_max;

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    itemp = Model::well(w).short_name().size() + 1;

    if(!w || itemp > well_name_size_max)
      //
      well_name_size_max = itemp;
    
    if(itemp < Model::log_precision + 7)
      //
      itemp = Model::log_precision + 7;

    well_name_size[w] = itemp;
  }

  if(well_name_size_max < 7)
    //
    well_name_size_max = 7;

  IO::log << IO::log_offset << std::setw(well_name_size_max) << "W\\W,P,E";

  if(chem_size)
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(well_name_size[w]) << Model::well(w).short_name();

  IO::log << std::setw(Model::log_precision + 7) << "Total";

  std::vector<int> bim_name_size(Model::bimolecular_size());
  
  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    itemp = Model::bimolecular(p).short_name().size() + 1;

    if(itemp < Model::log_precision + 7)
      //
      itemp = Model::log_precision + 7;

    bim_name_size[p] = itemp;
    
    IO::log << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
  }
  
  for(int e = 0; e < Model::escape_size(); ++e)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::escape_name(e);
    
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    dtemp = (Model::well(w).ground() - energy_reference()) / temperature();

    double fac;
      
    if(dtemp > -Limits::exp_pow_max()) {
      //
      fac = Model::well(w).weight(temperature()) / std::exp(dtemp);
    }
    else
      //
      continue;
      
    IO::log << IO::log_offset << std::setw(well_name_size_max) << Model::well(w).short_name();

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
      
	dtemp *= well(v).weight_sqrt() * well(w).weight_sqrt() * energy_step() / fac;

	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	diss -= dtemp;

	IO::log << std::setw(well_name_size[v]) << dtemp;
      }

    IO::log << std::setw(Model::log_precision + 7) << diss;
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(bim_name_size[p])
	//
	      << triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w), relax_lave, relax_size)
	//
	* energy_step()	* well(w).weight_sqrt() / fac;

    for(int e = 0; e < Model::escape_size(); ++e)
      //
      IO::log << std::setw(Model::log_precision + 7)
	//
	      << triple_product(&eigen_escape(chem_size, e), &eigen_pop(chem_size, w), relax_lave, relax_size)
	//
	* energy_step()	* well(w).weight_sqrt() / fac;
      
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

double MasterEquation::threshold_well_partition (Lapack::Matrix pop_chem, Partition& well_partition, Group& bimolecular_group) 
{
  const char funame [] = "MasterEquation::threshold_well_partition: ";

  IO::Marker funame_marker(funame);
  
  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  const int well_size = pop_chem.size1();

  const int chem_size = pop_chem.size2();
  
  if(Model::well_size() != well_size) {
    //
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1()
      //
	      << " differs from the number of wells = " << well_size << "\n";
    
    throw Error::Logic();
  }

  if(chem_size >= well_size) {
    //
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";
    
    throw Error::Logic();
  }

  IO::log << IO::log_offset << "well projection threshold  = " << well_projection_threshold << "\n";
  
  // projection-ordered wells
  //
  std::multimap<double, int> proj_well_map;
  
  bimolecular_group.clear();
  
  for(int w = 0; w < well_size; ++w) {
    //
    proj_well_map.insert(std::make_pair(Group(w).projection(pop_chem), w));
    
    bimolecular_group.insert(w);
  }

  // large projection wells map
  //
  std::vector<int> prim_well_map;

  std::multimap<double, int>::const_reverse_iterator pit;

  IO::log << IO::log_offset << "primary well projection/well:";
  
  for(pit = proj_well_map.rbegin(); pit != proj_well_map.rend(); ++pit) {
    //
    if(pit->first < well_projection_threshold && prim_well_map.size() >= chem_size)
      //
      break;

    IO::log << "   " << pit->first << "/" << pit->second;
    
    dtemp = pit->first;

    prim_well_map.push_back(pit->second);
    
    bimolecular_group.erase(pit->second);
  }

  IO::log << "\n";

  if(pit != proj_well_map.rend()) {
    //
    IO::log << IO::log_offset << "secondary well projection/well:";

    for(; pit != proj_well_map.rend(); ++pit)
      //
      IO::log << "   " << pit->first << "/" << pit->second;
  }

  IO::log << "\n";
    
  IO::log << IO::log_offset << "# of primary large projection wells / kinetically active species  = "
    
	  << prim_well_map.size() << " / " << chem_size << "\n";
  
  
  
  if(dtemp < well_projection_threshold)
    //
    IO::log << IO::log_offset << funame
      //
	    << "WARNING: primary well projection is smaller than the well projection threshold: "
      //
	    << dtemp << "\n";
  
  double pmax = -1.;
  
  for(PartitionGenerator pg(chem_size, prim_well_map.size()); !pg.end(); ++pg) {
    //
    // initialize new partition
    //
    Partition p(pg, prim_well_map);

    // partition projection
    //
    dtemp = p.projection(pop_chem);

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
    
    for(Group::const_iterator w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
      //
      for(int s = 0; s < chem_size; ++s) {
	//
	Group g = well_partition[s];

	g.insert(*w);
	
	dtemp = g.projection(pop_chem) - well_partition[s].projection(pop_chem);
	
	if(dtemp > pmax) {
	  //
	  pmax = dtemp;
	  
	  smax = s;
	  
	  wmax = *w;
	}
      }
    }
    
    if(pmax <= 0.)
      //
      break;

    well_partition[smax].insert(wmax);
    
    bimolecular_group.erase(wmax);
  }

  pmax = (double)chem_size - well_partition.projection(pop_chem);

  IO::log << IO::log_offset << "partition projection error = " << pmax << "\n";

  return pmax;
}

double MasterEquation::sort_well_partition (Lapack::Matrix pop_chem, Partition& well_partition, Group& bimolecular_group) 
{
  const char funame [] = "MasterEquation::sort_well_partition: ";
  
  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  const int well_size = pop_chem.size1();

  if(Model::well_size() != well_size) {
    //
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1() 
	      << " differs from the number of wells = " << well_size << "\n";
    
    throw Error::Logic();
  }

  if(chem_size >= well_size) {
    //
    std::cerr << funame << "number of kinetically active species should be less than the number of wells\n";

    throw Error::Logic();
  }

  std::multimap<double, std::pair<Partition, Group> > high_part;

  // partitioning the wells into equilibrated groups
  //
  int part_size = chem_size;
  
  for(PartitionGenerator pg(part_size, well_size); !pg.end(); ++pg) {
    //
    // initialize new partition
    //
    Partition part(pg);

    // partition projection
    //
    dtemp = part.projection(pop_chem);

    if(high_part.size() < red_out_num)
      //
      high_part.insert(std::make_pair(dtemp, std::make_pair(part, Group())));
    
    else if(dtemp > high_part.begin()->first) {
      //
      high_part.erase(high_part.begin());
      
      high_part.insert(std::make_pair(dtemp, std::make_pair(part, Group())));
    }
  }

  // the case when some of the wells equilibrate with bimolecular products
  // is solved by allowing the bimolecular group and then removing it
  //
  part_size = chem_size + 1;
  
  for(PartitionGenerator pg(part_size, well_size); !pg.end(); ++pg) {
    //
    // initialize new partition
    //
    Partition part(pg);
    
    for(int g = 0; g < part_size; ++g) {
      Partition q = part;
      Group bg = q[g];
      // removing the bimolecular group
      q.erase(q.begin() + g);
      // partition projection
      dtemp = q.projection(pop_chem);
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
  //
  IO::log << IO::log_offset << "single well projections onto chemical subspace:\n";
  
  IO::log << IO::log_offset << std::left << std::setw(6) << "well:" << std::right;
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << std::setw(6) << Model::well(w).short_name();
  
  IO::log << "\n";
  IO::log << IO::log_offset << std::left << std::setw(6) << "proj:" << std::right;
  
  for(int w = 0; w < Model::well_size(); ++w) {
    //
    Group g;
    
    g.insert(w);
    
    IO::log << std::setw(6) << g.projection(pop_chem);
  }
  
  IO::log << "\n";
  IO::log.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);

  IO::log << IO::log_offset << "closest well partitions:\n";
  std::multimap<double, std::pair<Partition, Group> >::const_reverse_iterator rit;
  IO::log << IO::log_offset 
	  << std::setw(Model::log_precision + 7) << "reduction error" 
	  << std::setw(50) << "equilibrated_groups"
	  << std::setw(50) << "wells_equilibrated_with_products"
	  << "\n";

  //IO::log << std::setprecision(6);
    
  for(rit = high_part.rbegin(); rit != high_part.rend(); ++rit) {
    //
    IO::log << IO::log_offset << std::setw(Model::log_precision + 7) << (double)chem_size - rit->first;
    
    // equilibrated bound species
    //
    stemp.clear();
    for(int g = 0; g < rit->second.first.size(); ++g) {
      //
      if(!g)
	//
	stemp += " ";
	
      for(Group::const_iterator w = rit->second.first[g].begin(); w != rit->second.first[g].end(); ++w) {
	//
	if(w != rit->second.first[g].begin())
	  //
	  stemp += "+";
	  
	stemp += Model::well(*w).short_name();
      }
    }

    if(stemp.size())
      //
      IO::log << std::setw(50) << stemp;
    else
      //
      IO::log << std::setw(50) << "---";

    // wells equilibrated with the products
    //
    if(rit->second.second.size()) {
      //
      stemp.clear();
	
      for(Group::const_iterator w = rit->second.second.begin(); w != rit->second.second.end(); ++w) {
	//
	if(w != rit->second.second.begin())
	  //
	  stemp += " ";
	  
	stemp += Model::well(*w).short_name();
      }

      IO::log << std::setw(50) << stemp;    
    }
    else
      //
      IO::log << std::setw(50) << "---";
      
    IO::log << "\n";
  }

  return (double)chem_size - high_part.rbegin()->first;
}

double MasterEquation::incremental_well_partition (Lapack::Matrix pop_chem, Partition& well_partition, Group& bimolecular_group) 
{
  const char funame [] = "MasterEquation::increment_well_partition: ";
  
  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  const int well_size = pop_chem.size1();

  if(Model::well_size() != well_size) {
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
    proj_well_map.insert(std::make_pair(Group(w).projection(pop_chem), w));
    bimolecular_group.insert(w);
  }

  well_partition.resize(chem_size);

  itemp = 0;
  for(std::multimap<double, int>::const_reverse_iterator i = proj_well_map.rbegin(); itemp < chem_size; ++i, ++itemp) {
    well_partition[itemp].insert(i->second);
    bimolecular_group.erase(i->second);
  }
    
  while(bimolecular_group.size()) {
    //
    int smax, wmax;
    
    double proj_max = -1.;
    
    for(int s = 0; s < chem_size; ++s) {
      //
      double sp = well_partition[s].projection(pop_chem);
      
      for(Group::const_iterator w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
	//
	Group g(well_partition[s]);

	g.insert(*w);
	
	dtemp = g.projection(pop_chem) - sp;
		      
	if(dtemp > proj_max) {
	  //
	  proj_max = dtemp;
	  
	  smax = s;
	  
	  wmax = *w;
	}
      }
    }
      
    if(proj_max < 0.)
      //
      break;

    well_partition[smax].insert(wmax);
    bimolecular_group.erase(wmax);
  }


  return (double)chem_size - well_partition.projection(pop_chem);
}

double MasterEquation::sequential_well_partition (Lapack::Matrix pop_chem, Partition& well_partition, Group& bimolecular_group) 
{
  const char funame [] = "MasterEquation::sequential_well_partition: ";
  
  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;
  std::string stemp;

  const int chem_size = pop_chem.size2();
  
  const int well_size = pop_chem.size1();

  if(Model::well_size() != well_size) {
    //
    std::cerr << funame << "first dimension of the pop_chem matrix = " << pop_chem.size1()
      //
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
    //
    double spec_proj = 0.;
    
    while(bimolecular_group.size()) {
      //
      int wmax;
      
      double proj_diff = -1.;
      
      for(Group::const_iterator w = bimolecular_group.begin(); w !=bimolecular_group.end(); ++w) {
	//
	Group g(well_partition[spec]);

	g.insert(*w);
	
	dtemp = g.projection(pop_chem) - spec_proj;
	
	if(dtemp > proj_diff) {
	  //
	  proj_diff = dtemp;
	  
	  wmax = *w;
	}
      }
      
      if(proj_diff < 0.)
	//
	break;

      well_partition[spec].insert(wmax);
      
      bimolecular_group.erase(wmax);
      
      spec_proj += proj_diff;
    }

    if(!well_partition[spec].size()) {
      //
      std::cerr << funame << "empty species\n";
      
      throw Error::Logic();
    }
  }
  
  return (double)chem_size - well_partition.projection(pop_chem);
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

std::vector<std::set<int> > MasterEquation::PartitionGenerator::partition (const std::vector<int>& map) const
{
  const char funame[] = "MasterEquation::PartitionGenerator::partition: ";

  std::vector<std::set<int> > res(partition_size());
  
  if(!map.size()) {
    //
    for(int i = 0; i < size(); ++i)
      //
      res[(*this)[i]].insert(i);
  }
  else if(map.size() == size()) {
    //
    for(int i = 0; i < size(); ++i)
      //
      res[(*this)[i]].insert(map[i]);
  }
  else {
    //
    std::cerr << funame << "map dimension mismatch\n";
    
    throw Error::Range();
  }

  return res;
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

  if(!map.size()) {
    //
    for(int i = 0; i < g.size(); ++i)
      //
      (*this)[g[i]].insert(i);
  }
  else if(map.size() == g.size()) {
    //
    for(int i = 0; i < g.size(); ++i)
      //
      (*this)[g[i]].insert(map[i]);
  }
  else {
    //
    std::cerr << funame << "map dimension mismatch\n";
    
    throw Error::Range();
  }
}

std::vector<double> MasterEquation::Partition::weight () const
{
  std::vector<double> res(size());

  for(const_iterator g = begin(); g != end(); ++g)
    //
    res[g - begin()] = g->weight();
  
  return res;
}

std::vector<double> MasterEquation::Partition::real_weight (double t) const
{
  std::vector<double> res(size());
  
  for(const_iterator g = begin(); g != end(); ++g)
    //
    res[g - begin()] = g->real_weight(t);
  
  return res;
}

std::vector<double> MasterEquation::Partition::real_ground () const
{
  std::vector<double> res(size());
  
  for(const_iterator g = begin(); g != end(); ++g)
    //
    res[g - begin()] = g->real_ground();
  
  return res;
}

std::vector<int> MasterEquation::Partition::group_index () const 
{
  std::vector<int>  res(size());
  
  for(const_iterator g = begin(); g != end(); ++g)
    //
    res[g - begin()] = g->group_index();
  
  return res;
}

double MasterEquation::Partition::projection (Lapack::Matrix pop_chem) const 
{
  double res = 0.;
  
  for(const_iterator g = begin(); g != end(); ++g)
    
    res += g->projection(pop_chem);

  return res;
}
	
Lapack::Matrix MasterEquation::Partition::basis () const 
{
  Lapack::Matrix res(Model::well_size(), size());
  
  for(const_iterator g = begin(); g != end(); ++g)
    //
    res.column(g - begin()) = g->basis_vector();
  
  return res;
}

/**************************************** Well Group *******************************************/

void MasterEquation::Group::_assert () const
{
  const char funame [] = " MasterEquation::Group::_assert: ";

  if(!size()) {
    //
    IO::log << IO::log_offset << funame << "group is empty" << std::endl;

    throw Error::Init();
  }

  if(*begin() < 0 || *rbegin() >= Model::well_size()) {
    //
    IO::log << IO::log_offset << funame << "indices out of range: " << *begin() << ", " << *rbegin() << std::endl;

    throw Error::Range();
  }
}

double MasterEquation::Group::weight () const
{
  _assert();
  
  double res = 0.;
  
  for(const_iterator w = begin(); w != end(); ++w)
    //
    res += well(*w).weight();
  
  return res;
}

double MasterEquation::Group::real_ground () const
{
  _assert();
  
  double res;
  
  for(const_iterator w = begin(); w != end(); ++w)
    //
    if(w == begin() || res > Model::well(*w).ground())
      //
      res = Model::well(*w).ground();
  
  return res;
}

double MasterEquation::Group::real_weight (double t) const
{
  const char funame [] = " MasterEquation::Group::real_weight: ";

  _assert();
  
  if(t <= 0.) {
    //
    IO::log << IO::log_offset << funame << "temperature out of range: " << t / Phys_const::kelv << std::endl;

    throw Error::Range();
  }
  
  double dtemp;

  double res = 0.;

  double g = real_ground();
  
  for(const_iterator w = begin(); w != end(); ++w) {
    //
    dtemp = (Model::well(*w).ground() - g) / t;

    if(dtemp < Limits::exp_pow_max())
      //
      res += Model::well(*w).weight(t) / std::exp(dtemp);
  }
  
  return res;
}

int MasterEquation::Group::group_index () const
{
  _assert();
  
  double dtemp, gmin;

  int res;
  
  for(const_iterator w = begin(); w != end(); ++w) {
    //
    dtemp =  Model::well(*w).ground();
    
    if(w == begin() || dtemp < gmin) {
      //
      gmin = dtemp;
      
      res = *w;
    }
  }
  
  return res;
}

Lapack::Vector MasterEquation::Group::basis_vector () const
{
  _assert();
  
  Lapack::Vector res(Model::well_size());
  
  res = 0.;

  if(size() == 1) {
    //
    res[*begin()] = 1.;
    
    return res;
  }
    
  double dtemp = std::sqrt(weight());

  for(const_iterator w = begin(); w != end(); ++w)
    //
    res[*w] = well(*w).weight_sqrt() / dtemp;
  
  return res;
}

bool MasterEquation::Group::insert (const Group& group)
{
  bool btemp;
  
  bool res = true;
  
  for(const_iterator w = group.begin(); w != group.end(); ++w) {
    //
    btemp = insert(*w);

    res = btemp && res;
  }
  
  return res;
}

bool MasterEquation::Group::erase (const Group& group)
{
  bool btemp;
  
  bool res = true;
  
  for(const_iterator g = group.begin(); g != group.end(); ++g) {
    //
    btemp = erase(*g);

    res = btemp && res;
  }
  
  return res;
}

bool MasterEquation::Group::contain (const Group& group) const
{
  for(const_iterator w = group.begin(); w != group.end(); ++w)
    //
    if(find(*w) == end())
      //
      return false;
  
  return true;
}

double MasterEquation::Group::projection (Lapack::Matrix pop_chem) const
{
  const char funame [] = "MasterEquation::Group::projection: ";

  _assert();
  
  double dtemp;

  const int well_size = pop_chem.size1();

  const int chem_size = pop_chem.size2();
  
  if(Model::well_size() != well_size) {
    //
    std::cerr << funame << "pop_chem: first dimension out of range: " << pop_chem.size1() << "\n";
    
    throw Error::Logic();
  }

  if(size() == 1)
    //
    return vdot(pop_chem.row(*begin()));

  double res = 0.;
  
  for(int l = 0; l < chem_size; ++l) {
    //
    dtemp = 0.;
    
    for(const_iterator w = begin(); w != end(); ++w)
      //
      dtemp += well(*w).weight_sqrt() * pop_chem(*w, l);
    
    res += dtemp * dtemp;
  }

  // normalization
  //
  res /= weight();

  return res;
}
