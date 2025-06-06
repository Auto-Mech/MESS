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

#include "new_mess.hh"
#include "units.hh"
#include "key.hh"
#include "io.hh"
#include "shared.hh"
#include "offload.hh"
#include "mpack.hh"

namespace MasterEquation {

  typedef dd_real float_t;
  
  double  temperature = -1.;

  double  pressure    = -1.;

  /********************************* USER DEFINED PARAMETERS ********************************/

  // use multi-precision for global matrix diagonalization
  //
  int                                                      use_mp = 0;
  
  // energy step
  //
  double                                                    energy_step_over_temperature = -1.;

  // highest energy
  //
  double                                                    excess_energy_over_temperature = -1.;

  // thermal distribution width multiplier
  //
  double                                                    dist_width_factor     = 5.;
  
  // well cutoff parameter
  //
  double                                                     well_cutoff           = 10.;
  
  // highest chemecial eigenvalue
  //
  double                                                     chemical_threshold     = -2.;
  
  // smallest chemical eigenvalue
  //
  double                                                     chemical_tolerance     = -1.;
  
  // fastest isomerization eigenvalue over collision frequency
  //
  double                                                     reduction_threshold    = 10.;

  // slowest isomerization eigenvalue over collision frequency
  //
  double                                                     reduction_tolerance    = .1;

  // minimal gap between thermal and active isomerization eigenvalues
  //
  double                                                     reduction_increment    = 100.;
  
  // microcanonical rate maximum
  //
  double                                                     micro_rate_max          = -1.;
  
  // pressure units
  //
  int                                                        pressure_unit = BAR;

  // well partition thresholds
  //
  double                                                     well_projection_threshold = 0.2;

  double                                                     well_projection_tolerance = -1.;
  
  // well extention parameter
  //
  double                                                     well_extension = -1.;

  // well extension correction
  //
  double                                                     well_ext_corr = 0.;

  // time propagation parameters
  //
  double                                                     time_step  = 0.1; // time step * fastest relaxation eigenvalue = coll_freq * reduction_threshold

  double                                                     time_limit = 10.; // time limit * collisional frequency

  double                                                     time_output = 1.; // output time interval * collisional frequency

  // bimolecular-to-bimolecular rates calculation method
  //
  int                                                        use_matrix_inversion = 1;

  // rate calculation management
  //
  double                                                     ener_gap_min = 10.;
 
  double                                                     ener_range_max = 20.;

  int                                                        use_thermal_map = 0;
  //
}// MasterEquation

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
    return -1;
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
 *********************************** REACTIVE COMPLEX ***************************************
 ********************************************************************************************/

void MasterEquation::ReactiveComplex::PesObject::assert(int e) const
{
  const char funame [] = " MasterEquation::ReactiveComplex::PesObject::assert: ";
  
  if(e >= size()) {
    //
    IO::log << IO::log_offset << funame << "index out of range\n";

    throw Error::Range();
  }
}

void MasterEquation::ReactiveComplex::PesObject::set ()
{
  const char funame [] = "MasterEquation::ReactiveComplex::PesObject::set: ";

  double dtemp;
  
  if(_weight > 0.) {
    //
    IO::log << IO::log_offset << funame << "has been set already\n";
    
    throw Error::Init();
  }

  // statistical weight
  //
  _weight = 0.;

  for(int e = 0; e < size(); ++e) {
    //
    dtemp = _states[e];

    if(dtemp <= 0.) {
      //
      IO::log << IO::log_offset << funame << "non-positive states: " << dtemp << "\n";

      throw Error::Range();
    }
    
    _weight += dtemp * thermal_factor(e);
  }
}

double MasterEquation::ReactiveComplex::PesObject::real_weight () const
{
  const char funame [] = "MasterEquation::ReactiveComplex::PesObject::real_weight: ";
      
  switch(_type) {
    //
  case Model::WELL:
    //
    return Model::well(_index).weight(temperature) * std::exp(-Model::well(_index).ground() / temperature);

  case Model::INNER:
    //
    return Model::inner_barrier(_index).weight(temperature) * std::exp(-Model::inner_barrier(_index).ground() / temperature);

  case Model::OUTER:
    //
    return Model::outer_barrier(_index).weight(temperature) * std::exp(-Model::outer_barrier(_index).ground() / temperature);
  }

  IO::log << IO::log_offset << funame << "should not be here\n";

  throw Error::Logic();
}
    
std::string MasterEquation::ReactiveComplex::PesObject::name () const
{
  const char funame [] = "MasterEquation::ReactiveComplex::PesObject::name: ";
      
  switch(_type) {
    //
  case Model::WELL:
    //
    return Model::well(_index).short_name();

  case Model::INNER:
    //
    return Model::inner_barrier(_index).short_name();

  case Model::OUTER:
    //
    return Model::outer_barrier(_index).short_name();
  }

  IO::log << IO::log_offset << funame << "should not be here\n";

  throw Error::Logic();
}
    
double  MasterEquation::ReactiveComplex::PesObject::weight () const
{
  const char funame [] = "MasterEquation::ReactiveComplex::PesObject::weight: ";
	
  if(_weight < 0.) {
    //
    IO::log << IO::log_offset << funame << name() << ": not initialized: " << _weight << "\n";

    throw Error::Init();
  }
	
  return _weight;
}

void MasterEquation::ReactiveComplex::Well::init (int w, const ReactiveComplex& rc)
{
  const char funame [] = "MasterEquation::ReactiveComplex::Well::init: ";

  double dtemp;
  int    itemp;
  
  PesObject::init(Model::WELL, w, rc);

  _kernel_fraction.resize(Model::buffer_size());
  
  _collision_frequency = 0.;

  for(int b = 0; b < Model::buffer_size(); ++b) {
    //
    dtemp = Model::well(w).collision(b)(temperature) * Model::buffer_fraction(b);
    
    _kernel_fraction[b] = dtemp;
    
    _collision_frequency += dtemp;
  }

  for(int b = 0; b < Model::buffer_size(); ++b)
    //
    _kernel_fraction[b] /= _collision_frequency;

  _collision_frequency *= pressure;
}

void MasterEquation::ReactiveComplex::Well::set ()
{
  const char funame [] = "MasterEquation::ReactiveComplex::Well::set: ";
  
  double  dtemp;
  int     itemp;
  bool    btemp;

  if(iset()) {
    //
    IO::log << IO::log_offset << funame << "has been set already\n";
    
    throw Error::Init();
  }
  
  PesObject::set();
  
  /***********************************************************************************
   **************************** ENERGY TRANSFER KERNEL *******************************
   ***********************************************************************************/
  
  _kernel.resize(size());
    
  _kernel = 0.;

  Lapack::Matrix tmp_kernel(size());

  double nfac;

  for(int b = 0; b < Model::buffer_size(); ++b) {
    //
    tmp_kernel = 0.;

    // collision energy transfer on the grid
    //
    itemp = (int)std::ceil(Model::well(index()).kernel(b).cutoff_energy(temperature) / energy_step());

    if(itemp < 2) {
      //
      IO::log << IO::log_offset << "no collision energy transfer\n";

      throw Error::Range();
    }
    
    std::vector<double> energy_transfer(itemp);
	
    for(int e = 0; e < energy_transfer.size(); ++e)
      //
      energy_transfer[e] = Model::well(index()).kernel(b)((double)e * energy_step(), temperature);

    if(!b || kernel_bandwidth < energy_transfer.size())
      //
      kernel_bandwidth = energy_transfer.size();

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

	tmp_kernel(itemp, e) = energy_transfer[i] * states(e) / states(itemp) / thermal_factor(i);

	nfac += tmp_kernel(e, itemp) + tmp_kernel(itemp, e);
      }

      nfac /= _kernel_fraction[b];
	  
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
      
    _kernel += tmp_kernel;
  }
}

double MasterEquation::ReactiveComplex::Barrier::states (int e) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::Barrier::states: ";
  
  double dtemp;

  const double res = PesObject::states(e);
  
  if(micro_rate_max < 0.)
    //
    return res;
  
  switch(type()) {
    //
  case Model::INNER:
    //
    dtemp = std::min(reactive_complex().well[connect.first ].states(e), reactive_complex().well[connect.second].states(e)) * micro_rate_max * 2. * M_PI;

    break;
    //
  case Model::OUTER:
    //
    dtemp = reactive_complex().well[connect.first].states(e) * micro_rate_max * 2. * M_PI;

    break;
    //
  default:
    //
    IO::log << IO::log_offset << "wrong type: " << type() << "\n";

    throw Error::Range();
  }

  return std::min(dtemp, res);
}

// reactive complex initialization
//
MasterEquation::ReactiveComplex::ReactiveComplex (const Model::ChemGraph& cg, int er, int flags)
  //
  : Model::ChemGraph(cg), _ener_index_begin(er), full_energy_range(0)
{
  const char funame [] = "MasterEquation::ReactiveComplex::ReactiveComplex: ";
  
  IO::Marker marker("ReactiveComplex");

  double  dtemp;
  int     itemp;
  bool    btemp;

  cg.assert();

  IO::log << IO::log_offset << "energy reference =      " << std::floor(energy_bin(0) / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol\n";

  IO::log << IO::log_offset << "well/ground[kcal/mol]:";

  for(std::set<int>::const_iterator w = well_set.begin(); w != well_set.end(); ++w)
    //
    IO::log << "  " << Model::well(*w).short_name() << "/" << std::floor(Model::well(*w).ground() / Phys_const::kcal * 10. + 0.5) / 10.;
      
  IO::log << "\n";

  if(inner_set.size()) {
    //
    IO::log << IO::log_offset << "inner/hight[kcal/mol]:";

    for(std::set<int>::const_iterator b = inner_set.begin(); b != inner_set.end(); ++b)
      //
      IO::log << "  " << Model::inner_barrier(*b).short_name() << "/" << std::floor(Model::inner_barrier(*b).real_ground() / Phys_const::kcal * 10. + 0.5) / 10.;
      
    IO::log << "\n";
  }

  if(outer_set.size()) {
    //
    IO::log << IO::log_offset << "outer/hight[kcal/mol]:";

    for(std::set<int>::const_iterator b = outer_set.begin(); b != outer_set.end(); ++b)
      //
      IO::log << "  " << Model::outer_barrier(*b).short_name() << "/" << std::floor(Model::outer_barrier(*b).real_ground() / Phys_const::kcal * 10. + 0.5) / 10.;
      
    IO::log << "\n";
  }
  else
    //
    IO::log << IO::log_offset << "isomerization only\n";

  IO::log << "\n";
  
  if(!is_connected()) {
    //
    IO::log << IO::log_offset << "graph is not connected\n";

    throw Error::Logic();
  }

  // well initialization
  //
  int count = 0;

  well.resize(well_set.size());

  _index_well_map.resize(well_set.size());
  
  for(std::set<int>::const_iterator w = well_set.begin(); w != well_set.end(); ++w, ++count) {
    //
    _well_index_map[*w] = count;

    _index_well_map[count] = *w;

    well[count].init(*w, *this);
  }

  // inner barrier initialization
  //
  count = 0;
  
  inner_barrier.resize(inner_set.size());
  
  for(std::set<int>::const_iterator b = inner_set.begin(); b != inner_set.end(); ++b, ++count) {
    //
    inner_barrier[count].init(Model::INNER, *b, *this);

    const int w1 = Model::inner_connect(*b).first;
    
    const int w2 = Model::inner_connect(*b).second;

    inner_barrier[count].connect.first  = well_to_index(w1);
    
    inner_barrier[count].connect.second = well_to_index(w2);
  }
  
  // outer barrier initialization
  //
  count = 0;
  
  outer_barrier.resize(outer_set.size());
  
  for(std::set<int>::const_iterator b = outer_set.begin(); b != outer_set.end(); ++b, ++count) {
    //
    outer_barrier[count].init(Model::OUTER, *b, *this);

    const int w = Model::outer_connect(*b).first;
    
    const int p = Model::outer_connect(*b).second;

    if(_bim_index_map.find(p) == _bim_index_map.end()) {
      //
      _bim_index_map[p] = _index_bim_map.size();

      _index_bim_map.push_back(p);
    }

    outer_barrier[count].connect.first = well_to_index(w);

    outer_barrier[count].connect.second = bim_to_index(p);
  }

  set_kinetic_matrix();
}

int MasterEquation::ReactiveComplex::KineticBasis::well_to_index (int w) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::KineticBasis::well_to_index: ";
  
  std::map<int, int>::const_iterator wit = well_index_map.find(w);

  if(wit != well_index_map.end())
    //
    return wit->second;

  IO::log << IO::log_offset << funame << "well index not in the map: " << w << "\n";
  
  throw Error::Logic();
}

std::set<int> MasterEquation::ReactiveComplex::KineticBasis::well_to_index (const std::set<int>& g) const
{
  std::set<int> res;
  
  for(std::set<int>::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res.insert(well_to_index(*w));

  return res;
}
  
int MasterEquation::ReactiveComplex::KineticBasis::index_to_well (int i) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::KineticBasis::index_to_well: ";
  
  if(i >= 0 && i < index_well_map.size())
    //
    return index_well_map[i];

  IO::log << IO::log_offset << funame << "index out of range: " << i << "\n";

  throw Error::Logic();
}

std::set<int> MasterEquation::ReactiveComplex::KineticBasis::index_to_well (const std::set<int>& g) const
{
  std::set<int> res;
  
  for(std::set<int>::const_iterator w = g.begin(); w != g.end(); ++w)
    //
    res.insert(index_to_well(*w));

  return res;
}
  
int MasterEquation::ReactiveComplex::well_to_index (int w) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::well_to_index: ";
  
  std::map<int, int>::const_iterator wit = _well_index_map.find(w);

  if(wit != _well_index_map.end())
    //
    return wit->second;

  IO::log << IO::log_offset << funame << "well index not in the map: " << w << "\n";
  
  throw Error::Logic();
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
  
int MasterEquation::ReactiveComplex::index_to_well (int i) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::index_to_well: ";
  
  if(i >= 0 && i < _index_well_map.size())
    //
    return _index_well_map[i];

  IO::log << IO::log_offset << funame << "index out of range: " << i << "\n";

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
  
std::vector<std::set<int> > MasterEquation::ReactiveComplex::index_to_well (const std::vector<std::set<int> >& p) const
{
  std::vector<std::set<int> > res(p.size());

  for(int g = 0; g < p.size(); ++g)
    //
    res[g] = index_to_well(p[g]);
  
  return res;
}
  
int MasterEquation::ReactiveComplex::bim_to_index (int p) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::bim_to_index: ";
  
  std::map<int, int>::const_iterator pit = _bim_index_map.find(p);

  if(pit != _bim_index_map.end())
    //
    return pit->second;

  IO::log << IO::log_offset << funame << "bimolecular index not in the map: " << p << "\n";

  throw Error::Logic();
}

int MasterEquation::ReactiveComplex::index_to_bim (int i) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::index_to_bim: ";
  
  if(i >= 0 && i < _index_bim_map.size())
    //
    return _index_bim_map[i];

  IO::log << IO::log_offset << funame << "index out of range: " << i << "\n";

  throw Error::Logic();
}

std::string MasterEquation::ReactiveComplex::name (const group_t& g) const
{
  _assert(g);

  std::string res;
      
  for(group_t::const_iterator w = g.begin(); w != g.end(); ++w) {
    //
    if(w != g.begin())
      //
      res += "+";

    res += well[*w].name();
  }
  
  return res;
}
      
std::string MasterEquation::ReactiveComplex::name (const partition_t& p) const
{
  std::string res;
      
  for(int g = 0; g < p.size(); ++g) {
    //
    if(g != 0)
      //
      res += "  ";

    res += name(p[g]);
  }
  
  return res;
}
      
void MasterEquation::ReactiveComplex::_assert (const group_t& g) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::_assert:";
  
  if(!g.size()) {
    //
    IO::log << IO::log_offset << "group is empty\n";

    throw Error::Init();
  }
	
  if(*g.begin() < 0 || *g.rbegin() >= well_size()) {
    //
    IO::log << IO::log_offset << funame << "well indices out of range\n";
	
    throw Error::Range();
  }
}

/***************************************************************************************
 ********************************** RATE PARTITIONING **********************************
 ***************************************************************************************/

// populations in different wells
//
Lapack::Vector MasterEquation::ReactiveComplex::population (double* state_vector) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::population: ";

  double dtemp;
	
  Lapack::Vector res(well.size(), 0.);

  for(int w = 0; w < well.size(); ++w) {
    
    for(int e = 0; e < well[w].size(); ++e) {
      //
      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	IO::log << IO::log_offset << funame << "well is not in the kinetic basis space\n";

	throw Error::Logic();
      }

      dtemp = 0.;
      
      for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	//
	dtemp += state_vector[well_shift[e] + k] * kinetic_basis[e].eigenvector(i->second, k);

      dtemp *= std::sqrt(well[w].states(e) * thermal_factor(e));

      res[w] += dtemp;
    }
  }
  
  return res;
}

// escape flux through the outer barrier
//
Lapack::Vector MasterEquation::ReactiveComplex::escape_flux (Lapack::Vector state) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::escape_flux: ";

  double dtemp;
  int    itemp;
  
  if(state.size() != global_size) {
    //
    IO::log << IO::log_offset << funame << "wrong state vector dimension: " << state.size() << "/" << global_size << "\n";

    throw Error::Logic();
  }

  Lapack::Vector res(bimolecular_size(), 0.);

  double dd_val;
    
  for(int b = 0; b < outer_barrier.size(); ++b) {
    //
    const int w = outer_barrier[b].connect.first;
      
    const int p = outer_barrier[b].connect.second;

    dd_val = 0.;
    
#pragma omp parallel for default(shared) private(dtemp) reduction(+: dd_val) schedule(dynamic)

    for(int e = 0; e < outer_barrier[b].size(); ++e) {
      //
      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	IO::log << IO::log_offset << funame << "well is not in the kinetic basis space\n";
	
	throw Error::Logic();
      }

      dtemp = 0.;
	
      for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	//
	dtemp += state[well_shift[e] + k] * kinetic_basis[e].eigenvector(i->second, k);
      }

      dtemp *= outer_barrier[b].states(e) * thermal_factor_sqrt(e) / well[w].states_sqrt(e);

      dd_val += dtemp;
    }

    dd_val /= 2. * M_PI;
    
    res[p] += dd_val;
  }

  for(int b = 0;  b < outer_barrier.size(); ++b) {
    //
    const int w = outer_barrier[b].connect.first;
		
    const int p = outer_barrier[b].connect.second;

    dd_val = 0.;
		
#pragma omp parallel for default(shared) private(dtemp, itemp) reduction(+: dd_val) schedule(dynamic)

    for(int e = 0; e < outer_barrier[b].size(); ++e) {

      std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);
		  
      if(i == kinetic_basis[e].well_index_map.end()) {
	//
	IO::log << IO::log_offset << funame << "high eigenstate correction: well w is not in the e map\n";

	throw Error::Logic();
      }
	
      for(int w1 = 0; w1 < well_size(); ++w1) {
	//
	std::map<int, int>::const_iterator i1 = kinetic_basis[e].well_index_map.find(w1);

	if(i1 == kinetic_basis[e].well_index_map.end())
	  //
	  continue;

	double proj = 0.;

	for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
	  //
	  proj += kinetic_basis[e].eigenvector(i1->second, k)
	    //
	    * kinetic_basis[e].eigenvector(i->second, k)
	    //
	    / kinetic_basis[e].eigenvalue[k];

	proj *= outer_barrier[b].states(e) / well[w].states_sqrt(e) * thermal_factor_sqrt(e);

	itemp = e + well[w1].kernel_bandwidth;

	const int e1_max = itemp < well[w1].size() ? itemp : well[w1].size();
	  
	itemp = e - well[w1].kernel_bandwidth + 1;

	const int e1_min = itemp > 0 ? itemp : 0;
	    
	for(int e1 = e1_min; e1 < e1_max; ++e1) {
	  //
	  std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);

	  if(i11 == kinetic_basis[e1].well_index_map.end()) {
	    //
	    IO::log << IO::log_offset << funame << "high eigenstate correction: well w1 is not in the e1 map\n";

	    throw Error::Logic();
	  }
	  
	  dtemp = 0.;
      
	  for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
	    //
	    dtemp += state[well_shift[e1] + k] * kinetic_basis[e1].eigenvector(i11->second, k);
	      
	  dtemp *= proj * well[w1].kernel(e1, e) * std::sqrt(well[w1].states(e1) * thermal_factor(e1 - e) / well[w1].states(e));
	  
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

Lapack::Vector MasterEquation::ReactiveComplex::propagate (Lapack::Vector state, Lapack::Vector escape) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::propagate: ";

  if(state.size() != global_size) {
    //
    IO::log << IO::log_offset << funame << "wrong state vector dimension: " << state.size() << "\n";

    throw Error::Init();
  }
  
  double         dtemp;
  int            itemp;
  Lapack::Vector vtemp;

  int      tout = time_output / time_step;
  
  int      tmax = time_limit  / time_step;

  double  tstep = time_step   / well[0].collision_frequency();

  if(bimolecular_size() && (!escape.isinit() || escape.size() != bimolecular_size())) {
    //
    escape.resize(bimolecular_size());

    escape = 0.;
  }
  
  if(tout > 0) {
    //
    IO::log << IO::log_offset << "branching fractions:\n"
      //
	    << IO::log_offset << std::setw(10) << "time";

    for(int w = 0; w < well.size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << well[w].name();

    for(int p = 0; p < bimolecular_size(); ++p)
      //
      IO::log << std::setw(Model::log_precision + 7) << bimolecular(p).name();
    
    IO::log << std::setw(Model::log_precision + 7) << "Total" << "\n";
  }
  
  for(int t = 0; t < tmax; ++t) {
    //
    // log output
    //
    if(tout > 0 && !(t % tout)) {
      //
      vtemp = population(state);

      IO::log << IO::log_offset << std::setw(10) << t * time_step;

      for(int w = 0; w < well.size(); ++w)
	//
	IO::log << std::setw(Model::log_precision + 7) << vtemp[w];

      for(int p = 0; p < bimolecular_size(); ++p)
	//
	IO::log << std::setw(Model::log_precision + 7) << escape[p];
    
      dtemp = sum(vtemp);

      if(bimolecular_size())
	//
	dtemp += sum(escape);

      IO::log << std::setw(Model::log_precision + 7) << dtemp << "\n";
    }
    
    // escape population contribution
    //
    if(bimolecular_size()) {
      //
      vtemp = escape_flux(state);

      vtemp *= tstep;

      escape += vtemp;
    }
    
    // second order time-propagation
    //
    vtemp = kin_mat * state;

    vtemp *= - tstep / 2.;

    vtemp += state;

    vtemp = kin_mat * vtemp;

    vtemp *= tstep;

    state -= vtemp;
  }
  
  vtemp = population(state);

  return vtemp;  
}

Lapack::Vector MasterEquation::ReactiveComplex::rate_partition ()
{
  const char funame [] = "MasterEquation::rate_partition: ";

  IO::Marker funame_marker(funame);
  
  double dtemp;
  int    itemp;
  bool   btemp;

  if(!kinetic_basis.size()) {
    //
    IO::log << IO::log_offset << "kinetic basis is not initialized\n";

    throw Error::Init();
  }
  
  Lapack::Vector state(global_size);

  state = 0.;

  // initial distribution: thermal one at reference energy
  //
  for(int k = 0; k < kinetic_basis[0].active_size; ++k) {
    //
    dtemp = 0.;
    
    for(int i = 0; i < kinetic_basis[0].size(); ++i) {
      //
      const int w = kinetic_basis[0].index_well_map[i];
      
      dtemp +=  kinetic_basis[0].eigenvector(i, k) * std::sqrt(well[w].states(0));
    }
    
    state[k] = dtemp;
  }
  
  Lapack::Vector res = propagate(state);

  dtemp = sum(res);

  res /= dtemp;
  
  double weight = 0.;

  for(int w = 0; w < well_size(); ++w)
    //
    weight += well[w].real_weight();

  IO::log << IO::log_offset << "real vs. thermal fractions:\n";

  IO::log << IO::log_offset
	  << std::setw(Model::log_precision + 7) << "well"
    	  << std::setw(Model::log_precision + 7) << "real"
	  << std::setw(Model::log_precision + 7) << "therm"
	  << "\n";

  for(int w = 0; w < well_size(); ++w)
    //
    IO::log << IO::log_offset 
	    << std::setw(Model::log_precision + 7) << well[w].name()
	    << std::setw(Model::log_precision + 7) << res[w]
	    << std::setw(Model::log_precision + 7) << well[w].real_weight() / weight
	    << "\n";

  return res;
}

/***********************************************************************************************
 ********************************** REACTION RATE MANAGEMENT ***********************************
 ***********************************************************************************************/

double MasterEquation::default_energy_reference (const Model::ChemGraph& cg)
{
  const char funame [] = "MasterEquation::default_energy_reference: ";

  double dtemp;
  int    itemp;
  bool   btemp;

  double e, s, c;
    
  btemp = true;

  double res;
  
  for(std::set<int>::const_iterator b = cg.inner_set.begin(); b != cg.inner_set.end(); ++b) {
    //
    Model::inner_barrier(*b).esc_parameters(temperature, e, s, c);
      
    e += Model::inner_barrier(*b).ground() + dist_width_factor * temperature * std::sqrt(c);
    
    if(btemp || res < e) {
      //
      btemp = false;
        
      res = e;
    }
  }
  
  for(std::set<int>::const_iterator b = cg.outer_set.begin(); b != cg.outer_set.end(); ++b) {
    //
    Model::outer_barrier(*b).esc_parameters(temperature, e, s, c);
      
    e += Model::outer_barrier(*b).ground() + dist_width_factor * temperature * std::sqrt(c);
    
    if(btemp || res < e) {
      //
      btemp = false;
        
      res = e;
    }
  }
  
  return res;
}

void MasterEquation::get_rate_data (std::list<RateData>& rate_data_list, Model::ChemGraph chem_graph, int eref, int flags)
{
  const char funame [] = "MasterEquation::get_rate_data:";

  double dtemp;
  bool   btemp;
  int    itemp;

  int root = 0;
  
  { IO::Marker marker("get_rate_data");

    IO::log << IO::log_offset << "pressure = ";
  
    switch(pressure_unit) {
      //
    case BAR:
      //
      IO::log << pressure / Phys_const::bar << " bar";
    
      break;
      //
    case TORR:
      //
      IO::log << pressure / Phys_const::tor << " torr";
    
      break;
      //
    case ATM:
      //
      IO::log << pressure / Phys_const::atm << " atm";
    
      break;
    }
  
    IO::log << "\t  temperature = " << (int)std::floor(temperature / Phys_const::kelv + 0.5) << " K\n";

    // initialization to full reactive complex
    //
    if(!chem_graph.well_set.size()) {
      //
      root = 1;
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	chem_graph.well_set.insert(w);
    
      for(int b = 0; b < Model::inner_barrier_size(); ++b)
	//
	chem_graph.inner_set.insert(b);
    
      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	chem_graph.outer_set.insert(b);

      eref = default_energy_reference(chem_graph) / energy_step();
    }

    if(!chem_graph.outer_set.size())
      //
      eref = default_energy_reference(chem_graph) / energy_step();

    ReactiveComplex rc(chem_graph, eref);
  
    // actual rate calculation
    //
    RateData rate_data;
	       
    rc.well_reduction_method(rate_data);

    // isomerization processes
    //
    if(!rc.outer_size()) {
      //
      if(rc.full_energy_range) {
	//
	if(rate_data.well_partition.size() > 1)
	  //
	  rate_data_list.push_back(rate_data);

	return;
      }
      
      if(rate_data.well_partition.size() > 1) {
	//
	rate_data_list.push_back(rate_data);

	for(int i = 0; i < rate_data.well_partition.size(); ++i)
	  //
	  if(rate_data.well_partition[i].size() > 1) {
	    //
	    IO::log << "\n" << IO::log_offset << "...continued..." << "\n\n";
	
	    get_rate_data(rate_data_list, Model::ChemGraph(rate_data.well_partition[i]));
	  }
      }
      else if(rc.xg.size() > 1) {
	//
	IO::log << "\n" << IO::log_offset << "...continued..." << "\n\n";
	
	get_rate_data(rate_data_list, Model::ChemGraph(rc.index_to_well(rc.xg)));
      }
      else if(!rc.xg.size()){
	//
	IO::log << IO::log_offset << "well partitioning error: there should be at least two well groups: "
	  //
		<< "increase WellCutoff parameter: current value: " << well_cutoff << "\n";

	throw Error::Run();
      }

      return;
    }

    rate_data.index_bim_map.resize(rc.index_bim_map().size());
    
    rate_data.index_bim_map = rc.index_bim_map();
    
    rate_data_list.push_back(rate_data);

    // all wells are included
    //
    if(rc.full_energy_range) {
      //
      if(root)
	//
	print_rate_data(rate_data_list);
      
      return;
    }

    // additional izomerization rates
    //
    for(int i = 0; i < rate_data.well_partition.size(); ++i)
      //
      if(rate_data.well_partition[i].size() > 1) {
	//
	IO::log << "\n" << IO::log_offset << "...continued..." << "\n\n";
	
	get_rate_data(rate_data_list, Model::ChemGraph(rate_data.well_partition[i]));
      }

    // no uncoupled low energy wells
    //
    if(!rc.xg.size()) {
      //
      if(root)
	//
	print_rate_data(rate_data_list);

      return;
    }
    
    for(int w = 0; w < rc.well_size(); ++w)
      //
      if(rc.xg.find(w) == rc.xg.end())
	//
	chem_graph.remove(std::make_pair(Model::WELL, rc.well[w].index()));

    eref -= rc.kinetic_basis.size();
    
    dtemp = eref * energy_step();

    for(int b = 0; b < rc.inner_size(); ++b) {
      //
      itemp = rc.inner_barrier[b].index();

      if(Model::inner_barrier(itemp).ground() >= dtemp)
	//
	chem_graph.inner_set.erase(itemp);
    }

    for(int b = 0; b < rc.outer_size(); ++b) {
      //
      itemp = rc.outer_barrier[b].index();

      if(Model::outer_barrier(itemp).ground() >= dtemp)
	//
	chem_graph.outer_set.erase(itemp);
    }

    std::list<Model::ChemGraph> chem_graph_list = chem_graph.factorize();

    for(std::list<Model::ChemGraph>::const_iterator g = chem_graph_list.begin(); g != chem_graph_list.end(); ++g) {
      //
      if(!g->outer_set.size()) {
	//
	IO::log << IO::log_offset << "no outer barrier in the graph\n";

	throw Error::Logic();
      }

      IO::log << "\n" << IO::log_offset << "...continued..." << "\n\n";
	
      get_rate_data(rate_data_list, *g, eref);
    }
  }

  if(root)
    //
    print_rate_data(rate_data_list);
}

void MasterEquation::check_rate_data (const std::list<RateData>& rate_data_list)
{
  const char funame [] = "MasterEquation::check_rate_data: ";
  
  for(std::list<RateData>::const_iterator r = rate_data_list.begin(); r != rate_data_list.end(); ++r) {
    //
    const int ws = r->well_partition.size();

    const int bs = r->index_bim_map.size();
      
    if(bs && (!r->bb_rate.isinit() || r->bb_rate.size() != bs)) {
      //
      IO::log << IO::log_offset << funame << "bimolecular-to-bimolecular rate data not initialized properly\n";

      throw Error::Init();
    }
      
    if(ws && (!r->ww_rate.isinit() || r->ww_rate.size() != ws)) {
      //
      IO::log << IO::log_offset << funame << "well-to-well rate data not initialized properly\n";

      throw Error::Init();

    }

    if(ws && bs && (!r->wb_rate.isinit() || r->wb_rate.size1() != ws || r->wb_rate.size2() != bs ||
		    //
		    !r->bw_rate.isinit() || r->bw_rate.size1() != bs || r->bw_rate.size2() != ws)) {
      //
      IO::log << IO::log_offset << funame << "well-to-bimolecular rate data not initialized properly\n";

      throw Error::Init();
    }
  }
}

void MasterEquation::print_rate_data (const std::list<RateData>& rate_data_list)
{
  const char funame [] = "MasterEquation::print_rate_data: ";

  check_rate_data(rate_data_list);

  int itemp;
  
  // local output
  //
  itemp = 0;
    
  for(std::list<RateData>::const_iterator rit = rate_data_list.begin(); rit != rate_data_list.end(); ++rit)
    //
    if(rit->well_partition.size() > 1) {
      //
      itemp = 1;

      break;
    }

  IO::log << "\n";
	       
  if(itemp) {
    //
    IO::log << IO::log_offset << "well-to-well rates:\n\n";

    for(std::list<RateData>::const_iterator rit = rate_data_list.begin(); rit != rate_data_list.end(); ++rit) {
      //
      if(rit->well_partition.size() > 1) {
	//
	// well name width
	//
	int welw = 0;

	for(int g = 0; g < rit->well_partition.size(); ++g) {
	  //
	  itemp = Model::group_name(rit->well_partition[g]).size();
	  
	  if(!g || itemp > welw)
	    //
	    welw = itemp;
	}

	welw += 2;

	itemp = Model::log_precision + 7;

	welw = welw > itemp ? welw: itemp;
	
	IO::log << IO::log_offset << std::setw(welw) << "W\\W";

	for(int g = 0; g < rit->well_partition.size(); ++g)
	  //
	  IO::log << std::setw(welw) << Model::group_name(rit->well_partition[g]);

	IO::log << "\n";

	for(int g = 0; g < rit->well_partition.size(); ++g) {
	  //
	  IO::log << IO::log_offset << std::setw(welw) << Model::group_name(rit->well_partition[g]);

	  for(int f = 0; f < rit->well_partition.size(); ++f)
	    //
	    IO::log << std::setw(welw) << rit->ww_rate(g, f);

	  IO::log << "\n";
	}//

	IO::log << "\n";
	//
      }// need output
      //
    }// data cycle
    //
  }// need well-to-well rates
    
  itemp = 0;
    
  for(std::list<RateData>::const_iterator rit = rate_data_list.begin(); rit != rate_data_list.end(); ++rit)
    //
    if(rit->well_partition.size() && rit->index_bim_map.size()) {
      //
      itemp = 1;

      break;
    }

  if(itemp) {
    //
    IO::log << IO::log_offset << "well-to-bimolecular and bimolecular-to-well rates:\n\n";
      
    for(std::list<RateData>::const_iterator rit = rate_data_list.begin(); rit != rate_data_list.end(); ++rit) {
      //
      if(rit->well_partition.size() && rit->index_bim_map.size()) {
	//
	// well name width
	//
	int welw = 0;

	for(int g = 0; g < rit->well_partition.size(); ++g) {
	  //
	  itemp = Model::group_name(rit->well_partition[g]).size();
	  
	  if(!g || itemp > welw)
	    //
	    welw = itemp;
	}

	welw += 2;

	itemp = Model::log_precision + 7;

	welw = welw > itemp ? welw: itemp;

	// bimolecular name width
	//
	int bimw = 0;

	for(int p = 0; p < rit->index_bim_map.size(); ++p) {
	  //
	  itemp = Model::bimolecular(rit->index_bim_map[p]).name().size();
	  
	  if(!p || itemp > bimw)
	    //
	    bimw = itemp;
	}
	
	bimw += 2;

	itemp = Model::log_precision + 7;

	bimw = bimw > itemp ? bimw: itemp;
	
	IO::log << IO::log_offset << std::setw(bimw) << "B\\W";

	for(int g = 0; g < rit->well_partition.size(); ++g)
	  //
	  IO::log << std::setw(welw) << Model::group_name(rit->well_partition[g]);

	IO::log << "\n";

	for(int p = 0; p < rit->index_bim_map.size(); ++p) {
	  //
	  IO::log << IO::log_offset << std::setw(bimw) << Model::bimolecular(rit->index_bim_map[p]).name();

	  for(int g = 0; g < rit->well_partition.size(); ++g)
	    //
	    IO::log << std::setw(welw) << rit->bw_rate(p, g);

	  IO::log << "\n";
	}

	IO::log << "\n" << IO::log_offset << std::setw(welw) << "W\\B";

	for(int p = 0; p < rit->index_bim_map.size(); ++p)
	  //
	  IO::log << std::setw(bimw) << Model::bimolecular(rit->index_bim_map[p]).name();

	IO::log << "\n";
	  
	for(int g = 0; g < rit->well_partition.size(); ++g) {
	  //
	  IO::log << IO::log_offset << std::setw(welw) << Model::group_name(rit->well_partition[g]);

	  for(int p = 0; p < rit->index_bim_map.size(); ++p)
	    //
	    IO::log << std::setw(bimw) << rit->wb_rate(g, p);
	    
	  IO::log << "\n";
	}
	  
	IO::log << "\n";
	//
      }// need data
      //
    }// data cycle
    //
  }// need bw data 
}

void MasterEquation::aggregate_rate_data (const std::list<RateData>& rate_data_list, RateData& rd)
{
  const char funame [] = "MasterEquation::aggregate_rate_data: ";

  int itemp;

  check_rate_data(rate_data_list);
  
  itemp = 0;
    
  for(std::list<RateData>::const_iterator r = rate_data_list.begin(); r != rate_data_list.end(); ++r) {
    //
    for(int g = 0; g < r->well_partition.size(); ++g) {
      //
      std::map<std::set<int>, int>::const_iterator i = rd.group_index_map.find(r->well_partition[g]);

      if(i != rd.group_index_map.end()) {
	//
	IO::log << IO::log_offset << funame << "duplicated group: " << Model::group_name(r->well_partition[g]) << "\n";

	throw Error::Logic();
      }
	
      rd.group_index_map[r->well_partition[g]] = itemp++;
    }
  }

  rd.ww_rate.resize(itemp);

  rd.ww_rate = 0.;
    
  rd.wb_rate.resize(itemp, Model::bimolecular_size());

  rd.wb_rate = 0.;
    
  rd.bw_rate.resize(Model::bimolecular_size(), itemp);

  rd.bw_rate = 0.;
    
  rd.bb_rate.resize(Model::bimolecular_size());

  rd.bb_rate = 0.;
    
  for(std::list<RateData>::const_iterator r = rate_data_list.begin(); r != rate_data_list.end(); ++r) {
    //
    const int ws = r->well_partition.size();

    const int bs = r->index_bim_map.size();
      
    for(int i = 0; i < ws; ++i)
      //
      for(int j = 0; j < ws; ++j)
	//    
	rd.ww_rate(rd.group_index_map[r->well_partition[i]], rd.group_index_map[r->well_partition[j]]) = r->ww_rate(i, j);

    for(int i = 0; i < ws; ++i)
      //
      for(int j = 0; j < bs; ++j) {
	//    
	rd.wb_rate(rd.group_index_map[r->well_partition[i]], r->index_bim_map[j]) = r->wb_rate(i, j);
	  
	rd.bw_rate(r->index_bim_map[j], rd.group_index_map[r->well_partition[i]]) = r->bw_rate(j, i);
      }

    for(int i = 0; i < bs; ++i)
      //
      for(int j = i; j < bs; ++j)
	//
	rd.bb_rate(r->index_bim_map[i], r->index_bim_map[j]) += r->bb_rate(i, j);
  }
}

/********************************************************************************************
 ********************************** SETTING KINETIC BASIS ***********************************
 ********************************************************************************************/

void MasterEquation::ReactiveComplex::set_kinetic_matrix ()
{    
  const char funame [] = "MasterEquation::ReactiveComplex::set_kinetic_matrix: ";

  //IO::Marker marker("setting_kinetic_matrix");

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  Lapack::Vector      vtemp;
  group_t             gtemp;

  for(group_t::const_iterator w = well_set.begin(); w != well_set.end(); ++w)
    //
    if(w == well_set.begin() || Model::well(*w).ground() < dtemp)
      //
      dtemp = Model::well(*w).ground();

  int ener_index_max = _ener_index_begin - (int)std::floor(dtemp / energy_step());
  
  for(int w = 0; w < well_size(); ++w)
    //
    well[w].resize(ener_index_max);
  
  for(int b = 0; b < inner_size(); ++b)
    //
    inner_barrier[b].resize(ener_index_max);
  
  for(int b = 0; b < outer_size(); ++b)
    //
    outer_barrier[b].resize(ener_index_max);

  kinetic_basis.resize(ener_index_max);

  /**********************************************************************************
   **************************** KINETICALLY ACTIVE BASIS ****************************
   **********************************************************************************/

  // equilibrated well group - highest energy map
  //
  typedef std::map<std::set<int>, int> gem_t;

  gem_t gem;

  int ener_index_active = -1;
  
  for(int e = 0; e < ener_index_max; ++e) {// energy cycle
    //
    const double ener = energy_bin(e);

    std::vector<int>   index_well;
     
    std::map<int, int> well_index;
      
    std::vector<float_t> weight;

    for(int w = 0; w < well_size(); ++w) {
      //
      if(e >= well[w].size()) 
	//
	continue;

      itemp = well[w].index();
      
      if(ener < Model::well(itemp).ground() + energy_step() / 2.) {
	//
	well[w].resize(e);

	continue;
      }

      dtemp = Model::well(itemp).states(ener);

      if(dtemp <= 0.) {
	//
	IO::log << IO::log_offset << "non-positive density of states\n";

	throw Error::Range();
      }
    
      well[w].states(e, dtemp);
      
      well_index[w] = index_well.size();
	  
      index_well.push_back(w);

      weight.push_back((float_t)dtemp);
    }
  
    const int index_size = index_well.size();

    if(!index_size) {
      //
      full_energy_range = 1;
      
      kinetic_basis.resize(e);

      break;
    }

    std::vector<float_t> weight_sqrt(index_size);

    for(int i = 0; i < index_size; ++i)
      //
      weight_sqrt[i] = sqrt(weight[i]);
    
    kinetic_basis[e].well_index_map = well_index;
      
    kinetic_basis[e].index_well_map = index_well;

    // microcanonical kinetic matrix
    //
    Matrix<float_t> km(index_size);
      
    km = (float_t)0.;
  
    // inner barrier contribution
    //
    for(int b = 0; b < inner_size(); ++b) {
      //
      if(e >= inner_barrier[b].size())
	//
	continue;

      itemp = inner_barrier[b].index();
      
      if(ener < Model::inner_barrier(itemp).ground()) {
	//
	inner_barrier[b].resize(e);

	continue;
      }

      const float_t ibs = Model::inner_barrier(itemp).states(ener);

      if(ibs <= 0.) {
	/*
	  IO::log << IO::log_offset << "WARNING: " << inner_barrier[b].name() << ": non-positive inner number of states = " << ibs
	  //
	  << " at energy = " << std::floor(energy_bin(e) / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol: truncating\n";
	*/
	inner_barrier[b].resize(e);

	continue;
      }

      inner_barrier[b].states(e, (double)ibs);
      
      const int w1 = inner_barrier[b].connect.first;
	  
      const int w2 = inner_barrier[b].connect.second;
	  
      std::map<int, int>::const_iterator p1 = well_index.find(w1);
	  
      std::map<int, int>::const_iterator p2 = well_index.find(w2);
      
      if(p1 == well_index.end()) {
	/*
	  IO::log << IO::log_offset << "WARNING: " << inner_barrier[b].name() << " inner barrier energy[kcal/mol] = "
	  //
	  << std::floor(ener / Phys_const::kcal * 10. + 0.5) / 10.
	  //
	  << " is lower than connected " << well[w1].name() << " well ground = "
	  //
	  << std::floor(well[w1].ground() / Phys_const::kcal * 10. + 0.5) / 10.
	  //
	  << ": truncating\n";
	*/
	inner_barrier[b].resize(e);
	
	continue;
      }
	
      if(p2 == well_index.end()) {
	/*
	  IO::log << IO::log_offset << "WARNING: " << inner_barrier[b].name() << " inner barrier energy[kcal/mol] = "
	  //
	  << std::floor(ener / Phys_const::kcal * 10. + 0.5) / 10.
	  //
	  << " is lower than connected " << well[w2].name() << " well ground = "
	  //
	  << std::floor(well[w2].ground() / Phys_const::kcal * 10. + 0.5) / 10.
	  //
	  << ": truncating\n";
	*/
	inner_barrier[b].resize(e);
	
	continue;
      }

      ibs /= 2. * M_PI;
      
      km(p1->second, p2->second) -= ibs / weight_sqrt[p1->second] / weight_sqrt[p2->second];

      km(p2->second, p1->second) = km(p1->second, p2->second);
      
      km(p1->second, p1->second) += ibs / weight[p1->second];
	  
      km(p2->second, p2->second) += ibs / weight[p2->second];
      //
    }// inner barrier
  
    // outer barrier contribution
    //
    for(int b = 0; b < outer_size(); ++b) {
      //
      if(e >= outer_barrier[b].size())
	//
	continue;

      itemp = outer_barrier[b].index();

      if(ener < Model::outer_barrier(itemp).ground()) {
	//
	outer_barrier[b].resize(e);

	continue;
      }

      const float_t obs = Model::outer_barrier(itemp).states(ener);

      if(obs <= 0.) {
	/*
	  IO::log << IO::log_offset << "WARNING: " << outer_barrier[b].name() << ": non-positive outer number of states = " << obs
	  //
	  << " at energy = " << std::floor(energy_bin(e) / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol: truncating\n";
	*/
	outer_barrier[b].resize(e);

	continue;
      }

      outer_barrier[b].states(e, (double)obs);
      
      const int w = outer_barrier[b].connect.first;

      std::map<int, int>::const_iterator p = well_index.find(w);

      if(p == well_index.end()) {
	/*
	  IO::log << IO::log_offset << "WARNING: " << outer_barrier[b].name() << " outer barrier energy[kcal/mol] = "
	  //
	  << std::floor(ener / Phys_const::kcal * 10. + 0.5) / 10.
	  //
	  << " is lower than connected " << well[w].name() << " well ground = "
	  //
	  << std::floor(well[w].ground() / Phys_const::kcal * 10. + 0.5) / 10.
	  //
	  << ": truncating\n";
	*/
	outer_barrier[b].resize(e);
	  
	continue;
      }

      obs /= 2. * M_PI;
      
      km(p->second, p->second) += obs / weight[p->second];
      //
    }// outer barrier
  
    // relaxation eigenvalues
    //
    Matrix<float_t> evec(index_size);

    Vector<float_t> eval(index_size);

    Mpack::eigenvalues<float_t>(km, eval, &evec);
    
    kinetic_basis[e].eigenvalue  = (Lapack::Vector)eval;
      
    kinetic_basis[e].eigenvector = (Lapack::Matrix)evec;
      
    // kinetically active subspace
    //
    for(itemp = 0; itemp < index_size; ++itemp)
      //
      if((double)eval[itemp] > well[0].collision_frequency() * reduction_threshold)
	//
	break;

    const int active_size = itemp;
    
    kinetic_basis[e].active_size = itemp;

    if(!e && (outer_size() && active_size || !outer_size() && active_size > 1)) {
      //
      IO::log << IO::log_offset << "WARNING: there are kinetically active processes at the highest energy: increase ThermalWidthFactor (default 5)\n";

      IO::log << IO::log_offset << "eigenvalues over collision_frequency:";

      for(int i = 0; i < index_size; ++i)
	//
	IO::log << "  " << (double)eval[i] / well[0].collision_frequency();

      IO::log << "\n\n";
    }

    if((outer_size() && active_size || !outer_size() && active_size > 1) && ener_index_active < 0) {
      //
      IO::log << IO::log_offset << "kinetically active processes first appear: energy = "
	//
	      << std::floor(energy_bin(e) / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol\n";
      
      IO::log << IO::log_offset << "eigenvalues over collision_frequency:";

      for(int i = 0; i < index_size; ++i)
	//
	IO::log << "  " << (double)eval[i] / well[0].collision_frequency();

      IO::log << "\n";

      ener_index_active = e;
    }
    
    // no kinetically active processes
    //
    if(outer_size() && !active_size || !outer_size() && active_size == 1) {
      //
      // we are at the top
      //
      if(ener_index_active < 0)
	//
	continue;

      if(index_well.size() == well_size()) {
	//
	IO::log << IO::log_offset << "the only reason for kinetically active process to disappear is that some well is gone at lower energy\n";

	throw Error::Logic();
      }
      
      for(int i = 0; i < index_size; ++i)
	//
	xg.insert(index_well[i]);
      
      IO::log << "\n" << IO::log_offset << "energy = " << std::floor(ener / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol\n";
      
      IO::log << IO::log_offset << "eigenvalues over collision_frequency:";

      for(int i = 0; i < index_size; ++i)
	//
	IO::log << "  " << (double)eval[i] / well[0].collision_frequency();

      IO::log << "\n";

      IO::log << "\n" << IO::log_offset << "all isomerization processes are fast => stop\n";

      kinetic_basis.resize(e);

      break;
    }

    // number of low eigenvalues
    //
    for(itemp = 0; itemp < index_size; ++itemp)
      //
      if((double)eval[itemp] > well[0].collision_frequency() * reduction_tolerance)
	//
	break;

    if(itemp != index_size && itemp != 0)
      //
      for(; itemp != 0; --itemp)
	//
	if(eval[itemp - 1] <= 0. || eval[itemp] / eval[itemp - 1] > (float_t)reduction_increment)
	  //
	  break;
    
    const int spec_size = itemp;

    if(outer_size() && !spec_size || !outer_size() && spec_size == 1)
      //
      continue;

    IO::log << "\n" << IO::log_offset << "................\n\n";
    
    IO::log << IO::log_offset << "energy = " << std::floor(ener / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol\n";

    IO::log << IO::log_offset << "wells:";

    for(int i = 0; i < index_size; ++i)
      //
      IO::log << " " << well[index_well[i]].name();

    IO::log << "\n";
    
    IO::log << IO::log_offset << "eigenvalue over collision frequency:";

    for(int i = 0; i < index_size; ++i)
      //
      IO::log << "  " << (double)eval[i] / well[0].collision_frequency();

    IO::log << "\n";

    // low eigenvectors
    //
    Matrix<float_t> pop_chem(index_size, spec_size);

    for(int i = 0; i < spec_size; ++i)
      //
      pop_chem.column(i) = evec.column(i);

    // partitioning the wells
    //
    group_t     bg;
    
    partition_t wp;
    
    WellSet<float_t> well_set(weight);

    itemp = 0;
    
    if(!outer_size())
      //
      itemp = WellSet<float_t>::INTERNAL;
    
    dtemp = well_set.threshold_well_partition(pop_chem, wp, bg, itemp);

    IO::log << IO::log_offset << "wells partition error  = " << dtemp << "\n";

    // check that partitioned wells are connected
    //
    for(int g = 0; g < wp.size(); ++g) {
      //
      Model::ChemGraph cg(index_to_well(kinetic_basis[e].index_to_well(wp[g])));
      
      if(!cg.is_connected()) {
	//
	IO::log << IO::log_offset << Model::group_name(cg.well_set) << " is not connected\n";

	typedef std::list<Model::ChemGraph> gl_t;

	gl_t gl = cg.factorize();

	typedef std::multimap<double, std::set<int> > wm_t;

	wm_t wm;
	
	IO::log << IO::log_offset << "factorized group/weight:";

	for(gl_t::const_iterator i = gl.begin(); i != gl.end(); ++i) {
	  //
	  gtemp = kinetic_basis[e].well_to_index(well_to_index(i->well_set));

	  dtemp = 0.;
	  
	  for(std::set<int>::const_iterator w = gtemp.begin(); w != gtemp.end(); ++w)
	    //
	    dtemp += (double)weight[*w];
	  
	  wm.insert(std::make_pair(dtemp, gtemp));
	  
	  IO::log << "  " << Model::group_name(i->well_set) << "/" << dtemp;
	}

	IO::log << "\n";
	
	throw Error::Run();
      }
	
	/*
	if(outer_size()) {
	  //
	  IO::log << "\nmoving the low weight associated wells to bimolecular group\n";

	  for(wm_t::const_reverse_iterator i = wm.rbegin(); i != wm.rend(); ++i)
	    //
	    if(i != wm.rbegin())
	      //
	      for(std::set<int>::const_iterator w = i->second.begin(); w != i->second.end(); ++w) {
		//
		wp[g].erase(*w);

		bg.insert(*w);
	      }
	}
      }
      */
    }
    
    IO::log << IO::log_offset << "wells partition size =   " << wp.size() << "\n";

    IO::log << IO::log_offset << "wells partition:       ";

    for(int g = 0; g < wp.size(); ++g)
      //
      IO::log << "  " << name(kinetic_basis[e].index_to_well(wp[g]));

    IO::log << "\n";

    if(outer_size()) {
      //
      IO::log << IO::log_offset << "bimolecular group size = " << bg.size() << "\n";
    
      if(bg.size())
	//
	IO::log << IO::log_offset << "bimolecular group:       " << name(kinetic_basis[e].index_to_well(bg)) << "\n";
    }
    
    // new bound groups
    //
    itemp = gem.size();
    
    for(int g = 0; g < wp.size(); ++g) {
      //
      btemp = true;
      
      gtemp = kinetic_basis[e].index_to_well(wp[g]);

      for(gem_t::const_iterator g = gem.begin(); g != gem.end(); ++g) {
	//
	if((gtemp & g->first).size()) {
	  //
	  btemp = false;

	  break;
	}
      }
      
      if(btemp) {
	//
	gem[gtemp] = e;
	
	IO::log << IO::log_offset << "new bound group found:   " << name(gtemp) << "\n";
      }
    }

    if(itemp != gem.size())
      //
      continue;
    
    // check if all bound groups have sufficient depth
    //
    double wd;
    
    for(gem_t::const_iterator g = gem.begin(); g != gem.end(); ++g) {
      //
      dtemp = (e - g->second) * energy_step_over_temperature;
      
      if(g == gem.begin() || dtemp < wd)
	//
	wd = dtemp;
    }

    IO::log << IO::log_offset << "minimal well depth     = " << wd << " (vs. " << well_cutoff << ")\n";

    if(wd < well_cutoff)
      //
      continue;

    if(!outer_size()) {
      //
      IO::log << "\n" << IO::log_offset << "required depth has been reached => stop\n";

      if(index_well.size() != well_size())
	//
	for(int i = 0; i < index_size; ++i)
	  //
	  xg.insert(index_well[i]);
      
      kinetic_basis.resize(e);

      break;
    }
    
    if(!bg.size()) {
      //
      IO::log << "\n" << IO::log_offset << "no bimolecular group => stop\n";

      kinetic_basis.resize(e);

      break;
    }
      
    Matrix<float_t> tkm(bg.size());

    itemp = 0;
    
    for(group_t::const_iterator i = bg.begin(); i != bg.end(); ++i, ++itemp) {
      //
      int count = 0;
      
      for(group_t::const_iterator j = bg.begin(); j != bg.end(); ++j, ++count)
	//
	tkm(itemp, count) = km(*i, *j);
    }

    Vector<float_t> teval(bg.size());

    Mpack::eigenvalues(tkm, teval);
    
    for(itemp = 0; itemp < bg.size();++itemp)
      //
      if(teval[itemp] > well[0].collision_frequency() * reduction_threshold)
	//
	break;
    
    IO::log << IO::log_offset << "active subspace size   = " << itemp << "\n";

    if(!itemp) {
      //
      IO::log << "\n" << IO::log_offset << "no kinetically active processes => stop\n";
      
      xg = kinetic_basis[e].index_to_well(bg);
      
      kinetic_basis.resize(e);
      
      break;
    }
    
  }//energy cycle

  if(kinetic_basis.size() == ener_index_max) {
    //
    full_energy_range = 1;
  }
  else
    //
    ener_index_max = kinetic_basis.size();

  if(full_energy_range) {
    //
    IO::log << IO::log_offset << "energy = " << std::floor(energy_bin(ener_index_max) / Phys_const::kcal * 10. + 0.5) / 10. << " kcal/mol\n";

    IO::log << "\n" << IO::log_offset << "full energy range has been reached => stop\n";
  }

  IO::log << "\n" << IO::log_offset << "****************************************\n\n";
  
  IO::log << IO::log_offset << "kinetic basis initialization done\n";

  for(int w = 0; w < well_size(); ++w) {
    //
    if(well[w].size() > ener_index_max)
      //
      well[w].resize(ener_index_max);

    well[w].set();

    if(!w || kernel_bandwidth_max < well[w].kernel_bandwidth)
      //
      kernel_bandwidth_max = well[w].kernel_bandwidth;
  }

  for(int b = 0; b < inner_size(); ++b) {
    //
    if(inner_barrier[b].size() > ener_index_max)
      //
      inner_barrier[b].resize(ener_index_max);

    inner_barrier[b].set();
  }
  
  for(int b = 0; b < outer_size(); ++b) {
    //
    if(outer_barrier[b].size() > ener_index_max)
      //
      outer_barrier[b].resize(ener_index_max);

    outer_barrier[b].set();
  }

  IO::log << "\n" << IO::log_offset << "energy range over temperature = "
    //
	  << std::floor((ener_index_max - ener_index_active) * energy_step() / temperature + 0.5) << "\n\n";
  
  /************************************************************************************* 
   ********************************* KINETIC MATRIX ************************************
   *************************************************************************************/

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
    itemp += well[w].size();
    
  IO::log << IO::log_offset << "original kinetic matrix size = " << itemp << "\n";
  
  IO::log << IO::log_offset << "reduced  kinetic matrix size = " << global_size << "\n";
  
  if(!global_size)
    //
    return;
  
  kin_mat.resize(global_size);
  
  kin_mat = 0.;
  
  // add (slow) isomerization relaxation
  //
  for(int e = 0; e < ener_index_max; ++e) {
    //
    for(int l = 0; l < kinetic_basis[e].active_size; ++l) {
      //
      itemp = well_shift[e] + l;
      
      kin_mat(itemp, itemp) = kinetic_basis[e].eigenvalue[l];
    }
  }
    
  // collision relaxation
  //
  int cycle_size = ener_index_max * ener_index_max;
  
#pragma omp parallel for default(shared) private(itemp, dtemp) schedule(static)

  for(int cycle = 0; cycle < cycle_size; ++cycle) {
    //
    int e1 = cycle / ener_index_max;

    int e2 = cycle % ener_index_max;

    if(e1 > e2 || e2 - e1 >= kernel_bandwidth_max)
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
	  
	for(int w = 0; w < well_size(); ++w) {
	  //
	  std::map<int, int>::const_iterator i1 = kinetic_basis[e1].well_index_map.find(w);

	  std::map<int, int>::const_iterator i2 = kinetic_basis[e2].well_index_map.find(w);
	    
	  itemp = e2 - e1;
	    
	  if(i1 != kinetic_basis[e1].well_index_map.end() &&
	     //
	     i2 != kinetic_basis[e2].well_index_map.end()) {
	    //
	    dtemp +=  well[w].kernel(e1, e2) * well[w].collision_frequency()
	      //
	      * std::sqrt(well[w].states(e1) / well[w].states(e2) / thermal_factor(e2 - e1))
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
  }// main cycle
}

/********************************************************************************************
 ********************************** WELL REDUCTION METHOD ***********************************
 ********************************************************************************************/

void MasterEquation::ReactiveComplex::well_reduction_method (RateData& rate_data, int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::well_reduction_method: ";

  IO::Marker marker("well_reduction_method");
  
  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  Lapack::Vector      vtemp;

  IO::log << IO::log_offset << "collision frequency = " << well[0].collision_frequency() / Phys_const::herz << " 1/sec\n";

  if(!kinetic_basis.size()) {
    //
    IO::log << IO::log_offset << "kinetic basis not initialized\n";

    throw Error::Init();
  }

  Lapack::SymmetricMatrix bb_rate;

  if(bimolecular_size())
    //
    bb_rate = fast_bb_rate();

  if(!global_size) {
    //
    IO::log << IO::log_offset << "only fast isomerization is available\n";

    if(bimolecular_size()) {
      //
      bb_rate /= std::exp(_ener_index_begin * energy_step_over_temperature);
  
      rate_data.bb_rate = bb_rate;
    }
    
    return;
  }
  
  /**************************************************************************************
   ************************* KINETIC MATRIX DIAGONALIZATION *****************************
   **************************************************************************************/

  Lapack::Vector eigenval;
  
  Lapack::Matrix global_eigen(global_size);

  { IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);
    //
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

  if(global_size > well_size()) {
    //
    dtemp = eigenval[well_size()];
  }
  else
    //
    dtemp = eigenval.back();
  
  const double relax_eval_min = dtemp;
  
  const double relax_eval_max = eigenval.back();

  IO::log << IO::log_offset << "minimal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_min / well[0].collision_frequency() << "\n";
  
  IO::log << IO::log_offset << "maximal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_max / well[0].collision_frequency() << "\n";

  /*****************************************************************************************************************
   ********************************************** SETTING GLOBAL MATRICES ******************************************
   *****************************************************************************************************************/
  
  // eigenvectors projection on thermal subspace;
  //
  Lapack::Matrix eigen_pop(global_size, well_size());

#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)
  
  for(int l = 0; l < global_size; ++l) {
    //
    for(int w = 0; w < well_size(); ++w) {
      //
      double prod = 0.;
      
      for(int e = 0; e < well[w].size(); ++e) {
	//
	std::map<int, int>::const_iterator wi = kinetic_basis[e].well_index_map.find(w);

	if(wi == kinetic_basis[e].well_index_map.end()) {
	  //
	  IO::log << IO::log_offset << "well is not in the kinetic basis space\n";

	  throw Error::Logic();
	}

	dtemp = 0.;
	
	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(wi->second, k);

	prod += dtemp * std::sqrt(well[w].states(e) * thermal_factor(e));
      }
      
      eigen_pop(l, w) = prod / well[w].weight_sqrt();
    }//
    //
  }// thermal
  
  // eigenvectors projection on the bimolecular subspace;
  //
  Lapack::Matrix eigen_bim;

  if(bimolecular_size()) {
    //
    eigen_bim.resize(global_size, bimolecular_size());

    eigen_bim = 0.;
  
#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)
  
    for(int l = 0; l < global_size; ++l) {
      //
      for(int b = 0; b < outer_barrier.size(); ++b) {
	//
	const int& w = outer_barrier[b].connect.first;
      
	const int& p = outer_barrier[b].connect.second;

	for(int e = 0; e < outer_barrier[b].size(); ++e) {
	  //
	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end())
	    //
	    continue;

	  dtemp = 0.;
	
	  for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	    //
	    dtemp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	  eigen_bim(l, p) += dtemp * outer_barrier[b].states(e) / 2. / M_PI
	    //
	    * std::sqrt(thermal_factor(e) / well[w].states(e));
	}//
	//
      }//
      //
    }//
    //
  }// bimolecular
  
  /***************************************************************
   **************** CHEMICAL SUBSPACE DIMENSION ******************
   ***************************************************************/

  std::vector<double> relaxation_projection(well_size(), 0.);

  if(global_size >=  well_size()) {
    //
    for(int l = 0; l < well_size(); ++l)
      //
      relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));


    // absolute eigenvalue threshold
    //
    if(chemical_threshold > 1.) {
      //
      for(itemp = 0; itemp < well_size(); ++itemp)
	//
	if(relax_eval_min / eigenval[itemp]  < chemical_threshold)
	  //
	  break;
    }
    // relaxation projection threshold
    //
    else if(chemical_threshold < 1. && chemical_threshold > 0.) {
      //
      for(itemp = 0; itemp < well_size(); ++itemp)
	//
	if(relaxation_projection[itemp] > chemical_threshold)
	  //
	  break;
    }
    // relative eigenvalue threshold
    //
    else if(chemical_threshold < -1. ) {
      //
      for(itemp = well_size(); itemp > 0; --itemp)
	//
	if(eigenval[itemp] / eigenval[itemp - 1] > - chemical_threshold)
	  //
	  break;
    }
    //
    else
      //
      itemp = well_size();
  }
  else if(bimolecular_size()) {
    //
    itemp = 0;
  }
  else
    //
    return;
	
  const int chem_size = itemp;

  /*******************************************************************
   ***************** CHEMICAL EIGENSTATES CORRECTION *****************
   *******************************************************************/

  for(int l = 0; l < chem_size; ++l) {
    //
    for(int w1 = 0; w1 < well_size(); ++w1) {
      //
      for(int e1 = 0; e1 < well[w1].size(); ++e1) {
	//
	std::map<int, int>::const_iterator w1i = kinetic_basis[e1].well_index_map.find(w1);
	  
	dtemp = 0.;
	
	for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
	  //
	  dtemp += global_eigen(well_shift[e1] + k, l) * kinetic_basis[e1].eigenvector(w1i->second, k);
	 
	const double proj = dtemp;

	itemp = e1 + well[w1].kernel_bandwidth;

	const int e2_max = itemp < well[w1].size() ? itemp : well[w1].size();
	  
	itemp = e1 - well[w1].kernel_bandwidth + 1;

	const int e2_min = itemp > 0 ? itemp : 0;
	    
	for(int e2 = e2_min; e2 < e2_max; ++e2) {
	  //
	  w1i = kinetic_basis[e2].well_index_map.find(w1);

	  // population correction
	  //
	  for(int w2 = 0; w2 < well_size(); ++w2) {
	    //
	    std::map<int, int>::const_iterator w2i = kinetic_basis[e2].well_index_map.find(w2);

	    if(w2i == kinetic_basis[e2].well_index_map.end())
	      //
	      continue;

	    dtemp = 0.;
	      
	    for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
	      //
	      dtemp += kinetic_basis[e2].eigenvector(w1i->second, k)
		//
		* kinetic_basis[e2].eigenvector(w2i->second, k)
		//
		/ kinetic_basis[e2].eigenvalue[k];

	    eigen_pop(l, w2) -= dtemp *  proj * well[w1].kernel(e1, e2) * well[w1].collision_frequency()
	      //
	      * std::sqrt(well[w1].states(e1) * thermal_factor(e1) * well[w2].states(e2) / well[w1].states(e2) / well[w2].weight());
	    //
	  }// population correction
	  //
	  // bimolecular correction
	  //
	  for(int b = 0; b < outer_barrier.size(); ++b) {
	    //
	    if(e2 >= outer_barrier[b].size())
	      //
	      continue;

	    const int& w2 = outer_barrier[b].connect.first;
		
	    const int& p  = outer_barrier[b].connect.second;
	    
	    std::map<int, int>::const_iterator w2i = kinetic_basis[e2].well_index_map.find(w2);
		  
	    if(w2i == kinetic_basis[e2].well_index_map.end())
	      //
	      continue;

	    dtemp = 0.;
	      
	    for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
	      //
	      dtemp += kinetic_basis[e2].eigenvector(w1i->second, k)
		//
		* kinetic_basis[e2].eigenvector(w2i->second, k)
		//
		/ kinetic_basis[e2].eigenvalue[k];

	    dtemp *= well[w1].kernel(e1, e2) * well[w1].collision_frequency();
	    
	    eigen_bim(l, p) -= dtemp * proj * outer_barrier[b].states(e2) / 2. / M_PI
	      //
	      * std::sqrt(well[w1].states(e1) * thermal_factor(e1) / well[w1].states(e2) / well[w2].states(e2));
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
    
  /************************************************************************
   ************************* LOW EIGENPAIR OUTPUT *************************
   ************************************************************************/

  if(global_size >= well_size()) {
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
      IO::log << std::setw(Model::log_precision + 7) << well[w].name();
  
    IO::log << "\n";

    for(int l = 0; l < well_size(); ++l) {
      //
      double nfac;
    
      for(int w = 0; w < well_size(); ++w) {
	//
	dtemp = std::fabs(eigen_pop(l, w) * well[w].weight_sqrt());

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
    
      for(int w = 0; w < well_size(); ++w) {
	//
	dtemp = eigen_pop(l, w) * well[w].weight_sqrt() / nfac;
      
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

    for(int w = 0; w < well_size(); ++w)
      //
      if(!w || dtemp < well[w].weight_sqrt())
	//
	dtemp = well[w].weight_sqrt();
  
  
    for(int w = 0; w < well_size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << well[w].weight_sqrt() / dtemp;
  
    IO::log << "\n";

    IO::log << IO::log_offset
      //
	    << "*R - eigenvalue over the relaxation limit\n"
      //
	    << IO::log_offset
      //
	    << "*P - eigenvector projection squared on the relaxation subspace (= 1 - F_ne)\n"
      //
	    << IO::log_offset
      //
	    << "*Z - well partition function square root (normalized)\n";
  }
  
  /************************************************************************************
   *************************** BIMOLECULAR RATE CONSTANTS *****************************
   *************************************************************************************/

  const int relax_size = global_size - chem_size;
  
  // bimolecular-to-bimolecular rate coefficients: slow isomerization contribution
  //
  if(bimolecular_size()) {
    //
    Lapack::Vector relax_lave(relax_size);
  
    for(int r = 0; r < relax_size; ++r) {
      //
      itemp = r + chem_size;

      relax_lave[r] = 1. / eigenval[itemp];
    }
      
    for(int p = 0; p < bimolecular_size(); ++p)
      //
      for(int q = p; q < bimolecular_size(); ++q) {
	//
	dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_bim(chem_size, q), relax_lave, relax_size) * energy_step();
	
	if(p != q) {
	  //
	  bb_rate(p, q) += dtemp;
	}
	else {
	  //
	  bb_rate(p, q) -= dtemp;
	}
      }
    
    bb_rate /= std::exp(_ener_index_begin * energy_step_over_temperature);
  
    rate_data.bb_rate = bb_rate; 
  }

  /*************************************************************************
   *************************** CHEMICAL SUBSPACE ***************************
   *************************************************************************/

  IO::log << IO::log_offset << "chemical subspace dimension = " << chem_size << "\n";

  if(chem_size) {
    //
    // chemical and thermal eigenvectors cross-products
    //
    Matrix<double> pop_chem(well_size(), chem_size);
    
    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

    bound_species_rates(pop_chem, eigenval, eigen_bim, rate_data, flags);
    //
  }//
  //
}// well-reduction method

Lapack::SymmetricMatrix MasterEquation::ReactiveComplex::fast_bb_rate () const
{
  const char funame [] = " MasterEquation::ReactiveComplex::fast_bb_rate: ";

  Lapack::SymmetricMatrix res;

  if(!bimolecular_size())
    //
    return res;

  Lapack::Vector vtemp(bimolecular_size());
  
  res.resize(bimolecular_size());

  res = 0.;
    
  std::set<int> low_outer_set;
    
  for(int e = 0; e < kinetic_basis.size(); ++e) {// energy index cycle
    //
    for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {// kinetic basis cycle
      //
      vtemp = 0.;
	
      for(int b = 0; b < outer_barrier.size(); ++b) {
	//
	const int& w = outer_barrier[b].connect.first;

	const int& p = outer_barrier[b].connect.second;

	if(e >= outer_barrier[b].size())
	  //
	  continue;

	if(low_outer_set.find(b) != low_outer_set.end())
	  //
	  continue;
	  
	std::map<int, int>::const_iterator wi = kinetic_basis[e].well_index_map.find(w);

	if(wi == kinetic_basis[e].well_index_map.end()) {
	  //
	  IO::log << IO::log_offset << "WARNING: " << outer_barrier[b].name() << " outer barrier energy[kcal/mol] = "
	    //
		  << outer_barrier[b].ground() / Phys_const::kcal << " is less than the connected " << well[w].name()
	    //
		  << " well ground energy = " << well[w].ground() / Phys_const::kcal << "\n";

	  low_outer_set.insert(b);
	    
	  continue;
	}
	  
	vtemp[p] += kinetic_basis[e].eigenvector(wi->second, k)
	  //
	  * outer_barrier[b].states(e) / 2. / M_PI / well[w].states_sqrt(e);
      }
	
      for(int p = 0; p < bimolecular_size(); ++p)
	//
	for(int q = p; q < bimolecular_size(); ++q)
	  //
	  res(p, q) += vtemp[p] * vtemp[q] * thermal_factor(e) / kinetic_basis[e].eigenvalue[k];

    }// kinetic basis cycle
    //
  }//energy index cycle

  for(int p = 0; p < bimolecular_size(); ++p)
    //
    res(p, p) = -res(p, p);
      
  for(int b = 0; b < outer_size(); ++b) {      
    //
    int p = outer_barrier[b].connect.second;

    // add capture rate
    //
    res(p, p) += outer_barrier[b].weight() / 2. / M_PI;
  }
    
  res *= energy_step();

  return res;
}

/********************************************************************************************
 ************ THE DIRECT DIAGONALIZATION OF THE GLOBAL KINETIC RELAXATION MATRIX ************
 ********************************************************************************************/

void MasterEquation::ReactiveComplex::direct_diagonalization_method (RateData& rate_data, int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::direct_diagonalization_method: ";

  IO::Marker funame_marker(funame);

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;
  Lapack::Vector      vtemp;

  IO::log << IO::log_offset << "Pressure = ";
  
  switch(pressure_unit) {
    //
  case BAR:
    //
    IO::log << pressure / Phys_const::bar << " bar";
    
    break;
    //
  case TORR:
    //
    IO::log << pressure / Phys_const::tor << " torr";
    
    break;
    //
  case ATM:
    //
    IO::log << pressure / Phys_const::atm << " atm";
    
    break;
  }

  IO::log << "\t Temperature = " << (int)std::floor(temperature / Phys_const::kelv + 0.5) << " K"
    //
	  << "\t Collision frequency = " << well[0].collision_frequency() / Phys_const::herz << " 1/sec\n";
  
  // total kinetic relaxation matrix dimension and index shifts for individual wells
  //
  well_shift.resize(well_size());
  
  itemp = 0;
  
  for(int w = 0; w < well_size(); itemp += well[w++].size())
    //
    well_shift[w] = itemp;

  global_size = itemp;

  /**************************************************************************************
   **************************** SETTING GLOBAL MATRICES *********************************
   **************************************************************************************/

  // kinetic relaxation matrix initialization
  //
  kin_mat.resize(global_size);
  
  kin_mat = 0.;

  // inner barrier isomerization contribution
  //
  for(int b = 0; b < inner_barrier.size(); ++b) {
    //
    const int& w1 = inner_barrier[b].connect.first;
      
    const int& w2 = inner_barrier[b].connect.second;
      
    for(int e = 0; e < inner_barrier[b].size(); ++e) {
      //
      if(e >= well[w1].size() || e >= well[w2].size()) {
	//
	IO::log << IO::log_offset << "WARNING: " << inner_barrier[b].name() << " inner barrier energy[kcal/mol] = "
	  //
		<< energy_bin(inner_barrier[b].size()) / Phys_const::kcal << " is lower than one of the connected wells ground energies: "
	  //
		<< well[w1].name() << " energy = " << energy_bin(well[w1].size()) / Phys_const::kcal << ", " << well[w2].name()
	  //
		<< " energy = " << energy_bin(well[w2].size()) / Phys_const::kcal << "\n";

	break;
      }
      
      kin_mat(e + well_shift[w1], e + well_shift[w2]) -= inner_barrier[b].states(e) / 2. / M_PI
	//
	/ std::sqrt(well[w1].states(e) * well[w2].states(e));
	
      kin_mat(e + well_shift[w1], e + well_shift[w1]) += inner_barrier[b].states(e) / 2. / M_PI
	//
	/ well[w1].states(e);
	
      kin_mat(e + well_shift[w2], e + well_shift[w2]) += inner_barrier[b].states(e) / 2. / M_PI
	//
	/ well[w2].states(e);
    }
  }

  // outer barrier isomerization contribution
  //
  for(int b = 0; b < outer_barrier.size(); ++b) {
    //
    const int& w = outer_barrier[b].connect.first;
      
    for(int e = 0; e < outer_barrier[b].size(); ++e) {
      //
      if(e >= well[w].size()) {
	//
	IO::log << IO::log_offset << "WARNING: " << outer_barrier[b].name() << " outer barrier energy[kcal/mol] = "
	  //
		<< energy_bin(outer_barrier[b].size()) / Phys_const::kcal << " is lower than the connected "
	  //
		<< well[w].name() << " well ground energy = " << energy_bin(well[w].size()) / Phys_const::kcal << "\n";

	break;
      }
      
      kin_mat(e + well_shift[w], e + well_shift[w]) += outer_barrier[b].states(e) / 2. / M_PI
	//
	/ well[w].states(e);
    }
  }

  // collision contribution
  //
  for(int w = 0; w < well_size(); ++w) {
    //
    for(int e = 0; e < well[w].size(); ++e) {
      //
      itemp = e + well[w].kernel_bandwidth;

      if(well[w].size() < itemp)
	//
	itemp = well[w].size();
      
      for(int f = e; f < itemp; ++f)
	//
	kin_mat(e + well_shift[w], f + well_shift[w]) +=  well[w].collision_frequency() * well[w].kernel(e, f)
	  //
	  * std::sqrt(well[w].states(e) / well[w].states(f) / thermal_factor(f - e));
    }
  }

  // bimolecular product vectors
  //
  Lapack::Matrix global_bim;
  
  if(bimolecular_size()) {
    //
    global_bim.resize(global_size,  bimolecular_size());
    
    global_bim = 0.;
    
    for(int b = 0; b < outer_barrier.size(); ++b) {
      //
      const int& w = outer_barrier[b].connect.first;
      
      const int& p = outer_barrier[b].connect.second;

      for(int e = 0; e < outer_barrier[b].size(); ++e) {
	//
	if(e >= well[w].size())
	  //
	  break;
	
	global_bim(e + well_shift[w], p) += outer_barrier[b].states(e) / 2. / M_PI
	  //
	  * std::sqrt(thermal_factor(e) / well[w].states(e));
      }
    }
  }
    
  // Boltzmann distributions
  //
  Lapack::Matrix global_pop(global_size, well_size());
  
  global_pop = 0.;

  for(int w = 0; w < well_size(); ++w) {
    //
    dtemp = well[w].weight_sqrt();
    
    for(int e = 0; e < well[w].size(); ++e)
      //
      global_pop(e + well_shift[w], w) = std::sqrt(well[w].states(e) * thermal_factor(e)) / dtemp;
  }
    
  /**************************************************************************************************
   ********************************* KINETIC MATRIX DIAGONALIZATION *********************************
   **************************************************************************************************/
    
  Lapack::Vector eigenval;
  
  Lapack::Matrix eigen_global(global_size);

  {
    IO::Marker solve_marker("diagonalizing kinetic matrix");

    IO::log << IO::log_offset << "global kinetic matrix size = " << global_size << "\n";
  
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
    //
  }// kinetic matrix diagonalization

  const double relax_eval_min = eigenval[well_size()];
  
  IO::log << IO::log_offset << "minimal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_min / well[0].collision_frequency() << "\n";
  
  IO::log << IO::log_offset << "maximal relaxation eigenvalue over collision frequency = "
    //
	  << eigenval.back() / well[0].collision_frequency() << "\n";

  Lapack::Matrix eigen_pop = eigen_global * global_pop;

  Lapack::Matrix eigen_bim;
  
  if(bimolecular_size())
    //
    eigen_bim = eigen_global * global_bim;
 
  /***********************************************************************
   ************************ LOW EIGENPAIR OUTPUT *************************
   ***********************************************************************/

  std::vector<double> relaxation_projection(well_size());
  
  for(int l = 0; l < well_size(); ++l)
    //
    relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));

  IO::log  << IO::log_offset << "eigenvector populations normalized:\n"
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
    IO::log << std::setw(Model::log_precision + 7) << well[w].name();
  
  IO::log << "\n";

  // maximal population
  //
  for(int l = 0; l < well_size(); ++l) {
    //
    double pos_pop = 0.;
    
    double neg_pop = 0.;
    
    for(int w = 0; w < well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well[w].weight_sqrt();
      
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
	    << std::setw(Model::log_precision + 7) << eigenval[l] / relax_eval_min
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * well[w].weight_sqrt() / max_pop;
      
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

  for(int w = 0; w < well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << well[w].name();
  
  IO::log << "\n";

  for(int l = 0; l < well_size(); ++l) {
    //
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(Model::log_precision + 7) << eigenval[l] / well[0].collision_frequency()
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < well_size(); ++w)
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
  
  for(int w = 0; w < well_size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << well[w].weight_sqrt();
  
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
  
  /***************************************************************************************
   *************************** CHEMICAL SUBSPACE DIMENSION *******************************
   ***************************************************************************************/

  // absolute eigenvalue threshold
  //
  if(chemical_threshold > 1.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(relax_eval_min / eigenval[itemp]  < chemical_threshold)
	//
	break;
  }
  // relaxation projection threshold
  //
  else if(chemical_threshold < 1. && chemical_threshold > 0.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(relaxation_projection[itemp] > chemical_threshold)
	//
	break;
  }
  // relative eigenvalue threshold
  //
  else if(chemical_threshold < -1. ) {
    //
    for(itemp = well_size(); itemp > 0; --itemp)
      //
      if(eigenval[itemp] / eigenval[itemp - 1] > - chemical_threshold)
	//
	break;
  }
  //
  else
    //
    itemp = well_size();
	
  const int chem_size = itemp;

  IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n";

  const int relax_size = global_size - chem_size;

  /**********************************************************************************************
   ******************************* BIMOLECULAR RATE COEFICIENTS *********************************
   **********************************************************************************************/

  if(bimolecular_size()) {
    //
    Lapack::SymmetricMatrix bb_rate(bimolecular_size());

    Lapack::Matrix kappa(well_size(), bimolecular_size());
    
    // use modified kinetic matrix approach
    //
    if(use_matrix_inversion) {
      //
      int cycle_size = global_size * global_size;
    
#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)

      // kinetic matrix modified
      //
      for(int cycle = 0; cycle < cycle_size; ++cycle) {
	//
	const int i = cycle / global_size;
      
	const int j = cycle % global_size;

	if(j < i)
	  //
	  continue;
      
	dtemp = 0.;
      
	for(int l = 0; l < chem_size; ++l)
	  //
	  dtemp += eigen_global(l, i) * eigen_global(l, j);
      
	kin_mat(i, j) += dtemp * well[0].collision_frequency();
      }

      for(int p = 0; p < bimolecular_size(); ++p)
	//
	for(int l = 0; l < chem_size; ++l)
	  //
	  parallel_orthogonalize(&global_bim(0, p), &eigen_global(l, 0), global_size, 1, global_size);
    
      Lapack::Matrix  global_bim_over_km = Lapack::Cholesky(kin_mat).invert(global_bim);

      for(int w = 0; w < well_size(); ++w)
	//
	for(int l = 0; l < chem_size; ++l)
	  //
	  parallel_orthogonalize(&global_pop(0, w), &eigen_global(l, 0), global_size, 1, global_size);
  
      // kappa matrix
      //
      kappa = global_pop.transpose() * global_bim_over_km;

      // bimolecular-to-bimolecular rate coefficients
      //
      bb_rate = Lapack::SymmetricMatrix(global_bim.transpose() * global_bim_over_km);
    }
    // using relaxational eigenvectors
    //
    else {
      //
      Lapack::Vector relax_lave(relax_size);
  
      for(int r = 0; r < relax_size; ++r) {
	//
	itemp = r + chem_size;

	relax_lave[r] = 1. / eigenval[itemp];
      }

      // kappa matrix
      //
      for(int w = 0; w < well_size(); ++w)
	//
	for(int p = 0; p < bimolecular_size(); ++p)
	  //
	  kappa(w, p) = triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w), relax_lave, relax_size);

      // bimolecular-to-bimolecular rate coefficients
      //
      for(int p = 0; p < bimolecular_size(); ++p)
	//
	for(int q = p; q < bimolecular_size(); ++q)
	  //
	  bb_rate(p, q) = triple_product(&eigen_bim(chem_size, p), &eigen_bim(chem_size, q), relax_lave, relax_size);
    }
    
    // kappa matrix output
    //
    IO::log << IO::log_offset << "isomers-to-bimolecular equilibrium coefficients (kappa matrix):\n"
      //
	    << IO::log_offset << std::setw(5) << "W\\P";
    
    for(int p = 0; p < bimolecular_size(); ++p)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(p).short_name();
    
    IO::log << "\n";

    for(int w = 0; w < well_size(); ++w) {
      //
      IO::log << IO::log_offset << std::setw(5) << well[w].name();
      
      for(int p = 0; p < bimolecular_size(); ++p) {
	//
	dtemp = kappa(w, p) / well[w].weight_sqrt();
	
	IO::log << std::setw(Model::log_precision + 7);
	
	if(dtemp < 0.1 && dtemp > -0.1) {
	  //
	  IO::log << "<.1";
	}
	else
	  //
	  IO::log << dtemp;
      }
      
      IO::log << "\n";
    }
    
    for(int p = 0; p < bimolecular_size(); ++p)
      //
      bb_rate(p, p) = -bb_rate(p, p);
    
    for(int b = 0; b < outer_size(); ++b) {
      //
      const int& p = outer_barrier[b].connect.second;

      bb_rate(p, p) += outer_barrier[b].weight() / 2. / M_PI;
    }

    bb_rate *= energy_step();

    rate_data.bb_rate = bb_rate;
  }

  /**************************************************************
   *********************** BOUND SPECIES ************************
   **************************************************************/

  if(chem_size) {
    //
    // projection of the chemical eigenvectors on the thermal subspace
    //
    Matrix<double> pop_chem(well_size(), chem_size);

    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

    bound_species_rates(pop_chem, eigenval, eigen_bim, rate_data, flags);
  }
  //
}// direct diagonalization method

// well-to-well, well-to-bimolecular, bimolecular-to-well rate coefficients and wells partition
//
void MasterEquation::ReactiveComplex::bound_species_rates (const Matrix<double>& pop_chem, Lapack::Vector eigenval, Lapack::Matrix eigen_bim,
							   //
							   RateData& rate_data, int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::bound_species_rates: ";

  double dtemp;
  
  int    itemp;
  
  const int chem_size = pop_chem.size2();

  partition_t wp;
  
  group_t     bg;

  std::vector<double> weight(well_size());

  for(int w = 0; w < well_size(); ++w)
    //
    weight[w] = well[w].weight();
  
  // partitioning wells into equilibrated groups
  //
  WellSet<double> ws(weight);

  itemp = 0;

  if(!bimolecular_size())
    //
    itemp = WellSet<double>::INTERNAL;
  
  dtemp = ws.threshold_well_partition(pop_chem, wp, bg, itemp);

  IO::log << IO::log_offset << "well partition:        " << name(wp) << "\n";

  if(bg.size())
    //
    IO::log << IO::log_offset << "bimolecular group:     " << name(bg) << "\n";
  
  IO::log << IO::log_offset << "well partition error = " << dtemp << "\n";

  Lapack::Matrix m_direct;
    
  if(chem_size != well_size()) {
    //
    // convert chemical eigenvectors in the new basis
    //
    m_direct = ws.basis(wp).transpose() * (Lapack::Matrix)pop_chem;

    weight = ws.weight(wp);
  }
  else
    //
    m_direct = (Lapack::Matrix)pop_chem;

  Lapack::Matrix m_inverse = m_direct.invert();
  
  // well-to-well rate coefficients
  //
  itemp = 0;
  
  if(!bimolecular_size())
    //
    itemp = 1;
  
  Lapack::Matrix ww_rate(chem_size);
    
  ww_rate = 0.;
    
  for(int i = 0; i < chem_size; ++i) {
    //
    for(int j = 0; j < chem_size; ++j) {
      //
      for(int l = itemp; l < chem_size; ++l) {
	//
	if(i != j) {
	  //
	  ww_rate(i, j) -= m_direct(j, l) * m_inverse(l, i) * eigenval[l];
	}
	else {
	  //
	  ww_rate(i, j) += m_direct(j, l) * m_inverse(l, i) * eigenval[l];
	}
      }
      
      ww_rate(i, j) *= std::sqrt(weight[i] * weight[j]) * energy_step();
    }
  }

  ww_rate /= std::exp(_ener_index_begin * energy_step_over_temperature);
  
  rate_data.ww_rate = ww_rate;
    
  // well-to-bimolecular and bimolecular-to-well rate coefficients
  //
  if(bimolecular_size()) {
    //
    Lapack::Matrix wb_rate(chem_size, bimolecular_size());
      
    Lapack::Matrix bw_rate(bimolecular_size(), chem_size);
      
    for(int w = 0; w < chem_size; ++w) {
      //
      dtemp = std::sqrt(weight[w]) * energy_step();
	
      for(int p = 0; p < bimolecular_size(); ++p) {
	//
	wb_rate(w, p) = vdot(m_inverse.column(w), &eigen_bim(0, p)) * dtemp;

	bw_rate(p, w) = vdot(m_direct.row(w), &eigen_bim(0, p))     * dtemp;
      }
    }
      
    wb_rate /= std::exp(_ener_index_begin * energy_step_over_temperature);
  
    bw_rate /= std::exp(_ener_index_begin * energy_step_over_temperature);
  
    rate_data.wb_rate = wb_rate;
      
    rate_data.bw_rate = bw_rate;
  }

  rate_data.well_partition    = index_to_well(wp);
}

/***********************************************************************************************
 ************************************** GENERATORS *********************************************
 ***********************************************************************************************/

/************** Generator of partitions of n items into m groups ****************/

MasterEquation::PartitionGenerator::PartitionGenerator (int m, int n)
  //
  : _group_index(n), _frame(m + 1), _end(false), _partition_size(m)
{
  const char funame [] = "MasterEquation::PartitionGenerator::PartitionGenerator: ";
  
  if(m > n) {
    //
    IO::log << IO::log_offset << "number of groups is bigger than the number of items\n";
    
    throw Error::Range();
  }

  _frame[m] = n;
  
  for(int i = 0; i < m; ++i)
    //
    _group_index[i]= _frame[i] = i;
}

void MasterEquation::PartitionGenerator::operator++ ()
{
  if(_end)
    //
    return;

  // iterate within a frame
  //
  for(int i = _frame[1] + 1, j = 1; i < _group_index.size(); ++i)
    //
    if(i == _frame[j + 1]) {
      //
      ++j;
    }
    else if(_group_index[i] == j) {
      //
      _group_index[i] = 0;
    }
    else {
      //
      ++_group_index[i];
      
      return;
    }
    
  // iterate over different frames
  //
  for(int j = 1; j < _partition_size; ++j)
    //
    if(_frame[j] + 1 == _frame[j + 1]) {
      //
      _frame[j] = _frame[j - 1] + 1;
    }
    else {
      //
      ++_frame[j];
      
      // initialize the index array
      //
      for(int i = 1, k = 1; i < _group_index.size(); ++i)
	//
	if(i == _frame[k]) {
	  //
	  _group_index[i] = k;
	  
	  ++k;
	}
	else
	  //
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
    IO::log << IO::log_offset << "map dimension mismatch\n";
    
    throw Error::Range();
  }

  return res;
}

/*************** Generator of group of m elements from the pool of n elements *******************/

MasterEquation::GroupGenerator::GroupGenerator (int m, int n)
  //
  : _index(m + 1), _end(false), _size(m)
{
  const char funame [] = "MasterEquation::GroupGenerator::GroupGenerator: ";
  
  if(m > n) {
    //
    IO::log << IO::log_offset << "group size is bigger than the pool size\n";
    
    throw Error::Range();
  }

  for(int i = 0; i < m; ++i)
    //
    _index[i] = i;
  
  _index[m] = n;
}

void MasterEquation::GroupGenerator::operator++ ()
{
  bool btemp;

  if(_end)
    //
    return;

  btemp = false;
  
  for(int i = 0; i < size(); ++i)
    //
    if(_index[i] + 1 == _index[i + 1]) {
      //
      if(i) {
	//
	_index[i] = _index[i - 1] + 1;
      }
      else {
	//
	_index[i] = 0;
      }
    }
    else {
      //
      ++_index[i];
      
      btemp = true;
      
      break;
    }

  if(btemp)
    //
    return;
  
  _end = true;
}

