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

#include "crossrate.hh"
#include "random.hh"
#include "read.hh"
#include "units.hh"
#include "key.hh"

#include <iomanip>
#include <cmath>
#include <set>

namespace CrossRate {

  bool _isinit = false;
  bool isinit () { return _isinit; }

  // random potential error flag; opposite assumes that potential error corresponds to large positive potential
  int rand_pot_err_flag = 0;

  // output
  int         raden_flag;
  std::string raden_file;
  int         amproj_flag;
  std::string amproj_file;

 // minimum number of surface samplings to find facets
  int    MultiArray::min_sur_size;
 // minimal importance samplings array size for the facet
  int    MultiArray::min_imp_size;
 // maximum importance samplings array size for the facet
  int    MultiArray::max_imp_size;
  // minimum number of facet samplings before the accuracy can be estimated
  int    MultiArray::min_pot_size;
  // maximum number of facet samplings
  int    MultiArray::max_pot_size;

  // tolerances
  double MultiArray::face_rel_tol;  // facet   statistical flux relative tolerance
  double MultiArray::spec_rel_tol;  // species statistical flux relative tolerance
  double MultiArray::reac_rel_tol;  // species    reactive flux relative tolerance
  double MultiArray::tran_rel_tol;  // reactive transition flux relative tolerance

  // reactive transition
  std::vector<int> MultiArray::reactive_transition;

  // output stream
  std::ofstream xout;

  job_t _job;
  job_t job () { return _job; }
   
  mode_t _mode;
  mode_t mode () { return _mode; }

  double _temperature;
  double temperature () { return _temperature; }

  double _tot_ener;
  double energy_value () { return _tot_ener; }

  double _tot_amom;
  double amom_value () { return _tot_amom; }

  double _norm_factor;
  double norm_factor () { return _norm_factor; }

  int _reactant;
  int reactant () { return _reactant; }

  double _reac_ener = 1.;
  double   reactive_energy ()         { return _reac_ener; }
  void set_reactive_energy (double e) { _reac_ener = e; }
  
  int _seed;
  int seed () { return _seed; }

  std::vector<double> traj_rel_tol;
  std::vector<double> traj_abs_tol;
}

void CrossRate::init (std::istream& from) 
{
    const char funame [] = "CrossRate::init: ";

    if(isinit()) {
      std::cerr << funame << "has been initialized already\n";
      throw Error::Init();
    }
    _isinit = true;

    if(!Structure::isinit()) {
      std::cerr << funame << "structure has not been initialized\n";
      throw Error::Init();
    }

    IO::Marker funame_marker(funame);

    int itemp;
    double dtemp;
    bool btemp;

    std::map<std::string, Read> input;
    std::map<std::string, Read>::iterator idit;

    double trt;
    std::string calc_mode, job_type;
    input ["JobType"                    ] = Read(job_type, "dynamical");
    input ["CalculationMode"            ] = Read(calc_mode, "canonical");
    input ["Reactant"                   ] = Read(_reactant, -1);
    input ["Temperature"                ] = Read(_temperature, -1);
    input ["Energy"                     ] = Read(_tot_ener, -100);
    input ["AngularMomentum"            ] = Read(_tot_amom, -1);
    input ["MinSurfaceSamplingNumber"   ] = Read(MultiArray::min_sur_size, 10000);
    input ["MinImportanceSamplingNumber"] = Read(MultiArray::min_imp_size, 100);
    input ["MaxImportanceSamplingNumber"] = Read(MultiArray::max_imp_size, 10000);
    input ["MinPotentialSamplingNumber" ] = Read(MultiArray::min_pot_size, 100);
    input ["MaxPotentialSamplingNumber" ] = Read(MultiArray::max_pot_size, 1000000);
    input ["FacetStatisticalTolerance"  ] = Read(MultiArray::face_rel_tol, 0.1);
    input ["SpeciesStatisticalTolerance"] = Read(MultiArray::spec_rel_tol, 0.1);
    input ["SpeciesReactiveTolerance"   ] = Read(MultiArray::reac_rel_tol, 0.1);
    input ["ReactiveTransitionTolerance"] = Read(MultiArray::tran_rel_tol, 0.1);
    input ["ReactiveTransition"         ] = Read(MultiArray::reactive_transition, std::vector<int>());
    input ["TrajectRelativeTolerance"   ] = Read(trt, 1.e-5);
    input ["RandomPotentialErrorFlag"   ] = Read(rand_pot_err_flag, 0);
    input ["RadialEnergyFlag"           ] = Read(raden_flag, 0);
    input ["RadialEnergyFile"           ] = Read(raden_file, "raden.out");
    input ["AngularMomentumProjectFlag" ] = Read(amproj_flag, 0);
    input ["AngularMomentumProjectFile" ] = Read(amproj_file, "amproj.out");
    input ["RandomSeed"                 ] = Read(_seed, 0);
    input ["MinAtomDistance"            ] = Read(Dynamic::Coordinates::min_atom_dist, 1.5);

    std::string key, comment;
    while(from >> key) {
	if(key == IO::end_key())
	    break;
	idit = input.find(key);
	if (idit == input.end()) {
	  std::cerr << funame << "WARNING: did not find the key " << key
		    << ", taking it as a comment\n";
	  std::getline(from, comment);
	}
	else
	    from >> idit->second;
    }

    if(!from) {
	std::cerr << funame << "input stream is corrupted\n";
	throw Error::Form();
    }

    // check if all parameters were initialized
    btemp = true;
    for(idit = input.begin(); idit != input.end(); ++idit)
	if(!idit->second.is_init()) {
	    std::cerr << funame << idit->first << " is not initialized\n";
	    btemp = false;
	}
    if(!btemp)
	throw Error::Init();

    // default values
    IO::log << IO::log_offset << "Default parameters:\n" << std::left;
    for(idit = input.begin(); idit != input.end(); ++idit)
      if(idit->second.is_default())
	IO::log << IO::log_offset << "   " << std::setw(20) << idit->first << " = " << idit->second << "\n";
    IO::log << std::right;

    if(job_type == "dynamical")
	_job = DYN_JOB;
    else if(job_type == "statistical")
	_job = STAT_JOB;
    else if(job_type == "test")
	_job = TEST_JOB;
    else {
      std::cerr << funame << "unknown job type: " << job_type
		<< "; possible jobs: dynamical, statistical, and test\n";
      throw Error::Init();
    }
    
    if(calc_mode == "canonical")
      _mode = T_MODE;
    else if(calc_mode == "microcanonical")
      _mode = E_MODE;
    else if(calc_mode == "e,j-resolved")
      _mode = J_MODE;
    else {
      std::cerr << funame << "unknown caculation mode: " << calc_mode
		<< "; possible modes: canonical, microcanonical, and e,j-resolved\n";
      throw Error::Init();
    }

    switch(mode()) {
    case T_MODE:
      if(temperature() <= 0.) {
	std::cerr << funame << "Temperature, " << temperature() / Phys_const::kelv << " K, should be positive\n";
	throw Error::Init();
      }
      break;
    case E_MODE:
      if(energy_value() < -1.) {
	std::cerr << funame << "Energy, " << energy_value() / Phys_const::kcal << " kcal/mol, is out of range\n";
	throw Error::Init();
      }
      break;
    case J_MODE:
      if(energy_value() <= -1.) {
	std::cerr << funame << "Energy, " << energy_value() / Phys_const::kcal << " kcal/mol, is out of range\n";
	throw Error::Init();
      }
      if(amom_value() < 0.) {
	std::cerr << funame << "Angular momentum, " << amom_value() << ", should not be negative\n";
	throw Error::Init();
      }
      break;
    }

    // trajectories tolerances
    traj_rel_tol.resize(Structure::dv_size(), trt);
    traj_abs_tol.resize(Structure::dv_size(), trt);

    double scale;
    if(mode() == T_MODE) // using the temperature as a scaling factor
      scale = std::sqrt(temperature());
    else if(energy_value() < 0.)// using total energy as a scaling factor
      scale = std::sqrt(-energy_value());
    else
      scale = std::sqrt(energy_value());

    for(int i = 0; i < 3; ++i) {
      itemp = Structure::pos_size() + Structure::orb_vel() + i;
      traj_abs_tol[itemp] *= scale/Structure::mass_sqrt();
    }
    for(int frag = 0; frag < 2; ++frag)
      for(int i = 0; i < Structure::fragment(frag).vel_size(); ++i) {
	itemp = Structure::pos_size() + Structure::ang_vel(frag) + i;
	traj_abs_tol[itemp] *= scale/Structure::fragment(frag).imom_sqrt(i);
      }

    // normalization factor
    _norm_factor = std::pow(Structure::mass(), 1.5) * 2. / M_PI / M_PI;
    for(int frag = 0; frag < 2; ++frag)
      switch (Structure::fragment(frag).type()) {
      case Molecule::MONOATOMIC:
	_norm_factor *= M_PI;
	break;
      case Molecule::LINEAR:
	_norm_factor *= 2. * M_PI * Structure::fragment(frag).imom(0);
      break;
      case Molecule::NONLINEAR:
	for (int i = 0; i < 3; ++i)
	  _norm_factor *= std::sqrt(2. * M_PI * Structure::fragment(frag).imom(i));
	break;
    }

    switch(mode()) {
    case T_MODE:
      _norm_factor *= std::pow(temperature(), double(Structure::tm_dof() - 1) / 2.);
      break;
    case E_MODE:
      _norm_factor /= tgamma(double(Structure::tm_dof() + 1) / 2.);
      break;
    case J_MODE:
      //_norm_factor *= amom_value() * amom_value() * std::sqrt(2. / M_PI) 
      // / tgamma(double(Structure::tm_dof() - 2) / 2.);
      _norm_factor *= std::sqrt(2. / M_PI) / tgamma(double(Structure::tm_dof() - 2) / 2.);
      break;
    }
}

/**************************************************************************
 *                           Sampling methods                             *
 **************************************************************************/


CrossRate::DynSmp::DynSmp (Potential::Wrap pot, const DivSur::MultiSur& surface, int prim, const Dynamic::Coordinates& dc)
  : Dynamic::Vars(dc)
{
  const char funame [] = "CrossRate::DynSmp::DynSmp: ";
  
  static const double max_exp = 100.;

  double dtemp;

  // potential energy
  _energy = pot(dc);
  
  // random variable
  _ranval = Random::flat();

  // statistical weight and random velocities
  switch(mode()) {
  case T_MODE: // canonical calculation
    _weight = surface.t_velocity(prim, temperature(), *this);

    dtemp = - potential_energy() / temperature();
    if(dtemp > max_exp) {
      std::cerr << funame << "potential energy, " << potential_energy() / Phys_const::kcal << " kcal/mol, is too low\n";
      throw Error::Range();
    }
    else if(dtemp > - max_exp)
      _weight *= std::exp(dtemp);
    else
      _weight = -1.;

    break;
  case E_MODE: // microcanonical calculation
    if(potential_energy() < energy_value())
      _weight = surface.e_velocity(prim, energy_value() - potential_energy(), *this);
    else
      _weight = -1.;
    break;
  case J_MODE: // E,J-resolved calculation
    if(potential_energy() < energy_value())
      _weight = surface.j_velocity(prim, energy_value() - potential_energy(), amom_value(), *this);
    else
      _weight = -1.;
    break;
  }

  back.init(new DynRes(pot, *this, BACKWARD));
  forw.init(new DynRes(pot, *this,  FORWARD));
}

// propagate trajectory both forward and backward
void CrossRate::DynSmp::run_traj (const DivSur::MultiSur& ms, const DivSur::face_t& face, Dynamic::CCP stop)
{
  const char funame [] = "CrossRate::DynSmp::run_traj: ";

  if(!isinit()) {
    std::cerr << funame << "crossrate environment has not yet been initialized\n";
    throw Error::Init();
  }

  if(is_run()) {
    std::cerr << funame << "trajectory has been run already\n";
    throw Error::Init();
  }

  if(reactant() >= 0 && face.first != reactant() && face.second != reactant())
    return;

  int    itemp;
  double dtemp;

  back->rel_tol = traj_rel_tol;
  back->abs_tol = traj_abs_tol;
  forw->rel_tol = traj_rel_tol;
  forw->abs_tol = traj_abs_tol;
  if(temperature() < 0.) {
    dtemp = std::sqrt(energy_value() - potential_energy());
    for(int i = 0; i < 3; ++i) {
      itemp = Structure::pos_size() + Structure::orb_vel() + i;
      back->abs_tol[itemp] *= dtemp / Structure::mass_sqrt();
      forw->abs_tol[itemp] *= dtemp / Structure::mass_sqrt();
    }
  
    for(int frag = 0; frag < 2; ++frag)
      for(int i = 0; i < Structure::fragment(frag).vel_size(); ++i) {
	itemp = Structure::pos_size() + Structure::ang_vel(frag) + i;
	back->abs_tol[itemp] *= dtemp / Structure::fragment(frag).imom_sqrt(i);
	forw->abs_tol[itemp] *= dtemp / Structure::fragment(frag).imom_sqrt(i);
      }
  }

  if(reactant() >= 0) {
    if(reactant() == face.first){
      try {
	back->run(stop | Dynamic::negate(Dynamic::CCP(new SpecCondition(ms, face.first))), ms);

	if(back->species() == face.first)
	  back->stat = DynRes::DIRECT;
	else {
	  back->stat = DynRes::RECROSS;
	  return;
	}
      }
      catch(Trajectory::PotentialFailure) {
	back->stat =  DynRes::POT_FAIL;
	return;
      }
      catch(Trajectory::RunFailure) {
	back->stat =  DynRes::RUN_FAIL;
	return;
      }
      catch(Trajectory::ExcludeRegionHit) {
	back->stat =  DynRes::EXCLUDE;
	return;
      }

      try {
	forw->run(stop, ms);
	forw->stat = DynRes::PASS;
      }
      catch(Trajectory::PotentialFailure) {
	forw->stat =  DynRes::POT_FAIL;
	return;
      }
      catch(Trajectory::RunFailure) {
	forw->stat =  DynRes::RUN_FAIL;
	return;
      }
      catch(Trajectory::ExcludeRegionHit) {
	forw->stat =  DynRes::EXCLUDE;
	return;
      }
    }
    else if(reactant() == face.second){
      try {
	forw->run(stop | Dynamic::negate(Dynamic::CCP(new SpecCondition(ms, face.second))), ms);

	if(forw->species() == face.second)
	  forw->stat = DynRes::DIRECT;
	else {
	  forw->stat = DynRes::RECROSS;
	  return;
	}
      }
      catch(Trajectory::PotentialFailure) {
	forw->stat =  DynRes::POT_FAIL;
	return;
      }
      catch(Trajectory::RunFailure) {
	forw->stat =  DynRes::RUN_FAIL;
	return;
      }
      catch(Trajectory::ExcludeRegionHit) {
	forw->stat =  DynRes::EXCLUDE;
	return;
      }

      try {
	back->run(stop, ms);
	back->stat = DynRes::PASS;
      }
      catch(Trajectory::PotentialFailure) {
	back->stat =  DynRes::POT_FAIL;
	return;
      }
      catch(Trajectory::RunFailure) {
	back->stat =  DynRes::RUN_FAIL;
	return;
      }
      catch(Trajectory::ExcludeRegionHit) {
	back->stat =  DynRes::EXCLUDE;
	return;
      }
    }
    return;
  }// reactant

  try {
    back->run(stop | Dynamic::negate(Dynamic::CCP(new SpecCondition(ms, face.first))), ms);

    if(back->species() == face.first)
      back->stat = DynRes::DIRECT;
    else
      back->stat = DynRes::RECROSS;
  }
  catch(Trajectory::PotentialFailure) {
    back->stat =  DynRes::POT_FAIL;
    return;
  }
  catch(Trajectory::RunFailure) {
    back->stat =  DynRes::RUN_FAIL;
    return;
  }
  catch(Trajectory::ExcludeRegionHit) {
    back->stat =  DynRes::EXCLUDE;
    return;
  }

  try {
    forw->run(stop | Dynamic::negate(Dynamic::CCP(new SpecCondition(ms, face.second))), ms);

    if(forw->species() == face.second)
      forw->stat = DynRes::DIRECT;
    else
      forw->stat = DynRes::RECROSS;
  }
  catch(Trajectory::PotentialFailure) {
    forw->stat =  DynRes::POT_FAIL;
    return;
  }
  catch(Trajectory::RunFailure) {
    forw->stat =  DynRes::RUN_FAIL;
    return;
  }
  catch(Trajectory::ExcludeRegionHit) {
    forw->stat =  DynRes::EXCLUDE;
    return;
  }

  // either forward and backward trajectories both recrossed or both finished
  if(forw->stat == back->stat)
    return;

  // run backward to the end
  if(back->stat == DynRes::RECROSS)
    try {
      back->run(stop, ms);
    }
    catch(Trajectory::PotentialFailure) {
      back->stat =  DynRes::POT_FAIL;
      return;
    }
    catch(Trajectory::RunFailure) {
      back->stat =  DynRes::RUN_FAIL;
      return;
    }
    catch(Trajectory::ExcludeRegionHit) {
      back->stat =  DynRes::EXCLUDE;
      return;
    }

  // run forward to the end
  if(forw->stat == DynRes::RECROSS)
    try {
      forw->run(stop, ms);
    }
    catch(Trajectory::PotentialFailure) {
      forw->stat =  DynRes::POT_FAIL;
      return;
    }
    catch(Trajectory::RunFailure) {
      forw->stat =  DynRes::RUN_FAIL;
      return;
    }
    catch(Trajectory::ExcludeRegionHit) {
      forw->stat =  DynRes::EXCLUDE;
      return;
    }
}

/************************************************************************
 *               Facet sampling data methods: FacetArray                *
 ************************************************************************/

int CrossRate::FacetArray::init_traj_num () const
{
  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit)
    if(!fit->is_run())
      ++res;
  return res;
}

int CrossRate::FacetArray::potfail_traj_num () const
{
  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit)
    if(fit->is_pot_fail())
      ++res;
  return res;
}

int CrossRate::FacetArray::runfail_traj_num () const
{
  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit)
    if(fit->is_run_fail())
      ++res;
  return res;
}

int CrossRate::FacetArray::exclude_traj_num () const
{
  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit)
    if(fit->is_exclude())
      ++res;
  return res;
}

int CrossRate::FacetArray::run_traj_num () const
{
  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit) {
    if(!fit->is_run() || fit->is_pot_fail() || fit->is_run_fail() || fit->is_exclude())
      continue;
    ++res;
  }
  return res;
}

int CrossRate::FacetArray::recross_traj_num (int dir) const
{
  const char funame [] = "CrossRate::FacetArray::recross_traj_num: ";

  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit) {
    if(!fit->is_run() || fit->is_pot_fail() || fit->is_run_fail() || fit->is_exclude())
      continue;

    switch(dir) {
    case BACKWARD:
      if(fit->forw->stat == DynRes::RECROSS)
	++res;
      break;
    case FORWARD:
      if(fit->back->stat == DynRes::RECROSS)
	++res;
      break;
    default:
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }
  }
  return res;
}

int CrossRate::FacetArray::reac_traj_num (int dir, int spec) const
{
  const char funame [] = "CrossRate::FacetArray::reac_traj_num: ";

  int res = 0;
  for(const_iterator fit = begin(); fit != end(); ++fit) {
    if(!fit->is_run() || fit->is_pot_fail() || fit->is_run_fail() || fit->is_exclude())
      continue;

    switch(dir) {
    case BACKWARD:
      if(fit->forw->stat == DynRes::DIRECT && fit->back->species() == spec)
	++res;
      break;
    case FORWARD:
      if(fit->back->stat == DynRes::DIRECT && fit->forw->species() == spec)
	++res;
      break;
    default:
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }
  }
  return res;
}

// propagate trajectories for facet samplings
void CrossRate::FacetArray::run_traj (const DivSur::MultiSur& ms, const DivSur::face_t& face, Dynamic::CCP stop)
{
  const char funame [] = "CrossRate::FacetArray::run_traj: ";
  
  if(reactant() >= 0 && face.first != reactant() && face.second != reactant())
    return;

  if(!size())
    return;

  int new_traj_num = init_traj_num();
  if(!new_traj_num)
    return;

#ifndef DEBUG

  int new_share, old_share = 0;

#endif

  int count = 0;
  for(iterator fit = begin(); fit != end(); ++fit) // sampling cycle
    if(!fit->is_run()) {
      fit->run_traj(ms, face, stop);
      ++count;

#ifdef DEBUG

      if(!fit->is_run_fail() &&  !fit->is_pot_fail() &&  !fit->is_exclude()) {
	IO::log << count << "-th trajectory:\n";
	IO::log << std::setw(13) << "time/a.u."
		  << std::setw(13) << "dist/bohr"
		  << std::setw(13) << "ener/kcal"
		  << std::setw(13) << "amom/a.u." << "\n";

	IO::log << std::setw(13) << "0"
		  << std::setw(13) << fit->interfragment_distance()
		  << std::setw(13) << fit->total_energy() / Phys_const::kcal;

	D3::Vector tam;
	fit->total_angular_momentum(tam);
	for(int i = 0; i < 3; ++i)
	  IO::log << std::setw(13) << tam[i];
	IO::log << "\n";

	if(fit->forw->stat != DynRes::INIT) {
	  IO::log << std::setw(13) << fit->forw->time()
		    << std::setw(13) << fit->forw->interfragment_distance() 
		    << std::setw(13) << fit->forw->total_energy() / Phys_const::kcal;

	  fit->forw->total_angular_momentum(tam);
	  for(int i = 0; i < 3; ++i)
	    IO::log << std::setw(13) << tam[i];
	  IO::log << "\n";
	}

	if(fit->back->stat != DynRes::INIT) {
	  IO::log << std::setw(13) << fit->back->time()
		    << std::setw(13) << fit->back->interfragment_distance() 
		    << std::setw(13) << fit->back->total_energy() / Phys_const::kcal;

	  fit->back->total_angular_momentum(tam);
	  for(int i = 0; i < 3; ++i)
	    IO::log << std::setw(13) << tam[i];
	  IO::log << "\n";
	}
	IO::log << "\n";
      }
      else {
	IO::log << "   failed\n\n";
      }

#else

      new_share =(int)((double)count / (double)new_traj_num * 100.);
      print_progress(old_share, new_share);

#endif

    }// sampling cycle

  //IO::log << "\n";

  // checking
  for(const_iterator fit = begin(); fit != end(); ++fit) { // sampling cycle
    if(fit->is_run_fail() || fit->is_pot_fail() || fit->is_exclude() || !fit->is_run())
      continue;

    if(fit->back->stat == DynRes::DIRECT && fit->back->species() != face.first) {
      std::cerr << funame << "WARNING: backward direct trajectory finished in the wrong well\n";
    }
	    
    if(fit->forw->stat == DynRes::DIRECT && fit->forw->species() != face.second) {
      std::cerr << funame << "WARNING: forward direct trajectory finished in the wrong well\n";
    }    
  }// sampling cycle
}

bool CrossRate::FacetArray::add_smp (Potential::Wrap pot, const DivSur::MultiSur& surface, int prim, const Dynamic::Coordinates& dc)
{
  const char funame [] = "CrossRate::FacetArray::add_smp: ";

  try {
    DynSmp smp(pot, surface, prim, dc);

    // update minimal energy
    if(!flux_num() || smp.potential_energy() < _min_ener) {
      _min_ener = smp.potential_energy();
      _min_geom = dc;
    }

    ++_flux_num;

    // zero weight sampling
    if(smp.weight() <= 0.)
      return true;

    // update flux
    _flux += smp.weight();
    _fvar += smp.weight() * smp.weight();

    // update importance samplings data
    if(!size()) {
      _max_weight = smp.weight();
      push_back(smp);
      return true;
    }

    if(smp.weight() > _max_weight) {
      // use new sampling weight as a reference weight
      _max_weight = smp.weight();

      // remove all non-quallifying samplings from the
      // importance samplings list
      iterator it = begin();
      while(it != end())
	if(it->ranval() * _max_weight >= it->weight())
	  it = erase(it);
	else
	  ++it;
    
      // add sampling to the importance samplings list
      push_back(smp);
    }
    else if(smp.weight() > smp.ranval() * _max_weight)
      push_back(smp);

    return true;
  }
  catch(Error::General) {
    //std::cerr << funame << "potential energy calculation failed\n";
    ++_fail_num;
    return false;
  }
}
  
inline int CrossRate::FacetArray::samp_num () const 
{ 
  if(rand_pot_err_flag)
    return  _flux_num;
  else
    return  _flux_num + _fail_num;
}

double CrossRate::FacetArray::flux_val () const
{
  if(_flux == 0.) 
    return 0.;
  
  return _flux / double(samp_num());
}

double CrossRate::FacetArray::flux_var () const
{
  const char funame [] = "CrossRate::FacetArray::flux_var: ";

  if(_flux == 0.)
    return 0.;

  double dtemp = flux_val();
  double res   = _fvar / double(samp_num()) - dtemp * dtemp;

  if(res < 0.) {
    std::cerr << funame << "WARNING: negative variance\n";
    return 0.;
  }

  return res;
}
  
double CrossRate::FacetArray::flux_rel_var () const
{
  const char funame [] = "CrossRate::FacetArray::flux_rel_var: ";

  if(_flux == 0.)
    return 0.;
  
  double res = _fvar / _flux / _flux - 1. / (double)(samp_num());

  if(res <= 0.) {
    IO::log << funame << "WARNING: negative variance\n";
    return 0.;
  }

  return res;
}

/****************************************************************************************
 *************************************** MULTIARRAY *************************************
 ****************************************************************************************/

std::ostream& operator<< (std::ostream& to, const DivSur::face_t& face) 
{
  std::ostringstream oss;
  
  oss << "(" << face.first << " -> " << face.second << ")";
  
  to << oss.str();
    
  return to;
}


namespace IO {
  //
  template<>
  LogOut& operator<< (LogOut& to, const DivSur::face_t& face) 
  {
    std::ostringstream oss;
  
    oss << "(" << face.first << " -> " << face.second << ")";

    to << oss.str();

    return to;
  }
}


// Satisfy requirements for the minimal potential energy facet samplings number,  
// for the minimal importance facet samplings  number,
// for the uniform relative facet flux error,  
// and for the uniform facet volume fraction relative error 

void CrossRate::MultiArray::run_traj (Dynamic::CCP stop) 
{
  // check that minimal potential enery is bigger than reactive energy
  for(const_iterator mit = begin(); mit != end(); ++mit)
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit)
      if(sit->second.size() && sit->second.min_ener() < reactive_energy()) {
	std::cerr << mit - begin() << "-th surface,  " << sit->first << " facet: " 
		  << "minimal energy[kcal/mol], "<< sit->second.min_ener() /Phys_const::kcal
		  << ", is smaller than the reactive energy: "
		  << "decrease the reactive energy or change the transition state dividing surface\n";
	throw Error::Run();
      }

  for(iterator mit = begin(); mit != end(); ++mit) {// primitives cycle
    IO::log << "   " << mit - begin() << "-th surface:\n";
    for(SurArray::iterator sit = mit->begin(); sit != mit->end(); ++sit)
      if(sit->second.init_traj_num()) {// facet cycle
	IO::log <<  "      " << sit->first << " facet:\n";
	sit->second.run_traj(_ms, sit->first, stop);
      }// facet cycle
  }// primitives cycle
}


int CrossRate::MultiArray::_sample (iterator mit, const std::set<DivSur::face_t >& face_work)
{
  const char funame [] = "CrossRate::MultiArray::_sample: ";

  static int  imp_warn = 1;
  static int samp_warn = 1;

  Dynamic::Coordinates new_conf;
  SurArray::iterator sit;

  const int sur = mit - begin();

  _ms.random_orient(sur, new_conf);
  DivSur::MultiSur::SmpRes smp_res = _ms.facet_test(sur, new_conf);

  switch(smp_res.stat) {
  case DivSur::MultiSur::SmpRes::CLOSE:
    mit->add_close();
    return 0;
  case DivSur::MultiSur::SmpRes::EXCLUDE:
    mit->add_excl();
    return 0;
  case DivSur::MultiSur::SmpRes::INNER:
    mit->add_inner();
    return 0;
  case DivSur::MultiSur::SmpRes::FAIL:
    mit->add_fail();
    return 0;
  case DivSur::MultiSur::SmpRes::FACET:
    // skip non-reactant facets
    if(reactant() >= 0 && smp_res.face.first != reactant() && smp_res.face.second != reactant()) {
      mit->add_skip();
      return 0;
    }

    sit = mit->find(smp_res.face);
    // new facet
    if(sit == mit->end() || !sit->second.face_num()) {
      IO::log << "      " << sur << "-th surface: new " << smp_res.face << " facet\n";
      (*mit)[smp_res.face].add_smp(_pot, _ms, sur, new_conf);
      return 1;
    }

    // maximum facet sampling number warning
    if(samp_warn && sit->second.samp_num() >= max_pot_size) {
      std::cerr << funame << sur << "-th surface, " << sit->first << " facet: "
		<< "WARNING: maximal facet samplings number has been reached\n";
      samp_warn = 0;
    } 

    // maximum importance sampling number warning
    if(imp_warn && sit->second.size() >= max_imp_size) {
      std::cerr << funame << sur << "-th surface, " << sit->first << " facet: "
		<< "WARNING: maximal importance samplings number has been reached\n";
      imp_warn = 0;
    } 

    // no potential calculation is needed
    if(face_work.find(smp_res.face) == face_work.end())
      sit->second.add_fake();
    // add facet samping
    else
      sit->second.add_smp(_pot, _ms, sur, new_conf);

    return 0;
  default:
    std::cerr << funame << "wrong case\n";
    throw Error::Logic();
  }
}

void CrossRate::MultiArray::_get_stat_flux(std::vector<double>& flux_mean, 
					   std::vector<double>& flux_rmsd, 
					   std::vector<double>& flux_variance) const
{
  flux_mean.clear();
  flux_mean.resize(_ms.species_size());
  flux_rmsd.clear();
  flux_rmsd.resize(_ms.species_size());
  flux_variance.clear();
  flux_variance.resize(_ms.species_size());

  double face_flux, face_rmsd, face_var;
  for(const_iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit) {// facet cycle
      if(reactant() >= 0 && reactant() != sit->first.first && reactant() != sit->first.second)
	continue;

      if(sit->second.size()) {
	face_flux = mit->face_flux(sit);
	face_rmsd = mit->vol_frac(sit) * std::sqrt(sit->second.flux_var());
	face_var  = face_flux * face_flux * (sit->second.flux_rel_var() + mit->vol_rel_var(sit));

	if(reactant() >= 0) {
	  flux_mean[reactant()]      += face_flux;
	  flux_rmsd[reactant()]      += face_rmsd;
	  flux_variance[reactant()]  += face_var;
	}
	else {
	  flux_mean[sit->first.first]      += face_flux;
	  flux_rmsd[sit->first.first]      += face_rmsd;
	  flux_variance[sit->first.first]  += face_var;

	  flux_mean[sit->first.second]      += face_flux;
	  flux_rmsd[sit->first.second]      += face_rmsd;
	  flux_variance[sit->first.second]  += face_var;
	}
      } 
    }// facet cycle
  }// surface cycle
}

void CrossRate::MultiArray::_get_reac_flux (std::map<DivSur::face_t, double>& flux_mean, 
					    std::map<DivSur::face_t, double>& flux_rmsd, 
					    std::map<DivSur::face_t, double>& flux_variance) const
{
  const char funame [] = "CrossRate::MultiArray::_get_reac_flux: ";

  flux_mean.clear();
  flux_rmsd.clear();
  flux_variance.clear();

  double dtemp;

  for(const_iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit)// facet cycle
      if(sit->second.run_traj_num()) {
	double face_flux = mit->face_flux(sit);
	for (int ward = 0; ward < 2; ++ward) {// ward cycle
	  int react;
	  switch(ward) {
	  case BACKWARD:
	    react = sit->first.second;
	    break;
	  case FORWARD:
	    react = sit->first.first;
	    break;
	  }
	  if(reactant() >= 0 && react != reactant())
	    continue;

	  for(int spec = 0; spec < _ms.species_size(); ++spec) {//product cycle
	    if(spec != react && sit->second.reac_traj_num(ward, spec)) {
	      // dynamical correction factor
	      double dyn_fac = (double)sit->second.reac_traj_num(ward, spec) / (double)sit->second.run_traj_num();
	      double dyn_var = dyn_fac * (1. - dyn_fac);

	      DivSur::face_t tran(react, spec);

	      // reactive flux
	      flux_mean[tran] += dyn_fac * face_flux;

	      // reactive flux standard deviation
	      flux_rmsd[tran] += std::sqrt(dyn_var) * face_flux;

	      // reactive flux variance
	      flux_variance[tran] += dyn_var * face_flux * face_flux / (double)sit->second.run_traj_num();
	    }
	  }// product cycle
	}// ward cycle
      }// facet cycle 
  }// surface cycle
} 

void CrossRate::MultiArray::_print_sampling_results () const
{
  IO::log << "Sampling results:\n";
  for(const_iterator cmit = begin(); cmit != end(); ++cmit) {// surface cycle
    if(!cmit->size())
      continue;
    IO::log << cmit - begin() << "-th surface:\n";
    if(cmit->tot_smp_num())
      IO::log << "     Number of all samplings:             " << cmit->tot_smp_num()  << "\n";
    if(cmit->fail_num())
      IO::log << "     Number of failed samplings:          " << cmit->fail_num()  << "\n";
    if(cmit->inner_num())
      IO::log << "     Number of inner region samplings:    " << cmit->inner_num() << "\n";
    if(cmit->excl_num())
      IO::log << "     Number of excluded region samplings: " << cmit->excl_num()  << "\n";
    if(cmit->close_num())
      IO::log << "     Number of close atoms samplings:     " << cmit->close_num() << "\n";
    if(cmit->skip_num())
      IO::log <<  "     Number of skipped facets samplings:  " << cmit->skip_num()  << "\n";

    for(SurArray::const_iterator cit = cmit->begin(); cit != cmit->end(); ++cit) {
      IO::log << "     " << cit->first << " facet:\n"
	      << "          Statistical flux:               "
	      << cmit->face_flux(cit) * norm_factor() <<  " (" 
	      << std::sqrt(cit->second.flux_rel_var() + cmit->vol_rel_var(cit)) * 100 << "%)\n"
	      << "          Number of importance samplings: " << cit->second.size() << "\n"
	      << "          Number of flux samplings:  " << cit->second.flux_num() << "\n";
      if(cit->second.fake_num())
	IO::log	<< "          Number of fake samplings:       " << cit->second.fake_num() << "\n";
      if(cit->second.fail_num())
	IO::log << "          Number of failed samplings:     " << cit->second.fail_num() << "\n";
      IO::log	<< "          Minimal energy configuration:\n"
		<< "             Energy, kcal/mol:       " 
		<< cit->second.min_ener() / Phys_const::kcal << "\n";
      cit->second.min_geom().print_geom(IO::log, "             ");
    }
    IO::log << std::endl;
  }

  std::vector<double> flux_mean, flux_rmsd, flux_variance;
  _get_stat_flux(flux_mean, flux_rmsd, flux_variance);

  // output
  int old_precision = IO::log.precision(3);
  IO::log << "Statistical flux:\n";
  IO::log << std::setw(13) << "Reactant" 
	    << std::setw(13) << "Flux"
	    << std::setw(13) << "Error, %"
	    << "\n";
  for(int spec = 0; spec < _ms.species_size(); ++spec)
    if(reactant() < 0 || spec == reactant()) {
      IO::log << std::setw(13) << spec
		<< std::setw(13) << flux_mean[spec] * norm_factor();
      if(flux_mean[spec] == 0.)
	IO::log << std::setw(13) << "***";
      else
	IO::log << std::setw(13) << 100. * std::sqrt(flux_variance[spec]) / flux_mean[spec];
      IO::log << "\n";
    }
  IO::log << "\n";
  IO::log.precision(old_precision);
}

void CrossRate::MultiArray::_print_traj_results () const
{
  int reac, prod;
  SharedPointer<DynRes> DynSmp::* region [2];

  double face_flux, dyn_fac;
  double dtemp;

  int old_precision = IO::log.precision(3);

  IO::log << "Trajectory results:\n"; 
  for(const_iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    if(mit->size())
      IO::log << mit - begin() << "-th surface:\n";
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit)// facet cycle
      if(sit->second.size()) {
	face_flux =  mit->face_flux(sit) * norm_factor();
	IO::log << "   " << sit->first <<  " facet:\n";
	//IO::log << "      " << "Volume fraction:  " << mit->vol_frac(sit) << "\n"; 
	//IO::log << "      " << "Statistical flux: " << face_flux << "\n"; 
	IO::log << "      " << "Number of trajectories run: " << sit->second.run_traj_num() << "\n";
	if(sit->second.potfail_traj_num())
	  IO::log << "      " << "Number of potential failed trajectories: " 
		    << sit->second.potfail_traj_num() << "\n";
	if(sit->second.runfail_traj_num())
	  IO::log << "      " << "Number of integrator failed trajectories: " 
		    << sit->second.runfail_traj_num() << "\n";
	if(sit->second.exclude_traj_num())
	  IO::log << "      " << "Number of exlcuded trajectories: " 
		    << sit->second.exclude_traj_num() << "\n";
	if(sit->second.init_traj_num())
	  IO::log << "      " << "Number of non running trajectories: " 
		    << sit->second.init_traj_num() << "\n";

	if(sit->second.run_traj_num()) 
	  for (int ward = 0; ward < 2; ++ward) {// ward cycle
	    switch(ward) {
	    case BACKWARD:
	      region[0] = &DynSmp::forw;
	      region[1] = &DynSmp::back;
	      reac = sit->first.second;
	      break;
	    case FORWARD:
	      region[0] = &DynSmp::back;
	      region[1] = &DynSmp::forw;
	      reac = sit->first.first;
	      break;
	    }

	    if(reactant() >= 0 && reac != reactant())
	      continue;

	    IO::log << "      " << reac << "-th reactant:\n";
	    IO::log << "         " << "Fraction (#) of recrossed trajectories: " 
		      << (double)sit->second.recross_traj_num(ward) / (double)sit->second.run_traj_num() 
		      << " (" << sit->second.recross_traj_num(ward) <<")\n";
	    IO::log << "         " << "Fraction (#) of returned  trajectories: " 
		      << (double)sit->second.reac_traj_num(ward, reac)  / (double)sit->second.run_traj_num()
		      << " (" << sit->second.reac_traj_num(ward, reac) <<")\n";

	    IO::log << "         "   << "Reactive trajectories:\n";
	    IO::log << std::setw(20) << "R -> P" 
		    << std::setw(15) << "Reactive Flux"
		    << std::setw(15) << "Fraction"
		    << std::setw(15) << "Number"
		    << std::setw(15) << "Error, %"
		    << "\n";

	    for(prod = 0; prod < _ms.species_size(); ++prod) {// product cycle
	      if(prod == reac)
		continue;
	      
	      std::ostringstream oss;
	      oss << reac << " -> " << prod;
	      dyn_fac = (double)sit->second.reac_traj_num(ward, prod) / (double)sit->second.run_traj_num();
	      IO::log << std::setw(20)  << oss.str()
			<< std::setw(15) << dyn_fac * face_flux
			<< std::setw(15) << dyn_fac
			<< std::setw(15)  << sit->second.reac_traj_num(ward, prod);
	      if(sit->second.reac_traj_num(ward, prod)) {
		dtemp = 1./(double)sit->second.reac_traj_num(ward, prod) - 1./(double)sit->second.run_traj_num();
		dtemp = dtemp > 0. ? std::sqrt(dtemp) * 100. : 0.;
		IO::log << std::setw(15) << dtemp;
	      }
	      else
		IO::log << std::setw(15) << "***";
	      IO::log << "\n";
	    }//product cycle
	    IO::log << "\n";

	    // energy and angular momentum conservation error output
	    double ee[2][2];
	    int ee_count;

	    IO::log << "         "   << "Energy conservation error (kcal/mol):\n"
		    << std::setw(20) << "R -> P" 
		    << std::setw(15) << "Direction"
		    << std::setw(15) << "Average"
		    << std::setw(15) << "Max"
		    << std::setw(15) << "Traj. #"
		    << "\n";

	    for(prod = 0; prod < _ms.species_size(); ++prod) {// product cycle
	      if(prod == reac)
		continue;

	      ee_count = 0;
	      for(int i = 0; i < 2; ++i)
		for(int ireg = 0; ireg < 2; ++ireg)
		  ee[ireg][i] = 0.;

	      for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit) {//sampling cycle
		if(!fit->is_run() || fit->is_run_fail() || fit->is_pot_fail() || fit->is_exclude())
		  continue;

		const DynSmp* s = fit.operator->();
		if((s->*region[0])->stat == DynRes::DIRECT && (s->*region[1])->species() == prod) {
		  ee_count++;
		  for(int ireg = 0; ireg < 2; ++ireg) { 
		    dtemp = (s->*region[ireg])->total_energy() - fit->total_energy();
		    dtemp = dtemp > 0. ? dtemp : -dtemp;
		    ee[ireg][0] += dtemp;
		    ee[ireg][1] = dtemp > ee[ireg][1] ? dtemp : ee[ireg][1];
		  }
		}
	      }// sampling cycle
	      
	      // normalization
	      if(ee_count)
		for(int ireg = 0; ireg < 2; ++ireg)
		  ee[ireg][0] /= (double)ee_count;
		
	      
	      // output
	      std::ostringstream sso;
	      sso << reac << " -> " << prod;

	      for(int ireg = 0; ireg < 2; ++ireg) {
		IO::log << std::setw(20) << sso.str();
		if(!ireg)
		  IO::log << std::setw(15) << "backward";
		else
		  IO::log << std::setw(15) << "forward";
		for(int i = 0; i < 2; ++i)
		  IO::log << std::setw(15) << ee[ireg][i] / Phys_const::kcal;
		IO::log << std::setw(15) << ee_count << "\n";
	      }
	    }// product cycle
	    IO::log << "\n";
	    
	    IO::log << "         "   << "Angular momentum conservation error (a.u.):\n"
		    << std::setw(20) << "R -> P" 
		    << std::setw(15) << "Direction"
		    << std::setw(15) << "Average"
		    << std::setw(15) << "Max"
		    << std::setw(15) << "Traj. #"
		    << "\n";

	    for(prod = 0; prod < _ms.species_size(); ++prod) {// product cycle
	      if(prod == reac)
		continue;

	      ee_count = 0;
	      for(int i = 0; i < 2; ++i)
		for(int ireg = 0; ireg < 2; ++ireg)
		  ee[ireg][i] = 0.;

	      for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit) {//sampling cycle
		if(!fit->is_run() || fit->is_run_fail() || fit->is_pot_fail() || fit->is_exclude())
		  continue;

		const DynSmp* s = fit.operator->();
		if((s->*region[0])->stat == DynRes::DIRECT && (s->*region[1])->species() == prod) {
		  ee_count++;
		  for(int ireg = 0; ireg < 2; ++ireg) { 
		    dtemp = ((s->*region[ireg])->total_angular_momentum() - fit->total_angular_momentum()).vlength();
		    ee[ireg][0] += dtemp;
		    ee[ireg][1] = dtemp > ee[ireg][1] ? dtemp : ee[ireg][1];
		  }
		}
	      }// sampling cycle
	      
	      // normalization
	      if(ee_count)
		for(int ireg = 0; ireg < 2; ++ireg)
		  ee[ireg][0] /= (double)ee_count;
	      
	      // output
	      std::ostringstream sso;
	      sso << reac << " -> " << prod;

	      for(int ireg = 0; ireg < 2; ++ireg) {
		IO::log << std::setw(20) << sso.str();
		if(!ireg)
		  IO::log << std::setw(15) << "backward";
		else
		  IO::log << std::setw(15) << "forward";
		for(int i = 0; i < 2; ++i)
		  IO::log << std::setw(15) << ee[ireg][i];
		IO::log << std::setw(15) << ee_count << "\n";
	      }
	    }// product cycle
	    IO::log << "\n";

	  }// ward cycle
      }// facet cycle 
  }// surface cycle
  IO::log << "\n";


  // crossrate output
  std::map<DivSur::face_t, double> flux_mean, flux_rmsd, flux_variance;
  _get_reac_flux(flux_mean, flux_rmsd, flux_variance);
  
  IO::log << "Reactive transition fluxes:\n"
	  << std::setw(13) << "Transition"
	  << std::setw(11) << "Flux"
	  << std::setw(9) << "Err, %"
	  << "\n";
  for(std::map<DivSur::face_t, double>::const_iterator it = flux_mean.begin(); it != flux_mean.end(); ++it) {
    if(reactant() >= 0 && it->first.first != reactant())
      continue; 

    IO::log << std::setw(13) << it->first 
	    << std::setw(11) << it->second * norm_factor();

    if(it->second == 0.)
      IO::log<< std::setw(9) << "***";
    else
      IO::log<< std::setw(9) << 100. * std::sqrt(flux_variance.find(it->first)->second) / it->second;

    IO::log << "\n";
  }
  IO::log << "\n";

  std::vector<double> reac_flux(_ms.species_size());
  std::vector<double> reac_var(_ms.species_size());
  for(std::map<DivSur::face_t, double>::const_iterator it = flux_mean.begin(); it != flux_mean.end(); ++it) {
    if(reactant() >= 0 && it->first.first != reactant())
      continue; 
    reac_flux[it->first.first] += it->second;
    reac_var[it->first.first]  += flux_variance.find(it->first)->second;
  }

  IO::log << "Cumulative reactive fluxes:\n"
	  << std::setw(3)  << "R"
	  << std::setw(11) << "Flux"
	  << std::setw(9)  << "Err, %"
	  << "\n";
  for(int r =  0; r < _ms.species_size(); ++r) {
    if(reactant() >= 0 && r != reactant())
      continue; 

    dtemp = reac_flux[r];
    IO::log << std::setw(3)  << r 
	    << std::setw(11) << dtemp * norm_factor();

    if( dtemp == 0.)
      IO::log << std::setw(9) << "***";
    else
      IO::log << std::setw(9) << 100. * std::sqrt(reac_var[r]) / dtemp;

    IO::log << "\n";
  }
  IO::log << "\n";

  IO::log.precision(old_precision);
  
}

void CrossRate::MultiArray::_analize () const
{
  double dtemp;
  DynRes::Stat stat;
  int reac, prod; 
  SharedPointer<DynRes> DynSmp::* reac_p;
  SharedPointer<DynRes> DynSmp::* prod_p;

  int old_precision;

  if(raden_flag || amproj_flag) {// auxiliary output
    std::ofstream amproj_stream;
    std::ofstream raden_stream;
    int m_proj_out;
    if(amproj_flag) {
      amproj_stream.open(amproj_file.c_str());
      amproj_stream << "Abbreviations:\n"
		    << "\t_reac, _tran, and _prod subscripts refer to the reactant, transition state, and product regions\n"
		    << "\tR -> P: reactant to product trajectory\n"
		    << "\tK: the projection of the total angular momentum on the interfragment axis\n"
		    << "\tJ: is the total angular momentum\n"
		    << "\td: interfragment distance\n"
		    << "\tj: the length of the vector sum of the angular momenta of the fragments\n"
		    << "\tall quantities are in atomic units\n"
		    << "\n";

      m_proj_out = 0;
      for(int frag = 0; frag < 2; ++frag)
	if(Structure::fragment(frag).type() == Molecule::NONLINEAR &&
	   (Structure::fragment(frag).top() == Molecule::PROLATE || 
	    Structure::fragment(frag).top() == Molecule::OBLATE)) {
	  m_proj_out = 1;
	  break;
	}
    }
    if(raden_flag)
      raden_stream.open(raden_file.c_str());    
    
    for(const_iterator mit = begin(); mit != end(); ++mit) {// surface cycle
      for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit)// facet cycle
	if(sit->second.size()) {
	  std::vector<std::multiset<double> > raden[2];

	  if(raden_flag)
	    for(int ward = 0; ward < 2; ++ward) {
	      raden[ward].resize(_ms.species_size());
	    }

	  if(amproj_flag) {// angular momentum projection output
	    old_precision = amproj_stream.precision(3);
	    amproj_stream << mit - begin() << "-th surface, " << sit->first << " facet:\n";
	    for(int ward = 0; ward < 2; ++ward) {// ward cycle
	      switch(ward) {
	      case BACKWARD:
		reac = sit->first.second;
		break;
	      case FORWARD:
		reac = sit->first.first;
		break;
	      }
	      if(reactant() >= 0 && reac != reactant())
		continue;
	      for(prod = 0; prod < _ms.species_size(); ++prod) {// product cycle
		if(prod == reac)
		  continue;
		  
		for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit) {//sampling cycle
		  if(fit->is_run_fail() || fit->is_pot_fail() || fit->is_exclude())
		    continue;
		  const DynSmp* s = fit.operator->();
		  switch(ward) {
		  case BACKWARD:
		    reac_p = &DynSmp::forw;
		    prod_p = &DynSmp::back;
		    break;
		  case FORWARD:
		    reac_p = &DynSmp::back;
		    prod_p = &DynSmp::forw;
		    break;
		  }
		  if((s->*reac_p)->stat == DynRes::DIRECT && (s->*prod_p)->species() == prod) {// output

		    amproj_stream  << std::setw(10) << "R -> P"
				   << std::setw(10) << "d_reac" 
				   << std::setw(10) << "d_tran" 
				   << std::setw(10) << "d_prod" 
				   << std::setw(10) << "J" 
				   << std::setw(14) << "E, kcal/mol"
				   << "\n";

		    std::ostringstream sso;
		    sso << (s->*reac_p)->species() << " -> " << (s->*prod_p)->species();
		    amproj_stream << std::setw(10) << sso.str()
				  << std::setw(10) << (s->*reac_p)->interfragment_distance()
				  << std::setw(10) << fit->interfragment_distance()
				  << std::setw(10) << (s->*prod_p)->interfragment_distance()
				  << std::setw(10) << fit->total_angular_momentum().vlength()
				  << std::setw(14) << fit->total_energy() / Phys_const::kcal
				  << "\n\n";

		    amproj_stream << std::setw(20) << "Fragment"
				  << std::setw(10) << "K_reac" 
				  << std::setw(10) << "K_tran" 
				  << std::setw(10) << "K_prod";
		    if(m_proj_out)
		      amproj_stream << std::setw(10) << "M_reac" 
				    << std::setw(10) << "M_tran" 
				    << std::setw(10) << "M_prod";
		    amproj_stream << "\n";

		    for(int frag = 0; frag < 2; ++frag) {
		      amproj_stream << std::setw(20) << frag
				    << std::setw(10) << (s->*reac_p)->angular_momentum_k_projection(frag)
				    << std::setw(10) << fit->angular_momentum_k_projection(frag)
				    << std::setw(10) << (s->*prod_p)->angular_momentum_k_projection(frag);
		      if(Structure::fragment(frag).type() == Molecule::NONLINEAR &&
			 (Structure::fragment(frag).top() == Molecule::PROLATE || 
			  Structure::fragment(frag).top() == Molecule::OBLATE))
			amproj_stream  << std::setw(10) << (s->*reac_p)->angular_momentum_m_projection(frag)
				       << std::setw(10) << fit->angular_momentum_m_projection(frag)
				       << std::setw(10) << (s->*prod_p)->angular_momentum_m_projection(frag);
		      amproj_stream << "\n";
		    }
		    amproj_stream << std::setw(20) << "Sum"
				  << std::setw(10) << (s->*reac_p)->angular_momentum_k_projection(0)
		      + (s->*reac_p)->angular_momentum_k_projection(1)
				  << std::setw(10) << fit->angular_momentum_k_projection(0) 
		      + fit->angular_momentum_k_projection(1)
				  << std::setw(10) << (s->*prod_p)->angular_momentum_k_projection(0)
		      + (s->*prod_p)->angular_momentum_k_projection(1)
				  << "\n\n";
		  }//output
		}// sampling cycle
		amproj_stream << "\n";
	      }// product cycle
	    }// ward cycle
	    amproj_stream.precision(old_precision);
	  }// angular momentum projection output

	  if(raden_flag) {// radial energy output
	    for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit) { // sampling cycle
	      if(fit->is_run_fail() || fit->is_pot_fail() || fit->is_exclude())
		continue;

	      if(fit->back->stat == DynRes::DIRECT || fit->forw->stat == DynRes::DIRECT)
		for(int ward = 0; ward < 2; ++ward) {
		  switch(ward) {
		  case BACKWARD:
		    stat  = fit->forw->stat;
		    prod  = fit->back->species();
		    break;
		  case FORWARD:
		    stat  = fit->back->stat;
		    prod  = fit->forw->species();
		    break;
		  }
		  if(stat == DynRes::DIRECT) {
		    raden[ward][prod].insert(fit->radial_kinetic_energy());
		  }
		}
	    }// sampling cycle

	    for (int ward = 0; ward < 2; ++ward) {// ward cycle
	      switch(ward) {
	      case BACKWARD:
		reac = sit->first.second;
		break;
	      case FORWARD:
		reac = sit->first.first;
		break;
	      }
	      if(reactant() >= 0 && reac != reactant())
		continue;
	      for(prod = 0; prod < _ms.species_size(); ++prod)// product cycle
		if(reac != prod) {
		  raden_stream << mit - begin() << "-th surface, "
			       << sit->first << " facet, "
			       << reac << " -> " << prod << " transition:\n"
			       << std::setw(15) << "E, kcal/mol" << std::setw(15) << "N(E)\n";
		  int i = 0;
		  for(std::multiset<double>::const_iterator rit = raden[ward][prod].begin(); 
		      rit != raden[ward][prod].end(); ++rit)
		    raden_stream << *rit / Phys_const::kcal << "   " 
				 << (double)(i++)/(double)raden[ward][prod].size() << "\n";
		  raden_stream << "\n";
		  
		}// product cycle 
	    }// ward cycle
	  }// radial energy output
	}// facet cycle 
    }// surface cycle
  }// auxiliary output
}

void CrossRate::print_progress (int& old_share, int new_share)
{
  static std::ios::pos_type p0;
  if(new_share > old_share) {
    if(old_share)
      IO::log.seekp(p0);
    //for(int i = 0; i < 106; ++i)
    //IO::log << "\b";
    else {
      p0 = IO::log.tellp();
    }
    old_share = new_share;
    for(int i = 0; i < 100; ++i)
      if(i < new_share)
	IO::log << ".";
      else
	IO::log << " ";
    IO::log << std::setw(5) << new_share << "%\n";
    IO::log.flush();
  }
}

bool CrossRate::MultiArray::_work (Dynamic::CCP stop) 
{
  const char funame [] = "CrossRate::MultiArray::_work: ";

  double dtemp;
  int itemp;
  bool btemp;

  int count, count_max;
  int old_share, new_share;

  SurArray::iterator sit;

  std::map<DivSur::face_t, int> proj_samp_num; // projected number of samplings per facet
  std::map<DivSur::face_t, int>::iterator nit;

  /**********************************************************************
   ********** SATISFY  MINIMAL NON-ZERO FLUX SAMPLINGS NUMBER ***********
   *********************************************************************/

  IO::log << "Sampling to satisfy minimal flux samplings number (flux samplings counted) ...\n";
  for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle

    count_max = 0;
    for(sit = mit->begin(); sit != mit->end(); ++sit) {
      itemp = min_pot_size - sit->second.samp_num();
      if(itemp > 0)
	count_max += itemp;
    }

    if(!count_max)
      continue;

    IO::log << "   " << mit - begin() << "-th surface:\n";
    old_share = 0;
    while(1) {// sampling loop
      std::set<DivSur::face_t > face_work;

      count = 0;
      for(sit = mit->begin(); sit != mit->end(); ++sit) {
	itemp = min_pot_size - sit->second.samp_num();
	if(itemp > 0) { 
	  face_work.insert(sit->first);
	  count += itemp;
	}
      }

      new_share = int((double)(count_max - count) / (double)count_max * 100.);
      print_progress(old_share, new_share);

      // exit condition
      if(!face_work.size())
	break;

      if(_sample(mit, face_work))
	return false;
    }// sampling loop
  }// surface cycle
  IO::log << "done\n\n";

  /**********************************************************************
   ********** SATISFY FACET STATISTICAL FLUX RELATIVE TOLERANCE *********
   *********************************************************************/

  IO::log << "Sampling to satisfy facet statistical flux relative tolerance (flux samplings counted) ...\n";
  for(iterator mit = begin(); mit != end(); ++mit) {// primitive cycle
    proj_samp_num.clear();

    for(sit = mit->begin(); sit != mit->end(); ++sit) {
      dtemp = face_rel_tol * sit->second.flux_val();
      if(dtemp != 0.)
	proj_samp_num[sit->first] = int(sit->second.flux_var() / dtemp / dtemp);
    }

    if(proj_samp_num.size()) {
      IO::log << "   " << mit - begin() << "-th surface:\n"
	      << std::setw(15) << "Facet" 
	      << std::setw(15) << "Projected" 
	      << std::setw(15) << "Current"  
	      << "\n";

      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	IO::log << std::setw(15) << nit->first 
		<< std::setw(15) << nit->second 
		<< std::setw(15) << mit->find(nit->first)->second.samp_num() 
		<< "\n";
    }

    count_max = 0;
    for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
      sit = mit->find(nit->first);
      itemp =  nit->second - sit->second.samp_num();
      if(itemp > 0) {
	count_max += itemp; 
      }
    }
    
    if(!count_max)
      continue;

    old_share = 0;
    while(1) {// sampling loop
      std::set<DivSur::face_t > face_work;

      count = 0;
      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	sit = mit->find(nit->first);
	itemp =  nit->second - sit->second.samp_num();
	if(itemp > 0) {
	  count += itemp; 
	  if(sit->second.samp_num() < max_pot_size)
	    face_work.insert(nit->first);
	}
      }

      new_share = int((double)(count_max - count) / (double)count_max * 100.);
      print_progress(old_share, new_share);

      // exit condition
      if(!face_work.size())
	break;

      if(_sample(mit, face_work))
	return false;
    }// sampling loop
  }// primitive cycle
  IO::log << "done\n\n";

  /**********************************************************************
   ******** SATISFY SPECIES STATISTICAL FLUX RELATIVE TOLERANCE *********
   *********************************************************************/

  std::vector<double> stat_flux_mean, stat_flux_rmsd, stat_flux_variance;
  _get_stat_flux(stat_flux_mean, stat_flux_rmsd, stat_flux_variance);

  IO::log << "Sampling to satisfy species statistical flux relative tolerance (flux samplings counted) ... \n";

  for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    proj_samp_num.clear();
    for(sit = mit->begin(); sit != mit->end(); ++sit) {
      int react;
      for(int ward = 0; ward < 2; ++ward) {// ward cycle
	switch(ward) {
	case BACKWARD:
	  react = sit->first.second;
	  break;
	case FORWARD:
	  react = sit->first.first;
	  break;
	}

	if(reactant() >= 0 && reactant() != react) 
	  continue;

	// projected number of samplings
	dtemp = spec_rel_tol * stat_flux_mean[react];
	if(dtemp != 0.) {
	  double face_rmsd = mit->vol_frac(sit) * std::sqrt(sit->second.flux_var());
	  itemp = (int)(face_rmsd * stat_flux_rmsd[react] / dtemp / dtemp);

	  nit = proj_samp_num.find(sit->first);
	  if(nit == proj_samp_num.end() || itemp > nit->second)
	    proj_samp_num[sit->first] = itemp;
	}
      }// ward cycle
    }// facet cycle
    
    if(proj_samp_num.size()) {
      IO::log << "   " << mit - begin() << "-th surface:\n"
	      << std::setw(15) << "Facet" 
	      << std::setw(15) << "Projected" 
	      << std::setw(15) << "Current"  
	      << "\n";

      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	IO::log << std::setw(15) << nit->first
		<< std::setw(15) << nit->second 
		<< std::setw(15) << mit->find(nit->first)->second.samp_num() 
		<< "\n";
    }

    count_max = 0;
    for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
      sit = mit->find(nit->first);
      itemp =  nit->second - sit->second.samp_num();
      if(itemp > 0) {
	count_max += itemp; 
      }
    }

    if(!count_max)
      continue;

    old_share = 0;
    while(1) {// sampling loop
      std::set<DivSur::face_t > face_work;
      count = 0;
      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	sit = mit->find(nit->first);
	itemp =  nit->second - sit->second.samp_num();
	if(itemp > 0) {
	  count += itemp; 
	  if(sit->second.samp_num() < max_pot_size)
	    face_work.insert(nit->first);
	}
      }
      
      new_share = int((double)(count_max - count) / (double)count_max * 100.);
      print_progress(old_share, new_share);

      if(!face_work.size())
	break;

      if(_sample(mit, face_work))
	return false;
    }// sampling loop
  }// surface cycle
  IO::log << "done\n\n";

  if(job() == DYN_JOB) {

    /****************************** DYNAMICAL CALCULATIONS ***************************/

    /*************************************************************
     ******* SATISFY MINIMAL IMPORTANCE SAMPLINGS NUMBER *********
     *************************************************************/

    IO::log << "Sampling to satisfy minimal importance samplings number (importance samplings counted) ...\n";
    for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle

      IO::log << "   " << mit - begin() << "-th surface:\n"
	      << std::setw(15) << "Facet" 
	      << std::setw(15) << "Projected" 
	      << std::setw(15) << "Current" << "\n";

      count_max = 0;
      for(sit = mit->begin(); sit != mit->end(); ++sit) {
	IO::log << std::setw(15) << sit->first
		<< std::setw(15) << min_imp_size
		<< std::setw(15) << sit->second.size()
		<< "\n";
	itemp = min_imp_size - sit->second.size();
	if(sit->second.size() && itemp > 0)
	  count_max += itemp;
	
      }

      if(!count_max)
	continue;
      
      old_share = 0;
      while(1) {// sampling loop
	// tests
	std::set<DivSur::face_t > face_work;
	count = 0;
	for(sit = mit->begin(); sit != mit->end(); ++sit)  {
	  itemp = min_imp_size - sit->second.size();
	  if(sit->second.size() && itemp > 0) {
	    if(sit->second.samp_num() >= max_pot_size) {
	      std::cerr << funame << mit - begin() << "-th surface, " << sit->first << " facet: "
			<< "cannot satisfy minimal importance sampling number: "
			<< "maximal number of flux samplings has been reached\n";
	      throw Error::Run();
	    }
	    count += itemp;
	    face_work.insert(sit->first);
	  }
	}

	new_share = int((double)(count_max - count) / (double)count_max * 100.);
	print_progress(old_share, new_share);

	// exit condition
	if(!face_work.size())
	  break;

	if(_sample(mit, face_work))
	  return false;
      }// sampling loop
    }// surface cycle
    IO::log << "done\n\n";
  
    _print_sampling_results();
    
    /*************************************************************
     *********************** RUN TRAJECTORIES ********************
     *************************************************************/

    IO::log << "Running trajectories ...\n";
    run_traj(stop);
    IO::log << "done\n\n";

    _print_traj_results();

    /*************************************************************
     * SATISFY SPECIES OUTGOING REACTIVE FLUX RELATIVE TOLERANCE *
     *************************************************************/

    std::map<DivSur::face_t, double> cross_flux_mean, cross_flux_rmsd, cross_flux_variance;
    _get_reac_flux(cross_flux_mean, cross_flux_rmsd, cross_flux_variance);

    std::vector<double> reac_flux(_ms.species_size());
    std::vector<double> reac_rmsd(_ms.species_size());
    for(std::map<DivSur::face_t, double>::const_iterator it = cross_flux_mean.begin();
	it != cross_flux_mean.end(); ++it) {
      reac_flux[it->first.first] += it->second;
      reac_rmsd[it->first.first] += cross_flux_rmsd.find(it->first)->second;
    }

    IO::log << "Sampling to satisfy species outgoing reactive flux relative tolerance (importance samplings counted) ...\n";
    for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
      // estimate necessary number of importance samplings for each facet on a given surface
      proj_samp_num.clear();
      int react;
      for(sit = mit->begin(); sit != mit->end(); ++sit) {// facet cycle
	if(sit->second.run_traj_num())
	  for(int ward = 0; ward < 2; ++ward) {
	    switch(ward) {
	    case BACKWARD:
	      react = sit->first.second;
	      break;
	    case FORWARD:
	      react = sit->first.first;
	      break;
	    }
	    
	    if(reactant() >= 0 && reactant() != react)
	      continue;
	    
	    for(int spec = 0; spec < _ms.species_size(); ++spec) {//product cycle
	      if(spec != react && (itemp = sit->second.reac_traj_num(ward, spec))) {

		// recrossing factor
		dtemp = (double)itemp / (double)sit->second.run_traj_num();

		// facet reactive flux standard deviation
		double face_rmsd = std::sqrt(dtemp * (1. - dtemp)) * mit->face_flux(sit);
		
		dtemp = reac_rel_tol * reac_flux[react];
		itemp = (int)(face_rmsd * reac_rmsd[react] / dtemp / dtemp);
		itemp = int((double)itemp * (double)sit->second.size() / (double)sit->second.run_traj_num());

		nit = proj_samp_num.find(sit->first);
		if(nit == proj_samp_num.end() || itemp > nit->second)
		  proj_samp_num[sit->first] = itemp;
	      }
	    }// product cycle
	  } // ward cycle
      }// facet cycle 
    
      if(proj_samp_num.size()) {
	IO::log << "   " << mit - begin() << "-th surface:\n"
		<< std::setw(15) << "Facet" 
		<< std::setw(15) << "Projected" 
		<< std::setw(15) << "Current" << "\n";

	for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	  IO::log << std::setw(15) << nit->first 
		  << std::setw(15) << nit->second 
		  << std::setw(15) << mit->find(nit->first)->second.size() 
		  << "\n";
      }

      count_max = 0;
      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	sit = mit->find(nit->first);
	itemp = nit->second - sit->second.size();
	if(itemp > 0) {
	  count_max += itemp;
	}
      }

      if(!count_max)
	continue;

      old_share = 0;
      while(1) {// sampling loop
	std::set<DivSur::face_t > face_work;
	count = 0;
	for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	  sit = mit->find(nit->first);
	  itemp = nit->second - sit->second.size();
	  if(itemp > 0) {
	    count += itemp;
	    if(sit->second.samp_num() < max_pot_size && sit->second.size() < max_imp_size)
	      face_work.insert(nit->first);
	  }
	}

	new_share = int((double)(count_max - count) / (double)count_max * 100.);
	print_progress(old_share, new_share);

	if(!face_work.size())
	  break;

	if(_sample(mit, face_work))
	  return false;
      }// sampling loop
    }// surface cycle
    IO::log << "done\n\n";
  
    /*******************************************************
     * SATISFY REACTIVE TRANSITION FLUX RELATIVE TOLERANCE *
     *******************************************************/

    if(reactive_transition.size()) {
      // checking
      if(reactive_transition.size() != 2 || reactive_transition[0] == reactive_transition[1] 
	 || reactive_transition[0] < 0 || reactive_transition[0] >= _ms.species_size()
	 || reactive_transition[1] < 0 || reactive_transition[1] >= _ms.species_size()) {
	std::cerr << funame << "reactive transition pair is not properly defined\n";
	throw Error::Init();
      }

      IO::log << "Sampling to satisfy reactive transition flux relative tolerance (importance samplings counted) ...\n";
      for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
	proj_samp_num.clear();
	for(sit = mit->begin(); sit != mit->end(); ++sit)
	  if(sit->second.run_traj_num()) {// facet cycle
	    int react, prod = -1;

	    for(int ward = 0; ward < 2; ++ward) {
	      switch(ward) {
	      case BACKWARD:
		react = sit->first.second;
		break;
	      case FORWARD:
		react = sit->first.first;
		break;
	      }
	    
	      if(reactant() >= 0 && reactant() != react)
		continue;
	    
	      for(int i = 0; i < 2; ++i)
		if(reactive_transition[i] == react)
		  prod = reactive_transition[1 - i];
	    
	      if(prod < 0)
		continue;

	      if(itemp = sit->second.reac_traj_num(ward, prod)) {
	      
		// recrossing factor
		dtemp = (double)itemp / (double)sit->second.run_traj_num();

		// facet reqctive flux standard deviation
		double face_rmsd = std::sqrt(dtemp * (1. - dtemp)) * mit->face_flux(sit);
		
		dtemp = tran_rel_tol  * cross_flux_mean[DivSur::face_t(react, prod)];
		itemp = int(face_rmsd * cross_flux_rmsd[DivSur::face_t(react, prod)] / dtemp / dtemp);
		itemp = int((double)itemp * (double)sit->second.size() / (double)sit->second.run_traj_num());

		nit = proj_samp_num.find(sit->first);
		if(nit == proj_samp_num.end() || itemp > nit->second)
		  proj_samp_num[sit->first] = itemp;
	      }
	    } // ward cycle
	  }// facet cycle
    
	if(proj_samp_num.size()) {
	  IO::log << "   " << mit - begin() << "-th surface:\n"
		  << std::setw(15) << "Facet" 
		  << std::setw(15) << "Projected" 
		  << std::setw(15) << "Current" << "\n";

	  for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	    IO::log << std::setw(15) << nit->first 
		    << std::setw(15) << nit->second 
		    << std::setw(15) << mit->find(nit->first)->second.size() 
		    << "\n";
	}

	count_max = 0;
	for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	  sit = mit->find(nit->first);
	  itemp = nit->second - sit->second.size();
	  if(itemp > 0) {
	    count_max += itemp;
	  }
	}

	if(!count_max)
	  continue;

	old_share = 0;
	while(1) {// sampling loop
	  std::set<DivSur::face_t > face_work;
	  count = 0;
	  for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	    sit = mit->find(nit->first);
	    itemp = nit->second - sit->second.size();
	    if(itemp > 0) {
	      count += itemp;
	      if(sit->second.samp_num() < max_pot_size && sit->second.size() < max_imp_size)
		face_work.insert(nit->first);
	    }
	  }

	  new_share = int((double)(count_max - count) / (double)count_max * 100.);
	  print_progress(old_share, new_share);

	  if(!face_work.size())
	    break;

	  if(_sample(mit, face_work))
	    return false;
	}// sampling loop
      }// surface cycle
      IO::log << "done\n\n";
    }

    /*************************************************************
     *********************** RUN TRAJECTORIES ********************
     *************************************************************/

    // run trajectories
    IO::log << "Running trajectories again ...\n";
    run_traj(stop);
    IO::log << "done\n\n";
  
  }

  _print_sampling_results();

  if(job() == DYN_JOB) {
    _print_traj_results();
    _analize();
  }

  return true;
} 

CrossRate::MultiArray::MultiArray(const DivSur::MultiSur& ms, Potential::Wrap pot, Dynamic::CCP stop)
   : std::vector<SurArray>(ms.primitive_size()), _ms(ms), _pot(pot)
{
  const char funame [] = "CrossRate::MultiArray::MultiArray: ";
 
  static const int max_fail = 10;

  SurArray::iterator sit;

  // initial sampling to check which facets are actually present
  IO::log << "Preliminary sampling with " << min_sur_size 
	  << " samplings per primitive surface ...\n";
  for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    IO::log << "   " << mit - begin() << "-th surface:\n";
    for(int s = 0; s < min_sur_size; ++s) {// sampling loop

      std::set<DivSur::face_t > face_work;
      for(sit = mit->begin(); sit != mit->end(); ++sit)
	face_work.insert(sit->first);
      _sample(mit, face_work);
    }// sampling loop
  }// surface cycle
  IO::log << "done\n\n";

  // check for failures
  for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    if(mit->fail_num() > max_fail) {
      std::cerr << funame << "too many failures for " << mit - begin() << "-th surface\n";
      throw Error::Run();
    }
    for(sit = mit->begin(); sit != mit->end(); ++sit)
      if(!sit->second.flux_num()) {
	std::cerr << funame << mit - begin() << "-th surface, " << sit->first << " facet: "
		  <<"no flux samplings\n";
	throw Error::Run();
      }
  }// surface cycle

  // sample the surface and run trajectories to satisfy predefined tolerances
  while(!_work(stop)) {}
}
