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
#include "limits.hh"

#include <iomanip>
#include <cmath>
#include <set>

namespace CrossRate {
  //
  bool _isinit = false;
  bool isinit () { return _isinit; }

  // random potential error flag; opposite assumes that potential error corresponds to large positive potential
  //
  int rand_pot_err_flag = 0;

  // output flags and files
  //
  int         raden_flag;
  std::string raden_file;
  int         amproj_flag;
  std::string amproj_file;

  // minimum number of surface samplings to find facets
  //
  int MultiArray::min_sur_size;

  // minimal importance samplings array size for the facet
  //
  int MultiArray::min_imp_size;

  // maximum importance samplings array size for the facet
  //
  int MultiArray::max_imp_size;

  // minimum number of facet samplings before the accuracy can be estimated
  //
  int MultiArray::min_pot_size;

  // number of samplings performed at once
  //
  int MultiArray::smp_set_size;

  // maximal primitive's sampling failurs
  //
  int max_fail_num;

  // tolerances
  //
  double MultiArray::face_rel_tol;  // facet   statistical flux relative tolerance
  double MultiArray::spec_rel_tol;  // species statistical flux relative tolerance
  double MultiArray::reac_rel_tol;  // species    reactive flux relative tolerance
  double MultiArray::tran_rel_tol;  // reactive transition flux relative tolerance

  // energy and angular momentum trajectory deviation thresholds
  //
  double MultiArray::amom_dev_max = -1.;
  double MultiArray::ener_dev_max = -1.;
  
  // reactive transition
  //
  std::vector<int> MultiArray::reactive_transition;

  // job type
  //
  job_t _job;
  job_t job () { return _job; }
   
  // calculation mode
  //
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
}

void CrossRate::init (std::istream& from) 
{
  const char funame [] = "CrossRate::init: ";

  if(isinit()) {
    //
    std::cerr << funame << "has been initialized already\n";

    throw Error::Init();
  }
  _isinit = true;

  if(!Structure::isinit()) {
    //
    std::cerr << funame << "structure should be initialized first\n";

    throw Error::Init();
  }

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;

  std::map<std::string, Read> input;
  std::map<std::string, Read>::iterator idit;

  double trt, exp_max;
  std::string calc_mode, job_type;

  input ["JobType"                ] = Read(job_type,  "dynamical");
  input ["CalculationMode"        ] = Read(calc_mode, "canonical");
  input ["Reactant"               ] = Read(_reactant, -1);
  input ["ReactiveTransition"     ] = Read(MultiArray::reactive_transition, std::vector<int>());
  input ["Temperature"            ] = Read(_temperature, -1.);
  input ["Energy"                 ] = Read(_tot_ener, -100.);
  input ["AngularMomentum"        ] = Read(_tot_amom, -1.);
  input ["SamplingSize"           ] = Read(MultiArray::smp_set_size, 1000);
  input ["SurfaceSamplingMin"     ] = Read(MultiArray::min_sur_size, 1000000);
  input ["PotentialSamplingMin"   ] = Read(MultiArray::min_pot_size, 100);
  input ["ImportanceSamplingMin"  ] = Read(MultiArray::min_imp_size, 100);
  input ["ImportanceSamplingMax"  ] = Read(MultiArray::max_imp_size, 10000);
  input ["FailSamplingMax"        ] = Read(max_fail_num, 10);
  input ["TrajEnerDeviationMax"   ] = Read(MultiArray::ener_dev_max, -1.);
  input ["TrajAmomDeviationMax"   ] = Read(MultiArray::amom_dev_max, -1.);
  input ["FaceStatTolerance"      ] = Read(MultiArray::face_rel_tol, 0.1);
  input ["SpecStatTolerance"      ] = Read(MultiArray::spec_rel_tol, 0.1);
  input ["SpecReacTolerance"      ] = Read(MultiArray::reac_rel_tol, 0.1);
  input ["ReacTranTolerance"      ] = Read(MultiArray::tran_rel_tol, 0.1);
  input ["TrajTolerance"          ] = Read(trt, 1.e-5);
  input ["AngLengthTolerance"     ] = Read(Dynamic::Coordinates::ang_len_tol, 2.);
  input ["TrajTimeStep[au]"       ] = Read(Trajectory::default_step, 100.);
  input ["RandomPotErrorFlag"     ] = Read(rand_pot_err_flag, 0);
  input ["RadialEnergyFlag"       ] = Read(raden_flag, 0);
  input ["RadialEnergyFile"       ] = Read(raden_file, "raden.out");
  input ["AmomProjectFlag"        ] = Read(amproj_flag, 0);
  input ["AmomProjectFile"        ] = Read(amproj_file, "amproj.out");
  input ["RandomSeed"             ] = Read(_seed, 0);
  input ["MinAtomDist[Bohr]"      ] = Read(Dynamic::Coordinates::min_atom_dist, 1.5);
  input ["ExpPowerMax"            ] = Read(exp_max, 200.);

  std::string token, key, comment, name;

  // maximal parameter name size
  //
  itemp = 0;
    
  for(idit = input.begin(); idit != input.end(); ++idit)
    //
    if(idit->first.size() > itemp)
      //
      itemp = idit->first.size();

  const int name_len_max = itemp;

  IO::log << IO::log_offset << "taken as a comment:\n";
  
  itemp = 0;

  while(from >> token) {
    //
    if(token == IO::end_key())
      //
      break;

    idit = input.find(token);

    if(idit == input.end()) {
      //
      ++itemp;

      IO::log << IO::log_offset << "   " << token << "\n";

      std::getline(from, comment);
    }
    else
      //
      from >> idit->second;
  }

  if(!itemp)
    //
    IO::log << IO::log_offset << "   " << "None\n";

  if(!from) {
    //
    std::cerr << funame << "input stream is corrupted\n";

    throw Error::Form();
  }

  // check if all parameters were initialized
  //
  btemp = false;

  for(idit = input.begin(); idit != input.end(); ++idit)
    //
    if(!idit->second.is_init()) {
      //
      std::cerr << funame << idit->first << " is not initialized\n";

      btemp = true;
    }

  if(btemp)
    //
    throw Error::Init();

  // default values
  //
  IO::log << IO::log_offset << "default parameters:\n" << std::left;

  for(idit = input.begin(); idit != input.end(); ++idit)
    //
    if(idit->second.is_default())
      //
      IO::log << IO::log_offset << "   " << std::setw(name_len_max) << idit->first << " = " << idit->second << "\n";

  IO::log << std::right;
  
  // maximal exponent power
  //
  if(exp_max <= 0.) {
    //
    std::cerr << funame << "ExpPowerMax out of range: " << exp_max << "\n";

    throw Error::Range();
  }

  Limits::set_exp_pow_max(exp_max);

  // job type
  //
  if(job_type == "dynamical") {
    //
    _job = DYN_JOB;
  }
  else if(job_type == "statistical") {
    //
    _job = STAT_JOB;
  }
  else if(job_type == "test") {
    //
    _job = TEST_JOB;
  }
  else {
    //
    std::cerr << funame << "unknown job type: " << job_type << "; possible jobs: dynamical, statistical, and test\n";

    throw Error::Init();
  }
    
  // calculation mode
  //
  if(calc_mode == "canonical") {
    //
    _mode = T_MODE;
  }
  else if(calc_mode == "microcanonical") {
    //
    _mode = E_MODE;
  }
  else if(calc_mode == "e,j-resolved") {
    //
    _mode = J_MODE;
  }
  else {
    //
    std::cerr << funame << "unknown caculation mode: " << calc_mode << "; possible modes: canonical, microcanonical, or e,j-resolved\n";

    throw Error::Init();
  }

  switch(mode()) {
    //
  case T_MODE:
    //
    if(temperature() <= 0.) {
      //
      std::cerr << funame << "Temperature, " << temperature() / Phys_const::kelv << " K, should be positive\n";

      throw Error::Init();
    }

    break;
    //
  case E_MODE:
    //
    if(energy_value() < -.01) {
      //
      std::cerr << funame << "Energy, " << energy_value() / Phys_const::kcal << " kcal/mol, is out of range\n";

      throw Error::Init();
    }

    break;
    //
  case J_MODE:
    //
    if(energy_value() <= -1.) {
      //
      std::cerr << funame << "Energy, " << energy_value() / Phys_const::kcal << " kcal/mol, is out of range\n";

      throw Error::Init();
    }

    if(amom_value() < 0.) {
      //
      std::cerr << funame << "Angular momentum, " << amom_value() << ", should not be negative\n";

      throw Error::Init();
    }
  }

  // trajectories tolerances
  //
  Trajectory::traj_rel_tol.resize(Structure::dv_size(), trt);
  Trajectory::traj_abs_tol.resize(Structure::dv_size(), trt);

  // scaling factor
  //
  double scale;

  // using the temperature as a scaling factor
  //
  if(mode() == T_MODE) {
    //
    //
    scale = std::sqrt(temperature());
  }
  // using total energy as a scaling factor
  //
  else if(energy_value() < 0.) {
    //
    scale = std::sqrt(-energy_value());
  }
  else
    //
    scale = std::sqrt(energy_value());

  for(int i = 0; i < 3; ++i) {
    //
    itemp = Structure::pos_size() + Structure::orb_vel() + i;

    Trajectory::traj_abs_tol[itemp] *= scale / Structure::mass_sqrt();
  }

  for(int frag = 0; frag < 2; ++frag)
    //
    for(int i = 0; i < Structure::fragment(frag).vel_size(); ++i) {
      //
      itemp = Structure::pos_size() + Structure::ang_vel(frag) + i;

      Trajectory::traj_abs_tol[itemp] *= scale / Structure::fragment(frag).imom_sqrt(i);
    }

  // normalization factor
  //
  _norm_factor = std::pow(Structure::mass(), 1.5) * 2. / M_PI / M_PI;

  for(int frag = 0; frag < 2; ++frag)
    //
    switch (Structure::fragment(frag).type()) {
      //
    case Molecule::MONOATOMIC:
      //
      _norm_factor *= M_PI;

      break;
      //
    case Molecule::LINEAR:
      //
      _norm_factor *= 2. * M_PI * Structure::fragment(frag).imom(0);

      break;
      //
    case Molecule::NONLINEAR:
      //
      for(int i = 0; i < 3; ++i)
	//
	_norm_factor *= std::sqrt(2. * M_PI * Structure::fragment(frag).imom(i));
    }

  switch(mode()) {
    //
  case T_MODE:
    //
    _norm_factor *= std::pow(temperature(), double(Structure::tm_dof() - 1) / 2.);

    break;
    //
  case E_MODE:
    //
    _norm_factor /= tgamma(double(Structure::tm_dof() + 1) / 2.);

    break;
    //
  case J_MODE:
    //
    _norm_factor *= std::sqrt(2. / M_PI) / tgamma(double(Structure::tm_dof() - 2) / 2.) * amom_value() * amom_value();
  }
}

/**************************************************************************
 ************************ DYNAMIC SAMPLING METHODS ************************
 **************************************************************************/

void CrossRate::DynSmp::init (const Potential::Wrap& pot, const DivSur::MultiSur& ms, int prim, const Dynamic::Coordinates& dc)
{
  const char funame [] = "CrossRate::DynSmp::init: ";
  
  double dtemp;

  if(back) {
    //
    //std::cerr << funame << "already initialized\n";

    throw InitErr();
  }

  // generalized coordinates
  //
  (Dynamic::Coordinates&)*this = dc;

  // potential energy
  //
  _energy = pot(dc);

  //  
  if(_energy < reactive_energy()) {
    //
    throw LowEnerErr();
  }

  // random variable
  //
  _ranval = Random::flat();

  // statistical weight and random velocities
  //
  switch(mode()) {
    //
  case T_MODE: // canonical calculation
    //
    _weight = ms.t_velocity(prim, temperature(), *this);

    dtemp = -potential_energy() / temperature();
    
    // high energy
    //
    if(dtemp < -Limits::exp_pow_max()) {
      //
      _weight = -1.;
    }
    else
      //
      _weight *= std::exp(dtemp);

    break;
    //
  case E_MODE: // microcanonical calculation
    //
    dtemp = energy_value() - potential_energy();
    
    if(dtemp > 0.) {
      //
      _weight = ms.e_velocity(prim, dtemp, *this);
    }
    else
      //
      _weight = -1.;
    
    break;
    //
  case J_MODE: // E,J-resolved calculation
    //
    dtemp = energy_value() - potential_energy();
    
    if(dtemp > 0.) {
      //
      _weight = ms.j_velocity(prim, dtemp, amom_value(), *this);
    }
    else
      //
      _weight = -1.;
  }

  back.init(new DynRes(pot, *this, BACKWARD));
  forw.init(new DynRes(pot, *this,  FORWARD));
}

// propagate trajectory both forward and backward
//
void CrossRate::DynSmp::run_traj (const DivSur::MultiSur& ms, const DivSur::face_t& face, const Dynamic::CCP& stop)
{
  const char funame [] = "CrossRate::DynSmp::run_traj: ";

  if(!isinit()) {
    //
    std::cerr << funame << "crossrate environment has not yet been initialized\n";

    throw Error::Init();
  }

  if(is_run()) {
    //
    std::cerr << funame << "trajectory has been run already\n";

    throw Error::Init();
  }

  int    itemp;
  double dtemp;

  // trajectory is directed away from the reactant
  //
  if(reactant() == face.first){
    //
    try {
      //
      back->run(stop, ms, Dynamic::CCP(new SpecCondition(ms, face.first)));

      if(back->species() == face.first) {
	//
	back->stat = DynRes::DIRECT;
      }
      else {
	//
	back->stat = DynRes::RECROSS;

	return;
      }
    }
    catch(Trajectory::PotentialFailure) {
      //
      back->stat =  DynRes::POT_FAIL;

      return;
    }
    catch(Trajectory::RunFailure) {
      //
      back->stat =  DynRes::RUN_FAIL;

      return;
    }

    try {
      //
      forw->run(stop, ms);
	
      forw->stat = DynRes::PASS;

      return;
    }
    catch(Trajectory::PotentialFailure) {
      //
      forw->stat =  DynRes::POT_FAIL;

      return;
    }
    catch(Trajectory::RunFailure) {
      //
      forw->stat =  DynRes::RUN_FAIL;

      return;
    }
  }

  // trajectory is in the direction of reactant
  //
  if(reactant() == face.second){
    //
    try {
      //
      forw->run(stop, ms, Dynamic::CCP(new SpecCondition(ms, face.second)));

      if(forw->species() == face.second) {
	//
	forw->stat = DynRes::DIRECT;
      }
      else {
	//
	forw->stat = DynRes::RECROSS;

	return;
      }
    }
    catch(Trajectory::PotentialFailure) {
      //
      forw->stat =  DynRes::POT_FAIL;

      return;
    }
    catch(Trajectory::RunFailure) {
      //
      forw->stat =  DynRes::RUN_FAIL;

      return;
    }

    try {
      //
      back->run(stop, ms);

      back->stat = DynRes::PASS;

      return;
    }
    catch(Trajectory::PotentialFailure) {
      //
      back->stat =  DynRes::POT_FAIL;

      return;
    }
    catch(Trajectory::RunFailure) {
      //
      back->stat =  DynRes::RUN_FAIL;

      return;
    }
  }

  try {
    //
    back->run(stop, ms, Dynamic::CCP(new SpecCondition(ms, face.first)));

    if(back->species() == face.first) {
      //
      back->stat = DynRes::DIRECT;
    }
    else
      //
      back->stat = DynRes::RECROSS;
  }
  catch(Trajectory::PotentialFailure) {
    //
    back->stat =  DynRes::POT_FAIL;

    return;
  }
  catch(Trajectory::RunFailure) {
    //
    back->stat =  DynRes::RUN_FAIL;

    return;
  }

  try {
    //
    forw->run(stop, ms, Dynamic::CCP(new SpecCondition(ms, face.second)));

    if(forw->species() == face.second) {
      //
      forw->stat = DynRes::DIRECT;
    }
    else
      //
      forw->stat = DynRes::RECROSS;
  }
  catch(Trajectory::PotentialFailure) {
    //
    forw->stat =  DynRes::POT_FAIL;

    return;
  }
  catch(Trajectory::RunFailure) {
    //
    forw->stat =  DynRes::RUN_FAIL;

    return;
  }

  // either forward and backward trajectories both recrossed or both finished
  //
  if(forw->stat == back->stat)
    //
    return;

  // run backward to the end
  //
  if(back->stat == DynRes::RECROSS) {
    //
    try {
      //
      back->run(stop, ms);
    }
    catch(Trajectory::PotentialFailure) {
      //
      back->stat =  DynRes::POT_FAIL;

      return;
    }
    catch(Trajectory::RunFailure) {
      //
      back->stat =  DynRes::RUN_FAIL;

      return;
    }
  }

  // run forward to the end
  //
  if(forw->stat == DynRes::RECROSS) {
    //
    try {
      //
      forw->run(stop, ms);
    }
    catch(Trajectory::PotentialFailure) {
      //
      forw->stat =  DynRes::POT_FAIL;

      return;
    }
    catch(Trajectory::RunFailure) {
      //
      forw->stat =  DynRes::RUN_FAIL;

      return;
    }
  }
}

void CrossRate::DynSmp::print_traj_results () const
{
  int old_precision = IO::log.precision(3);
  
  IO::log << std::setw(20) << "time/a.u."
	  << std::setw(13) << "dist/bohr"
	  << std::setw(13) << "ener/kcal"
	  << std::setw(13) << "amom/a.u."
	  << std::setw(13) << "amom/a.u."
	  << std::setw(13) << "amom/a.u."

	  << "\n";

  IO::log << std::setw(20) << "0"
	  << std::setw(13) << orb_len()
	  << std::setw(13) << total_energy() / Phys_const::kcal;

  D3::Vector tam;
      
  total_angular_momentum(tam);
      
  for(int i = 0; i < 3; ++i)
    //
    IO::log << std::setw(13) << tam[i];
      
  IO::log << "\n";

  if(forw->stat != DynRes::INIT) {
    //
    IO::log << std::setw(20) << forw->time()
	    << std::setw(13) << forw->orb_len() 
	    << std::setw(13) << forw->final_energy() / Phys_const::kcal;

    forw->total_angular_momentum(tam);

    for(int i = 0; i < 3; ++i)
      //
      IO::log << std::setw(13) << tam[i];
	
    IO::log << "\n";
  }

  if(back->stat != DynRes::INIT) {
    //
    IO::log << std::setw(20) << back->time()
	    << std::setw(13) << back->orb_len() 
	    << std::setw(13) << back->final_energy() / Phys_const::kcal;

    back->total_angular_momentum(tam);
	
    for(int i = 0; i < 3; ++i)
      //
      IO::log << std::setw(13) << tam[i];
	
    IO::log << "\n";
  }

  IO::log.precision(old_precision);
}

/************************************************************************
 ************************* FACET ARRAY METHODS **************************
 ************************************************************************/

int CrossRate::FacetArray::init_traj_num () const
{
  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(!fit->is_run())
      //
      ++res;
  
  return res;
}

int CrossRate::FacetArray::potfail_traj_num () const
{
  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(fit->is_pot_fail())
      //
      ++res;
  
  return res;
}

int CrossRate::FacetArray::runfail_traj_num () const
{
  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(fit->is_run_fail())
      //
      ++res;
  
  return res;
}

int CrossRate::FacetArray::exclude_traj_num () const
{
  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(fit->is_exclude(0) || fit->is_exclude(1))
      //
      ++res;
  
  return res;
}

int CrossRate::FacetArray::run_traj_num () const
{
  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(fit->is_run() && !fit->is_pot_fail() && !fit->is_run_fail())
      //
      ++res;

  return res;
}

int CrossRate::FacetArray::direct_traj_num () const
{
  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(!fit->is_pot_fail() && !fit->is_run_fail() &&
	//
       fit->forw->stat == DynRes::DIRECT && fit->back->stat == DynRes::DIRECT)
      //
      ++res;

  return res;
}

int CrossRate::FacetArray::recross_traj_num (int dir) const
{
  const char funame [] = "CrossRate::FacetArray::recross_traj_num: ";

  int res = 0;
  
  for(const_iterator fit = begin(); fit != end(); ++fit)
    //
    if(!fit->is_pot_fail() && !fit->is_run_fail() && 
	//
	fit->dyn_res(dir)->stat == DynRes::RECROSS)
      //
      ++res;
      
  return res;
}

int CrossRate::FacetArray::direct_traj_num (int dir, int* exclude) const
{
  const char funame [] = "CrossRate::FacetArray::direct_traj_num: ";

  int res = 0;
  
  if(exclude)
    //
    *exclude = 0;

  for(const_iterator fit = begin(); fit != end(); ++fit) {
    //
    if(!fit->is_pot_fail() && !fit->is_run_fail() &&
      //
      fit->dyn_res(dir)->stat == DynRes::DIRECT) {
      //
      ++res;

      if(exclude && fit->is_exclude(dir))
	//
	(*exclude)++;
    }
  }

  return res;
}

int CrossRate::FacetArray::reac_traj_num (int dir, int spec, int* back, int* forw) const
{
  const char funame [] = "CrossRate::FacetArray::reac_traj_num: ";

  int res = 0;
  
  if(forw)
    //
    *forw = 0;

  if(back)
    //
    *back = 0;


  for(const_iterator fit = begin(); fit != end(); ++fit) {
    //
    if(!fit->is_pot_fail() && !fit->is_run_fail() &&
       //
       fit->dyn_res(!dir)->stat == DynRes::DIRECT && 
       //
       fit->dyn_res(dir)->species() == spec) {
	 //
	 ++res;
      
	if(forw && fit->is_exclude(dir))
	  //
	  (*forw)++;

	if(back && fit->is_exclude(!dir))
	  //
	  (*back)++;
      }
  }

  return res;
}

// propagate trajectories for facet samplings
//
void CrossRate::FacetArray::run_traj (const DivSur::MultiSur& ms, const DivSur::face_t& face, Dynamic::CCP stop)
{
  const char funame [] = "CrossRate::FacetArray::run_traj: ";
  
  if(reactant() >= 0 && face.first != reactant() && face.second != reactant())
    //
    return;

  if(!size())
    //
    return;

  const int count_max = init_traj_num();
  
  if(!count_max)
    //
    return;

  int new_share, old_share = 0;

  int count = 0;

  // run trajectory cycle
  //

#ifdef OPENMP

#pragma omp parallel for default(shared) schedule(static)

  for(int traj = 0; traj < size(); ++traj) {
    //
    iterator fit = begin();

    std::advance(fit, traj);

#else

  for(iterator fit = begin(); fit != end(); ++fit) {

#endif

    if(fit->is_run())
      //
      continue;
    
    fit->run_traj(ms, face, stop);
    
    ++count;

#ifndef OPENMP
#ifdef  DEBUG

    if(fit->is_run_fail() || fit->is_pot_fail()) {
      //
      IO::log << "   failed\n\n";

      continue;
    }

    IO::log << count << "-th trajectory:\n";
      
    fit->print_traj_results();

    IO::log << "\n";

#else

    new_share =(int)((double)count / (double)count_max * 100.);
    print_progress(old_share, new_share);

#endif
#endif

  }// run trajectory cycle
  
  // checking
  //
  for(const_iterator fit = begin(); fit != end(); ++fit) {
    //
    if(fit->is_run_fail() || fit->is_pot_fail() || !fit->is_run())
      //
      continue;

    if(fit->back->stat == DynRes::DIRECT && fit->back->species() != face.first) {
      //
      std::cerr << "WARNING: backward direct trajectory finished in the wrong well\n";
    }
	    
    if(fit->forw->stat == DynRes::DIRECT && fit->forw->species() != face.second) {
      //
      std::cerr << "WARNING: forward direct trajectory finished in the wrong well\n";
    }    
  }
}

void CrossRate::FacetArray::merge (const FacetArray& a)
{
  _flux_num += a._flux_num;
  _fail_num += a._fail_num;
  _fake_num += a._fake_num;

  _flux += a._flux;
  _fvar += a._fvar;

  if(a._min_ener < _min_ener) {
    //
    _min_ener = a._min_ener;
    _min_geom = a._min_geom;
  }

  if(a._max_weight > _max_weight) {
    //
    _max_weight = a._max_weight;

    // remove all non-quallifying samplings from the importance samplings list
    //
    iterator it = begin();

    while(it != end()) {
      //
      if(it->weight() < it->ranval() * _max_weight) {
	//
	it = erase(it);
      }
      else
	//
	++it;
    }
    
    for(const_iterator cit = a.begin(); cit != a.end(); ++cit)
      //
      push_back(*cit);
  }
  else
    //
    for(const_iterator cit = a.begin(); cit != a.end(); ++cit)
      //
      if(cit->weight() > cit->ranval() * _max_weight)
	//
	push_back(*cit);
} 
 
void CrossRate::FacetArray::add_smp (const DynSmp& smp)
{
  const char funame [] = "CrossRate::FacetArray::add_smp: ";

  // update minimal energy
  //
  if(!flux_num() || smp.potential_energy() < _min_ener) {
    //
    _min_ener = smp.potential_energy();
    _min_geom = smp;
  }

  ++_flux_num;

  // zero weight sampling
  //
  if(smp.weight() <= 0.)
    //
    return;
    
  // update flux
  //
  _flux += smp.weight();
  _fvar += smp.weight() * smp.weight();

  // update importance samplings data
  //
  if(!size()) {
    //
    _max_weight = smp.weight();
    push_back(smp);

    return;
  }

  if(smp.weight() > _max_weight) {
    //
    // use new sampling weight as a reference weight
    //
    _max_weight = smp.weight();

    // remove all non-quallifying samplings from the importance samplings list
    //
    iterator it = begin();

    while(it != end())
      //
      if(it->weight() < it->ranval() * _max_weight) {
	//
	it = erase(it);
      }
      else
	//
	++it;
    
    // add sampling to the importance samplings list
    //
    push_back(smp);
  }
  else if(smp.weight() > smp.ranval() * _max_weight)
    //
    push_back(smp);

  return;
}

long CrossRate::FacetArray::samp_num () const 
{ 
  if(rand_pot_err_flag) {
    //
    return  _flux_num;
  }
  else
    //
    return  _flux_num + _fail_num;
}

double CrossRate::FacetArray::flux_value () const
{
  if(!samp_num()) 
    //
    return 0.;
  
  return _flux / double(samp_num());
}

double CrossRate::FacetArray::flux_var () const
{
  const char funame [] = "CrossRate::FacetArray::flux_var: ";

  if(!samp_num())
    //
    return 0.;

  double dtemp = flux_value();

  double res   = _fvar / double(samp_num()) - dtemp * dtemp;

  if(res < 0.) {
    //
    std::cerr << funame << "WARNING: negative variance\n";

    return 0.;
  }

  return res;
}
  
double CrossRate::FacetArray::flux_rel_var () const
{
  const char funame [] = "CrossRate::FacetArray::flux_rel_var: ";

  if(!samp_num() || _flux == 0.)
    //
    return 0.;
  
  double res = _fvar / _flux / _flux - 1. / double(samp_num());

  if(res <= 0.) {
    //
    IO::log << funame << "WARNING: negative variance\n";

    return 0.;
  }

  return res;
}

/**************************************************************************************
 *************************************** SurArray *************************************
 **************************************************************************************/

int CrossRate::SurArray::merge (const SurArray& a)
{
  int res = 0;

  _fail    += a._fail;
  _inner   += a._inner;
  _exclude += a._exclude;
  _close   += a._close;


  for(const_iterator sit = a.begin(); sit != a.end(); ++sit) {
    //
    if(find(sit->first) == end()) {
      //
      (*this)[sit->first] = sit->second;

      res = 1;
    }
    else
      //
      (*this)[sit->first].merge(sit->second);
  }

  return res;
}


/****************************************************************************************
 *************************************** MultiArray *************************************
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

void CrossRate::MultiArray::run_traj ()
{
  for(iterator mit = begin(); mit != end(); ++mit) {
    //
    bool btemp = true;

    for(SurArray::iterator sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      // check that minimal potential enery is bigger than reactive energy
      //
      if(sit->second.size() && sit->second.min_ener() < reactive_energy()) {
	//
	std::cerr << mit - begin() << "-th surface,  " << sit->first << " facet: " 
		  << "minimal energy[kcal/mol], "<< sit->second.min_ener() /Phys_const::kcal
		  << ", is smaller than the reactive energy: "
		  << "decrease the reactive energy or change the transition state dividing surface\n";

	IO::log.flush();

	throw Error::Run();
      }

      if(sit->second.init_traj_num()) {
	//
	if(btemp) {
	  //
	  IO::log << "   " << mit - begin() << "-th surface:\n";

	  btemp = false;
        }


#ifdef OPENMP

	IO::log <<  "      " << sit->first << " facet: ... ";

#else

	IO::log <<  "      " << sit->first << " facet:\n";

#endif

	sit->second.run_traj(_ms, sit->first, _stop);

#ifdef OPENMP

	IO::log << "done\n";
#endif

      }//
      //
    }//
    //
  }//
}//

int CrossRate::MultiArray::_sample (iterator mit, const std::set<DivSur::face_t >& face_work)
{
  const char funame [] = "CrossRate::MultiArray::_sample: ";

  int    itemp;
  double dtemp;

  std::set<int> ierr;

  const int sur = mit - begin();

  std::vector<SmpRes> smp_res(smp_set_size);

#pragma omp parallel for default(shared) schedule(static)
  //
  for(int i = 0; i < smp_set_size; ++i) {
    //
    Dynamic::Coordinates dc;

    _ms.random_orient(sur, dc);
  
    (DivSur::SmpRes&)smp_res[i] = _ms.surface_test(sur, dc);

    if(smp_res[i].stat != DivSur::FACET)
      //
      continue;

    // facet is not of interest
    //
    if(reactant() >= 0 && smp_res[i].face.first != reactant() && smp_res[i].face.second != reactant()) {
      //
      smp_res[i].smp_stat = SmpRes::FAKE;

      continue;
    }

    // no potential calculation is needed
    //
    if(mit->find(smp_res[i].face) != mit->end() && face_work.find(smp_res[i].face) == face_work.end()) {
      //
      smp_res[i].smp_stat = SmpRes::FAKE;

      continue;
    }

    try {
      //
      smp_res[i].init(_pot, _ms, sur, dc);

      smp_res[i].smp_stat = SmpRes::SUCCESS;
    }
    catch(Potential::Exception) {
      //
      smp_res[i].smp_stat = SmpRes::FAIL;
    }
    catch(DynSmp::InitErr) {
    //
      smp_res[i].error = DynSmp::INIT_ERR;
    }
    catch(DynSmp::LowEnerErr) {
      //
      smp_res[i].error = DynSmp::LOW_ENER_ERR;
    }
  }

  itemp = 0;
  dtemp = 0.;

  for(int i = 0; i < smp_set_size; ++i) {
    //
    itemp |= smp_res[i].error;

    if(smp_res[i].error == DynSmp::LOW_ENER_ERR && smp_res[i].potential_energy() < dtemp)
      //
      dtemp = smp_res[i].potential_energy();
  }

  if(itemp & DynSmp::INIT_ERR)
    //
    std::cerr << sur << "-th surface: some samplings have been already initialized\n";

  if(itemp & DynSmp::LOW_ENER_ERR)
    //
    std::cerr << sur << "-th surface: sampling potential energy, " << dtemp / Phys_const::kcal 
	//
	      << " kcal/mol, is lower than the reactive energy, " << reactive_energy() / Phys_const::kcal << "\n";

  if(itemp) {
    //
    IO::log.flush();

    throw Error::Run();
  }

  int isnew = 0;
  
  for(std::vector<SmpRes>::const_iterator rit = smp_res.begin(); rit != smp_res.end(); ++rit) {
    //
    switch(rit->stat) {
      //
    case DivSur::INNER:
      //
      mit->add_inner();

      continue;
      //
    case DivSur::CLOSE:
      //
      mit->add_close();

      continue;
      //
    case DivSur::EXCLUDE:
      //
      mit->add_excl();

      continue;
      //
    case DivSur::FAIL:
      //
      mit->add_fail();

      continue;
    }

    // new facet
    //
    if(mit->find(rit->face) == mit->end()) {
      //
      IO::log << "      new " << rit->face << " facet\n";
      
      isnew = 1;
    }

    switch(rit->smp_stat) {
      //
    case SmpRes::FAKE:
      //
      (*mit)[rit->face].add_fake();

      continue;
      //
    case SmpRes::FAIL:
      //
      (*mit)[rit->face].add_fail();

      continue;
    }

    (*mit)[rit->face].add_smp(*rit);
  }
  
  return isnew;
}

void CrossRate::MultiArray::_get_stat_flux(std::vector<double>& flux_mean, std::vector<double>& flux_rmsd, std::vector<double>& flux_variance) const
{
  flux_mean.clear();
  flux_mean.resize(_ms.species_size());
  flux_rmsd.clear();
  flux_rmsd.resize(_ms.species_size());
  flux_variance.clear();
  flux_variance.resize(_ms.species_size());

  double face_flux, face_rmsd, face_var;

  for(const_iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    //
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit) {// facet cycle
      //
      if(reactant() >= 0 && reactant() != sit->first.first && reactant() != sit->first.second)
	//
	continue;

      if(sit->second.size()) {
	//
	face_flux = mit->face_flux(sit);
	face_rmsd = mit->vol_frac(sit) * std::sqrt(sit->second.flux_var());
	face_var  = face_flux * face_flux * (sit->second.flux_rel_var() + mit->vol_rel_var(sit));

	if(reactant() >= 0) {
	  //
	  flux_mean[reactant()]      += face_flux;
	  flux_rmsd[reactant()]      += face_rmsd;
	  flux_variance[reactant()]  += face_var;
	}
	else {
	  //
	  flux_mean[sit->first.first]      += face_flux;
	  flux_rmsd[sit->first.first]      += face_rmsd;
	  flux_variance[sit->first.first]  += face_var;

	  flux_mean[sit->first.second]      += face_flux;
	  flux_rmsd[sit->first.second]      += face_rmsd;
	  flux_variance[sit->first.second]  += face_var;
	}//
	//
      }//
      //
    }// facet cycle
    //
  }// surface cycle
  //
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

  // primitive surface cycle
  //
  for(const_iterator mit = begin(); mit != end(); ++mit) {
    //
    // facet cycle
    //
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      if(!sit->second.run_traj_num())
	//
	continue;

      double face_flux = mit->face_flux(sit);

      // ward cycle
      //
      for (int ward = 0; ward < 2; ++ward) {
	//
	int react;
	
	switch(ward) {
	  //
	case BACKWARD:
	  //
	  react = sit->first.second;

	  break;
	  //
	case FORWARD:
	  //
	  react = sit->first.first;
	}

	if(reactant() >= 0 && react != reactant())
	  //
	  continue;

	//product cycle
	//
	for(int prod = 0; prod < _ms.species_size(); ++prod) {
	  //
	  if(prod != react && sit->second.reac_traj_num(ward, prod)) {
	    //
	    // dynamical correction factor
	    //
	    double dyn_fac = (double)sit->second.reac_traj_num(ward, prod) / (double)sit->second.run_traj_num();

	    double dyn_var = dyn_fac * (1. - dyn_fac);

	    DivSur::face_t tran(react, prod);

	    // reactive flux
	    //
	    flux_mean[tran] += dyn_fac * face_flux;

	    // reactive flux standard deviation
	    //
	    flux_rmsd[tran] += std::sqrt(dyn_var) * face_flux;

	    // reactive flux variance
	    //
	    flux_variance[tran] += dyn_var * face_flux * face_flux / (double)sit->second.run_traj_num();
	  }//
	  //
	}// product cycle
	//
      }// ward cycle
      //
    }// facet cycle 
    //
  }// surface cycle
} 

void CrossRate::MultiArray::_print_sampling_results () const
{
  IO::log << "sampling results:\n";

  int old_prec = IO::log.precision(3);

  // surface cycle
  //
  for(const_iterator cmit = begin(); cmit != end(); ++cmit) {
    //
    if(!cmit->size())
      //
      continue;
    
    IO::log << cmit - begin() << "-th surface:\n";
    
    if(cmit->tot_smp_num())
      //
      IO::log << "     number of all samplings:             " << cmit->tot_smp_num()  << "\n";
    
    if(cmit->fail_num())
      //
      IO::log << "     number of failed samplings:          " << cmit->fail_num()  << "\n";
    
    if(cmit->inner_num())
      //
      IO::log << "     number of inner region samplings:    " << cmit->inner_num() << "\n";
    
    if(cmit->excl_num())
      //
      IO::log << "     number of excluded region samplings: " << cmit->excl_num()  << "\n";
    
    if(cmit->close_num())
      //
      IO::log << "     number of close atoms samplings:     " << cmit->close_num() << "\n";
    
    for(SurArray::const_iterator cit = cmit->begin(); cit != cmit->end(); ++cit) {
      //
      IO::log << "     " << cit->first << " facet:\n"
	//
	      << "          statistical flux:               " << cmit->face_flux(cit) * norm_factor()
	//
	      <<  " (" << std::sqrt(cit->second.flux_rel_var() + cmit->vol_rel_var(cit)) * 100 << " %)\n"
	//
	      << "          number of importance samplings: " << cit->second.size() << "\n"
	//
	      << "          number of flux samplings:       " << cit->second.flux_num()
	//
	      << "\n";
      
      if(cit->second.fake_num())
	//
	IO::log	<< "          number of fake samplings:       " << cit->second.fake_num() << "\n";
      
      if(cit->second.fail_num())
	//
	IO::log << "          number of failed samplings:     " << cit->second.fail_num() << "\n";
      
      IO::log   << "          minimal energy configuration:\n"
	//
		<< "             energy, kcal/mol:            "
	//
		<< cit->second.min_ener() / Phys_const::kcal << "\n";
      
      cit->second.min_geom().print_geom(IO::log, "             ");
    }
    IO::log << std::endl;
  }

  std::vector<double> flux_mean, flux_rmsd, flux_variance;
  _get_stat_flux(flux_mean, flux_rmsd, flux_variance);

  // output
  //
  IO::log << "outgoing statistical flux:\n";
  
  IO::log << std::setw(10) << "reactant" 
	  << std::setw(10) << "flux"
	  << std::setw(10) << "err, %"
	  << "\n";
  
  for(int spec = 0; spec < _ms.species_size(); ++spec) {
    //
    if(reactant() >= 0 && spec != reactant())
      //
      continue;

    // statistical flux
    //
    IO::log << std::setw(10) << spec
      //
	    << std::setw(10) << flux_mean[spec] * norm_factor();

    // flux rms
    //
    if(flux_mean[spec] == 0.) {
      //
      IO::log << std::setw(10) << "***";
    }
    else
      //
      IO::log << std::setw(10) << 100. * std::sqrt(flux_variance[spec]) / flux_mean[spec];
    
    IO::log << "\n";
  }
  IO::log << "\n";
  
  IO::log.precision(old_prec);
}

void CrossRate::MultiArray::_print_traj_results () const
{
  double dtemp;
  int    itemp;

  int reac, prod;

  double rms_dev, max_dev;
	
  int count, smp_index;

  double face_flux, dyn_fac;

  std::ostringstream oss;
	
  int old_precision = IO::log.precision(3);

  IO::log << "trajectory results:\n"; 
  
  for(const_iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    //
    if(!mit->size())
      //
      continue;
    
    IO::log << mit - begin() << "-th surface:\n";

    // facet cycle
    //
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      // no importance samplings
      //
      if(!sit->second.size())
	//
	continue;

      face_flux =  mit->face_flux(sit) * norm_factor();

      IO::log << "   " << sit->first <<  " facet:   " << sit->second.run_traj_num() << "\n";

      //IO::log << "      " << "number of trajectories run:               " << sit->second.run_traj_num() << "\n";

      IO::log << "      " << "number of direct trajectories:            " << sit->second.direct_traj_num() << "\n";
      
      //IO::log << "      " << "number of exlcuded region hits:           " << sit->second.exclude_traj_num() << "\n";

      if(sit->second.potfail_traj_num())
	//
	IO::log << "      " << "number of potential-failed trajectories:  " << sit->second.potfail_traj_num() << "\n";

      if(sit->second.runfail_traj_num())
	//
	IO::log << "      " << "number of integrator-failed trajectories: " << sit->second.runfail_traj_num() << "\n";

      if(sit->second.init_traj_num())
	//
	IO::log << "      " << "number of non-run trajectories:           " << sit->second.init_traj_num() << "\n";
      
      // trajectories with large energy or angular momentum deviations
      //
      std::set<std::pair<int, int> > big_dev_traj;
      
      //sampling cycle
      //
      count = 0;

      if(ener_dev_max > 0. || amom_dev_max > 0.) {
	//
	for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit, ++count) {
	  //
	  if(fit->is_run_fail() || fit->is_pot_fail())
	    //
	    continue;

	  for(int dir = 0; dir < 2; ++dir) {
	    //
	    if(fit->dyn_res(dir)->stat == DynRes::INIT)
	      //
	      continue;
	    
	    if(ener_dev_max > 0.) {
	      //
	      dtemp = fit->ener_deviation(dir);

	      if(dtemp > ener_dev_max)
		//
		big_dev_traj.insert(std::make_pair(dir, count));
	    }
	    
	    if(amom_dev_max > 0.) {
	      //
	      dtemp = fit->amom_deviation(dir);

	      if(dtemp > amom_dev_max)
		//
		big_dev_traj.insert(std::make_pair(dir, count));
	    }
	  }
	}
      }
      
      if(big_dev_traj.size()) {
	//
	IO::log << "      " << "number of large deviation trajectories:   " <<  big_dev_traj.size() << "      // ";


	if(ener_dev_max > 0.) {
	  //
	  IO::log << "energy deviation threshold = " << ener_dev_max / Phys_const::kcal << " kcal/mol";

	  if(amom_dev_max > 0.)
	    //
	    IO::log << ", ";
	}

	if(amom_dev_max > 0.)
	  //
	  IO::log << "angular momentum deviation threshold = " << amom_dev_max << " au";

	IO::log << "\n";
      }

      // print large deviation trajectories to the standard error
      //
      count = 0;

      if(big_dev_traj.size()) {
	//
	std::cerr << mit - begin() << "-th surface: " << sit->first << "-th facet: re-running " << big_dev_traj.size() << " troubled trajectories:\n";
	
	
	for(std::set<std::pair<int, int> >::const_iterator bit = big_dev_traj.begin(); bit != big_dev_traj.end(); ++bit, ++count) {
	  //
	  std::cerr << "running " << count << "-th trajectory: starts\n";
	
	  FacetArray::const_iterator fit = sit->second.begin();

	  std::advance(fit, bit->second);

	  int dir = bit->first;

	  Trajectory::Propagator prop(_pot, *fit, dir, std::cerr);

	  prop.run(_stop, _ms);

	  std::cerr << "running " << count << "-th trajectory: done\n\n";
	}
      }
      
      // no trajectories
      //
      if(!sit->second.run_traj_num()) 
	//
	continue;

      // direction cycle
      //
      for (int ward = 0; ward < 2; ++ward) {
	//
	switch(ward) {
	  //
	case BACKWARD:
	  //
	  reac = sit->first.second;

	  break;
	  //
	case FORWARD:
	  //
	  reac = sit->first.first;
	}

	if(reactant() >= 0 && reac != reactant())
	  //
	  continue;

	IO::log << "      " << reac << "-th reactant:\n";

	// recrossed trajectories
	//
	IO::log << "         " << "recrossed trajectories fraction (#): "
	  //
		<< std::setw(16) << (double)sit->second.recross_traj_num(!ward) / (double)sit->second.run_traj_num()
	  //
		<< " (" << sit->second.recross_traj_num(!ward) <<")\n";

	// backward direct trajectories & their excluded region hits part
	//
	//IO::log << "         " << "backward  trajectories / excluded region hits: "
	//
	//	<< std::setw(6) << sit->second.direct_traj_num(!ward, &itemp) << "/" << itemp << "\n";


	// reactive trajectories
	//
	IO::log << "         "   << "reactive  trajectories / excluded region hits:\n";
	IO::log << std::setw(15) << "R -> P" 
		<< std::setw(15) << "reactive flux"
		<< std::setw(15) << "fraction"
		<< std::setw(15) << "number"
		<< std::setw(15) << "error, %"
		<< "\n";

	for(prod = 0; prod < _ms.species_size(); ++prod) {
	  //
	  oss.str("");
	  oss << reac << " -> " << prod;
	  
	  dyn_fac = (double)sit->second.reac_traj_num(ward, prod) / (double)sit->second.run_traj_num();
	  
	  IO::log << std::setw(15)  << oss.str()
		  << std::setw(15) << dyn_fac * face_flux
		  << std::setw(15) << dyn_fac;

	  // reactive trajectories number & its excluded region hits number part
	  //
	  oss.str("");
	  
	  int back, forw;

	  oss << sit->second.reac_traj_num(ward, prod, &back, &forw) << "(" << back << "," << forw << ")";

	  IO::log << std::setw(15)  << oss.str();
	  
	  if(sit->second.reac_traj_num(ward, prod)) {
	    //
	    dtemp = 1./(double)sit->second.reac_traj_num(ward, prod) - 1./(double)sit->second.run_traj_num();
	    
	    dtemp = dtemp > 0. ? std::sqrt(dtemp) * 100. : 0.;
	    
	    IO::log << std::setw(15) << dtemp;
	  }
	  else
	    //
	    IO::log << std::setw(15) << "***";
	  
	  IO::log << "\n";
	  //
	}//
	
	IO::log << "\n";

	// energy conservation error output
	//
	IO::log << "         "   << "energy conservation error [kcal/mol]:\n"
		<< std::setw(15) << "R -> P" 
		<< std::setw(15) << "direction"
		<< std::setw(15) << "average"
		<< std::setw(15) << "max"
		<< std::setw(15) << "max traj. #"
		<< "\n";

	// product cycle
	//
	for(prod = 0; prod < _ms.species_size(); ++prod) {
	  //
	  for(int dir = 0; dir < 2; ++dir) {
	    //
	    rms_dev = max_dev = 0.;

	    count = smp_index = 0;

	    int max_dev_smp = -1;
	    
	    //sampling cycle
	    //
	    for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit, ++smp_index) {
	      //
	      if(!fit->is_run_fail() && !fit->is_pot_fail() &&
		 //
		 fit->dyn_res(!ward)->stat == DynRes::DIRECT &&
		 //
		 fit->dyn_res(ward)->species() == prod) {
		//
		++count;

		dtemp = fit->ener_deviation(1 - (ward + dir) % 2);

		rms_dev += dtemp * dtemp;

		if(dtemp > max_dev) {
		  //
		  max_dev = dtemp;

		  max_dev_smp = smp_index;
		}
		//
	      }//
	      //
	    }// sampling cycle

	    if(count) {
	      //
	      rms_dev /= (double)count;

	      rms_dev = std::sqrt(rms_dev);
	    }
	      
	    // output
	    //
	    oss.str("");
	  
	    oss << reac << " -> " << prod;

	    IO::log << std::setw(15) << oss.str();

	    switch(dir) {
	      //
	    case BACKWARD:
	      //
	      IO::log << std::setw(15) << "backward";

	      break;
	      //
	    case FORWARD:
	      //
	      IO::log << std::setw(15) << "forward";
	    }

	    IO::log << std::setw(15) << rms_dev / Phys_const::kcal
	      //
		    << std::setw(15) << max_dev / Phys_const::kcal
	      //
		    << std::setw(15) << max_dev_smp << "\n";
	    //
	  }// direction cycle
	  //
	}// product cycle

	IO::log << "\n";
	
	// angular momentum conservation error output
	//
	IO::log << "         "   << "angular momentum conservation error [a.u.]:\n"
		<< std::setw(15) << "R -> P" 
		<< std::setw(15) << "direction"
		<< std::setw(15) << "average"
		<< std::setw(15) << "max"
		<< std::setw(15) << "max traj. #"
		<< "\n";
	
	// product cycle
	//
	for(prod = 0; prod < _ms.species_size(); ++prod) {
	  //
	  for(int dir = 0; dir < 2; ++dir) {
	    //
	    rms_dev = max_dev = 0.;

	    count = smp_index = 0;

	    int max_dev_smp = -1;
	    
	    //sampling cycle
	    //
	    for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit, ++smp_index) {
	      //
	      if(!fit->is_run_fail() && !fit->is_pot_fail() &&
		 //
		 fit->dyn_res(!ward)->stat == DynRes::DIRECT &&
		 //
		 fit->dyn_res(ward)->species() == prod) {
		//
		++count;

		dtemp = fit->amom_deviation(1 - (ward + dir) % 2);

		rms_dev += dtemp * dtemp;

		if(dtemp > max_dev) {
		  //
		  max_dev = dtemp;

		  max_dev_smp = smp_index;
		}//
		//
	      }//
	      //
	    }// sampling cycle

	    // normalization
	    //
	    if(count) {
	      //
	      rms_dev /= (double)count;

	      rms_dev = std::sqrt(rms_dev);
	    }
	    
	    // output
	    //
	    oss.str("");
	  
	    oss << reac << " -> " << prod;

	    IO::log << std::setw(15) << oss.str();

	    switch(dir) {
	      //
	    case BACKWARD:
	      //
	      IO::log << std::setw(15) << "backward";

	      break;
	      //
	    case FORWARD:
	      //
	      IO::log << std::setw(15) << "forward";
	    }
	    
	    IO::log << std::setw(15) << rms_dev
	      //
		    << std::setw(15) << max_dev
	      //
		    << std::setw(15) << max_dev_smp << "\n";
	    
	  }// direction cycle
	  //
	}// product cycle

	IO::log << "\n";
	//
      }// ward cycle
      //
    }// facet cycle
    //
  }// surface cycle
  
  // crossrate output
  //
  std::map<DivSur::face_t, double> flux_mean, flux_rmsd, flux_variance;
  
  _get_reac_flux(flux_mean, flux_rmsd, flux_variance);
  
  IO::log << "transition reactive flux:\n"
    //
	  << std::setw(10) << "transition"
	  << std::setw(10) << "flux"
	  << std::setw(10)  << "err, %"
	  << "\n";
  
  for(std::map<DivSur::face_t, double>::const_iterator it = flux_mean.begin(); it != flux_mean.end(); ++it) {
    //
    if(reactant() >= 0 && it->first.first != reactant())
      //
      continue; 

    IO::log << std::setw(10) << it->first
      //
	    << std::setw(10) << it->second * norm_factor();

    if(it->second == 0.) {
      //
      IO::log<< std::setw(10) << "***";
    }
    else
      //
      IO::log<< std::setw(10) << 100. * std::sqrt(flux_variance[it->first]) / it->second;

    IO::log << "\n";
  }
  IO::log << "\n";

  std::vector<double> reac_flux(_ms.species_size(), 0.);
  std::vector<double> reac_var(_ms.species_size(), 0.);
  
  for(std::map<DivSur::face_t, double>::const_iterator it = flux_mean.begin(); it != flux_mean.end(); ++it) {
    //
    if(reactant() >= 0 && it->first.first != reactant())
      //
      continue;
    
    reac_flux[it->first.first] += it->second;
    reac_var[it->first.first]  += flux_variance[it->first];
  }

  IO::log << "outgoing reactive flux:\n"
    //
	  << std::setw(10) << "reactant"
	  << std::setw(10) << "flux"
	  << std::setw(10) << "err, %"
	  << "\n";
  
  for(int r =  0; r < _ms.species_size(); ++r) {
    //
    if(reactant() >= 0 && r != reactant())
      //
      continue; 

    dtemp = reac_flux[r];
    
    IO::log << std::setw(10) << r
      //
	    << std::setw(10) << dtemp * norm_factor();

    if( dtemp == 0.) {
      //
      IO::log << std::setw(10) << "***";
    }
    else
      //
      IO::log << std::setw(10) << 100. * std::sqrt(reac_var[r]) / dtemp;

    IO::log << "\n";
  }
  IO::log << "\n";

  IO::log.precision(old_precision);
}

void CrossRate::MultiArray::_analize () const
{
  double dtemp;
  int stat, reac, prod; 
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

		//sampling cycle
		//
		for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit) {
		  if(fit->is_run_fail() || fit->is_pot_fail())
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
				  << std::setw(10) << (s->*reac_p)->orb_len()
				  << std::setw(10) << fit->orb_len()
				  << std::setw(10) << (s->*prod_p)->orb_len()
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
	    //
	    // sampling cycle
	    //
	    for(FacetArray::const_iterator fit = sit->second.begin(); fit != sit->second.end(); ++fit) {
	      //
	      if(fit->is_run_fail() || fit->is_pot_fail())
		//
		continue;

	      if(fit->back->stat == DynRes::DIRECT || fit->forw->stat == DynRes::DIRECT) {
		//
		for(int ward = 0; ward < 2; ++ward) {
		  //
		  switch(ward) {
		    //
		  case BACKWARD:
		    //
		    stat  = fit->forw->stat;
		    prod  = fit->back->species();
		    
		    break;
		    //
		  case FORWARD:
		    //
		    stat  = fit->back->stat;
		    prod  = fit->forw->species();
		  }
		  
		  if(stat == DynRes::DIRECT) {
		    //
		    raden[ward][prod].insert(fit->radial_kinetic_energy());
		  }
		}
	      }
	    }// sampling cycle

	    // ward cycle
	    //
	    for (int ward = 0; ward < 2; ++ward) {
	      //
	      switch(ward) {
		//
	      case BACKWARD:
		//
		reac = sit->first.second;
		
		break;
		//
	      case FORWARD:
		//
		reac = sit->first.first;
	      }
	      
	      if(reactant() >= 0 && reac != reactant())
		//
		continue;

	      // product cycle
	      //
	      for(prod = 0; prod < _ms.species_size(); ++prod)
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

bool CrossRate::MultiArray::_work () 
{
  const char funame [] = "CrossRate::MultiArray::_work: ";

  double dtemp;
  long   itemp;
  bool   btemp;

  long count, count_max;
  int old_share, new_share;

  SurArray::iterator sit;

  std::map<DivSur::face_t, long> proj_samp_num; // projected number of samplings per facet
  std::map<DivSur::face_t, long>::iterator nit;

  /**********************************************************************
   ********** SATISFY  MINIMAL NON-ZERO FLUX SAMPLINGS NUMBER ***********
   *********************************************************************/

  IO::log << "sampling to satisfy minimal flux samplings number (flux samplings counted) ..." << std::endl;
  
  for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
    //
    count_max = 0;

    for(sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      itemp = min_pot_size - sit->second.samp_num();

      if(itemp > 0)
        //
	count_max += itemp;
    }

    if(!count_max)
      //
      continue;

    IO::log << "   " << mit - begin() << "-th surface:\n";

    old_share = 0;

    // sampling loop
    //
    while(1) {
      //
      std::set<DivSur::face_t > face_work;

      count = 0;

      for(sit = mit->begin(); sit != mit->end(); ++sit) {
	//
	itemp = min_pot_size - sit->second.samp_num();

	if(itemp > 0) { 
	  //
	  face_work.insert(sit->first);

	  count += itemp;
	}
      }

      new_share = int((double)(count_max - count) / (double)count_max * 100.);
      print_progress(old_share, new_share);

      // exit condition
      //
      if(!face_work.size())
	//
	break;

      if(_sample(mit, face_work))
	//
	return false;
      //
    }// sampling loop
    //
  }// surface cycle

  IO::log << "done\n\n";

  /**********************************************************************
   ********** SATISFY FACET STATISTICAL FLUX RELATIVE TOLERANCE *********
   *********************************************************************/

  IO::log << "sampling to satisfy facet statistical flux relative tolerance (flux samplings counted) ..." << std::endl;

  for(iterator mit = begin(); mit != end(); ++mit) {// primitive cycle
    //
    proj_samp_num.clear();

    for(sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      dtemp = face_rel_tol * sit->second.flux_value();

      if(dtemp != 0.)
	//
	proj_samp_num[sit->first] = std::floor(sit->second.flux_var() / dtemp / dtemp + 0.5);
    }

    if(proj_samp_num.size()) {
      //
      IO::log << "   " << mit - begin() << "-th surface:\n"
	      << std::setw(15) << "facet" 
	      << std::setw(15) << "projected" 
	      << std::setw(15) << "current"  
	      << "\n";

      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	//
	IO::log << std::setw(15) << nit->first 
		<< std::setw(15) << nit->second 
		<< std::setw(15) << mit->find(nit->first)->second.samp_num() 
		<< "\n";
    }

    count_max = 0;

    for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
      //
      sit = mit->find(nit->first);

      itemp =  nit->second - sit->second.samp_num();

      if(itemp > 0) {
	//
	count_max += itemp; 
      }
    }
    
    if(!count_max)
      //
      continue;

    old_share = 0;

    // sampling loop
    //
    while(1) {
      //
      std::set<DivSur::face_t > face_work;

      count = 0;

      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	//
	sit = mit->find(nit->first);

	itemp =  nit->second - sit->second.samp_num();

	if(itemp > 0) {
	  //
	  count += itemp; 

	  face_work.insert(nit->first);
	}
      }

      new_share = int((double)(count_max - count) / (double)count_max * 100.);

      print_progress(old_share, new_share);

      // exit condition
      //
      if(!face_work.size())
	//
	break;

      if(_sample(mit, face_work))
	//
	return false;
      //
    }// sampling loop
    //
  }// primitive cycle

  IO::log << "done\n\n";

  /**********************************************************************
   ******** SATISFY SPECIES STATISTICAL FLUX RELATIVE TOLERANCE *********
   *********************************************************************/

  IO::log << "sampling to satisfy outgoing statistical flux relative tolerance (flux samplings counted) ..." << std::endl;

  std::vector<double> stat_flux_mean, stat_flux_rmsd, stat_flux_variance;

  _get_stat_flux(stat_flux_mean, stat_flux_rmsd, stat_flux_variance);

#ifdef DEBUG
  /*
  IO::log << std::setw(7) << "Species"  << std::setw(15) << "stat_flux_mean" 
	//
 	  << std::setw(15) << "stat_flux_rmsd" << std::setw(15) << "stat_flux_var" << std::endl;

  for(int s = 0; s < _ms.species_size(); ++s)
    //
    IO::log << std::setw(7) << s << std::setw(15) << stat_flux_mean[s] << std::setw(15) << stat_flux_rmsd[s] 
	//
	    << std::setw(15) << stat_flux_variance[s] << "\n";
*/
#endif

  // surface cycle
  //
  for(iterator mit = begin(); mit != end(); ++mit) {
    //
    proj_samp_num.clear();

    for(sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      int react;

      for(int ward = 0; ward < 2; ++ward) {// ward cycle
	//
	switch(ward) {
	  //
	case BACKWARD:
	  //
	  react = sit->first.second;

	  break;
	  //
	case FORWARD:
	  //
	  react = sit->first.first;
	}

	if(reactant() >= 0 && reactant() != react) 
	  //
	  continue;

	// projected number of samplings
	//
	dtemp = spec_rel_tol * stat_flux_mean[react];

	if(dtemp != 0.) {
	  //
	  double face_rmsd = mit->vol_frac(sit) * std::sqrt(sit->second.flux_var());

	  itemp = std::floor(face_rmsd * stat_flux_rmsd[react] / dtemp / dtemp + 0.5);

	  nit = proj_samp_num.find(sit->first);

	  if(nit == proj_samp_num.end() || itemp > nit->second)
	    //
	    proj_samp_num[sit->first] = itemp;
	}
	//
      }// ward cycle
      //
    }// facet cycle
    
    if(proj_samp_num.size()) {
      //
      IO::log << "   " << mit - begin() << "-th surface:\n"
	      << std::setw(15) << "facet" 
	      << std::setw(15) << "projected" 
	      << std::setw(15) << "current"  
	      << "\n";

      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	//
	IO::log << std::setw(15) << nit->first
		<< std::setw(15) << nit->second 
		<< std::setw(15) << mit->find(nit->first)->second.samp_num() 
		<< "\n";
    }

    count_max = 0;

    for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
      //
      sit = mit->find(nit->first);

      itemp =  nit->second - sit->second.samp_num();

      if(itemp > 0) {
	//
	count_max += itemp; 
      }
    }

    if(!count_max)
      //
      continue;

    old_share = 0;

    // sampling loop
    //
    while(1) {
      //
      std::set<DivSur::face_t > face_work;

      count = 0;
      //
      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	//
	sit = mit->find(nit->first);

	itemp =  nit->second - sit->second.samp_num();

	if(itemp > 0) {
	  //
	  count += itemp; 

	  face_work.insert(nit->first);
	}
      }
      
      new_share = int((double)(count_max - count) / (double)count_max * 100.);
      print_progress(old_share, new_share);

      if(!face_work.size())
	//
	break;

      if(_sample(mit, face_work))
	//
	return false;
      //
    }// sampling loop
    //
  }// surface cycle

  IO::log << "done\n\n";

  if(job() != DYN_JOB) {
    //
    _stat_pot_hit = _pot.reset_count();
  }
  else {

    /****************************** DYNAMICAL CALCULATIONS ***************************/

    /*************************************************************
     ******* SATISFY MINIMAL IMPORTANCE SAMPLINGS NUMBER *********
     *************************************************************/

    IO::log << "sampling to satisfy minimal importance samplings number (importance samplings counted) ..." << std::endl;

    for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
      //
      IO::log << "   " << mit - begin() << "-th surface:\n"
	      << std::setw(15) << "facet" 
	      << std::setw(15) << "projected" 
	      << std::setw(15) << "current" << "\n";

      count_max = 0;
      for(sit = mit->begin(); sit != mit->end(); ++sit) {
	//
	IO::log << std::setw(15) << sit->first
		<< std::setw(15) << min_imp_size
		<< std::setw(15) << sit->second.size()
		<< "\n";

	itemp = min_imp_size - sit->second.size();

	if(sit->second.size() && itemp > 0)
	  //
	  count_max += itemp;
      }

      if(!count_max)
	//
	continue;
      
      old_share = 0;

      // sampling loop
      //
      while(1) {
	//
	// tests
	//
	std::set<DivSur::face_t > face_work;

	count = 0;

	for(sit = mit->begin(); sit != mit->end(); ++sit)  {
	  //
	  itemp = min_imp_size - sit->second.size();

	  if(sit->second.size() && itemp > 0) {
	    //
	    count += itemp;

	    face_work.insert(sit->first);
	  }
	}

	new_share = int((double)(count_max - count) / (double)count_max * 100.);
	print_progress(old_share, new_share);

	// exit condition
	//
	if(!face_work.size())
	  //
	  break;

	if(_sample(mit, face_work))
	  //
	  return false;
	//
      }// sampling loop
      //
    }// surface cycle

    IO::log << "done\n\n";
  
    /*************************************************************
     *********************** RUN TRAJECTORIES ********************
     *************************************************************/

    _stat_pot_hit = _pot.reset_count();

    IO::log << "running trajectories ...\n";
    
    run_traj();
    
    IO::log << "done\n\n";

    _traj_pot_hit = _pot.reset_count();

    /*************************************************************
     * SATISFY SPECIES OUTGOING REACTIVE FLUX RELATIVE TOLERANCE *
     *************************************************************/

    bool run_again = false;

    IO::log << "sampling to satisfy outgoing reactive flux relative tolerance (importance samplings counted) ..." << std::endl;

    std::map<DivSur::face_t, double> cross_flux_mean, cross_flux_rmsd, cross_flux_variance;

    _get_reac_flux(cross_flux_mean, cross_flux_rmsd, cross_flux_variance);

    std::vector<double> reac_flux(_ms.species_size(), 0.);
    std::vector<double> reac_rmsd(_ms.species_size(), 0.);

    for(std::map<DivSur::face_t, double>::const_iterator it = cross_flux_mean.begin(); it != cross_flux_mean.end(); ++it) {
      //
      reac_flux[it->first.first] += it->second;
      reac_rmsd[it->first.first] += cross_flux_rmsd.find(it->first)->second;
    }

    // surface cycle
    //
    for(iterator mit = begin(); mit != end(); ++mit) {
      //
      // estimate necessary number of importance samplings for each facet on a given surface
      //
      proj_samp_num.clear();
      int react;

      // facet cycle
      //
      for(sit = mit->begin(); sit != mit->end(); ++sit) {
	//
	if(sit->second.run_traj_num()) {
	  //
	  for(int ward = 0; ward < 2; ++ward) {
	    //
	    switch(ward) {
	      //
	    case BACKWARD:
	      //
	      react = sit->first.second;

	      break;
	      //
	    case FORWARD:
	      //
	      react = sit->first.first;
	    }
	    
	    if(reactant() >= 0 && reactant() != react)
	      //
	      continue;
	    
	    // product cycle
	    //
	    for(int prod = 0; prod < _ms.species_size(); ++prod) {
	      //
	      if(prod != react && (itemp = sit->second.reac_traj_num(ward, prod))) {
		//
		// recrossing factor
		//
		dtemp = (double)itemp / (double)sit->second.run_traj_num();

		// facet reactive flux standard deviation
		//
		double face_rmsd = std::sqrt(dtemp * (1. - dtemp)) * mit->face_flux(sit);
		
		dtemp = reac_rel_tol * reac_flux[react];

		dtemp = face_rmsd * reac_rmsd[react] / dtemp / dtemp * (double)sit->second.size() / (double)sit->second.run_traj_num();

		if(dtemp > (double)max_imp_size) {
		  //
		  IO::log << "\n   WARNING: " << mit - begin() << "-th surface: " << sit->first << " facet: " 
			//
			  << react << "-th reactant: requested importance samplings #, " << std::floor(dtemp) << ",\n"
			//
			  << std::setw(12) << " " << "exceeds maximal facet importance samplings #, " << max_imp_size 
			//
			  << ": truncating" << "\n\n";

		    itemp = max_imp_size;
		}
		else
		  //
		  itemp = (long)dtemp;

		nit = proj_samp_num.find(sit->first);

		if(nit == proj_samp_num.end() || itemp > nit->second)
		  //
		  proj_samp_num[sit->first] = itemp;
		//
	      }//
	      //
	    }// product cycle
	    //
	  }// ward cycle
	  //
	}//
	//
      }// facet cycle 
    
      if(proj_samp_num.size()) {
	//
	IO::log << "   " << mit - begin() << "-th surface:\n"
		<< std::setw(15) << "facet" 
		<< std::setw(15) << "projected" 
		<< std::setw(15) << "current" << "\n";

	for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	  //
	  IO::log << std::setw(15) << nit->first 
		  << std::setw(15) << nit->second 
		  << std::setw(15) << mit->find(nit->first)->second.size() 
		  << "\n";
      }

      count_max = 0;

      for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	//
	sit = mit->find(nit->first);

	itemp = nit->second - sit->second.size();

	if(itemp > 0) {
	  //
	  count_max += itemp;
	}
      }

      if(!count_max)
	//
	continue;

      old_share = 0;

      // sampling loop
      //
      while(1) {
	//
	std::set<DivSur::face_t > face_work;

	count = 0;

	for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	  //
	  sit = mit->find(nit->first);

	  itemp = nit->second - sit->second.size();

	  if(itemp > 0) {
	    //
	    count += itemp;

	    face_work.insert(nit->first);
	  }
	}

	new_share = int((double)(count_max - count) / (double)count_max * 100.);

	print_progress(old_share, new_share);

	if(!face_work.size()) {
	  //
	  break;
	}
	else
	  //
	  run_again = true;


	if(_sample(mit, face_work))
	  //
	  return false;
	//
      }// sampling loop
      //
    }// surface cycle

    IO::log << "done\n\n";
  
    /*******************************************************
     * SATISFY REACTIVE TRANSITION FLUX RELATIVE TOLERANCE *
     *******************************************************/

    if(reactive_transition.size()) {
      //
      IO::log << "Sampling to satisfy reactive transition flux relative tolerance (importance samplings counted) ..." << std::endl;

      // checking
      //
      if(reactive_transition.size() != 2 || reactive_transition[0] == reactive_transition[1] ||
	 //
	 reactive_transition[0] < 0      || reactive_transition[0] >= _ms.species_size() ||
	 //
	  reactive_transition[1] < 0     || reactive_transition[1] >= _ms.species_size()) {
	//
	std::cerr << funame << "reactive transition pair is not properly defined\n";

        IO::log.flush();

	throw Error::Init();
      }

      // primitive's cycle
      //
      for(iterator mit = begin(); mit != end(); ++mit) {// surface cycle
	//
	proj_samp_num.clear();

	// facet cycle
	//
	for(sit = mit->begin(); sit != mit->end(); ++sit)
	  //
	  if(sit->second.run_traj_num()) {
	    //
	    int react, prod = -1;

	    for(int ward = 0; ward < 2; ++ward) {
	      //
	      switch(ward) {
		//
	      case BACKWARD:
		//
		react = sit->first.second;

		break;
		//
	      case FORWARD:
		//
		react = sit->first.first;
	      }
	    
	      if(reactant() >= 0 && reactant() != react)
		//
		continue;
	    
	      for(int i = 0; i < 2; ++i)
		//
		if(reactive_transition[i] == react)
		  //
		  prod = reactive_transition[1 - i];
	    
	      if(prod < 0)
		//
		continue;

	      if(itemp = sit->second.reac_traj_num(ward, prod)) {
		//
		// recrossing factor
		//
		dtemp = (double)itemp / (double)sit->second.run_traj_num();

		// facet reqctive flux standard deviation
		//
		double face_rmsd = std::sqrt(dtemp * (1. - dtemp)) * mit->face_flux(sit);
		
		dtemp = tran_rel_tol  * cross_flux_mean[DivSur::face_t(react, prod)];

		dtemp = face_rmsd * cross_flux_rmsd[DivSur::face_t(react, prod)] / dtemp / dtemp 
		//
		 * (double)sit->second.size() / (double)sit->second.run_traj_num();

		if(dtemp > (double)max_imp_size) {
		  //
		  IO::log << mit - begin() << "-th surface, " << sit->first << " facet: WARNING: requested importance samplings #, " 
			//
			  << std::floor(dtemp) << ", exceeds maximal facet importance samplings #, " << max_imp_size << ": truncating\n";

		    itemp = max_imp_size;
		}
		else
		  //
		  itemp = (long)dtemp;

		nit = proj_samp_num.find(sit->first);

		if(nit == proj_samp_num.end() || itemp > nit->second)
		  //
		  proj_samp_num[sit->first] = itemp;
	      }//
	      //
	    }// ward cycle
	    //
	  }// facet cycle
    
	if(proj_samp_num.size()) {
	  //
	  IO::log << "   " << mit - begin() << "-th surface:\n"
		  << std::setw(15) << "facet" 
		  << std::setw(15) << "projected" 
		  << std::setw(15) << "current" << "\n";

	  for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit)
	    //
	    IO::log << std::setw(15) << nit->first 
		    << std::setw(15) << nit->second 
		    << std::setw(15) << mit->find(nit->first)->second.size() 
		    << "\n";
	}

	count_max = 0;

	for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	  //
	  sit = mit->find(nit->first);

	  itemp = nit->second - sit->second.size();

	  if(itemp > 0) {
	    //
	    count_max += itemp;
	  }
	}

	if(!count_max)
	  //
	  continue;

	old_share = 0;

	// sampling loop
	//
	while(1) {
	  //
	  std::set<DivSur::face_t > face_work;

	  count = 0;

	  for(nit = proj_samp_num.begin(); nit != proj_samp_num.end(); ++nit) {
	    //
	    sit = mit->find(nit->first);

	    itemp = nit->second - sit->second.size();

	    if(itemp > 0) {
	      //
	      count += itemp;

	      face_work.insert(nit->first);
	    }
	  }

	  new_share = int((double)(count_max - count) / (double)count_max * 100.);
	  print_progress(old_share, new_share);

	  if(!face_work.size()) {
	    //
	    break;
	  }
	  else
	    //
	    run_again = true;

	  if(_sample(mit, face_work))
	    //
	    return false;
	  //
	}// sampling loop
	//
      }// surface cycle

      IO::log << "done\n\n";
    }

    /*************************************************************
     *********************** RUN TRAJECTORIES ********************
     *************************************************************/

    // run trajectories again
    //
    if(run_again) {
      //
      _stat_pot_hit += _pot.reset_count();

      IO::log << "running trajectories again ...\n";

      run_traj();

      IO::log << "done\n\n";

      _traj_pot_hit += _pot.reset_count();
    }
  }

  _print_sampling_results();

  if(job() == DYN_JOB) {
    //
    _print_traj_results();
    //_analize();
  }

  return true;
} 

CrossRate::MultiArray::MultiArray(const DivSur::MultiSur& ms, Potential::Wrap pot, Dynamic::CCP stop)
  //
  : std::vector<SurArray>(ms.primitive_size()), _ms(ms), _pot(pot), _stop(stop), _stat_pot_hit(0), _traj_pot_hit(0)
{
  const char funame [] = "CrossRate::MultiArray::MultiArray: ";
 
  // initial sampling to check which facets are actually present
  //
  IO::log << "preliminary sampling with " << min_sur_size << " samplings per primitive surface:\n";

  // surface cycle
  //
  for(iterator mit = begin(); mit != end(); ++mit) {
    //
    IO::log << "   " << mit - begin() << "-th surface:\n";

    // sampling loop
    //
    const int n = min_sur_size / smp_set_size;

    for(int s = 0; s < n; ++s) {
      //
      std::set<DivSur::face_t > face_work;

      for(SurArray::iterator sit = mit->begin(); sit != mit->end(); ++sit)
	//
	face_work.insert(sit->first);

      _sample(mit, face_work);
      //
    }// sampling loop

    // check for failures
    //
    if(mit->fail_num() > max_fail_num) {
      //
      std::cerr << funame << "too many failures for " << mit - begin() << "-th surface\n";

      IO::log.flush();

      throw Error::Run();
    }

    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      if(!sit->second.flux_num()) {
	//
	std::cerr << funame << mit - begin() << "-th surface, " << sit->first << " facet: " <<"no flux samplings\n";

	IO::log.flush();

	throw Error::Run();
      }//
      //
    }//
    //
  }// surface cycle

  IO::log << "done\n\n";

  IO::log << "preliminary sampling results:\n";

  for(const_iterator mit = begin(); mit != end(); ++mit) {
    //
    IO::log << "   " << mit - begin() << "-th surface:\n"
      //
	    << std::setw(15) << "facet"
	    << std::setw(15) << "inner"
	    << std::setw(15) << "exclude"
	    << std::setw(15) << "close atoms"
	    << std::setw(15) << "logic fails"
	    << "\n"
	    << std::setw(15) << "---"
	    << std::setw(15) << mit->inner_num()
	    << std::setw(15) << mit->excl_num()
	    << std::setw(15) << mit->close_num()
	    << std::setw(15) << mit->fail_num()
	    << "\n\n";
    
    IO::log << std::setw(15) << "facet"
	    << std::setw(15) << "flux number"
	    << std::setw(15) << "fail number"
	    << std::setw(15) << "fake number"
	    << std::setw(15) << "imp. number"
	    << std::setw(15) << "flux value"
	    << std::setw(15) << "emin, kcal"
	    << "\n";
      
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit)
      //
      IO::log << std::setw(15) << sit->first 
	      << std::setw(15) << sit->second.flux_num()
	      << std::setw(15) << sit->second.fail_num()
	      << std::setw(15) << sit->second.fake_num()
	      << std::setw(15) << sit->second.size()
	      << std::setw(15) << sit->second.flux_value()
	      << std::setw(15) << sit->second.min_ener() / Phys_const::kcal
	      << "\n";

    IO::log << "\n";
  }
  
  for(const_iterator mit = begin(); mit != end(); ++mit) {
    //
    for(SurArray::const_iterator sit = mit->begin(); sit != mit->end(); ++sit) {
      //
      if(!sit->second.flux_num()) {
	//
	std::cerr << mit - begin() << "-th surface,  " << sit->first << " facet: no potential samplings\n";

	IO::log.flush();

	throw Error::Run();
      }

      if(sit->second.min_ener() < reactive_energy()) {
	//
	std::cerr << mit - begin() << "-th surface,  " << sit->first << " facet: " 
		//
		  << "minimal energy = " << std::floor(sit->second.min_ener() /Phys_const::kcal * 10. + 0.5) / 10.
		//
		  << " kcal/mol is less than the reactive energy = " << std::floor(reactive_energy() / Phys_const::kcal * 10. + 0.5) / 10.
		//
		  << " kcal/mol\n";

	IO::log.flush();

	throw Error::Run();
      }
    }
  }

  // sample the surface and run trajectories to satisfy predefined tolerances
  //
  while(!_work()) {}

  IO::log << "potential sampling   calculations # = " << _stat_pot_hit << "\n"
	//
	  << "potential trajectory calculations # = " << _traj_pot_hit << "\n\n";
}
