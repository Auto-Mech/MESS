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

#include <ios>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <list>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <cstdarg>

#include "atom.hh"
#include "model.hh"
#include "units.hh"
#include "key.hh"
#include "io.hh"
#include "slatec.h"

// sparse matrix eigenvalue solver
//#include "feast.h"
//#include "feast_sparse.h"

namespace Model {

  SharedPointer<TimeEvolution> time_evolution;

  /********************************************************************************************
   **************************************** CONTROLS ******************************************
   ********************************************************************************************/

  // minimal interatomic distance
  double atom_dist_min = 1.4;
  // maximal tunneling action
  double Tunnel::_action_max = 100.;

  /********************************************************************************************
   **************************************** GLOBAL OBJECTS ************************************
   ********************************************************************************************/

  bool _isinit = false;
  bool isinit () { return _isinit; }

  bool _no_run = false;
  bool no_run () { return _no_run; }

  // global collision frequency model
  std::vector<SharedPointer<Collision> > _default_collision;
  //ConstSharedPointer<Collision> collision (int i) { return _collision[i]; }

  // energy transfer kernel
  std::vector<SharedPointer<Kernel> > _default_kernel;
  ConstSharedPointer<Kernel> default_kernel(int i) { return _default_kernel[i]; }

  // energy transfer kernel flags
  int Kernel::_flags = 0;

  // buffer fraction
  std::vector<double> _buffer_fraction;
  double buffer_fraction (int i) { return _buffer_fraction[i]; }
  int buffer_size () { return _buffer_fraction.size(); }

  std::vector<int> _escape_well_index;
  int escape_size       ()      { return _escape_well_index.size(); }
  int escape_well_index (int i) { return _escape_well_index[i];     }

  // bound species
  std::vector<Well> _well;
  // bimolecular products
  std::vector<SharedPointer<Bimolecular> > _bimolecular;
  // well-to-well barriers
  std::vector<SharedPointer<Species> >   _inner_barrier;
  // well-to-well connection scheme
  std::vector<std::pair<int ,int> >      _inner_connect;
  // well-to-bimolecular barriers
  std::vector<SharedPointer<Species> >   _outer_barrier;
  // well-to-bimolecular connection scheme
  std::vector<std::pair<int ,int> >      _outer_connect;

  int            well_size () { return            _well.size(); }
  int     bimolecular_size () { return     _bimolecular.size(); }
  int   inner_barrier_size () { return   _inner_barrier.size(); }
  int   outer_barrier_size () { return   _outer_barrier.size(); }

  const Well&                                   well (int w) { return            _well[w]; }
  const Bimolecular&                     bimolecular (int p) { return    *_bimolecular[p]; }
  const Species&                       inner_barrier (int b) { return  *_inner_barrier[b]; }
  const Species&                       outer_barrier (int b) { return  *_outer_barrier[b]; }
  const std::pair<int, int>&           inner_connect (int b) { return   _inner_connect[b]; }
  const std::pair<int, int>&           outer_connect (int b) { return   _outer_connect[b]; }

  double  _maximum_barrier_height;
  double  maximum_barrier_height() { return _maximum_barrier_height; }

  // bimolecular product to be used as a reference
  std::string reactant;
  double _energy_shift = 0.;
  double energy_shift () { return _energy_shift; }

  /*************************************** ENERGY LIMIT **************************************/

  bool   _is_energy_limit = false;
  double _energy_limit;
  void set_energy_limit (double e) { _is_energy_limit = true; _energy_limit = e; }
  bool is_energy_limit () { return _is_energy_limit; }

  double  energy_limit ()  
  {
    const char funame [] = "Model::energy_limit: ";

    if(!_is_energy_limit) {
      std::cerr << funame << "not defined\n";
      IO::log << std::flush;
      throw Error::Init();
    }
    return _energy_limit;
  }

  /***************************************** HELPERS ******************************************/

  std::vector<D3::Vector> adjusted_normal_mode (const std::vector<Atom>&, const std::map<int, D3::Vector>&,
						Lapack::Vector* =0);
  double vibrational_sum (double ener, const Lapack::Vector& freq, double, int = 0);
  std::vector<std::vector<int> > population (double ener, const std::vector<double>& freq, int = 0); 
  int rotation_matrix_element (int m, int n, int k, double& fac); // <m|f(k)|n>, f(k)=sin, cos
}

Model::TimeEvolution::TimeEvolution (IO::KeyBufferStream& from) 
  : _start(-1.), _finish(-1.), _size(-1), _reactant(-1), _excess(-1.), _temperature(-1.)
{
  const char funame [] = "Model::TimeEvolution::TimeEvolution: ";

  int         itemp;
  double      dtemp;
  std::string stemp;

  KeyGroup TimeEvolutionModel;

  Key start_key("Start[s]" );
  Key   fin_key("Finish[s]" );
  Key  size_key("Size"      );
  Key  reac_key("Reactant"  );
  Key   exc_key("ExcessReactantConcentration[molecule/cm^3]");
  Key   ext_key("EffectiveTemperature[K]");
  Key   out_key("TimeOutput");
  
  std::string token, comment;

  while(from >> token) {
    // end input 
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // excess reactant concentration
    else if(exc_key == token) {
      if(_excess > 0. || _temperature > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> _excess)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_excess <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      
      _excess /= Phys_const::cm * Phys_const::cm * Phys_const::cm;
    }
    // effective  temperature
    else if(ext_key == token) {
      if(_excess > 0. || _temperature > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> _temperature)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_temperature <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      
      _temperature *= Phys_const::kelv;
    }
    // start time
    else if(start_key == token) {
      if(_start >= 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> _start)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_start <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      
      _start /= Phys_const::herz;
    }
    // finish time
    else if(fin_key == token) {
      if(_finish >= 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> _finish)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_finish <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      
      _finish /= Phys_const::herz;
    }
    // time grid size
    else if(size_key == token) {
      if(_size > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> _size)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_size < 2) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // reactant
    else if(reac_key == token) {
      if(_reactant_name.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> _reactant_name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
    }
    // output stream
    else if(out_key == token) {
      if(out.is_open()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      
      out.open(stemp.c_str());

      if(!out) {
	std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	throw Error::Open();
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

  if(!_reactant_name.size()) {
    std::cerr << funame << "reactant name not defined\n";
    throw Error::Init();
  }

  if(_excess < 0. && _temperature < 0.) {
    std::cerr << funame << "excess reactant concentration not initialized\n";
    throw Error::Init();
  }
  if(_start <= 0.) {
    std::cerr << funame << "time grid start not initialized\n";
    throw Error::Init();
  }
  if(_finish <= 0.) {
    std::cerr << funame << "time grid finish not initialized\n";
    throw Error::Init();
  }
  if(_size < 2) {
    std::cerr << funame << "time grid size not initialized\n";
    throw Error::Init();
  }
  if(_finish <= _start) {
    std::cerr << funame << "time grid finish should be bigger than start\n";
    throw Error::Init();
  }
  if(!out.is_open()) {
    std::cerr << funame << "output stream not initialized\n";
    throw Error::Init();
  }

  _step = std::pow(_finish / _start, 1. / double(_size - 1));

}

void Model::TimeEvolution::set_reactant () const
{
  const char funame [] = "Model::TimeEvolution::set_reactant: ";

  if(_excess > 0.) {
    for(_reactant = 0; _reactant < bimolecular_size(); ++_reactant)
      if(bimolecular(_reactant).name() == _reactant_name)
	return;
    std::cerr << funame << "unknown bimolecular: " << _reactant_name <<"\n";
    throw Error::Range();
  }
  
  
  for(_reactant = 0; _reactant < well_size(); ++_reactant)
    if(well(_reactant).name() == _reactant_name)
      return;

  std::cerr << funame << "unknown well: " << _reactant_name <<"\n";
  throw Error::Range();
}

/********************************************************************************************
 ******************************** MODEL INITIALIZATION **************************************
 ********************************************************************************************/

void Model::init (IO::KeyBufferStream& from) 
{
  const char funame [] = "Model::init: ";

  if(isinit()) {
    std::cerr << funame << "allready initialized\n";
    throw Error::Init();
  }

  if(!IO::out.is_open()) {
    std::cerr << funame  << "output stream is not open\n";
    throw Error::Init();
  }      

  if(!IO::log.is_open()) {
    std::cerr << funame  << "log stream is not open\n";
    throw Error::Init();
  }      

  IO::Marker funame_marker(funame);

  _isinit = true;

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  KeyGroup ModelInit;

  Key     kflag_key("EnergyRelaxationFlag"    );
  Key      freq_key("CollisionFrequency"      );
  Key       cer_key("EnergyRelaxation"        );
  Key      buff_key("BufferFraction"          );
  Key      well_key("Well"                    );
  Key      barr_key("Barrier"                 );
  Key     bimol_key("Bimolecular"             );
  Key     tstep_key("OutputTemperatureStep[K]");
  Key      tmin_key("OutputTemperatureMin[K]" );
  Key     tsize_key("OutputTemperatureSize"   );
  Key     eref_key("OutputReferenceEnergy[kcal/mol]");
  Key     wout_key("ThermodynamicDataOutput");
  Key    tincr_key("RelativeTemperatureIncrement");

  /************************************* INPUT *********************************************/

  double eref = 0.;
  int tstep = 100;
  int tsize = 20;
  int tmin = -1;

  std::string token, name, comment;
  std::pair<std::string, std::string> species_pair;
  std::vector<std::pair<std::string, std::string> > connect_verbal;
  typedef std::vector<std::pair<std::string, std::string> >::const_iterator It;
  std::vector<SharedPointer<Species> > barrier;   // barrier pool
  std::set<std::string> species_name;  // all wells, barriers, and bimolecular products names
  std::map<std::string, int>        well_index;
  std::map<std::string, int> bimolecular_index;

  std::string wout_file;
  double temp_rel_incr = 0.001;


  while(from >> token) {
    // end input 
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // relative temperature increment
    else if(tincr_key == token) {
      if(!(from >> temp_rel_incr)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      if(temp_rel_incr <= 0. || temp_rel_incr >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // weight output
    else if(wout_key == token) {
      if(!(from >> wout_file)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // output reference energy
    else if(eref_key == token) {
      if(!(from >> eref)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      eref *= Phys_const::kcal;
    }
    // output temperature step
    else if(tstep_key == token) {
      if(!(from >> tstep)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(tstep <= 0) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
	
    }
    // output temperature start
    else if(tmin_key == token) {
      if(!(from >> tmin)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(tmin <= 0) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
	
    }
    // output temperature size
    else if(tsize_key == token) {
      if(!(from >> tsize)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(tsize <= 0) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }	
    }
    // collision frequency model
    else if(freq_key == token) {
      _default_collision.push_back(new_collision(from));
    }
    // energy relaxation kernel
    else if(cer_key == token) {
      _default_kernel.push_back(new_kernel(from));
    }
    // buffer gas fraction
    else if(buff_key == token) {
      if(_buffer_fraction.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }

      IO::LineInput lin(from);
      while(lin >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}
	_buffer_fraction.push_back(dtemp);
      }

      if(!_buffer_fraction.size()) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Init();
      }

      dtemp = 0.;
      for(int i = 0; i < _buffer_fraction.size(); ++i)
	dtemp += _buffer_fraction[i];

      for(int i = 0; i < _buffer_fraction.size(); ++i)
	_buffer_fraction[i] /= dtemp;
    }
    // energy relaxation kernel flags
    else if(kflag_key == token) {
      from >> stemp;
      std::getline(from, comment);
      if(stemp == "up")
	Kernel::add_flag(Kernel::UP);
      else if(stemp == "density")
	Kernel::add_flag(Kernel::DENSITY);
      else if(stemp == "notruncation")
	Kernel::add_flag(Kernel::NOTRUN);
      else {
	std::cerr << funame << token << ": unknown key: " << stemp 
		  << " possible keys: up, density, notruncation\n";
	throw Error::Range();
      }
    }
    // new well
    else if(well_key == token) {
      if(!(from >> name)) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(species_name.find(name) != species_name.end()) {
	std::cerr << funame << token << ": name " << name << " already in use\n";
	throw Error::Logic();
      }
      species_name.insert(name);
      well_index[name] = _well.size();
      IO::log << IO::log_offset << "WELL: " << name << "\n";
      _well.push_back(Well(from, name));
      //_well.rbegin()->set_name(name);
    }
    // new bimolecular
    else if(bimol_key == token) {
      if(!(from >> name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(species_name.find(name) != species_name.end()) {
	std::cerr << funame << token << ": name " << name << " already in use\n";
	throw Error::Logic();
      }
      species_name.insert(name);
      bimolecular_index[name] = _bimolecular.size();
      IO::log << IO::log_offset << "BIMOLECULAR: " << name << "\n";
      _bimolecular.push_back(new_bimolecular(from, name));
      //_bimolecular.back()->set_name(name);
    }
    // new barrier
    else if(barr_key == token) {
      if(!(from >> name >> species_pair.first >> species_pair.second)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(species_name.find(name) != species_name.end()) {
	std::cerr << funame << token << ": name " << name << " already in use\n";
	throw Error::Logic();
      }
      if(species_pair.first == species_pair.second) {
	std::cerr << funame << name << "barrier connects " 
		  << species_pair.first <<" well with itself\n";
	throw Error::Logic();
      }
      for(It b = connect_verbal.begin(); b != connect_verbal.end(); ++b) {
	itemp = b - connect_verbal.begin();
	if(b->first == species_pair.second && b->second == species_pair.first ||
	   *b == species_pair) {
	  std::cerr << funame << name << " and " << barrier[itemp]->name()
		    << " barriers connect the same pair of species\n";
	  throw Error::Logic();
	}
      } 
      species_name.insert(name);
      connect_verbal.push_back(species_pair);
      IO::log << IO::log_offset << "BARRIER: " << name << "\n";
      barrier.push_back(new_species(from, name, NUMBER));
      //barrier.back()->set_name(name);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
  
  /************************************* CHECKING *********************************************/

  {
    IO::Marker  input_marker("checking input", IO::Marker::ONE_LINE| IO::Marker::NOTIME);

    if(!from) {
      std::cerr << funame << "input stream is currupted\n";
      throw Error::Input();
    }

    if(!_default_collision.size()) {
      std::cerr << funame  << "default collision model is not initialized\n";
      throw Error::Init();
    }

    if(!_default_kernel.size()) {
      std::cerr << funame  << "default energy transfer kernel is not defined\n";
      throw Error::Init();
    }

    if(_default_collision.size() != _default_kernel.size()) {
      std::cerr << funame  << "number of energy transfer kernels and collision frequency models mismatch\n";
      throw Error::Init();
    }

    if(_default_collision.size() == 1) {
      _buffer_fraction.resize(1);
      _buffer_fraction[0] = 1.;
    }
    
    if(_default_collision.size() != _buffer_fraction.size()) {
      std::cerr << funame  << "number of collision frequency models and buffer gas fractions mismatch\n";
      throw Error::Init();
    }

    if(!well_size()) {
      std::cerr << funame << "no wells\n";
      throw Error::Init();
    }

    for(int w = 0; w < well_size(); ++w)
      if(well(w).escape())
	_escape_well_index.push_back(w);
  }

  /****************************** SETTING CONNECTIVITY SCHEME **********************************/

  {
    IO::Marker  set_marker("setting connectivity", IO::Marker::ONE_LINE| IO::Marker::NOTIME);

    for(It b = connect_verbal.begin(); b != connect_verbal.end(); ++b) {
      itemp = b - connect_verbal.begin();
      // inner barrier
      if(well_index.find(b->first)  != well_index.end() &&
	 well_index.find(b->second) != well_index.end()) {
	_inner_barrier.push_back(barrier[itemp]);
	_inner_connect.push_back(std::make_pair(well_index[b->first], well_index[b->second]));
      }
      // outer barrier
      else if(well_index.find(b->first) != well_index.end() &&
	      bimolecular_index.find(b->second)  != bimolecular_index.end()) {
	_outer_barrier.push_back(barrier[itemp]);
	_outer_connect.push_back(std::make_pair(well_index[b->first], 
						bimolecular_index[b->second]));
      }
      // outer barrier
      else if(well_index.find(b->second) != well_index.end() &&
	      bimolecular_index.find(b->first)  != bimolecular_index.end()) {
	_outer_barrier.push_back(barrier[itemp]);
	_outer_connect.push_back(std::make_pair(well_index[b->second], 
						bimolecular_index[b->first]));
      }
      // wrong connection
      else {
	std::cerr << funame << barrier[itemp]->name() << " barrier connects wrong species: "
		  << b->first << " and " << b->second << "\n";
	throw Error::Logic();
      }
    }

  }

  /****************************** CONNECTIVITY CHECK ***********************************/

  {
    IO::Marker  check_marker("checking connectivity", IO::Marker::ONE_LINE| IO::Marker::NOTIME);

    std::vector<std::set<int> > pool;
    typedef std::vector<std::set<int> >::iterator Pit;
    // well connectivity
    for(int b = 0; b < inner_barrier_size(); ++b) {// inner barrier cycle
      Pit p1, p2;
      for(p1 = pool.begin(); p1 != pool.end(); ++p1)
	if(p1->find(inner_connect(b).first) != p1->end())
	  break;
      for(p2 = pool.begin(); p2 != pool.end(); ++p2)
	if(p2->find(inner_connect(b).second) != p2->end())
	  break;


      if(p1 == pool.end() && p2 == pool.end()) {
	pool.push_back(std::set<int>());
	pool.rbegin()->insert(inner_connect(b).first);
	pool.rbegin()->insert(inner_connect(b).second);
      }
      else if(p1 == pool.end())
	p2->insert(inner_connect(b).first);
      else if(p2 == pool.end())
	p1->insert(inner_connect(b).second);
      else if(p1 != p2) {
	p1->insert(p2->begin(), p2->end());
	pool.erase(p2);
      }
    }// inner barrier cycle

    if(!pool.size()) {
      if(well_size() > 1) {
	std::cerr << funame << "wells are not connected\n";
	throw Error::Init();
      }
    }
    else if(pool.size() > 1) {
      std::cerr << funame << "there are " << pool.size() << " unconnected groups of wells: ";
      for(Pit p = pool.begin(); p != pool.end(); ++p) {
	if(p != pool.begin())
	  std::cerr << ", ";
	for(std::set<int>::iterator s = p->begin(); s != p->end(); ++s) {
	  if(s != p->begin())
	    std::cerr << "+";
	  std::cerr << well(*s).name();
	}
      }
      std::cerr << std::endl;
      throw Error::Input();
    }
    else if(pool.begin()->size() != well_size()) {
      std::cerr << funame << "there are unconnected wells:";
      for(int w = 0; w < well_size(); ++w)
	if(pool.begin()->find(w) == pool.begin()->end()) 
	  std::cerr << " " << well(w).name();
      std::cerr << std::endl;
      throw Error::Input();
    }
  
    // product connectivity
    std::set<int> product_pool;
    for(int b = 0; b < outer_barrier_size(); ++b)
      product_pool.insert(outer_connect(b).second);
    if(product_pool.size() != bimolecular_size()) {
      std::cerr << funame << "there are unconnected bimolecular products:";
      for(int p = 0; p < bimolecular_size(); ++p)
	if(product_pool.find(p) == product_pool.end()) 
	  std::cerr << " " << bimolecular(p).name();
      std::cerr << std::endl;
      throw Error::Input();
    }
  }
 
  /****************************** SHIFTING ENERGY ***********************************/

  // bimolecular reactant index
  if(bimolecular_index.find(reactant) != bimolecular_index.end()) {
    IO::Marker zero_marker("shifting energy zero", IO::Marker::ONE_LINE | IO::Marker::NOTIME);

    _energy_shift = -bimolecular(bimolecular_index[reactant]).ground();
    for(int w = 0; w < well_size(); ++w)
      _well[w].shift_ground(_energy_shift);
    for(int b = 0; b < inner_barrier_size(); ++b)
      _inner_barrier[b]->shift_ground(_energy_shift);
    for(int b = 0; b < outer_barrier_size(); ++b)
      _outer_barrier[b]->shift_ground(_energy_shift);
    for(int p = 0; p < bimolecular_size(); ++p)
      _bimolecular[p]->shift_ground(_energy_shift);
  }

  /****************************** DISSOCIATION LIMIT ***********************************/

  {
    IO::Marker diss_marker("setting dissociation limit", IO::Marker::ONE_LINE | IO::Marker::NOTIME);

    for(int w = 0; w < well_size(); ++w) {
      btemp = true;
      for(int b = 0; b < outer_barrier_size(); ++b) {
	dtemp =  outer_barrier(b).real_ground();
	if(outer_connect(b).first == w && (btemp || dtemp < _well[w].dissociation_limit)) {
	  btemp = false;
	  _well[w].dissociation_limit = dtemp;
	}
      }

      for(int b = 0; b < inner_barrier_size(); ++b) {
	dtemp =  inner_barrier(b).real_ground();
	if((inner_connect(b).first == w || inner_connect(b).second == w) && 
	   (btemp || dtemp < _well[w].dissociation_limit)) {
	  btemp = false;
	  _well[w].dissociation_limit = dtemp;
	}
      }
      if(btemp) {
	std::cerr << funame << "no barrier associated with " << well(w).name() << " well found\n";
	throw Error::Init();
      }
    }
  }

  /****************************** MAXIMUM BARIER HEIGHT ***********************************/

  {
    IO::Marker height_marker("setting maximum barrier height", IO::Marker::ONE_LINE | IO::Marker::NOTIME);

    btemp = true;
    for(int b = 0; b < outer_barrier_size(); ++b) {
      dtemp = outer_barrier(b).real_ground();
      if(btemp || _maximum_barrier_height < dtemp) {
	btemp = false;
	_maximum_barrier_height = dtemp;
      }
    }
    for(int b = 0; b < inner_barrier_size(); ++b) {
      dtemp = inner_barrier(b).real_ground();
      if(btemp || _maximum_barrier_height < dtemp) {
	btemp = false;
	_maximum_barrier_height = dtemp;
      }
    }
  }
      
  /************************************** OUTPUT ***************************************/

  IO::out.precision(1);
  IO::out << std::fixed;

  IO::out << "Wells (G - ground energy, D - dissociation limit, kcal/mol):\n"
	  << IO::first_offset 
	  << std::setw(5) << "Name"
	  << std::setw(9) << "G"
	  << std::setw(9) << "D"
	  << "\n";
  for(int w = 0; w < well_size(); ++w)
    IO::out << IO::first_offset 
	    << std::setw(5) << well(w).name()
	    << std::setw(9) << well(w).ground() / Phys_const::kcal
	    << std::setw(9) << well(w).dissociation_limit / Phys_const::kcal
	    << "\n";
  IO::out << "\n";

  if(bimolecular_size()) {
    IO::out << "Bimolecular Products (G - ground energy, kcal/mol):\n"
	    << IO::first_offset 
	    << std::setw(5) << "Name"
	    << std::setw(9) << "G"
	    << "\n";
    for(int p = 0; p < bimolecular_size(); ++p)
      if(!bimolecular(p).dummy())
	IO::out << IO::first_offset 
		<< std::setw(5) << bimolecular(p).name()
		<< std::setw(9) << bimolecular(p).ground() / Phys_const::kcal
		<< "\n";
    IO::out << "\n";
  }

  if(outer_barrier_size()) {
    IO::out << "Well-to-Bimolecular Barriers (H/G - barrier height/well depth, kcal/mol)\n"
	    << IO::first_offset 
	    << std::setw(5) << "Name"
	    << std::setw(9) << "H"
	    << std::setw(9) << "Well"
	    << std::setw(9) << "G"
	    << std::setw(9) << "Product"
	    << "\n";
    for(int b = 0; b < outer_barrier_size(); ++b) {
      int w = outer_connect(b).first;
      int p = outer_connect(b).second;
      IO::out << IO::first_offset 
	      << std::setw(5) << outer_barrier(b).name()
	      << std::setw(9) << outer_barrier(b).real_ground() / Phys_const::kcal
	      << std::setw(9) << well(w).name()
	      << std::setw(9) << well(w).ground() / Phys_const::kcal
	      << std::setw(9) << bimolecular(p).name()
	      << "\n";
    }
    IO::out << "\n";
  }
  
  if(inner_barrier_size()) {
    IO::out << "Well-to-Well Barriers (H/G - barrier height/well depth, kcal/mol):\n"
	    << IO::first_offset 
	    << std::setw(5) << "Name"
	    << std::setw(9) << "H"
	    << std::setw(9) << "Well"
	    << std::setw(9) << "G1"
	    << std::setw(9) << "Well"
	    << std::setw(9) << "G2"
	    << "\n";
    for(int b = 0; b < inner_barrier_size(); ++b) {
      int w1 = inner_connect(b).first;
      int w2 = inner_connect(b).second;
      IO::out << IO::first_offset 
	      << std::setw(5) << inner_barrier(b).name()
	      << std::setw(9) << inner_barrier(b).real_ground() / Phys_const::kcal
	      << std::setw(9) << well(w1).name()
	      << std::setw(9) << well(w1).ground() / Phys_const::kcal
	      << std::setw(9) << well(w2).name()
	      << std::setw(9) << well(w2).ground() / Phys_const::kcal
	      << "\n";
    }
    IO::out << "\n";
  }

  IO::out.precision(6);
  IO::out << std::resetiosflags(std::ios_base::floatfield);

  IO::out << "________________________________________________________________________________________\n\n";

  // partition functions output

  // bimolecular partition function units
  const double bpu = Phys_const::cm * Phys_const::cm * Phys_const::cm;


  if(wout_file.size()) {

    std::ofstream wout(wout_file.c_str());
 
    wout << "partition function logs (1/cm^3, relative to the ground level) and their derivatives:\n"
	 << std::setw(7) << "T\\Q";

    // wells
    for(int w = 0; w < well_size(); ++w)
      wout << std::setw(13) << well(w).name()
	   << std::setw(13) << "first"
	   << std::setw(13) << "second";

    // bimolecular products
    for(int p = 0; p < bimolecular_size(); ++p)
      if(!bimolecular(p).dummy())
	for(int f = 0; f < 2; ++f)
	  wout << std::setw(13) << bimolecular(p).fragment_name(f)
	       << std::setw(13) << "first"
	       << std::setw(13) << "second";

    wout << "\n";

    // temperature cycle
    for(int t = 0; t <= tsize; ++t) {

      if(t < tsize) {
	if(tmin > 0)
	  itemp = tmin;
	else
	  itemp = tstep;
	itemp += t * tstep;
	dtemp = (double)itemp;
      }
      else
	dtemp = 298.15;

      wout << std::setw(7) << dtemp;

      const double tval = dtemp * Phys_const::kelv;

      double temp_incr = tval * temp_rel_incr;

      double tt[3];
      tt[0] = tval;
      tt[1] = tval - temp_incr;
      tt[2] = tval + temp_incr;

      temp_incr /= Phys_const::kelv;

      double zz[3];

      // well partition functions and derivatives
      for(int w = 0; w < well_size(); ++w) {
	for(int i = 0; i < 3; ++i)
	  zz[i] = std::log(well(w).weight(tt[i]) * std::pow(well(w).mass() * tt[i] / 2. / M_PI, 1.5) * bpu);

	wout << std::setw(13) << zz[0]
	     << std::setw(13) << (zz[2] - zz[1]) / 2. / temp_incr
	     << std::setw(13) << (zz[2] + zz[1] - 2. * zz[0]) / temp_incr / temp_incr;

      }

      // bimolecular product partition functions
      for(int p = 0; p < bimolecular_size(); ++p) 
	if(!bimolecular(p).dummy())
	  for(int f = 0; f < 2; ++f) {
	    for(int i = 0; i < 3; ++i)
	      zz[i] = std::log(bimolecular(p).fragment_weight(f, tt[i]) * bpu);
	  
	    wout << std::setw(13) << zz[0]
		 << std::setw(13) << (zz[2] - zz[1]) / 2. / temp_incr
		 << std::setw(13) << (zz[2] + zz[1] - 2. * zz[0]) / temp_incr / temp_incr;
	  }
      wout << "\n";
    }// temperature cycle
  }

  IO::log << IO::log_offset << "partition functions (relative to the ground level, units - 1/cm3):\n"
	  << IO::log_offset << std::setw(5) << "T\\Q";

  // wells
  for(int w = 0; w < well_size(); ++w)
    IO::log << std::setw(13) << well(w).name();
  // bimolecular products
  for(int p = 0; p < bimolecular_size(); ++p)
    if(!bimolecular(p).dummy())
      for(int f = 0; f < 2; ++f)
	IO::log << std::setw(13) << bimolecular(p).fragment_name(f);
  // inner barriers
  for(int b = 0; b < inner_barrier_size(); ++b)
    IO::log << std::setw(13) << inner_barrier(b).name();
  // outer barriers
  for(int b = 0; b < outer_barrier_size(); ++b)
    IO::log << std::setw(13) << outer_barrier(b).name();
  IO::log << "\n";

  // temperature cycle
  for(int t = 0; t < tsize; ++t) {
    if(tmin > 0)
      itemp = tmin;
    else
      itemp = tstep;
    itemp += t * tstep;
    const double tval = itemp * Phys_const::kelv;

    IO::log << IO::log_offset << std::setw(5) << itemp;

    // well partition functions and derivatives
    for(int w = 0; w < well_size(); ++w)
      IO::log << std::setw(13) << well(w).weight(tval) * std::pow(well(w).mass() * tval / 2. / M_PI, 1.5) * bpu;
    // bimolecular product partition functions
    for(int p = 0; p < bimolecular_size(); ++p) 
      if(!bimolecular(p).dummy())
	for(int f = 0; f < 2; ++f)
	  IO::log << std::setw(13) << bimolecular(p).fragment_weight(f, tval) * bpu;
    // inner barrier partition functions
    for(int b = 0; b < inner_barrier_size(); ++b)
      IO::log << std::setw(13) << inner_barrier(b).weight(tval) 
	* std::exp((inner_barrier(b).real_ground() - inner_barrier(b).ground())/ tval)
	* std::pow(inner_barrier(b).mass() * tval / 2. / M_PI, 1.5) * bpu;

    // outer barrier partition functions
    for(int b = 0; b < outer_barrier_size(); ++b)
      IO::log << std::setw(13) << outer_barrier(b).weight(tval) 
	* std::exp((outer_barrier(b).real_ground() - outer_barrier(b).ground())/ tval)
	* std::pow(outer_barrier(b).mass() * tval / 2. / M_PI, 1.5) * bpu;
    IO::log << "\n";
  }// temperature cycle

  IO::log << IO::log_offset << "partition functions (relative to " << eref / Phys_const::kcal
	  << " kcal/mol, units - 1/cm3):\n"
	  << IO::log_offset << std::setw(5) << "T\\Q";

  // wells
  for(int w = 0; w < well_size(); ++w)
    IO::log << std::setw(13) << well(w).name();
  // inner barriers
  for(int b = 0; b < inner_barrier_size(); ++b)
    IO::log << std::setw(13) << inner_barrier(b).name();
  // outer barriers
  for(int b = 0; b < outer_barrier_size(); ++b)
    IO::log << std::setw(13) << outer_barrier(b).name();
  IO::log << "\n";

  // temperature cycle
  for(int t = 0; t < tsize; ++t) {
    if(tmin > 0)
      itemp = tmin;
    else
      itemp = tstep;
    itemp += t * tstep;
    const double tval = itemp * Phys_const::kelv;

    IO::log << IO::log_offset << std::setw(5) << itemp;

    // well partition functions
    for(int w = 0; w < well_size(); ++w)
      IO::log << std::setw(13) << well(w).weight(tval) 
	* std::exp((eref - well(w).ground()) / tval)
	* std::pow(well(w).mass() * tval / 2. / M_PI, 1.5) * bpu;

    // inner barrier partition functions
    for(int b = 0; b < inner_barrier_size(); ++b)
      IO::log << std::setw(13) << inner_barrier(b).weight(tval) 
	* std::exp((eref - inner_barrier(b).ground())/ tval)
	* std::pow(inner_barrier(b).mass() * tval / 2. / M_PI, 1.5) * bpu;

    // outer barrier partition functions
    for(int b = 0; b < outer_barrier_size(); ++b)
      IO::log << std::setw(13) << outer_barrier(b).weight(tval) 
	* std::exp((eref - outer_barrier(b).ground())/ tval)
	* std::pow(outer_barrier(b).mass() * tval / 2. / M_PI, 1.5) * bpu;
    IO::log << "\n";
  }// temperature cycle

  IO::log << IO::log_offset 
	  << "tunneling partition function correction factors &"
    " high pressure effective energy shifts (D), kcal/mol:\n"
	  << IO::log_offset << std::setw(5) << "T\\B";
  // inner barrier
  for(int b = 0; b < inner_barrier_size(); ++b)
    IO::log << std::setw(13) << inner_barrier(b).name()
	    << std::setw(13) << "D";
  // outer barrier
  for(int b = 0; b < outer_barrier_size(); ++b)
    IO::log << std::setw(13) << outer_barrier(b).name()
	    << std::setw(13) << "D";
  IO::log << "\n";

  // temperature cycle
  for(int t = 0; t < tsize; ++t) {
    if(tmin > 0)
      itemp = tmin;
    else
      itemp = tstep;
    itemp += t * tstep;
    const double tval = itemp * Phys_const::kelv;

    IO::log << IO::log_offset << std::setw(5) << itemp;

    // inner barrier
    for(int b = 0; b < inner_barrier_size(); ++b) {
      dtemp = std::log(inner_barrier(b).tunnel_weight(tval)) * tval 
	+ inner_barrier(b).real_ground() - inner_barrier(b).ground();

      IO::log << std::setw(13) << std::exp(dtemp / tval)
	      << std::setw(13) << -dtemp / Phys_const::kcal;
    }
    // outer barrier
    for(int b = 0; b < outer_barrier_size(); ++b) {
      dtemp = std::log(outer_barrier(b).tunnel_weight(tval)) * tval 
	+ outer_barrier(b).real_ground() - outer_barrier(b).ground();
      
      IO::log << std::setw(13) << std::exp(dtemp / tval)
	      << std::setw(13) << -dtemp / Phys_const::kcal;
    }
    IO::log << "\n";
  }// temperature cycle

}// model initialization

/********************************************************************************************
 ************************************* COMMON FUNCTONS **************************************
 ********************************************************************************************/

// <m|f(k)|n> for m, n, and k harmonics
//
int Model::rotation_matrix_element (int m, int n, int p, double& fac)
{
  if(!p && n == m)
    //
    return 0;
  
  if(!m && n == p || !n && m == p) {
    //
    fac /= M_SQRT2;

    return 0;
  }

  // all three fourier harmonics are cosines
  //
  if(m % 2 + n % 2 + p % 2 == 0 && (m + n == p || m + p == n || n + p == m)) {
    //
    fac /= 2.;

    return 0;
  }

  // two fourier harmonics are sines and  one is cosine
  //
  if(m % 2 + n % 2 + p % 2 == 2) {
    //
    if(!(m % 2)) {
      //
      std::swap(m, p);
    }

    if(!(n % 2)) {
      //
      std::swap(n, p);
    }

    if(m == n + p || n == m + p) {
      //
      fac /= 2.;

      return 0;
    }

    if(p == m + n + 2) {
      //
      fac /= -2.;

      return 0;
    }
  }
  return 1;
}

// adjust 1D motion so that the angular and linear momenta dissapear
// center of mass is assumed to be at zero
std::vector<D3::Vector> Model::adjusted_normal_mode 
(const std::vector<Atom>& atom, const std::map<int, D3::Vector>& motion,  Lapack::Vector* avp)
{
  const char funame [] = "Model::adjusted_normal_mode: ";

  if(motion.begin()->first < 0 || motion.rbegin()->first >= atom.size()) {
    std::cerr << funame << "motion index out of range\n";
    throw Error::Range();
  }

  double dtemp;

  // linear and angular momenta
  D3::Vector am; // angular momentum
  D3::Vector lm;// linear  momentum
  std::map<int, D3::Vector>::const_iterator mit;
  for(mit = motion.begin(); mit != motion.end(); ++mit) {
    am += atom[mit->first].mass() * D3::vprod(atom[mit->first], mit->second);
    lm += atom[mit->first].mass() * mit->second;
  }

  // total mass
  dtemp = 0.;
  for(int a = 0; a < atom.size(); ++a)
    dtemp += atom[a].mass();

  // linear velocity
  lm /= dtemp;

  // inverted inertia moment matrix; center of mass is assumed to be at zero
  Lapack::SymmetricMatrix im = inertia_moment_matrix(atom).invert();
      
  // angular velocity      
  Lapack::Vector av = im * am;
  if(avp)
    *avp = av;

  // normal mode adjusted
  std::vector<D3::Vector> res(atom.size());
  for(int a = 0; a < atom.size(); ++a)
    res[a] = D3::vprod(atom[a], av) - lm;
  for(mit = motion.begin(); mit != motion.end(); ++mit)
    res[mit->first] += mit->second;

  return res;
}

Lapack::SymmetricMatrix Model::inertia_moment_matrix(const std::vector<Atom>& atom)
{
  double dtemp;
  Lapack::SymmetricMatrix res(3);
  res = 0.;
  for(int i = 0; i < 3; ++i)
    for(int j = i; j < 3; ++j)
      for(std::vector<Atom>::const_iterator at = atom.begin(); at != atom.end(); ++at)
	res(i, j) -= at->mass() * (*at)[i] * (*at)[j];
      
  dtemp = 0.;
  for(int i = 0; i < 3; ++i)
    dtemp += res(i, i);

  for(int i = 0; i < 3; ++i)
    res(i, i) -= dtemp;

  return res;
}

void Model::shift_cm_to_zero(std::vector<Atom>& atom)
{
  D3::Vector vtemp;

  // shift the zero to the center of mass
  double mass = 0.;
  for(std::vector<Atom>::const_iterator at = atom.begin(); at != atom.end(); ++at) {
    vtemp += at->mass() * (*at);
    mass  += at->mass();
  }
  vtemp /= mass;
  for(std::vector<Atom>::iterator at = atom.begin(); at != atom.end(); ++at)
    *at -= vtemp;
}

// vibrational population inside given energy
std::vector<std::vector<int> > Model::population (double ener, const std::vector<double>& freq, int findex) 
{
  const char funame [] = "Model::population: ";

  if(findex < 0 || findex > freq.size()) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

  if(ener <= 0.) {
    std::cerr << funame << "negative energy\n";
    throw Error::Range();
  }

  if(findex == freq.size())
      return std::vector<std::vector<int> >(1); 

  double f = freq[findex];
  if(f <= 0.) {
    std::cerr << funame << findex << "-th frequency negative\n";
    throw Error::Range();
  }

  ++findex;
  std::vector<std::vector<int> > res;
  std::vector<std::vector<int> > vv;
  std::vector<std::vector<int> >::iterator vit;
  int n = 0;
  while(ener > 0.) {
    vv = population(ener, freq, findex);
    for(vit = vv.begin(); vit != vv.end() ; ++vit)
      vit->insert(vit->begin(), n);
    res.insert(res.end(), vv.begin(), vv.end());
    ener -= f;
    n++;
  }
  return res;
}

// \sum (e - \sum f_i * n_i)^p
double Model::vibrational_sum (double ener, const Lapack::Vector& freq, double power, int findex)
{
  const char funame [] = "Model::vibrational_sum: ";

  if(findex < 0 || findex > freq.size()) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

  if(ener <= 0.) {
    std::cerr << funame << "negative energy\n";
    throw Error::Range();
  }

  if(findex == freq.size())
    return std::pow(ener, power);

  double f = freq[findex];
  if(f <= 0.) {
    std::cerr << funame << findex << "-th frequency negative\n";
    throw Error::Range();
  }

  ++findex;
  double res = 0.;
  while(ener > 0.) {
    res += vibrational_sum(ener, freq, power, findex);
    ener -= f;
  }
  return res;
}

void Model::read_geometry (IO::KeyBufferStream& from, std::vector<Atom>& atom, int units)
{
  const char funame [] = "Model::read_geometry: ";

  int itemp;
  double dtemp;

  // number of atoms
  IO::LineInput size_input(from);
  size_input >> itemp;      
  if(!size_input) {
    std::cerr << funame << "cannot read number of atoms\n";
    throw Error::Input();
  }
  atom.resize(itemp);
  // read atoms
  for(std::vector<Atom>::iterator at = atom.begin(); at != atom.end(); ++at) {
    IO::LineInput data_input(from);
    if(!(data_input >> *at)) {
      std::cerr << funame << at - atom.begin() << "-th atom input failed\n";
      throw Error::Input();
    } 

    switch(units) {
    case ANGSTROM: 
      *at *= Phys_const::angstrom;
      break;
    case BOHR:
      break;
    default:
      std::cerr << funame << "unknown units\n";
      throw Error::Logic();
    }
  }

  // check interatomic distances
  for(int i = 0; i < atom.size(); ++i)
    for(int j = i + 1; j < atom.size(); ++j) {
      dtemp = vdistance(atom[i], atom[j]);
      if(dtemp < atom_dist_min) {
	std::cerr << funame << "the distance between the " 
		  << i + 1 << " and " << j + 1 << "-th atoms = " 
		  << dtemp << " bohr is less then " << atom_dist_min << " bohr: check the geometry\n";
	throw Error::Range();
      }
    }
  // shift center of mass to zero
  shift_cm_to_zero(atom);
}

/********************************************************************************************
 ******************************************** READERS ***************************************
 ********************************************************************************************/

SharedPointer<Model::Escape> Model::new_escape(IO::KeyBufferStream& from)
{
  const char funame [] = "Model::new_escape: ";

  KeyGroup NewEscapeModel;

  Key const_key("Constant");
  Key   fit_key("Fit");
  
  std::string token, comment;
  while (from >> token) {
    // energy independent model
    if(const_key == token) {
      std::getline(from, comment);
      return SharedPointer<Escape>(new ConstEscape(from)); 
    }   
    // spline fit
    if(fit_key == token) {
      std::getline(from, comment);
      return SharedPointer<Escape>(new FitEscape(from)); 
    }   
    // unknown keyword
    if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      IO::log << std::flush;
      throw Error::Init();
    }
  }  

  std::cerr << funame << "corrupted\n";
  IO::log << std::flush;
  throw Error::Input();
}


SharedPointer<Model::Collision> Model::new_collision(IO::KeyBufferStream& from)
{
  const char funame [] = "Model::new_collision: ";

  KeyGroup NewCollisionModel;

  Key lj_key("LennardJones");
  
  std::string token, comment;
  while (from >> token) {
    // Lennard-Jones model
    if(lj_key == token) {
      std::getline(from, comment);
      return SharedPointer<Collision>(new LennardJonesCollision(from)); 
    }   
    // unknown keyword
    if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      IO::log << std::flush;
      throw Error::Init();
    }
  }  

  std::cerr << funame << "corrupted\n";
  IO::log << std::flush;
  throw Error::Input();
}


SharedPointer<Model::Kernel> Model::new_kernel(IO::KeyBufferStream& from)
{
  const char funame [] = "Model::new_kernel: ";

  KeyGroup NewKernelModel;

  Key exp_key("Exponential");
  
  std::string token, comment;
  while(from >> token) {
    // single exponential model
    if(exp_key == token) {
      std::getline(from, comment);
      return SharedPointer<Kernel>(new ExponentialKernel(from));
    }
    // unknown keyword
    if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      IO::log << std::flush;
      throw Error::Init();
    }
  }

  std::cerr << funame << "corrupted\n";
  IO::log << std::flush;
  throw Error::Input();
}

SharedPointer<Model::Rotor> Model::new_rotor(IO::KeyBufferStream& from, const std::vector<Atom>& atom) 
{
  const char funame [] = "Model::new_rotor: ";

  KeyGroup NewRotorModel;

  Key hind_key("Hindered");
  Key free_key("Free"    );
  Key umbr_key("Umbrella");

  std::string token, comment;
  while(from >> token) {
    // hindered rotor model
    if(hind_key == token) {
      std::getline(from, comment);
      return SharedPointer<Rotor>(new HinderedRotor(from, atom));
    }
    // free rotor model
    if(free_key == token) {
      std::getline(from, comment);
      return SharedPointer<Rotor>(new FreeRotor(from, atom));
    }
    // umbrella mode model
    if(umbr_key == token) {
      std::getline(from, comment);
      return SharedPointer<Rotor>(new Umbrella(from, atom));
    }
  // unknown keyword
    if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      IO::log << std::flush;
      throw Error::Init();
    }
  }

  std::cerr << funame << "corrupted\n";
  IO::log << std::flush;
  throw Error::Input();
}

SharedPointer<Model::Tunnel> Model::new_tunnel(IO::KeyBufferStream& from) 
{
  const char funame [] = "Model::new_tunnel: ";

  KeyGroup NewTunnelModel;

  Key harmon_key("Harmonic");
  Key eckart_key("Eckart"  );
  Key  quart_key("Quartic" );
  Key   read_key("Read"    );

  std::string token, comment;
  while(from >> token) {
  // parabolic barrier tunneling
    if(harmon_key == token) {
      std::getline(from, comment);
      return SharedPointer<Tunnel>(new HarmonicTunnel(from));
    }
  // Eckart barrier tunneling
    if(eckart_key == token) {
      std::getline(from, comment);
      return SharedPointer<Tunnel>(new   EckartTunnel(from));
    }
    // Quartic potential barrier tunneling
    if(quart_key == token) {
      std::getline(from, comment);
      return SharedPointer<Tunnel>(new  QuarticTunnel(from));
    }
    if(read_key == token) {
      std::getline(from, comment);
      return SharedPointer<Tunnel>(new  ReadTunnel(from));
    }
    // unknown keyword
    if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      IO::log << std::flush;
      throw Error::Init();
    }
  }

  std::cerr << funame << "corrupted\n";
  IO::log << std::flush;
  throw Error::Input();  
}

SharedPointer<Model::Core> Model::new_core(IO::KeyBufferStream& from, const std::vector<Atom>& atom,  int mode) 
{
  const char funame [] = "Model::new_core: ";

  IO::Marker funame_marker(funame);
  
  KeyGroup NewCoreModel;

  Key rigid_key("RigidRotor"      );
  Key multi_key("MultiRotor"      );
  Key  rotd_key("Rotd"            );
  Key   pst_key("PhaseSpaceTheory");

  std::string token, comment;
  while(from >> token) {
  // rigid rotor
    if(rigid_key == token) {
      std::getline(from, comment);
      return SharedPointer<Core>(new RigidRotor(from, atom, mode));
    }
    // coupled internal rotations and vibrations
    if(multi_key == token) {
      std::getline(from, comment);
      return SharedPointer<Core>(new MultiRotor(from, atom, mode));
    }
    // transition modes number of states from Rotd
    if(rotd_key == token) {
      std::getline(from, comment);
      return SharedPointer<Core>(new Rotd(from, mode));
    }
    // phase space theory
    if(pst_key == token) {
      std::getline(from, comment);
      return SharedPointer<Core>(new PhaseSpaceTheory(from));
    }
    // unknown keyword
    if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      IO::log << std::flush;
      throw Error::Init();
    }
  }

  std::cerr << funame << "corrupted\n";
  IO::log << std::flush;
  throw Error::Input();
}

SharedPointer<Model::Species> Model::new_species(IO::KeyBufferStream& from, const std::string& name, int mode)
{
  const char funame [] = "Model::new_species: ";

  KeyGroup NewSpeciesModel;

  Key  rrho_key("RRHO"       );
  Key union_key("Union"      );
  Key   var_key("Variational");
  Key  read_key("Read"       );
  Key  atom_key("Atom"       );
  Key   mc_key("MonteCarlo"  );
  Key  mcd_key("MonteCarloWithDummyAtoms");
  
  std::string token, comment;

  while(from >> token)
    //
    if(IO::skip_comment(token, from))
      //
      break;
  
  if(!from) {
    //
    std::cerr << funame << "end of file\n";
    
    throw Error::Input();
  }
  
  // rigid-rotor-harmonic-oscillator model
  //
  if(rrho_key == token) {
    std::getline(from, comment);
    return SharedPointer<Species>(new RRHO(from, name, mode));
  }
  // read states from file
  //
  if(read_key == token) {
    std::getline(from, comment);
    return SharedPointer<Species>(new ReadSpecies(from, name, mode));
  }
  // union of species
  //
  if(union_key == token) {
    std::getline(from, comment);
    return SharedPointer<Species>(new UnionSpecies(from, name, mode));
  }
  // variational barrier model
  //
  if(var_key == token) {
    if(mode == DENSITY) {
      std::cerr << funame << token << ": wrong mode\n";
      throw Error::Init();
    }
    std::getline(from, comment);
    return SharedPointer<Species>(new VarBarrier(from, name));
  }
  // atomic fragment
  //
  if(atom_key == token) {
    if(mode != NOSTATES) {
      std::cerr << funame << token << ": wrong mode\n";
      throw Error::Init();
    }
    std::getline(from, comment);
    return SharedPointer<Species>(new AtomicSpecies(from, name));
  }
  // crude Monte Carlo
  //
  if(mc_key == token) {
    //
    if(mode != NOSTATES) {
      //
      std::cerr << funame << token << ": wrong mode\n";
      
      throw Error::Init();
    }
    
    std::getline(from, comment);
    
    return SharedPointer<Species>(new MonteCarlo(from, name, mode));
  }
  // crude Monte Carlo with Dummy atoms
  //
  if(mcd_key == token) {
    //
    if(mode != NOSTATES) {
      //
      std::cerr << funame << token << ": wrong mode\n";
      
      throw Error::Init();
    }
    
    std::getline(from, comment);
    
    return SharedPointer<Species>(new MonteCarloWithDummy(from, name, mode));
  }
  // no species
  //
  if(IO::end_key() == token) {
    std::getline(from, comment);
    return SharedPointer<Species>(0);
  }
  
  // unknown keyword
  //
  std::cerr << funame << "unknown keyword " << token << "\n";
  
  Key::show_all(std::cerr);

  std::cerr << "\n";
  
  IO::log << std::flush;
  
  throw Error::Init();
}

SharedPointer<Model::Bimolecular> Model::new_bimolecular(IO::KeyBufferStream& from, const std::string& name)
{
  return SharedPointer<Bimolecular>(new Bimolecular(from, name));
}


/********************************************************************************************
 ************************************* COLLISION MODEL **************************************
 ********************************************************************************************/

Model::Collision::~Collision ()
{
  //std::cout << "Model::Collision destroyed\n";
}

/********************************************************************************************
 **************************** LENNARD-JONES COLLISION MODEL *********************************
 ********************************************************************************************/

Model::LennardJonesCollision::LennardJonesCollision(IO::KeyBufferStream& from) 
  : _epsilon(-1.)
{
  const char funame [] = "Model::LennardJonesCollision::LennardJonesCollision: ";

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE | IO::Marker::NOTIME);

  KeyGroup LennardJonesCollisionModel;

  Key epsw_key("Epsilons[1/cm]"  );
  Key epsk_key("Epsilons[K]"     );
  Key  sig_key("Sigmas[angstrom]");
  Key mass_key("Masses[amu]"     );

  double  mass = -1.;
  double sigma = -1.;

  double dtemp, dtemp1, dtemp2;
  std::string token, comment;

  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    //  epsilons
    else if(epsw_key == token || epsk_key == token) {
      if(_epsilon > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> dtemp1 >> dtemp2)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp1 <= 0. || dtemp2 <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      _epsilon = std::sqrt(dtemp1 * dtemp2);
      // units
      if(epsw_key == token)
	_epsilon *= Phys_const::incm;
      else
	_epsilon *= Phys_const::kelv;

    }
    //  sigmas
    else if(sig_key == token) {
      if(sigma > 0.) {
	std::cerr << funame << token << ": alread initialized\n";
	throw Error::Init();
      }
      if(!(from >> dtemp1 >> dtemp2)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp1 <= 0. || dtemp2 <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      sigma = (dtemp1 + dtemp2) / 2. * Phys_const::angstrom;
    }
    //  masses
    else if(mass_key == token) {
      if(mass > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> dtemp1 >> dtemp2)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp1 <= 0. || dtemp2 <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      mass = Phys_const::amu / (1./dtemp1 + 1./dtemp2);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
  
  // checking
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_epsilon < 0.) {
    std::cerr << funame << "epsilons not initialized\n";
    throw Error::Init();
  }

  if(sigma < 0.) {
    std::cerr << funame << "sigmas not initialized\n";
    throw Error::Init();
  }

  if(mass < 0.) {
    std::cerr << funame << "masses not initialized\n";
    throw Error::Init();
  }

  _frequency_factor = std::sqrt(8. * M_PI / mass) * sigma * sigma;

}// Lennard-Jones Collision

Model::LennardJonesCollision::~LennardJonesCollision ()
{
  //std::cout << "Model::LennardJonesCollision destroyed\n";
}

double Model::LennardJonesCollision::_omega_22_star (double t) const
{
  t /= _epsilon;
  return 1.16145 / std::pow(t, 0.14874) 
    + 0.52487 / std::exp(0.7732 * t) 
    + 2.16178 / std::exp(2.437887  * t);

}

double Model::LennardJonesCollision::operator () (double t) const
{
  return _frequency_factor / std::sqrt(t) * _omega_22_star(t);
}


/********************************************************************************************
 ***************************** COLLISION RELAXATION KERNEL MODEL ****************************
 ********************************************************************************************/

Model::Kernel::~Kernel ()
{
  //std::cout << "Model::Kernel destroyed\n";
}

/********************************************************************************************
 *********************************** EXPONENTIAL KERNEL *************************************
 ********************************************************************************************/

Model::ExponentialKernel::ExponentialKernel(IO::KeyBufferStream& from) 
  : _cutoff(10.)
{
  const char funame [] = "Model::ExponentialKernel::ExponentialKernel: ";

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE | IO::Marker::NOTIME);

  KeyGroup ExponentialKernelModel;

  Key factor_key("Factor[1/cm]"  );
  Key  power_key("Power"         );
  Key  fract_key("Fraction"      );
  Key cutoff_key("ExponentCutoff");

  double dtemp;

  std::string token, comment;
  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // delta energy down factor
    else if(factor_key == token) {
      if(_factor.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      IO::LineInput lin(from);
      while(lin >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": out of range\n";
	  throw Error::Range();
	}
	_factor.push_back(dtemp * Phys_const::incm);
      }
    }
    // delta energy down power
    else if(power_key == token) {
      if(_power.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      IO::LineInput lin(from);
      while(lin >> dtemp) {
	if(dtemp < 0.) {
	  std::cerr << funame << token << ": out of range\n";
	  throw Error::Range();
	}
	_power.push_back(dtemp);
      }
    }
    // single exponential fraction
    else if(fract_key == token) {
      if(_fraction.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      IO::LineInput lin(from);
      while(lin >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": out of range\n";
	  throw Error::Range();
	}
	_fraction.push_back(dtemp);
      }
    }
    // energy cutoff parameter
    else if(cutoff_key == token) {
      if(!(from >> _cutoff)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_cutoff <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
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
  
  // checking
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_factor.size() == 1 && !_fraction.size())
    _fraction.push_back(1.);

  if(!_power.size() || _power.size() != _factor.size() || _power.size() != _fraction.size()) {
    std::cerr << funame << "initialization error\n";
    throw Error::Init();
  }

  dtemp = 0.;
  for(int i = 0; i < _fraction.size(); ++i)
    dtemp += _fraction[i];

  for(int i = 0; i < _fraction.size(); ++i)
    _fraction[i] = std::log(dtemp / _fraction[i]);


}// Exponential Kernel

Model::ExponentialKernel::~ExponentialKernel ()
{
  //std::cout << "Model::ExponentialKernel destroyed\n";
}

double Model::ExponentialKernel::_energy_down (int f, double temperature) const
{ 
  static const double normal_temperature = 300. * Phys_const::kelv;

  double res = _factor[f] * std::pow(temperature / normal_temperature, _power[f]);

  if(Kernel::flags() & UP) {
    return temperature * res / (temperature + res);
  }
  else
    return res;
}

double Model::ExponentialKernel::operator () (double energy, double temperature) const
{
  const char funame [] = "Model::ExponentialKernel::operator(): ";

  double dtemp;
 
  if(energy < 0.) {
    std::cerr << funame << "negative energy\n";
    throw Error::Range();
  }
  
  double res = 0.;

  for(int f = 0; f < _fraction.size(); ++f) {
    dtemp = energy / _energy_down(f, temperature) + _fraction[f];
    if(dtemp < _cutoff)
      res += std::exp(-dtemp);
  }

  return res;
}

double Model::ExponentialKernel::cutoff_energy (double temperature) const 
{ 
  double dtemp;

  double res = 0.;
  for(int f = 0; f < _fraction.size(); ++f) {
    dtemp = _energy_down(f, temperature) * (_cutoff - _fraction[f]);
    if(dtemp > res)
      res = dtemp; 
  }

  return res;
}

/********************************************************************************************
 ***************************************** TUNNELING ****************************************
 ********************************************************************************************/

Model::Tunnel::Tunnel (IO::KeyBufferStream& from) 
  : _cutoff(-1.), _wtol(1.e-2), _freq(-1.)
{
  const char funame [] = "Model::Tunnel::Tunnel : ";

  KeyGroup TunnelModel;

  Key kcal_cut_key("CutoffEnergy[kcal/mol]"    );
  Key incm_cut_key("CutoffEnergy[1/cm]"        );
  Key   kj_cut_key("CutoffEnergy[kJ/mol]"      );
  Key     wtol_key("StatisticalWeightTolerance");
  Key     freq_key("ImaginaryFrequency[1/cm]"  );

  std::string token, comment;
  while(from >> token) {
    // energy cutoff
    if(kcal_cut_key == token || incm_cut_key == token || kj_cut_key == token) {
      if(!(from >> _cutoff)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }    
      std::getline(from, comment);

      if(_cutoff < 0.) {
	std::cerr << funame << token << ": should not be negative\n";
      throw Error::Range();
      }
      if(kcal_cut_key == token)
	_cutoff *= Phys_const::kcal;
      if(incm_cut_key == token)
	_cutoff *= Phys_const::incm;
      if(kj_cut_key == token)
	_cutoff *= Phys_const::kjoul;
    }
    // imaginary frequency
    else if(freq_key == token) {
      if(_freq > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> _freq)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_freq <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      _freq *= Phys_const::incm;
    }
    // statistical weight tolerance
    else if(wtol_key == token) {
      if(!(from >> _wtol)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }    
      std::getline(from, comment);

      if(_wtol <= 0. || _wtol >= 1.) {
	std::cerr << funame << token << ": out of range\n";
      throw Error::Range();
      }
    }
    // unknown key
    else if(IO::skip_comment(token, from)) {
      //IO::last_key = token;
      from.put_back(token);
      break;
    }
  }

  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_freq < 0.) {
    std::cerr << funame << "imaginary frequency not initialized\n";
    throw Error::Init();
  }

}

Model::Tunnel::~Tunnel ()
{
  //std::cout << "Model::Tunnel destroyed\n";
}

// convolution of the number of states with the tunneling density
void Model::Tunnel::convolute(Array<double>& stat, double step) const
{
  const char funame [] = "Model::Tunnel::convolute: ";

  if(cutoff() < 0.) {
    std::cerr << funame << "cutoff energy is not initialized\n";
    throw Error::Init();
  }

  double dtemp;

  // tunneling density of states
  Array<double> td((int)std::ceil(2. * cutoff() / step) + 1);
  double ener = - cutoff();
  double fac = 0.;
  for(int i = 0; i < td.size(); ++i, ener += step) {
    dtemp = density(ener);
    td[i] = dtemp;
    fac  += dtemp;
  }

  td /= fac;

  Array<double> new_stat(stat.size());
  
#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic, 1)

  for(int j = 0; j < stat.size(); ++j) {
    dtemp = 0.;
    for(int i = 0; i < td.size(); ++i) {
      if(i > j)
	break;
      dtemp += stat[j - i] * td[i];
    }
    new_stat[j] = dtemp;
  }
  
  stat = new_stat;
}

// statistical weight relative to cutoff energy
double Model::Tunnel::weight (double temperature) const
{
  const char funame [] = "Model::Tunnel::weight: ";

  if(cutoff() < 0.) {
    std::cerr << funame << "cutoff energy is not initialized\n";
    throw Error::Init();
  }

  double dtemp;
  int    itemp;

  double estep = temperature * _wtol; // discretization energy step
  itemp = (int)std::floor(cutoff() / estep);
  int imax = 2 * itemp + 1;
  double ener = -double(itemp) * estep;

  double fac = 0.;
  double res = 0.;
  for(int i = 0; i < imax; ++i, ener += estep) {
    dtemp = density(ener);
    res += dtemp * std::exp(-ener / temperature);
    fac += dtemp;
  }

  res /= fac * std::exp(cutoff() / temperature);

  return res;
}

double Model::Tunnel::density (double ener) const
{
  double dtemp = action(ener, 0);
  if(dtemp > action_max() || dtemp < - action_max())
    return 0.;

  dtemp = std::cosh(dtemp / 2);
  return action(ener, 1) / 4. / dtemp / dtemp;
}

double Model::Tunnel::factor (double ener) const
{
  double dtemp = action(ener, 0);

  if(dtemp > action_max()) 
    return 1.;
  if(dtemp < - action_max())
    return 0.;
  return 1. / (1. + std::exp(-dtemp));
}

/**************************************************************************************
 ****************************** READ BARRIER TUNNELING ********************************
 **************************************************************************************/

Model::ReadTunnel::ReadTunnel (IO::KeyBufferStream& from) 
  : Tunnel(from)
{
  const char funame [] = "Model::ReadTunnel::ReadTunnel: ";

  int    itemp;
  double dtemp;

  IO::Marker funame_marker(funame);

  KeyGroup ReadTunnelModel;

  Key file_key("File");

  std::string file_name;

  std::string comment, token;

  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // data file name
    else if(file_key == token) {
      if(file_name.size()) {
	std::cerr << funame << token <<  ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> file_name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
 
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(!file_name.size()) {
    std::cerr << funame << "data file is not initialized\n";
    throw Error::Init();
  }

  /************************************* ACTION INPUT **********************************/

  std::ifstream fin(file_name.c_str());
  if(!fin) {
    std::cerr << funame << "cannot open rotd input file " << file_name << "\n";
    throw Error::Open();
  }

  std::map<double, double> ener_action_map;

  double ener;
  while(fin >> ener) {
    ener *= Phys_const::incm; // energy

    if(!(fin >> dtemp)) { // tunneling action
      std::cerr << funame << "reading tunneling action failed\n";
      throw Error::Input();
    }

    ener_action_map[ener] = -dtemp; // action is negative at negative energies
  }

  if(ener_action_map.begin()->first >= 0.) {
    std::cerr << funame << "no energies below the barrier\n";
    throw Error::Range();
  }

  if(ener_action_map.rbegin()->first < 0.)
    ener_action_map[0.] = 0.;

  Array<double> x((int)ener_action_map.size());
  Array<double> y((int)ener_action_map.size());

  itemp = 0;
  for(std::map<double, double>::const_iterator it = ener_action_map.begin(); it != ener_action_map.end(); ++it, ++itemp) {
    x[itemp] = it->first;
    y[itemp] = it->second;
  }

  _action.init(x, y, itemp);

  if(_action.arg_max() == 0.)
    IO::log << IO::log_offset << "no data for tunneling action above the barrier: mirrow immage will be used\n";
}

double Model::ReadTunnel::action (double ener, int n) const
{
  const char funame [] = "Model::ReadTunnel::action: ";

  int    itemp;
  double dtemp;

  double fac = 1.;
  if(ener > 0. && _action.arg_max() == 0.) {
    ener = -ener;
    fac = -1.;
  }
  
  if(ener < _action.arg_min()) {
    //
    dtemp = _action(_action.arg_min(), 1);
    
    switch(n) {
      //
    case 0:
      return fac * (_action.fun_min() + (ener - _action.arg_min()) * dtemp);

    case 1:
      return dtemp;

    default:
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }
  }
  else if(ener <= _action.arg_max()) {
    //
    switch(n) {
      //
    case 0:
      return fac * _action(ener, 0);

    case 1:
      return _action(ener, 1);

    default:
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }
  }
  else {
    //
    dtemp = _action(_action.arg_max(), 1);
    
    switch(n) {
      //
    case 0:
      return _action.fun_max() + (ener - _action.arg_max()) * dtemp;

    case 1:
      return dtemp;

    default:
      std::cerr << funame << "wrong case\n";
      throw Error::Logic();
    }
  }
}

Model::ReadTunnel::~ReadTunnel ()
{
}

/********************************************************************************************
 ******************************** PARABOLIC BARRIER TUNNELING *******************************
 ********************************************************************************************/

Model::HarmonicTunnel::HarmonicTunnel (IO::KeyBufferStream& from) 
  : Tunnel(from)
{
  const char funame [] = "Model::HarmonicTunnel::HarmonicTunnel: ";

  IO::Marker funame_marker(funame);

  KeyGroup HarmonicTunnelModel;

  std::string token, comment;
  //std::string token = IO::last_key, comment;
  //IO::last_key.clear();
  while(from >> token) {
    // end
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  } //while(from >> token);
  
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(cutoff() < 0.) {
    std::cerr << funame << "cutoff energy is not initialized\n";
    throw Error::Init();
  }

  IO::log << IO::log_offset << "cutoff energy = " << -cutoff() / Phys_const::kcal 
	  << " kcal/mol   tunneling factor = "
	  << factor(-cutoff()) << "\n";

}

Model::HarmonicTunnel::~HarmonicTunnel ()
{
  //std::cout << "Model::HarmonicTunnel destroyed\n";
}

double Model::HarmonicTunnel::action (double ener, int n) const
{
  const char funame [] = "Model::HarmonicTunnel::action: ";

  switch(n) {
  case 0: 
    return 2. * M_PI * ener / _freq;
  case 1: 
    return 2. * M_PI  / _freq;
  default:
    std::cerr << funame << "should not be here\n";
    throw Error::Logic();
  }
}

/********************************************************************************************
 **************************************** ECKART BARRIER ************************************
 ********************************************************************************************/

Model::EckartTunnel::EckartTunnel (IO::KeyBufferStream& from) 
  : Tunnel(from)
{
  const char funame [] = "Model::EckartTunnel::EckartTunnel: ";

#ifdef DEBUG

  IO::Marker funame_marker(funame);

#else

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE | IO::Marker::NOTIME);

#endif

  KeyGroup EckartTunnelModel;

  Key kcal_well_key("WellDepth[kcal/mol]");
  Key incm_well_key("WellDepth[1/cm]"    );
  Key   kj_well_key("WellDepth[kJ/mol]"  );

  double dtemp;
  int    itemp;

  std::string token, comment;
  //std::string token = IO::last_key, comment;
  //IO::last_key.clear();
  while(from >> token) {
    // end
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // well depth
    else if(kcal_well_key == token || incm_well_key == token || kj_well_key == token) {
      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      if(kcal_well_key == token)
	dtemp *= Phys_const::kcal;
      if(incm_well_key == token)
	dtemp *= Phys_const::incm;
      if(kj_well_key == token)
	dtemp *= Phys_const::kjoul;

      if(cutoff() >=0. && dtemp < cutoff()) {
	std::cerr << funame << token << ": should not be smaller than cutoff energy\n";
	throw Error::Range();
      }

      _depth.push_back(dtemp);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  } //while(from >> token);
  
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_depth.size() != 2) {
    std::cerr << funame << "wrong number of wells\n";
    throw Error::Init();
  }

  if(_depth[0] > _depth[1])
    std::swap(_depth[0], _depth[1]);

  // set up cutoff energy
  if(cutoff() < 0.)
    _cutoff = _depth[0];

  for(int w = 0; w < 2; ++w)
    _depth[w] /= _freq;

  dtemp = 0.;
  for(int w = 0; w < 2; ++w)
    dtemp += 1. / std::sqrt(_depth[w]);

  _factor = 4. * M_PI / dtemp;

#ifdef DEBUG

  IO::log << IO::log_offset << "Tunneling action S versus energy E[kcal/mol]:\n"
	  << IO::log_offset << std::setw(2) << "E" 
	  << std::setw(13) << "S_harm" 
	  << std::setw(13) << "S_full" 
	  << "\n";

  int emax = (int)std::ceil(cutoff() / Phys_const::kcal);
  for(int e = 1; e < emax; ++e)
    IO::log << IO::log_offset << std::setw(2) << e 
	    << std::setw(13) << 2. * M_PI * (double)e * Phys_const::kcal / _freq
	    << std::setw(13) << -action(-(double)e * Phys_const::kcal) << "\n";

  IO::log << IO::log_offset << "cutoff energy = " << -cutoff() / Phys_const::kcal 
	  << " kcal/mol   tunneling factor = "
	  << factor(-cutoff()) << "\n";

#endif

}

Model::EckartTunnel::~EckartTunnel ()
{
  //std::cout << "Model::EckartTunnel destroyed\n";
}

double Model::EckartTunnel::action (double ener, int der) const
{
  const char funame [] = "Model::EckartTunnel::action: ";
  
  double dtemp;
  double res = 0.;

  ener /= _freq;
  for(int w = 0; w < 2; ++w) {
    dtemp = ener + _depth[w];
    if(dtemp < 0.) 
      dtemp = dtemp < 0. ? 0. : dtemp;
    switch(der) {
    case 0:
      res += std::sqrt(dtemp) - std::sqrt(_depth[w]); 
      break;
    case 1:
      if(dtemp > 0.)
	res += 1./ std::sqrt(dtemp);
      break;
    default:
      std::cerr << funame << "should not be here\n";
      throw Error::Logic();
    }
  }
  res *= _factor;

  if(der)    
    res /= 2. * _freq;  

  return res;
}

  /**************************************************************************************
   ************************** QUARTIC BARRIER TUNNELING ********************************
   **************************************************************************************/

double Model::QuarticTunnel::XratioSearch::operator() (double x, int n) const
{
  const char funame [] = "Model::QuarticTunnel::Xratio::operator(): ";

  double dtemp;
  switch(n) {
  case 0:
    return (x + 2.) * x * x * x / (2. * x + 1.) - _vratio;
  case 1:
    dtemp = x * (x + 1.) / (2. * x + 1.);
    return 6. * dtemp * dtemp;
  default:
    std::cerr << funame << "wrong case\n";
    throw Error::Logic();
  }
}

Model::QuarticTunnel::QuarticTunnel(IO::KeyBufferStream& from) 
  : Tunnel(from)
{
  const char funame [] = "Model::QuarticTunnel::QuarticTunnel: ";

  double xtol   = 1.e-10;// potential minima coordinates ratio search tolerance
  double xstep  = 1.e-3; // coordinate discretization step
  int ener_size = 50;    // interpolation size
  std::vector<double> depth;

  IO::Marker funame_marker(funame);

  KeyGroup QuarticTunnelModel;

  Key kcal_well_key("WellDepth[kcal/mol]"         );
  Key incm_well_key("WellDepth[1/cm]"             );
  Key   kj_well_key("WellDepth[kJ/mol]"           );
  Key      xtol_key("CoordinateTolerance"         );
  Key      step_key("CoordinateDiscretizationStep");
  Key      ener_key("ActionInterpolationSize"     );

  double dtemp;
  int    itemp;

  std::string token, comment;
  //std::string token = IO::last_key, comment;
  //IO::last_key.clear();
  while(from >> token) {
    // end
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // well depth
    else if(kcal_well_key == token || incm_well_key == token || kj_well_key == token) {
      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      if(kcal_well_key == token)
	dtemp *= Phys_const::kcal;
      if(incm_well_key == token)
	dtemp *= Phys_const::incm;
      if(kj_well_key == token)
	dtemp *= Phys_const::kjoul;

      if(cutoff() >= 0. && dtemp < cutoff()) {
	std::cerr << funame << token << ": should be bigger than cutoff energy\n";
	throw Error::Range();
      }
      depth.push_back(dtemp);
    }
    // potential minima coordinates ratio search tolerance
    else if(xtol_key == token) {
      if(!(from >> xtol)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(xtol <= 0. || xtol >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();	
      }
    }
    // coordinate discretization step
    else if(step_key == token) {
      if(!(from >> xstep)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(xstep <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();	
      }
    }
    // energy interpolation size
    else if(ener_key == token) {
      if(!(from >> ener_size)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(ener_size <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();	
      }
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  } //while(from >> token);
  
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  double xmin, xmax;
  switch(depth.size()) {
  case 1: // cubic potential
    IO::log << IO::log_offset << "assuming cubic potential\n";
    xmin = std::sqrt(6.);
    _v3 = -1. / xmin / 3.;
    _v4 = 0.;    
    break;
  case 2: // quartic potential
    if(depth[0] > depth[1])
      std::swap(depth[0], depth[1]);

    // ratio of potential minima 
    dtemp = depth[1] / depth[0];
    // ratio of potential minima coordinates
    xmax = std::pow(dtemp, 1./3.);// initial guess
    XratioSearch(dtemp, xtol).find(xmax);

    // potential minimal coordinates 
    xmin = std::sqrt(12. * xmax / (2. * xmax + 1.));
    xmax *= xmin;

    // potential coefficients
    _v3 = (1./xmax - 1./xmin) / 3.;
    _v4 = -0.25 / xmax / xmin;

    IO::log << IO::log_offset << "small reduced potential minimum         = " << _potential(xmin)   << "\n"
	    << IO::log_offset << "large reduced potential minimum         = " << _potential(-xmax)  << "\n"
	    << IO::log_offset << "real potential minima ratio             = " <<  depth[1]/depth[0] << "\n";
    break;
  default:
    std::cerr << funame << "wrong number of wells\n";
    throw Error::Init();
  }

  // lowest tunneling energy
  _vmin = depth[0];

  // setting up cutoff energy
  if(cutoff() < 0.)
    _cutoff = _vmin;

  // action quadratic correction coefficient in terms of potential anharmonicities
  double quad_corr = 3.75 * _v3 * _v3 - 1.5 * _v4;
  IO::log << IO::log_offset << "quadratic energy term coefficient (in action) = " <<  quad_corr << "\n";

  // coordinate discretization size and step
  itemp = (int)std::ceil(xmin / xstep);
  xstep = xmin / (double)itemp;
  int xsize = 2 * itemp + 1;

  // energy interpolation step
  double ener_step = 1. / (double)ener_size;

  // energy and action on the grid
  Array<double>   ener_grid(ener_size + 1);
  Array<double> action_grid(ener_size + 1);
  ener_grid[0]   = 0.;
  action_grid[0] = 0.;

  double x;
  double ener = 1.;
  for(int e = ener_size; e > 0; --e, ener -= ener_step) {
    ener_grid[e] = ener;
    x = xmin;
    action_grid[e] = 0.;
    for(int i = 0; i < xsize; ++i, x -= xstep) {
      dtemp = ener - _potential(x);
      if(dtemp > 0.)
	action_grid[e] += std::sqrt(dtemp);
    }
    action_grid[e] *= xstep * M_SQRT2 * 2.;
  }
  
  // action spline interpolation
  _action.init(ener_grid, action_grid, ener_grid.size());
  
  IO::log << IO::log_offset << "Tunneling action versus energy E[kcal/mol]:\n"
	  << IO::log_offset << std::setw(2) << "E" 
	  << std::setw(13) << "S_harm" 
	  << std::setw(13) << "S_corr" 
	  << std::setw(13) << "S_full" 
	  << "\n";

  int emax = (int)std::ceil(cutoff() / Phys_const::kcal);
  for(int e = 1; e < emax; ++e) {
    ener = (double)e * Phys_const::kcal;
    IO::log << IO::log_offset << std::setw(2) << e 
	    << std::setw(13) << 2. * M_PI * ener / _freq 
	    << std::setw(13) << 2. * M_PI * ener / _freq * (1. + quad_corr * ener / _vmin)
	    << std::setw(13) << -action(-ener) 
	    << "\n";
  }

  IO::log << IO::log_offset << "cutoff energy = " << -cutoff() / Phys_const::kcal 
	  << " kcal/mol   tunneling factor = "
	  << factor(-cutoff()) << "\n";

}

Model::QuarticTunnel::~QuarticTunnel ()
{
  //std::cout << "Model::QuarticTunnel destroyed\n";
}

double Model::QuarticTunnel::action (double ener, int der) const 
{
  const char funame [] = "Model::QuarticTunnel::action: ";
  
  double dtemp;
  if(ener >= 0.)
    switch(der) {
    case 0:
      return 2. * M_PI * ener / _freq;
    case 1:
      return 2. * M_PI / _freq;
    default:
      std::cerr << funame << "wrong derivative\n";
      throw Error::Range();
    }

  ener /= -_vmin;

  if(ener < 1.)
    switch(der) {
    case 0:
      return - _vmin * _action(ener, 0) / _freq;
    case 1:
      return _action(ener, 1) / _freq;
    default:
      std::cerr << funame << "wrong derivative\n";
      throw Error::Range();
    }
  else
    switch(der) {
    case 0:
      return - _vmin * _action.fun_max() / _freq;
    case 1:
      return 0.;
    default:
      std::cerr << funame << "wrong derivative\n";
      throw Error::Logic();
    }
}

/********************************************************************************************
 ******************************** INTERNAL ROTATION DEFINITION ******************************
 ********************************************************************************************/

Model::InternalRotationBase::InternalRotationBase (IO::KeyBufferStream& from) 
  : _axis(-1, -1), _symmetry(1), _isinit(true)
{
  const char funame [] = "Model::InternalRotationBase::InternalRotationBase: ";

  KeyGroup InternalRotationBaseDefinition;

  Key group_key("Group"   );
  Key  axis_key("Axis"    );  
  Key  symm_key("Symmetry");

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  std::string token, line, comment;
  //if(IO::last_key.size()) {
  //  token = IO::last_key;
  //  IO::last_key.clear();
  //}
  //else if(!(from >> token)) {
  //  std::cerr << funame << "stream is corrupted\n";
  //  throw Error::Input();
  //}  
  
  while(from >> token) {
    // moving group definition
    if(group_key == token) {
      if(_group.size()) {
	std::cerr << funame << token << ": group is already defined\n";
	throw Error::Init();
      }
      IO::LineInput group_input(from);
      while(group_input >> itemp) {
	// using Fortran type indexing in the input
	if(itemp < 1) {
	  std::cerr << funame << token << ": atomic index should be positive\n";
	  throw Error::Range();
	}
	if(!_group.insert(itemp - 1).second) {
	  std::cerr << funame << token << ": " << itemp 
		    << "-th atom already has been used in the internal motion defininition\n";
	  throw Error::Init();
	}
      }
    }
    // internal rotation axis definition
    else if(axis_key == token) {
      if(_axis.first >= 0) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }

      IO::LineInput axis_input(from);
      for(int i = 0; i < 2; ++i) {
	if(!(axis_input >> itemp)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}

	if(itemp < 1) {// using Fortran type indexing in the input
	  std::cerr << funame << token << ": atomic index should be positive\n";
	  throw Error::Range();
	}

	--itemp;

	if(i)
	  _axis.second = itemp;
	else
	  _axis.first  = itemp;
      }
    
      if(_axis.first == _axis.second) {
	std::cerr << funame << token << ": atomic indices should be different\n";
	throw Error::Range();
      }
    }
    // symmetry
    else if(symm_key == token) {
      if(!(from >> _symmetry)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_symmetry < 1) {
	std::cerr << funame << token << ": symmetry should be positive\n";
	throw Error::Range();
      }
    }
    // unknown key
    else if(IO::skip_comment(token, from)) {
      from.put_back(token);
      //IO::last_key = token;
      break;
    }
  } // while(from >> token);

  if(!from) {
    std::cerr << funame << "input is corrupted\n";
    throw Error::Input();
  }

  if(!_group.size()) {
    std::cerr << funame  << "rotor group not defined\n";
    throw Error::Init();
  }

  if(_axis.first < 0) {
    std::cerr << funame  << "axis not defined\n";
    throw Error::Init();
  }


  if(_group.find(_axis.first ) != _group.end() ||
     _group.find(_axis.second) != _group.end()) {
    std::cerr << funame << "rotor group should not have common atoms with the axis\n";
    throw Error::Init();
  }

  // maximal atomic index
  _imax = *_group.rbegin();
  _imax = _axis.first  > _imax ? _axis.first  : _imax;
  _imax = _axis.second > _imax ? _axis.second : _imax;

}// Internal Rotation Base

Model::InternalRotationBase::~InternalRotationBase ()
{
  //IO::log << IO::log_offset << "Model::InternalRotationBase destroyed\n";
}

// rotating internal group around internal axis 
std::vector<Atom> Model::InternalRotationBase::rotate (const std::vector<Atom>& atom, double angle) const
{
  const char funame [] = "Model::InternalRotationBase::rotate: ";

  if(!_isinit) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(_imax >= atom.size() || _group.size() + 2 >= atom.size()) {
    std::cerr << funame << "atomic indices are inconsistent with the number of atoms\n";
    throw Error::Logic();
  }

  std::vector<Atom> res(atom);

  int        itemp;
  double     dtemp;
  D3::Vector vtemp;

  D3::Vector ref = atom[_axis.first];
  D3::Vector nz  = atom[_axis.second] - ref;
  nz.normalize();

  D3::Vector nx, ny;

  // new geometry
  std::set<int>::const_iterator cit;
  for(cit = _group.begin(); cit != _group.end(); ++cit) {
    vtemp = res[*cit] - ref;
    ny = D3::vprod(nz, vtemp);
    dtemp = ny.normalize();
    nx = D3::vprod(ny, nz);
    nx *= dtemp * std::cos(angle);
    ny *= dtemp * std::sin(angle);
    dtemp = vdot(nz, vtemp);
    vtemp = dtemp * nz;
    res[*cit] = ref + nx + ny + vtemp;
  }
  return res;
}

// center of mass is assumed to be at zero
std::vector<D3::Vector> Model::InternalRotationBase::normal_mode 
(const std::vector<Atom>& atom, Lapack::Vector* avp) const
{
  const char funame [] = "Model::InternalRotationBase::normal_mode: ";

  if(!_isinit) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(_imax >= atom.size() || _group.size() + 2 >= atom.size()) {
    std::cerr << funame << "atomic indices are inconsistent with the number of atoms\n";
    throw Error::Logic();
  }

  int        itemp;
  double     dtemp;
  D3::Vector vtemp;

  // rotation axis
  D3::Vector dir = atom[_axis.second] - atom[_axis.first];
  dir.normalize();

  // rotational motion
  std::map<int, D3::Vector> motion;
  for(std::set<int>::const_iterator cit = _group.begin(); cit != _group.end(); ++cit)
    motion[*cit] = D3::vprod(dir, atom[*cit] - atom[_axis.first]);

  // normal mode 
  return adjusted_normal_mode(atom, motion, avp);
}

/********************************************************************************************
 *********************** 1D UMBRELLA MODE/ ROTOR CLASS FAMILY  ******************************
 ********************************************************************************************/

Model::Rotor::Rotor () :_ham_size_min(99), _ham_size_max(999), _grid_size(1000), _therm_pow_max(10.) {}

Model::Rotor::Rotor (IO::KeyBufferStream& from, const std::vector<Atom>& atom) 
  :_ham_size_min(99), _ham_size_max(999), _grid_size(1000), _therm_pow_max(10.), _atom(atom)
{
  const char funame [] = "Model::Rotor::Rotor: ";

  KeyGroup RotorModel;

  Key ang_geom_key("Geometry[angstrom]");
  Key bor_geom_key("Geometry[au]"      );
  Key     hmax_key("HamiltonSizeMax"   );
  Key     hmin_key("HamiltonSizeMin"   );
  Key     grid_key("GridSize"          );
  Key     pmax_key("ThermalPowerMax"   );

  int         itemp;
  double      dtemp;
  std::string stemp;

  std::string token, line, comment;
  while(from >> token) {
    //
    // hindered rotor geometry
    //
    if(ang_geom_key == token || bor_geom_key == token) {
      //
      if(_atom.size()) {
	//
	IO::log << IO::log_offset << "WARNING: geometry has been changed\n";

	_atom.clear();
      }
      
      if(ang_geom_key == token)
	//
	read_geometry(from, _atom, ANGSTROM);
      
      if(bor_geom_key == token)
	//
	read_geometry(from, _atom, BOHR);
    }
    // Hamiltonian maximum size
    //
    else if(hmax_key == token) {
      //
      if(!(from >> _ham_size_max)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }
      
      std::getline(from, comment);

      if(_ham_size_max < 1) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }
      
      if(!(_ham_size_max % 2))
	//
	++_ham_size_max;
    }
    // Hamiltonian minimum size
    //
    else if(hmin_key == token) {
      if(!(from >> _ham_size_min)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_ham_size_min < 1) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      if(!(_ham_size_min % 2))
	++_ham_size_min;
    }
    // potential discretization size
    else if(grid_key == token) {
      if(!(from >> _grid_size)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_grid_size < 1) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
    }
    // thermal exponent power maximum
    else if(pmax_key == token) {
      if(!(from >> _therm_pow_max)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_therm_pow_max <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();	
      }
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      from.put_back(token);
      //IO::last_key = token;
      break;
    }
  } //while(from >> token);
  
  // stream state checking
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  // geometry checking
  if(!_atom.size()) {
    std::cerr << funame << "geometry is not defined\n";
    throw Error::Init();
  }
}

Model::Rotor::~Rotor ()
{
  //std::cout << "Model::Rotor destroyed\n";
}

void Model::Rotor::convolute(Array<double>& stat_grid, double ener_quant) const
{
  const char funame [] = "Model::Rotor::convolute: ";
  
  int    itemp;
  double dtemp;

  Array<double> new_stat_grid = stat_grid;

  std::map<int, int> shift;
  
  for(int n = 1; n < level_size(); ++n) {
    //
    dtemp = std::ceil(energy_level(n) / ener_quant);

    if(dtemp < 0.) {
      //
      std::cerr << funame << n << "-th energy level is negative: "
		<< energy_level(n) / Phys_const::kcal << " kcal/mol\n";

      throw Error::Range();
    }
    
    if(dtemp < (double)stat_grid.size()) {
      //
      itemp = (int)dtemp;
      
      shift[itemp]++;
    }
    else {
      //
      IO::log << IO::log_offset << funame << "WARNING: "
	      << n << "-th energy level is too high: "
	      << energy_level(n) / Phys_const::kcal << " kcal/mol\n";
    }
  }

#pragma omp parallel for default(shared) schedule(static, 10)

  for(int i = 0; i < stat_grid.size(); ++i)
    //
    for(std::map<int, int>::const_iterator it = shift.begin(); it != shift.end(); ++it)
      //
      if(i >= it->first)
	//
	new_stat_grid[i] += (double)it->second * stat_grid[i - it->first];

  stat_grid = new_stat_grid; 
}

/********************************************************************************************
 ********************************* 1-D HINDERED ROTOR MODEL *********************************
 ********************************************************************************************/

Model::RotorBase::RotorBase (IO::KeyBufferStream& from, const std::vector<Atom>& atom) 
  :  Rotor(from, atom), InternalRotationBase(from)
{
  const char funame [] = "Model::RotorBase::RotorBase: ";

  double dtemp;

  // normal mode
  //
  std::vector<D3::Vector> nm = normal_mode(_atom);
  
  // effective mass
  //
  dtemp = 0.;
  
  for(int a = 0; a < _atom.size(); ++a)
    //
    dtemp += _atom[a].mass() * vdot(nm[a]);

  // rotational constant
  //
  _rotational_constant = 0.5 / dtemp;

}// 1D Rotor base

Model::RotorBase::~RotorBase ()
{
  //std::cout << "Model::RotorBase destroyed\n";
}

/********************************************************************************************
 **************************************** FREE ROTOR ****************************************
 ********************************************************************************************/

Model::FreeRotor::FreeRotor (IO::KeyBufferStream& from, const std::vector<Atom>& atom) 
  : RotorBase(from, atom), _level_size(1)
{
  const char funame [] = "Model::FreeRotor::FreeRotor: ";

  IO::Marker funame_marker(funame);

  IO::log << IO::log_offset << "effective rotational constant = "
	    << rotational_constant() / Phys_const::incm << " 1/cm\n";
      
  KeyGroup FreeRotorModel;

  int    itemp;
  double dtemp;

  std::string token, line, comment;
  //std::string token = IO::last_key, line, comment;
  //IO::last_key.clear();
  while(from >> token) {
    // end input
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  } //while(from >> token);
  
  /******************************************* Checking *************************************************/
  // stream state
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

}// Free 1D rotor

Model::FreeRotor::~FreeRotor ()
{
  //std::cout << "Model::FreeRotor destroyed\n";
}

double Model::FreeRotor::ground () const
{
  return 0.;
}

double Model::FreeRotor::energy_level (int level) const
{
  if(!level)
    return 0.;

  int itemp = (level + 1) / 2 * symmetry();
  return rotational_constant() * double(itemp * itemp);
}

int Model::FreeRotor::level_size () const
{
  return _level_size;
}

double Model::FreeRotor::weight (double temperature) const
{
  /*
  static const double max_exp_pow = 35.;

  double dtemp;

  double res = 1.;
  for(int l = 1; l < level_size(); ++l) {
    dtemp = energy_level(l) / temperature;
    if(dtemp > max_exp_pow)
      return res;
    res += std::exp(-dtemp);
  }

  return res;
  */

  return std::sqrt(M_PI * temperature / rotational_constant()) / (double)symmetry();
}

void  Model::FreeRotor::set (double ener_max)
{
  if(ener_max <= 0.) {
    _level_size = 1;
    return;
  }

  _level_size = int(std::sqrt(ener_max / rotational_constant())) / symmetry() * 2 + 1;
  return;
}

/********************************************************************************************
 ********************************** QUANTUM HINDERED ROTOR **********************************
 ********************************************************************************************/

Model::HinderedRotor::HinderedRotor (const std::map<int, double>& p, double r, int s)
  : RotorBase(r, s), _pot_four(p), _weight_output_temperature_step(100), _weight_output_temperature_max(1000),
    _weight_output_temperature_min(-1), _use_quantum_weight(false)
{
  const char funame [] = "Model::HinderedRotor::HinderedRotor: ";

  IO::Marker funame_marker(funame);

  _init();
}

Model::HinderedRotor::HinderedRotor (IO::KeyBufferStream& from, const std::vector<Atom>& atom) 
  : RotorBase(from, atom), _weight_output_temperature_step(100), _weight_output_temperature_max(1000),
    _weight_output_temperature_min(-1), _use_quantum_weight(false)
{
  const char funame [] = "Model::HinderedRotor::HinderedRotor: ";

  IO::Marker funame_marker(funame);

  _read(from);
  _init();
}

void Model::HinderedRotor::_read(IO::KeyBufferStream& from) 
{
  const char funame [] = "Model::HinderedRotor::_read: ";

  int    itemp;
  double dtemp;

  IO::Marker funame_marker(funame);

  KeyGroup HinderedRotorModel;

  Key  kcal_pot_key("Potential[kcal/mol]"       );
  Key  incm_pot_key("Potential[1/cm]"           );
  Key    kj_pot_key("Potential[kJ/mol]"         );
  Key kcal_four_key("FourierExpansion[kcal/mol]");
  Key incm_four_key("FourierExpansion[1/cm]"    );
  Key   kj_four_key("FourierExpansion[kJ/mol]"  );
  Key     tstep_key("OutputTemperatureStep[K]"  );
  Key      tmax_key("OutputTemperatureMax[K]"   );
  Key      tmin_key("OutputTemperatureMin[K]"   );
  Key  use_quan_key("UseQuantumWeight"          );

  std::string token, line, comment;
  
  while(from >> token) {
    //
    // input end
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);

      break;
    }
    // potential on the grid
    //
    else if(kcal_pot_key == token || incm_pot_key == token || kj_pot_key == token) {
      //
      if(_pot_four.size()) {
	//
	std::cerr << funame << "potential fourier expansion has been already defined\n";

	throw Error::Init();
      }

      // potential sampling size
      //
      IO::LineInput size_input(from);
      
      if(!(size_input >> itemp)) {
	//
	std::cerr << funame << token << ": potential sampling size unreadable\n";

	throw Error::Input();
      }
      
      if(itemp < 2) {
	//
	std::cerr << funame << token << ": potential sampling size = " << itemp << " too small\n";

	throw Error::Range();
      }

      // potential sampling
      //
      std::vector<double> pval(itemp);

      for(int i = 0; i < pval.size(); ++i) {
	//
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": cannot read potential\n";

	  throw Error::Input();
	}

	if(kcal_pot_key == token)
	  //
	  dtemp *= Phys_const::kcal;
	
	if(incm_pot_key == token)
	  //
	  dtemp *= Phys_const::incm;

	if(kj_pot_key == token)
	  //
	  dtemp *= Phys_const::kjoul;

	pval[i] = dtemp;
      }
      
      std::getline(from, comment);
      
      // fourier transform
      //
      itemp = pval.size() % 2 ? pval.size() : pval.size() + 1;
      
      for(int i = 0; i < itemp; ++i)
	//
	for(int j = 0; j < pval.size(); ++j)
	  //
	  if(i % 2) {
	    //
	    _pot_four[i] += pval[j] * std::sin(M_PI * double((i + 1)  * j) / double(pval.size()));
	  }
	  else
	    //
	    _pot_four[i] += pval[j] * std::cos(M_PI * double(i * j) / double(pval.size()));

      // normalization
      //
      for(int i = 0; i < _pot_four.size(); ++i)
	//
	if(!i || i == pval.size()) {
	  //
	  _pot_four[i] /= double(pval.size());
	}
	else
	  //
	  _pot_four[i] /= double(pval.size()) / 2.;

      IO::log << IO::log_offset << "Fourier Expansion Coefficients(kcal/mol):\n";
      
      for(int i = 0; i < _pot_four.size(); ++i)
	//
	IO::log << IO::log_offset << std::setw(3) << i << std::setw(15) <<  _pot_four[i] / Phys_const::kcal << "\n";

#ifdef DEBUG

      IO::log << IO::log_offset << funame << "Hindered Rotor Potential:\n";
      
      IO::log << IO::log_offset 
	      << std::setw(18) << "Angle[degrees]" 
	      << std::setw(18) << "Energy[kcal/mol]"
	      << "\n";
      
      for(int i = 0; i < pval.size(); ++i) {
	//
	dtemp = double(i) * 360. / double(pval.size() * symmetry());
	
	IO::log << IO::log_offset 
		<< std::setw(18) << dtemp
		<< std::setw(18) << potential(dtemp * M_PI / 180.) / Phys_const::kcal
		<< "\n";
      }

#endif
      
    }
    // potential fourier expansion
    //
    else if(kcal_four_key == token || incm_four_key == token || kj_four_key == token) {
      //
      if(_pot_four.size()) {
	//
	std::cerr << funame << "potential fourier expansion has been already defined\n";

	throw Error::Init();
      }
      
      // fourier expansion size
      //
      IO::LineInput size_input(from);
      
      int esize;
      
      if(!(size_input >> esize)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      if(esize < 1) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }
      
      // potential fourier expansion data
      //
      for(int i = 0; i < esize; ++i) {
	//
	IO::LineInput data_input(from);
	
	if(!(data_input >> itemp)) {
	  //
	  std::cerr << funame << token << i << "-th index is unreadable\n";
	  
	  throw Error::Input();
	}
	
	if(itemp < 0) {
	  //
	  std::cerr << funame << token << ": negative index\n";
	  
	  throw Error::Init();
	}
	
	if(_pot_four.find(itemp) != _pot_four.end()) {
	  //
	  std::cerr << funame << token << ": identical indices\n";
	  
	  throw Error::Init();
	}
	
	if(!(data_input >> dtemp)) {
	  //
	  std::cerr << funame << token << i << "-th fourier coeficient is unreadable\n";

	  throw Error::Input();
	}

	if(kcal_four_key == token)
	  //
	  dtemp *= Phys_const::kcal;
	
	if(incm_four_key == token)
	  //
	  dtemp *= Phys_const::incm;
	
	if(kj_four_key == token)
	  //
	  dtemp *= Phys_const::kjoul;

	_pot_four[itemp] = dtemp;
      }
    }
    // output temperature step
    //
    else if(tstep_key == token) {
      //
      if(!(from >> itemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      std::getline(from, comment);

      if(itemp <= 0) {
	//
	std::cerr << funame << token << ": should be positive\n";
	
	throw Error::Range();
      }
      
      _weight_output_temperature_step = itemp;
    }
    // output temperature max
    //
    else if(tmax_key == token) {
      //
      if(!(from >> itemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }
      
      std::getline(from, comment);

      if(itemp <= 0) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }
      
      _weight_output_temperature_max = itemp;
    }
    // output temperature min
    //
    else if(tmin_key == token) {
      //
      if(!(from >> itemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      std::getline(from, comment);

      if(itemp <= 0) {
	//
	std::cerr << funame << token << ": should be positive\n";
	
	throw Error::Range();
      }
      
      _weight_output_temperature_min = itemp;
    }
    // use quantum weight
    //
    else if(use_quan_key == token) {
      //
      std::getline(from, comment);
      
      _use_quantum_weight = true;
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
    //
  }// input cycle

  // stream state checking
  //
  if(!from) {
    //
    std::cerr << funame << "input stream corrupted\n";
    
    throw Error::Input();
  }
}

void Model::HinderedRotor::_init ()
{
  const char funame [] = "Model::HinderedRotor::_init: ";

  int    itemp;
  double dtemp;

  // fourier expansion checking
  if(!_pot_four.size()) {
    std::cerr << funame << "fourier expansion is not initialized\n";
    throw Error::Init();
  }
 
  // hamiltonian size adjusting
  itemp = _pot_four.rbegin()->first;
  if(itemp % 2)
    itemp += 4;
  else
    itemp += 3;

  if(itemp > _ham_size_min)
    _ham_size_min = itemp;
      

  if(_ham_size_max < _ham_size_min)
    _ham_size_max = _ham_size_min;

  // minimum and maximum potential energies and potential discretization
  _grid_step = 2. * M_PI / double(symmetry() * _grid_size);
  _pot_grid.resize(_grid_size);
  _freq_grid.resize(_grid_size);
  // potential
  _pot_grid[0] = _pot_min = _pot_max = potential(0.);
  // frequency
  dtemp = potential(0., 2);
  if(dtemp < 0.)
    dtemp = -std::sqrt(-2. * dtemp * rotational_constant());
  else
    dtemp = std::sqrt(2. * dtemp * rotational_constant());
  _freq_grid[0] = _freq_min = _freq_max = dtemp;

  int imax, imin;
  imin = imax = 0;
  double a = _grid_step;
  for(int i = 1; i < _grid_size; ++i, a += _grid_step) {
    // potential
    dtemp = potential(a);
    _pot_grid[i] = dtemp; 
    if(dtemp < _pot_min) {
      _pot_min = dtemp;
      imin = i;
    }
    if(dtemp > _pot_max) {
      _pot_max = dtemp;
      imax = i;
    }
    // frequency
    dtemp = potential(a, 2);
    if(dtemp < 0.)
      dtemp = -std::sqrt(-2. * dtemp * rotational_constant());
    else
      dtemp = std::sqrt(2. * dtemp * rotational_constant());
    _freq_grid[i] = dtemp;
    if(dtemp < _freq_min) {
      _freq_min = dtemp;
    }
    if(dtemp > _freq_max) {
      _freq_max = dtemp;
    }
  }

  // maximum energy correction
  dtemp = _pot_grid[imax == _pot_grid.size() - 1 ? 0: imax + 1] -
    _pot_grid[imax == 0 ? _pot_grid.size() - 1 : imax - 1];
  
  a = _pot_grid[imax == _pot_grid.size() - 1 ? 0: imax + 1] +
    _pot_grid[imax == 0 ? _pot_grid.size() - 1 : imax - 1] - 2. * _pot_grid[imax];

  if(a < -1.e-5)
    _pot_max -= dtemp * dtemp / a / 8.;

  // minimum energy correction
  dtemp = _pot_grid[imin == _pot_grid.size() - 1 ? 0: imin + 1] -
    _pot_grid[imin == 0 ? _pot_grid.size() - 1 : imin - 1];
  
  a = _pot_grid[imin == _pot_grid.size() - 1 ? 0: imin + 1] +
    _pot_grid[imin == 0 ? _pot_grid.size() - 1 : imin - 1] - 2. * _pot_grid[imin];

  if(a > 1.e-5)
    _pot_min -= dtemp * dtemp / a / 8.;

  // numerical vibrational frequency
  double freq_num = 0.;
  if(a > 0.)
    freq_num = std::sqrt(a * 2. * rotational_constant()) / _grid_step;

  // analytical vibrational frequency at the minimum
  double freq_anal;
  dtemp = potential((double)imin * _grid_step, 2);
  if(dtemp > 0.)
    freq_anal = std::sqrt(dtemp * 2. * rotational_constant());
  else {
    IO::log << IO::log_offset << funame 
	    << "WARNING: potential second derivative at the minimum is negative, assuming zero\n";
    freq_anal = 0.;
  }

  _harm_ground = _pot_min + freq_anal / 2.;
  
  IO::log << IO::log_offset << "effective rotational constant[1/cm]  = " 
	  << rotational_constant() / Phys_const::incm << "\n";
  IO::log << IO::log_offset << "analytic  frequency at minimum[1/cm] = " 
	  << freq_anal / Phys_const::incm << "\n";
  IO::log << IO::log_offset << "numerical frequency at minimum[1/cm] = " 
	  << freq_num  / Phys_const::incm << "\n";    
  IO::log << IO::log_offset << "minimum energy[kcal/mol]             = "
	  << _pot_min  / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "maximum energy[kcal/mol]             = "
	  << _pot_max  / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "maximum frequency[1/cm]              = "
	  << _freq_max / Phys_const::incm << "\n";    
  IO::log << IO::log_offset << "maximum imaginary frequency[1/cm]    = "
	  <<-_freq_min / Phys_const::incm << "\n";    

  // calculating energy levels in momentum space
  _set_energy_levels(_ham_size_min);

  IO::log << IO::log_offset << "ground energy [kcal/mol]             = "
	  << ground() / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "highest energy level [kcal/mol]      = "
	  << (_energy_level.back() + ground()) / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "number of levels                     = "
	  << level_size()  << "\n";

  // energy levels output
  itemp = level_size() < 10 ? level_size() : 10;
  IO::log << IO::log_offset << itemp << " lowest excited states [kcal/mol] relative to the ground:" 
	  << std::setprecision(3);
  for(int l = 1; l < itemp; ++l)
    IO::log  << " " << energy_level(l) / Phys_const::kcal;
  IO::log << std::setprecision(6) << "\n";

  // statistical weight
  IO::log << IO::log_offset << "Statistical Weight (*** - deep tunneling regime):\n";
  IO::log << IO::log_offset 
	  << std::setw(5) << "T, K" 
	  << std::setw(15) << "Quantum"
	  << std::setw(15) << "Classical"
	  << std::setw(15) << "Semiclassical"
	  << "  ***\n";

  if(_weight_output_temperature_min < 0)
    _weight_output_temperature_min = _weight_output_temperature_step;

  for(int t = _weight_output_temperature_min; t <= _weight_output_temperature_max; t += _weight_output_temperature_step) {
    double tval = (double)t * Phys_const::kelv;
    double cw, sw;
    itemp = get_semiclassical_weight(tval, cw, sw);
    IO::log << IO::log_offset 
	    << std::setw(5) << t
	    << std::setw(15) << quantum_weight(tval)
	    << std::setw(15) << cw
	    << std::setw(15) << sw;
    if(itemp)
      IO::log << "  ***";
    IO::log << "\n";
  }

#ifdef DEBUG

  // calculating energy levels in real space
  Lapack::Vector rl = real_space_energy_levels();

  IO::log << IO::log_offset << "Energy Levels [kcal/mol]:\n";
  IO::log << IO::log_offset 
	  << std::setw(3)  << "#"
	  << std::setw(5)  << "*N"
	  << std::setw(13) << "*M"
	  << std::setw(13) << "*R"
	  << "\n";

  itemp = rl.size() < level_size() ? rl.size() : level_size();
  for(int l = 0; l < itemp; ++l)
    IO::log << IO::log_offset
	    << std::setw(3) << l
	    << std::setw(5) << semiclassical_states_number(energy_level(l) + ground())
	    << std::setw(13) << (energy_level(l) + ground()) / Phys_const::kcal
	    << std::setw(13) << rl[l] / Phys_const::kcal
	    << "\n";
  IO::log << IO::log_offset << "*N  - semiclassical number of states\n"
	  << IO::log_offset << "*M  - momentum space energy levels\n"
	  << IO::log_offset << "*R  - real space energy levels\n";
  
#endif

}// 1D Hindered Rotor

Model::HinderedRotor::~HinderedRotor ()
{
  //std::cout << "Model::HinderedRotor destroyed\n";
}

void  Model::HinderedRotor::set (double ener_max)
{
  const char funame [] = "Model::HinderedRotor::set: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;

  /******************************** setting hamiltonian size *********************************/

  int hsize = 0;

  dtemp = ener_max + ground() - _pot_min;

  if(dtemp > 0.)
    hsize = 2 * (int)std::ceil(std::sqrt(dtemp / rotational_constant()) / (double)symmetry()) + 1;
  
  if(hsize <= _ham_size_min)
    return;

  if(hsize > _ham_size_max) {
    IO::log << IO::log_offset << "WARNING: requested Hamiltonian size = " << hsize 
	    << " exceeds the current limit = " << _ham_size_max << "=> truncating\n";
    hsize = _ham_size_max;
  }
  
  IO::log << IO::log_offset << "hamiltonian size                = " << hsize << "\n";

  // setting quantum levels
  _set_energy_levels(hsize);

  IO::log << IO::log_offset << "ground energy [kcal/mol]        = " 
	  << ground() / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "highest level energy [kcal/mol] = " 
	  << _energy_level.back() / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "number of levels                = "
	  << level_size()  << "\n";


#ifdef DEBUG

  // energy levels
  IO::log << IO::log_offset << "energy levels [kcal/mol]:\n";
  IO::log << IO::log_offset 
	  << std::setw(5) << "#" 
	  << std::setw(13) << "E-E0" 
	  << "\n";
  for(int l = 1; l < level_size(); ++l)
    IO::log << IO::log_offset 
	    << std::setw(5)  << l
	    << std::setw(13) << energy_level(l) / Phys_const::kcal
	    << "\n";  
  
#endif

}

double Model::HinderedRotor::ground () const
{
  return _ground;
}

double Model::HinderedRotor::energy_level (int n) const
{
  return _energy_level[n];
}

int Model::HinderedRotor::level_size () const
{
  return _energy_level.size();
}

double Model::HinderedRotor::weight (double temperature) const
{
  if(_use_quantum_weight)
    return quantum_weight(temperature);

  double cw, sw;
  get_semiclassical_weight(temperature, cw, sw);
  return sw;
}

double Model::HinderedRotor::potential (double angle, int der) const
{
  double dtemp;

  double res = 0.;
  for(std::map<int, double>::const_iterator pit = _pot_four.begin(); pit != _pot_four.end(); ++pit)
    if(pit->first) {
      int n;
      if(pit->first % 2)
	n =  (pit->first + 1) / 2  * symmetry();
      else
	n =  pit->first / 2  * symmetry();

      int fac = 1;
      for(int i = 0; i < der; ++i)
	fac *= n;

      double (*f)(double);
      if((pit->first + der) % 2)
	f = std::sin;
      else
	f = std::cos;

      dtemp = (double)fac * f((double)n * angle) * pit->second;
      
      if((der + 1 - pit->first % 2) / 2 % 2)
	res -= dtemp;
      else
	res += dtemp;
    }
    else if(!der)
      res += pit->second;

  return res;
}

void Model::HinderedRotor::_set_energy_levels (int hsize) 
{
  const char funame [] = "Model::HinderedRotor::_set_energy_levels: ";

  int    itemp;
  double dtemp;


  /************************************ setting Hamiltonian ************************************/

  Lapack::SymmetricMatrix ham(hsize);
  ham = 0.;

  /*
  //kinetic energy contribution
  for(int ml = 1; ml < hsize; ++ml) 
    for(int nl = ml; nl < hsize; ++nl) {// ket state cycle

      int ifac = ((ml + 1) / 2) * ((nl + 1) / 2) * symmetry() * symmetry();
      int m = ml;
      if(ml % 2)// sine
	++m;
      else {// cosine
	--m;
	ifac = -ifac;
      }
      int n = nl;
      if(nl % 2)// sine
	++n;
      else {// cosine
	--n;
	ifac = -ifac;
      }
	
      // fourier expansion cycle
      for(std::map<int, double>::const_iterator  mit = _mobility.begin(); mit != _mobility.end(); ++mit) {
	dtemp = mit->second * (double)ifac;
	if(!rotation_matrix_element(m, n, mit->first, dtemp))
	  ham(m, n) += dtemp;
      }
    }// ket state cycle
  */

  // kinetic energy contribution
  for(int i = 0; i < hsize; ++i) {
    dtemp = double((i + 1) / 2 * symmetry());
    ham(i, i) = rotational_constant() * dtemp * dtemp;
  }

  // potential contribution
  for(int m = 0; m < hsize; ++m)
    for(int n = m; n < hsize; ++n)
      for(std::map<int, double>::const_iterator pit = _pot_four.begin(); pit != _pot_four.end(); ++pit) {
	dtemp = pit->second;
	if(!rotation_matrix_element(m, n, pit->first, dtemp))
	  ham(m, n) += dtemp;
      }

  Lapack::Vector el = ham.eigenvalues();

  // ground state energy
  _ground = el[0];

  // relative excited state energies
  _energy_level.clear();
  _energy_level.reserve(hsize);
  itemp = hsize / 2 * symmetry();
  dtemp = rotational_constant() * double(itemp * itemp) + _pot_min;
  for(int i = 0; i < hsize; ++i) {
    if(el[i] > dtemp)
      break;
    _energy_level.push_back(el[i] - _ground);
  }
}

int Model::HinderedRotor::get_semiclassical_weight (double temperature, double& cw, double& sw) const
{
  static const double eps = 0.01;
  static const double amin = 0.1 - M_PI;

  double dtemp;

  int res = 0;
  if(_freq_min / temperature / 2. < amin)
    res = 1;

  cw = sw = 0.;
  double fac;
  for(int i = 0; i < _pot_grid.size(); ++i) {
    dtemp = _freq_grid[i] / temperature / 2.;
    if(dtemp > eps)
      fac = dtemp / std::sinh(dtemp);
    else if(dtemp < amin)
      fac = -1.;
    else if(dtemp < -eps)
      fac = dtemp / std::sin(dtemp);
    else
      fac = 1.;

    dtemp = std::exp(-_pot_grid[i] / temperature);
    if(fac > 0.)
      sw += dtemp * fac;
    cw += dtemp;
  }

  dtemp = _grid_step * std::sqrt(temperature / rotational_constant() / 4. / M_PI) 
    * std::exp(_ground / temperature);
    //* std::exp(_harm_ground / temperature);
  cw *= dtemp;
  sw *= dtemp;

  return res;
}

double Model::HinderedRotor::quantum_weight (double temperature) const
{
  double dtemp;
  double res = 1.;
  for(int l = 1; l < level_size(); ++l) {
    dtemp = energy_level(l) / temperature;
    if(dtemp > _therm_pow_max)
      break;
    res += std::exp(-dtemp);
  }

  return res;
}

int Model::HinderedRotor::semiclassical_states_number (double ener) const
{
  const char funame [] = "Model::HinderedRotor::semiclassical_states_number: ";

  static const double eps = 1.e-6;

  const double fac = _grid_step / 1.5 / std::sqrt(rotational_constant()) / M_PI;

  double dtemp;

  double y1, y2, z1, z2;

  y1 = ener - *_pot_grid.rbegin();
  bool incomplete = false;
  if(y1 > 0.) {
    z1 = y1 * std::sqrt(y1);
    incomplete = true;
  }

  double inact;
  bool isfirst = true;
  int res = 0;
  double action = 0.;
  for(std::vector<double>::const_iterator it = _pot_grid.begin(); it != _pot_grid.end(); ++it) {
    y2 = ener - *it;
    if(y2 > 0.)
      z2 = y2 * std::sqrt(y2);

    if(y1 > 0. && y2 <= 0.) {
      action += z1  / (y1 - y2);
      if(incomplete && isfirst) 
	inact = action;
      else
	res += int(action * fac + 0.5);
      isfirst = false;
    }
    else if(y2 > 0. && y1 <= 0.) {
      action = z2  / (y2 - y1);
    }
    else if(y1 > 0. && y2 > 0.) {
      dtemp = (y2 - y1) / y1;
      if(dtemp < eps && dtemp > -eps)
	action += 1.5 * z1 / y1;
      else
	action += (z2 - z1) / (y2 - y1);
    }

    y1 = y2;
    z1 = z2;
  }

  if(incomplete && isfirst)
    return 2 * int(action * fac / 2.) + 1;

  if(incomplete)
    res += int((action + inact) * fac + 0.5);

  return res;
}

Lapack::Vector Model::HinderedRotor::real_space_energy_levels () const
{
  double dtemp;

  Lapack::SymmetricMatrix ham(_pot_grid.size());
  dtemp = 2. * rotational_constant() / _grid_step / _grid_step;
  ham = dtemp;
  dtemp /= -2.;
  for(int i = 0; i < ham.size(); ++i) {
    ham(i, i) += _pot_grid[i];
    if(i)
      ham(i - 1, i) = dtemp;
    else
      ham(0, ham.size() - 1) = dtemp;
  }
  return ham.eigenvalues();
}

void Model::HinderedRotor::integrate (Array<double>& states, double ener_step) const
{
  int itemp;
  std::map<int, int> shift;
  for(int g = 0; g < _pot_grid.size(); ++g) {
    itemp = (int)round((_pot_grid[g] - _pot_min) / ener_step);
    shift[itemp]++;
  }

  Array<double> new_states(states.size(), 0.);

#pragma omp parallel for default(shared) schedule(static, 10)

  for(int i = 0; i < states.size(); ++i) {
    for(std::map<int, int>::const_iterator sit = shift.begin(); sit != shift.end(); ++sit)
      if(i >= sit->first)
	new_states[i] += (double)sit->second * states[i - sit->first];
  }

  states = new_states;

  states *= _grid_step;
}

/********************************************************************************************
 ***************************************** UMBRELLA MODE ************************************
 ********************************************************************************************/

Model::Umbrella::Umbrella (IO::KeyBufferStream& from, const std::vector<Atom>& atom) 
  : Rotor(from, atom)
{
  const char funame [] = "Model::Umbrella::Umbrella: ";

  IO::Marker funame_marker(funame);

  KeyGroup UmbrellaModeModel;

  Key    group_key("Group"              );
  Key      ref_key("ReferenceAtom"      );
  Key     plan_key("Plane"              );
  Key      dir_key("Direction"          );
  Key kcal_pot_key("Potential[kcal/mol]");
  Key   kj_pot_key("Potential[kJ/mol]"  );
  Key incm_pot_key("Potential[1/cm]"    );

  int         itemp;
  double      dtemp;
  std::string stemp;
  D3::Vector  vtemp;

  std::set<int> set_temp;

  std::string token, line, comment;
  //std::string token = IO::last_key, line, comment;
  //IO::last_key.clear();

  std::set<int> group_def;// moving group
  int ref_def = -1; // reference atom
  std::vector<int> plan_def;// umbrella plane
  Array<double> dir; // direction of motion
  double xrange = -1.; // coordinate variation range

  while(from >> token) {
    // end input
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // moving group definition
    else if(group_key == token) {
      if(group_def.size()) {
	std::cerr << funame << token <<  ": already defined\n";
	throw Error::Init();
      }
      IO::LineInput group_input(from);
      while(group_input >> itemp) {
	if(itemp < 1) {// Fortran indexing
	  std::cerr << funame << token << ": index out of range\n";
	  throw Error::Range();
	}
	if(!group_def.insert(itemp - 1).second) { 
	  std::cerr << funame << token << ": " << itemp 
		    << "-th atom already has been used in umbrella group defininition\n";
	  throw Error::Init();
	}
      }
    }
    // umbrella plane definition
    else if(plan_key == token) {
      if(dir.size() || plan_def.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      IO::LineInput plan_input(from);
      set_temp.clear();
      for(int i = 0; i < 3; ++i) {
	if(!(plan_input >> itemp))  {
	  std::cerr << funame << token << ": corrupted";
	  throw Error::Input();
	}
	if(itemp < 1) {// Fortran indexing
	  std::cerr << funame << token << ": index out of range\n";
	  throw Error::Range();
	}
	if(!set_temp.insert(itemp).second) { 
	  std::cerr << funame << token << ": " << itemp 
		    << "-th atom already has been used in the plane defininition\n";
	  throw Error::Init();
	}
	plan_def.push_back(itemp - 1);
      }
    }
    // umbrella motion direction
    else if(dir_key == token) {
      if(dir.size() || plan_def.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      dir.resize(3);
      IO::LineInput dir_input(from);
      for(int i = 0; i < 3; ++i) {
	if(!(dir_input >> dir[i]))  {
	  std::cerr << funame << token << ": corrupted";
	  throw Error::Input();
	}
      }
      normalize(dir, 3);
    }
    // reference atom
    else if(ref_key == token) {
      if(ref_def >= 0) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      if(!(from >> ref_def)) {
	std::cerr << funame << token << ": corrupted";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(ref_def-- < 1) {// Fortran indexing
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // potential
    else if(kcal_pot_key == token || kj_pot_key == token || incm_pot_key == token) {
      // sampling size
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": sampling size is corrupted";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(itemp < 3) {
	std::cerr << funame << token << ": not enough sampling points\n";
	throw Error::Range();
      }

      // sampling data
      std::map<double, double> pot_data;
      for(int i = 0; i < itemp; ++i) {
	IO::LineInput data_input(from);
	double xval;

	if(!(data_input >> xval)) {
	  std::cerr << funame << token << ": cannot read " << i << "-th coordinate value\n";
	  throw Error::Input();
	}

	if(!(data_input >> dtemp)) {
	  std::cerr << funame << token << ": cannot read " << i << "-th potential value\n";
	  throw Error::Input();
	}

	if(kcal_pot_key == token)
	  dtemp *= Phys_const::kcal;
	if(incm_pot_key == token)
	  dtemp *= Phys_const::incm;
	if(kj_pot_key == token)
	  dtemp *= Phys_const::kjoul;

	pot_data[xval] = dtemp;
      }

      // coordinate variation range
      xrange = pot_data.rbegin()->first - pot_data.begin()->first;

      // convert potential samplings into power expansion coefficients
      _pot_coef.resize(pot_data.size());
      Lapack::Vector xval(pot_data.size()); 
      std::map<double, double>::const_iterator pit;
      for(pit = pot_data.begin(), itemp = 0; pit != pot_data.end(); ++pit, ++itemp) {	
	xval[itemp] = (pit->first  - pot_data.begin()->first) / xrange;
	_pot_coef[itemp] = pit->second;
      }

      // conversion matrix
      Lapack::Matrix xmat(xval.size(), xval.size());
      for(int i = 0; i < xval.size(); ++i) {
	dtemp = 1.;
	for(int j = 0; j < xval.size(); ++j, dtemp *= xval[i])
	  xmat(i, j) = dtemp;
      }

      //_pot_coef = xmat.invert() * _pot_coef; // does not work
      //_pot_coef = Lapack::LU(xmat).invert(_pot_coef);
      _pot_coef = xmat.invert(_pot_coef);

      // potential at the sampling points
      IO::log << IO::log_offset << "potential[kcal/mol] at the sampling points:\n"
	      << IO::log_offset
	      << std::setw(7) << "X"
	      << std::setw(7) << "V"
	      << "\n" << std::setprecision(3) << std::fixed;
      for(int i = 0; i < xval.size(); ++i)
	IO::log << IO::log_offset 
		<< std::setw(7) << xval[i] * xrange
		<< std::setw(7) << potential(xval[i]) / Phys_const::kcal 
		<< "\n";
      IO::log << std::setprecision(6) << std::resetiosflags(std::ios_base::floatfield);
    }    
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  } //while(from >> token);
  
  /******************************************* Checking *************************************************/

  // stream state
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(!_pot_coef.size()) {
    std::cerr << funame << "potential is not defined\n";
    throw Error::Init();
  }

  if(!atom.size()) {
    std::cerr << funame << "geometry is not defined\n";
    throw Error::Init();
  }

  if(!dir.size() && !plan_def.size()) {
    std::cerr << funame << "motion direction is not defined\n";
    throw Error::Init();
  }

  if(!group_def.size()) {
    std::cerr << funame << "moving group is not defined\n";
    throw Error::Init();
  }

  if(*group_def.rbegin() >= atom.size()) {
    std::cerr << funame << "atom index in moving group definition is out of range\n";
    throw Error::Range();
  }

  if(group_def.find(ref_def) != group_def.end()) {
    std::cerr << funame << "reference atom should not belong to moving group\n";
    throw Error::Logic();
  }

  if(ref_def >= atom.size()) {
    std::cerr << funame << "refernce atom index  is out of range\n";
    throw Error::Range();
  }

  if(plan_def.size())
    for(int i = 0; i < 3; ++i)
      if(plan_def[i] >= atom.size()) {
	std::cerr << funame << i << "-th atom index in the umbrella plane definition is out of range\n";
	throw Error::Range();
      }

  /******************************************* Effective mass *************************************************/

  // umbrella motion direction
  if(plan_def.size()) {
    dir.resize(3);
    dir = D3::vprod(atom[plan_def[1]] - atom[plan_def[0]], atom[plan_def[2]] - atom[plan_def[0]]);
    normalize(dir, 3);
  }

  // umbrella motion
  std::map<int, D3::Vector> motion;
  if(ref_def < 0)
    for(std::set<int>::const_iterator  cit = group_def.begin(); cit != group_def.end(); ++cit)
      motion[*cit] = dir;
  else
    for(std::set<int>::const_iterator  cit = group_def.begin(); cit != group_def.end(); ++cit) {
      vtemp = dir;
      dtemp = vdistance(atom[*cit], atom[ref_def], 3);
      vtemp *= dtemp;
      motion[*cit] = vtemp;
    }

  // adjust normal mode so that angular, and linear momenta disapear
  std::vector<D3::Vector> normal_mode = adjusted_normal_mode(atom, motion);  

  // effective mass 
  _mass = 0.;
  for(int a = 0; a < atom.size(); ++a)
    _mass += atom[a].mass() * vdot(normal_mode[a]);

  /******************************************* POTENTIAL *************************************************/

  if(ref_def < 0) {
    IO::log << IO::log_offset << "coordinate in the potential definition assumed to be in angstrom\n";
    xrange *= Phys_const::angstrom;
  }
  else {
    IO::log << IO::log_offset << "coordinate in the potential definition assumed to be in angular degrees\n";
    xrange *= M_PI / 180.;
  }

  // mass renormalization
  _mass *= xrange * xrange;

  // minimum potential energy and potential discretization
  _astep = 1. / double(_grid_size - 1);
  _pot_grid.resize(_grid_size);
  _freq_grid.resize(_grid_size);

  int imin;
  double a = 0;
  for(int i = 0; i < _grid_size; ++i, a += _astep) {
    // potential
    dtemp = potential(a);
    _pot_grid[i] = dtemp; 
    if(!i || dtemp < _pot_min) {
      _pot_min = dtemp;
      imin = i;
    }
    // frequency
    dtemp = potential(a, 2);
    if(dtemp < 0.)
      dtemp = -std::sqrt(- dtemp / _mass);
    else
      dtemp = std::sqrt(dtemp / _mass);
    _freq_grid[i] = dtemp;
  }

  // numerical vibrational frequency
  if(imin != 0 && imin != _grid_size - 1) {  
    dtemp = _pot_grid[imin + 1] + _pot_grid[imin - 1] - 2. * _pot_min;
    if(dtemp >= 0.)
      dtemp = std::sqrt(dtemp / _mass) / _astep;
    else
      dtemp = -std::sqrt(- dtemp / _mass) / _astep;
    IO::log << IO::log_offset << "numerical frequency at minimum [1/cm] = " << dtemp  / Phys_const::incm << "\n";    
  }

  // analytical vibrational frequency at the minimum
  dtemp = potential((double)imin * _astep, 2);
  if(dtemp >= 0.)
    dtemp = std::sqrt(dtemp / _mass);
  else
    dtemp = -std::sqrt(- dtemp / _mass);
  IO::log << IO::log_offset << "analytic frequency at minimum [1/cm]  = " << dtemp / Phys_const::incm << "\n";

  // energy levels in momentum space
  _set_energy_levels(_ham_size_min);

  IO::log << IO::log_offset << "minimum energy [kcal/mol]             = " 
	  << _pot_min  / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "number of quantum levels              = " 
	  << level_size()  << "\n";
  IO::log << IO::log_offset << "ground energy [kcal/mol]              = "
	  << ground() / Phys_const::kcal << "\n";
  IO::log << IO::log_offset << "highest level energy [kcal/mol]       = "
	  << (_energy_level.back() + ground()) / Phys_const::kcal << "\n";

  // statistical weight output
  IO::log << IO::log_offset << "statistical weight (*** - deep tunneling regime):\n";
  IO::log << IO::log_offset 
	  << std::setw(5) << "T, K" 
	  << std::setw(15) << "Quantum"
	  << std::setw(15) << "Classical"
	  << std::setw(15) << "Semiclassical"
	  << "  ***\n";
  for(int t = 100; t <= 1000 ; t+= 100) {
    double tval = (double)t * Phys_const::kelv;
    double cw, sw;
    itemp = get_semiclassical_weight(tval, cw, sw);
    IO::log << IO::log_offset 
	    << std::setw(5) << t
	    << std::setw(15) << quantum_weight(tval)
	    << std::setw(15) << cw
	    << std::setw(15) << sw;
    if(itemp)
      IO::log << "  ***";
    IO::log << "\n";
  }

  // energy levels output
  itemp = level_size() < 10 ? level_size() : 10;
  IO::log << IO::log_offset << itemp << " lowest excited states[kcal/mol]:" << std::setprecision(3);
  for(int l = 1; l < itemp; ++l)
    IO::log  << " " << (energy_level(l) + ground()) / Phys_const::kcal;
  IO::log << std::setprecision(6) << "\n";

}// Umbrella

Model::Umbrella::~Umbrella ()
{
  //std::cout << "Model::Umbrella destroyed\n";
}

void Model::Umbrella::set (double ener_max)
{
  const char funame [] = "Model::Umbrella::set: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;

  int hsize = 0;
  
  dtemp = ener_max + ground() - _pot_min;

  if(dtemp >= 0.)
    hsize = (int)std::ceil(std::sqrt(2. * dtemp * _mass)/ M_PI);
  
  if(hsize <= _ham_size_min)
    return;

  if(hsize > _ham_size_max) {
    IO::log << IO::log_offset << "WARNING: requested Hamiltonian size = " << hsize 
	    << " exceeds the current limit = " << _ham_size_max << "=> truncating\n";
    hsize = _ham_size_max;
  }
  
  IO::log << IO::log_offset << "Hamiltonian Size = " << hsize << "\n";
  
  // setting quantum levels
  _set_energy_levels(hsize);

}


void Model::Umbrella::_set_energy_levels (int hsize) 
{
  const char funame [] = "Model::Umbrella::_set_energy_levels: ";

  int    itemp;
  double dtemp;

  // setting  hamiltonian
  Lapack::SymmetricMatrix ham(hsize);

  for(int m = 0; m < hsize; ++m)
    for(int n = m; n < hsize; ++n) {
      dtemp = 0.;
      for(int p = 1; p < _pot_coef.size(); ++p)
	dtemp += _pot_coef[p] * (_integral(p, n - m) - _integral(p, m + n + 2));
      ham(m, n) = dtemp;
    }

  for(int m = 0; m < hsize; ++m) {
    dtemp = M_PI * double(m + 1);
    ham(m, m) += _pot_coef[0] + dtemp * dtemp / 2. / _mass;
  }

  Lapack::Vector el = ham.eigenvalues();

  // updating energy levels
  _ground  = el.front();
  _energy_level.clear();
  _energy_level.reserve(hsize);
  dtemp = (double)hsize * M_PI;
  dtemp = dtemp * dtemp / 2. / _mass + _pot_min;
  for(int i = 0; i < hsize; ++i) {
    if(el[i] > dtemp)
      break;
    _energy_level.push_back(el[i] - _ground);
  }
}

double Model::Umbrella::potential(double x, int der) const 
{
  const char funame [] = "Model::Umbrella::potential: ";

  if(der < 0) {
    std::cerr << funame << "negative derivative power requested\n";
    throw Error::Range();
  }

  int    itemp;
  double dtemp;

  double xfac = 1.;
  double res = 0.;
  for(int i = der; i < _pot_coef.size(); ++i, xfac *= x) {
    dtemp = _pot_coef[i] * xfac;
    if(der) {
      itemp = 1;
      for(int j = i - der + 1; j <= i; ++j)
	itemp *= j;
      dtemp *= (double)itemp;
    }
    res += dtemp;
  }
    
  return res;
}

double Model::Umbrella::_integral(int p, int n)
{
  if(!n)
    return 1. / (1. + (double)p);

  if(!p)
    return 0.;

  n = n < 0 ? -n : n;

  int v = n % 2 ? -1 : 1;
  double dtemp = M_PI * (double)n;
  dtemp *= dtemp;

  if(p == 1)
    return double(v - 1) / dtemp;

  return ((double)v - double(p - 1) * _integral(p - 2, n)) * double(p) / dtemp;
}

double Model::Umbrella::ground () const
{
  return _ground;
}

double Model::Umbrella::energy_level (int i) const
{
  return _energy_level[i];
}

int Model::Umbrella::level_size () const
{
  return _energy_level.size();
}

double Model::Umbrella::quantum_weight (double temperature) const
{
  double dtemp;

  double res = 1.;
  for(int l = 1; l < _energy_level.size(); ++l) {
    dtemp = _energy_level[l] / temperature;
    if(dtemp > _therm_pow_max)
      break;
    res += std::exp(-dtemp);
  }

  return res;
}

double Model::Umbrella::weight (double temperature) const
{
  return quantum_weight(temperature);

  //  double cw, sw;
  //get_semiclassical_weight(temperature, cw, sw);
  //return sw;
}

int Model::Umbrella::get_semiclassical_weight (double temperature, double& cw, double& sw) const
{
  static const double eps = 0.01;
  static const double amin = 0.1 - M_PI;

  double dtemp;

  int res = 0;

  cw = sw = 0.;
  double fac;
  for(int i = 0; i < _pot_grid.size(); ++i) {
    dtemp = _freq_grid[i] / temperature / 2.;
    if(dtemp > eps)
      fac = dtemp / std::sinh(dtemp);
    else if(dtemp < amin) {
      fac = -1.;
      res = 1;
    }
    else if(dtemp < -eps)
      fac = dtemp / std::sin(dtemp);
    else
      fac = 1.;

    dtemp = std::exp(-_pot_grid[i] / temperature);
    if(fac > 0.)
      sw += dtemp * fac;
    cw += dtemp;
  }

  dtemp = _astep * std::sqrt(temperature * _mass / 2. / M_PI) 
    * std::exp(_ground / temperature);
  cw *= dtemp;
  sw *= dtemp;

  return res;
}

/********************************************************************************************
 *****************************************  CORE ********************************************
 ********************************************************************************************/

Model::Core::~Core ()
{
  //std::cout << "Model::Core destroyed\n";
}

/********************************************************************************************
 *************************** PHASE SPACE THEORY NUMBER OF STATES ****************************
 ********************************************************************************************/

Model::PhaseSpaceTheory::PhaseSpaceTheory (IO::KeyBufferStream& from) 
  : Core(NUMBER), _states_factor(1.)
{
  const char funame [] = "Model::PhaseSpaceTheory::PhaseSpaceTheory: ";

  IO::Marker funame_marker(funame);

  double pot_exp = -1.; // power exponent n in potential expression V(R) = V_0 / R^n;
  double pot_fac = -1.; // potential prefactor V_0

  int    itemp;
  double dtemp;

  KeyGroup PhaseSpaceTheoryModel;

  Key ang_geom_key("FragmentGeometry[angstrom]");
  Key bor_geom_key("FragmentGeometry[au]"      );
  Key     mass_key("FragmentMasses[amu]"       );
  Key  wn_rcon_key("RotationalConstants[1/cm]" );
  Key ghz_rcon_key("RotationalConstants[GHz]"  );
  Key     symm_key("SymmetryFactor"            );
  Key     pfac_key("PotentialPrefactor[au]"    );
  Key     pexp_key("PotentialPowerExponent"    );
  
  int gcount = 0; // number of read geometries
  std::vector<double> frag_mass; // fragments masses
  std::vector<double> frag_rcon; // fragments rotational constants
  

  std::string token, comment;
  while(from >> token) {
    // end input
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // fragment geometry
    else if(ang_geom_key == token || bor_geom_key == token) {
      if(!gcount && (frag_mass.size() || frag_rcon.size())) {
	std::cerr << funame << token << ": no geometry if fragments masses or rotational constants are used\n";
	throw Error::Init();
      }
      gcount++;

      std::vector<Atom> atom;
      if(ang_geom_key == token)
	read_geometry(from, atom, ANGSTROM);
      if(bor_geom_key == token)
	read_geometry(from, atom, BOHR);

      dtemp = 0.;
      for(std::vector<Atom>::iterator at = atom.begin(); at != atom.end(); ++at)
	dtemp += at->mass();
      
      frag_mass.push_back(dtemp);

      // rotational constants
      if(atom.size() != 1) {
	IO::log << IO::log_offset << "molecular geometry will be used to estimate rotational constants\n";
	Lapack::Vector imom = inertia_moment_matrix(atom).eigenvalues();
	if(imom[0] / imom[1] < 1.e-5) {// linear configuration
	  dtemp = 0.5 / imom[1];
	  IO::log << IO::log_offset << "rotational configuration: linear\n"
		  << IO::log_offset << "rotational constant[1/cm]: "
		  << std::setw(15) << dtemp / Phys_const::incm
		  << "\n";
	  frag_rcon.push_back(dtemp);
	  frag_rcon.push_back(dtemp);
	}
	else {// nonlinear configuration
	  IO::log << IO::log_offset << "rotational configuration: nonlinear\n"
		  << IO::log_offset << "rotational constants[1/cm]: ";
	  for(int i = 0; i < 3; ++i) {
	    dtemp = 0.5 / imom[i];
	    frag_rcon.push_back(dtemp);
	    IO::log << std::setw(15) << dtemp / Phys_const::incm;
	  }
	  IO::log << "\n";
	}
      }
    }
    // fragment masses
    else if(mass_key == token) {
      if(gcount || frag_mass.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      IO::LineInput mass_input(from);
      for(int i = 0; i < 2; ++i) {
	if(!(mass_input >> dtemp)) {
	  std::cerr << funame << token << ": corrupted\n";
	  throw Error::Input();
	}
      
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}
	frag_mass.push_back(dtemp * Phys_const::amu);
      }
    }
    // rotational constants
    else if(wn_rcon_key == token || ghz_rcon_key == token) {
      if(gcount || frag_rcon.size()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      IO::LineInput rcon_input(from);
      while(rcon_input >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}

	if(wn_rcon_key == token)
	  dtemp *= Phys_const::incm;
	else
	  dtemp *= Phys_const::herz * 2.e+9 * M_PI;

	frag_rcon.push_back(dtemp);
      }
    }
    // symmetry factor (number of symmetry operations)
    else if(symm_key == token) {
      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      _states_factor = 1./dtemp;
    }
    // potential prefactor V_0
    else if(pfac_key == token) {
      if(!(from >> pot_fac)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(pot_fac <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
    }
    // potential distance dependence power exponent
    else if(pexp_key == token) {
      if(!(from >> pot_exp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(pot_exp <= 1.) {
	std::cerr << funame << token << ": should be more than 1\n";
	throw Error::Range();
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
 
  if(!from) {
    std::cerr << funame << "input corrupted\n";
    throw Error::Input();
  }
  
  if(frag_mass.size() != 2) {
    std::cerr << funame << "wrong number of fragments\n";
    throw Error::Init();
  }

  // numerical factor for the number of states
  switch(frag_rcon.size()) {
  case 2:// 2 + 0 
    break;
  case 3:// 3 + 0
    _states_factor *= 16. / 15.;
    break;
  case 4:// 2 + 2
    _states_factor *= 1. / 3.;
    break;
  case 5:// 2 + 3
    _states_factor *= 32. / 105.;
    break;
  case 6:// 3 + 3
    _states_factor *= M_PI / 12.;
    break;
  default:
    std::cerr << funame << "wrong number of rotational constants: " << frag_rcon.size() << "\n";
    throw Error::Init();
  }
  // effective mass factor
  dtemp = 0.;
  for(int i = 0; i < 2; ++i)
    dtemp += 1. / frag_mass[i];
  _states_factor /= dtemp;

  // rotational constants factor;
  dtemp = 1.;
  for(int i  = 0; i < frag_rcon.size(); ++i)
    dtemp *= frag_rcon[i];
  _states_factor /= std::sqrt(dtemp);

  // power parameter
  itemp = frag_rcon.size() + 2;
  dtemp = (double)(itemp) * pot_exp / 4. - 1.;
  _states_factor *= std::pow(dtemp * pot_fac, 2. / pot_exp) * std::pow(1. + 1. / dtemp, (double)itemp / 2.); 
  _power = (double)itemp / 2. - 2. / pot_exp;

  _weight_factor = _states_factor * tgamma(_power + 1.);

}// Phase Space Theory

Model::PhaseSpaceTheory::~PhaseSpaceTheory ()
{
  //std::cout << "Model::PhaseSpaceTheory destroyed\n";
}

double Model::PhaseSpaceTheory::ground       () const
{
  return 0.;
}

double Model::PhaseSpaceTheory::weight (double temperature) const
{
  return  _weight_factor * std::pow(temperature, _power);  
}

double Model::PhaseSpaceTheory::states (double ener) const
{
  return _states_factor * std::pow(ener, _power);
}

/********************************************************************************************
 *********************************** RIGID ROTOR CORE MODEL *********************************
 ********************************************************************************************/

Model::RigidRotor::RigidRotor(IO::KeyBufferStream& from, const std::vector<Atom>& atom,  int m)
   : Core(m), _factor(1.), _rdim(0), _rofactor(0.), _ground(0.), _emax(-1.) 
{
  const char funame [] = "Model::RigidRotor::RigidRotor: ";

  IO::Marker funame_marker(funame);

  int         itemp;
  double      dtemp;
  std::string stemp;
  std::vector<double> vtemp;

  KeyGroup RigidRotorModel;

  Key      rcon_key("RotationalConstants[1/cm]"       );
  Key      symm_key("SymmetryFactor"                  );
  Key      freq_key("Frequencies[1/cm]"               );
  Key    fscale_key("FrequencyScalingFactor"          );
  Key     degen_key("FrequencyDegeneracies"           );
  Key      harm_key("AreFrequenciesHarmonic"          );
  Key    anharm_key("Anharmonicities[1/cm]"           );
  Key     rovib_key("RovibrationalCouplings[1/cm]"    );
  Key      dist_key("RotationalDistortion[1/cm]"      );
  Key      emax_key("InterpolationEnergyMax[kcal/mol]");
  Key     estep_key("InterpolationEnergyStep[1/cm]"   );
  Key     extra_key("ExtrapolationStep"               );
  Key      zero_key("ZeroPointEnergy[1/cm]"           );
  
  bool issym  = false;
  bool isharm = false;
  bool iszero = false;
  
  double extra_step = 0.1;
  double ener_quant = Phys_const::incm;
  double fscale     = -1.;

  Lapack::SymmetricMatrix distort;
  std::vector<double> rcon;
  
  std::string token, comment;

  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // rotational dimension and constants
    else if(rcon_key == token) {
      if(_rdim > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      IO::LineInput data_input(from);
      while(data_input >> dtemp) {
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": should be positive\n";
	  throw Error::Range();
	}
	rcon.push_back(dtemp * Phys_const::incm);
      }

      switch(rcon.size()) {
      case 1:
	_rdim = 2;
	_factor /= rcon[0];
	break;
      case 2:
	if(rcon[0] != rcon[1])
	  IO::log << IO::log_offset << funame 
		  << "WARNING: inertia moments for linear molecule read from input file differ\n";
	_rdim = 2;
	_factor /= rcon[0];
	break;
      case 3:
	_rdim = 3;
	dtemp = 1.;
	for(int i = 0; i < 3; ++i)
	  dtemp *= rcon[i];
	_factor /= std::sqrt(dtemp);
	break;
      default:
	std::cerr << funame << token << ": wrong size\n";
	throw Error::Range();
      }
    }
    //  frequency quantum
    //
    else if(estep_key == token) {
      if(!(from >> ener_quant)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(ener_quant <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      ener_quant *= Phys_const::incm;
    }
    // interpolation maximal energy
    //
    else if(emax_key == token) {
      if(!(from >> _emax)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_emax <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      _emax *= Phys_const::kcal;
    }
    // extrapolation step
    //
    else if(extra_key == token) {
      if(!(from >> extra_step)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(extra_step <= 0. || extra_step >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // zero-point energy
    //
    else if(zero_key == token) {
      //
      if(!(from >> _ground)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_ground < 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      
      _ground *= Phys_const::incm;
      iszero = true;
    }
    // frequencies
    //
    else if(freq_key == token) {
      //
      if(_frequency.size()) {
	std::cerr << funame << token << ": have been initialized already\n";
	throw Error::Init();
      }
      
      // number of frequencies
      //
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": cannot read number of frequencies\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(itemp < 1) {
	std::cerr << funame << token << ": number of frequencies should be positive\n";
	throw Error::Range();
      }
      _frequency.resize(itemp);
      
      // read frequencies
      //
      for(int i = 0; i < _frequency.size(); ++i) {
	if(!(from >> dtemp)) {
	  std::cerr << funame << token << ": cannot read " << i << "-th frequency\n";
	  throw Error::Input();
	}
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": " << i << "-th frequency: should be positive\n";
	  throw Error::Range();
	}
	_frequency[i] = dtemp * Phys_const::incm;
      }
      
      if(!from) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // frequency scaling factor
    //
    else if(fscale_key == token) {
      //
      if(fscale > 0.) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": already initialized";
      }
	  
      if(!(from >> fscale)) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": corrupted";
      }
      
      if(fscale <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // are frequencies harmonic
    //
    else if(harm_key == token) {
      //
      std::getline(from, comment);
      
      isharm = true;
    }
    // frequency degeneracies
    //
    else if(degen_key == token) {
      //
      if(_fdegen.size()) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": already initialized";
      }
      
      // number of degeneracies
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": cannot read number of frequency degeneracies";
      }
      std::getline(from, comment);

      if(itemp < 1) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": number of frequency degeneracies should be positive";
      }
      _fdegen.resize(itemp);
      
      // read frequency degeneracies
      //
      for(int i = 0; i < _fdegen.size(); ++i) {
	//
	if(!(from >> itemp)) {
	  //
	  ErrOut err_out;
	  
	  err_out << funame << token << ": cannot read " << i << "-th degeneracy";
	}
	
	if(itemp <= 0) {
	  //
	  ErrOut err_out;
	  
	  err_out << funame << token << ": " << i << "-th degeneracy: should be positive";
	}
	_fdegen[i] = itemp;
      }

      if(!from) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": corrupted";
      }
      
      std::getline(from, comment);
    }
    // anharmonicities
    //
    else if(anharm_key == token) {
      //
      if(!_frequency.size()) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": frequencies should be initialized first";
      }
	  
      std::getline(from, comment);
      
      _anharm.resize(_frequency.size());
      //
      _anharm = 0.;
      
      for(int i = 0; i < _anharm.size(); ++i)
	//
	for(int j = 0; j <= i; ++j)
	  //
	  if(!(from >> _anharm(i, j))) {
	    //
	    ErrOut err_out;
	    
	    err_out << funame << token << ": cannot read anharmonicities";
	  }
      
      _anharm *= Phys_const::incm;
      
    }
    // rovibrational couplings
    //
    else if(rovib_key == token) {
      //
      if(!_frequency.size()) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": frequencies should be initialized first";
      }
      std::getline(from, comment);
      
      _rvc.resize(_frequency.size());

      for(int f = 0; f < _rvc.size(); ++f) {
	//
	IO::LineInput lin(from);
	
	while(lin >> dtemp)
	  //
	  _rvc[f].push_back(dtemp * Phys_const::incm);
      }
    }
    // centrifugal distortion
    //
    else if(dist_key == token) {
      //
      if(distort.isinit()) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": already initialized";
      }
	  
      distort.resize(3);
      
      distort = 0.;
      
      std::getline(from, stemp);
      
      while(1) {
	//
	IO::LineInput lin(from);
	
	lin >> stemp;
	
	if(stemp == IO::end_key())
	  //
	  break;

	if(!(lin >> dtemp)) {
	  //
	  ErrOut err_out;
	  
	  err_out << funame << token << ": corrupted";
	}
	
	dtemp *= Phys_const::incm;

	if(stemp == "aaaa") {
	  //
	  distort(0,0) += 0.75 * dtemp;
	}
	else if(stemp == "bbbb") {
	  //
	  distort(1,1) += 0.75 * dtemp;
	}
	else if(stemp == "cccc") {
	  //
	  distort(2,2) += 0.75 * dtemp;
	}
	else if(stemp == "bbaa") {
	  //
	  distort(0,1) += 0.5 * dtemp;
	}
	else if(stemp == "ccbb") {
	  //
	  distort(1,2) += 0.5 * dtemp;
	}
	else if(stemp == "ccaa") {
	  //
	  distort(0,2) += 0.5 * dtemp;
	}
	else if(stemp == "baba") {
	  //
	  distort(0,1) += dtemp;
	}
	else if(stemp == "cbcb") {
	  //
	  distort(1,2) += dtemp;
	}
	else if(stemp == "caca") {
	  //
	  distort(0,2) += dtemp;
	}
	else {
	  //
	  ErrOut err_out;
	  
	  err_out << funame << token << ": unknown index: " << stemp
		  << ": available indices: aaaa, bbbb, cccc, bbaa, ccaa, ccbb, baba, cbcb, caca";
	}
      }
    }
    // symmetry factor (number of symmetry operations)
    //
    else if(symm_key == token) {
      //
      if(issym) {
	//
	ErrOut err_out;
	
	err_out << funame << "symmetry number has been initialized already";
      }
      
      issym = true;
      
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": corrupted";
      }
      
      std::getline(from, comment);

      if(dtemp <= 0.) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": should be positive";
      }
      
      _factor /= dtemp;
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      ErrOut err_out;
      
      err_out << funame << "unknown keyword " << token << "\n" << Key::show_all() << "\n";
    }
  }
 
  if(!from) {
    //
    ErrOut err_out;
    
    err_out << funame << "input corrupted";
  }

  // rotational constants from geometry
  if(!_rdim) {
    IO::log << IO::log_offset << "molecular geometry will be used to estimate rotational constants\n";
    if(!atom.size()) {
      std::cerr << funame << "geometry is not available\n";
      throw Error::Init();
    }
    Lapack::Vector imom = inertia_moment_matrix(atom).eigenvalues();
    if(imom[0] / imom[1] < 1.e-5) {// linear configuration
      _rdim = 2;
      dtemp = 0.5 / imom[1];
      rcon.push_back(dtemp);
      rcon.push_back(dtemp);
      
      _factor *= 2. * imom[1];
      
      IO::log << IO::log_offset << "rotational configuration: linear\n";
      IO::log << IO::log_offset << "rotational constant[1/cm]: "
	      << std::setw(15) << 0.5 / imom[1] / Phys_const::incm
	      << "\n";
    }
    else {// nonlinear configuration
      _rdim = 3;
      for(int i = 0; i < 3; ++i) {
	rcon.push_back(0.5 / imom[i]);
	_factor *= std::sqrt(2. * imom[i]);
      }
      
      IO::log << IO::log_offset << "rotational configuration: nonlinear\n";
      IO::log << IO::log_offset << "rotational constants[1/cm]: ";
      for(int i = 0; i < 3; ++i)
	IO::log << std::setw(15) << 0.5 / imom[i] / Phys_const::incm;
      IO::log << "\n";

      if(distort.isinit()) {
	for(int i = 0; i < 3; ++i)
	  for(int j = i; j < 3; ++j)
	    _rofactor += distort(i, j) * imom[i] * imom[j];
	
	IO::log << IO::log_offset << "centrifugal distortion factor (rho factor, 1/K): " << _rofactor * Phys_const::kelv << "\n";
      }
    }
  } // rotational constants from geometry
  // centrifugal distortion factor
  else if(_rdim == 3 && distort.isinit()) {
    for(int i = 0; i < 3; ++i)
      for(int j = i; j < 3; ++j)
	_rofactor += distort(i, j) / rcon[i] / rcon[j];
    
    _rofactor /= 4.;
    
    IO::log << IO::log_offset << "centrifugal distortion factor (rho factor, 1/K): " << _rofactor * Phys_const::kelv << "\n";
  }

  // anharmonicities, rovibrational couplings, etc.
  if(_frequency.size()) {
    if(!iszero) {
      std::cerr << funame << "zero-point energy is not defined\n";
      throw Error::Init();
    }
    
    if(fscale > 0.)
      for(int f = 0; f < _frequency.size(); ++f)
	_frequency[f] *= fscale;
    
    if(_fdegen.size() > _frequency.size()) {
      std::cerr << funame << "number of frequency degeneracies exceeds number of frequencies\n";
      throw Error::Init();
    }

    for(int i = _fdegen.size(); i < _frequency.size(); ++i)
      _fdegen.push_back(1);

    if(isharm && _anharm.size()) {
      for(int i = 0; i < _frequency.size(); ++i)
	for(int j = 0; j < _frequency.size(); ++j)
	  if(j != i)
	    _frequency[i] += double(_fdegen[j]) * _anharm(i, j) / 2.;
	  else
	    _frequency[i] += double (_fdegen[i] + 1) * _anharm(i, i);

      IO::log << IO::log_offset << "fundamentals(1/cm):";
      for(int i = 0; i < _frequency.size(); ++i)
	IO::log << "   " << _frequency[i] / Phys_const::incm;
      IO::log << "\n";
    }
      
    if(_emax < 0.)
      _emax = energy_limit();
    itemp = (int)std::ceil(_emax / ener_quant);

    if(mode() != NOSTATES) {
      Array<double> ener_grid(itemp);
      Array<double> stat_grid(itemp);
      ener_grid[0] = 0.;
      stat_grid[0] = 0.;

      dtemp = ener_quant;
      for(int i = 1; i < ener_grid.size(); ++i, dtemp += ener_quant) {
	// energy grid
	ener_grid[i] = dtemp;
	// core density / number of states
	stat_grid[i] = _core_states(dtemp);
      }

      for(int f = 0; f < _frequency.size(); ++f) {
	itemp = (int)round(_frequency[f] / ener_quant);
	if(itemp < 0) {
	  std::cerr << funame << "negative frequency\n";
	  throw Error::Range();
	}
	for(int i = 0; i < _fdegen[f]; ++i)
	  for(int e = itemp + 1; e < ener_grid.size(); ++e)
	    stat_grid[e] += stat_grid[e - itemp];
      }

      _states.init(ener_grid, stat_grid, ener_grid.size());
      dtemp = _states.fun_max() / _states(_states.arg_max() * (1. - extra_step));
      _nmax = std::log(dtemp) / std::log(1. / (1. - extra_step));
    }
    
    // rovibrational coupling
    if(_rvc.size()) {
      // check the dimensionality
      for(int f = 0; f < _frequency.size(); ++f)
	if(_rvc[f].size() != _rdim) {
	  std::cerr << funame << f + 1 << "-th frequency: rovibrational coupling size " << _rvc[f].size()
		    << " inconsistent with rotational dimension " << _rdim << "\n";
	  throw Error::Range();
	}

      // zero-point rotational constant correction
      for(int i = 0; i < _rdim; ++i) {
	dtemp = 0.;
	for(int f = 0; f < _frequency.size(); ++f)
	  dtemp += _rvc[f][i] * _fdegen[f] / rcon[i];
	
	dtemp = 1. - 0.5 * dtemp;
	
	if(dtemp <= 0.) {
	  std::cerr << funame << "negative zero-point correction to " << i + 1 << "-th rotational constant\n";
	  throw Error::Range();
	}
	_factor /= std::sqrt(dtemp);
	
	rcon[i] *= dtemp;

	for(int f = 0; f < _frequency.size(); ++f)
	  _rvc[f][i] /= rcon[i];
	  
      }
    } // rovibrational coupling
  }
}// Rigid Rotor

double Model::RigidRotor::_rovib (int rank ...) const
{
  const char funame [] = " Model::RigidRotor::_rovib: ";

  double dtemp;
  
  if(!_rvc.size()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(rank < 1) {
    std::cerr << funame << "rank out of range\n";
    throw Error::Range();
  }

  std::vector<int> ff(rank);

  va_list ap;
  va_start(ap, rank);

  for(int i = 0; i < rank; ++i) {
    ff[i] = va_arg(ap, int);
    
    if(ff[i] < 0 || ff[i] >= _frequency.size()) {
      va_end(ap);
      std::cerr << "index out of range\n";
      throw Error::Range();
    }
  }
  va_end(ap);
  
  double res = 0.;
  for(int i = 0; i < _rdim; ++i) {
    dtemp = 1.;
    for(int r = 0; r < rank; ++r)
      dtemp *= _rvc[ff[r]][i];  
    res += dtemp;
  }

  return res / 2. / double(rank);  
}

Model::RigidRotor::~RigidRotor ()
{
  //std::cout << "Model::RigidRotor destroyed\n";
}

double Model::RigidRotor::ground () const
{
  return _ground;
}

double Model::RigidRotor::weight (double temperature) const
{
  const char funame [] = "Model::RigidRotor::weight: ";
  
  double dtemp;
  int    itemp;
  
  //if(_potex.size())
  //_graph_perturbation_theory_correction(temperature);

  // core contribution
  double res = _core_weight(temperature);

  if(!_frequency.size())
    return res;
  
  // vibrational contribution
  std::vector<double> rfactor(_frequency.size()), sfactor(_frequency.size());
  
  for(int  f = 0; f < _frequency.size(); ++f) {
    rfactor[f] = std::exp(- _frequency[f] / temperature);
    sfactor[f] = 1. / (1. - rfactor[f]);
    for(int i = 0; i < _fdegen[f]; ++i)
      res *= sfactor[f];
  }

  double acl;
  
  // vibrational anharmonic contributions
  //
  if(_anharm.isinit()) {
    //
    // first order correction
    //
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += _anharm(i, i) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i] * sfactor[i], 2.);

      for(int j = 0; j < i; ++j)
	acl += _anharm(i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * sfactor[i] * rfactor[j] * sfactor[j];
    }
    acl /= temperature;
    
    //std::cerr << funame << "temperature = " << temperature / Phys_const::kelv << "   acl1 = " << acl << std::endl;

    res *= std::exp(-acl);
    
    // second order correction
    //
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += 2. * std::pow(_anharm(i, i), 2.) * double(_fdegen[i] * (_fdegen[i] + 1))	* std::pow(rfactor[i], 2.)
	* (1. + double(2 * (_fdegen[i] + 1)) * rfactor[i]) * std::pow(sfactor[i], 4.);
    
      for(int j = 0; j < i; ++j)
	acl += std::pow(_anharm(i, j), 2.) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j]
	  * (1. + (double)_fdegen[i] * rfactor[i] + (double)_fdegen[j] * rfactor[j])
	  * std::pow(sfactor[i] * sfactor[j], 2.);
    
      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  if(i != j && i != k)
	    acl += 2. * _anharm(i, j) * _anharm(i, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k])
	      * rfactor[i] * rfactor[j] * rfactor[k] * std::pow(sfactor[i], 2.) * sfactor[j] * sfactor[k];
	
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 4. * _anharm(i, i) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) *_fdegen[j])
	    * std::pow(rfactor[i], 2.) * rfactor[j] * std::pow(sfactor[i], 3.) * sfactor[j];
    }
    acl /= 2. * temperature * temperature;

    //std::cerr << funame << "temperature = " << temperature / Phys_const::kelv << "   acl2 = " << acl << std::endl;

    res *= std::exp(acl);
    
    // third order correction
    //
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += 4. * std::pow(_anharm(i, i), 3.) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i], 2.) *
	(1. + double(8 * _fdegen[i] + 12) * rfactor[i] + double(8 * _fdegen[i] * _fdegen[i] + 22 * _fdegen[i] + 15) * std::pow(rfactor[i], 2.) +
	 double(2 * _fdegen[i] * _fdegen[i] + 4 * _fdegen[i] + 2) * std::pow(rfactor[i], 3.)) * std::pow(sfactor[i], 6.);

      for(int j = 0; j < i; ++j)
	acl += std::pow(_anharm(i, j), 3.) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] *
	  (1. + double(3 * _fdegen[i] + 1) * rfactor[i] + double(3 * _fdegen[j] + 1) * rfactor[j] +
	   std::pow(double(_fdegen[i]) * rfactor[i], 2.) * (1. + rfactor[j]) +
	   std::pow(double(_fdegen[j]) * rfactor[j], 2.) * (1. + rfactor[i]) +
	   double(6 * _fdegen[i] * _fdegen[j] + 3 * _fdegen[i] + 3 * _fdegen[j] + 1) * rfactor[i] * rfactor[j]) *
	  std::pow(sfactor[i] * sfactor[j], 3.);
      
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 12. * std::pow(_anharm(i, i), 2.) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j]) *
	    std::pow(rfactor[i], 2.) * rfactor[j] *
	    (1. + double(3 * _fdegen[i] + 4) * rfactor[i] + double(_fdegen[i] + 1) * std::pow(rfactor[i], 2.)) *
	    std::pow(sfactor[i], 5.) * sfactor[j]
	    + 6. * _anharm(i, i) * std::pow(_anharm(i, j), 2.) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j]) *
	    std::pow(rfactor[i], 2.) * rfactor[j] *
	    (2. + double(2 * _fdegen[i] + 1) * rfactor[i] + double(2 * _fdegen[j]) * rfactor[j] + double(_fdegen[j]) * rfactor[i] * rfactor[j]) *
	    std::pow(sfactor[i], 4.) * std::pow(sfactor[j], 2.);

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  if(i != j && i != k && k != j)
	    acl += 3. * std::pow(_anharm(i, j), 2.) * _anharm(i, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) *
	      rfactor[i] * rfactor[j] * rfactor[k] *
	      (1. + double(2 * _fdegen[i] + 1) * rfactor[i] + double(_fdegen[j]) * rfactor[j] * (1. + rfactor[i])) *
	      std::pow(sfactor[i] * sfactor[j], 2.) * sfactor[k];

      for(int j = 0; j < i; ++j)
	acl += 24. * _anharm(i, i) * _anharm(j, j) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j] * (_fdegen[j] + 1)) *
	  std::pow(rfactor[i] * rfactor[j], 2.) * std::pow(sfactor[i] * sfactor[j], 3.);

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  if(i != j && i != k)
	    acl += 12. * _anharm(i, i) * _anharm(i, j) * _anharm(i, k) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j] * _fdegen[k]) *
	      std::pow(rfactor[i], 2.) * rfactor[j] * rfactor[k] * (2. + rfactor[i]) * std::pow(sfactor[i], 4.) * sfactor[j] * sfactor[k];

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  if(i != j && i != k && k != j)
	    acl += 12. * _anharm(i, i) * _anharm(i, j) * _anharm(j, k) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j] * _fdegen[k]) *
	      std::pow(rfactor[i], 2.) * rfactor[j] * rfactor[k] * std::pow(sfactor[i], 3.) * std::pow(sfactor[j], 2.) * sfactor[k];

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  for(int l = 0; l < k; ++l)
	    if(i != j && i != k && i != l)
	      acl += 6. * _anharm(i, j) * _anharm(i, k) * _anharm(i, l) * double(_fdegen[i] * _fdegen[j] * _fdegen[k] * _fdegen[l]) *
		rfactor[i] * rfactor[j] * rfactor[k] * rfactor[l] * (1. + rfactor[i]) *
		std::pow(sfactor[i], 3.) * sfactor[j] * sfactor[k] * sfactor[l];

      // takes into account i->j, k->l symmetry
      for(int j = 0; j < i; ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  for(int l = 0; l < _frequency.size(); ++l)
	    if(i != k && i != l && j != k && j != l && k != l)
	      acl += 6. * _anharm(i, j) * _anharm(i, k) * _anharm(j, l) * double(_fdegen[i] * _fdegen[j] * _fdegen[k] * _fdegen[l]) *
		rfactor[i] * rfactor[j] * rfactor[k] * rfactor[l] * std::pow(sfactor[i] * sfactor[j], 2.) * sfactor[k] * sfactor[l];

      for(int j = 0; j < i; ++j)
	for(int k = 0; k < j; ++k)
	  acl += 6. * _anharm(i, j) * _anharm(i, k) * _anharm(j, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) *
	    rfactor[i] * rfactor[j] * rfactor[k] * (1. + _fdegen[i] * rfactor[i] + _fdegen[j] * rfactor[j] + _fdegen[k] * rfactor[k]) *
	    std::pow(sfactor[i] * sfactor[j] * sfactor[k], 2.);
    }
    acl /= 6. * std::pow(temperature, 3.);

    //std::cerr << funame << "temperature = " << temperature / Phys_const::kelv << "   acl3 = " << acl << std::endl;
												       
    res *= std::exp(-acl);
  }

  // rotational-vibrational contribution
  //
  if(_rvc.size()) {
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += _rovib(1, i) * double(_fdegen[i]) * rfactor[i] * sfactor[i];
      
      acl += _rovib(2, i, i) * double(_fdegen[i]) * rfactor[i] * std::pow(sfactor[i], 2.) * (1. + double(_fdegen[i]) * rfactor[i]);
    
      // a(i, j) = 2 * _rovib(2, i, j)
      for(int j = 0; j < i; ++j)
	acl += 2. * _rovib(2, i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] * sfactor[i] * sfactor[j];

      acl += 0.5 * std::pow(_rovib(1, i), 2.) * double(_fdegen[i]) * rfactor[i] * std::pow(sfactor[i], 2.);

      acl += _rovib(1, i) * _rovib(2, i, i) * double(_fdegen[i]) * rfactor[i] * (1. + double(1 + 2 * _fdegen[i]) * rfactor[i]) *
	std::pow(sfactor[i], 3.);

      // a(i, j) = 2 * _rovib(2, i, j)
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 2. * _rovib(1, i) * _rovib(2, i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] *
	    std::pow(sfactor[i], 2.) * sfactor[j];

      acl += std::pow(_rovib(1, i), 3.) * double(_fdegen[i]) * rfactor[i] * (1. + rfactor[i]) * std::pow(sfactor[i], 3.) / 6.;

      acl += _rovib(3, i, i, i) * double(_fdegen[i]) * rfactor[i] *
	(1. + double(1 + 3 * _fdegen[i]) * rfactor[i] + double(_fdegen[i] * _fdegen[i]) * rfactor[i] * rfactor[i]) * std::pow(sfactor[i], 3.);

      // a(i, i, j) = 3 * _rovib(3, i, i, j)
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 3. * _rovib(3, i, i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] * (1. + double(_fdegen[i]) * rfactor[i]) *
	    std::pow(sfactor[i], 2.) * sfactor[j];

      // a(i, j, k) = 6 * _rovib(3, i, j, k)
      for(int j = 0; j < i; ++j)
	for(int k = 0; k < j; ++k)
	  acl += 6. * _rovib(3, i, j, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) * rfactor[i] * rfactor[j] * rfactor[k] *
	    sfactor[i] * sfactor[j] * sfactor[k];
    }
    res *= std::exp(acl);
  }

  // mixing rotational-vibrational and anharmonic contribution
  //
  if(_rvc.size() && _anharm.isinit()) {
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += 2. * _rovib(1, i) * _anharm(i, i) * double(_fdegen[i] * (_fdegen[i] + 1))
	* std::pow(rfactor[i], 2.) * std::pow(sfactor[i], 3.);

      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += _rovib(1, i) * _anharm(i, j) * double(_fdegen[i] * _fdegen[j])
	    * rfactor[i] * rfactor[j] * std::pow(sfactor[i], 2.) * sfactor[j];

      acl += _rovib(2, i, i) * _anharm(i, i) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i], 2.) *
	(4. + double(4 * _fdegen[i] + 2) * rfactor[i]) * std::pow(sfactor[i], 4.);

      // a(i, j) = 2 * _rovib(2, i, j)
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 4. * _rovib(2, i, j) * _anharm(i, i) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i], 2.) * rfactor[j] *
	    std::pow(sfactor[i], 3.) * sfactor[j]
	    + _rovib(2, i, i) * _anharm(i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] *
	    (1. + double(2 * _fdegen[i] + 1) * rfactor[i]) * std::pow(sfactor[i], 3.) * sfactor[j];
      
      // a(i, j) = 2 * _rovib(2, i, j)
      for(int j = 0; j < i; ++j)
	acl += 2. * _rovib(2, i, j) * _anharm(i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] *
	  (1. + double(_fdegen[i]) * rfactor[i] + double(_fdegen[j]) * rfactor[j]) * std::pow(sfactor[i] * sfactor[j], 2.);

      // a(i, j) = 2 * _rovib(2, i, j)
      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  if(i != j && i != k && j != k)
	    acl += 2. * _rovib(2, i, j) * _anharm(k, j) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) * rfactor[i] * rfactor[j] * rfactor[k] *
	      sfactor[i] * std::pow(sfactor[j], 2.) * sfactor[k];
    }
    acl /= temperature;
    res *= std::exp(-acl);

    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += 4. * _rovib(1, i) * std::pow(_anharm(i, i), 2.) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i], 2.) *
	(1. + double(3 * _fdegen[i] + 4) * rfactor[i] + double(_fdegen[i] + 1) * std::pow(rfactor[i], 2.) +
	 double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i], 3.)) * std::pow(sfactor[i], 5.);

      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += _rovib(1, i) * std::pow(_anharm(i, j), 2.) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] *
	    (1. + double(2 * _fdegen[i] + 1) * rfactor[i] + double(_fdegen[j]) * rfactor[j] * (1. + rfactor[i])) *
	    std::pow(sfactor[i], 3.) * std::pow(sfactor[j], 2.)
	    + 4. * _rovib(1, i) * _anharm(i, i) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j]) *
	    std::pow(rfactor[i], 2.) * rfactor[j] * (2. + rfactor[i]) * std::pow(sfactor[i], 4.) * sfactor[j] // misprint in the paper ?
	    + 4. * _rovib(1, j) * _anharm(i, i) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j]) *
	    std::pow(rfactor[i], 2.) * rfactor[j] * std::pow(sfactor[i], 3.) * std::pow(sfactor[j], 2.);

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  if(i != j && i != k)
	    acl += 2. * _rovib(1, i) * _anharm(i, j) * _anharm(i, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) *
	      rfactor[i] * rfactor[j] * rfactor[k] * (1. + rfactor[i]) * std::pow(sfactor[i], 3.) * sfactor[j] * sfactor[k];
      
      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  if(i != j && i != k && j != k)
	    acl += 2. * _rovib(1, j) * _anharm(i, j) * _anharm(i, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) *
	      rfactor[i] * rfactor[j] * rfactor[k] * std::pow(sfactor[i] * sfactor[j], 2.) * sfactor[k];
    }
    acl /= 2. * temperature * temperature;
    res *= std::exp(acl);
  }

  return res;
}

double Model::RigidRotor::_core_weight (double temperature) const
{
  static const double nfac = std::sqrt(M_PI);
  
  double res = _factor * temperature * (1. + _rofactor * temperature);
  if(_rdim == 3)
    res *= std::sqrt(temperature * M_PI);
    
  return res;
}

double Model::RigidRotor::states (double ener) const
{
  const char funame [] = "Model::RigidRotor::states: ";

  if(mode() == NOSTATES) {
    std::cerr << funame << "wrong case\n";
    throw Error::Logic();
  }
  
  if(!_frequency.size())
    return _core_states(ener);

  ener -= ground();

  if(ener <= 0.)
    return 0.;

  if(ener >= _states.arg_max()) {
    if(ener >= _states.arg_max() * 2.)
      IO::log << IO::log_offset << funame << "WARNING: energy far beyond interpolation range\n";
    return _states.fun_max() * std::pow(ener / _states.arg_max(), _nmax);
  }

  return _states(ener);
}

double Model::RigidRotor::_core_states (double ener) const
{
  const char funame [] = "Model::RigidRotor::states: ";

  if(ener <= 0.)
    return 0.;

  switch(mode()) {
  case DENSITY:
    if(_rdim == 2)
      return _factor;
    else
      return _factor * 2. * std::sqrt(ener);
  case NUMBER:
    if(_rdim == 2)
      return _factor * ener;
    else
      return _factor * 1.3333333333333 * std::sqrt(ener) * ener;
  default:
    std::cerr << funame << "wrong case\n";
    throw Error::Logic();
  }
}

/*******************************************************************************************
 *************************** TRANSITIONAL MODES NUMBER OF STATES  **************************
 *******************************************************************************************/
  
Model::Rotd::Rotd(IO::KeyBufferStream& from, int m) 
  : Core(m)
{
  const char funame [] = "Model::Rotd::Rotd : ";

  IO::Marker funame_marker(funame);

  static const double rotd_etol = 1.e-10;
  static const double rotd_ntol = 1.;

  int    itemp;
  double dtemp;

  KeyGroup RotdModel;

  Key rotd_key("File"          );
  Key symm_key("SymmetryFactor");

  std::string rotd_name;
  double factor = 1.;

  bool issym = false;

  std::string comment, token;
  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // symmetry factor (number of symmetry operations)
    else if(symm_key == token) {
      if(issym) {
	std::cerr << funame << "symmetry number has been initialized already\n";
	throw Error::Init();
      }

      issym = true;

      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(dtemp <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      factor /= dtemp;
    }
    // rotd output file name
    else if(rotd_key == token) {
      if(rotd_name.size()) {
	std::cerr << funame << token <<  ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> rotd_name)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
 
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(!rotd_name.size()) {
    std::cerr << funame << "rotd output file name is not initialized\n";
    throw Error::Init();
  }

  /************************************* ROTD INPUT **********************************/

  std::ifstream rotd_in(rotd_name.c_str());
  if(!rotd_in) {
    std::cerr << funame << "cannot open rotd input file " << rotd_name << "\n";
    throw Error::Open();
  }

  std::map<double, double> read_nos;

  while(rotd_in >> dtemp) {
    dtemp *= Phys_const::incm; // energy
    rotd_in >> read_nos[dtemp];             // number of states
    if(!rotd_in) {
      std::cerr << funame << "reading transitional modes density failed\n";
      throw Error::Input();
    }
  }

  // find zero density energy
  typedef std::map<double, double>::const_reverse_iterator It;
  It izero;
  for( izero = read_nos.rbegin(), itemp = 1; izero != read_nos.rend(); ++izero, ++itemp)
    if(izero->second <= rotd_ntol)
      break;

  if(itemp < 3) {
    std::cerr << funame << "not enough data\n";
    throw Error::Init();
  }

  if(izero != read_nos.rend())
    _ground = izero->first;
  else {
    std::cerr << funame << "zero flux energy is not reached\n";
    throw Error::Input();
  }

  _rotd_ener.resize(itemp);
  _rotd_nos.resize(itemp);

  _rotd_ener[0] = 0.;
  _rotd_nos[0]  = 0.;

  IO::log << IO::log_offset << "grid size = " << itemp << "\n";

  --itemp;
  for(It it = read_nos.rbegin(); it != izero; ++it, --itemp) {
    _rotd_ener[itemp] = it->first - _ground;
    _rotd_nos[itemp]  = it->second * factor;
  }
  
  Array<double> x(_rotd_ener.size() - 1);
  Array<double> y(_rotd_ener.size() - 1);

  for(int i = 1; i < _rotd_ener.size(); ++i) {
    itemp = i - 1;
    x[itemp] = std::log(_rotd_ener[i]);
    y[itemp] = std::log( _rotd_nos[i]);
  }

  _rotd_spline.init(x, y, x.size());

  _rotd_emin = _rotd_ener[1]     * (1. + rotd_etol);
  _rotd_emax = _rotd_ener.back() * (1. - rotd_etol);

  _rotd_nmin = (y[1] - y[0]) / (x[1] - x[0]);
  _rotd_amin = std::exp(y[0] - x[0] * _rotd_nmin);

  int l1 = x.size() - 1;
  int l2 = x.size() - 2;
  _rotd_nmax = (y[l1] - y[l2]) / (x[l1] - x[l2]);
  _rotd_amax = std::exp(y[l1] - x[l1] * _rotd_nmax);

  IO::log << IO::log_offset << "effective power exponent at " 
	  << _rotd_emax / Phys_const::kcal << " kcal/mol = " << _rotd_nmax << "\n"
	  << IO::log_offset << "ground energy [kcal/mol] = " << _ground / Phys_const::kcal << "\n";

}// Rotd Core

Model::Rotd::~Rotd ()
{
  //std::cout << "Model::Rotd destroyed\n";
}

double Model::Rotd::ground () const
{
  return _ground;
}

double Model::Rotd::states (double ener) const 
{
  const char funame [] = "Model::Rotd::states: ";
  double dtemp;

  if(ener <= 0.)
    return 0.;

  switch(mode()) {
  case NUMBER:
    if(ener <= _rotd_emin)
      return _rotd_amin * std::pow(ener, _rotd_nmin);
    else if(ener >= _rotd_emax)
      return _rotd_amax * std::pow(ener, _rotd_nmax);
    else
      return std::exp(_rotd_spline(std::log(ener)));
  case DENSITY:
    if(ener <= _rotd_emin)
      return _rotd_amin * _rotd_nmin * std::pow(ener, _rotd_nmin - 1.);
    else if(ener >= _rotd_emax)
      return _rotd_amax * _rotd_nmax * std::pow(ener, _rotd_nmax - 1.);
    else {
      dtemp = std::log(ener);
      return std::exp(_rotd_spline(dtemp, 0)) * _rotd_spline(dtemp, 1) / ener;
    }
  default:
    std::cerr << funame << "wrong mode\n";
    throw Error::Logic();
  }
}

double Model::Rotd::weight (double temperature) const
{
  const char funame [] = "Model::Rotd::weight: ";

  const double eps = 1.e-3;

  double dtemp;

  Array<double> term(_rotd_nos.size());
  term[0] = 0.;
  for(int i = 1; i < _rotd_nos.size(); ++i)
    term[i] = _rotd_nos[i] / std::exp(_rotd_ener[i] / temperature);

  int info;
  double res;
  davint_(_rotd_ener, term, term.size(), _rotd_ener.front(), _rotd_ener.back(), res, info); 
  if (info != 1) {
    std::cerr << funame  << "davint integration error\n";
    throw Error::Logic();
  }
  res /= temperature;

  dtemp = _rotd_emax / temperature;
  if(dtemp <= _rotd_nmax) {
    IO::log << IO::log_offset << funame << "WARNING: " 
	    << "integration cutoff energy is less than weight maximum energy\n";
    return res;
  }

  dtemp = _rotd_amax * std::pow(_rotd_emax, _rotd_nmax) / std::exp(dtemp) / (1. - _rotd_nmax / dtemp);
  if(dtemp / res > eps)
    IO::log << IO::log_offset << funame << "WARNING: integration cutoff error = " << dtemp / res << "\n";    
  res += dtemp;

  return res;
}

/********************************************************************************************
 ******************************* INTERNAL ROTATION FOR MULTIROTOR ***************************
 ********************************************************************************************/

Model::InternalRotation::InternalRotation (IO::KeyBufferStream& from) 
  : InternalRotationBase(from), _msize(0), _wsize(0), _psize(0), _qmin(0), _qmax(0)
{
  const char funame [] = "Model::InternalRotation::InternalRotation: ";

  IO::Marker funame_marker("InternalRotation", IO::Marker::NOTIME);

  KeyGroup InternalRotationDefinition;

  Key msize_key("MassExpansionSize"     );
  Key psize_key("PotentialExpansionSize");
  Key  qmin_key("HamiltonSizeMin"       );
  Key  qmax_key("HamiltonSizeMax"       );
  Key wsize_key("GridSize"              );

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  std::string token, line, comment;

  while(from >> token) {
    //
    // input end
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);

      break;
    }
    // mass fourier expansion size
    //
    else if(msize_key == token) {
      //
      if(!(from >> _msize)) {
	//
	std::cerr << funame << token << ": corruped\n";

	throw Error::Input();
      }
      std::getline(from, comment);

      if(_msize < 1) {
	//
	std::cerr << funame << token << ": should be positive and odd\n";

	throw Error::Range();
      }

      if(!(_msize % 2)) {
	//
	++_msize;

	IO::log << IO::log_offset << token << ": WARNING: should be odd: changing to " << _msize << "\n";
      }
    }
    // potential fourier expansion size
    //
    else if(psize_key == token) {
      //
      if(_psize) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Init();
      }
      if(!(from >> _psize)) {
	//
	std::cerr << funame << token << ": corruped\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_psize < 1 || !(_psize % 2)) {
	//
	std::cerr << funame << token << ": should be positive and odd\n";

	throw Error::Range();
      }

      if(!(_psize % 2)) {
	//
	++_psize;

	IO::log << IO::log_offset << token << ": WARNING: should be odd: changing to " << _psize << "\n";
      }
    }
    // minimal quantum state size
    //
    else if(qmin_key == token) {
      //
      if(_qmin) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Init();
      }

      if(!(from >> _qmin)) {
	//
	std::cerr << funame << token << ": corruped\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_qmin < 1) {
	//
	std::cerr << funame << token << ": should be positive and odd\n";

	throw Error::Range();
      }

      if(!(_qmin % 2)) {
	//
	++_qmin;

	IO::log << IO::log_offset << token << ": WARNING: should be odd: changing to " << _qmin << "\n";
      }
    }
    // maximal quantum state size
    //
    else if(qmax_key == token) {
      //
      if(_qmax) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Init();
      }

      if(!(from >> _qmax)) {
	//
	std::cerr << funame << token << ": corruped\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_qmax < 1) {
	//
	std::cerr << funame << token << ": should be positive and odd\n";

	throw Error::Range();
      }

      if(!(_qmax % 2)) {
	//
	++_qmax;

	IO::log << IO::log_offset << token << ": WARNING: should be odd: changing to " << _qmax << "\n";
      }
    }
    // grid size
    //
    else if(wsize_key == token) {
      //
      if(!(from >> _wsize)) {
	//
	std::cerr << funame << token << ": corruped\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_wsize < 1) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
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
    }//
    //
  }// input cycle
  
  /************************************ CHECKING ***************************************/

  if(!from) {
    //
    std::cerr << funame << "input is corrupted\n";

    throw Error::Input();
  }

  if(!_msize) {
    //
    std::cerr << funame << "mass expansion size not initialized\n";
    
    throw Error::Init();
  }

  if(!_wsize) {
    //
    std::cerr << funame << "angular grid size not initialized\n";

    throw Error::Init();
  }

  //  if(!_qmin || !_qmax || _qmin > _qmax) {
  //
  //std::cerr << funame << "harmonic expansion size limits should be positive, odd, and ordered: " << _qmin << ", " << _qmax << "\n";
  
  //throw Error::Range();
  //
  //}//
  //
}// InternalRotation

Model::InternalRotation::~InternalRotation ()
{
  //std::cout << "Model::InternalRotation destroyed\n";
}

/********************************************************************************************
 ***************************** MULTIPLE COUPLED ROTORS MODEL ********************************
 ********************************************************************************************/

Model::MultiRotor::MultiRotor(IO::KeyBufferStream& from, const std::vector<Atom>& atom, int mm)  

  : Core(mm), _extra_ener(-1.), _extra_step(0.1), _ener_quant(Phys_const::incm), _level_ener_max(-1.),

    _with_ext_rot(true), _full_quantum_treatment(false), _amom_max(0), _ptol(-1.), _vtol(-1.), _mtol(-1.),

    _pot_shift(0.), _external_symmetry(1.), _with_ctf(false)
{
  const char funame [] = "Model::MultiRotor::MultiRotor: ";

  IO::Marker funame_marker(funame);

  KeyGroup MultiRotorModel;

  Key      irot_key("InternalRotation"                );
  Key      symm_key("SymmetryFactor"                  );
  Key  kcal_pes_key("PotentialEnergySurface[kcal/mol]");
  Key  incm_pes_key("PotentialEnergySurface[1/cm]"    );
  Key    kj_pes_key("PotentialEnergySurface[kJ/mol]"  );
  Key kcal_qmax_key("QuantumLevelEnergyMax[kcal/mol]" );  
  Key incm_qmax_key("QuantumLevelEnergyMax[1/cm]"     );  
  Key   kj_qmax_key("QuantumLevelEnergyMax[kJ/mol]"   );  
  Key kcal_emax_key("InterpolationEnergyMax[kcal/mol]");  
  Key incm_emax_key("InterpolationEnergyMax[1/cm]"    );  
  Key   kj_emax_key("InterpolationEnergyMax[kJ/mol]"  );  
  Key     estep_key("InterpolationEnergyStep[1/cm]"   );  
  Key     extra_key("ExtrapolationStep"               );  
  Key      erot_key("WithoutExternalRotation"         );  
  Key      full_key("FullQuantumTreatment"            );
  Key      ptol_key("PotentialTolerance"              );
  Key      vtol_key("FrequencyTolerance"              );
  Key      mtol_key("MobilityTolerance"               );
  Key      amom_key("AngularMomentumMax"              );
  Key     force_key("ForceQFactor"                    );
  Key       ctf_key("WithCurvlinearFactor"            );

  bool  force_qfactor = false;
  
  int             itemp;
  double          dtemp;
  bool            btemp;
  std::string     stemp;
  Lapack::Vector  vtemp;
  Lapack::Matrix  mtemp;

  std::vector<int>    ivec;
  std::vector<double> dvec;

  std::string token, line, comment;

  std::vector<std::vector<double> > vibration_sampling; // vibrational sampling data

  while(from >> token) {
    //
    // input end
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);

      break;
    }
    // internal rotation factor included
    //
    else if(ctf_key == token) {
      //
      std::getline(from, comment);

      _with_ctf = true;
    }
    // force qfactor
    //
    else if(force_key == token) {
      //
      std::getline(from, comment);

      force_qfactor = true;
    }
    // internal motion definition
    //
    else if(irot_key == token) {
      //
      std::getline(from, comment);

      IO::log << IO::log_offset << _internal_rotation.size() + 1 << "-th internal rotation definition:\n";

      _internal_rotation.push_back(InternalRotation(from));
    }
    // potential and vibrational samplings
    //
    else if(kcal_pes_key == token || incm_pes_key == token || kj_pes_key == token) {
      //
      if(_pot_index.rank()) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Init();
      }

      // file reading
      //
      IO::LineInput name_input(from);

      if(!(name_input >> stemp)) {
	//
	std::cerr << funame << token << ": cannot read potential energy surface file name\n";

	throw Error::Input();
      }

      std::ifstream file_input(stemp.c_str());

      if(!file_input) {
	//
	std::cerr << funame << token << ": cannot open potential energy surface file " << stemp << "\n";

	throw Error::Input();
      }

      // sampling dimensions      
      //
      ivec.clear();

      IO::LineInput sampling_dimensions_input(file_input);

      while(sampling_dimensions_input >> itemp) {
	//
	if(itemp < 1) {
	  //
	  std::cerr << funame << token << ": the sampling dimensions should be positive\n";

	  throw Error::Range();
	}

	ivec.push_back(itemp);
      }

      _pot_index.resize(ivec);

      _pot_real.resize(_pot_index.size());

      // vibrational frequencies  subset
      //
      std::set<int> vib_set;

      IO::LineInput vib_set_input(file_input);

      while(vib_set_input >> itemp) {
	//
	if(itemp < 1) {
	  //
	  std::cerr << funame << token << ": vibration index should be positive\n";

	  throw Error::Init();
	}

	if(!vib_set.insert(itemp).second) {
	  //
	  std::cerr << funame << token << ": identical vibration indices\n";

	  throw Error::Init();
	}
      }
      
      if(vib_set.size()) {
	//
	_vib_four.resize(vib_set.size());

	_eff_vib_four.resize(vib_set.size());

	vibration_sampling.resize(_pot_index.size());
      }

      // headline
      //
      std::getline(file_input, comment);

      std::set<int> pl_pool;

      // sampling reading
      //
      for(int pl = 0; pl < _pot_index.size(); ++pl) {
	//
	IO::LineInput sampling_input(file_input);

	// reading potential sampling index
	//
	for(int i = 0; i < ivec.size(); ++i) {
	  //
	  if(!(sampling_input >> itemp)) {
	    //
	    std::cerr << funame << token 
		      << ": potential data corrupted: line input format: i_1 ... i_n energy f_1 ... f_m\n";

	    throw Error::Input();
	  }

	  // Fortran style indexing in input
	  //
	  if(itemp < 1 || itemp > _pot_index.size(i)) {
	    //
	    std::cerr << funame << token 
		      << ": " << i + 1 
		      << "-th potential sampling index (Fortran style indexing) out of range\n";

	    throw Error::Range();
	  }

	  ivec[i] = itemp - 1;
	}

	// potential sampling linear index
	//
	int pl_curr = _pot_index(ivec);

	if(!pl_pool.insert(pl_curr).second) {
	  //
	  std::cerr << funame << token << ": identical sampling: (";

	  for(int i = 0; i < ivec.size(); ++i) {
	    //
	    if(i)
	      //
	      std::cerr << ", ";

	    std::cerr << ivec[i] + 1;	
	  }
	  std::cerr << ")\n";

	  throw Error::Init();
	}

	// reading potential energy
	//
	if(!(sampling_input >> dtemp)) {
	  //
	  std::cerr << funame << token 
		    << ": potential data corrupted: line input format: i_1 ... i_n energy f_1 .. f_m\n";

	  throw Error::Input();
	}

	if(kcal_pes_key == token)
	  //
	  dtemp *= Phys_const::kcal;

	if(incm_pes_key == token)
	  //
	  dtemp *= Phys_const::incm;

	if(kj_pes_key == token)
	  //
	  dtemp *= Phys_const::kjoul;
	
	_pot_real[pl_curr] = dtemp;

	// reading vibrational frequencies
	//
	if(_vib_four.size()) {
	  //
	  itemp = 0;

	  while(sampling_input >> dtemp) {
	    //
	    dtemp *= Phys_const::incm;

	    ++itemp;

	    if(vib_set.find(itemp) != vib_set.end()) {
	      //
	      if(dtemp <= 0.) {
		//
		std::cerr << funame << token << ": ("; 

		for(int i = 0; i < ivec.size(); ++i) {
		  //
		  if(i)
		    //
		    std::cerr << ", ";

		  std::cerr << ivec[i] + 1;
		}
		std::cerr << ")-th sampling: " << itemp << "-th frequency negative\n";

		throw Error::Range();
	      }

	      vibration_sampling[pl_curr].push_back(dtemp);
	    }	  
	  }

	  // checking
	  //
	  if(vibration_sampling[pl_curr].size() != _vib_four.size()) {
	    //
	    std::cerr << funame << token << ": (";

	    for(int i = 0; i < ivec.size(); ++i) {
	      //
	      if(i)
		//
		std::cerr << ", ";

	      std::cerr << ivec[i] + 1;
	    }
	    std::cerr << ")-th sampling: wrong frequencies # = "
		      << vibration_sampling[pl_curr].size() 
		      << ", total read freq. # = " << itemp << "\n";

	    throw Error::Range();
	  }
	  //
	}// reading vibrational frequencies
	//
      }// reading sampling
    }
    // external rotational symmetry
    //
    else if(symm_key == token) {
      //
      if(!(from >> _external_symmetry)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }
      std::getline(from, comment);

      if(_external_symmetry <= 0.) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }
    }
    // include external rotation
    //
    else if(erot_key == token) {
      //
      _with_ext_rot = false;
      
      std::getline(from, comment);
    }
    // full quantum treatment
    //
    else if(full_key == token) {
      //
      _full_quantum_treatment = true;

      std::getline(from, comment);
    }
    // maximum quantum level energy relative to the potential minimum
    //
    else if(kcal_qmax_key == token || incm_qmax_key == token || kj_qmax_key == token) {
      //
      if(!(from >> _level_ener_max)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(kcal_qmax_key == token)
	//
	_level_ener_max *= Phys_const::kcal;

      if(incm_qmax_key == token)
	//
	_level_ener_max *= Phys_const::incm;

      if(kj_qmax_key == token)
	//
	_level_ener_max *= Phys_const::kjoul;
    }
    // interpolation maximal energy
    //
    else if(kcal_emax_key == token || incm_emax_key == token || kj_emax_key == token) {
      //
      if(!(from >> _extra_ener)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_extra_ener <= 0.) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }

      if(kcal_emax_key == token)
	//
	_extra_ener *= Phys_const::kcal;

      if(incm_emax_key == token)
	//
	_extra_ener *= Phys_const::incm;

      if(kj_emax_key == token)
	//
	_extra_ener *= Phys_const::kjoul;
    }
    // extrapolation step
    //
    else if(extra_key == token) {
      //
      if(!(from >> _extra_step)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_extra_step <= 0. || _extra_step >= 1.) {
	//
	std::cerr << funame << token << ": out of range\n";

	throw Error::Range();
      }
    }
    // frequency quantum
    //
    else if(estep_key == token) {
      //
      if(!(from >> _ener_quant)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_ener_quant <= 0.) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }

      _ener_quant *= Phys_const::incm;
    }
    // angular momentum quantum number maximum
    //
    else if(amom_key == token) {
      //
      if(!(from >> _amom_max)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_amom_max <= 0) {
	//
	std::cerr << funame << token << ": should be positive\n";

	throw Error::Range();
      }
    }
    // potential matrix element tolerance
    //
    else if(ptol_key == token) {
      //
      if(!(from >> _ptol)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_ptol <= 0. || _ptol >= 0.5) {
	//
	std::cerr << funame << token << ": out of range\n";

	throw Error::Range();
      }
    }
    // vibrational frequency matrix element tolerance
    //
    else if(vtol_key == token) {
      //
      if(!(from >> _vtol)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_vtol <= 0. || _vtol >= 0.5) {
	//
	std::cerr << funame << token << ": out of range\n";

	throw Error::Range();
      }
    }
    // mass matrix element tolerance
    //
    else if(mtol_key == token) {
      //
      if(!(from >> _mtol)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(_mtol <= 0. || _mtol >= 0.5) {
	//
	std::cerr << funame << token << ": out of range\n";

	throw Error::Range();
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

  /***************************************** CHECKING ****************************************/

  if(!from) {
    //
    std::cerr << funame << "input stream corrupted\n";

    throw Error::Input();
  }

  // number of internal motions
  //
  if(!internal_size()) {
    //
    std::cerr << funame << "internal motions are not defined\n";

    throw Error::Init();
  }
    
  // initial geometry
  //
  if(!atom.size()) {
    //
    std::cerr << funame << "geometry is not initialized\n";

    throw Error::Init();
  }

  // potential sampling dimensionality
  //
  if(_pot_index.rank() != internal_size()) {
    //
    std::cerr << funame << "potential dimensionality is inconsistent with the number of internal motions\n";

    throw Error::Init();
  }

  if(_extra_ener < 0. && (mode() != NOSTATES || force_qfactor)) {
    //
    std::cerr << funame << "Maximum interpolation energy not initialized\n";

    throw Error::Init();
  }

  /******************************* INITIAL SETTINGS *******************************/
  
  // zero-energy shift
  //
  for(int v = 0; v < _vib_four.size(); ++v)
    //
    for(int g = 0; g < _pot_index.size(); ++g)
      //
      _pot_real[g] += vibration_sampling[g][v] / 2.;

  if(_ptol > 0. && _vtol < 0.)
    //
    _vtol = _ptol;
  

  /******************************* MASS FOURIER EXPANSION *******************************/
  
  if(1) {
    //
    IO::Marker work_marker("mass fourier expansion");

    // mass fourier expansion dimensions
    //
    ivec.resize(internal_size());

    for(int r = 0; r < internal_size(); ++r) 
      //
      ivec[r] = _internal_rotation[r].mass_fourier_size();

    _mass_index.resize(ivec);

    IO::log << IO::log_offset << "mass fourier expansion size = " << _mass_index.size() << std::endl;

    for(int h = 0; h < _mass_index.size(); ++h) {
      //
      _imm_four[h].resize(internal_size());

      _imm_four[h] = 0.;

      if(_with_ctf) {
	//
	_eff_imm_four[h].resize(internal_size());

	_eff_imm_four[h] = 0.;
      }
    }

    // external rotation quantized  
    //
    if(_amom_max) {
      //
      _internal_mobility_real.resize(_mass_index.size());

      _external_mobility_real.resize(_mass_index.size());

      _coriolis_coupling_real.resize(_mass_index.size());

      _ctf_real.resize(_mass_index.size());
    }

    // sampling geometries
    //
    for(int g = 0; g < _mass_index.size(); ++g) {// grid cycle
      //
      std::vector<int> gv = _mass_index(g);

      // new geometry
      //
      std::vector<Atom> atom_current = atom;

      if(g)
	//
	for(int r = 0; r < internal_size(); ++r) {
	  //
	  dtemp =  2. * M_PI * double(gv[r]) / double(symmetry(r) * _mass_index.size(r));// rotation angle

	  atom_current = _internal_rotation[r].rotate(atom_current, dtemp);
	}

      // shift center of mass to zero
      //
      shift_cm_to_zero(atom_current);

      // coriolis coupling matrix initialization
      //
      if(_amom_max)
	//
	_coriolis_coupling_real[g].resize(3, internal_size());

      // normal modes
      //
      std::vector<std::vector<D3::Vector> > nm(internal_size());

      for(int r = 0; r < internal_size(); ++r) {
	//
	// normal mode and coriolis coupling
	//
	if(_amom_max) {
	  //
	  nm[r] = _internal_rotation[r].normal_mode(atom_current, &vtemp);

	  _coriolis_coupling_real[g].column(r) = vtemp;
	}    
	// normal mode only
	//
	else {
	  //
	  nm[r] = _internal_rotation[r].normal_mode(atom_current);
	}
      }

      // generalized mass matrix
      //
      Lapack::SymmetricMatrix gmm(internal_size());

      for(int i = 0; i < internal_size(); ++i) {
	//
	// generalized mass matrix
	//
	for(int j = i; j < internal_size(); ++j) {
	  //
	  gmm(i, j) = 0.;

	  for(int a = 0; a < atom.size(); ++a)
	    //
	    gmm(i, j) += atom[a].mass() * vdot(nm[i][a], nm[j][a]);
	}
      }

      // output: one-dimensional internal mobility for original geometry
      //
      if(!g) {
	//
	IO::log << IO::log_offset << "1D internal mobilities for original geometry, 1/cm:";

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(12) << 0.5 / gmm(r, r) / Phys_const::incm;

	IO::log << "\n";

	IO::log << IO::log_offset << "generalized mass matrix [amu * bohr * bohr] for original geometry:\n";

	IO::log << IO::log_offset << std::setw(3) << "#\\#";

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(12) << r + 1;

	IO::log << "\n";

	for(int r = 0; r < internal_size(); ++r) {
	  //
	  IO::log << IO::log_offset 
		  << std::left <<  std::setw(3) << r + 1 << std::right;

	  for(int s = 0; s <= r; ++s)
	    //
	    IO::log << std::setw(12) << gmm(r, s) / Phys_const::amu;

	  IO::log << "\n";
	}

	Lapack::Vector pim = inertia_moment_matrix(atom_current).eigenvalues();

	_rotational_constant.resize(pim.size());

	IO::log << IO::log_offset << "external rotational constants for original geometry [1/cm]:";

	for(int i = 0; i < pim.size(); ++i) {
	  //
	  IO::log << "   ";

	  if(pim[i] > 0.) {
	    //
	    dtemp = 0.5 / pim[i];

	    _rotational_constant[i] = dtemp;

	    IO::log << dtemp / Phys_const::incm;
	  }
	  else {
	    //
	    std::cerr << funame << "molecule is linear?\n";

	    throw Error::Range();
	  }
	}

	IO::log << "\n";
      }

      // Cholesky representation
      //
      Lapack::Cholesky gmm_chol(gmm);

      // external rotation factor - sqrt(inertia moments product)
      //
      const double  erf = Lapack::Cholesky(inertia_moment_matrix(atom_current)).det_sqrt();

      // internal rotation factor
      //
      const double irf = gmm_chol.det_sqrt();

      // curvlinear transformation factor
      //
      const double ctf = irf * erf;

      // internal mobility matrix
      //
      Lapack::SymmetricMatrix imm = gmm_chol.invert();

      imm /= 2.;

      Lapack::SymmetricMatrix eff_imm;
      
      if(_with_ctf) {
	//
	eff_imm = imm.copy();
	
	eff_imm *= ctf;
      }
      
      // output: internal mobility matrix for original geometry
      //
      if(!g) {
	//
	IO::log << IO::log_offset << "internal mobility matrix for original geometry, 1/cm:\n";

	IO::log << IO::log_offset << std::setw(3) << "#\\#";

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(12) << r + 1;

	IO::log << "\n";

	for(int r = 0; r < internal_size(); ++r) {
	  //
	  IO::log << IO::log_offset 
		  << std::left <<  std::setw(3) << r + 1 << std::right;

	  for(int s = 0; s <= r; ++s)
	    //
	    IO::log << std::setw(12) << imm(r, s) / Phys_const::incm;

	  IO::log << "\n";
	}
      }

      if(internal_size() == 1) {
	//
	if(!g)
	  //
	  IO::log << IO::log_offset << "(single rotor) internal mobility, 1/ cm:\n";

	IO::log << IO::log_offset << std::setw(3) << g << std::setw(15) << imm(0, 0) / Phys_const::incm << "\n";
      }

      // coriolis coupling with external (overall) rotations
      //
      if(_amom_max) {
	//
	_internal_mobility_real[g] = imm.copy();

	Lapack::Matrix cc = _coriolis_coupling_real[g] * _internal_mobility_real[g];

	_external_mobility_real[g] = inertia_moment_matrix(atom_current).invert();
	
	_external_mobility_real[g] /= 2.;

	_external_mobility_real[g] += Lapack::SymmetricMatrix(cc * _coriolis_coupling_real[g].transpose());

	_coriolis_coupling_real[g] = cc;

	if(_with_ctf) {
	  //
	  _internal_mobility_real[g] *= ctf;

	  _external_mobility_real[g] *= ctf;

	  _coriolis_coupling_real[g] *= ctf;

	  _ctf_real[g] = ctf;
	}
      }

      for(int h = 0; h < _mass_index.size(); ++h) {// fourier expansion cycle
	//
	std::vector<int> hv = _mass_index(h);
	
	double dfac = 1.;
       
	for(int r = 0; r < internal_size(); ++r) {
	  //
	  if(!hv[r])
	    //
	    continue;
	  
	  if(hv[r] % 2) {
	    //
	    dfac *= std::sin(M_PI * double((hv[r] + 1) * gv[r])/ double(_mass_index.size(r)));
	  }
	  else
	    //
	    dfac *= std::cos(M_PI * double(hv[r] * gv[r])/ double(_mass_index.size(r)));
	}

	for(int i = 0; i < internal_size(); ++i) {
	  //
	  for(int j = i; j < internal_size(); ++j) {
	    //
	    //#pragma omp atomic
	    //
	    _imm_four[h](i, j) += imm(i, j) * dfac;

	    if(_with_ctf) {
	      //
	      _eff_imm_four[h](i, j) += eff_imm(i, j) * dfac;
	      //
	    }//
	    //
	  }//
	  //
	}//

	
	if(_with_ctf)
	  //
	  _eff_erf_four[h] += erf * ctf * dfac;

	_ctf_four[h] += ctf * dfac;
	
	_irf_four[h] += irf * dfac;

	_erf_four[h] += erf * dfac;
	
	//
      }// fourier expansion cycle
      //
    }// grid cycle    
  
    // normalization
    //
    for(int h = 0; h < _mass_index.size(); ++h) {// fourier expansion cycle
      //
      std::vector<int> hv = _mass_index(h);

      // normalization factor
      //
      double dfac = double(_mass_index.size());

      for(int r = 0; r < internal_size(); ++r)
	//
	if(hv[r])
	  //
	  dfac /= 2.;
    
      _imm_four[h] /= dfac;	  

      if(_with_ctf) {
	//
	_eff_imm_four[h] /= dfac;

	_eff_erf_four[h] /= dfac;
      }

      _ctf_four[h] /= dfac;

      _erf_four[h] /= dfac;
      
      _irf_four[h] /= dfac;
      //
    }// fourier expansion cycle

    // output: averaged internal rotational constant matrix
    //
    IO::log << IO::log_offset << "averaged internal mobility matrix, 1/cm:\n";

    IO::log << IO::log_offset 
	    << std::setw(3) << "#\\#";

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << std::setw(12) << r + 1;

    IO::log << "\n";

    for(int r = 0; r < internal_size(); ++r) {
      //
      IO::log << IO::log_offset 
	      << std::left <<  std::setw(3) << r + 1 << std::right;

      for(int s = 0; s <= r; ++s)
	//
	IO::log << std::setw(12) << _imm_four[0](r, s) / Phys_const::incm;

      IO::log << "\n";
    }

    IO::log << std::flush;

    // internal mobility matrix fourier expansion pruning
    //
    if(_mtol > 0.) {
      //
      IO::Marker  prune_marker("pruning internal mobility matrix fourier expansion");
	
      IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

      IO::log << IO::log_offset << "original internal mobility matrix fourier expansion size = " 
	      << _imm_four.size() << "\n";

      typedef std::map<int, Lapack::SymmetricMatrix>::iterator mit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;
      
      std::multimap<double, int>::const_iterator otit;

      for(mit_t mit = _imm_four.begin(); mit != _imm_four.end(); ++mit) {// fourier expnasion cycle
	//
	dtemp = 0.;

	for(int i = 0; i < internal_size(); ++i) {
	  //
	  for(int j = i; j < internal_size(); ++j) {
	    //
	    if(j != i) {
	      //
	      dtemp += 2. * mit->second(i, j) * mit->second(i, j);
	    }
	    else
	      //
	      dtemp += mit->second(i, j) * mit->second(i, j);
	  }
	}

	std::vector<int> gv = _mass_index(mit->first);
	  
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(gv[i])
	    //
	    itemp *= 2;
	  
	dtemp /= (double)itemp;

	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, mit->first));
	//
      }// fourier expansion cycle

      if(threshold == 0.) {
	//
	std::cerr << funame << "no kinetic energy term\n";

	throw Error::Logic();
      }

      threshold *= _mtol;
      
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;

	if(dtemp > threshold)
	  //
	  break;

	_imm_four.erase(otit->second);
      }
      
      IO::log << IO::log_offset << "pruned   internal mobility matrix fourier expansion size = " 
	      << _imm_four.size() << "\n";
      //
    }// internal mobility matrix fourier expansion pruning

    // effective internal mobility matrix fourier expansion pruning
    //
    if(_mtol > 0. && _with_ctf) {
      //
      IO::Marker  prune_marker("pruning effective internal mobility matrix fourier expansion");
	
      IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

      IO::log << IO::log_offset << "original effective internal mobility matrix fourier expansion size = " 
	      << _eff_imm_four.size() << "\n";

      typedef std::map<int, Lapack::SymmetricMatrix>::iterator mit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;
      
      std::multimap<double, int>::const_iterator otit;

      for(mit_t mit = _eff_imm_four.begin(); mit != _eff_imm_four.end(); ++mit) {// fourier expnasion cycle
	//
	dtemp = 0.;

	for(int i = 0; i < internal_size(); ++i) {
	  //
	  for(int j = i; j < internal_size(); ++j) {
	    //
	    if(j != i) {
	      //
	      dtemp += 2. * mit->second(i, j) * mit->second(i, j);
	    }
	    else
	      //
	      dtemp += mit->second(i, j) * mit->second(i, j);
	  }
	}

	std::vector<int> gv = _mass_index(mit->first);
	  
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(gv[i])
	    //
	    itemp *= 2;
	  
	dtemp /= (double)itemp;

	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, mit->first));
	//
      }// fourier expansion cycle

      if(threshold == 0.) {
	//
	std::cerr << funame << "no kinetic energy term\n";

	throw Error::Logic();
      }

      threshold *= _mtol;
      
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;

	if(dtemp > threshold)
	  //
	  break;

	_eff_imm_four.erase(otit->second);
      }
      
      IO::log << IO::log_offset << "pruned   effective internal mobility matrix fourier expansion size = " 
	      << _eff_imm_four.size() << "\n";
      //
    }// effective internal mobility matrix fourier expansion pruning

    // external rotational factor pruning
    //
    if(_mtol > 0.) {
      //
      IO::Marker  prune_marker("pruning external rotation factor fourier expansion");
      
      IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

      IO::log << IO::log_offset << "original external rotation factor fourier expansion size = " 
	      << _erf_four.size() << "\n";

      typedef std::map<int, double>::iterator pit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;
      
      std::multimap<double, int>::const_iterator otit;

      for(pit_t pit = _erf_four.begin(); pit != _erf_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _mass_index(pit->first);
	  
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(pv[i])
	    //
	    itemp *= 2;
	  
	dtemp = pit->second * pit->second / (double)itemp;
	
	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, pit->first));

      }// fourier expansion cycle

      if(threshold == 0.) {
	//
	std::cerr << funame << "no external rotation\n";

	throw Error::Logic();
      }

      threshold *= _mtol;
      
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;

	if(dtemp > threshold)
	  //
	  break;

	_erf_four.erase(otit->second);
      }

      IO::log << IO::log_offset << "pruned   external rotation factor fourier expansion size = " 
	      << _erf_four.size() << "\n";
      //
    }// external rotation factor pruning

    // effective external rotational factor pruning
    //
    if(_mtol > 0. && _with_ctf) {
      //
      IO::Marker  prune_marker("pruning effective external rotation factor fourier expansion");
      
      IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

      IO::log << IO::log_offset << "original effective external rotation factor fourier expansion size = " 
	      << _eff_erf_four.size() << "\n";

      typedef std::map<int, double>::iterator pit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;
      
      std::multimap<double, int>::const_iterator otit;

      for(pit_t pit = _eff_erf_four.begin(); pit != _eff_erf_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _mass_index(pit->first);
	  
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(pv[i])
	    //
	    itemp *= 2;
	  
	dtemp = pit->second * pit->second / (double)itemp;
	
	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, pit->first));

      }// fourier expansion cycle

      if(threshold == 0.) {
	//
	std::cerr << funame << "no external rotation\n";

	throw Error::Logic();
      }

      threshold *= _mtol;
      
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;

	if(dtemp > threshold)
	  //
	  break;

	_eff_erf_four.erase(otit->second);
      }

      IO::log << IO::log_offset << "pruned   effective external rotation factor fourier expansion size = " 
	      << _eff_erf_four.size() << "\n";
      //
    }// effective external rotation factor pruning

    // internal rotational factor pruning
    //
    if(_mtol > 0.) {
      //
      IO::Marker  prune_marker("pruning internal rotation factor fourier expansion");
      
      IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

      IO::log << IO::log_offset << "original internal rotation factor fourier expansion size = " 
	      << _irf_four.size() << "\n";

      typedef std::map<int, double>::iterator pit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;
      
      std::multimap<double, int>::const_iterator otit;

      for(pit_t pit = _irf_four.begin(); pit != _irf_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _mass_index(pit->first);
	  
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(pv[i])
	    //
	    itemp *= 2;
	  
	dtemp = pit->second * pit->second / (double)itemp;
	
	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, pit->first));

      }// fourier expansion cycle

      if(threshold == 0.) {
	//
	std::cerr << funame << "no external rotation\n";

	throw Error::Logic();
      }

      threshold *= _mtol;
      
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;

	if(dtemp > threshold)
	  //
	  break;

	_irf_four.erase(otit->second);
      }

      IO::log << IO::log_offset << "pruned   internal rotation factor fourier expansion size = " 
	      << _irf_four.size() << "\n";
      //
    }// internal rotation factor pruning

    // curvlinear transformation factor pruning
    //
    if(_mtol > 0.) {
      //
      IO::Marker  prune_marker("pruning curvlinear transformation factor fourier expansion");
      
      IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

      IO::log << IO::log_offset << "original curvlinear transformation factor fourier expansion size = " 
	      << _ctf_four.size() << "\n";

      typedef std::map<int, double>::iterator pit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;
      
      std::multimap<double, int>::const_iterator otit;

      for(pit_t pit = _ctf_four.begin(); pit != _ctf_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _mass_index(pit->first);
	  
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(pv[i])
	    //
	    itemp *= 2;
	  
	dtemp = pit->second * pit->second / (double)itemp;
	
	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, pit->first));

      }// fourier expansion cycle

      if(threshold == 0.) {
	//
	std::cerr << funame << "no external rotation\n";

	throw Error::Logic();
      }

      threshold *= _mtol;
      
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;

	if(dtemp > threshold)
	  //
	  break;

	_ctf_four.erase(otit->second);
      }

      IO::log << IO::log_offset << "pruned   curvlinear transformation factor fourier expansion size = " 
	      << _ctf_four.size() << "\n";
      //
    }// curvlinear transformation factor pruning

    // mass-related complex fourier expansions for rovibrational states quantization
    //
    if(_amom_max) {
      //
      IO::Marker complex_mass_fourier_marker("mass-related complex fourier expansions for rovibrational states quantization");

      vtemp.resize(_mass_index.size());

      Lapack::ComplexVector ctemp;

      typedef std::map<int, Lapack::ComplexMatrix>::iterator pit_t;

      for(int g = 0; g < _mass_index.size(); ++g)
	//
	if(_mass_index.conjugate(g) >= g) {
	  //
	  _internal_mobility_fourier[g].resize(internal_size());

	  _coriolis_coupling_fourier[g].resize(3, internal_size());

	  _external_mobility_fourier[g].resize(3);
	}

      // internal rotation fourier expansion initialization
      //
      for(int i = 0; i < internal_size(); ++i) {
	//
	for(int j = i; j < internal_size(); ++j) {
	  //
	  for(int l = 0; l < _mass_index.size(); ++l)
	    //
	    vtemp[l] = _internal_mobility_real[l](i, j);

	  ctemp = Lapack::fourier_transform(vtemp, _mass_index);

	  for(pit_t pit = _internal_mobility_fourier.begin(); pit != _internal_mobility_fourier.end(); ++pit) {
	    //
	    pit->second(i, j) = ctemp[pit->first];

	    if(i != j)
	      //
	      pit->second(j, i) = ctemp[pit->first];
	  }
	}
      }

      // external rotation fourier expansion initialization
      //
      for(int i = 0; i < 3; ++i) {
	//
	for(int j = i; j < 3; ++j) {
	  //
	  for(int l = 0; l < _mass_index.size(); ++l)
	    //
	    vtemp[l] = _external_mobility_real[l](i, j);

	  ctemp = Lapack::fourier_transform(vtemp, _mass_index);

	  for(pit_t pit = _external_mobility_fourier.begin(); pit != _external_mobility_fourier.end(); ++pit) {
	    //
	    pit->second(i, j) = ctemp[pit->first];

	    if(i != j)
	      //
	      pit->second(j, i) = ctemp[pit->first];
	  }
	}
      }

      // coriolis coupling fourier expansion initialization
      //
      for(int i = 0; i < 3; ++i) {
	//
	for(int j = 0; j < internal_size(); ++j) {
	  //
	  for(int l = 0; l < _mass_index.size(); ++l)
	    //
	    vtemp[l] = _coriolis_coupling_real[l](i, j);

	  ctemp = Lapack::fourier_transform(vtemp, _mass_index);

	  for(pit_t pit = _coriolis_coupling_fourier.begin(); pit != _coriolis_coupling_fourier.end(); ++pit)
	    //
	    pit->second(i, j) = ctemp[pit->first];
	}
      }
      
      // curvlinear transformation factor initialization
      //
      if(_with_ctf) {
	//
	ctemp = Lapack::fourier_transform(_ctf_real, _mass_index);

	for(int g = 0; g < _mass_index.size(); ++g)
	  //
	  if(_mass_index.conjugate(g) >= g)
	    //
	    _ctf_complex_fourier[g] = ctemp[g];
      }

      // mass-related complex fourier expansions pruning
      //
      if(_mtol > 0.) {
	//
	IO::Marker  prune_marker("pruning mass-related complex fourier expansions");
	
	double threshold;

	std::multimap<double, int> order_term;
	
	std::multimap<double, int>::const_iterator otit;

	IO::log << IO::log_offset << "tolerance[%] = " << _mtol * 100. << "\n";

	// internal mobility
	//
	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "original effective internal mobility expansion size = " 
		  << _internal_mobility_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "original internal mobility expansion size = " 
		  << _internal_mobility_fourier.size() << "\n";

	threshold = 0.;

	order_term.clear();

	for(pit_t pit = _internal_mobility_fourier.begin(); pit != _internal_mobility_fourier.end(); ++pit) {// fouricer expansion cycle
	  //
	  dtemp = 0.;
	  
	  for(int i = 0; i < internal_size(); ++i) {
	    //
	    for(int j = i; j < internal_size(); ++j) {
	      //
	      if(i != j) {
		//
		dtemp += 2. * std::norm(pit->second(i, j));
	      }
	      else 
		//
		dtemp += std::norm(pit->second(i, j));
	    }//
	    //
	  }//

	  threshold += dtemp;

	  order_term.insert(std::make_pair(dtemp, pit->first));
			    
	}// fourier expansion cycle

	if(threshold == 0.) {
	  //
	  std::cerr << funame << "no internal mobility\n";

	  throw Error::Logic();
	}

	threshold *= _mtol;
	
	dtemp = 0.;
      
	for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	  //
	  dtemp += otit->first;
	  
	  if(dtemp > threshold)
	    //
	    break;

	  _internal_mobility_fourier.erase(otit->second);
	}

	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "pruned   effective internal mobility expansion size = " 
		  << _internal_mobility_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "pruned   internal mobility expansion size = " 
		  << _internal_mobility_fourier.size() << "\n";

	// external mobility
	//
	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "original effective external mobility expansion size = " 
		  << _external_mobility_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "original external mobility expansion size = " 
		  << _external_mobility_fourier.size() << "\n";

	threshold = 0.;

	order_term.clear();

	for(pit_t pit = _external_mobility_fourier.begin(); pit != _external_mobility_fourier.end(); ++pit) {// fourier expansion cycle
	  //
	  dtemp = 0.;
	  
	  for(int i = 0; i < 3; ++i) {
	    //
	    for(int j = i; j < 3; ++j) {
	      //
	      if(i != j) {
		//
		dtemp += 2. * std::norm(pit->second(i, j));
	      }
	      else 
		//
		dtemp += std::norm(pit->second(i, j));
	      //
	    }//
	    //
	  }//

	  threshold += dtemp;

	  order_term.insert(std::make_pair(dtemp, pit->first));

	}// fourier expansion cycle

	if(threshold == 0.) {
	  //
	  std::cerr << funame << "no external mobility\n";
	  
	  throw Error::Logic();
	}

	threshold *= _mtol;
	
	dtemp = 0.;
      
	for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	  //
	  dtemp += otit->first;
	  
	  if(dtemp > threshold)
	    //
	    break;

	  _external_mobility_fourier.erase(otit->second);
	}

	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "pruned   effective external mobility expansion size = " 
		  << _external_mobility_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "pruned   external mobility expansion size = " 
		  << _external_mobility_fourier.size() << "\n";

	// coriolis coupling
	//
	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "original effective coriolis coupling expansion size = " 
		  << _coriolis_coupling_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "original coriolis coupling expansion size = " 
		  << _coriolis_coupling_fourier.size() << "\n";

	threshold = 0.;

	order_term.clear();

	for(pit_t pit = _coriolis_coupling_fourier.begin(); pit != _coriolis_coupling_fourier.end(); ++pit) {// fourier expansion cycle
	  //
	  dtemp = 0.;
	  
	  for(int i = 0; i < 3; ++i)
	    //
	    for(int j = 0; j < internal_size(); ++j)
	      //
	      dtemp += std::norm(pit->second(i, j));

	  threshold += dtemp;
	  
	  order_term.insert(std::make_pair(dtemp, pit->first));
	  
	}// fourier expansion cycle
	
	threshold *= _mtol;
	
	dtemp = 0.;
      
	for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	  //
	  dtemp += otit->first;
	  
	  if(dtemp > threshold)
	    //
	    break;

	  _coriolis_coupling_fourier.erase(otit->second);
	}
	
	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "pruned   effective coriolis coupling expansion size = " 
		  << _coriolis_coupling_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "pruned   coriolis coupling expansion size = " 
		  << _coriolis_coupling_fourier.size() << "\n";

	// curvlinear transformation factor pruning
	//
	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "original curvlinear transformation factor expansion size = " 
		  << _ctf_complex_fourier.size() << "\n";

	  threshold = 0.;

	  order_term.clear();

	  typedef std::map<int, Lapack::complex>::const_iterator pit_t;

	  for(pit_t pit = _ctf_complex_fourier.begin(); pit != _ctf_complex_fourier.end(); ++pit) {// fourier expansion cycle
	    //
	    dtemp = std::norm(pit->second);

	    threshold += dtemp;
	  
	    order_term.insert(std::make_pair(dtemp, pit->first));
	    //
	  }// fourier expansion cycle
	
	  threshold *= _mtol;
	
	  dtemp = 0.;
      
	  for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	    //
	    dtemp += otit->first;
	  
	    if(dtemp > threshold)
	      //
	      break;

	    _ctf_complex_fourier.erase(otit->second);
	  }
	  
	  IO::log << IO::log_offset << "pruned   curvlinear transformation factor expansion size = " 
		  << _coriolis_coupling_fourier.size() << "\n";
	  //
	}//
	//
      }// mass-related complex fourier expansions pruning

      // add conjugated points for internal mobility 
      //
      itemp = 0;

      ivec.resize(_internal_mobility_fourier.size());

      for(pit_t pit = _internal_mobility_fourier.begin(); pit != _internal_mobility_fourier.end(); ++pit, ++itemp)
	//
	ivec[itemp] = pit->first;
	
      for(std::vector<int>::const_iterator pit = ivec.begin(); pit != ivec.end(); ++pit) { 
	//
	itemp = _mass_index.conjugate(*pit);

	if(itemp != *pit) {
	  //
	  _internal_mobility_fourier[itemp].resize(internal_size());

	  for(int i = 0; i < internal_size(); ++i)
	    //
	    for(int j = 0; j < internal_size(); ++j)		
	      //
	      _internal_mobility_fourier[itemp](i, j) = std::conj(_internal_mobility_fourier[*pit](i, j));
	}
      }

      // add conjugated points for external mobility
      //
      itemp = 0;

      ivec.resize(_external_mobility_fourier.size());
      
      for(pit_t pit = _external_mobility_fourier.begin(); pit != _external_mobility_fourier.end(); ++pit, ++itemp)
	//
	ivec[itemp] = pit->first;
	
      for(std::vector<int>::const_iterator pit = ivec.begin(); pit != ivec.end(); ++pit) { 
	//
	itemp = _mass_index.conjugate(*pit);

	if(itemp != *pit) {
	  //
	  _external_mobility_fourier[itemp].resize(3);

	  for(int i = 0; i < 3; ++i)
	    //
	    for(int j = 0; j < 3; ++j)		
	      //
	      _external_mobility_fourier[itemp](i, j) = std::conj(_external_mobility_fourier[*pit](i, j));
	}
      }

      // add conjugated points for coriolis coupling
      //
      itemp = 0;

      ivec.resize(_coriolis_coupling_fourier.size());

      for(pit_t pit = _coriolis_coupling_fourier.begin(); pit != _coriolis_coupling_fourier.end(); ++pit, ++itemp)
	//
	ivec[itemp] = pit->first;
	
      for(std::vector<int>::const_iterator pit = ivec.begin(); pit != ivec.end(); ++pit) { 
	//
	itemp = _mass_index.conjugate(*pit);

	if(itemp != *pit) {
	  //
	  _coriolis_coupling_fourier[itemp].resize(3, internal_size());

	  for(int i = 0; i < 3; ++i)
	    //
	    for(int j = 0; j < internal_size(); ++j)		
	      //
	      _coriolis_coupling_fourier[itemp](i, j) = std::conj(_coriolis_coupling_fourier[*pit](i, j));
	}
      }
    
      // add conjugated points for curvlinear transformation factor
      //
      if(_with_ctf) {
	//
	typedef std::map<int, Lapack::complex>::const_iterator pit_t;
	

	itemp = 0;

	ivec.resize(_ctf_complex_fourier.size());

	for(pit_t  pit = _ctf_complex_fourier.begin(); pit != _ctf_complex_fourier.end(); ++pit, ++itemp)
	  //
	  ivec[itemp] = pit->first;
	
	for(std::vector<int>::const_iterator pit = ivec.begin(); pit != ivec.end(); ++pit) { 
	  //
	  itemp = _mass_index.conjugate(*pit);

	  if(itemp != *pit) {
	    //
	    _ctf_complex_fourier[itemp] =std::conj(_ctf_complex_fourier[*pit]);
	  }
	}
      }

      // checking
      /*
	std::map<int, Lapack::complex> cval;
	IO::log << IO::log_offset << "internal mobility test:\n";
	for(int i = 0; i < internal_size(); ++i)
	for(int j = 0; j < internal_size(); ++j) {		
	for(pit_t pit = _internal_mobility_fourier.begin(); pit != _internal_mobility_fourier.end();	++pit)
	cval[pit->first] = pit->second(i, j);
	ctemp = Lapack::fourier_transform(cval, _mass_index);
	double dif_val = 0.;
	double ref_val = 0.;
	for(int g = 0; g < _mass_index.size(); ++g) {
	dtemp = _internal_mobility_real[g](i, j);
	ref_val += dtemp * dtemp;
	dtemp -= std::real(ctemp[g]);
	dif_val += dtemp * dtemp;
	}
	if(ref_val != 0.)
	IO::log << IO::log_offset << "i = " << i << " j = " << j
	<< std::setw(13) << dif_val / ref_val << "\n";
	}
	IO::log << IO::log_offset << "external mobility test:\n";
	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j) {		
	for(pit_t pit = _external_mobility_fourier.begin(); pit != _external_mobility_fourier.end();	++pit)
	cval[pit->first] = pit->second(i, j);
	ctemp = Lapack::fourier_transform(cval, _mass_index);
	double dif_val = 0.;
	double ref_val = 0.;
	for(int g = 0; g < _mass_index.size(); ++g) {
	dtemp = _external_mobility_real[g](i, j);
	ref_val += dtemp * dtemp;
	dtemp -= std::real(ctemp[g]);
	dif_val += dtemp * dtemp;
	}
	if(ref_val != 0.)
	IO::log << IO::log_offset << "i = " << i << " j = " << j
	<< std::setw(13) << dif_val / ref_val << "\n";
	}
	IO::log << IO::log_offset << "coriolis coupling test:\n";
	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < internal_size(); ++j) {		
	for(pit_t pit = _coriolis_coupling_fourier.begin(); pit != _coriolis_coupling_fourier.end(); ++pit)
	cval[pit->first] = pit->second(i, j);
	ctemp = Lapack::fourier_transform(cval, _mass_index);
	double dif_val = 0.;
	double ref_val = 0.;
	for(int g = 0; g < _mass_index.size(); ++g) {
	dtemp = _coriolis_coupling_real[g](i, j);
	ref_val += dtemp * dtemp;
	dtemp -= std::real(ctemp[g]);
	dif_val += dtemp * dtemp;
	}
	if(ref_val != 0.)
	IO::log << IO::log_offset << "i = " << i << " j = " << j
	<< std::setw(13) << dif_val / ref_val << "\n";
	}
      */

      // output
      //
      if(_with_ctf) {
	//
	IO::log << IO::log_offset << "effective internal mobility expansion size (including conjugated terms) = " 
		<< _internal_mobility_fourier.size() << "\n";

	IO::log << IO::log_offset << "effective external mobility expansion size (including conjugated terms) = " 
		<< _external_mobility_fourier.size() << "\n";
	
	IO::log << IO::log_offset << "effective coriolis coupling expansion size (including conjugated terms) = " 
		<< _coriolis_coupling_fourier.size() << "\n";

	IO::log << IO::log_offset << "curvlinear factor expansion size (including conjugated terms) = " 
		<< _ctf_complex_fourier.size() << "\n";
      }
      else {
	//
	IO::log << IO::log_offset << "internal mobility expansion size (including conjugated terms) = " 
		<< _internal_mobility_fourier.size() << "\n";

	IO::log << IO::log_offset << "external mobility expansion size (including conjugated terms) = " 
		<< _external_mobility_fourier.size() << "\n";
	
	IO::log << IO::log_offset << "coriolis coupling expansion size (including conjugated terms) = " 
		<< _coriolis_coupling_fourier.size() << "\n";
      }
      //
    }// mass-related complex fourier expansions for rovibrational states quantization
    //
  }// mass fourier expansion
  
  /********************************* POTENTIAL FOURIER EXPANSION *************************************/

  if(1) {
    //
    IO::Marker work_marker("potential fourier expansion");

    // potential fourier expansion dimensions
    //
    ivec.resize(internal_size());

    for(int r = 0; r < internal_size(); ++r) {
      //
      itemp = _pot_index.size(r);

      itemp = itemp % 2 ? itemp : itemp + 1;

      if(_internal_rotation[r].potential_fourier_size()) {
	//
	if(_internal_rotation[r].potential_fourier_size() > itemp) {
	  //
	  IO::log << IO::log_offset << "WARNING: requested potential fourier expansion size should not be bigger than sampling size\n";
	}
	else 
	  //
	  itemp =  _internal_rotation[r].potential_fourier_size();
      }

      ivec[r] = itemp; 
    }

    _pot_four_index.resize(ivec);

    IO::log << IO::log_offset << "potential fourier expansion size = " << _pot_four_index.size() << std::endl;

    // potential fourier expansion
    //
    for(int g = 0; g < _pot_index.size(); ++g) {// grid cycle
      //
      std::vector<int> gv = _pot_index(g);

      if(_with_ctf) {
	//
	// sampling configuration
	//
	std::vector<double> angle(internal_size());

	for(int r = 0; r < internal_size(); ++r)
	  //
	  angle[r] = 2. * M_PI * double(gv[r]) / double(symmetry(r) * _pot_index.size(r));

	dtemp = curvlinear_factor(angle);
      }

      const double ctf = dtemp;
      
      for(int h = 0; h < _pot_four_index.size(); ++h) {// fourier expansion cycle
	//
	std::vector<int> hv = _pot_four_index(h);

	// fourier transform
	//
	double dfac = 1.;

	for(int r = 0; r < internal_size(); ++r) {
	  //
	  if(!hv[r])
	    //
	    continue;

	  if(hv[r] % 2) {
	    //
	    dfac *= std::sin(M_PI * double((hv[r] + 1)  * gv[r]) / double(_pot_index.size(r)));
	  }
	  else
	    //
	    dfac *= std::cos(M_PI * double(hv[r] * gv[r]) / double(_pot_index.size(r)));
	}

	// potential
	//
	_pot_four[h] += dfac * _pot_real[g];

	// vibrational frequencies
	//
	for(int v = 0; v < _vib_four.size(); ++v) {
	  //
	  _vib_four[v][h] += dfac * vibration_sampling[g][v];
	}//

	if(_with_ctf) {
	  //
	  _eff_pot_four[h] += dfac * _pot_real[g] * ctf;

	  // vibrational frequencies
	  //
	  for(int v = 0; v < _vib_four.size(); ++v) {
	    //
	    _eff_vib_four[v][h] += dfac * vibration_sampling[g][v] * ctf;
	  }//
	  //
	}//
	//
      }// fourier expansion cycle
      //
    }// grid cycle

    // normalization
    //
    for(int h = 0; h < _pot_four_index.size(); ++h) {// fourier expansion cycle
      //
      std::vector<int> hv = _pot_four_index(h);

      double dfac = double(_pot_index.size());

      for(int r = 0; r < internal_size(); ++r)
	//
	if(hv[r] && hv[r] != _pot_index.size(r))
	  //
	  dfac /= 2.;

      _pot_four[h] /= dfac;	  

      for(int v = 0; v < _vib_four.size(); ++v)
	//
	_vib_four[v][h] /= dfac;

      if(_with_ctf) {
	//
	_eff_pot_four[h] /= dfac;	  

	for(int v = 0; v < _vib_four.size(); ++v)
	  //
	  _eff_vib_four[v][h] /= dfac;
      }//
      //
    }// fourier expansion cycle 

    _pot_shift = _pot_four.begin()->second;

    _pot_four.erase(_pot_four.begin());
    
    //IO::log << IO::log_offset << "potential constant part[kcal/mol] = " << _pot_shift / Phys_const::kcal << "\n";
    
    // potential correlation analysis
    //
    if(_pot_four.size()) {
      //
      typedef std::map<int, double>::const_iterator pit_t;

      // different order correlation terms contributions
      //
      dvec.resize(internal_size());

      for(int i = 0;  i < dvec.size(); ++i)
	//
	dvec[i] = 0.;

      
      for(pit_t pit = _pot_four.begin(); pit != _pot_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _pot_four_index(pit->first);

	itemp = 0;

	for(int r = 0; r < internal_size(); ++r)
	  //
	  if(pv[r])
	    //
	    ++itemp;

	dvec[itemp - 1] += pit->second * pit->second / std::pow(2., (double)itemp);
	//
      }// fourier expansion cycle

      dtemp = 0.;

      for(int i = 0; i < dvec.size(); ++i)
	//
	dtemp += dvec[i];

      if(dtemp == 0.) {
	//
	IO::log << IO::log_offset << "WARNING: potential is a constant\n";
      }
      else {
	//
	IO::log << IO::log_offset 
		<< "different order potential correlation terms contributions, %:\n"
		<< IO::log_offset << std::setw(3) << "#";

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(13) << r + 1; 

	IO::log  << "\n";

	IO::log << IO::log_offset << std::setw(3) << "%"; 

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(13) << 100. * dvec[r] / dtemp;

	IO::log << "\n";
      }
    
      // first order correlation contributions from different rotations
      //
      dvec.resize(internal_size());

      for(int i = 0;  i < dvec.size(); ++i)
	//
	dvec[i] = 0.;

      for(pit_t pit = _pot_four.begin(); pit != _pot_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _pot_four_index(pit->first);

	itemp = 0;

	int rr;

	for(int r = 0; r < internal_size(); ++r) {
	  //
	  if(pv[r]) {
	    //
	    rr = r;

	    ++itemp;
	  }
	}

	if(itemp == 1)
	  //
	  dvec[rr] += pit->second * pit->second / 2.;
	//
      }// fourier expansion cycle

      dtemp = 0.;

      for(int i = 0; i < dvec.size(); ++i)
	//
	dtemp += dvec[i];

      if(dtemp > 0.) {
	//
	IO::log << IO::log_offset 
		<< "first order correlation contributions to the potential from different rotations, %:\n"
		<< IO::log_offset << std::setw(3) << "R";

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(13) << r; 

	IO::log  << "\n";
	
	IO::log << IO::log_offset << std::setw(3) << "%"; 

	for(int r = 0; r < internal_size(); ++r)
	  //
	  IO::log << std::setw(13) << 100. * dvec[r] / dtemp;

	IO::log << "\n";
	//
      }//

      // high order harmonics potential contributions
      //
      IO::log << IO::log_offset << "high order harmonics potential contributions, %:\n";

      IO::log << std::setprecision(3);

      itemp = 0;

      for(int r = 0; r < internal_size(); ++r)
	//
	if(_pot_four_index.size(r) > itemp)
	  //
	  itemp = _pot_four_index.size(r);

      itemp /= 2;

      IO::log << IO::log_offset << std::setw(3) << "R\\#";
    
      for(int i = 0; i < itemp; ++i)
	//
	IO::log << std::setw(9) << i + 1;

      IO::log << "\n";
    
      for(int r = 0;  r < internal_size(); ++r) {
	//
	dvec.resize(_pot_four_index.size(r) / 2 + 1);

	for(int i = 0;  i < dvec.size(); ++i)
	  //
	  dvec[i] = 0.;

	for(pit_t pit = _pot_four.begin(); pit != _pot_four.end(); ++pit) {// fourier expansion cycle
	  //
	  std::vector<int> pv = _pot_four_index(pit->first);

	  itemp = 1;

	  for(int i = 0; i < internal_size(); ++i)
	    //
	    if(pv[i])
	      //
	      itemp *= 2;

	  dvec[(pv[r] + 1) / 2] += pit->second * pit->second / (double)itemp;
	  //
	}//fourier expansion cycle

	dtemp = 0.;

	for(int i = 1; i < dvec.size(); ++i)
	  //
	  dtemp += dvec[i];

	IO::log << IO::log_offset << std::setw(3) << r;

	if(dtemp > 0.)
	  //
	  for(int i = 1; i < dvec.size(); ++i)
	    //
	    IO::log << std::setw(9) << 100. * dvec[i] / dtemp;

	IO::log << "\n";
      }

      IO::log << std::setprecision(6);
      //
    }// potential correlation analysis
    
    // vibrational frequencies correlation analysis
    //
    if(_pot_four_index.size() > 1 && _vib_four.size()) {
      //
      typedef std::map<int, double>::const_iterator pit_t;

      dvec.resize(internal_size());

      for(int v = 0; v < _vib_four.size(); ++v) {// vibrational frequencies cycle
	//
	for(int i = 0;  i < dvec.size(); ++i)
	  //
	  dvec[i] = 0.;

	for(pit_t pit = _vib_four[v].begin(); pit != _vib_four[v].end(); ++pit) {//fourier expansion cycle
	  //
	  if(pit != _vib_four[v].begin()) {
	    //
	    std::vector<int> pv = _pot_four_index(pit->first);

	    itemp = 0;

	    for(int r = 0; r < internal_size(); ++r)
	      //
	      if(pv[r])
		//
		++itemp;

	    dvec[itemp - 1] += pit->second * pit->second / std::pow(2., (double)itemp);
	  }
	}

	dtemp = 0.;

	for(int i = 0; i < dvec.size(); ++i)
	  //
	  dtemp += dvec[i];

	if(dtemp <= 0.) {
	  //
	  IO::log << IO::log_offset << "WARNING: " << v << "-th vibrational frequency is a constant\n";
	}
	else {
	  //
	  IO::log << IO::log_offset << std::setw(3) << std::left << v + 1 << std::right; 
	  
	  for(int r = 0; r < internal_size(); ++r)
	    //
	    IO::log << std::setw(13) << 100. * dvec[r] / dtemp;

	  IO::log << "\n";
	  //
	}//
	//
      }// vibrational frequencies cycle
      //
    }// vibrational frequencies correlation analysis

    std::vector<double> angle(internal_size());

    // output: potential check
    //
    IO::log << IO::log_offset << "Model potential:\n";
    
    IO::log << IO::log_offset;

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << std::setw(14) << "phi_" << r + 1;

    IO::log << std::setw(15) << "V, kcal/mol";

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << std::setw(14) << "f_" << r + 1;

    IO::log << "\n";
      
    for(int g = 0; g < _pot_index.size(); ++g) {
      //
      std::vector<int> gv = _pot_index(g);

      IO::log << IO::log_offset;

      for(int r = 0; r < internal_size(); ++r) {
	//
	angle[r] = 2. * M_PI * double(gv[r]) / double(symmetry(r) * _pot_index.size(r));

	IO::log << std::setw(15) << angle[r] * 180. / M_PI;
      }

      IO::log << std::setw(15) << potential(angle) / Phys_const::kcal;

      vtemp = frequencies(angle);

      for(int r = 0; r < internal_size(); ++r)
	//
	IO::log << std::setw(15) << vtemp[r] / Phys_const::incm;

      IO::log << "\n";
    }

    IO::log << "\n";

    // potential pruning
    //
    if(_ptol > 0.) {
      //
      IO::Marker  prune_marker("pruning potential fourier expansion");
      
      IO::log << IO::log_offset << "tolerance[%] = " << _ptol << "\n";

      IO::log << IO::log_offset << "original potential fourier expansion size = " 
	      << _pot_four.size() << "\n";

      typedef std::map<int, double>::iterator pit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;

      std::multimap<double, int>::const_iterator otit;

      for(pit_t pit = _pot_four.begin(); pit != _pot_four.end(); ++pit) {// fourier expansion cycle
	//
	std::vector<int> pv = _pot_four_index(pit->first);
	    
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(pv[i])
	    //
	    itemp *= 2;
	    
	dtemp = pit->second * pit->second / (double)itemp;

	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, pit->first));

      }// fourier expansion cycle

      threshold *= _ptol;
	
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;
	
	if(dtemp > threshold)
	  //
	  break;

	_pot_four.erase(otit->second);
      }
      
      IO::log << IO::log_offset << "pruned   potential fourier expansion size = " 
	      << _pot_four.size() << "\n";
      //
    }// potential pruning

    // effective potential pruning
    //
    if(_ptol > 0. && _with_ctf) {
      //
      IO::Marker  prune_marker("pruning effective potential fourier expansion");
      
      IO::log << IO::log_offset << "tolerance[%] = " << _ptol << "\n";

      IO::log << IO::log_offset << "original effective potential fourier expansion size = " 
	      << _eff_pot_four.size() << "\n";

      typedef std::map<int, double>::iterator pit_t;

      double threshold = 0.;

      std::multimap<double, int> order_term;

      std::multimap<double, int>::const_iterator otit;

      for(pit_t pit = _eff_pot_four.begin(); pit != _eff_pot_four.end(); ++pit) {// fourier expansion cycle
	//
	if(pit == _eff_pot_four.begin())
	  //
	  continue;

	std::vector<int> pv = _pot_four_index(pit->first);
	    
	itemp = 1;

	for(int i = 0; i < internal_size(); ++i)
	  //
	  if(pv[i])
	    //
	    itemp *= 2;
	    
	dtemp = pit->second * pit->second / (double)itemp;

	threshold += dtemp;

	order_term.insert(std::make_pair(dtemp, pit->first));

      }// fourier expansion cycle

      threshold *= _ptol;
	
      dtemp = 0.;
      
      for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	//
	dtemp += otit->first;
	
	if(dtemp > threshold)
	  //
	  break;

	_eff_pot_four.erase(otit->second);
      }
      
      IO::log << IO::log_offset << "pruned   effective potential fourier expansion size = " 
	      << _eff_pot_four.size() << "\n";
      //
    }// effective potential pruning

    // vibrational frequencies pruning
    //
    if(_vtol > 0. && _vib_four.size()) {
      //
      IO::Marker  prune_marker("pruning vibrational frequencies fourier expansion");

      IO::log << IO::log_offset << "tolerance[%] = " << _vtol << "\n";

      typedef std::map<int, double>::iterator pit_t;
      
      for(int v = 0; v < _vib_four.size(); ++v) {// vibrational frequencies cycle
	//
	IO::log << IO::log_offset << "original " << std::setw(2) << v << "-th vibrational frequency fourier expansion size = "
		<< _vib_four[v].size() << "\n";

	double threshold = 0.;

	std::multimap<double, int> order_term;

	std::multimap<double, int>::const_iterator otit;

	for(pit_t pit = _vib_four[v].begin(); pit != _vib_four[v].end(); ++pit) {// fourier expansion cycle
	  //
	  std::vector<int> pv = _pot_four_index(pit->first);
	    
	  itemp = 1;

	  for(int i = 0; i < internal_size(); ++i)
	    //
	    if(pv[i])
	      //
	      itemp *= 2;

	  dtemp = pit->second * pit->second / (double)itemp;;

	  threshold += dtemp;

	  order_term.insert(std::make_pair(dtemp, pit->first));

	}// fourier expansion cycle

	threshold *= _vtol;
	
	dtemp = 0.;
      
	for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	  //
	  dtemp += otit->first;
	  
	  if(dtemp > threshold)
	    //
	    break;

	  _vib_four[v].erase(otit->second);
	}

	IO::log << IO::log_offset << "pruned   " << std::setw(2) << v << "-th vibrational frequency fourier expansion size = "
		<< _vib_four[v].size() << "\n";
	//
      }// vibrational frequencies cycle
      //
    }// vibrational frequencies pruning

    // effective vibrational frequencies pruning
    //
    if(_vtol > 0. && _eff_vib_four.size() && _with_ctf) {
      //
      IO::Marker  prune_marker("pruning effective vibrational frequencies fourier expansion");

      IO::log << IO::log_offset << "tolerance[%] = " << _vtol << "\n";

      typedef std::map<int, double>::iterator pit_t;
      
      for(int v = 0; v < _eff_vib_four.size(); ++v) {// vibrational frequencies cycle
	//
	IO::log << IO::log_offset << "original " << std::setw(2) << v << "-th effective vibrational frequency fourier expansion size = "
		<< _eff_vib_four[v].size() << "\n";

	double threshold = 0.;

	std::multimap<double, int> order_term;

	std::multimap<double, int>::const_iterator otit;

	for(pit_t pit = _eff_vib_four[v].begin(); pit != _eff_vib_four[v].end(); ++pit) {// fourier expansion cycle
	  //
	  std::vector<int> pv = _pot_four_index(pit->first);
	    
	  itemp = 1;

	  for(int i = 0; i < internal_size(); ++i)
	    //
	    if(pv[i])
	      //
	      itemp *= 2;

	  dtemp = pit->second * pit->second / (double)itemp;;

	  threshold += dtemp;

	  order_term.insert(std::make_pair(dtemp, pit->first));

	}// fourier expansion cycle

	threshold *= _vtol;
	
	dtemp = 0.;
      
	for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	  //
	  dtemp += otit->first;
	  
	  if(dtemp > threshold)
	    //
	    break;

	  _eff_vib_four[v].erase(otit->second);
	}

	IO::log << IO::log_offset << "pruned   " << std::setw(2) << v << "-th effective vibrational frequency fourier expansion size = "
		<< _eff_vib_four[v].size() << "\n";
	//
      }// vibrational frequencies cycle
      //
    }// effective vibrational frequencies pruning

    // effective potential complex fourier expansion and pruning for rovibrational states quantization
    //
    if(_amom_max) {
      //
      IO::Marker complex_pot_four_marker("potential complex fourier expansion for rovibrational states quantization");

      vtemp.resize(_pot_index.size());

      for(int g = 0; g < _pot_index.size(); ++g) {
	//
	std::vector<int> gv = _pot_index(g);

	if(_with_ctf) {
	  //
	  // sampling configuration
	  //
	  std::vector<double> angle(internal_size());

	  for(int r = 0; r < internal_size(); ++r)
	    //
	    angle[r] = 2. * M_PI * double(gv[r]) / double(symmetry(r) * _pot_index.size(r));

	  vtemp[g] = _pot_real[g] * curvlinear_factor(angle);
	}
	else
	  //
	  vtemp[g] = _pot_real[g];
      }

      Lapack::ComplexVector ctemp = Lapack::fourier_transform(vtemp, _pot_index);

      typedef std::map<int, Lapack::complex>::iterator pit_t;
	
      // initialization
      //
      for(int g = 1; g < _pot_index.size(); ++g) {
	//
	std::vector<int> gv = _pot_index(g);

	btemp = true;
	
	for(int r = 0; r < internal_size(); ++r) {
	  //
	  itemp = _internal_rotation[r].potential_fourier_size() ? _internal_rotation[r].potential_fourier_size() / 2 : -1;

	  if(itemp >= 0 && itemp < gv[r] && itemp < _pot_index.size(r) - gv[r]) {
	    //
	    btemp = false;

	    break;
	  }
	}

	if(btemp && _pot_index.conjugate(g) >= g)
	  //
	  _pot_complex_fourier[g] = ctemp[g];
	//
      }// initialization

      // potential complex fourier expansion pruning
      //
      if(_ptol > 0.) {
	//
	IO::Marker  prune_marker("pruning potential complex fourier expansion");
	
	IO::log << IO::log_offset << "tolerance[%] = " << _ptol * 100. << "\n";

	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "original effective potential fourier expansion size = " << _pot_complex_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "original potential fourier expansion size = " << _pot_complex_fourier.size() << "\n";

	double threshold = 0.;

	std::multimap<double, int> order_term;

	std::multimap<double, int>::const_iterator otit;

	for(pit_t pit = _pot_complex_fourier.begin(); pit != _pot_complex_fourier.end(); ++pit) {// fourier expansion cycle
	  //
	  dtemp = std::norm(pit->second);

	  threshold += dtemp;

	  order_term.insert(std::make_pair(dtemp, pit->first));
	  
	}// fourier expansion cycle

	threshold *= _ptol;
	
	dtemp = 0.;
      
	for(otit = order_term.begin(); otit != order_term.end(); ++otit) {
	  //
	  dtemp += otit->first;
	  
	  if(dtemp > threshold)
	    //
	    break;

	  _pot_complex_fourier.erase(otit->second);
	}

	if(_with_ctf) {
	  //
	  IO::log << IO::log_offset << "pruned   effective potential fourier expansion size = " << _pot_complex_fourier.size() << "\n";
	}
	else
	  //
	  IO::log << IO::log_offset << "pruned   potential fourier expansion size = " << _pot_complex_fourier.size() << "\n";

	//
      }// potential fourier expansion pruning

      _pot_complex_fourier[0] = ctemp[0];
      
      // add conjugated points
      //
      itemp = 0;

      ivec.resize(_pot_complex_fourier.size());

      for(pit_t pit = _pot_complex_fourier.begin(); pit != _pot_complex_fourier.end(); ++pit, ++itemp)
	//
	ivec[itemp] = pit->first;
      
      for(std::vector<int>::const_iterator pit = ivec.begin(); pit != ivec.end(); ++pit) { 
	//
	itemp = _pot_index.conjugate(*pit);

	if(itemp != *pit)
	  //
	  _pot_complex_fourier[itemp] = std::conj(_pot_complex_fourier[*pit]);
      }
      
      // checking
      /*
	ctemp = Lapack::fourier_transform(_pot_complex_fourier, _pot_index);
	IO::log << IO::log_offset << "potential checking\n";
	for(int g = 0; g < _pot_index.size(); ++g)
	IO::log << IO::log_offset 
	<< std::setw(13) << _pot_real[g]
	<< std::setw(13) << std::real(ctemp[g])
	<< "\n";
      */

      if(_with_ctf) {
	//
	IO::log << IO::log_offset << "effective potential complex fourier expansion size (including conjugated terms) = " 
		<< _pot_complex_fourier.size() << "\n";
      }
      else
	//
	IO::log << IO::log_offset << "potential complex fourier expansion size (including conjugated terms) = " 
		<< _pot_complex_fourier.size() << "\n";

      //
    }// potential complex fourier expansion and pruning
    //
  }// potential expansion

  /*************************************** DESCRETIZATION *******************************************/

  if(1) {
    //
    IO::Marker work_marker("descretization");

    double pot_max;

    std::vector<double> angle_min(internal_size()), angle_max(internal_size());

    std::vector<double> vib_min, vib_max;

    // angular integration grid dimensions
    //
    ivec.resize(internal_size());

    for(int r = 0; r < internal_size(); ++r) 
      //
      ivec[r] = _internal_rotation[r].weight_sampling_size();

    _grid_index.resize(ivec);

    IO::log << IO::log_offset << "angular integration grid size = " 
	    << _grid_index.size() << std::endl;

    // angular step and volume
    //
    _angle_grid_step.resize(internal_size());

    _angle_grid_cell = 1.;

    for(int r = 0; r < internal_size(); ++r) { 
      //
      dtemp = 2. * M_PI / double(symmetry(r) * _grid_index.size(r)); 

      _angle_grid_step[r] = dtemp;

      _angle_grid_cell *= dtemp;
    }
  
    _pot_grid.resize(_grid_index.size());

    _freq_grid.resize(_grid_index.size());

    _irf_grid.resize(_grid_index.size());

    if(_with_ext_rot)
      //
      _erf_grid.resize(_grid_index.size());
      
    if(_vib_four.size()) {
      //
      _vib_grid.resize(_grid_index.size());

      vib_min.resize(_vib_four.size());

      vib_max.resize(_vib_four.size());
    }

#pragma omp parallel for default(shared) private(itemp, dtemp) schedule(static)

    // grid cycle
    //
    for(int g = 0; g < _grid_index.size(); ++g) {
      //
      // current angle
      //
      std::vector<int> gv = _grid_index(g);

      std::vector<double> angle(internal_size());

      for(int r = 0; r < internal_size(); ++r)
	//
	angle[r] = double(gv[r]) * _angle_grid_step[r];

      // potential
      //
      _pot_grid[g]  = potential(angle);

      // mass factor
      //
      try {
	//
	_irf_grid[g] =  internal_rotation_factor(angle);
      } 
      catch(Error::General) {
	//
	std::cerr << funame 
		  << "interpolated mass matrix is not positive definite: "
	  "check your internal rotation definitions and/or increase mass expansion sizes\n";
	throw;
      }

      // external rotation mass factor correction
      //
      if(_with_ext_rot) {
	//
	dtemp = external_rotation_factor(angle);

	if(dtemp < 0.) {
	  //
	  std::cerr << funame << "negative external rotation factor at (";

	  for(int r = 0; r < internal_size(); ++r) {
	    //
	    if(r) {
	      //
	      std::cerr << ", ";
	    }
	    std::cerr << angle[r] * 180. / M_PI;
	  }
	  std::cerr << ") grid point\n";

	  throw Error::Range();
	}

	_erf_grid[g] = dtemp;
      }

      // internal rotation local frequencies
      //
      _freq_grid[g] = frequencies(angle);

      // vibrational frequencies
      //
      if(_vib_four.size()) {
	//
	_vib_grid[g] = vibration(angle);
	
	for(int v = 0; v < _vib_four.size(); ++v) {
	  //
	  if(_vib_grid[g][v] <= 0.) {
	    //
	    std::cerr << funame << v + 1 << "-th negative vibrational frequency at (";

	    for(int r = 0; r < internal_size(); ++r) {
	      //
	      if(r) {
		//
		std::cerr << ", ";
	      }
	      std::cerr << angle[r] * 180. / M_PI;
	    }
	    std::cerr << ") grid point\n";

	    throw Error::Range();
	  }
	}
      }
    }// grid cycle

    // potential global minima and maxima
    //
    int gmin, gmax;

    for(int g = 0; g < _grid_index.size(); ++g) {// grid cycle
      //
      dtemp = _pot_grid[g];

      if(!g || dtemp < _pot_global_min) {
	//
	_pot_global_min   = dtemp;

	gmin = g;
      }

      if(!g || dtemp > pot_max) {
	//
	pot_max   = dtemp;

	gmax = g;
      }
      
      // vibrational frequencies minimum and maximum
      //
      for(int v = 0; v < _vib_four.size(); ++v) {
	//
	if(!g || _vib_grid[g][v] < vib_min[v]) {
	  //
	  vib_min[v] = _vib_grid[g][v];
	}
	if(!g || _vib_grid[g][v] > vib_max[v]) {
	  //
	  vib_max[v] = _vib_grid[g][v];
	}
      }
    }// grid cycle

    // potential minimum energy correction
    //
    ivec = _grid_index(gmin);
    //
    for(int r = 0; r < internal_size(); ++r)
      //
      angle_min[r] = double(ivec[r]) * _angle_grid_step[r];

    try {
      //
      vtemp = Lapack::Cholesky(force_constant_matrix(angle_min)).invert(potential_gradient(angle_min));

      for(int i = 0; i < internal_size(); ++i)
	//
	angle_min[i] -= vtemp[i];
      
      dtemp = potential(angle_min) - _pot_global_min;

      _pot_global_min += dtemp;

      IO::log << IO::log_offset << "potential minimum energy correction [kcal/mol] = "
	      << dtemp / Phys_const::kcal << "\n";
    }
    catch(Error::General) {
      //
      IO::log << IO::log_offset << "WARNING: force constant matrix at minimum is not positive definite: "
	"cannot get potential minimum energy correction\n";
    }

    // potential maximum energy correction
    //
    ivec = _grid_index(gmax);

    for(int r = 0; r < internal_size(); ++r)
      //
      angle_max[r] = double(ivec[r]) * _angle_grid_step[r];
    
    try {
      //
      Lapack::SymmetricMatrix fcm = force_constant_matrix(angle_max);

      fcm *= -1.;

      vtemp = Lapack::Cholesky(fcm).invert(potential_gradient(angle_max));

      for(int i = 0; i < internal_size(); ++i)
	//
	angle_max[i] += vtemp[i];
      
      dtemp = potential(angle_max) - pot_max;

      pot_max += dtemp;

      IO::log << IO::log_offset << "potential maximum energy correction [kcal/mol] = "
	      << dtemp / Phys_const::kcal << "\n";
      
    }
    catch(Error::General) {
      //
      IO::log << IO::log_offset  << "WARNING: force constant matrix at maximum is not negative definite: "
	"cannot get potential maximum energy correction\n";
    }

    // output: potential minimum and maximum
    //
    IO::log << IO::log_offset << "potential energy surface, kcal/mol:\n"
	    << IO::log_offset << std::setw(5) << "" << std::setw(13) << "E, kcal/mol";

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << std::setw(12) << "phi_" << r + 1 ;

    IO::log << "\n";

    IO::log << IO::log_offset << std::setw(5) << "min" << std::setw(13) << _pot_global_min / Phys_const::kcal;

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << std::setw(13) << angle_min[r] / M_PI * 180.;

    IO::log << "\n";

    IO::log << IO::log_offset << std::setw(5) << "max" << std::setw(13) << pot_max / Phys_const::kcal;

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << std::setw(13) << angle_max[r] / M_PI * 180.;

    IO::log << "\n";

    if(_vib_four.size()) {
      //
      IO::log << IO::log_offset << "vibrational frequencies, 1/cm:\n"
	      << IO::log_offset 
	      << std::left << std::setw(2) << "#" << std::right
	      << std::setw(13) << "min"
	      << std::setw(13) << "max"
	      << std::setw(13) << "mean"
	      << "\n";

      for(int v = 0; v < _vib_four.size(); ++v)
	//
	IO::log << IO::log_offset
		<< std::left << std::setw(2) << v + 1 << std::right
		<< std::setw(13) << vib_min[v]     / Phys_const::incm
		<< std::setw(13) << vib_max[v]     / Phys_const::incm
		<< std::setw(13) << _vib_four[v][0] / Phys_const::incm
		<< "\n";
    }

    // frequencies at minimum
    //
    Lapack::Vector freq_min = frequencies(angle_min);

    IO::log << IO::log_offset << "hindered rotation frequencies f at minimum, 1/cm:\n";

    IO::log << IO::log_offset 
	    << std::left << std::setw(2) << "#" << std::right
	    << std::setw(12) << "f"
	    << "\n";

    for(int f = 0; f < internal_size(); ++f)
      //
      IO::log << IO::log_offset 
	      << std::left << std::setw(2) << f + 1 << std::right
	      << std::setw(12) << freq_min[f] / Phys_const::incm
	      << "\n";

    // mobility parameter
    //
    _mobility_parameter.resize(internal_size());

    _mobility_min = mass(angle_min).invert();

    _mobility_min /= 2.;

    if(internal_size() == 1) {
      //
      _mobility_parameter[0] = _mobility_min(0, 0);

      if(_mobility_parameter[0] <= 0.) {
	//
	std::cerr << funame << "negative mobility parameter\n";

	throw Error::Range();
      }
    }
    else {
      //
      Lapack::SymmetricMatrix bmat(internal_size() - 1);

      Lapack::Vector          bvec(internal_size() - 1);

      for(int r = 0; r < internal_size(); ++r) {
	//
	for(int i = 0; i < bmat.size(); ++i) {
	  //
	  for(int j = i; j < bmat.size(); ++j) {
	    //
	    int im = i < r ? i : i + 1;

	    int jm = j < r ? j : j + 1;

	    bmat(i, j) = _mobility_min(im, jm);
	  }
	}

	for(int i = 0; i < bvec.size(); ++i) {
	  //
	  int im = i < r ? i : i + 1;

	  bvec[i] = _mobility_min(r, im);
	}

	_mobility_parameter[r] = (_mobility_min(r, r) - bvec * (bmat.invert() * bvec));

	if(_mobility_parameter[r] <= 0.) {
	  //
	  std::cerr << funame << "negative mobility parameter\n";

	  throw Error::Range();
	}
      }
    }

    IO::log << IO::log_offset << "effective one-dimensional internal mobilities[1/cm]:";

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << "   " << _mobility_parameter[r] / Phys_const::incm;

    IO::log << "\n";

    // ground state energy in the harmonic approximation
    //
    _ground = _pot_global_min;

    for(int f = 0; f < internal_size(); ++f)
      //
      _ground += freq_min[f] / 2.;

    IO::log << IO::log_offset << "ground energy in the harmonic oscilator approximation = "
	    << _ground / Phys_const::kcal << " kcal/mol\n";

    /*
      #ifdef DEBUG

      // internal rotation frequencies
      IO::log << IO::log_offset << "Internal frequencies on the sampling grid, 1/cm:\n";
      IO::log << IO::log_offset;
      for(int r = 0; r < internal_size(); ++r)
      IO::log << std::setw(4) << "i_" << r + 1;
      for(int r = 0; r < internal_size(); ++r)
      IO::log << std::setw(14) << "f_" << r + 1;
      IO::log << "\n";

      std::vector<double> angle(internal_size());

      for(int g = 0; g < _pot_index.size(); ++g) {
      std::vector<int> gv = _pot_index(g);
      IO::log << IO::log_offset;
      for(int r = 0; r < internal_size(); ++r) {
      angle[r] = 2. * M_PI * double(gv[r]) / double(symmetry(r) * _pot_index.size(r));
      IO::log << std::setw(5) << gv[r] + 1;
      }
      Lapack::Vector ivf = frequencies(angle);
      for(int r = 0; r < internal_size(); ++r)
      IO::log << std::setw(15) << ivf[r] / Phys_const::incm;
      IO::log << "\n";
      }
      IO::log << "\n";

      #endif
    */

  }// desretization

  /********************************* ONE-DIMENSIONAL ROTORS MODEL *****************************************/

  // one-dimensional rotors model ground energy
  //
  double one_ground;

  std::vector<SharedPointer<HinderedRotor> > one_rotor;

  if(_pot_four_index.size() > 1) {
    //
    IO::Marker work_marker("one-dimensional uncoupled rotors model initialization");
 
    typedef std::map<int, double>::const_iterator pit_t;

    // effective one-dimensional rotor potential
    //
    for(int r = 0; r < internal_size(); ++r) {
      //
      std::map<int, double> one_pot_four;

      std::vector<int> pv(internal_size(), 0);

      for(int n = 1; n < _pot_four_index.size(r); ++n) {
	//
	pv[r] = n;

	int pl = _pot_four_index(pv);

	pit_t pit = _pot_four.find(pl);

	if(pit != _pot_four.end())
	  //
	  one_pot_four[n] = pit->second;
      }
      
      // initializing one-dimensional rotor
      //
      one_rotor.push_back(SharedPointer<HinderedRotor>(new HinderedRotor(one_pot_four, _imm_four[0](r, r), symmetry(r))));

      (*one_rotor.rbegin())->set(_extra_ener);
    }

    // ground energy for one-dimensional rotors model
    //
    one_ground = 0.;

    for(int r = 0; r < internal_size(); ++r)
      //
      one_ground += one_rotor[r]->ground() - one_rotor[r]->potential_minimum();

    one_ground += _pot_global_min;

    IO::log << IO::log_offset << "one-dimensional uncoupled rotors model ground energy [kcal/mol] = "
	    << one_ground / Phys_const::kcal << "\n";
  }

  /************************************* INTERPOLATION ********************************************/

  Slatec::Spline one_qstates; // one-dimensional rotors model quantum   density/number of states

  Slatec::Spline one_cstates; // one_dimensional rotors model classical density/number of states
 
  // quantum correction factor
  //
  if(_level_ener_max > 0. && (mode() != NOSTATES || force_qfactor)) {
    //
    _set_qfactor();
  }
  else {
    //
    IO::log << IO::log_offset << "WARNING: harmonic approximation will be used for ground state energy evaluation\n";
  }

  // external rotation quantum treatment
  //
  if(_amom_max)
    //
    rotational_energy_levels();

  // classical number of states, one-dimensional rotors model
  //
  if(mode() != NOSTATES || force_qfactor) {
    //
    IO::Marker interpol_marker("interpolating states number/density");
    
    /************************* CLASSICAL DENSITY/NUMBER OF STATES *****************************/

    typedef std::map<int, double>::const_iterator pit_t;

    itemp = (int)std::ceil(_extra_ener / _ener_quant);

    Array<double> ener_grid(itemp);

    Array<double> stat_base(itemp);

    Array<double> stat_grid(itemp, 0.);

    // energy grid
    //
    double ener = 0.;

    for(int i = 0; i < ener_grid.size(); ++i, ener += _ener_quant)
      //
      ener_grid[i] = ener;

    // base number of states
    //
    _set_states_base(stat_base);

    //angular integration
    //
    IO::log << IO::log_offset << "angular integration ... ";

    if(_vib_four.size()) {

#pragma omp parallel for default(shared) private(itemp, dtemp) schedule(dynamic, 1)      

      for(int g = 0; g < _grid_index.size(); ++g) {
	//
	// local vibrations contribution
	//
	Array<double> stat_freq = stat_base;

	for(int v = 0; v < _vib_four.size(); ++v) {
	  //
	  itemp = (int)round(_vib_grid[g][v] / _ener_quant);

	  for(int i = itemp; i < ener_grid.size(); ++i)
	    //
	    stat_freq[i] += stat_freq[i - itemp];
	}
	
	// potential energy shift and mass factor
	//
	itemp = (int)round((_pot_grid[g] - _pot_global_min) / _ener_quant);

	for(int i = itemp; i < ener_grid.size(); ++i) {
	  //
	  dtemp = stat_freq[i - itemp] * _irf_grid[g];

	  if(_with_ext_rot)
	    //
	    dtemp *= _erf_grid[g];

#pragma omp atomic

	  stat_grid[i] += dtemp;
	}
      }
    }
    // no vibrations
    //
    else {
      //
      // potential energy shift and mass factor
      //
      std::map<int, double> shift_factor;

      for(int g = 0; g < _grid_index.size(); ++g) {
	//
	itemp = (int)round((_pot_grid[g] - _pot_global_min) / _ener_quant);

	dtemp = _irf_grid[g];

	if(_with_ext_rot)
	  //
	  dtemp *= _erf_grid[g];

	shift_factor[itemp] += dtemp;
      }

#pragma omp parallel for default(shared) private(itemp, dtemp)  schedule(static)

      for(int i = 0; i < ener_grid.size(); ++i) {
	//
	dtemp = 0.;
    
	for(pit_t pit = shift_factor.begin(); pit != shift_factor.end(); ++pit) {
	  //
	  itemp = i - pit->first;

	  if(itemp >= 0)
	    //
	    dtemp  += stat_base[itemp] * pit->second;
	}

	stat_grid[i] = dtemp;
      }//
      //
    }// no vibrations
    
    // angular integration done
    //
    IO::log << "done\n";

    // classical density/number of states interpolation
    //
    IO::log << IO::log_offset << "classical states interpolation ... ";

    _cstates.init(ener_grid, stat_grid, ener_grid.size());

    // classical density/number of states extrapolation
    //
    dtemp = _cstates.arg_max() * (1. - _extra_step);

    dtemp = _cstates.fun_max() / _cstates(dtemp);

    _cstates_pow = std::log(dtemp) / std::log(1. / (1. - _extra_step));

    IO::log << "done\n";
    
    IO::log << IO::log_offset << "effective power exponent at " 
	    << _cstates.arg_max() / Phys_const::kcal << " kcal/mol = "<< _cstates_pow << "\n";

    /***************** ONE-DIMENSIONAL UNCOUPLED ROTORS MODEL DENSITY/NUMBER OF STATES ********************/

    // rotational density/number of states
    //
    if(_with_ext_rot) {
      //
      dtemp = 4. * M_SQRT2 * _erf_four[0] / external_symmetry();

      if(mode() == NUMBER)
	//
	dtemp *= 2. / 3.;

      ener = 0.;

      for(int i = 0; i < ener_grid.size(); ++i, ener += _ener_quant) {
	//
	stat_grid[i] = dtemp * std::sqrt(ener);

	if(mode() == NUMBER)
	  //
	  stat_grid[i] *= ener;
      }
    }
    // no external rotation
    //
    else {
      //
      switch(mode()) {
	//
      case NUMBER:
	//
	stat_grid = 1.;

	break;

      default:
	//
	stat_grid = 0.;

	stat_grid[0] = 1. / _ener_quant;
      }
    }

    // convolution with one-dimensional uncoupled rotors
    //
    for(int r = 0; r < internal_size(); ++r)
      //
      one_rotor[r]->convolute(stat_grid, _ener_quant);

    // convolution with average vibrations
    //
    for(int v = 0; v < _vib_four.size(); ++v) {
      //
      itemp = (int)round(_vib_four[v][0] / _ener_quant);

      for(int i = itemp; i < ener_grid.size(); ++i)
	//
	stat_grid[i] += stat_grid[i - itemp];
    }

    // effective one-dimensional rotors quantum density/number of states interpolation
    //
    one_qstates.init(ener_grid, stat_grid, ener_grid.size());

    // one-dimensional rotors model classical density/number of states (relative to potential energy minimum)
    //
    stat_grid = stat_base;

    dtemp = 1.;

    for(int r = 0; r < internal_size(); ++r)
      //
      one_rotor[r]->integrate(stat_grid, _ener_quant);

    // normalization
    //
    dtemp = 1.;

    for(int r = 0; r < internal_size(); ++r)
      //
      dtemp /= 2. * one_rotor[r]->rotational_constant();

    dtemp = std::sqrt(dtemp) / _angle_grid_cell;

    if(_with_ext_rot)
      //
      dtemp *= _erf_four[0];  
  
    stat_grid *= dtemp;

    // convolution with average vibrations
    //
    for(int v = 0; v < _vib_four.size(); ++v) {
      //
      itemp = (int)round(_vib_four[v][0] / _ener_quant);

      for(int i = itemp; i < ener_grid.size(); ++i)
	//
	stat_grid[i] += stat_grid[i - itemp];
    }

    // one-dimensional rotors classical density/number of states interpolation
    //
    one_cstates.init(ener_grid, stat_grid, ener_grid.size());

    // density/number of states output
    //
    if(_qfactor.size()) {
      //
      IO::log << IO::log_offset;

      switch(mode()) {
	//
      case NUMBER:
	//
	IO::log << "number of states";

	break;

      default:
	//
	IO::log << "density of states [cm]";
      }

      IO::log << " (E - energy in kcal/mol relative to the ground level):\n";

      IO::log << IO::log_offset << std::setw(5)  << "E" 
	      << std::setw(13) << "Classical"
	      << std::setw(13) << "Q-corrected"
	      << std::setw(13) << "Q-Factor"
	      << std::setw(13) << "1D-Classical"
	      << std::setw(13) << "1D-Quantum"
	      << std::setw(13) << "1D-Q-Factor"
	      << "\n";

      double stat_unit = 1.;

      if(mode() != NUMBER)
	//
	stat_unit = Phys_const::incm;

      double estep = 0.1 * Phys_const::kcal;

      for(ener = estep; ener < _qfactor.arg_max(); ener += estep) {
	//
	// classical density/number of states
	//
	IO::log << IO::log_offset 
		<< std::setw(5) << ener / Phys_const::kcal
		<< std::setw(13) << _cstates(ener + _ground - _pot_global_min) * stat_unit;

	// quantum-factor-corrected density/number of states
	//
	IO::log << std::setw(13) << states(ener) * stat_unit
		<< std::setw(13) << _qfactor(ener);

	// one-dimensional rotors model classical density/number of states
	//
	IO::log << std::setw(13) << one_cstates(ener + _ground - _pot_global_min) * stat_unit;

	// one-dimensional rotors model quantum density/number of states
	//
	dtemp = ener + _ground - one_ground;

	if(dtemp > 0.) {
	  //
	  IO::log << std::setw(13) << one_qstates(dtemp) * stat_unit
		  << std::setw(13) << one_qstates(dtemp) / one_cstates(ener + _ground - _pot_global_min);
	}
	else {
	  //
	  IO::log << std::setw(13) << 0
		  << std::setw(13) << 0;
	}

	IO::log << "\n";
      }
    }
  }// number of states output

  // statistical weight output
  //
  IO::log << IO::log_offset << "statistical weight (*** - deep tunneling regime):\n";

  IO::log << IO::log_offset 
	  << std::setw(5) << "T, K" 
	  << std::setw(13) << "Classical";

  if(_full_quantum_treatment && _level_ener_max > 0. && (mode() != NOSTATES || force_qfactor))
    //
    IO::log << std::setw(13) << "Quantum";

  IO::log << std::setw(13) << "PathIntegral" << "\n";

  for(int t = 50; t <= 1000 ; t+= 50) {
    //
    double tval = (double)t * Phys_const::kelv;

    double cw, sw;

    itemp = get_semiclassical_weight(tval, cw, sw);

    IO::log << IO::log_offset 
	    << std::setw(5) << t
	    << std::setw(13) << cw;

    if(_full_quantum_treatment && _level_ener_max > 0. && (mode() != NOSTATES || force_qfactor))
      //
      IO::log << std::setw(13) << quantum_weight(tval);

    IO::log << std::setw(13) << sw;

    if(itemp)
      //
      IO::log << "  ***";

    IO::log << "\n";  
    //
  }// temperature cycle
  //
}// MultiRotor

Model::MultiRotor::~MultiRotor ()
{
  //std::cout << "Model::MultiRotor destroyed\n";
}

void Model::MultiRotor::_set_qfactor ()
{
  const char funame [] = "Model::MultiRotor::_set_qfactor: ";

  IO::Marker funame_marker(funame);

  typedef std::map<int, double>::const_iterator                  pit_t;

  typedef std::map<int, Lapack::SymmetricMatrix>::const_iterator mit_t;

  int            itemp;
  double         dtemp;
  bool           btemp;
  Lapack::Vector vtemp;

  std::vector<int>    ivec;
  std::vector<double> dvec;

  std::clock_t  start_time;

  IO::log << IO::log_offset << "potential fourier expansion size = " << _pot_four.size()   << "\n";

  if(_vib_four.size()) {
    //
    IO::log << IO::log_offset << "vibrational frequencies fourier expansion sizes:";

    for(int v = 0; v < _vib_four.size(); ++v)
      //
      IO::log << "   " << _vib_four[v].size();

    IO::log << "\n";
  }

  IO::log << IO::log_offset << "internal mobility matrix fourier expansion size = " 
	  << _imm_four.size() << "\n";

  IO::log << IO::log_offset << "external rotation factor fourier expansion size = " 
	  << _erf_four.size() << "\n";

  // quantum state dimensions
  //
  ivec.resize(internal_size());

  for(int r = 0; r < internal_size(); ++r) {
    //
    itemp = int(std::sqrt(_level_ener_max / _mobility_parameter[r])) / symmetry(r) * 2 + 1; 

    if(_internal_rotation[r].quantum_size_max() && itemp > _internal_rotation[r].quantum_size_max()) {
      //
      ivec[r] = _internal_rotation[r].quantum_size_max();
    }
    else if(itemp < _internal_rotation[r].quantum_size_min()) {
      //
      ivec[r] = _internal_rotation[r].quantum_size_min();
    }
    else
      //
      ivec[r] = itemp;
  }

  MultiIndexConvert quantum_index(ivec);

  IO::log << IO::log_offset << "quantum phase space dimensions:";

  for(int i = 0; i < ivec.size(); ++i)
    //
    IO::log << "   " << ivec[i];

  IO::log << "\n";

  IO::log << IO::log_offset << "reference hamiltonian size = " << quantum_index.size() << std::endl;

  // vibrational population vector
  //
  std::vector<std::vector<int> > vib_pop;

  if(_vib_four.size()) {
    //
    dvec.resize(_vib_four.size());

    for(int v = 0; v < _vib_four.size(); ++v) {
      //
      dvec[v] = _vib_four[v][0];
    }

    vib_pop = population(_level_ener_max, dvec);

    for(int v = 0; v < _vib_four.size(); ++v) {
      //
      if(vib_pop[0][v] != 0) {
	//
	std::cerr << funame << "first population vector should describe the ground vibrational state\n";

	throw Error::Logic();
      }
    }

    IO::log << IO::log_offset << "number of vibrational population states = " << vib_pop.size() << std::endl;
  }

  /*
  // pruned potential on the grid
  Lapack::Vector pruned_pot_grid(_grid_index.size());
  double pruned_pot_global_min;

  if(_ptol > 0.) {

  #pragma omp parallel for default(shared) schedule(static)

  for(int g = 0; g < _grid_index.size(); ++g) {// grid cycle
    
  // current angle
  std::vector<int> gv = _grid_index(g);
  std::vector<double> angle(internal_size());
  for(int r = 0; r < internal_size(); ++r)
  angle[r] = double(gv[r]) * _angle_grid_step[r];
    
  // potential
  pruned_pot_grid[g]  = potential(angle);
  }// grid cycle

  // pruned potential global minimum
  int gmin = 0;
  for(int g = 0; g < _grid_index.size(); ++g) {
  dtemp =  pruned_pot_grid[g];
  if(!g || dtemp < pruned_pot_global_min) {
  pruned_pot_global_min = dtemp;
  gmin = g;
  }
  }

  // correction to the potential
  ivec = _grid_index(gmin);
  std::vector<double> angle_min(internal_size());
  for(int r = 0; r < internal_size(); ++r)
  angle_min[r] = double(ivec[r]) * _angle_grid_step[r];

  try {
  vtemp = Lapack::Cholesky(force_constant_matrix(angle_min)).invert(potential_gradient(angle_min));
  for(int i = 0; i < internal_size(); ++i)
  angle_min[i] -= vtemp[i];

  dtemp = potential(angle_min) - pruned_pot_global_min;
  pruned_pot_global_min += dtemp;
  IO::log << IO::log_offset << "pruned potential minimum energy correction [kcal/mol] = "
  << dtemp / Phys_const::kcal << "\n";
  }
  catch(Error::Math) {
  IO::log << IO::log_offset << "WARNING: force constant matrix at minimum is not positive definite: "
  "cannot get potential minimum energy correction\n";
  }

  IO::log << IO::log_offset
  << "pruned potential global minimum [kcal/mol] = " 
  << pruned_pot_global_min / Phys_const::kcal << "\n";
  }
  else {

  #pragma omp parallel for default(shared) schedule(static)

  for(int g = 0; g < _grid_index.size(); ++g) {
  pruned_pot_grid[g] = _pot_grid[g];
  }

  pruned_pot_global_min = _pot_global_min;
  }

  std::vector<double> pruned_irf_grid(_grid_index.size());
  std::vector<double> pruned_erf_grid;
  if(_with_ext_rot)
  pruned_erf_grid.resize(_grid_index.size());

  if(_mtol > 0.) {

  #pragma omp parallel for default(shared) private(itemp, dtemp) schedule(static)

  for(int g = 0; g < _grid_index.size(); ++g) {// grid cycle

  // current angle
  std::vector<int> gv = _grid_index(g);
  std::vector<double> angle(internal_size());
  for(int r = 0; r < internal_size(); ++r)
  angle[r] = double(gv[r]) * _angle_grid_step[r];

  // mass factor
  try {
  pruned_irf_grid[g] =  Lapack::Cholesky(mass(angle)).det_sqrt();
  } 
  catch(Error::General) {
  std::cerr << funame 
  << "interpolated mass matrix is not positive definite: decrease mass tolerance\n";
  throw;
  }

  // external rotation mass factor correction
  if(_with_ext_rot) {
  dtemp = external_rotation_factor(angle);
  if(dtemp < 0.) {
  std::cerr << funame << "negative external rotation factor at (";
  for(int r = 0; r < internal_size(); ++r) {
  if(r)
  std::cerr << ", ";
  std::cerr << angle[r] * 180. / M_PI;
  }
  std::cerr << ") grid point\n";
  throw Error::Range();
  }
  pruned_erf_grid[g] = dtemp;
  }
  }// grid cycle
  }
  else {

  #pragma omp parallel for default(shared) schedule(static)

  for(int g = 0; g < _grid_index.size(); ++g) {
  pruned_irf_grid[g] = _irf_grid[g];
  if(_with_ext_rot)
  pruned_erf_grid[g]  = _erf_grid[g];
  }
  }
  */

  // classical number of states without external rotations at maximum energy
  //
  itemp = (int)std::ceil(_level_ener_max / _ener_quant);

  Array<double> stat_base(itemp);

  _set_states_base(stat_base, 1);

  dtemp = 0.;

  for(int g = 0; g < _grid_index.size(); ++g) {
    //
    itemp = stat_base.size() - 1 - (int)round((_pot_grid[g] - _pot_global_min) / _ener_quant);

    if(itemp >= 0)
      //
      dtemp += stat_base[itemp] * _irf_grid[g];
  }

  const int nos_max = (int)dtemp;
  
  IO::log << IO::log_offset << "classically estimated number of energy levels = " << nos_max << "\n";

  // ground vibrational state hamiltonian
  //
  Lapack::SymmetricMatrix ham_base(quantum_index.size());

  if(1) {
    //
    IO::Marker basic_marker("setting reference hamiltonian", IO::Marker::ONE_LINE);

    if(_with_ctf) {
      //
      ham_base = 0.;
    }
    else
      //
      ham_base = _pot_shift;
      
    const std::map<int, double>*                  eff_pot = &_pot_four;

    const std::map<int, Lapack::SymmetricMatrix>* eff_imm = &_imm_four;

    if(_with_ctf) {
      //
      eff_pot = &_eff_pot_four;
      
      eff_imm = &_eff_imm_four;
    }
    
#pragma omp parallel for default(shared) private(btemp, dtemp, itemp) schedule(dynamic, 1)

    for(int ml = 0; ml < quantum_index.size(); ++ml) {// left quantum state cycle
      //
      const std::vector<int> mv = quantum_index(ml);

      for(int nl = ml; nl < quantum_index.size(); ++nl) {// right quantum state cycle
	//
	const std::vector<int> nv = quantum_index(nl);

	// potential contribution
	//
	for(pit_t pit = eff_pot->begin(); pit != eff_pot->end(); ++pit) {// fourier expansion cycle
	  //
	  const std::vector<int> pv = _pot_four_index(pit->first);

	  btemp = true;
	  
	  dtemp = pit->second;

	  for(int r = 0; r < internal_size(); ++r) {
	    //
	    if(rotation_matrix_element(mv[r], nv[r], pv[r], dtemp)) {
	      //
	      btemp = false;

	      break;
	    }
	  }

	  if(btemp) 
	    //
	    ham_base(ml, nl) +=  dtemp;
	  //
	}// fourier expansion cycle

	//kinetic energy contribution
	//
	for(int mi = 0; mi < internal_size(); ++mi) {// mobility matrix first index cycle
	  //
	  if(!mv[mi])
	    //
	    continue;

	  int mfac = (mv[mi] + 1) / 2 * symmetry(mi);

	  std::vector<int> mw = mv;

	  // sin -> cos
	  //
	  if(mw[mi] % 2) {
	    //
	    ++mw[mi];
	  }
	  // cos -> sin
	  //
	  else {
	    //
	    --mw[mi];

	    mfac = -mfac;
	  }

	  for(int ni = 0; ni < internal_size(); ++ni) {// mobility matrix second index cycle
	    //
	    if(!nv[ni])
	      //
	      continue;

	    int nfac = (nv[ni] + 1) / 2 * symmetry(ni);

	    std::vector<int> nw = nv;

	    // sin -> cos
	    //
	    if(nw[ni] % 2) {
	      //
	      ++nw[ni];
	    }
	    // cos -> sin
	    //
	    else {
	      //
	      --nw[ni];

	      nfac = -nfac;
	    }
	     
	    nfac *= mfac;
	    
	    for(mit_t mit = eff_imm->begin(); mit != eff_imm->end(); ++mit) {// fourier expansion cycle
	      //
	      std::vector<int> kv = _mass_index(mit->first);

	      btemp = true;

	      dtemp = (double)nfac * mit->second(mi, ni);

	      for(int r = 0; r < internal_size(); ++r) {// internal rotation cycle
		//
		if(rotation_matrix_element(mw[r], nw[r], kv[r], dtemp)) {
		  //
		  btemp = false;

		  break;
		  //
		}//
		//
	      }// internal rotation cycle

	      if(btemp)
		//
		ham_base(ml, nl) += dtemp;
	      //
	    }// fourier expansion cycle      
	    //
	  }// mobility matrix second index cycle
	  //
	}// mobility matrix first index cycle
	//
      }// bra vector cycle
      //
    }// ket vector cycle
    //
  }// ground vibrational state hamiltonian

  // curvlinear transformation factor matrix
  //
  Lapack::SymmetricMatrix ctf_base;

  if(_with_ctf) {
    //
    IO::Marker ext_rot_marker("setting curvlinear transformation factor matrix elements", IO::Marker::ONE_LINE);

    ctf_base.resize(quantum_index.size());

    ctf_base = 0.;

#pragma omp parallel for default(shared) private(btemp, dtemp) schedule(dynamic, 1)

    for(int ml = 0; ml < quantum_index.size(); ++ml) {// bra cycle
      //
      std::vector<int> mv = quantum_index(ml);

      for(int nl = ml; nl < quantum_index.size(); ++nl) {// ket cycle
	//
	std::vector<int> nv = quantum_index(nl);

	for(pit_t pit = _ctf_four.begin(); pit != _ctf_four.end(); ++pit) {// fourier expansion cycle
	  //
	  std::vector<int> pv = _mass_index(pit->first);

	  btemp = true;

	  dtemp = pit->second;

	  for(int r = 0; r < internal_size(); ++r) {// internal rotation cycle
	    //
	    if(rotation_matrix_element(mv[r], nv[r], pv[r], dtemp)) {
	      //
	      btemp = false;

	      break;
	    }
	    //
	  }// internal rotation cycle

	  if(btemp)
	    //
	    ctf_base(ml, nl) +=  dtemp;
	  //
	}// fourier expansion cycle
	//
      }// ket cycle
      //
    }// bra cycle
    //
  }// internal rotation factor matrix elements

  // external rotation factor matrix
  //
  Lapack::SymmetricMatrix erf_base;

  if(_with_ext_rot) {
    //
    IO::Marker ext_rot_marker("setting external rotation factor matrix elements", IO::Marker::ONE_LINE);

    erf_base.resize(quantum_index.size());

    erf_base = 0.;

    const std::map<int, double>* eff_erf = &_erf_four;

    if(_with_ctf)
      //
      eff_erf = &_eff_erf_four;
      
#pragma omp parallel for default(shared) private(btemp, dtemp) schedule(dynamic, 1)

    for(int ml = 0; ml < quantum_index.size(); ++ml) {// bra cycle
      //
      std::vector<int> mv = quantum_index(ml);

      for(int nl = ml; nl < quantum_index.size(); ++nl) {// ket cycle
	//
	std::vector<int> nv = quantum_index(nl);

	for(pit_t pit = eff_erf->begin(); pit != eff_erf->end(); ++pit) {// fourier expansion cycle
	  //
	  std::vector<int> pv = _mass_index(pit->first);

	  btemp = true;

	  dtemp = pit->second;

	  for(int r = 0; r < internal_size(); ++r) {// internal rotation cycle
	    //
	    if(rotation_matrix_element(mv[r], nv[r], pv[r], dtemp)) {
	      //
	      btemp = false;

	      break;
	    }
	    //
	  }// internal rotation cycle

	  if(btemp)
	    //
	    erf_base(ml, nl) +=  dtemp;
	  //
	}// fourier expansion cycle
	//
      }// ket cycle
      //
    }// bra cycle
    //
  }// external rotation factor matrix elements

  // vibrationally shifted energy levels
  //
  if(_vib_four.size() && _full_quantum_treatment) {
    //
    _energy_level.resize(vib_pop.size());

    _mean_erf.resize(vib_pop.size());
  }
  else {
    //
    _energy_level.resize(1);

    _mean_erf.resize(1);
  }

  // setting energy levels
  //
  for(int vp = 0; vp < _energy_level.size(); ++vp) {// vibrational population cycle
    //
    if(_vib_four.size() && _full_quantum_treatment) {
      //
      // output: vibration population vector
      //
      IO::log << IO::log_offset << "vibrational population vector = [";

      for(int v = 0; v < _vib_four.size(); ++v) { 
	//
	if(v)
	  //
	  IO::log << ", ";

	IO::log << vib_pop[vp][v];
      }

      IO::log << "]\n";

      IO::log_offset.increase();
    }

    // maximum level energy shifted by vibrational energy
    //
    dtemp = _level_ener_max;

    for(int v = 0; v < _vib_four.size(); ++v)
      //
      dtemp -= _vib_four[v][0] * double(vib_pop[vp][v]);

    const double ener_max = dtemp;

    if(ener_max <= 0.)
      //
      continue;

    // stripping high harmonics
    //
    std::vector<int> s2u;// stripped-to-basic index converter

    dvec.resize(internal_size());

    for(int ml = 0; ml < quantum_index.size(); ++ml) {
      //
      ivec = quantum_index(ml);

      btemp = false;

      for(int i = 0; i < internal_size(); ++i) {// internal rotation cycle
	//
	if(ivec[i] < _internal_rotation[i].quantum_size_min()) {
	  //
	  btemp = true;

	  break;
	}

	dvec[i] = double((ivec[i] + 1) / 2 * symmetry(i));
	//
      }// internal rotation cycle

      if(btemp) {
	//
	s2u.push_back(ml);
      }
      else {
	//
	dtemp = 0.;

	for(int i = 0; i < internal_size(); ++i) {
	  //
	  for(int j = i; j < internal_size(); ++j) {
	    //
	    if(i == j) {
	      //
	      dtemp += _mobility_min(i, i) * dvec[i] * dvec[i];
	    }
	    else
	      //
	      dtemp += 2. * _mobility_min(i, j) * dvec[i] * dvec[j];
	  }
	}

	if(dtemp < ener_max)
	  //
	  s2u.push_back(ml);
	//
      }//
      //
    }//
    
    IO::log << IO::log_offset << "stripped hamiltonian size = " << s2u.size() << std::endl;

    // setting hamiltonian
    //
    Lapack::SymmetricMatrix ham(s2u.size());

#pragma omp parallel for default(shared) schedule(static)
      
    for(int i = 0; i < s2u.size(); ++i) {
      //
      for(int j = i; j < s2u.size(); ++j) {
	//
	ham(i, j) = ham_base(s2u[i], s2u[j]);
	//
      }//
      //
    }//

    // setting scalar product matrix
    //
    Lapack::SymmetricMatrix ctf_mat;

    if(_with_ctf) {
      //
      ctf_mat.resize(s2u.size());

#pragma omp parallel for default(shared) schedule(static)
      
      for(int i = 0; i < s2u.size(); ++i) {
	//
	for(int j = i; j < s2u.size(); ++j) {
	  //
	  ctf_mat(i, j) = ctf_base(s2u[i], s2u[j]);
	  //
	}//
	//
      }//
      //
    }//
    
    // excited vibrational states correction
    //
    if(vp) {
      //
      IO::Marker ham_marker("setting hamiltonian", IO::Marker::ONE_LINE);

      const std::vector<std::map<int, double> >* eff_vib = &_vib_four;

      if(_with_ctf)
	//
	eff_vib = &_eff_vib_four;
      
#pragma omp parallel for default(shared) private(btemp, dtemp) schedule(dynamic, 1)
      
      for(int ml = 0; ml < s2u.size(); ++ml) {// ket quantum state cycle
	//
	std::vector<int> mv = quantum_index(s2u[ml]);
	
	for(int nl = ml; nl < s2u.size(); ++nl) {// bra quantum state cycle
	  //
	  std::vector<int> nv = quantum_index(s2u[nl]);

	  for(int v = 0; v < eff_vib->size(); ++v) {// vibrational frequencies cycle
	    //
	    for(pit_t pit = (*eff_vib)[v].begin(); pit != (*eff_vib)[v].end(); ++pit) {// fourier expansion cycle
	      //
	      std::vector<int> pv = _pot_four_index(pit->first);
	      
	      btemp = true;

	      dtemp =  pit->second * double(vib_pop[vp][v]);

	      for(int r = 0; r < internal_size(); ++r) {
		//
		if(rotation_matrix_element(mv[r], nv[r], pv[r], dtemp)) {// internal rotation cycle
		  //
		  btemp = false;

		  break;
		}
	      }

	      if(btemp)
		//
		ham(ml, nl) += dtemp;
	      //
	    }// fourier expansion cycle 
	    //
	  }// vibration frequencies cycle
	  //
	}// bra quantum state cycle  
	//
      }// ket quantum state cycle
      //
    }// vibrational excited state correction

    // diagonalizing hamiltonian
    //
    Lapack::Vector eval;

    Lapack::Matrix evec;

    if(_with_ctf) {
      //
      IO::Marker diag_marker("diagonalizing hamiltonian and ctf matrix simultaneously", IO::Marker::ONE_LINE);

      if(_with_ext_rot) {
	//
	eval = Lapack::diagonalize(ham, ctf_mat, &evec);
      }
      else {
	//
	eval = Lapack::diagonalize(ham, ctf_mat);
      }
    }
    else {
      //
      IO::Marker diag_marker("diagonalizing hamiltonian", IO::Marker::ONE_LINE);

      if(_with_ext_rot) {
	//
	eval = ham.eigenvalues(&evec);
      }
      else
	//
	eval = ham.eigenvalues();
    }
    
    // ground energy level
    //
    if(!vp)
      //
      _ground = eval[0];
      
    // quantum states initialization
    //
    for(int l = 0; l < ham.size(); ++l) {
      //
      if(eval[l] > _level_ener_max + _pot_global_min)
	//
	break;

      _energy_level[vp].push_back(eval[l] - _ground);
    }
	
    // external rotation factor quantum state average
    //
    Lapack::SymmetricMatrix erf_mat;

    if(_with_ext_rot) {
      //
      IO::Marker erf_marker("averaging external rotation factor over quantum state", IO::Marker::ONE_LINE);
      
      erf_mat.resize(s2u.size());
      
#pragma omp parallel for default(shared) schedule(static)
      
      for(int m = 0; m < s2u.size(); ++m) {
	//
	for(int n = m; n < s2u.size(); ++n)
	  //
	  erf_mat(m, n) = erf_base(s2u[m], s2u[n]);
      }

      _mean_erf[vp].resize(_energy_level[vp].size());

      for(int l = 0; l < _mean_erf[vp].size(); ++l) {
	//
	vtemp = erf_mat * &evec(0,l);

	_mean_erf[vp][l] = parallel_vdot(vtemp, &evec(0, l), s2u.size());
	//
      }//
      //
    }// external rotation quantum state average

    IO::log << IO::log_offset << "quantum states  # = " << _energy_level[vp].size() << std::endl;

    // output: energy levels
    //
    itemp = _energy_level[vp].size() < 10 ? _energy_level[vp].size() : 10;

    if(itemp) {
      //
      IO::log << IO::log_offset << "first " << itemp << " energy levels, kcal/mol:";

      for(int i = 0; i < itemp; ++i)
	//
	IO::log  << " " << (_ground + _energy_level[vp][i])/ Phys_const::kcal;

      IO::log << std::endl;
    }

    if(_vib_four.size() && _full_quantum_treatment)
      //
      IO::log_offset.decrease();
    //
  }// vibrational population cycle

  IO::log << IO::log_offset << "ground state energy      [kcal/mol] = " << _ground / Phys_const::kcal << std::endl;

  IO::log << IO::log_offset << "potential energy minimum [kcal/mol] = " << _pot_global_min / Phys_const::kcal << std::endl;

  /**************** QUANTUM DENSITY/NUMBER OF STATES ***************/

  itemp = stat_base.size() - (int)std::ceil((_ground - _pot_global_min) / _ener_quant);

  if(itemp < 3) {
    //
    IO::log << IO::log_offset << "WARNING: zero-point energy is bigger than quantum state energy maximum => no quantum correction" << std::endl;
    
    return;
  }

  Array<double> qstat_grid(itemp);

  quantum_states(qstat_grid, _ener_quant, 1);

  // convolution with average vibrational frequencies
  //
  for(int v = 0; v < _vib_four.size(); ++v) {
    //
    itemp = (int)round(_vib_four[v][0] / _ener_quant);

    for(int i = itemp; i < qstat_grid.size(); ++i)
      //
      qstat_grid[i] += qstat_grid[i - itemp];
  }

  /*************** CLASSICAL DENSITY/NUMBER OF STATES *******************/

  // reset base number of states
  //
  _set_states_base(stat_base);

  std::map<int, double> shift_factor;
  //
  for(int g = 0; g < _grid_index.size(); ++g) {
    //
    itemp = (int)round((_pot_grid[g] - _pot_global_min) / _ener_quant);

    dtemp = _irf_grid[g];

    if(_with_ext_rot)
      //
      dtemp *= _erf_grid[g];

    shift_factor[itemp] += dtemp;
  }

  Array<double> stat_grid(stat_base.size());

#pragma omp parallel for default(shared) private(itemp, dtemp)  schedule(static)

  for(int i = 0; i < stat_grid.size(); ++i) {
    //
    dtemp = 0.;

    for(pit_t pit = shift_factor.begin(); pit != shift_factor.end(); ++pit) {
      //
      itemp = i - pit->first;

      if(itemp >= 0)
	//
	dtemp  += stat_base[itemp] * pit->second;
    }

    stat_grid[i] = dtemp;
  }

  // convolution with average vibrational frequencies
  //
  for(int v = 0; v < _vib_four.size(); ++v) {
    //
    itemp = (int)round(_vib_four[v][0] / _ener_quant);

    for(int i = itemp; i < stat_grid.size(); ++i)
      //
      stat_grid[i] += stat_grid[i - itemp];
  }

  // energy grid
  //
  Array<double> ener_grid(stat_grid.size());

  double ener = 0.;

  for(int i = 0; i < ener_grid.size(); ++i, ener += _ener_quant)
    //
    ener_grid[i] = ener;

  // classical density/number of states interpolation
  //
  Slatec::Spline eff_cstates(ener_grid, stat_grid, ener_grid.size());

  /**************** QUANTUM CORRECTION FACTOR  ***************/

  // quantum correction factor
  //
  for(int i = 0; i < qstat_grid.size(); ++i) {
    //
    dtemp = ener_grid[i] + _ground - _pot_global_min;

    if(dtemp < eff_cstates.arg_max()) {
      //
      qstat_grid[i] /= eff_cstates(dtemp);
    }
    else
      //
      qstat_grid[i] = 1.;
  }

  // quantum correction factor interpolation
  //
  _qfactor.init(ener_grid, qstat_grid, qstat_grid.size());
}

void Model::MultiRotor::_set_states_base (Array<double>& base, int flag) const
{
  const char funame [] = "Model::MultiRotor::_set_states_base: ";

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE);

  bool btemp = _with_ext_rot;

  int mm = mode();

  if(flag) {
    //
    btemp = false;

    mm = NUMBER;
  }
  
  double dim, nfac;

  if(btemp) {
    //
    switch(mm) {
      //
    case NUMBER:
      //
      dim = double(internal_size() + 3) / 2.;

      nfac = 8. * M_PI * M_PI / std::pow(2. * M_PI, dim) 
	/ tgamma(dim + 1.) / external_symmetry();

      break;

    default:
      //
      dim = double(internal_size() + 1) / 2.;

      nfac = 8. * M_PI * M_PI / std::pow(2. * M_PI, dim + 1.) 
	/ tgamma(dim + 1.) / external_symmetry();
    }
  }
  else {
    //
    switch(mm) {
      //
    case NUMBER:
      //
      dim = double(internal_size()) / 2.;

      nfac = 1. / std::pow(2. * M_PI, dim) / tgamma(dim + 1.);

      break;

    default:
      //
      dim = double(internal_size() - 2) / 2.;

      nfac = 1. / std::pow(2. * M_PI, dim + 1.) / tgamma(dim + 1.);
    }
  }

  double ener = _ener_quant;

  base[0] = 0.;

  for(int i = 1; i < base.size(); ++i, ener += _ener_quant) {
    //
    base[i] = std::pow(ener, dim) * nfac * _angle_grid_cell;    
  }
}

// generalized mass matrix
//
Lapack::SymmetricMatrix Model::MultiRotor::mass (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::mass: ";

  typedef std::map<int, Lapack::SymmetricMatrix>::const_iterator mit_t;

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle size\n";

    throw Error::Range();
  }

  double dtemp;

  Lapack::SymmetricMatrix res(internal_size());

  res = 0.;

  for(mit_t mit = _imm_four.begin(); mit != _imm_four.end(); ++mit) {// fourier expansion cycle
    //
    std::vector<int> hvec = _mass_index(mit->first);

    for(int i = 0; i < internal_size(); ++i) {
      //
      for(int j = i; j < internal_size(); ++j) {
	//
	dtemp = mit->second(i, j);

	for(int k = 0; k < internal_size(); ++k) {
	  //
	  if(!hvec[k])
	    //
	    continue;

	  if(hvec[k] % 2) {
	    //
	    dtemp *= std::sin(double((hvec[k] + 1) / 2 * symmetry(k)) * angle[k]);
	  }
	  else
	    //
	    dtemp *= std::cos(double(hvec[k] / 2 * symmetry(k)) * angle[k]);
	}

	res(i, j) += dtemp;
      }//
      //
    }//
    //
  }//fourier expansion size

  res = res.invert();

  res /= 2.;

  return res;
}

// external rotation inertia moments product
//
double Model::MultiRotor::external_rotation_factor (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::external_rotation_factor: ";

  typedef std::map<int, double>::const_iterator pit_t;

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle size\n";

    throw Error::Range();
  }

  double dtemp;

  double res = 0.;

  for(pit_t pit = _erf_four.begin(); pit != _erf_four.end(); ++pit) {
    //
    std::vector<int> hvec = _mass_index(pit->first);

    dtemp = pit->second;

    for(int k = 0; k < internal_size(); ++k) {
      //
      if(!hvec[k])
	//
	continue;
      
      if(hvec[k] % 2) {
	//
	dtemp *= std::sin(double((hvec[k] + 1) / 2 * symmetry(k)) * angle[k]);
      }
      else
	//
	dtemp *= std::cos(double(hvec[k] / 2 * symmetry(k)) * angle[k]);
    }

    res += dtemp;
  }
  
  return res;
}

// internal rotation inertia moments product
//
double Model::MultiRotor::internal_rotation_factor (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::internal_rotation_factor: ";

  typedef std::map<int, double>::const_iterator pit_t;

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle size\n";

    throw Error::Range();
  }

  double dtemp;

  double res = 0.;

  for(pit_t pit = _irf_four.begin(); pit != _irf_four.end(); ++pit) {
    //
    std::vector<int> hvec = _mass_index(pit->first);

    dtemp = pit->second;

    for(int k = 0; k < internal_size(); ++k) {
      //
      if(!hvec[k])
	//
	continue;
      
      if(hvec[k] % 2) {
	//
	dtemp *= std::sin(double((hvec[k] + 1) / 2 * symmetry(k)) * angle[k]);
      }
      else
	//
	dtemp *= std::cos(double(hvec[k] / 2 * symmetry(k)) * angle[k]);
    }

    res += dtemp;
  }
  
  return res;
}

// curvlinear transformation factor (dx/dx' determinant)
//
double Model::MultiRotor::curvlinear_factor (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::curvlinear_factor: ";

  typedef std::map<int, double>::const_iterator pit_t;

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle size\n";

    throw Error::Range();
  }

  double dtemp;

  double res = 0.;

  for(pit_t pit = _ctf_four.begin(); pit != _ctf_four.end(); ++pit) {
    //
    std::vector<int> hvec = _mass_index(pit->first);

    dtemp = pit->second;

    for(int k = 0; k < internal_size(); ++k) {
      //
      if(!hvec[k])
	//
	continue;
      
      if(hvec[k] % 2) {
	//
	dtemp *= std::sin(double((hvec[k] + 1) / 2 * symmetry(k)) * angle[k]);
      }
      else
	//
	dtemp *= std::cos(double(hvec[k] / 2 * symmetry(k)) * angle[k]);
    }

    res += dtemp;
  }
  
  return res;
}

// vibrational frequencies
//
Lapack::Vector  Model::MultiRotor::vibration (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::vibration: ";

  typedef std::map<int, double>::const_iterator pit_t;

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle size\n";

    throw Error::Range();
  }

  if(!_vib_four.size()) {
    //
    std::cerr << funame << "no vibrations\n";

    throw Error::Logic();
  }

  double dtemp;

  Lapack::Vector res(_vib_four.size());

  res = 0.;

  for(int v = 0; v < _vib_four.size(); ++v) {// vibrational frequencies cycle
    //
    for(pit_t pit = _vib_four[v].begin(); pit != _vib_four[v].end(); ++pit) {// fourier expansion cycle
      //
      std::vector<int> hvec = _pot_four_index(pit->first);

      dtemp = 1.;

      for(int k = 0; k < internal_size(); ++k) {
	//
	if(!hvec[k])
	  //
	  continue;

	if(hvec[k] % 2) {
	  //
	  dtemp *= std::sin(double((hvec[k] + 1) / 2 * symmetry(k)) * angle[k]);
	}
	else {
	  //
	  dtemp *= std::cos(double(hvec[k] / 2 * symmetry(k)) * angle[k]);
	}//
	//
      }//

      res[v] += dtemp * pit->second;
      //
    }// fourier expansion cycle
    //
  }// vibrational frequencies cycle

  return res;
}

// derivatives of the potential
//
double Model::MultiRotor::potential (const std::vector<double>& angle, const std::map<int, int>& der) const
{
  const char funame [] = "Model::MultiRotor::potential: ";

  typedef std::map<int, double>::const_iterator pit_t;

  double dtemp;

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle dimension\n";

    throw Error::Logic();
  }

  std::map<int, int>::const_iterator dit;

  for(dit = der.begin(); dit != der.end(); ++dit) {
    //
    if(dit->second <= 0) {
      //
      std::cerr << funame << "derivative order out of range: " << dit->second << "\n";

      throw Error::Range();
    }
  }

  if(der.size() && (der.begin()->first < 0 || der.rbegin()->first >= internal_size())) {
    //
    std::cerr << funame << "variable index out of range\n";

    throw Error::Range();
  }

  double res = 0.;

  for(pit_t pit = _pot_four.begin(); pit != _pot_four.end(); ++pit) {
    //
    std::vector<int> hvec = _pot_four_index(pit->first);

    double der_value = pit->second;

    for(int k = 0; k < internal_size(); ++k) {
      //
      int der_order = 0;

      dit = der.find(k);

      if(dit != der.end())
	//
	der_order = dit->second;

      if(hvec[k]) {
	//
	int n =  (hvec[k] + 1) / 2  * symmetry(k);

	int ifac = 1;

	for(int i = 0; i < der_order; ++i)
	  //
	  ifac *= n;

	if((hvec[k] + der_order) % 2) {
	  //
	  der_value *= (double)ifac * std::sin((double)n * angle[k]);
	}
	else
	  //
	  der_value *= (double)ifac * std::cos((double)n * angle[k]);

    
	if((der_order + 1 - hvec[k] % 2) / 2 % 2)
	  //
	  der_value = -der_value;
      }
      else if(der_order) {
	//
	der_value = 0.;

	break;
      }
    }
    
    res += der_value;
  }
  
  if(!der.size())
    //
    res += _pot_shift;

  return res;
}

Lapack::SymmetricMatrix Model::MultiRotor::force_constant_matrix (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::force_constant_matrix: ";

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle dimension\n";

    throw Error::Logic();
  }

  Lapack::SymmetricMatrix res(internal_size());

  for(int i = 0; i < internal_size(); ++i) {
    //
    for(int j = i; j < internal_size(); ++j) {
      //
      std::map<int, int> der;

      if(i == j) {
	//
	der[i] = 2;
      }
      else
	//
	der[i] = der[j] = 1;

      res(i, j) = potential(angle, der);
      //
    }//
    //
  }//

  return res;
}

Lapack::Vector Model::MultiRotor::potential_gradient(const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::potential_gradient: ";

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle dimension\n";

    throw Error::Logic();
  }

  // potential gradient
  //
  Lapack::Vector res(internal_size());

  for(int i = 0; i < internal_size(); ++i) {
    //
    std::map<int, int> der;

    der[i] = 1;

    res[i] = potential(angle, der);
  }

  return res;
}

// local frequencies
//
Lapack::Vector Model::MultiRotor::frequencies (const std::vector<double>& angle) const
{
  const char funame [] = "Model::MultiRotor::frequencies: ";

  if(angle.size() != internal_size()) {
    //
    std::cerr << funame << "wrong angle size\n";

    throw Error::Range();
  }

  double dtemp;

  // frequencies
  //
  Lapack::Vector res = Lapack::diagonalize(force_constant_matrix(angle), mass(angle));

  for(int r = 0; r < internal_size(); ++r) {
    //
    dtemp = res[r];

    if(dtemp < 0.) {
      //
      res[r] = -std::sqrt(-dtemp);
    }
    else
      //
      res[r] = std::sqrt(dtemp);
  }
  
  return res;
}

double Model::MultiRotor::quantum_weight (double temperature) const
{
  static const double pi_fac = 2. * std::sqrt(2. * M_PI);

  double res = 0.;
  double dtemp;

  for(int p = 0; p < _energy_level.size(); ++p)
    for(int l = 0; l < _energy_level[p].size(); ++l) {
      dtemp = std::exp(-_energy_level[p][l] / temperature);
      if(_with_ext_rot) 
	res += dtemp * _mean_erf[p][l];
      else 
	res += dtemp;
    }

  if(_with_ext_rot)
    res *=  pi_fac * temperature * std::sqrt(temperature) / external_symmetry();

  return res;
}

int Model::MultiRotor::get_semiclassical_weight (double temperature, double& classical_weight, double& quantum_weight) const
{
  static const double eps = 1.e-5;
  
  static const double pi_fac = 2. * std::sqrt(2. * M_PI);

  int    itemp;
  double dtemp;

  int res = 0;
  double cw = 0.;
  double qw = 0.;

#pragma omp parallel for default(shared) private(itemp, dtemp) reduction(+: cw, qw) schedule(static)

  for(int g = 0; g < _grid_index.size(); ++g) {// grid cycle
    //
    // quantum correction factor
    //
    double qfac = 1.;
    
    for(int r = 0; r < internal_size(); ++r) {
      //
      dtemp = _freq_grid[g][r] / temperature / 2.;

      if(dtemp > eps) {
	//
	qfac *= dtemp / std::sinh(dtemp);
      }
      else if(dtemp < eps - M_PI) {
	//
	qfac = -1.;

	res = 1;

	break;
      }
      else if(dtemp < -eps)
	//
	qfac *= dtemp / std::sin(dtemp);
    }

    // classical partition  function (internal rotations & vibrations part)
    //
    dtemp = std::exp(-_pot_grid[g] / temperature) * _irf_grid[g];
    
    if(_with_ext_rot)
      //
      dtemp *= _erf_grid[g];

    for(int v = 0; v < _vib_four.size(); ++v)
      //
      dtemp /= 1. - std::exp(-_vib_grid[g][v] / temperature);
   
    cw += dtemp;
    
    if(qfac > 0.)
      //
      qw += qfac * dtemp;
    //
  }// grid cycle
    
  // normalization
  //
  dtemp = std::pow(temperature / 2. / M_PI, double(internal_size()) / 2.)
    //
    * _angle_grid_cell * std::exp(_ground / temperature);

  if(_with_ext_rot)
    //
    dtemp *= pi_fac * temperature * std::sqrt(temperature) / external_symmetry();

  classical_weight = dtemp * cw;
  
  quantum_weight   = dtemp * qw;

  return res;
}

// number/density of quantum states
//
void Model::MultiRotor::quantum_states (Array<double>& stat_grid, double ener_step, int flag) const
{
  const char funame [] = "Model::MultiRotor::quantum_states: ";

  if(!_energy_level.size()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  int    itemp;
  double dtemp;

  stat_grid = 0.;

  const int pmax = flag ? 1 : _energy_level.size();

  if(_with_ext_rot) {
    //
    for(int p = 0; p < pmax; ++p) {
      //
      for(int l = 0; l < _energy_level[p].size(); ++l)  {
	//
	itemp = (int)std::ceil(_energy_level[p][l] / ener_step);

	for(int i = itemp; i < stat_grid.size(); ++i) {
	  //
	  double rot_ener = (double)i * ener_step - _energy_level[p][l];

	  if(rot_ener <= 0.)
	    //
	    continue;

	  dtemp = std::sqrt(rot_ener) * _mean_erf[p][l];

	  if(mode() == NUMBER)
	    //
	    dtemp *= rot_ener;

	  stat_grid[i] += dtemp;
	}
      }
    }

    // normalization
    //
    dtemp = 4. * M_SQRT2 / external_symmetry();
    
    if(mode() == NUMBER)
      //
      dtemp *= 2. / 3.;
    
    stat_grid *= dtemp;
    //
  }//
  // without external rotation
  //
  else {
    //
    for(int p = 0; p < pmax; ++p) {
      //
      for(int l = 0; l < _energy_level[p].size(); ++l) {
	//
	itemp = (int)std::ceil(_energy_level[p][l] / ener_step);

	if(itemp >= stat_grid.size())
	  //
	  continue;
	
	switch(mode()) {
	  //
	case NUMBER:
	  //
	  for(int i = itemp; i < stat_grid.size(); ++i)
	    //
	    stat_grid[i] += 1.;

	  break;

	default:
	  //
	  stat_grid[itemp] += 1.;
	}
      }
    }

    // normalization
    //
    if(mode() != NUMBER)
      //
      stat_grid /= ener_step;

  }// without external rotation

}

double Model::MultiRotor::states (double ener) const
{
  const char funame [] = "Model::MultiRotor::states: ";

  if(ener <= 0.)
    return 0.;

  double dtemp = ener + _ground - _pot_global_min;

  double res;
  if(dtemp < _cstates.arg_max())
    res = _cstates(dtemp);
  else
    res = _cstates.fun_max() * std::pow(dtemp / _cstates.arg_max(), _cstates_pow);

  if(_qfactor.size() && ener < _qfactor.arg_max())
    res *= _qfactor(ener);
  
  return res;
}
  
double Model::MultiRotor::ground () const 
{
  return _ground;
}

double Model::MultiRotor::weight (double temperature) const
{
  double cw, qw;
  get_semiclassical_weight(temperature, cw, qw);
  return qw;
}

void Model::MultiRotor::rotational_energy_levels () const
{
  const char funame [] = "Model::MultiRotor::rotational_energy_levels: ";

  if(_level_ener_max <= 0.) 
    //
    return;

  IO::Marker funame_marker(funame);

  double dtemp;
  int    itemp;
  bool   btemp;

  const Lapack::complex imaginary_unit(0., 1.);

  std::vector<int>    ivec(internal_size());
  std::vector<double> dvec(internal_size());

  IO::log << IO::log_offset << "(external) rotational constants for original geometry [1/cm]:";

  for(int i = 0; i < _rotational_constant.size(); ++i)
    //
    IO::log << "   " << _rotational_constant[i] / Phys_const::incm;

  IO::log << "\n";

  itemp = (int)std::sqrt(_level_ener_max / _rotational_constant.back());

  IO::log << IO::log_offset << "estimated maximum angular momentum needed  = " << itemp << "\n";

  int amom_max = itemp < _amom_max ? itemp : _amom_max;

  std::vector<Lapack::Vector> eigenvalue(amom_max);

  // angular momentum cycle
  //
  for(int amom = 0; amom < amom_max; ++amom) {
    //
    IO::Marker work_marker("fixed angular momentum calculation");

    const int multiplicity = 2 * amom + 1;

    IO::log << IO::log_offset << "J = " << amom << "\n";

    // internal rotation quantum  state dimensions
    //
    const double ener_max = _level_ener_max - double(amom * amom + amom) * _rotational_constant.back();
    
    for(int r = 0; r < internal_size(); ++r) {
      //
      itemp = int(std::sqrt(ener_max / _mobility_parameter[r])) / symmetry(r) * 2 + 1; 

      IO::log << IO::log_offset << r 
	      << "-th internal rotation: suggested optimal internal state dimension = " 
	      <<  itemp << "\n"; 

      if(_internal_rotation[r].quantum_size_max() && itemp > _internal_rotation[r].quantum_size_max()) {
	//
	ivec[r] = _internal_rotation[r].quantum_size_max();
      }
      else if(itemp < _internal_rotation[r].quantum_size_min()) {
	//
	ivec[r] = _internal_rotation[r].quantum_size_min();
      }
      else
	//
	ivec[r] = itemp;
    }

    IO::log << IO::log_offset << "internal phase space dimensions:";

    for(int r = 0; r < internal_size(); ++r)
      //
      IO::log << "  " << ivec[r];

    IO::log << "\n";

    MultiIndexConvert internal_state_index(ivec);

    itemp = multiplicity * internal_state_index.size();

    IO::log << IO::log_offset << "internal   size  = " << internal_state_index.size() << "\n"
	    << IO::log_offset << "multiplicity     = " << multiplicity                << "\n"
	    << IO::log_offset << "hamiltonian size = " << itemp                       << std::endl;

    Lapack::HermitianMatrix ham(itemp);

    // setting hamiltonian
    //
    if(1) {
      //
      IO::Marker work_marker("setting hamiltonian", IO::Marker::ONE_LINE);

      ham = 0.;

#pragma omp parallel for default(shared) private(itemp, dtemp, btemp) schedule(dynamic, 1)

      for(int ml = 0; ml < internal_state_index.size(); ++ml) {
	//
	const std::vector<int> mv = internal_state_index(ml);

	std::vector<int> ivec(internal_size());

	for(int i = 0; i < internal_size(); ++i)
	  //
	  ivec[i] = (mv[i] - internal_state_index.size(i) / 2) * symmetry(i);

	const std::vector<int> mw = ivec;

	for(int nl = 0; nl < internal_state_index.size(); ++nl) {

	  const std::vector<int> nv = internal_state_index(nl);

	  for(int i = 0; i < internal_size(); ++i)
	    //
	    ivec[i] = (nv[i] - internal_state_index.size(i) / 2) * symmetry(i);

	  const std::vector<int> nw = ivec;

	  btemp = false;

	  int ifac = 1;

	  for(int i = 0; i < internal_size(); ++i) {
	    //
	    itemp = nv[i] - mv[i];

	    if(itemp > _mass_index.size(i) / 2 || -itemp > _mass_index.size(i) / 2) {
	      //
	      btemp = true;

	      break;
	    }

	    if(itemp < 0)
	      //
	      itemp += _mass_index.size(i);

	    ivec[i] = itemp;

	    if(_mass_index.size(i) == 2 * itemp)
	      //
	      ifac *= 2;
	  }

	  if(btemp)
	    //
	    continue;

	  itemp = _mass_index(ivec);	  

	  std::map<int, Lapack::ComplexMatrix>::const_iterator imp = _internal_mobility_fourier.find(itemp);
	  std::map<int, Lapack::ComplexMatrix>::const_iterator emp = _external_mobility_fourier.find(itemp);
	  std::map<int, Lapack::ComplexMatrix>::const_iterator ccp = _coriolis_coupling_fourier.find(itemp);
	  
	  for(int mx = 0; mx < multiplicity; ++mx) {
	    //
	    itemp = mx + 3;

	    const int nx_max = itemp < multiplicity ? itemp : multiplicity;

	    const int mproj = mx - amom;

	    for(int nx = mx; nx < nx_max; ++nx) {
	      //
	      const int mtot = ml + mx * internal_state_index.size();

	      const int ntot = nl + nx * internal_state_index.size();
	
	      if(ntot < mtot)
		//
		continue;

	      Lapack::complex matel = 0.;// matrix element

	      // internal mobility
	      //
	      if(imp != _internal_mobility_fourier.end() && nx == mx) {
		//
		// p_i * p_j term
		//
		for(int i = 0; i < internal_size(); ++i)
		  //
		  for(int j = 0; j < internal_size(); ++j)
		    //
		    matel += double(mw[i] * nw[j]) * imp->second(i, j);
	      }

	      if(!amom) {
		//
		ham(mtot, ntot) = matel / (double)ifac;

		continue;
	      }

	      // external mobility
	      //
	      if(emp != _external_mobility_fourier.end()) {
		//
		switch(nx - mx) {
		  //
		case 0:
		  //
		  dtemp = double(amom * amom + amom - mproj * mproj) / 2.;
		  
		  // M_z * M_z term
		  //
		  matel += double(mproj * mproj)
		    //
		    * emp->second(2, 2);

		  // M_x * M_x term
		  //
		  matel += dtemp 
		    //
		    * emp->second(0, 0);

		  // M_y * M_y term
		  //
		  matel += dtemp 
		    //
		    * emp->second(1, 1);

		  break;

		case 1:
		  //
		  dtemp = std::sqrt(double((amom + mproj + 1) * (amom - mproj))) / 2.;

		  // M_x * M_z term
		  //
		  matel += dtemp * double(2 * mproj + 1) 
		    //
		    * emp->second(0, 2);

		  // M_y * M_z term
		  //
		  matel += imaginary_unit * dtemp * double(2 * mproj + 1)
		    //
		    * emp->second(1, 2);

		  break;

		case 2:
		  //
		  dtemp = std::sqrt(double((amom + mproj + 1) * (amom - mproj) * (amom + mproj + 2) *
					   //
					   (amom - mproj - 1))) / 4.;

		  // M_x * M_x term
		  //
		  matel += dtemp 
		    //
		    * emp->second(0, 0);

		  // M_y * M_y term
		  //
		  matel -= dtemp 
		    //
		    * emp->second(1,1);

		  // M_x * M_y term
		  //
		  matel += imaginary_unit * 2. * dtemp 
		    //
		    * emp->second(0, 1);

		  break;
		    
		}// external mobility
		//
	      }//

	      // coriolis coupling
	      //
	      if(ccp != _coriolis_coupling_fourier.end()) {
		//
		for(int i = 0; i < internal_size(); ++i) {
		  //
		  switch(nx - mx) {
		    //
		  case 0:
		    //
		    // M_z * p_i term
		    //
		    matel -= double(mproj * (mw[i] + nw[i])) 
		      //
		      * ccp->second(2, i);
		    
		    break;

		  case 1:
		    //
		    dtemp = std::sqrt(double((amom + mproj + 1) * (amom - mproj))) / 2.;

		    // M_x * p_i term
		    //
		    matel -= dtemp * double(mw[i] + nw[i])
		      //
		      * ccp->second(0, i);
		  
		    // M_y * p_i term
		    //
		    matel -= imaginary_unit * dtemp * double(mw[i] + nw[i])
		      //
		      * ccp->second(1, i);

		    break;
		    //
		  }// coriolis coupling
		}
	      }

	      ham(mtot, ntot) = matel / (double)ifac;
	      //
	    }// nx cycle
	    //
	  }// mx cycle
	  //
	}// nl cycle
	//
      }// ml cycle

      // potential contribution

#pragma omp parallel for default(shared) private(itemp, dtemp, btemp) schedule(dynamic, 1)

      for(int ml = 0; ml < internal_state_index.size(); ++ml) {
	//
	const std::vector<int> mv = internal_state_index(ml);

	std::vector<int> ivec(internal_size());
	//
	for(int nl = ml; nl < internal_state_index.size(); ++nl) {
	  //
	  const std::vector<int> nv = internal_state_index(nl);

	  int ifac = 1;

	  btemp = false;
	  //
	  for(int i = 0; i < internal_size(); ++i) {
	    //
	    itemp = nv[i] - mv[i];

	    if(itemp > _pot_index.size(i) / 2 || -itemp > _pot_index.size(i) / 2) {
	      //
	      btemp = true;

	      break;
	    }
	    if(itemp < 0)
	      //
	      itemp += _pot_index.size(i);

	    if(_pot_index.size(i) == 2 * itemp)
	      //
	      ifac *= 2;

	    ivec[i] = itemp;
	  }

	  if(btemp)
	    continue;
	    
	  itemp = _pot_index(ivec);

	  std::map<int, Lapack::complex>::const_iterator pp = _pot_complex_fourier.find(itemp);

	  if(pp != _pot_complex_fourier.end()) {
	    //
	    for(int mx = 0; mx < multiplicity; ++mx) {
	      //
	      const int mtot = ml + mx * internal_state_index.size();

	      const int ntot = nl + mx * internal_state_index.size();

	      ham(mtot, ntot) += pp->second / (double)ifac;
	      //
	    }// mx cycle
	    //
	  }//
	  //
	}//nl cycle
	//
      }//ml cycle
      //
    }//

    Lapack::HermitianMatrix ctf_mat;

    if(_with_ctf) {// basis scalar product matrix
      //
      IO::Marker work_marker("setting orthogonality matrix for basis states", IO::Marker::ONE_LINE);

      ctf_mat.resize(multiplicity * internal_state_index.size());

      ctf_mat = 0.;

      // potential contribution

#pragma omp parallel for default(shared) private(itemp, dtemp, btemp) schedule(dynamic, 1)

      for(int ml = 0; ml < internal_state_index.size(); ++ml) {
	//
	const std::vector<int> mv = internal_state_index(ml);

	std::vector<int> ivec(internal_size());
	//
	for(int nl = ml; nl < internal_state_index.size(); ++nl) {
	  //
	  const std::vector<int> nv = internal_state_index(nl);

	  int ifac = 1;

	  btemp = false;
	  //
	  for(int i = 0; i < internal_size(); ++i) {
	    //
	    itemp = nv[i] - mv[i];

	    if(itemp > _mass_index.size(i) / 2 || -itemp > _mass_index.size(i) / 2) {
	      //
	      btemp = true;

	      break;
	    }
	    if(itemp < 0)
	      //
	      itemp += _mass_index.size(i);

	    if(_mass_index.size(i) == 2 * itemp)
	      //
	      ifac *= 2;

	    ivec[i] = itemp;
	  }

	  if(btemp)
	    continue;
	    
	  itemp = _mass_index(ivec);

	  std::map<int, Lapack::complex>::const_iterator pp = _ctf_complex_fourier.find(itemp);

	  if(pp != _ctf_complex_fourier.end()) {
	    //
	    for(int mx = 0; mx < multiplicity; ++mx) {
	      //
	      const int mtot = ml + mx * internal_state_index.size();

	      const int ntot = nl + mx * internal_state_index.size();

	      ctf_mat(mtot, ntot) += pp->second / (double)ifac;
	      //
	    }// mx cycle
	    //
	  }//
	  //
	}// nl cycle
	//
      }// ml cycle
      //
    }// basis scalar product matrix

    if(_with_ctf) {
      //
      IO::Marker work_marker("diagonalizing hamiltonian and ctf matrix simultaneously", IO::Marker::ONE_LINE);

      eigenvalue[amom] = Lapack::diagonalize(ham, ctf_mat);
    }
    else {
      //
      IO::Marker work_marker("diagonalizing hamiltonian", IO::Marker::ONE_LINE);

      eigenvalue[amom] = ham.eigenvalues();
    }//
    //
  }//amom cycle

  // rotational energy levels
  //
  std::vector<Lapack::Vector> roten(amom_max);

  roten[0].resize(1);
  
  roten[0][0] = 0.;
  
  for(int j = 1; j < amom_max; ++j) {
    //
    const int multiplicity = 2 * j + 1;

    // setting the hamiltonian
    //
    Lapack::SymmetricMatrix ham(multiplicity);

    ham = 0.;

    for(int mi = 0; mi < multiplicity; ++mi) {
      //
      const int m = mi - j;

      ham(mi, mi) = _rotational_constant[2] * double(m * m)
	+ (_rotational_constant[0] + _rotational_constant[1]) * double(j * j + j - m * m) / 2.;

      if(mi > 1)
	//
	ham(mi - 2, mi) = (_rotational_constant[0] - _rotational_constant[1])
	  * std::sqrt((j + m) * (j - m + 1) * (j + m - 1) * (j - m + 2)) / 4.;
    }

    // diagonalization
    //
    roten[j] = ham.eigenvalues();
  }
  
  // output
  //
  const int lmax = 2 * amom_max - 1;

  const int nmax = 10.;

  IO::log << "\n" << IO::log_offset << "rotational energy levels [1/cm] for original geometry:\n\n";

  IO::log << IO::log_offset << std::setw(3) << "J";

  IO::log << std::setw(13) << "ref";

  for(int j = 0; j < amom_max; ++j)
    //
    IO::log << std::setw(13) << j;

  IO::log << "\n";

  for(int l = 0; l < lmax; ++l) {
    //
    IO::log << IO::log_offset << std::setw(3) << " " << std::setw(13);

    if(!l) {
      //
      IO::log << "-" << std::setw(13) << "0";
    }
    else
      //
      IO::log << " " << std::setw(13) << " ";
    
    
    for(int j = 1; j < amom_max; ++j) {
      //
      IO::log << std::setw(13);
      
      if(l < 2 * j + 1) {
	//
	IO::log << roten[j][l] / Phys_const::incm;
      }
      else
	//
	IO::log << " ";
    }

    IO::log << "\n";
  }
  
  IO::log << "\n" << IO::log_offset << "rovibrational energy levels [1/cm]:\n\n";

  IO::log << IO::log_offset << std::setw(3) << "N\\J";

  IO::log << std::setw(13) << "ref";

  for(int j = 0; j < amom_max; ++j)
    //
    IO::log << std::setw(13) << j;

  IO::log << "\n";

  for(int n = 0; n < nmax; ++n) {
    //
    for(int l = 0; l < lmax; ++l) {
      //
      IO::log << IO::log_offset << std::setw(3);

      if(!l) {
	//
	IO::log << n << std::setw(13);

	if(_energy_level.size()) {
	  //
	  if(n < _energy_level[0].size()) {
	    //
	    IO::log << _energy_level[0][n] / Phys_const::incm;
	  }
	  else
	    //
	    IO::log << "***";
	}
	else if(n < eigenvalue[0].size()) {
	  //
	  IO::log << (eigenvalue[0][n] - _ground) / Phys_const::incm;
	}
	else
	  //
	  IO::log << "***";
      }
      else
	//
	IO::log << "" << std::setw(13) << "";
	
      for(int j = 0; j < amom_max; ++j) {
	//
	const int multiplicity = 2 * j + 1;
	
	const int lindex = multiplicity * n + l;

	IO::log << std::setw(13);

	if(l < multiplicity) {
	  //
	  if(!j && _energy_level.size()) {
	    //
	    if(n < _energy_level[0].size() && lindex < eigenvalue[j].size()) {
	      //
	      dtemp = (eigenvalue[j][lindex] - _energy_level[0][n] - _ground) / Phys_const::incm;

	      if(dtemp < 1.e-7 && dtemp > -1.e-7) {
		//
		IO::log << "< 1.e-7";
	      }
	      else
		//
		IO::log << dtemp;
	    }
	    else
	      //
	      IO::log << "***";
	  }
	  else {
	    //
	    if(n < eigenvalue[0].size() && lindex < eigenvalue[j].size()) {
	      //
	      IO::log << (eigenvalue[j][lindex] - eigenvalue[0][n]) / Phys_const::incm;
	    }
	    else
	      // 
	      IO::log << "***";
	  }
	}
	else
	  //
	  IO::log << "";
      }

      IO::log << "\n";
    }
  }

  IO::log << "\n";
}

/********************************************************************************************
 *********** ABSTRACT CLASS REPRESENTING WELL, BARRIER, AND BIMOLECULAR FRAGMENT ************
 ********************************************************************************************/

Model::Species::Species (IO::KeyBufferStream& from, const std::string& n, int m) 
  : _name(n), _mode(m), _mass(-1.), _print_step(-1.), _ground(0.)
{
  const char funame [] = "Model::Species::Species: ";

  IO::Marker funame_marker(funame);
  
  KeyGroup SpeciesModel;

  Key ang_geom_key("Geometry[angstrom]");
  Key bor_geom_key("Geometry[au]"      );
  Key     mass_key("Mass[amu]"         );
  Key     stoi_key("Stoichiometry"     );
  Key    print_key("Print[1/cm]"       );

  int    itemp;
  double dtemp;

  std::string token, stemp, word, comment, line;

  while(from >> token) {
    // geometry
    if(ang_geom_key == token || bor_geom_key == token) {
      if(_atom.size() || _mass > 0.) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      if(ang_geom_key == token)
	read_geometry(from, _atom, ANGSTROM);
      if(bor_geom_key == token)
	read_geometry(from, _atom, BOHR);
    }
    // print parameters
    else if(print_key == token) {
      IO::LineInput lin(from);
      if(!(lin >> _print_min >> _print_max >> _print_step)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(_print_step <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // mass
    else if(mass_key == token) {
      if(_atom.size() || _mass > 0.) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      if(!(from >> _mass)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      };
      std::getline(from, comment);

      if(_mass <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      };
	
      _mass *= Phys_const::amu;

    }
    // stoichiometry
    else if(stoi_key == token) {

      token += "(stoichiometry format example: C1H3O2): ";

      if(_atom.size() || _mass > 0.) {
	std::cerr << funame << token << "already defined\n";
	throw Error::Init();
      }

      if(!(from >> word)) {
	std::cerr << funame << token << "corrupted\n";
	throw Error::Input();
      };
      std::getline(from, comment);

      std::string name, number;
      for(int i = 0; i < word.size(); ++i) {
	char c = word[i];
	if(std::isdigit(c))
	  number.push_back(c);
	else if(std::isalpha(c)) {
	  if(number.size()) {
	    itemp = std::atoi(number.c_str());
	    if(itemp <= 0) {
	      std::cerr << funame << token << "wrong number of atoms: " << itemp << "\n";
	      throw Error::Range();
	    }
	    _mass += AtomBase(name).mass() * double(itemp);
	    name.clear();
	    number.clear();
	  }
	  name.push_back(c);
	}
	else {
	  std::cerr << funame << token << "unrecognized character: " << c << "\n";
	  throw Error::Input();
	}
      }

      if(!name.size() || !number.size()) {
	std::cerr << funame << token << "wrong format: " << word << "\n";
	throw Error::Form();
      }

      itemp = std::atoi(number.c_str());
      if(itemp <= 0) {
	std::cerr << funame << token << "wrong number of atoms: " << itemp << "\n";
	throw Error::Range();
      }
      _mass += AtomBase(name).mass() * double(itemp);	
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      from.put_back(token);

      break;
    }
  }

  /******************************************* Checking *************************************************/

  if(mode() != DENSITY && mode() != NUMBER && mode() !=  NOSTATES) {
    std::cerr << funame << "wrong mode\n";
    throw Error::Logic();
  }

  // stream state
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_atom.size()) {
    _mass = 0.;
    for(int a = 0; a < _atom.size(); ++a)
      _mass += _atom[a].mass();
  }

}// Species

Model::Species::~Species ()
{
  //std::cout << "Model::Species destroyed\n";
}

void Model::Species::_print () const
{
  if(_print_step < 0. || mode() == NOSTATES)
    return;

  std::ofstream sout;

  if(mode() == NUMBER)
    sout.open((name() + "_nos.out").c_str());
  else
    sout.open((name() + "_dos.out").c_str());
  
  for(double ener = _print_min; ener <= _print_max; ener += _print_step) {
    sout << std::setw(15) << ener << std::setw(15);
    
    if(mode() == NUMBER)
      sout << states(ener * Phys_const::incm);
    else
      sout << states(ener * Phys_const::incm) * Phys_const::incm;
    
    sout << "\n";
  }
}

double Model::Species::mass () const 
{
  const char funame [] = "Model::Species::mass: ";

  if(_mass <= 0.) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  return _mass;
}

double Model::Species::tunnel_weight (double) const { return 1.; }

// radiational transitions
double Model::Species::oscillator_frequency (int) const
{
  const char funame [] = "Model::Species::oscillator_frequency: ";

  std::cerr << funame << "should not be here\n";
  throw Error::Logic();
}

double  Model::Species::infrared_intensity (double, int) const
{
  const char funame [] = "Model::Species::infrared_intensity: ";

  std::cerr << funame << "should not be here\n";
  throw Error::Logic();
}

int Model::Species::oscillator_size () const
{
  return 0;
}

/*******************************************************************************************
 ************************* READ STATES FROM THE FILE AND INTERPOLATE ***********************
 *******************************************************************************************/
  
Model::ReadSpecies::ReadSpecies (std::istream& from, const std::string& n, int m) 
  : Species(n, m), _etol(1.e-10), _dtol(1.)
{
  const char funame [] = "Model::ReadSpecies::ReadSpecies:";

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE);
  //IO::Marker funame_marker(funame);

  KeyGroup ReadModel;

  Key  mass_key("Mass[amu]"  );
  Key  unit_key("EnergyUnits");
  Key  file_key("File"       );

  Key  incm_ener_key("GroundEnergy[1/cm]");
  Key  kcal_ener_key("GroundEnergy[kcal/mol]");
  Key    kj_ener_key("GroundEnergy[kJ/mol]");
  Key    ev_ener_key("GroundEnergy[eV]");
  Key    au_ener_key("GroundEnergy[au]");
  Key   zero_nos_key("GroundStateThreshold");

  double energy_unit  = Phys_const::incm;
  double ground_shift = 0.;

  if(mode() == NOSTATES) {
    std::cerr << funame << "wrong calculation mode\n";
    throw Error::Logic();
  }
  
  int    itemp;
  double dtemp;

  std::map<double, double> read_states;

  std::string token, comment, name, stemp;

  while(from >> token) {
    // end input
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // energy shift
    else if(incm_ener_key == token || 
	    kcal_ener_key == token || 
	    kj_ener_key == token || 
	    ev_ener_key == token ||
	    au_ener_key == token) {

      if(!(from >> ground_shift)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      
      if(incm_ener_key == token) 
	ground_shift *= Phys_const::incm;
      else if(kcal_ener_key == token)
	ground_shift *= Phys_const::kcal;
      else if(kj_ener_key == token)
	ground_shift *= Phys_const::kjoul;
      else if(ev_ener_key == token)
	ground_shift *= Phys_const::ev;
    }
    // units
    else if(unit_key == token) {
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      
      if(stemp == "kcal/mol")
	energy_unit = Phys_const::kcal;
      else if(stemp == "1/cm")
	energy_unit = Phys_const::incm;
      else if(stemp == "eV")
	energy_unit = Phys_const::ev;
      else if(stemp == "au")
	energy_unit = 1.;
      else if(stemp == "kJ/mol")
	energy_unit = Phys_const::kjoul;
      else {
	std::cerr << funame << token << ": unknown energy unit: " << stemp << ": available energy units: kcal/mol, 1/cm, eV, au, kJ/mol\n";
	throw Error::Range();
      }
    }
    // ground state threshold
    else if(zero_nos_key == token) {
      if(!(from >> _dtol)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(_dtol <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // mass
    else if(mass_key == token) {
      if(_mass > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      if(!(from >> _mass)) {
	std::cerr << funame << token << ": unreadable\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      
      if(_mass <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      _mass *= Phys_const::amu;
    }
    // reading density/number of states from file
    else if(file_key == token) {
      if(!(from >> name)) {
	std::cerr << funame << token << ": corrupted\n";
	IO::log << std::flush;
	throw Error::Input();
      }
      std::getline(from, comment);
      std::ifstream file(name.c_str());
      if(!file) {
	std::cerr << funame << "cannot open file " << name << " for reading\n";
	IO::log << std::flush;
	throw Error::File();
      }

      double eval, dval;
      while(file >> eval) {
	file >> dval;        // density of states
	if(!file) {
	  std::cerr << funame << token << ": reading number/density failed\n";
	  throw Error::Input();
	}
	read_states[eval] = dval;
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

  if(mode() == DENSITY)
    for(std::map<double, double>::iterator i = read_states.begin(); i != read_states.end(); ++i)
      i->second /= energy_unit;

  // find zero density energy
  typedef std::map<double, double>::const_reverse_iterator It;
  It izero;
  for( izero = read_states.rbegin(), itemp = 1; izero != read_states.rend(); ++izero, ++itemp)
    if(izero->second <= _dtol)
      break;

  if(itemp < 3) {
    std::cerr << funame << "not enough data\n";
    throw Error::Init();
  }

  _ener.resize(itemp);
  _states.resize(itemp);

  _ener[0] = 0.;
  _states[0]  = 0.;
  
  --itemp;
  for(It it = read_states.rbegin(); it != izero; ++it, --itemp) {
    _ener[itemp]    = it->first * energy_unit;
    _states[itemp]  = it->second;
  }

  //IO::log << IO::log_offset << "original # of points = " << read_states.size() << "\n";

  if(izero != read_states.rend())
    _ground = izero->first * energy_unit;
  else
    _ground = 2. * _ener[1] - _ener[2];  

  for(int i = 1; i < _ener.size(); ++i)
    _ener[i] -= _ground;

  _ground += ground_shift;

  Array<double> x(_ener.size() - 1);
  Array<double> y(_ener.size() - 1);

  for(int i = 1; i < _ener.size(); ++i) {
    itemp = i - 1;
    x[itemp] = std::log(_ener[i]);
    y[itemp] = std::log( _states[i]);
  }

  _spline.init(x, y, x.size());

  _emin = _ener[1]     * (1. + _etol);
  _emax = _ener.back() * (1. - _etol);

  _nmin = (y[1] - y[0]) / (x[1] - x[0]);
  _amin = std::exp(y[0] - x[0] * _nmin);

  int l1 = x.size() - 1;
  int l2 = x.size() - 2;
  _nmax = (y[l1] - y[l2]) / (x[l1] - x[l2]);
  _amax = std::exp(y[l1] - x[l1] * _nmax);

  _print();
}// Read Model

Model::ReadSpecies::~ReadSpecies ()
{
  //std::cout << "Model::ReadSpecies destroyed\n";
}

double Model::ReadSpecies::states (double en) const
{
  en -= _ground;

  if(en <= 0.)
    return 0.;

  if(en <= _emin)
    return _amin * std::pow(en, _nmin);

  if(en >= _emax)
    return _amax * std::pow(en, _nmax);

  return std::exp(_spline(std::log(en)));
}

double Model::ReadSpecies::weight (double temperature) const
{
  const char funame [] = "Model::ReadSpecies::weight: ";

  Array<double> term(_states.size());
  term[0] = 0.;
  for(int i = 1; i < _states.size(); ++i)
    term[i] = _states[i] / std::exp(_ener[i] / temperature);

  int info;
  double res;
  davint_(_ener, term, term.size(), _ener.front(), _ener.back(), res, info); 
  if (info != 1) {
    std::cerr << funame  << "davint integration error\n";
    throw Error::Logic();
  }

  if(mode() == NUMBER)
    res /= temperature;

  return res;
}

/***********************************************************************************************************
 ************************************* CRUDE MONTE-CARLO SAMPLING ******************************************
 ***********************************************************************************************************/

double Model::IntMod::increment = 1.e-2;

Model::IntMod::IntMod (int molec_size, IO::KeyBufferStream& from)
{
  const char funame [] = "Model::IntMod::IntMod: ";

  IO::Marker funame_marker(funame);
  
  int    itemp;
  double dtemp;

  if(molec_size < 2) {
    //
    std::cerr << funame << "molecule size is out of range: " << molec_size << "\n";

    throw Error::Range();
  }
  
  KeyGroup InternalModeDefinition;

  Key  atom_key("AtomIndices");
  
  std::string token, comment, stemp;

  while(from >> token) {
    //
    // atomic indices
    //
    if(atom_key == token) {
      //
      if(_atoms.size()) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Input();
      }
      
      // indices should be on the same line with the keyword, start from 1 (not from 0)
      //      
      IO::LineInput lin(from);

      while(lin >> itemp) {
	//
	if(itemp < 1 || itemp > molec_size) {
	  //
	  std::cerr << funame << token << ": atomic index is out of range: " << itemp << "\n";

	  throw Error::Range();
	}

	if(!_pool.insert(--itemp).second) {
	  //
	  std::cerr << funame << token << ": identical atomic indices: " << itemp + 1 << "\n";

	  throw Error::Input();
	}

	_atoms.push_back(itemp);
      }

      if(_atoms.size() < 2 || _atoms.size() > 4) {
	//
	std::cerr << funame << token << ": number of atomic indices is out of range: " << _atoms.size() << "\n";

	throw Error::Range();
      }
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      from.put_back(token);

      break;
    }
  }

  // checking
  //
  if(!_atoms.size()) {
    //
    std::cerr << funame << "not defined\n";

    throw Error::Init();
  }
}

double Model::IntMod::evaluate (Lapack::Vector          cart_pos, // atomic cartesian coordinates
				const std::vector<int>& sign      // derivative signature
				) const
{
  const char funame [] = "Model::IntMod::evaluate: ";
  
  double dtemp;
  
  int    itemp;
  
  // check if there are atoms in the derivative signature absent in the definition
  //
 for(int i = 0; i < sign.size(); ++i)
    //
    if(_pool.find(sign[i] / 3) == _pool.end())
      //
      return 0.;
      
  if(sign.size()) {
    //
    std::vector<int> new_sign = sign;

    new_sign.pop_back();

    Lapack::Vector new_pos = cart_pos.copy();

    new_pos[sign.back()] -= increment;

    dtemp = evaluate(new_pos, new_sign);

    new_pos[sign.back()] += 2. * increment;

    dtemp = evaluate(new_pos, new_sign) - dtemp;

    // correct for 360 jump for dihedral angle
    //
    if(!new_sign.size() && type() == DIHEDRAL) {
      //
      if(dtemp > 180.)
	//
	dtemp -= 360.;

      if(dtemp < -180.)
	//
	dtemp += 360.;
    }
    
    return dtemp / 2. / increment;
  }

  std::vector<D3::Vector> atom_pos(_atoms.size());

  for(int a = 0; a < _atoms.size(); ++a)
    //
    for(int ai = 0; ai < 3; ++ai)
      //
      atom_pos[a][ai] = cart_pos[_atoms[a] * 3 + ai];

  D3::Vector v1, v2, n;
  
  switch(type()) {
    //
    // interatomic distance in angstrom
    //
  case DISTANCE:
    //
    return (atom_pos[1] - atom_pos[0]).vlength() / Phys_const::angstrom;

    // plane angle in degrees
    //
  case ANGLE:
    //
    v1 = atom_pos[0] - atom_pos[1];
    
    v2 = atom_pos[2] - atom_pos[1];
    
    dtemp = vdot(v1, v2) / v1.vlength() / v2.vlength();
    
    if(dtemp < -1.)
      //
      dtemp = -1.;
    
    if(dtemp > 1.)
      //
      dtemp = 1.;

    return std::acos(dtemp) * 180. / M_PI;

    // dihedral angle in degrees
    //
  case DIHEDRAL:
    //
    n  = atom_pos[2] - atom_pos[1];
    
    v1 = atom_pos[0] - atom_pos[1];
    
    v2 = atom_pos[3] - atom_pos[2];

    v1.orthogonalize(n);
    
    v2.orthogonalize(n);
    
    dtemp = v1.vlength() * v2.vlength();
    
    if(dtemp < 1.e-8)
      //
      return 0.;

    dtemp = vdot(v1, v2) / dtemp;
    
    if(dtemp < -1.)
      //
      dtemp = -1.;
  
    if(dtemp > 1.)
      //
      dtemp = 1.;

    dtemp = std::acos(dtemp) * 180. / M_PI;

    if(volume(n, v1, v2) > 0.) {
      //
      return dtemp;
    }
    else
      //
      return 360. - dtemp;

    // wrong case
    //
  default:

    std::cerr << funame << "should not be here\n";

    throw Error::Logic();
  }
}

Model::Fluxional::Fluxional (int molec_size, IO::KeyBufferStream& from) : IntMod(molec_size, from), _span(-1.)
{
  const char funame [] = "Model::Fluxional::Fluxional: ";

  IO::Marker funame_marker(funame);
  
  int    itemp;
  double dtemp;

  KeyGroup FluxionalModeDefinition;

  Key  span_key("Span");
  
  std::string token, comment, stemp;

  while(from >> token) {
    //
    // end input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // span
    //
    else if(span_key == token) {
      //
      if(_span > 0.) {
	//
	std::cerr << funame << token << ": already defined\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);

      if(!(lin >> _span)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      if(_span <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << _span << "\n";
	
	throw Error::Range();
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

  if(_span < 0.) {
    //
    switch(type()) {
      //
    case DISTANCE:
      //
      IO::log << IO::log_offset << "ERROR: for interatomic distances span should be explicitely defined\n";

      std::cerr << funame << "for interatomic distances span should be explicitely defined\n";

      throw Error::Init();

    case ANGLE:
      //
      IO::log << IO::log_offset << "WARNING: span is not defined, default span of 180 for plane angle is assumed" << std::endl;

      _span = 180.;

      break;

    case DIHEDRAL:
      //
      IO::log << IO::log_offset << "WARNING: span is not defined, default span of 360 for dihedral angle is assumed" << std::endl;

      _span = 360.;

      break;

    default:

      std::cerr << funame << "should not be here\n";

      throw Error::Logic();
    }
  }
}

Model::Constrain::Constrain (int molec_size, IO::KeyBufferStream& from) : IntMod(molec_size, from), _value(-1.)
{
  const char funame [] = "Model::Constrain::Constrain: ";

  IO::Marker funame_marker(funame);
  
  int    itemp;
  double dtemp;

  KeyGroup ConstrainDefinition;

  Key  val_key("Value");
  
  std::string token, comment, stemp;

  while(from >> token) {
    //
    // end input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // value
    //
    else if(val_key == token) {
      //
      if(_value >= 0.) {
	//
	std::cerr << funame << token << ": already defined\n";
	
	throw Error::Init();
      }

      IO::LineInput lin(from);

      if(!(lin >> _value)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }

      switch(type()) {
	//
	// bond distance
	//
      case DISTANCE:
	//
	if(_value <= 0.) {
	  //
	  std::cerr << funame << token << ": out of range: " << _value << "\n";
	
	  throw Error::Range();
	}

	break;
	
	// bond angle
	//
      case ANGLE:
	//
	if(_value < 0. || _value > 180.) {
	  //
	  std::cerr << funame << token << ": out of range: " << _value << "\n";
	
	  throw Error::Range();
	}
	  
	break;

	// dihedral angle
	//
      case DIHEDRAL:
	//
	if(_value < 0. || _value > 360.) {
	  //
	  std::cerr << funame << token << ": out of range: " << _value << "\n";
	
	  throw Error::Range();
	}

	break;

	// wrong case
	//
      default:
	//
	std::cerr << funame << "should not be here\n";

	throw Error::Logic();
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

  if(_value < 0.) { 
    //
    std::cerr << funame << "value not defined\n";

    throw Error::Init();
  }      
}

double Model::MonteCarlo::_RefPot::operator() (const double* pos) const
{
  const char funame [] = "Model::MonteCarlo::_RefPot::operator(): ";

  if(!_pot) {
    //
    std::cerr << funame << "Potential not initialized\n";

    throw Error::Init();
  }

  double ener;
  
  int ifail;
  
  _pot(ifail, ener, pos);

  if(ifail) {
    //
    std::cerr << funame << "reference potential failed with the code: " << ifail << "\n";

    throw PotFail();
  }

  return ener;
}

double Model::MonteCarlo::_RefPot::weight (double temperature) const
{
  const char funame [] = "Model::MonteCarlo::_RefPot::weight: ";

  if(!_weight) {
    //
    std::cerr << funame << "not initialized\n";

    throw Error::Init();
  }

  if(temperature <= 0.) {
    //
    std::cerr << funame << "temperature out of range: " << temperature / Phys_const::kelv << "\n";

    throw Error::Range();
  }

  double res;
  
  _weight(res, temperature);

  return res;
}

void Model::MonteCarlo::_RefPot::init(std::istream& from)
{
  const char funame [] = "Model::MonteCarlo::_RefPot::init: ";

  int    itemp;
  double dtemp;

  std::string token, comment, stemp;

  KeyGroup ReferencePotentialModel;

  Key     lib_key("Library"  );
  Key     pot_key("Potential");
  Key  weight_key("Weight"   );
  
  while(from >> token) {
    //
    // end input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // library
    //
    else if(lib_key == token) {
      //
      if(_lib.isopen()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      IO::LineInput lin(from);

      if(!(lin >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      _lib.open(stemp);
    }
    // potential symbol
    //
    else if(pot_key == token) {
      //
      //
      if(_pot) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!_lib.isopen()) {
	//
	std::cerr << funame << token << ": library should be initialized first\n";

	throw Error::Init();
      }
	
      IO::LineInput lin(from);

      if(!(lin >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      _pot = (refp_t)_lib.member(stemp);
    }
    // statistical symbol
    //
    else if(weight_key == token) {
      //
      //
      if(_weight) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!_lib.isopen()) {
	//
	std::cerr << funame << token << ": library should be initialized first\n";

	throw Error::Init();
      }
	
      IO::LineInput lin(from);

      if(!(lin >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      _weight = (refw_t)_lib.member(stemp);
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
  // some checking
  //
  if(!from) {
    //
    std::cerr << funame << "input stream is currupted\n";

    throw Error::Input();
  }

  if(!_pot) {
    //
    std::cerr << funame << "potential not initialized\n";

    throw Error::Init();
  }

  if(!_weight) {
    //
    std::cerr << funame << "statistical weight not initialized\n";

    throw Error::Init();
  }
}

Model::MonteCarlo::~MonteCarlo() {}

Model::MonteCarlo::MonteCarlo(IO::KeyBufferStream& from, const std::string& n, int m)
  : Species(from, n, m), _symm_fac(1.), _noqf(false), _nohess(false), _nocurv(false),
    _cmshift(false), _ists(false), _ref_tem(-1.)
{
  const char funame [] = "Model::MonteCarlo::MonteCarlo: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  
  double dtemp;

  std::string token, comment, stemp;

  KeyGroup MonteCarloModel;

  Key    atom_key("MoleculeSpecification"        );
  Key    flux_key("FluxionalMode"                );
  Key    incr_key("NumericIncrement[Bohr]"       );
  Key    symm_key("SymmetryFactor"               );
  Key    data_key("DataFile"                     );
  Key    rpot_key("ReferencePotential"           );
  Key    rtem_key("ReferenceTemperature[K]"      );
  Key    noqf_key("NoQuantumCorrection"          );
  Key  nohess_key("NoHessian"                    );
  Key  nocurv_key("NoCurvlinearCorrection"       );
  Key cmshift_key("UseCMShift"                   );
  Key    freq_key("NonFluxionalFrequencies[1/cm]");
  Key      ts_key("TransitionState");
  
  Key incm_refen_key("ReferenceEnergy[1/cm]"          );
  Key kcal_refen_key("ReferenceEnergy[kcal/mol]"      );
  Key   kj_refen_key("ReferenceEnergy[kJ/mol]"        );
  Key   ev_refen_key("ReferenceEnergy[eV]"            );
  Key   au_refen_key("ReferenceEnergy[au]"            );

  Key incm_ground_key("GroundEnergy[1/cm]"          );
  Key kcal_ground_key("GroundEnergy[kcal/mol]"      );
  Key   kj_ground_key("GroundEnergy[kJ/mol]"        );
  Key   ev_ground_key("GroundEnergy[eV]"            );
  Key   au_ground_key("GroundEnergy[au]"            );

  Key incm_elev_key("ElectronicLevels[1/cm]"          );
  Key kcal_elev_key("ElectronicLevels[kcal/mol]"      );
  Key   kj_elev_key("ElectronicLevels[kJ/mol]"        );
  Key   ev_elev_key("ElectronicLevels[eV]"            );
  Key   au_elev_key("ElectronicLevels[au]"            );

  bool isrefen = false;

  bool isground = false;
  
  while(from >> token) {
    //
    // end input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // exclude quantum correction
    //
    else if(noqf_key == token) {
      //
      std::getline(from, comment);

      _noqf = true;
    }
    // no hessian data
    //
    else if(nohess_key == token) {
      //
      std::getline(from, comment);

      _nohess = true;
    }
    // no curvlinear correction to hessian
    //
    else if(nocurv_key == token) {
      //
      std::getline(from, comment);

      _nocurv = true;
    }
    // use center-of-mass shift
    //
    else if(cmshift_key == token) {
      //
      std::getline(from, comment);

      _cmshift = true;
    }
    // transition state calculation
    //
    else if(ts_key == token) {
      //
      std::getline(from, comment);

      _ists = true;
    }
    // electronic energy levels & degeneracies
    //
    else if(incm_elev_key == token || au_elev_key == token ||
	    kcal_elev_key == token || kj_elev_key == token ||
	    ev_elev_key == token) {
      //
      int num;
      
      if(!(from >> num)) {
	std::cerr << funame << token << ": levels number unreadable\n";
	throw Error::Input();
      }

      if(num < 1) {
	std::cerr << funame << token << ": levels number should be positive\n";
	throw Error::Range();
      }

      std::getline(from, comment);

      for(int l = 0; l < num; ++l) {
	//
	IO::LineInput level_input(from);
	
	if(!(level_input >> dtemp >> itemp)) {
	  //
	  std::cerr << funame << token << ": format: energy degeneracy(>=1)\n";

	  throw Error::Input();
	}

	if(incm_elev_key == token)
	  //
	  dtemp *= Phys_const::incm;
	
	if(kcal_elev_key == token)
	  //
	  dtemp *= Phys_const::kcal;
	
	if(kj_elev_key == token)
	  //
	  dtemp *= Phys_const::kjoul;
	
	if(ev_elev_key == token)
	  //
	  dtemp *= Phys_const::ev;

	if(_elevel.find(dtemp) != _elevel.end()) {
	  //
	  std::cerr << funame << token << ": identical energy levels\n";
	  
	  throw Error::Range();
	}

	if(itemp < 1) {
	  //
	  std::cerr << funame << token << ": degeneracy should be positive\n";
	  
	  throw Error::Range();
	}

	_elevel[dtemp] = itemp;
      }
    }
    // non-fluxional modes frequencies
    //
    else if(freq_key == token) {
      //
      if(_nm_freq.size()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!(from >> itemp)) {
	//
	std::cerr << funame << token << ": cannot read frequencies #\n";

	throw Error::Input();
      }

      if(itemp < 1) {
	//
	std::cerr << funame << token << ": frequencies # out of range: " << itemp << "\n";

	throw Error::Range();
      }

      _nm_freq.resize(itemp);

      for(int f = 0; f < _nm_freq.size(); ++f) {
	//
	if(!(from >> dtemp)) {
	  //
	  std::cerr << funame << token << ": cannot read " << f + 1 << "-th frequency\n";

	  throw Error::Input();
	}

	if(dtemp < nm_freq_min) {
	  //
	  std::cerr << funame << token << ": " << f + 1 << "-th non-fluxional mode frequency is below the threshold ("
		    << nm_freq_min << "): " << dtemp << "\n";

	  throw Error::Range();
	}

	_nm_freq[f] = dtemp * Phys_const::incm;
      }
    }
    //  reference energy
    //
    else if(token == incm_refen_key ||
	    token == kcal_refen_key ||
	    token ==   kj_refen_key || 
	    token ==   au_refen_key ||
	    token ==   ev_refen_key) {
      //
      if(isrefen) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Init();
      }
      
      isrefen = true;
      
      if(!(from >> _refen)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      std::getline(from, comment);

      if(incm_refen_key == token) {
	//
	_refen *= Phys_const::incm;
      }
      if(kcal_refen_key == token) {
	//
	_refen *= Phys_const::kcal;
      }
      if(kj_refen_key == token) {
	//
	_refen *= Phys_const::kjoul;
      }
      if(ev_refen_key == token) {
	//
	_refen *= Phys_const::ev;
      }
    }
    //  ground energy
    //
    else if(token == incm_ground_key ||
	    token == kcal_ground_key ||
	    token ==   kj_ground_key || 
	    token ==   au_ground_key ||
	    token ==   ev_ground_key) {
      //
      if(isground) {
	//
	std::cerr << funame << token << ": already initialized\n";
	
	throw Error::Init();
      }
      
      isground = true;
      
      if(!(from >> _ground)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      std::getline(from, comment);

      if(incm_ground_key == token) {
	//
	_ground *= Phys_const::incm;
      }
      if(kcal_ground_key == token) {
	//
	_ground *= Phys_const::kcal;
      }
      if(kj_ground_key == token) {
	//
	_ground *= Phys_const::kjoul;
      }
      if(ev_ground_key == token) {
	//
	_ground *= Phys_const::ev;
      }
    }
    // atomic masses
    //
    else if(atom_key == token) {
      //
      if(_mass_sqrt.size()) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Init();
      }

      IO::LineInput lin(from);

      _mass = 0.;
      
      while(lin >> stemp) {
	//
	// dot position if any
	//
	std::string::size_type dot_pos = stemp.find('.');

	int isot = -1;

	// isotope specification
	//
	if(dot_pos < stemp.size()) {
	  //
	  if(dot_pos == stemp.size() - 1 || dot_pos == 0) {
	    //
	    std::cerr << funame << token << ": dot position out of range: " << stemp << "\n";

	    throw Error::Range();
	  }

	  // isotope number starts at
	  //
	  itemp = dot_pos + 1;
      
	  isot = (int)IO::String(stemp.substr(itemp));

	  if(isot <= 0) {
	    //
	    std::cerr << funame << token << ": isotope number out of range: " << stemp << "\n";

	    throw Error::Range();
	  }

	  stemp = stemp.substr(0, dot_pos);
	}

	AtomBase a;

	if(isot > 0) {
	  //
	  a.set(stemp, isot);
	}
	else
	  //
	  a.set(stemp);

	dtemp = a.mass();

	_mass += dtemp;

	_mass_sqrt.push_back(std::sqrt(dtemp));
      }

      _total_mass_sqrt = std::sqrt(_mass);
    }
    // fluxional mode definition
    //
    else if(flux_key == token) {
      //
      std::getline(from, comment);

      if(!atom_size()) {
	//
	std::cerr << funame << token << ": molecular specification should go first\n";

	throw Error::Init();
      }

      _fluxional.push_back(Fluxional(atom_size(), from));
    }
    // numerical differentiation increment
    //
    else if(incr_key == token) {
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

      IntMod::increment = dtemp;
    }
    // reference potential
    //
    else if(rpot_key == token) {
      //
      std::getline(from, comment);

      if(_ref_pot) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      _ref_pot.init(from);
    }
    // reference temperature
    //
    else if(rtem_key == token) {
      //
      if(_ref_tem > 0.) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> _ref_tem)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(_ref_tem <= 0.) {
	//
	std::cerr << funame << token << ": out of range\n";

	throw Error::Range();
      }

      _ref_tem *= Phys_const::kelv;
    }
    // symmetry factor
    //
    else if(symm_key == token) {
      //
      IO::LineInput lin(from);

      if(!(lin >> dtemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(dtemp <= 0) {
	//
	std::cerr << funame << token << ": out of range: " << dtemp << "\n";

	throw Error::Range();
      }

      _symm_fac = dtemp;
    }
    // data file
    //
    else if(data_key == token) {
      //
      IO::LineInput lin(from);

      if(_data_file.size()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }
      
      if(!(lin >> _data_file)) {
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

  // some checking
  //
  if(!from) {
    //
    std::cerr << funame << "input stream is currupted\n";

    throw Error::Input();
  }

  if(!atom_size()) {
    //
    std::cerr << funame << "no atoms specified\n";

    throw Error::Init();
  }

  if(!_fluxional.size()) {
    //
    std::cerr << funame << "no fluxional modes specified\n";

    throw Error::Init();
  }

  if(_fluxional.size() + 6 > atom_size() * 3) {
    //
    std::cerr << funame << "number of fluxional modes (" << _fluxional.size()
	      << ") is inconsistent with the number of atoms, " << atom_size() << "\n";

    throw Error::Range();
  }

  if(_nohess && _nm_freq.size() + _fluxional.size() + 6 != atom_size() * 3) {
    //
    std::cerr << funame << "non-fluxional modes # (" << _nm_freq.size() << ") and fluxional modes # (" << _fluxional.size()
	      << ") mismatch with the atoms #: " << atom_size() << "\n";

    throw Error::Range();
  }
  
  if(!_data_file.size()) {
    //
    std::cerr << funame << "no data file is provided\n";

    throw Error::Init();
  }

  if(!isrefen && isground) {
    //
    IO::log << IO::log_offset << "WARNING: the reference energy set to the ground state energy: " << _ground  / Phys_const::kcal << " kcal/mol\n";

    _refen = _ground;
  }
  
  if(isrefen && !isground) {
    //
    IO::log << IO::log_offset << "WARNING: the ground state energy set to the reference energy: " << _refen  / Phys_const::kcal << " kcal/mol\n";

    _ground = _refen;
  }
  
  if(!isground && !isrefen) {
    //
    _set_reference_energy();

    _ground = _refen;
    
    IO::log << IO::log_offset << "WARNING: ground and reference energies initialized to the minimal sampling energy (with ZPE): "
	    << _refen / Phys_const::kcal << " kcal/mol\n";
  }

  if(_ref_pot && _ref_tem < 0. || !_ref_pot && _ref_tem > 0.) {
    //
    std::cerr << funame << "the reference potential and the reference temperature should be defined simultaneously\n";

    throw Error::Init();
  }

  if(!_elevel.size())
    //
    _elevel[0.] = 1;
}

double Model::MonteCarlo::states (double ener) const
{
  const char funame [] = "Model::MonteCarlo::states: ";

  std::cerr << funame << "not implemented yet, sorry\n";

  throw Error::General();

  return 0.;
}

double Model::MonteCarlo::weight_with_error (double temperature, double& werr) const
{
  const char funame [] = "Model::MonteCarlo::weight_with_error: ";

  int    itemp;
  
  double dtemp;

  IO::Marker funame_marker(funame);
  
  std::ifstream from(_data_file.c_str());

  if(!from) {
    //
    std::cerr << funame << "cannot open data file " << _data_file << "\n";

    throw Error::File();
  }

  std::string comment;

  std::getline(from, comment);
  
  double res = 0.;

  double variance = 0.;

  double ener;

  // cartesian coordinates
  //
  Lapack::Vector          cart_pos(atom_size() * 3);

  // cartesian gradient
  //
  Lapack::Vector          cart_grad(atom_size() * 3);

  // cartesian force constant matrix
  //
  Lapack::SymmetricMatrix cart_fc(atom_size() * 3);

  int count = 0, samp = 0;
  
  while(_read(from, ener, cart_pos, cart_grad, cart_fc)) {
    //
    ++samp;
    
    IO::log << IO::log_offset << "Sampling " << samp << "\n";
    
    dtemp = _local_weight(ener, cart_pos, cart_grad, cart_fc, temperature);

    if(dtemp < 0.) {
      //
      std::cerr << funame << samp << "-th sampling failed\n";

      continue;
    }
    
    ++count;
    
    res += dtemp;

    variance += dtemp * dtemp;
  }

  if(!count) {
    //
    std::cerr << funame << "no data\n";

    throw Error::Input();
  }
  
  // normalize
  //
  res /= (double)count;

  variance /= (double)count * res * res;

  variance -= 1.;

  variance /= (double)count;

  if(variance < 0.)
    //
    variance = 0.;

  werr = std::sqrt(variance) * 100.;
  
  // span factor
  //
  for(int f = 0; f < _fluxional.size(); ++f)
    //
    res *= _fluxional[f].span();

  // reference potential statistical weight
  //
  if(_ref_pot)
    //
    res *= _ref_pot.weight(_ref_tem);
  
  // pi factors for fluxional and external rotation modes
  //
  res *= 2. / std::pow(2. * M_PI, double(_fluxional.size() - 1) / 2.);
  
  // fluxional and external rotation temperature factor
  //
  res *= std::pow(temperature, double(_fluxional.size() + 3) / 2.);

  // non-fluxional modes temperature factor
  //
  itemp = atom_size() * 3 - _fluxional.size() - 6;

  if(_ists)
    //
    --itemp;
  
  res *= std::pow(temperature, double(itemp));

  // symmetry factor
  //
  res /= _symm_fac;

  // electronic energy levels contribution
  //
  dtemp = 0.;
  
  for(std::map<double, int>::const_iterator cit = _elevel.begin(); cit != _elevel.end(); ++cit)
    //
    dtemp += (double)cit->second * std::exp(-cit->first / temperature);

  res *= dtemp;
  
  // non-fluxional mode frequencies factor
  //
  if(_nohess)
    //
    for(int f = 0; f < _nm_freq.size(); ++f)
      //
      res /= _nm_freq[f];
  

  IO::log << IO::log_offset
	  << std::setw(7)  << "T, K"
	  << std::setw(13) << "Z"
	  << std::setw(13) << "Error, %"
	  << "\n"
	  << std::setw(7)  << temperature / Phys_const::kelv
	  << std::setw(13) << res
	  << std::setw(13) << werr
	  << "\n";
  
  return res;
}

// read data from file
//
bool Model::MonteCarlo::_read (std::istream&           from,       // data stream
			       double&                 ener,       // energy
			       Lapack::Vector          cart_pos,   // cartesian coordinates    
			       Lapack::Vector          cart_grad,  // energy gradient in cartesian coordinates
			       Lapack::SymmetricMatrix cart_fc     // cartesian force constant matrix
			       ) const
{
  const char funame [] = "Model::MonteCarlo::_read: ";

  int itemp;

  double dtemp;

  std::string stemp, token, comment;

  // sampling point label
  //
  if(!std::getline(from, comment)) {
    //
    // end of file
    //
    if(from.eof())
      //
      return false;

    std::cerr << funame << "corrupted\n";

    throw Error::Input();
  }

  // energy label
  //
  token = "";
  
  if(!(from >> token) || token != "Energy") {
    //
    std::cerr << funame << "cannot read energy label: " << token << "\n";

    throw Error::Input();
  }

  // energy value in kcal/mol
  //
  if(!(from >> ener)) {
    //
    std::cerr << funame << "cannot read energy\n";

    throw Error::Input();
  }

  ener *= Phys_const::kcal;

  // geometry label
  //
  token = "";
  
  if(!(from >> token) || token != "Geometry") {
    //
    std::cerr << funame << "cannot read geometry label: " << token << "\n";

    throw Error::Input();
  }

  // # of atoms
  //
  itemp = 0;
  
  if(!(from >> itemp) || itemp != atom_size()) {
    //
    std::cerr << funame << "wrong atoms #: " << itemp << "\n";

    throw Error::Input();
  }

  // atom specification
  //
  for(int a = 0; a < atom_size(); ++a) {
    //
    // atom name
    //
    if(!(from >> stemp)) {
      //
      std::cerr << funame << "cannot read " << a + 1 << "-th atom name\n";

      throw Error::Input();
    }

    // atom cartesian coordinates in angstrom
    //
    for(int i = 0; i < 3; ++i) {
      //
      if(!(from >> dtemp)) {
	//
	std::cerr << funame << "cannot read " << a + 1 << "-th atom " << i + 1 << "-th coordinate\n";

	throw Error::Input();
      }

      cart_pos[3 * a + i] = dtemp * Phys_const::angstrom;
    }
  }

  if(_nohess) {
    //
    std::getline(from, comment);

    // uncomment the next line if there is an empty line before next sampling point
    //
    //std::getline(from, comment);
    
    return true;
  }
  // gradient label
  //
  token = "";
  
  if(!(from >> token) || token != "Gradient") {
    //
    std::cerr << funame << "cannot read energy gradient label: " << token << "\n";

    throw Error::Input();
  }

  // gradient components
  //
  for(int i = 0; i < cart_grad.size(); ++i)
    //
    if(!(from >> cart_grad[i])) {
      //
      std::cerr << funame << "cannot read " << i + 1 << "-th energy gradient component\n";

      throw Error::Input();
    }

  // hessian label
  //
  token = "";
  
  if(!(from >> token) || token != "Hessian") {
    //
    std::cerr << funame << "cannot read hessian label: " << token << "\n";

    throw Error::Input();
  }

  // hessian components
  //
  for(int i = 0; i < cart_fc.size(); ++i)
    //
    for(int j = 0; j < cart_fc.size(); ++j) {
      //
      if(!(from >> dtemp)) {
      //
	std::cerr << funame << "cannot read (" << i + 1 << ", " << j + 1 << ")-th hessian component\n";

	throw Error::Input();
      }

      if(i < j)
	//
	continue;

      cart_fc(i, j) = dtemp;
    }
  
  std::getline(from, comment);
  
  // empty line
  //
  std::getline(from, comment);
    
  return true;
}

void Model::MonteCarlo::_set_reference_energy ()
{
  const char funame [] = "Model::MonteCarlo::_set_reference_energy: ";

  int    itemp;
  
  double dtemp;

  if(!_data_file.size()) {
    //
    std::cerr << funame << "no data file\n";

    throw Error::Init();
  }
  
  std::ifstream from(_data_file.c_str());

  if(!from) {
    //
    std::cerr << funame << "cannot open data file " << _data_file << "\n";

    throw Error::File();
  }

  std::string comment;

  std::getline(from, comment);
  
  double res = 0.;

  double ener;

  const int cart_size = atom_size() * 3;
  
  // cartesian coordinates
  //
  Lapack::Vector          cart_pos(cart_size);

  // cartesian gradient
  //
  Lapack::Vector          cart_grad(cart_size);

  // cartesian force constant matrix
  //
  Lapack::SymmetricMatrix cart_fc(cart_size);

  int count = 0;
  
  while(_read(from, ener, cart_pos, cart_grad, cart_fc)) {
    //
    if(!_noqf && !_nohess) {
      //
      for(int c = 0; c < cart_size; ++c)
	//
	for(int d = c; d < cart_size; ++d)
	  //
	  cart_fc(c, d) /= _mass_sqrt[c / 3] * _mass_sqrt[d / 3];
      
    
      Lapack::Vector eval = cart_fc.eigenvalues();

      dtemp = 0.;
    
      for(int c = 0; c < cart_size; ++c)
	//
	if(eval[c] > 0.)
	  //
	  dtemp += std::sqrt(eval[c]);

      // zero-point energy correction
      //
      ener += dtemp / 2.;
    }

    if(!count++ || ener < _refen)
      //
      _refen = ener;
  }

  // non-fluxional modes zero-point energy correction
  //
  if(_nohess && !_noqf)
    //
    for(int f = 0; f < _nm_freq.size(); ++f)
      //
      _refen += _nm_freq[f] / 2.;

  // rounding for better output
  //
  _refen = std::floor(_refen * 10. / Phys_const::kcal) / 10. * Phys_const::kcal;
}

// non-fluxional mode minimal frequency in 1/cm
//
double Model::MonteCarlo::nm_freq_min = 10.;

Lapack::SymmetricMatrix Model::MonteCarlo::_inertia_matrix (Lapack::Vector pos) const
{
  const char funame [] = "Model::MonteCarlo::_inertia_matrix: ";

  double dtemp;

  if(pos.size() != atom_size() * 3) {
    //
    std::cerr << funame << "sizes mismatch: " << pos.size() << " vs. " << atom_size() * 3 << "\n";

    throw Error::Range();
  }
  
  Lapack::SymmetricMatrix res(3);

  dtemp = 0.;
  
  for(int a = 0; a < atom_size(); ++a)
    //
    dtemp += _atom_mass(a) * vdot(pos + 3 * a, 3);

  
  res = dtemp;

  for(int a = 0; a < atom_size(); ++a)
    //
    for(int i = 0; i < 3; ++i)
      //
      for(int j = i; j < 3; ++j)
	//
	res(i, j) -= _atom_mass(a) * pos[3 * a + i] * pos[3 * a + j]; 
  
  return res;
}

void Model::MonteCarlo::_make_cm_shift (Lapack::Vector pos) const
{
  const char funame [] = "Model::MonteCarlo::_make_cm_shift: ";

  double dtemp;

  if(pos.size() != atom_size() * 3) {
    //
    std::cerr << funame << "sizes mismatch: " << pos.size() << " vs. " << atom_size() * 3 << "\n";

    throw Error::Range();
  }
  
  for(int i =  0; i < 3; ++i) {
    //
    double shift = 0.;

    for(int a = 0; a < atom_size(); ++a)
      //
      shift += _atom_mass(a) * pos[a * 3 + i];

    shift /= mass();
    
    for(int a = 0; a < atom_size(); ++a)
      //
      pos[a * 3 + i] -= shift;
  }
}
// local weight contribution for the well
//
double Model::MonteCarlo::_local_weight (double                  ener,       // energy of the sampling 
					 Lapack::Vector          cart_pos,   // cartesian coordinates    
					 Lapack::Vector          cart_grad,  // energy gradient in cartesian coordinates
					 Lapack::SymmetricMatrix cart_fc,    // cartesian force constant matrix
					 double                  temperature // temprature
					 ) const
{
  const char funame [] = "Model::MonteCarlo::_local_weight: ";

  // high frequency threshold x = freq / 2 / T.
  //
  static const double high_freq_thres = 5.;

  // maximal exponent argument
  //
  static const double exp_arg_max = 50.;
  
  /*****************************************************************************
   **************************** CALCULATION BEGINS *****************************
   *****************************************************************************/
  
  int    itemp;
  
  double dtemp;

  IO::Marker funame_marker(funame);
  
  // center of mass shift
  //
  if(_cmshift)
    //
    _make_cm_shift(cart_pos);
  
  // fluxional modes values output
  //
  IO::log << IO::log_offset << "fluxional modes values:";
  
  for(int f = 0; f < _fluxional.size(); ++f)
    //
    IO::log << std::setw(15) << _fluxional[f].evaluate(cart_pos);

  IO::log << "\n";
  
  // fluxional modes first derivatives
  //
  const int cart_size = _mass_sqrt.size() * 3;

  Lapack::Matrix fmfd(cart_size, _fluxional.size());

  for(int c = 0; c < cart_size; ++c) {
    //
    // first derivative signature
    //
    std::vector<int> sign(1);

    sign[0] = c;
    
    for(int f = 0; f < _fluxional.size(); ++f)
      //
      fmfd(c, f) = _fluxional[f].evaluate(cart_pos, sign);
  }

  // force constant curvlinear correction
  //
  if(!_nohess && !_nocurv && !_ists) {
    //
    // potential energy gradient over fluxional modes coordinates
    //
    double residue;
  
    Lapack::Vector flux_grad = Lapack::svd_solve(fmfd, cart_grad, &residue);

    IO::log << IO::log_offset << "energy, kcal/mol, over fluxional modes coordinates gradient:";

    for(int f = 0; f < _fluxional.size(); ++f)
      //
      IO::log << std::setw(15) << flux_grad[f] / Phys_const::kcal;

    IO::log << "\n";
  
    IO::log << IO::log_offset << "cartesian gradient length,  kcal/mol/Bohr: " << std::setw(15)
	    << std::sqrt(cart_grad.vdot()) / Phys_const::kcal << "\n";
  
    IO::log << IO::log_offset << "cartesian gradient residue, kcal/mol/Bohr: " << std::setw(15)
	    << residue / Phys_const::kcal << "\n";

    // check residue: => OK
    //
    /*
      Lapack::Vector appr_grad = fmfd * flux_grad;

      appr_grad -= cart_grad;

      IO::log << IO::log_offset << "test residue,               kcal/mol/Bohr: " << std::setw(15)
      << std::sqrt(appr_grad.vdot()) / Phys_const::kcal << "\n";
    */
    
    for(int c = 0; c < cart_size; ++c) {
      //
      for(int d = c; d < cart_size; ++d) {
	//
	std::vector<int> sign(2);

	sign[0] = c;
	sign[1] = d;

	for(int i = 0; i < _fluxional.size(); ++i)
	  //
	  cart_fc(c, d) -= flux_grad[i] * _fluxional[i].evaluate(cart_pos, sign);
      }
    }
  }

  if(!_nohess)
    //
    for(int c = 0; c < cart_size; ++c)
      //
      for(int d = c; d < cart_size; ++d)
	//
	cart_fc(c, d) /= _mass_sqrt[c / 3] * _mass_sqrt[d / 3];

  // basis orthogonal to external rotations and translations
  //
  Lapack::Matrix basis(cart_size);

  for(int c = 0; c < cart_size; ++c) {
    //
    // atomic index
    //
    const int a = c / 3;

    const int b = c % 3;
    
    // current mode index
    //
    int m = 0;
    
    // overall rotations
    //
    for(int i = 0; i < 3; ++i, ++m) {
      //
      const int i1 = (i + 1) % 3;

      const int i2 = (i + 2) % 3;

      if(b == i1) {
	//
	basis(c, m) =  _mass_sqrt[a] * cart_pos[a * 3 + i2];
      }
      else if(b == i2) {
	//
	basis(c, m) = -_mass_sqrt[a] * cart_pos[a * 3 + i1];
      }
      else
	//
	basis(c, m) = 0.;
    }
  
    // overall translations (normalized)
    //
    for(int i = 0; i < 3; ++i, ++m) {
      //
      if(b == i) {
	//
	basis(c, m) = _mass_sqrt[a] / _total_mass_sqrt;
      }
      else
	//
	basis(c, m) = 0.;
      //      
    }//
    //
  }//

  double wfac = Lapack::orthogonalize(basis, 6);

  // internal configurational space size (without overall rotations and translations)
  //
  const int in_size = cart_size - 6;

  // non-fluxional (tight) modes size
  //
  const int nm_size = in_size - _fluxional.size();

  // Hessian without rotations and translations
  //
  Lapack::SymmetricMatrix fc(in_size);

  if(!_nohess)
    //
    for(int i = 0; i < in_size; ++i)
      //
      for(int j = i; j < in_size; ++j)
	//
	fc(i, j) = cart_fc * &basis(0, i + 6) * &basis(0, j + 6);

  Lapack::Matrix evec(in_size);
    
  Lapack::Vector eval;
  
  bool deep_tunnel = false;
  
  // quantum correction factor
  //
  if(!_noqf && !_nohess) {
    //
    eval = fc.eigenvalues(&evec);

    for(int f = 0; f < in_size; ++f) {
      //
      if(_ists && !f && eval[f] >= 0.) {
	//
	std::cerr << funame << "not a saddle point: lowest frequency = " << std::sqrt(eval[0]) / Phys_const::incm << " 1/cm\n";

	throw Error::Range();
      }
      
      if(eval[f] < 0.) {
	//
	dtemp = std::sqrt(-eval[f]) / temperature / 2.;

	if(dtemp < 0.9 * M_PI) {
	  //
	  wfac *= dtemp / std::sin(dtemp);
	}
	else if(!deep_tunnel) {
	  //
	  deep_tunnel = true;
	
	  std::cerr << funame << "WARNING: the system is in the deep tunneling regime, check the log file\n";

	  IO::log << IO::log_offset << "WARNING: the system is in the deep tunneling regime" << std::endl;
	}
      }
      else if(eval[f] > 0.) {
	//
	dtemp = std::sqrt(eval[f]) / temperature / 2.;
      
	if(dtemp > high_freq_thres) {
	  //
	  ener += dtemp * temperature;
	
	  wfac *= 2. * dtemp;
	}
	else
	  //
	  wfac *= dtemp / std::sinh(dtemp);
      }
    }

    // deep tunneling regime output
    //
    if(deep_tunnel) {
      //
      IO::log << IO::log_offset << "Deep tunneling regime:\n";
      
      IO::log << IO::log_offset << "Energy (including zero-point energy correction) = "
	      << (ener - _refen) / Phys_const::kcal << " kcal/mol" << std::endl;

      IO::log << IO::log_offset << "Frequencies, 1/cm:";

      for(int f = 0; f < in_size; ++f) {
	//
	IO::log << "   ";
    
	if(eval[f] < 0.) {
	  //
	  IO::log << -std::sqrt(-eval[f]) / Phys_const::incm;
	}
	else
	  //
	  IO::log << std::sqrt(eval[f]) / Phys_const::incm;
      }

      IO::log << std::endl;
      //
    }// deep tunneling

  }// quantum correction factor
  
  // transition state calculation
  //
  if(_ists) {
    //
    // fluxional modes in mass-waited coordinates
    //
    Lapack::Matrix mtemp(_fluxional.size(), cart_size);
    
    for(int c = 0; c < cart_size; ++c)
      //
      for(int f = 0; f < _fluxional.size(); ++f)
	//
	mtemp(f, c) = fmfd(c, f) / _mass_sqrt[c / 3];

    // fluxional modes without overal rotations and translations
    //
    mtemp = mtemp * basis;

    Lapack::Matrix fspace(_fluxional.size(), in_size);

    for(int i = 0; i < in_size; ++i)
      //
      for(int f = 0; f < _fluxional.size(); ++f)
	//
	fspace(f, i) = mtemp(f, i + 6);

    // fluxional modes in eigenvector coordinates
    //
    fspace = fspace * evec;

    // potential gradient in mass-waited coordinates
    //
    for(int c = 0; c < cart_size; ++c)
      //
      cart_grad[c] /= _mass_sqrt[c / 3];

    // gradient without rotations and translations
    //
    Lapack::Vector vtemp = cart_grad * basis;

    Lapack::Vector grad(in_size);

    grad = &vtemp[6];
    
    // gradient in eigenvector coordinates
    //
    grad = grad * evec;
    
    // generalized mass factor
    //
    Lapack::SymmetricMatrix smtemp(_fluxional.size());

    smtemp = 0.;
    
    for(int f = 0; f < _fluxional.size(); ++f)
      //
      for(int g = f; g < _fluxional.size(); ++g)
	//
	for(int i = 0; i < in_size; ++i)
	  //
	  smtemp(f, g) += fspace(f, i) * fspace(g, i) / eval[i] / eval[i];

    wfac *= Lapack::Cholesky(smtemp).det_sqrt();

    smtemp = 0.;
    
    for(int f = 0; f < _fluxional.size(); ++f)
      //
      for(int g = f; g < _fluxional.size(); ++g)
	//
	for(int i = 0; i < in_size; ++i)
	  //
	  smtemp(f, g) += fspace(f, i) * fspace(g, i) / eval[i];

    wfac /= fabs(Lapack::SymLU(smtemp).det());
    
    // true fluxional modes complimentary space
    //
    basis.resize(in_size);

    for(int i = 0; i < in_size; ++i)
      //
      for(int f = 0; f < _fluxional.size(); ++f)
	//
	basis(i, f) = fspace(f, i) / eval[i];

    Lapack::orthogonalize(basis, _fluxional.size());

    // gradient in the complimentary space
    //
    vtemp = grad * basis;

    grad.resize(nm_size);

    grad = &vtemp[_fluxional.size()];
    
    // Hessian in complimentary space
    //
    fc.resize(nm_size);

    fc = 0.;
    
    for(int m = 0; m < nm_size; ++m)
      //
      for(int n = m; n < nm_size; ++n)
	//
	for(int i = 0; i < in_size; ++i)
	  //
	  fc(m, n) += eval[i] * basis(i, m + _fluxional.size()) * basis(i, n + _fluxional.size());

    evec.resize(nm_size);
    
    eval = fc.eigenvalues(&evec);

    IO::log << IO::log_offset << "projected frequencies, 1/cm:";
    
    for(int f = 0; f < nm_size; ++f) {
      //
      double freq = eval[f] >= 0. ? std::sqrt(eval[f]) : -std::sqrt(-eval[f]);

      IO::log << "   " << freq / Phys_const::incm;
      
      if(!f) {
	//
	if(freq >= -nm_freq_min) {
	  //
	  std::cerr << funame << "out of range: imaginary frequency = " << freq /Phys_const::incm << "\n";
      
	  throw Error::Logic();
	}

	continue;
      }
      
      if(freq < nm_freq_min) {
	//
	std::cerr << funame << f + 1 << "-th projected frequency is too low: " << freq / Phys_const::incm << " 1/cm\n";

	throw Error::Logic();
      }

      wfac /= freq;
    }

    IO::log << "\n";
    
    // gradient in the eigenvector complimentary space
    //
    grad = grad * evec;

    // energy correction
    //
    dtemp = 0.;
    
    for(int i = 0; i < nm_size; ++i)
      //
      dtemp -= grad[i] * grad[i] / eval[i];

    dtemp /= 2.;

    IO::log << IO::log_offset << "Energy correction = " << dtemp / Phys_const::kcal << " kcal/mol\n";

    ener += dtemp;
    //
  }// transition state calculation
  //
  // stable species calculation
  //
  else {
    //
    for(int c = 0; c < cart_size; ++c) {
      //
      // fluxional modes start index 
      //
      itemp = 6;
    
      // fluxional modes in mass-weighted coordinates
      //
      for(int i = 0; i < _fluxional.size(); ++i, ++itemp)
	//
	basis(c, itemp) = fmfd(c, i) / _mass_sqrt[c / 3];
    }
  
    // square root of fluxional modes masses 
    //
    wfac /= Lapack::orthogonalize(basis, 6 + _fluxional.size());

    // non-fluxional modes frequency factor
    //
    if(!_nohess) {
      //
      Lapack::SymmetricMatrix nm_fc(nm_size);

      for(int i = 0; i < nm_size; ++i)
	//
	for(int j = i; j < nm_size; ++j)
	  //
	  nm_fc(i, j) = cart_fc * &basis(0, i + 6 + _fluxional.size()) * &basis(0, j + 6 + _fluxional.size());

      Lapack::Vector eval = nm_fc.eigenvalues();

      // after testing uncomment next clause
      //
      // if(deep_tunnel) {
      //
      IO::log << IO::log_offset << "non-fluxional modes (" << nm_size << ") frequencies, 1/cm:";

      for(int f = 0; f < nm_size; ++f) {
	//
	IO::log << "   ";
    
	if(eval[f] < 0.) {
	  //
	  IO::log << -std::sqrt(-eval[f]) / Phys_const::incm;
	}
	else
	  //
	  IO::log << std::sqrt(eval[f]) / Phys_const::incm;
      }

      IO::log << std::endl;
      //
      // }
  
      // check for soft non-fluxional modes
      //
      if(eval[0] < 0.) {
	//
	dtemp = -std::sqrt(-eval[0]);
      }
      else
	//
	dtemp = std::sqrt(eval[0]);

      dtemp /= Phys_const::incm;
    
      if(dtemp < nm_freq_min) {
	//
	std::cerr << funame << "non-fluxional mode frequency is too low: " << dtemp  << " 1/cm \n";

	return -1.;
      }
  
      // classical weight factor for non-fluxional modes
      //
      for(int f = 0; f < nm_size; ++f) {
	//
	wfac /= std::sqrt(eval[f]);
      }
    }

    if(_nohess && !_noqf) {
      //
      for(int f = 0; f < _nm_freq.size(); ++f) {
	//
	dtemp = _nm_freq[f] / temperature / 2.;
      
	if(dtemp > high_freq_thres) {
	  //
	  ener += dtemp * temperature;
	
	  wfac *= 2. * dtemp;
	}
	else
	  //
	  wfac *= dtemp / std::sinh(dtemp);
      }
    }
  }// stable species
  
  // Boltzmann factor, including zero-point energy of high-frequency modes relative to the reference energy
  //
  double bpow = (ener - _refen) / temperature;

  // reference potential correction
  //
  if(_ref_pot) {
    //
    Lapack::Vector flux_pos(_fluxional.size());

    for(int f = 0; f < _fluxional.size(); ++f)
      //
      flux_pos[f] = _fluxional[f].evaluate(cart_pos);
  
    bpow -= _ref_pot(flux_pos) / _ref_tem;
  }

  if(bpow > exp_arg_max)
    //
    return 0.;

  if(bpow < -exp_arg_max) {
    //
    std::cerr << funame << "energy (relative to the reference energy) is too low: "
	      << (ener - _refen) / Phys_const::kcal
	      << " kcal/mol. Check the reference energy\n";

    throw Error::Range();
  }
  
  return wfac * std::exp(-bpow);
}

/***********************************************************************************************************
 ****************************** CRUDE MONTE-CARLO SAMPLING WITH DUMMY ATOMS ********************************
 ***********************************************************************************************************/

Model::MonteCarloWithDummy::~MonteCarloWithDummy() {}

Model::MonteCarloWithDummy::MonteCarloWithDummy(IO::KeyBufferStream& from, const std::string& n, int m) : Species(from, n, m)
{
  const char funame [] = "Model::MonteCarloWithDummy::MonteCarloWithDummy: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  
  double dtemp;

  std::string token, comment, stemp;

  KeyGroup MonteCarloWithDummyModel;

  Key  atom_key("MoleculeSpecification" );
  Key  flux_key("FluxionalMode"         );
  Key   con_key("Constrain"             );
  Key  incr_key("NumericIncrement[Bohr]");
  Key  data_key("DataFile"              );
  
  while(from >> token) {
    //
    // end input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // molecule specification
    //
    else if(atom_key == token) {
      //
      if(_atom_array.size()) {
	//
	std::cerr << funame << token << ": already defined\n";

	throw Error::Init();
      }

      IO::LineInput lin(from);

      int na;
      
      if(!(lin >> na)) {
	//
	std::cerr << funame << token << ": cannot read number of atoms\n";

	throw Error::Input();
      }

      if(na <= 0) {
	//
	std::cerr << funame << token << ": number of atoms out of range: " << na << "\n";

	throw Error::Range();
      }

      for(int a = 0; a < na; ++a) {
	//
	_atom_array.push_back(Atom(from));

	if(_atom_array.back().number()) {
	  //
	  _real_atom.push_back(a);
	}
	else
	  //
	  _dummy_atom.push_back(a);
      }
    }
    // fluxional mode specification
    //
    else if(flux_key == token) {
      //
      std::getline(from, comment);

      if(!_atom_array.size()) {
	//
	std::cerr << funame << token << ": molecule specification should go first\n";

	throw Error::Init();
      }

      _fluxional.push_back(Fluxional(_atom_array.size(), from));
    }
    // dummy constrain specification
    //
    else if(con_key == token) {
      //
      std::getline(from, comment);

      if(!_atom_array.size()) {
	//
	std::cerr << funame << token << ": molecule specification should go first\n";

	throw Error::Init();
      }

      _constrain.push_back(Constrain(_atom_array.size(), from));
    }
    // numerical differentiation increment
    //
    else if(incr_key == token) {
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

      IntMod::increment = dtemp;
    }
    // data file
    //
    else if(data_key == token) {
      //
      IO::LineInput lin(from);

      if(_data_file.size()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }
      
      if(!(lin >> _data_file)) {
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

  // some checking
  //
  if(!from) {
    //
    std::cerr << funame << "input stream is currupted\n";

    throw Error::Input();
  }

  if(!_atom_array.size()) {
    //
    std::cerr << funame << "no molecule specification\n";

    throw Error::Init();
  }

  if(!_fluxional.size()) {
    //
    std::cerr << funame << "no fluxional modes specified\n";

    throw Error::Init();
  }

  if(_fluxional.size() + 6 > _real_atom.size() * 3) {
    //
    std::cerr << funame << "number of fluxional modes, " << _fluxional.size()
	      << ", is inconsistent with the number of real atoms, " << _real_atom.size() << "\n";

    throw Error::Init();
  }

  if(_constrain.size() != _dummy_atom.size() * 3) {
    //
    std::cerr << funame << "number of dummy constrains, " << _constrain.size() << ", is inconsistent with the number of dummy atoms, "
	      << _dummy_atom.size() << "\n";

    throw Error::Init();
  }

  if(!_data_file.size()) {
    //
    std::cerr << funame << "no data file is provided\n";

    throw Error::Init();
  }
}

void Model::MonteCarloWithDummy::_assert(ConstSlice<double> v) const
{
  const char funame [] = "Model::MonteCarloWithDummy::_assert: ";

  if(v.size() == _atom_array.size() * 3)
    //
    return;

  std::cerr << funame << "vector size out of range: " << v.size();
  
  throw Error::Range();
}

double Model::MonteCarloWithDummy::_vdot(ConstSlice<double> v1, ConstSlice<double> v2) const
{
  const char funame [] = "Model::MonteCarloWithDummy::_vdot: ";

  double dtemp;
  
  int    itemp;

  _assert(v1);

  _assert(v2);

  double res = 0.;

  for(int r = 0; r < _real_atom.size(); ++r) {
    //
    const int a = _real_atom[r];

    dtemp = 0.;

    itemp = a * 3;
    
    for(int i = 0; i < 3; ++i)
      //
      dtemp += v1[itemp + i] * v2[itemp + i];

    res += dtemp * _atom_array[a].mass();
  }

  return res;
}

double Model::MonteCarloWithDummy::_normalize(Slice<double> a) const
{
  const char funame [] = "Model::MonteCarloWithDummy::_normalize: ";

  static const double eps = 1.e-20;
  
  double aa = _vdot(a, a);

  if(aa < eps) {
    //
    std::cerr << funame << "vector has zero length\n";

    throw Error::Range();
  }

  aa = std::sqrt(aa);
  
  for(Slice<double>::iterator it = a.begin(); it != a.end(); ++it)
    //
    *it /= aa;

  return aa;
}
    
void Model::MonteCarloWithDummy::_orthogonalize (Slice<double> a, ConstSlice<double> b) const
{
  double ab = _vdot(a, b);
  
  double bb = -vdot(b, b);

  for(int i = 0; i < a.size(); ++i)
    //
    a[i] -= b[i] * ab / bb;
}
    
void Model::MonteCarloWithDummy::_orthogonalize (Lapack::Matrix m) const
{
  for(int c = 0; c < m.size2(); ++c) {
    //
    for(int d  = 0; d < c; ++d)
      //
      _orthogonalize(m.column(c), m.column(d));

    _normalize(m.column(c));
  }
}
    
double Model::MonteCarloWithDummy::_local_weight (double                  ener,           // energy of the sampling 
						  Lapack::Vector          cart_pos,       // cartesian coordinates of all atoms (real and dummy ones)   
						  Lapack::Vector          real_cart_grad, // energy gradient in cartesian coordinates
						  Lapack::SymmetricMatrix real_cart_fc,   // cartesian force constant matrix
						  double                  temperature     // temprature
						  ) const
{
  const char funame [] = "Model::MonteCarloWithDummy::_local_weight: ";

  // non-fluxional mode minimal frequency in 1/cm
  //
  static const double nm_freq_min = 10.;

  // high frequency threshold x = freq / 2 / T.
  //
  static const double high_freq_thres = 5.;

  // maximal exponent argument
  //
  static const double exp_arg_max = 50.;
  
  /*****************************************************************************
   **************************** CALCULATION BEGINS *****************************
   *****************************************************************************/
  
  int    itemp;
  
  double dtemp;

  IO::Marker funame_marker(funame);

  // full cartesian size including dummy atoms
  //
  const int cart_size = _atom_array.size() * 3;

  // number of curvlinear modes, fluxional ones and constrains
  //
  const int curv_size = _fluxional.size() + _constrain.size();
  
  // internal (fluxional and constrain) modes first derivatives
  //
  Lapack::Matrix imfd(cart_size, curv_size);

  for(int c = 0; c < cart_size; ++c) {
    //
    // first derivative signature
    //
    std::vector<int> sign(1, c);

    itemp = 0;
    
    for(int i = 0; i < _constrain.size(); ++i, ++itemp)
      //
      imfd(c, itemp) = _constrain[i].evaluate(cart_pos, sign);

    for(int i = 0; i < _fluxional.size(); ++i, ++itemp)
      //
      imfd(c, itemp) = _fluxional[i].evaluate(cart_pos, sign);
  }

  // potential energy gradient over all cartesian coordinates (real and dummy ones)
  //
  Lapack::Vector cart_grad(cart_size);

  cart_grad = 0.;

  for(int r = 0; r < _real_atom.size(); ++r)
    //
    for(int ri = 0; ri < 3; ++ri)
      //
      cart_grad[_real_atom[r] * 3 + ri] = real_cart_grad[r * 3 + ri];

  // energy gradient over internal coordinates (fluxional and constrained ones)
  //
  double residue;
  
  Lapack::Vector curv_grad = Lapack::svd_solve(imfd, cart_grad, &residue);
  
  // force constant matrix over all cartesian coordinates (real and dummy ones)
  //
  Lapack::SymmetricMatrix cart_fc(cart_size);

  cart_fc = 0.;

  for(int r = 0; r < _real_atom.size(); ++r) {
    //
    for(int ri = 0; ri < 3; ++ri) {
      //
      for(int s = r; s < _real_atom.size(); ++s) {
	//
	for(int si = 0; si < 3; ++si) {
	  //
	  if(s == r && si < ri)
	    //
	    continue;

	  cart_fc(_real_atom[r] * 3 + ri, _real_atom[s] * 3 + si) = real_cart_fc(r * 3 + ri, s * 3 + si);
	}
      }
    }
  }

  // modified cartesian force constant matrix
  //
  for(int c = 0; c < cart_size; ++c)
    //
    for(int d = c; d < cart_size; ++d) {
      //
      std::vector<int> sign(2);

      sign[0] = c;
      
      sign[1] = d;

      itemp = 0;
      
      for(int i = 0; i < _constrain.size(); ++i, ++itemp)
	//
	cart_fc(c, d) -= curv_grad[itemp] * _constrain[i].evaluate(cart_pos, sign);

      for(int i = 0; i < _fluxional.size(); ++i, ++itemp)
	//
	cart_fc(c, d) -= curv_grad[itemp] * _fluxional[i].evaluate(cart_pos, sign);
    }
  
  /****************** FLUXIONAL MODES MASS FACTOR ************************/
  //

  // derivatives of dummy constrain over dummy cartesian coordinates
  //
  Lapack::Matrix dqdx(_constrain.size());

  for(int d = 0; d < _dummy_atom.size(); ++d)
    //
    for(int di = 0; di < 3; ++di)
      //
      for(int c = 0; c < _constrain.size(); ++c)
	//
	dqdx(c, d * 3 + di) = imfd(_dummy_atom[d] * 3 + di, c);

  // mass-weighted fluxional modes constrained derivatives
  //
  Lapack::Matrix fmcd(_real_atom.size() * 3, _fluxional.size());

  for(int r = 0; r < _real_atom.size(); ++r) {
    //
    double mass_sqrt = std::sqrt(_atom_array[_real_atom[r]].mass());
    //
    for(int ri = 0; ri < 3; ++ri) {

      // dummy constrain derivative over real atom cartesian coordinate
      //
      Lapack::Vector dqdr(_constrain.size());

      for(int c = 0; c < _constrain.size(); ++c)
	//
	dqdr[c] = imfd(_real_atom[r] * 3 + ri, c);

      // dummy atom cartesian coordinate constrained derivative over real atom cartesian coordinate
      //
      Lapack::Vector dxdr = Lapack::svd_solve(dqdx, dqdr);

      for(int f = 0; f < _fluxional.size(); ++f) {
	//
	// fluxional mode first derivative over real atom cartesian coordinate
	//
	dtemp = imfd(_real_atom[r] * 3 + ri, _constrain.size() + f);

	for(int d = 0; d < _dummy_atom.size(); ++d)
	  //
	  for(int di = 0; di < 3; ++di)
	    //
	    dtemp -= imfd(_dummy_atom[d] * 3 + di, _constrain.size() + f) * dxdr[d * 3 + di];

	fmcd(r * 3 + ri, f) = dtemp / mass_sqrt;
      }
    }
  }
  
  // orthogonalize mass-weighted fluxional modes constrained cartesian derivatives
  //
  double wfac = 1. / fmcd.orthogonalize();

  // dummy constrain matrix for full frequency calculations
  //
  Lapack::Matrix dcm(_constrain.size() + 6, cart_size);

  dcm = 0.;

  // translation
  //
  for(int i = 0; i < 3; ++i) {
    //
    for(int r = 0; r < _real_atom.size(); ++r) {
      //
      const int a = _real_atom[r];
      
      dcm(i, a * 3 + i) = _atom_array[a].mass();
    }
  }
  
  // overall rotations
  //
  for(int i = 0; i < 3; ++i) {
    //
    const int i1 = (i + 1) % 3;

    const int i2 = (i + 2) % 3;

    for(int r = 0; r < _real_atom.size(); ++r) {
      //
      const int a = _real_atom[r];
      
      dcm(i + 3, a * 3 + i1) =  _atom_array[a].mass() * cart_pos[a * 3 + i2];

      dcm(i + 3, a * 3 + i2) = -_atom_array[a].mass() * cart_pos[a * 3 + i1];
    }
  }

  // constrain modes
  //
  for(int i = 0; i < _constrain.size(); ++i)
    //
    for(int c = 0; c < cart_size; ++c)
      //
      dcm(i + 6, c) = imfd(c, i);

  // constrained subspace basis
  //
  Lapack::Matrix basis = dcm.kernel();

  _orthogonalize(basis);

  // internal (fluxional and non-fluxional) modes frequencies
  //
  const int in_size = basis.size2(); // should be 3 * "real atom #" - 6
  
  Lapack::SymmetricMatrix in_fc(in_size);

  for(int i = 0; i < in_size; ++i)
    //
    for(int j = i; j < in_size; ++j)
      //
      in_fc(i, j) = cart_fc * basis.column(i) * basis.column(j);
  
  // eigenvalues
  //
  Lapack::Vector eval = in_fc.eigenvalues();

  // quantum correction factor
  //
  bool deep_tunnel = false;
  
  for(int f = 0; f < in_size; ++f) {
    //
    if(eval[f] < 0.) {
      //
      dtemp = std::sqrt(-eval[f]) / temperature / 2.;

      if(dtemp > 0.9 * M_PI && !deep_tunnel) {
	//
	deep_tunnel = true;
	
	std::cerr << funame << "WARNING: the system is in the deep tunneling regime, check the log file\n";

	IO::log << IO::log_offset << "WARNING: the system is in the deep tunneling regime" << std::endl;
      }
      else {
	//
	wfac *= dtemp / std::sin(dtemp);
      }
    }
    else if(eval[f] > 0.) {
      //
      dtemp = std::sqrt(eval[f]) / temperature / 2.;
      
      if(dtemp > high_freq_thres) {
	//
	ener += dtemp * temperature;
	
	wfac *= 2. * dtemp;
      }
      else
	//
	wfac *= dtemp / std::sinh(dtemp);
    }
  }

  // internal mode frequencies output
  
  // after testing uncomment next clause
  //
  // if(deep_tunnel) {
  //
  IO::log << IO::log_offset << "internal modes (" << in_size << ") frequencies, 1/cm:";

  for(int f = 0; f < in_size; ++f) {
    //
    IO::log << "   ";
    
    if(eval[f] < 0.) {
      //
      IO::log << -std::sqrt(-eval[f]) / Phys_const::incm;
    }
    else
      //
      IO::log << std::sqrt(eval[f]) / Phys_const::incm;
  }

  IO::log << std::endl;

  // }
  
  // non-fluxional modes frequencies
  //
  const int nm_size = in_size - _fluxional.size();

  // full constrain matrix, with both dummy and fluxional constrains
  //
  Lapack::Matrix fcm(_constrain.size() + _fluxional.size() + 6, cart_size);

  for(int c = 0; c < cart_size; ++c) {
    //
    for(int i = 0; i < dcm.size1(); ++i)
      //
      fcm(i, c) = dcm(i, c);

    for(int i = 0; i < _fluxional.size(); ++i)
      //
      fcm(i + dcm.size1(), c) = imfd(c, i + _constrain.size());
  }
  
  basis = fcm.kernel();

  _orthogonalize(basis);
  
  // non-fluxional modes force constant matrix
  //
  Lapack::SymmetricMatrix nm_fc(nm_size);

  for(int i = 0; i < nm_size; ++i)
    //
    for(int j = i; j < nm_size; ++j)
      //
      nm_fc(i, j) = cart_fc * basis.column(i) * basis.column(j);

  // eigenvalues
  //
  eval = nm_fc.eigenvalues();

  // after testing uncomment next clause
  //
  // if(deep_tunnel) {
  //
  IO::log << IO::log_offset << "non-fluxional modes (" << nm_size << ") frequencies, 1/cm:";

  for(int f = 0; f < nm_size; ++f) {
    //
    IO::log << "   ";
    
    if(eval[f] < 0.) {
      //
      IO::log << -std::sqrt(-eval[f]) / Phys_const::incm;
    }
    else
      //
      IO::log << std::sqrt(eval[f]) / Phys_const::incm;
  }

  IO::log << std::endl;
  //
  // }
  
  // check for soft non-fluxional modes
  //
  if(eval[0] < 0.) {
    //
    dtemp = -std::sqrt(-eval[0]);
  }
  else
    //
    dtemp = std::sqrt(eval[0]);

  dtemp /= Phys_const::incm;
    
  if(dtemp < nm_freq_min) {
    //
    std::cerr << funame << "non-fluxional mode frequency is too low: " << dtemp  << " 1/cm \n";

    throw Error::Range();
  }
  
  // classical weight factor for non-fluxional modes
  //
  for(int f = 0; f < nm_size; ++f) {
    //
    wfac /= std::sqrt(eval[f]);
  }

  // deep tunneling regime energy output
  if(deep_tunnel) {
    //
    IO::log << IO::log_offset << "energy (including zero-point energy correction), kcal/mol = "
	    << (ener - ground()) / Phys_const::kcal << std::endl;
  }
  
  // Boltzmann factor, including zero-point energy of high-frequency modes
  //
  dtemp = (ener - ground()) / temperature;

  if(dtemp > exp_arg_max)
    //
    return 0.;

  if(dtemp < -exp_arg_max) {
    //
    std::cerr << funame << "energy is too low: " << (ener - ground()) / Phys_const::kcal << " kcal/mol. Check the ground energy\n";

    throw Error::Range();
  }
  
  return wfac * std::exp(-dtemp);
}

double Model::MonteCarloWithDummy::states (double ener) const
{
  const char funame [] = "Model::MonteCarloWithDummy::states: ";

  std::cerr << funame << "not implemented yet, sorry\n";

  throw Error::General();

  return 0.;
}

double Model::MonteCarloWithDummy::weight (double temperature) const
{
  const char funame [] = "Model::MonteCarloWithDummy::weight: ";

  static bool first_time = true;

  int    itemp;
  
  double dtemp;
  
  std::ifstream from(_data_file.c_str());

  if(!from) {
    //
    std::cerr << funame << "cannot open data file " << _data_file << "\n";

    throw Error::File();
  }

  double res = 0.;

  Lapack::Vector          cart_pos(_atom_array.size() * 3);
    
  Lapack::Vector          cart_grad(_real_atom.size() * 3);
    
  Lapack::SymmetricMatrix cart_fc  (_real_atom.size() * 3);

  double ener;

  int count = 0;

  std::vector<std::pair<double, double> > flimits(_fluxional.size());

  std::vector<double> con_rms(_constrain.size());
  
  while(from >> ener) {
    //
    // read energies, atom cartesian positions, gradients, and force constants
    // ...

    if(!from) {
      //
      std::cerr << funame << "data file " << _data_file << " is corrupted\n";

      throw Error::Input();
    }
    
    res += _local_weight(ener, cart_pos, cart_grad, cart_fc, temperature);

    if(first_time) {
      //
      for(int f = 0; f < _fluxional.size(); ++f) {
	//
	dtemp = _fluxional[f].evaluate(cart_pos);

	if(!count || dtemp < flimits[f].first)
	  //
	  flimits[f].first = dtemp;

	if(!count || dtemp > flimits[f].second)
	  //
	  flimits[f].second = dtemp;
      }

      for(int c = 0; c < _constrain.size(); ++c) {
	//
	dtemp = _constrain[c].evaluate(cart_pos) - _constrain[c].value();

	con_rms[c] += dtemp * dtemp;
      }
    }
    
    ++count;
  }

  if(!count) {
    //
    std::cerr << funame << "no data\n";

    throw Error::Input();
  }

  if(first_time) {
    //
    IO::log << IO::log_offset << funame << "\n";
    
    IO::log << IO::log_offset << "Constrain rms deviation:\n";
    
    for(int c = 0; c < _constrain.size(); ++c)
      //
      IO::log << IO::log_offset
	      << std::setw(3) << c
	      << std::setw(15) << std::sqrt(con_rms[c] / (double)count)
	      << "\n";

    IO::log << IO::log_offset << "Fluxional modes real span versus assumed one:\n";
    
    for(int f = 0; f < _fluxional.size(); ++f)
      //
      IO::log << IO::log_offset
	      << std::setw(3) << f
	      << std::setw(15) << flimits[f].second - flimits[f].first
	      << std::setw(15) << _fluxional[f].span()
	      << "\n";

    first_time = false;
  }
  
  // normalize
  //
  res /= double(count);

  // span factor
  //
  for(int f = 0; f < _fluxional.size(); ++f)
    //
    res *= _fluxional[f].span();

  // pi factors for fluxional and external rotation modes
  //
  res *= 2. / std::pow(2. * M_PI, double(_fluxional.size() - 1) / 2.);
  
  // fluxional and external rotation temperature factor
  //
  res *= std::pow(temperature, double(_fluxional.size() + 3) / 2.);

  // non-fluxional modes temperature factor
  //
  res *= std::pow(temperature, double(_real_atom.size() * 3 - _fluxional.size() - 6));
  
  return res;
}

/************************* RIGID ROTOR HARMONIC OSCILLATOR MODEL **************************/

Model::RRHO::RRHO(IO::KeyBufferStream& from, const std::string& n, int m) 
  :Species(from, n, m),  _emax(-1.), _sym_num(1.)
{
  const char funame [] = "Model::RRHO::RRHO: ";

  IO::Marker funame_marker(funame);

  KeyGroup RRHO_Model;

  Key incm_ener_key("ElectronicEnergy[1/cm]"          );
  Key kcal_ener_key("ElectronicEnergy[kcal/mol]"      );
  Key   kj_ener_key("ElectronicEnergy[kJ/mol]"        );
  Key   ev_ener_key("ElectronicEnergy[eV]"            );
  Key   au_ener_key("ElectronicEnergy[au]"            );
  Key incm_zero_key("ZeroEnergy[1/cm]"                );
  Key kcal_zero_key("ZeroEnergy[kcal/mol]"            );
  Key   kj_zero_key("ZeroEnergy[kJ/mol]"              );
  Key   ev_zero_key("ZeroEnergy[eV]"                  );
  Key   au_zero_key("ZeroEnergy[au]"                  );
  Key incm_elev_key("ElectronicLevels[1/cm]"          );
  Key kcal_elev_key("ElectronicLevels[kcal/mol]"      );
  Key   kj_elev_key("ElectronicLevels[kJ/mol]"        );
  Key   ev_elev_key("ElectronicLevels[eV]"            );
  Key   au_elev_key("ElectronicLevels[au]"            );
  Key      emax_key("InterpolationEnergyMax[kcal/mol]");
  Key     estep_key("InterpolationEnergyStep[1/cm]"   );
  Key     extra_key("ExtrapolationStep"               );
  Key      freq_key("Frequencies[1/cm]"               );
  Key    fscale_key("FrequencyScalingFactor"          );
  Key     degen_key("FrequencyDegeneracies"           );
  Key      harm_key("AreFrequenciesHarmonic"          );
  Key    anharm_key("Anharmonicities[1/cm]"           );
  Key      hrot_key("Rotor"                           );
  Key      core_key("Core"                            );
  Key      tunn_key("Tunneling"                       );
  Key      irin_key("InfraredIntensities[km/mol]"     );
  Key     graph_key("GraphPerturbationTheory"         );
  Key       sym_key("SymmetryFactor"                  );
 
  double extra_step = 0.1;
  double ener_quant = Phys_const::incm;
  double fscale     = -1.;

  int         itemp;
  double      dtemp;
  std::string stemp;

  bool isener = false;
  bool iszero = false;
  bool isharm = false;

  std::map<double, int> elevel_map;

  std::string token, comment;
  //
  while(from >> token) {
    // end input
    //
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // frequency scaling factor
    //
    else if(fscale_key == token) {
      //
      if(fscale > 0.) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": already initialized";
      }
	  
      if(!(from >> fscale)) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": corrupted";
      }
      
      if(fscale <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // are frequencies harmonic
    //
    else if(harm_key == token) {
      //
      std::getline(from, comment);
      
      isharm = true;
    }
    // frequency degeneracies
    //
    else if(degen_key == token) {
      //
      if(_fdegen.size()) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": already initialized";
      }
      
      // number of degeneracies
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": cannot read number of frequency degeneracies";
      }
      std::getline(from, comment);

      if(itemp < 1) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": number of frequency degeneracies should be positive";
      }
      _fdegen.resize(itemp);
      
      // read frequency degeneracies
      //
      for(int i = 0; i < _fdegen.size(); ++i) {
	//
	if(!(from >> itemp)) {
	  //
	  ErrOut err_out;
	  
	  err_out << funame << token << ": cannot read " << i << "-th degeneracy";
	}
	
	if(itemp <= 0) {
	  //
	  ErrOut err_out;
	  
	  err_out << funame << token << ": " << i << "-th degeneracy: should be positive";
	}
	_fdegen[i] = itemp;
      }

      if(!from) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": corrupted";
      }
      
      std::getline(from, comment);
    }
    // anharmonicities
    //
    else if(anharm_key == token) {
      //
      if(!_frequency.size()) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": frequencies should be initialized first";
      }
	  
      std::getline(from, comment);
      
      _anharm.resize(_frequency.size());
      //
      _anharm = 0.;
      
      for(int i = 0; i < _anharm.size(); ++i)
	//
	for(int j = 0; j <= i; ++j)
	  //
	  if(!(from >> _anharm(i, j))) {
	    //
	    ErrOut err_out;
	    
	    err_out << funame << token << ": cannot read anharmonicities";
	  }
      
      _anharm *= Phys_const::incm;
      
    }
    // symmetry number
    //
    else if(sym_key == token) {
      if(!(from >> _sym_num)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
   
      if(_sym_num <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // graph perturbation theory
    //
    // ONLY WORK WITH NO DEGENERACIES
    //
    else if(graph_key == token) {      
      if(!_frequency.size()) {
	std::cerr << funame << token << ": frequencies should be initialized first\n";
	throw Error::Init();
      }
      
      std::getline(from, comment);

      _init_graphex(from);
    }
    // tunneling
    else if(tunn_key == token) {
      if(_tunnel) {
	std::cerr << funame << token << ": already intialized\n";
	throw Error::Init();
      }
      if(mode() == DENSITY) {
	std::cerr << funame << token << ": only for barriers\n";
	throw Error::Init();
      }
      _tunnel = new_tunnel(from);
    }
    // RRHO core
    else if(core_key == token) {
      if(_core) {
	std::cerr << funame << token << ": already intialized\n";
	throw Error::Init();
      }
      _core = new_core(from, geometry(), mode());
    }
    // hindered rotor
    else if(hrot_key == token) {
      IO::log << IO::log_offset << _rotor.size() + 1 << "-th ROTOR:\n";
      _rotor.push_back(new_rotor(from, geometry()));
    }
    // frequency scaling factor
    else if(fscale_key == token) {
      if(fscale > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      if(!(from >> fscale)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(fscale <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // frequencies
    else if(freq_key == token) {
      if(_frequency.size()) {
	std::cerr << funame << token << ": have been initialized already\n";
	throw Error::Init();
      }
      // number of frequencies
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": cannot read number of frequencies\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(itemp < 1) {
	std::cerr << funame << token << ": number of frequencies should be positive\n";
	throw Error::Range();
      }
      _frequency.resize(itemp);
      // read frequencies
      for(int i = 0; i < _frequency.size(); ++i) {
	if(!(from >> dtemp)) {
	  std::cerr << funame << token << ": cannot read " << i << "-th frequency\n";
	  throw Error::Input();
	}
	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": " << i << "-th frequency: should be positive\n";
	  throw Error::Range();
	}
	_frequency[i] = dtemp * Phys_const::incm;
      }
      if(!from) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // infrared intensities
    else if(irin_key == token) {
      if(mode() != DENSITY) {
	std::cerr << funame << token << "only for wells\n";
	throw Error::Logic();
      }
      if(_osc_int.size()) {
	std::cerr << funame << token << ": have been initialized already\n";
	throw Error::Init();
      }
      // number of frequencies
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": cannot read number of oscillators\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(itemp < 1) {
	std::cerr << funame << token << ": number of oscillators should be positive\n";
	throw Error::Range();
      }
      _osc_int.resize(itemp);
      // read intensities
      for(int i = 0; i < _osc_int.size(); ++i) {
	if(!(from >> dtemp)) {
	  std::cerr << funame << token << ": cannot read " << i << "-th intensity\n";
	  throw Error::Input();
	}
	if(dtemp < 0.) {
	  std::cerr << funame << token << ": " << i << "-th intensity: should be positive\n";
	  throw Error::Range();
	}
	_osc_int[i] = dtemp * 1.e5 * Phys_const::cm / Phys_const::avogadro;
      }
      if(!from) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    //  frequency quantum
    else if(estep_key == token) {
      if(!(from >> ener_quant)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(ener_quant <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      ener_quant *= Phys_const::incm;
    }
    // interpolation maximal energy
    else if(emax_key == token) {
      if(!(from >> _emax)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_emax <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      _emax *= Phys_const::kcal;
    }
    // extrapolation step
    else if(extra_key == token) {
      if(!(from >> extra_step)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(extra_step <= 0. || extra_step >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    //  Electronic energy
    else if(incm_ener_key == token || kj_ener_key == token || 
	    kcal_ener_key == token || au_ener_key == token || ev_ener_key == token) {
      if(isener || iszero) {
	std::cerr << funame << "ground energy has been initialized already\n";
	throw Error::Init();
      }
      isener = true;
      if(!(from >> _ground)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(incm_ener_key == token)
	_ground *= Phys_const::incm;
      if(kcal_ener_key == token)
	_ground *= Phys_const::kcal;
      if(kj_ener_key == token)
	_ground *= Phys_const::kjoul;
      if(ev_ener_key == token)
	_ground *= Phys_const::ev;
    }
    //  Electronic energy plus zero-point energy
    else if(incm_zero_key == token || ev_zero_key == token ||
	    kcal_zero_key == token || kj_zero_key == token || au_zero_key == token) {
      if(isener || iszero) {
	std::cerr << funame << "ground energy has been initialized already\n";
	throw Error::Init();
      }
      iszero = true;
      if(!(from >> _ground)) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(incm_zero_key == token)
	_ground *= Phys_const::incm;
      if(kcal_zero_key == token)
	_ground *= Phys_const::kcal;
      if(kj_zero_key == token)
	_ground *= Phys_const::kjoul;
      if(ev_zero_key == token)
	_ground *= Phys_const::ev;
    }
    // electronic energy levels & degeneracies
    else if(incm_elev_key == token || au_elev_key == token ||
	    kcal_elev_key == token || kj_elev_key == token || ev_elev_key == token) {
      int num;
      if(!(from >> num)) {
	std::cerr << funame << token << ": levels number unreadable\n";
	throw Error::Input();
      }

      if(num < 1) {
	std::cerr << funame << token << ": levels number should be positive\n";
	throw Error::Range();
      }

      std::getline(from, comment);

      for(int l = 0; l < num; ++l) {
	IO::LineInput level_input(from);
	if(!(level_input >> dtemp >> itemp)) {
	  std::cerr << funame << token << ": format: energy(1/cm) degeneracy(>=1)\n";
	  throw Error::Input();
	}

	if(incm_elev_key == token)
	  dtemp *= Phys_const::incm;
	if(kcal_elev_key == token)
	  dtemp *= Phys_const::kcal;
	if(kj_elev_key == token)
	  dtemp *= Phys_const::kjoul;
	if(ev_elev_key == token)
	  dtemp *= Phys_const::ev;

	if(elevel_map.find(dtemp) != elevel_map.end()) {
	  std::cerr << funame << token << ": identical energy levels\n";
	  throw Error::Range();
	}

	if(itemp < 1) {
	  std::cerr << funame << token << ": degeneracy should be positive\n";
	  throw Error::Range();
	}

	elevel_map[dtemp] = itemp;
      }
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }// while(from >> token);
 
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_frequency.size()) {
    //
    if(fscale > 0.)
      //
      for(int f = 0; f < _frequency.size(); ++f)
	//
	_frequency[f] *= fscale;
    
    if(_fdegen.size() > _frequency.size()) {
      std::cerr << funame << "number of frequency degeneracies exceeds number of frequencies\n";
      throw Error::Init();
    }

    for(int i = _fdegen.size(); i < _frequency.size(); ++i)
      _fdegen.push_back(1);

    if(isharm && _anharm.size()) {
      for(int i = 0; i < _frequency.size(); ++i)
	for(int j = 0; j < _frequency.size(); ++j)
	  if(j != i)
	    _frequency[i] += double(_fdegen[j]) * _anharm(i, j) / 2.;
	  else
	    _frequency[i] += double (_fdegen[i] + 1) * _anharm(i, i);

      IO::log << IO::log_offset << "fundamentals(1/cm):";
      for(int i = 0; i < _frequency.size(); ++i)
	IO::log << "   " << _frequency[i] / Phys_const::incm;
      IO::log << "\n";
    }
  }      
  
  /************************************* CHECKING *************************************/ 

  // core
  if(!_core)
    //
    _no_run = true;
  //std::cerr << funame << "core has not been initialized\n";
  //throw Error::Init();
  //}

  // zero energy
  if(!isener && !iszero) {
    std::cerr << funame << "energy reference has not been initialized\n";
    throw Error::Init();
  }

  if(_osc_int.size() && _osc_int.size() != _frequency.size()) {
    std::cerr << funame << "numbers of vibrational frequencies and infrared intensities mismatch\n";
    throw Error::Init();
  }

  /************************************* SETTING *************************************/ 

  // electronic level ground energy correction
  if(!elevel_map.size()) {
    _elevel.resize(1, 0.);
    _edegen.resize(1, 1);
  }
  else {
    _ground += elevel_map.begin()->first;
    std::map<double, int>::const_iterator it;
    for(it = elevel_map.begin(); it != elevel_map.end(); ++it) {
      _elevel.push_back(it->first - elevel_map.begin()->first);
      _edegen.push_back(it->second);
    }
  }

  // vibrational zero-point energy correction
  if(!iszero) {
    for(int f = 0; f < _frequency.size(); ++f)
      //
      _ground += _frequency[f] / 2. * _fdegen[f];  

    // core ground energy correction
    if(_core)
      _ground += _core->ground();

    // hindered rotor energy correction
    for(int r = 0; r < _rotor.size(); ++r)
      _ground += _rotor[r]->ground();
  }

  _real_ground = _ground;

  IO::log << IO::log_offset << "ground state energy, kcal/mol:   "
	  << _real_ground / Phys_const::kcal << "\n";
  
  // tunneling correction
  //
  if(_tunnel)
    //
    _ground -= _tunnel->cutoff();
  

  // interpolating states density/number
  //
  if(mode() != NOSTATES) {
    //
    IO::Marker interpol_marker("interpolating states number/density");

    if(_emax > 0.) {
      //
      dtemp = _emax;
    }
    else
      //
      dtemp = energy_limit() - ground();
    
    itemp = (int)std::ceil(dtemp / ener_quant);

    // energy grid
    //
    Array<double> ener_grid(itemp);
    
    dtemp = 0.;
    
    for(int i = 0; i < ener_grid.size(); ++i, dtemp += ener_quant)
      //
      ener_grid[i] = dtemp;
    
    Array<double> stat_grid(itemp);

    // core states
    //
    if(_core) {
      //
      IO::Marker core_marker("core state contribution", IO::Marker::ONE_LINE);

      stat_grid[0] = 0.;
      
      dtemp = ener_quant;
      
      for(int i = 1; i < ener_grid.size(); ++i, dtemp += ener_quant)
	//
	stat_grid[i] = _core->states(dtemp);
    }
    else {
      //
      if(mode() == NUMBER) {
	//
	stat_grid = 1./ _sym_num;
      }
      else {
	//
	stat_grid = 0.;
	stat_grid[0] = 1. / ener_quant / _sym_num;
      }
    }

    Array<double> new_stat_grid(ener_grid.size());

    // electronic states contribution
    //
    if(_elevel.size() != 1) {
      //
      IO::Marker elev_marker("electronic states contribution", IO::Marker::ONE_LINE);

      new_stat_grid = 0.;
      for(int l = 0; l < _elevel.size(); ++l) {
	itemp = (int)round(_elevel[l] / ener_quant);
	for(int i = itemp; i < ener_grid.size(); ++i)
	  new_stat_grid[i] += stat_grid[i - itemp] * double(_edegen[l]);
      }
      stat_grid = new_stat_grid;
    }
    //
    else if(_edegen[0] != 1) {
      //
      IO::Marker elev_marker("electronic states contribution", IO::Marker::ONE_LINE);

      dtemp = double(_edegen[0]);
      for(int i = 1; i < ener_grid.size(); ++i)
	stat_grid[i] *= dtemp;
    }

    // vibrational frequency iteration
    //
    if(_frequency.size()) {
      //
      IO::Marker vib_marker("vibrational modes contribution", IO::Marker::ONE_LINE);

      for(int f = 0; f < _frequency.size(); ++f)
	for(int d = 0; d < _fdegen[f]; ++d) {
	  itemp = (int)round(_frequency[f] / ener_quant);
	  if(itemp < 0) {
	    std::cerr << funame << "negative frequency\n";
	    throw Error::Range();
	  }
	  for(int e = itemp; e < ener_grid.size(); ++e)
	    stat_grid[e] += stat_grid[e - itemp];
	}
    }
  

    // hindered rotor contribution
    //
    if(_rotor.size()) {
      //
      IO::Marker rotor_marker("hindered rotors contribution");
      
      for(int r = 0; r < _rotor.size(); ++r) {// 1D rotor cycle
	//
	_rotor[r]->set(energy_limit() - ground());
	
	_rotor[r]->convolute(stat_grid, ener_quant);
      }// 1D rotor cycle
    }


    // tunneling
    if(_tunnel) {
      IO::Marker tunnel_marker("tunneling contribution", IO::Marker::ONE_LINE);

      // convolute the number of states with the tunneling density
      _tunnel->convolute(stat_grid, ener_quant);
    }

    // occupation numbers for radiative transitions
    
    // ONLY WORK WITH NO DEGENERACIES
    //
    if(_osc_int.size()) {
      _occ_num.resize(_osc_int.size());
      _occ_num_der.resize(_osc_int.size());
      for(int f = 0; f < _frequency.size(); ++f) {
	new_stat_grid = stat_grid;
	itemp = (int)round(_frequency[f] / ener_quant);

	for(int e = itemp; e < ener_grid.size(); ++e)
	  new_stat_grid[e] += new_stat_grid[e - itemp];

	for(int e = itemp; e < ener_grid.size(); ++e)
	  if(stat_grid[e] != 0.) {
	    new_stat_grid[e] /= stat_grid[e];
	    new_stat_grid[e] -= 1.;
	  }
	new_stat_grid[itemp] = 0.;
	
	_occ_num[f].init(&ener_grid[itemp], &new_stat_grid[itemp], ener_grid.size() - itemp);
	dtemp = _occ_num[f].arg_max() * extra_step;
	_occ_num_der[f] = (_occ_num[f].fun_max() - _occ_num[f](_occ_num[f].arg_max() - dtemp)) / dtemp;
      }
    }

    // interpolation
    _states.init(ener_grid, stat_grid, ener_grid.size());
    dtemp = _states.fun_max() / _states(_states.arg_max() * (1. - extra_step));
    _nmax = std::log(dtemp) / std::log(1. / (1. - extra_step));

    IO::log << IO::log_offset << "effective power exponent at " 
	    << _states.arg_max() / Phys_const::kcal << " kcal/mol = "<< _nmax << "\n";

    // checking the number of states
    /*
      double tt = Phys_const::kelv * 2000.;
      double ww = 0.;
      for(int i = 1; i < ener_grid.size(); ++i) { 
      dtemp = ener_grid[i] / tt;
      if(dtemp > 100.)
      break;
      ww += stat_grid[i] * std::exp(-dtemp);      
      }
      ww *= ener_quant;
      if(mode() == NUMBER)
      ww /= tt;
      IO::log << IO::log_offset << "statistical weight at 2000K: "
      << std::setw(13) << ww 
      << std::setw(13) << weight(tt)
      << "\n";
    */

  }// interpolating states number/density

  _print();

  //graph_perturbation_theory_correction();
  
}// RRHO

// ONLY WORK WITH NO DEGENERACIES
//
void Model::RRHO::_init_graphex (std::istream& from)
{
  const char funame [] = "Model::RRHO::_init_graphex: ";

  int    itemp;
  double dtemp;
  
  KeyGroup InitGraphexGroup;

  Key     potex_key("PotentialExpansion[1/cm]"        );
  Key      bond_key("BondNumberMax"                   );
  Key      ftol_key("FrequencyTolerance"              );
  Key      keep_key("KeepPermutedGraphs"              );
  Key      scut_key("FourierSumCutoff"                );
  Key      redt_key("ReductionThreshold"              );

  // potential expansion
  //
  typedef std::map<std::multiset<int>, double> potex_t;
  //
  potex_t potex;
  
  std::string token, comment, line, name;
  //
  while(from >> token) {
    //
    // end of input
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);

      break;
    }
    // potential expansion
    //
    else if(potex_key == token) {
      //
      if(!_frequency.size()) { 
	//
	ErrOut err_out;

	err_out << funame << token << ": frequencies should be initialized before graph expansion initialization\n";
      }
      
      std::getline(from, comment);
      
      Graph::read_potex(_frequency, from, potex);
    }
    // maximum number of graph bonds (edges)
    //
    else if(bond_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      if(itemp <= 0) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< itemp;
      }

      Graph::bond_max = itemp;
      
      std::getline(from, comment);
    }
    // frequency reduction tolerance
    //
    else if(ftol_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      if(dtemp > 1. || dtemp < 0.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: " << dtemp;
      }

      Graph::Expansion::freq_tol = dtemp;

      std::getline(from, comment);
    }	
    // keep permuted graphs in the databases
    //
    else if(keep_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << "corrupted";
      }

      if(itemp) {
	//
	Graph::Expansion::mod_flag |= Graph::Expansion::KEEP_PERM;
      }
      else if(Graph::Expansion::mod_flag & Graph::Expansion::KEEP_PERM) {
	//
	Graph::Expansion::mod_flag ^= Graph::Expansion::KEEP_PERM;
      }

      std::getline(from, comment);
    }
    // low temperature / high frequency graph reduction threshold
    //
    else if(redt_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      if(dtemp <= 1.) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< dtemp;
      }
      
      Graph::FreqGraph::red_thresh = dtemp;

      std::getline(from, comment);
    }	
    // fourier sum graph evaluation cutoff
    //
    else if(scut_key == token) {
      //
      if(!(from >> itemp)) {
	//
	ErrOut err_out;

	err_out << funame << token << ": corrupted";
      }

      if(itemp < 2) {
	//
	ErrOut err_out;

	err_out << funame << token << ": out of range: "<< itemp;
      }

      Graph::FreqGraph::four_cut = itemp;

      std::getline(from, comment);
    }	
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      // 
      ErrOut err_out;

      err_out << funame << "unknown keyword: " << token << "\n" << Key::show_all() << "\n";
    }
  }

  if(!potex.size()) {
    //
    ErrOut err_out;

    err_out << funame << "potential expansion not initialized\n";
  }
  
  itemp = 0;
  //
  for(potex_t::const_iterator pit = potex.begin(); pit != potex.end(); ++pit)
    //
    if(pit->first.size() > itemp)
      //
      itemp = pit->first.size();

  Graph::potex_max = itemp;

  Graph::init();

  _graphex.init(_frequency, potex);
}

Model::RRHO::~RRHO ()
{
  //std::cout << "Model::RRHO destroyed\n";
}

double Model::RRHO::states (double ener) const
{
  const char funame [] = "Model::RRHO::states: ";
  
  ener -= ground();

  if(ener <= 0.)
    return 0.;

  if(ener >= _states.arg_max()) {
    if(ener >= _states.arg_max() * 2.)
      IO::log << IO::log_offset << funame << "WARNING: energy far beyond interpolation range\n";
    return _states.fun_max() * std::pow(ener / _states.arg_max(), _nmax);
  }

  return _states(ener);
}

double Model::RRHO::weight (double temperature) const
{
  double dtemp;
  int    itemp;

  // electronic level contribution
  double res = 0.;
  for(int l = 0; l < _elevel.size(); ++l)
    res += std::exp(-_elevel[l] / temperature) * double(_edegen[l]);

  // core contribution
  if(_core)
    res *= _core->weight(temperature);
  else
    res /= _sym_num;

  // vibrational contribution
  std::vector<double> rfactor(_frequency.size()), sfactor(_frequency.size());
  
  for(int  f = 0; f < _frequency.size(); ++f) {
    rfactor[f] = std::exp(- _frequency[f] / temperature);
    sfactor[f] = 1. / (1. - rfactor[f]);
    for(int i = 0; i < _fdegen[f]; ++i)
      res *= sfactor[f];
  }

  double acl;
  
  // vibrational anharmonic contributions
  //
  if(_anharm.isinit()) {
    //
    // first order correction
    //
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += _anharm(i, i) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i] * sfactor[i], 2.);

      for(int j = 0; j < i; ++j)
	acl += _anharm(i, j) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * sfactor[i] * rfactor[j] * sfactor[j];
    }
    acl /= temperature;
    
    //std::cerr << funame << "temperature = " << temperature / Phys_const::kelv << "   acl1 = " << acl << std::endl;

    res *= std::exp(-acl);
    
    // second order correction
    //
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += 2. * std::pow(_anharm(i, i), 2.) * double(_fdegen[i] * (_fdegen[i] + 1))	* std::pow(rfactor[i], 2.)
	* (1. + double(2 * (_fdegen[i] + 1)) * rfactor[i]) * std::pow(sfactor[i], 4.);
    
      for(int j = 0; j < i; ++j)
	acl += std::pow(_anharm(i, j), 2.) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j]
	  * (1. + (double)_fdegen[i] * rfactor[i] + (double)_fdegen[j] * rfactor[j])
	  * std::pow(sfactor[i] * sfactor[j], 2.);
    
      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  if(i != j && i != k)
	    acl += 2. * _anharm(i, j) * _anharm(i, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k])
	      * rfactor[i] * rfactor[j] * rfactor[k] * std::pow(sfactor[i], 2.) * sfactor[j] * sfactor[k];
	
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 4. * _anharm(i, i) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) *_fdegen[j])
	    * std::pow(rfactor[i], 2.) * rfactor[j] * std::pow(sfactor[i], 3.) * sfactor[j];
    }
    acl /= 2. * temperature * temperature;

    //std::cerr << funame << "temperature = " << temperature / Phys_const::kelv << "   acl2 = " << acl << std::endl;

    res *= std::exp(acl);
    
    // third order correction
    //
    acl = 0.;
    for(int i = 0; i < _frequency.size(); ++i) {
      acl += 4. * std::pow(_anharm(i, i), 3.) * double(_fdegen[i] * (_fdegen[i] + 1)) * std::pow(rfactor[i], 2.) *
	(1. + double(8 * _fdegen[i] + 12) * rfactor[i] + double(8 * _fdegen[i] * _fdegen[i] + 22 * _fdegen[i] + 15) * std::pow(rfactor[i], 2.) +
	 double(2 * _fdegen[i] * _fdegen[i] + 4 * _fdegen[i] + 2) * std::pow(rfactor[i], 3.)) * std::pow(sfactor[i], 6.);

      for(int j = 0; j < i; ++j)
	acl += std::pow(_anharm(i, j), 3.) * double(_fdegen[i] * _fdegen[j]) * rfactor[i] * rfactor[j] *
	  (1. + double(3 * _fdegen[i] + 1) * rfactor[i] + double(3 * _fdegen[j] + 1) * rfactor[j] +
	   std::pow(double(_fdegen[i]) * rfactor[i], 2.) * (1. + rfactor[j]) +
	   std::pow(double(_fdegen[j]) * rfactor[j], 2.) * (1. + rfactor[i]) +
	   double(6 * _fdegen[i] * _fdegen[j] + 3 * _fdegen[i] + 3 * _fdegen[j] + 1) * rfactor[i] * rfactor[j]) *
	  std::pow(sfactor[i] * sfactor[j], 3.);
      
      for(int j = 0; j < _frequency.size(); ++j)
	if(i != j)
	  acl += 12. * std::pow(_anharm(i, i), 2.) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j]) *
	    std::pow(rfactor[i], 2.) * rfactor[j] *
	    (1. + double(3 * _fdegen[i] + 4) * rfactor[i] + double(_fdegen[i] + 1) * std::pow(rfactor[i], 2.)) *
	    std::pow(sfactor[i], 5.) * sfactor[j]
	    + 6. * _anharm(i, i) * std::pow(_anharm(i, j), 2.) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j]) *
	    std::pow(rfactor[i], 2.) * rfactor[j] *
	    (2. + double(2 * _fdegen[i] + 1) * rfactor[i] + double(2 * _fdegen[j]) * rfactor[j] + double(_fdegen[j]) * rfactor[i] * rfactor[j]) *
	    std::pow(sfactor[i], 4.) * std::pow(sfactor[j], 2.);

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  if(i != j && i != k && k != j)
	    acl += 3. * std::pow(_anharm(i, j), 2.) * _anharm(i, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) *
	      rfactor[i] * rfactor[j] * rfactor[k] *
	      (1. + double(2 * _fdegen[i] + 1) * rfactor[i] + double(_fdegen[j]) * rfactor[j] * (1. + rfactor[i])) *
	      std::pow(sfactor[i] * sfactor[j], 2.) * sfactor[k];

      for(int j = 0; j < i; ++j)
	acl += 24. * _anharm(i, i) * _anharm(j, j) * _anharm(i, j) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j] * (_fdegen[j] + 1)) *
	  std::pow(rfactor[i] * rfactor[j], 2.) * std::pow(sfactor[i] * sfactor[j], 3.);

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  if(i != j && i != k)
	    acl += 12. * _anharm(i, i) * _anharm(i, j) * _anharm(i, k) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j] * _fdegen[k]) *
	      std::pow(rfactor[i], 2.) * rfactor[j] * rfactor[k] * (2. + rfactor[i]) * std::pow(sfactor[i], 4.) * sfactor[j] * sfactor[k];

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  if(i != j && i != k && k != j)
	    acl += 12. * _anharm(i, i) * _anharm(i, j) * _anharm(j, k) * double(_fdegen[i] * (_fdegen[i] + 1) * _fdegen[j] * _fdegen[k]) *
	      std::pow(rfactor[i], 2.) * rfactor[j] * rfactor[k] * std::pow(sfactor[i], 3.) * std::pow(sfactor[j], 2.) * sfactor[k];

      for(int j = 0; j < _frequency.size(); ++j)
	for(int k = 0; k < j; ++k)
	  for(int l = 0; l < k; ++l)
	    if(i != j && i != k && i != l)
	      acl += 6. * _anharm(i, j) * _anharm(i, k) * _anharm(i, l) * double(_fdegen[i] * _fdegen[j] * _fdegen[k] * _fdegen[l]) *
		rfactor[i] * rfactor[j] * rfactor[k] * rfactor[l] * (1. + rfactor[i]) *
		std::pow(sfactor[i], 3.) * sfactor[j] * sfactor[k] * sfactor[l];

      // takes into account i->j, k->l symmetry
      for(int j = 0; j < i; ++j)
	for(int k = 0; k < _frequency.size(); ++k)
	  for(int l = 0; l < _frequency.size(); ++l)
	    if(i != k && i != l && j != k && j != l && k != l)
	      acl += 6. * _anharm(i, j) * _anharm(i, k) * _anharm(j, l) * double(_fdegen[i] * _fdegen[j] * _fdegen[k] * _fdegen[l]) *
		rfactor[i] * rfactor[j] * rfactor[k] * rfactor[l] * std::pow(sfactor[i] * sfactor[j], 2.) * sfactor[k] * sfactor[l];

      for(int j = 0; j < i; ++j)
	for(int k = 0; k < j; ++k)
	  acl += 6. * _anharm(i, j) * _anharm(i, k) * _anharm(j, k) * double(_fdegen[i] * _fdegen[j] * _fdegen[k]) *
	    rfactor[i] * rfactor[j] * rfactor[k] * (1. + _fdegen[i] * rfactor[i] + _fdegen[j] * rfactor[j] + _fdegen[k] * rfactor[k]) *
	    std::pow(sfactor[i] * sfactor[j] * sfactor[k], 2.);
    }
    acl /= 6. * std::pow(temperature, 3.);

    //std::cerr << funame << "temperature = " << temperature / Phys_const::kelv << "   acl3 = " << acl << std::endl;
												       
    res *= std::exp(-acl);
  }

  // hindered 1D rotor contribution
  for(int r = 0; r < _rotor.size(); ++r)
    res *= _rotor[r]->weight(temperature);

  // tunneling contribution
  if(_tunnel)
    res *= _tunnel->weight(temperature);

  return res;
}

double Model::RRHO::tunnel_weight (double temperature) const 
{
  if(_tunnel)
    return _tunnel->weight(temperature);
 
  return 1.; 
}

// radiational transitions
int Model::RRHO::oscillator_size () const
{
  return _osc_int.size(); 
}

double Model::RRHO::oscillator_frequency (int f) const
{
  const char funame [] = "Model::RRHO::oscillator_frequency: ";
  
  if(f < 0 || f >= _osc_int.size()) {
    std::cerr << funame << "oscillator index out of range\n";
    throw Error::Range();
  }
  
  return _frequency[f];
}

double Model::RRHO::infrared_intensity (double ener, int f) const
{
  const char funame [] = "Model::RRHO::occupation_number: ";
  
  if(f < 0 || f >= _osc_int.size()) {
    std::cerr << funame << "oscillator index out of range\n";
    throw Error::Range();
  }

  if(_osc_int[f] == 0.)
    return 0.;

  ener -= ground();

  if(ener <= _occ_num[f].arg_min())
    return 0.;

  double res = _osc_int[f];

  if(ener >= _occ_num[f].arg_max()) {

    if(ener >= _occ_num[f].arg_max() * 2.)
      IO::log << IO::log_offset << funame << "WARNING: energy far beyond interpolation range\n";

    res *= _occ_num[f].fun_max() + (ener - _occ_num[f].arg_max()) * _occ_num_der[f];
  }
  else
    res *= _occ_num[f](ener);

  return res;
}

/*************************************************************************************
 *********************************** UNION OF WELLS **********************************
 *************************************************************************************/

Model::UnionSpecies::UnionSpecies (IO::KeyBufferStream& from, const std::string& n, int m) 
  : Species(n, m)
{
  const char funame [] = "Model::UnionSpecies::UnionSpecies: ";

  IO::Marker funame_marker(funame);

  SharedPointer<Species> p;
  while(p = new_species(from, n, mode()))
    _species.push_back(p);
  
  if(!_species.size()) {
    std::cerr << funame << "no species found\n";
    throw Error::Init();
  }

  for(int i = 0; i < _species.size(); ++i) {
    try {
      _mass = _species[i]->mass();
    }
    catch(Error::General) {
      IO::log << IO::log_offset << funame << "WARNING: mass of the " << i + 1 << "-th species not defined\n";
      continue;
    }
    break;
  }

  for(_Cit w = _species.begin(); w != _species.end(); ++w)
    if(_species.begin() == w || (*w)->ground() < _ground)
      _ground = (*w)->ground();

  for(_Cit w = _species.begin(); w != _species.end(); ++w)
    if(_species.begin() == w || (*w)->real_ground() < _real_ground)
      _real_ground = (*w)->real_ground();

  // radiational transitions
  if(mode() == DENSITY)
    for(int w = 0; w < _species.size(); ++w) {
      _osc_shift.push_back(_osc_spec_index.size());
      for(int i = 0; i < _species[w]->oscillator_size(); ++i)
	_osc_spec_index.push_back(w);
    }

  _print();
}// Union Species

Model::UnionSpecies::~UnionSpecies ()
{
  //std::cout << "Model::UnionSpecies destroyed\n";
}

double Model::UnionSpecies::states (double energy) const 
{
  if(energy <= _ground)
    return 0.;

  double res = 0.;
  for(_Cit w = _species.begin(); w != _species.end(); ++w)
    res += (*w)->states(energy);
  return res;
}

double Model::UnionSpecies::weight (double temperature) const
{
  double res = 0.;
  for(_Cit w = _species.begin(); w != _species.end(); ++w)
    res += (*w)->weight(temperature) * std::exp((ground() - (*w)->ground()) / temperature);
  return res;
}

void Model::UnionSpecies::shift_ground (double e)
{
  _ground += e;
  _real_ground += e;
  for(_Cit w = _species.begin(); w != _species.end(); ++w)
    (*w)->shift_ground(e);
}

// radiational transitions
double Model::UnionSpecies::oscillator_frequency (int num) const
{
  const char funame [] = "Model::UnionSpecies::oscillator_frequency: ";

  if(num >= _osc_spec_index.size() || num < 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
    
  int itemp = _osc_spec_index[num];

  return _species[itemp]->oscillator_frequency(num - _osc_shift[itemp]);
}

double  Model::UnionSpecies::infrared_intensity (double ener, int num) const
{
  const char funame [] = "Model::UnionSpecies::infrared_intensity: ";

  if(num >= _osc_spec_index.size() || num < 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
    
  double dtemp = states(ener);
  if(dtemp <= 0.)
    return 0.;

  int itemp = _osc_spec_index[num];

  return _species[itemp]->infrared_intensity(ener, num - _osc_shift[itemp]) 
    * _species[itemp]->states(ener) / dtemp;
}

int Model::UnionSpecies::oscillator_size () const
{
  return _osc_spec_index.size();
}

/***********************************************************************************************
 ********************************* VARIATIONAL BARRIER MODEL ***********************************
 ***********************************************************************************************/

Model::VarBarrier::VarBarrier(IO::KeyBufferStream& from, const std::string& n) 
  : Species(n, NUMBER), _tunnel(0), _ener_quant(Phys_const::incm), _emax(-1.), _tts_method(STATISTICAL)
{
  const char funame [] = "Model::VarBarrier::VarBarrier: ";

  IO::Marker funame_marker(funame);

  KeyGroup VariationalBarrierModel;

  Key estep_key("InterpolationEnergyStep[1/cm]"   );
  Key  emax_key("InterpolationEnergyMax[kcal/mol]");
  Key extra_key("ExtrapolationStep"               );
  Key  tunn_key("Tunneling"                       );
  Key  rrho_key("RRHO"                            );
  Key outer_key("OuterRRHO"                       );
  Key   tts_key("2TSMethod"                       );
  Key  tmin_key("OutputTemperatureMin[K]"         );
  Key  tmax_key("OutputTemperatureMax[K]"         );
  Key tstep_key("OutputTemperatureStep[K]"        );


  int tmin = -1, tmax = 2000, tstep = 100;

  double extra_step = 0.1;

  int    itemp;
  double dtemp;

  std::string token, comment, stemp;
  while(from >> token) {
    // end input
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // tunneling
    else if(tunn_key == token) {
      if(_tunnel) {
	std::cerr << funame << token << ": already intialized\n";
	throw Error::Init();
      }
      _tunnel = new_tunnel(from);
    }
    // RRHO
    else if(rrho_key == token) {
      std::getline(from, comment);
      IO::log << IO::log_offset << _rrho.size() + 1 << "-th RRHO:\n";
      _rrho.push_back(SharedPointer<RRHO>(new RRHO(from, n, NUMBER)));
    }
    // outer barrier
    else if(outer_key == token) {
      if(_outer) {
	std::cerr << funame << token << ": already intialized\n";
	throw Error::Init();
      }
      std::getline(from, comment);
      IO::log << IO::log_offset << "Outer RRHO:\n";
      _outer = SharedPointer<RRHO>(new RRHO(from, n, NUMBER));
    }
    //  two-transition-states model calculation method
    else if(tts_key == token) {
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(stemp == "statistical")
	_tts_method = STATISTICAL;
      else if(stemp == "dynamical")
	_tts_method = DYNAMICAL;
      else {
	std::cerr << funame << token << ": unknown method: " << stemp << "\n";
	throw Error::Input();
      }
    }
    //  energy step
    else if(estep_key == token) {
      if(!(from >> _ener_quant)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_ener_quant <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      _ener_quant *= Phys_const::incm;
    }
    // extrapolation step
    else if(extra_key == token) {
      if(!(from >> extra_step)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(extra_step <= 0. || extra_step >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // interpolation maximal energy
    else if(emax_key == token) {
      if(!(from >> _emax)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(_emax <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }
      _emax *= Phys_const::kcal;
    }
    // output temperature min
    else if(tmin_key == token) {
      if(!(from >> tmin)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(tmin <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // output temperature max
    else if(tmax_key == token) {
      if(!(from >> tmax)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(tmax <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // output temperature step
    else if(tstep_key == token) {
      if(!(from >> tstep)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(tstep <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
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
 
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  // rrho checking
  if(!_rrho.size()) {
    std::cerr << funame << "RRHOs not initialized\n";
    throw Error::Init();
  }

  // check for tunneling inside of RRHOs
  for(int i = 0; i < _rrho.size(); ++i)
    if(_rrho[i]->istunnel()) {
      std::cerr << funame << "tunneling inside " << i + 1 << "-th RRHO\n";
      throw Error::Init();
    }

  if(_outer && _outer->istunnel()) {
    std::cerr << funame << "tunneling inside " << "outer RRHO\n";
    throw Error::Init();
  }

  // mass
  for(int i = 0; i < _rrho.size(); ++i) {
    try { 
      _mass = _rrho[i]->mass();
    } catch(Error::General) {
      IO::log << IO::log_offset << funame << "WARNING: mass of " << i + 1 << "-th configuration not defined\n"; 
      continue;
    }
    break;
  }

  // ground setting
  for(int v = 0; v < _rrho.size(); ++v)
    if(!v || _rrho[v]->ground() > _ground)
      _ground = _rrho[v]->ground();

  // interpolating states density/number
  {
    IO::Marker interpol_marker("interpolating states number/density", IO::Marker::ONE_LINE);

    if(_emax > 0.)
      dtemp = _emax;
    else
      dtemp = energy_limit() - _ground;
    itemp = (int)std::ceil(dtemp / _ener_quant);

    Array<double> ener_grid(itemp);
    _stat_grid.resize(itemp);
    ener_grid[0]  = 0.;
    _stat_grid[0] = 0.;

    double ener = _ener_quant;
    for(int i = 1; i < ener_grid.size(); ++i, ener += _ener_quant) {
      // energy grid
      ener_grid[i] = ener;
      // density / number of states
      for(int v = 0; v < _rrho.size(); ++v) {
	dtemp = _rrho[v]->states(ener + _ground);
	if(!v || dtemp < _stat_grid[i])
	  _stat_grid[i] = dtemp;
      }
    }

    _real_ground = _ground;

    // tunneling
    if(_tunnel) {
      _ground -= _tunnel->cutoff();
      // convolute number of states with tunneling density
      _tunnel->convolute(_stat_grid, _ener_quant);
    }

    // two-transition-states model
    if(_outer) {
      int iground = 0;
      dtemp = _outer->ground() - _ground;

      // reset energy grid
      if(dtemp > 0.) {
	iground = (int)std::floor(dtemp / _ener_quant);
	dtemp -= (double)iground * _ener_quant;
	ener_grid.resize(ener_grid.size() - iground);

	for(int i = 1; i < ener_grid.size(); ++i)
	  ener_grid[i] -= dtemp;

	_ground = _outer->ground();
      }

      if(_outer->ground() > _real_ground)
	_real_ground = _outer->ground();

      // two-transition-states model number of states
      switch(_tts_method) {
      case DYNAMICAL:
	for(int i = 1; i < ener_grid.size(); ++i) {
	  double inn_stat = _stat_grid[i + iground];
	  double out_stat = _outer->states(ener_grid[i] + _ground);
	  if(out_stat < inn_stat)
	    _stat_grid[i] = out_stat;
	  else
	    _stat_grid[i] = inn_stat;
	}
	break;

      case STATISTICAL:
	for(int i = 1; i < ener_grid.size(); ++i) {
	  double inn_stat = _stat_grid[i + iground];
	  double out_stat = _outer->states(ener_grid[i] + _ground);
	  if(inn_stat != 0. && out_stat != 0.)
	    _stat_grid[i] = inn_stat * out_stat / (inn_stat + out_stat);
	  else
	    _stat_grid[i] = 0.;
	}
	break;

      default:
	std::cerr << funame << "wrong two-transition-states model case\n";
	throw Error::Logic();
      }

      _stat_grid.resize(ener_grid.size());
    }

    // interpolation
    _states.init(ener_grid, _stat_grid, ener_grid.size());

    // extrapolation
    dtemp = _states.fun_max() / _states(_states.arg_max() * (1. - extra_step));
    _nmax = std::log(dtemp) / std::log(1. / (1. - extra_step));
  }  

  IO::log << IO::log_offset << "effective power exponent at " 
	  << _states.arg_max() / Phys_const::kcal << " kcal/mol = " << _nmax << "\n";

  IO::log << IO::log_offset << "Constituent transition states partition functions\n";
  if(tmin < 0)
    tmin = tstep;
  IO::log << IO::log_offset << std::setw(5) << "T, K";
  for(int v = 0; v <= _rrho.size(); ++v)
    if(v < _rrho.size())
      IO::log << std::setw(13) << v + 1;
    else
      IO::log << std::setw(13) << "Min";

  IO::log << "\n";

  for(int t = tmin; t <= tmax; t += tstep) {
    IO::log << IO::log_offset << std::setw(5) << t;
    double tval = (double)t * Phys_const::kelv;
    for(int v = 0; v <= _rrho.size(); ++v)
      if(v < _rrho.size())
	IO::log << std::setw(13) << _rrho[v]->weight(tval) * std::exp((_ground - _rrho[v]->ground()) / tval);
      else {
	dtemp = weight(tval);
	if(_tunnel)
	  dtemp /= _tunnel->weight(tval) * std::exp(_tunnel->cutoff() / tval);
	IO::log << std::setw(13) << dtemp;
      }
    IO::log << "\n";
  }
  IO::log << "\n";

  _print();
}// Variational barrier

Model::VarBarrier::~VarBarrier ()
{
  //std::cout << "Model::VarBarrier destroyed\n";
}

double Model::VarBarrier::states (double ener) const
{
  const char funame [] = "Model::VarBarrier::state_number: ";

  ener -= _ground;

  if(ener <= 0.)
    return 0.;

  if(ener >= _states.arg_max()) {

    if(ener > _states.arg_max() * 2.)
      IO::log << IO::log_offset << funame << "WARNING: energy far beyond interpolation range\n";

    return _states.fun_max() * std::pow(ener / _states.arg_max(), _nmax);
  }

  return _states(ener);
}

double Model::VarBarrier::weight (double temperature) const
{
  const char funame [] = "Model::VarBarrier::weight: ";

  static const double eps = 1.e-3;

  double dtemp;
  
  double res = 0.;
  double ener = _ener_quant;
  for(int i = 1; i < _stat_grid.size(); ++i, ener += _ener_quant)
    res += _stat_grid[i] / std::exp(ener / temperature);
  res *= _ener_quant / temperature;

  dtemp = _states.arg_max() / temperature;
  if(dtemp <= _nmax) {
    IO::log << IO::log_offset << funame 
	    << "WARNING: integration cutoff energy is less than the distribution maximum energy\n";
    return res;
  }

  dtemp = _states.fun_max() / std::exp(dtemp) / (1. - _nmax / dtemp);
  if(dtemp / res > eps)
    IO::log << IO::log_offset << funame << "WARNING: integration cutoff error = " << dtemp / res << "\n";
  res += dtemp;
  return res;
}

double Model::VarBarrier::tunnel_weight (double temperature) const 
{
  if(_tunnel)
    return _tunnel->weight(temperature);
 
  return 1.; 
}

/********************************************************************************************
 ************************************* ATOMIC FRAGMENT **************************************
 ********************************************************************************************/

Model::AtomicSpecies::AtomicSpecies (IO::KeyBufferStream& from, const std::string& n) 
  : Species(n, NOSTATES)
{
  const char funame [] = "Model::AtomicFragment::AtomicFragment: ";

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE | IO::Marker::NOTIME);

  int    itemp;
  double dtemp;

  KeyGroup AtomicSpeciesModel;

  Key mass_key("Mass[amu]"                      );
  Key name_key("Name"                           );
  Key incm_elev_key("ElectronicLevels[1/cm]"    );
  Key kcal_elev_key("ElectronicLevels[kcal/mol]");
  Key   kj_elev_key("ElectronicLevels[kJ/mol]"  );
  Key   ev_elev_key("ElectronicLevels[eV]"      );
  Key   au_elev_key("ElectronicLevels[au]"      );

  std::map<double, int> elevel_map;

  std::string token, comment;
  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // mass
    else if(mass_key == token) {
      if(_mass > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      
      if(!(from >> _mass)) {
	std::cerr << funame << token << ": unreadable\n";
	throw Error::Input();
      }
      std::getline(from, comment);
      
      if(_mass <= 0.) {
	std::cerr << funame << token << ": should be positive\n";
	throw Error::Range();
      }

      _mass *= Phys_const::amu;
    }
    // atom by name
    else if(name_key == token) {
      if(_mass > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      _mass = Atom(from).mass();
    }      
    // electronic energy levels & degeneracies
    else if(incm_elev_key == token || au_elev_key == token ||
	    kcal_elev_key == token || kj_elev_key == token || ev_elev_key == token) {
      int num;
      if(!(from >> num)) {
	std::cerr << funame << token << ": levels number unreadable\n";
	throw Error::Input();
      }

      if(num < 1) {
	std::cerr << funame << token << ": levels number should be positive\n";
	throw Error::Range();
      }

      std::getline(from, comment);

      for(int l = 0; l < num; ++l) {
	IO::LineInput level_input(from);
	if(!(level_input >> dtemp >> itemp)) {
	  std::cerr << funame << token << ": format: energy[1/cm] degeneracy(>=1)\n";
	  throw Error::Input();
	}

	if(incm_elev_key == token)
	  dtemp *= Phys_const::incm;
	if(kcal_elev_key == token)
	  dtemp *= Phys_const::kcal;
	if(kj_elev_key == token)
	  dtemp *= Phys_const::kjoul;
	if(ev_elev_key == token)
	  dtemp *= Phys_const::ev;

	if(elevel_map.find(dtemp) != elevel_map.end()) {
	  std::cerr << funame << token << ": identical energy levels\n";
	  throw Error::Range();
	}

	if(itemp < 1) {
	  std::cerr << funame << token << ": degeneracy should be positive\n";
	  throw Error::Range();
	}

	elevel_map[dtemp] = itemp;
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
 
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  if(_mass < 0.) {
    std::cerr << funame << "mass is not initialized\n";
    throw Error::Init();
  }

  if(!elevel_map.size()) {
    _ground = 0.;
    _elevel.resize(1, 0.);
    _edegen.resize(1, 1);
  }
  else {
    _ground = elevel_map.begin()->first;
    std::map<double, int>::const_iterator it;
    for(it = elevel_map.begin(); it != elevel_map.end(); ++it) {
      _elevel.push_back(it->first - _ground);
      _edegen.push_back(it->second);
    }
  }

  _print();
}

Model::AtomicSpecies::~AtomicSpecies ()
{
  //std::cout << "Model::AtomicSpecies destroyed\n";
}

double Model::AtomicSpecies::weight (double temperature) const
{
  double res = 0.;
  
  for(int l = 0; l < _elevel.size(); ++l)
    //
    res += std::exp(-_elevel[l] / temperature) * double(_edegen[l]);

  return res;
}

double Model::AtomicSpecies::states (double) const
{
  const char funame [] = "Model::AtomicSpecies::states: ";

  return 0.;
}

void Model::AtomicSpecies::shift_ground (double e)
{
  const char funame [] = "Model::AtomicSpecies::shift_ground: ";

  _ground += e;
}

/********************************************************************************************
 ************************************ BIMOLECULAR MODEL *************************************
 ********************************************************************************************/

Model::Bimolecular::Bimolecular(IO::KeyBufferStream& from, const std::string& n) 
  : _name(n), _dummy(false)
{
  const char funame [] = "Model::Bimolecular::Bimolecular: ";

  IO::Marker funame_marker(funame);

  KeyGroup BimolecularModel;

  Key incm_ener_key("GroundEnergy[1/cm]"    );
  Key kcal_ener_key("GroundEnergy[kcal/mol]");
  Key   au_ener_key("GroundEnergy[au]"      );
  Key   kj_ener_key("GroundEnergy[kJ/mol]"  );
  Key      frag_key("Fragment"              );
  Key     dummy_key("Dummy"                 );
  Key      dens_key("StatesDensity"         );
  
  int    itemp;
  double dtemp;

  bool isener = false;
  bool isdensity = false;

  std::string token, comment, stemp;
  while(from >> token) {
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // total energy
    else if(incm_ener_key == token || au_ener_key == token ||
	    kcal_ener_key == token || kj_ener_key == token) {
      if(isener) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      isener = true;

      if(!(from >> _ground)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(incm_ener_key == token)
	_ground *= Phys_const::incm;
      if(kcal_ener_key == token)
	_ground *= Phys_const::kcal;
      if(kj_ener_key == token)
	_ground *= Phys_const::kjoul;      
    }
    // set density of states
    else if(dens_key == token) {
      if(_fragment.size()) {
	std::cerr << funame << token << "should be before any fragment definition\n";
	throw Error::Init();
      }
      std::getline(from, comment);
      isdensity = true;
    }
    // fragment
    else if(frag_key == token) {
      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }      
      std::getline(from, comment);
      if(isdensity)
	_fragment.push_back(new_species(from, stemp, DENSITY));
      else
	_fragment.push_back(new_species(from, stemp, NOSTATES));
    }
    // dummy
    else if(dummy_key == token) {
      if(_fragment.size()) {
	std::cerr << funame << token << "should be before any fragment definition\n";
	throw Error::Init();
      }
      _dummy = true;
      std::getline(from, comment);
      return;
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
 
  if(!from) {
    std::cerr << funame << "input stream corrupted\n";
    throw Error::Input();
  }

  /************************************* CHECKING *************************************/ 

  if(!isener) {
    std::cerr << funame << "ground energy has not been initialized\n";
    throw Error::Init();
  }

  if(_fragment.size() != 2) {
    std::cerr << funame << "wrong number of fragments\n";
    throw Error::Init();
  }

  /************************************* SETTING *************************************/ 

  // translational partition function & ground energy
  dtemp = 0.;
  for(int i = 0; i < 2; ++i) {
    dtemp += 1. / _fragment[i]->mass();
    _ground += _fragment[i]->ground();
    _fragment[i]->shift_ground(-_fragment[i]->ground());
  }
  
  dtemp *= 2. * M_PI;
  _weight_fac = 1. / dtemp / std::sqrt(dtemp);

}// Bimolecular

Model::Bimolecular::~Bimolecular ()
{
  //std::cout << "Model::Bimolecular destroyed\n";
}

double Model::Bimolecular::fragment_weight (int i, double temperature) const 
{ 
  double dtemp = _fragment[i]->mass() * temperature / 2. / M_PI;

  return dtemp * std::sqrt(dtemp) * _fragment[i]->weight(temperature);
}

double Model::Bimolecular::ground () const 
{ 
  if(_dummy)
    return 0.;

  return _ground;
}

double Model::Bimolecular::weight (double temperature) const
{
  if(_dummy)
    return -1.;

  double res = _weight_fac * temperature * std::sqrt(temperature);
  for(int i = 0; i < 2; ++i)
    res *= _fragment[i]->weight(temperature);

  return res;
}

void Model::Bimolecular::shift_ground (double e)
{
  _ground += e;
}

/********************************************************************************************
 ************************************** ESCAPE MODEL ****************************************
 ********************************************************************************************/

Model::ConstEscape::ConstEscape(IO::KeyBufferStream& from) 
  : _rate(-1.)
{
  const char funame [] = "Model::ConstEscape::ConstEscape: ";

  KeyGroup ConstEscapeModel;
  Key first_key("PseudoFirstOrderRateConstant[1/sec]");

  double dtemp;
  int    itemp;

  std::string token, comment;
  
  // actual input
  //
  while(from >> token) {
    //
    // end key
    //
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // rate value
    //
    if(first_key == token) {
      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      
      std::getline(from, comment);
      if(dtemp <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      _rate = dtemp * Phys_const::herz;
    }
    // unknown key
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }
  
  if(!from) {
    std::cerr << funame << "stream is corrupted\n";
    throw Error::Input(); 
  }

  if(_rate <= 0.) {
    //
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
}

Model::FitEscape::FitEscape(IO::KeyBufferStream& from) 
  : _ground(0.)
{
  const char funame [] = "Model::FitEscape::FitEscape: ";

  IO::Marker funame_marker(funame, IO::Marker::ONE_LINE | IO::Marker::NOTIME);

  int         itemp;
  double      dtemp;
  std::string stemp;

  // input parameters
  std::ifstream rate_in;

  // input keys
  KeyGroup FitEscapeModel;
  Key file_key("RateDataFile[kcal/mol,1/sec]");

  // actual input
  std::string comment, token;
  while(from >> token) {
    //
    // end key
    //
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // rate data file
    else if(file_key == token) {

      if(rate_in.is_open()) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      std::getline(from, comment);

      rate_in.open(stemp.c_str());
      if(!rate_in) {
	std::cerr << funame << token << ": cannot open  " << stemp << " file\n";
	throw Error::Open();
      }
    }
    // unknown key
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  if(!from) {
    std::cerr << funame << "corrupted\n";
    throw Error::Input(); 
  }

  if(!rate_in) {
    std::cerr << funame << "rate data not initialized\n";
    throw Error::Init();
  }

  std::map<double, double> rate_data;

  double eval, rval;
  while(rate_in >> eval) {
    eval *= Phys_const::kcal; // energy

    if(!(rate_in >> rval)) {
      std::cerr << funame << "reading rate data failed\n";
      throw Error::Input();
    }
    
    if(rval < 0.) {
      std::cerr << funame << "negative rate\n";
      throw Error::Range();
    }

    rate_data[eval] = Phys_const::herz * rval;
  }

  if(!rate_data.size()) {
    std::cerr << funame << "not initialized";
    throw Error::Init();
  }
  
  Array<double> x((int)rate_data.size());
  Array<double> y((int)rate_data.size());

  itemp = 0;
  for(std::map<double, double>::const_iterator it = rate_data.begin(); it != rate_data.end(); ++it, ++itemp) {
    x[itemp] = it->first;
    y[itemp] = it->second;
  }
  
  _rate.init(x, y, rate_data.size());
}

double Model::FitEscape::rate (double ener) const
{
  ener -= _ground;

  double dtemp;
  
  if(ener >= _rate.arg_max()) {
    //
    dtemp = _rate.fun_max();
  }
  else if(ener <= _rate.arg_min()) {
    //
    dtemp = _rate.fun_min();
  }
  else {
    //
    dtemp = _rate(ener);
  }

  if(dtemp < 0.)
    //
    return 0.;

  return dtemp;
}

/********************************************************************************************
 **************************************** WELL MODEL ****************************************
 ********************************************************************************************/

Model::Well::Well(IO::KeyBufferStream& from, const std::string& n) 
  : _extension(-1.)
{
  const char funame [] = "Model::Well::Well: ";

  KeyGroup WellModel;
  Key  relax_key("EnergyRelaxation");
  Key   spec_key("Species");
  Key escape_key("Escape");
  Key    ext_key("Extension");
  Key   freq_key("CollisionFrequency");

  std::string token, comment;

  while(from >> token) {// input cycle
    // end input 
    if(IO::end_key() == token) {
      std::getline(from, comment);
      break;
    }
    // collisional energy relaxation kernel
    else if(relax_key == token) {
      _kernel.push_back(new_kernel(from));
    }
    // collision frequency model
    else if(freq_key == token) {
      _collision.push_back(new_collision(from));
    }
    // species specification
    else if(spec_key == token) {
      if(_species) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      _species = new_species(from, n, DENSITY);
    }
    // escape specification
    //
    else if(escape_key == token) {
      if(_escape) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      _escape = new_escape(from);
    }
    // well extension
    //
    else if(ext_key == token) {
      if(_extension > 0.) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      if(!(from >> _extension)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(_extension <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // unknown key
    else if(IO::skip_comment(token, from)) {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }// input cycle

  if(!_default_kernel.size()) {
    std::cerr << funame << "default kernel(s) has not been defined yet\n";
    throw Error::Init();
  }

  if(!_kernel.size()) {
    _kernel = _default_kernel;
  }
  else if(_kernel.size() != _default_kernel.size()) {
    std::cerr << funame << "number of kernels mismatch with the default\n";
    throw Error::Init();
  }

  if(!_default_collision.size()) {
    std::cerr << funame << "default collision model has not been defined yet\n";
    throw Error::Init();
  }

  if(!_collision.size()) {
    _collision = _default_collision;
  }
  else if(_collision.size() != _default_collision.size()) {
    std::cerr << funame << "number of collision models mismatch with the default\n";
    throw Error::Init();
  }

  if(!_species) {
    std::cerr << funame << "species not inititalized\n";
    throw Error::Init();
  }
}

// radiational transition down probability 
double Model::Well::transition_probability (double ener, double temperature, int num) const
{
  const char funame [] = "Model::Well::transition_probability: ";

  static const double nfac = 2. / M_PI / Phys_const::light_speed;

  if(!_species) {
    std::cerr << funame << "species not initialized\n";
    throw Error::Init();
  }

  if(num < 0 || num >= oscillator_size()) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }

  if(temperature <= 0.) {
    std::cerr << funame << "negative temperature\n";
    throw Error::Range();
  }

  double res = _species->infrared_intensity(ener, num);
  if(res == 0.)
    return res;

  double dtemp = _species->oscillator_frequency(num);
  
  res *=  nfac * dtemp * dtemp;

  dtemp /= temperature;

  if(dtemp < 30.)
    res /= 1. - std::exp(-dtemp);

  return res;
}

/********************************************************************************************
 **************************************** THERMOCHEMISTRY ************************************
 ********************************************************************************************/

