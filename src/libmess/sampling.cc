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

#include "structure.hh"
#include "read.hh"
#include "configuration.hh"
#include "random.hh"
#include "key.hh"

#include <list>
#include <fstream>
#include <csignal>
#include <ctime>
#include <omp.h>

std::string state_data_file;
std::list<Configuration::State> vertex;
typedef std::list<Configuration::State>::const_iterator Vit;

void save_state_data ()
{
  if(!vertex.size())
    return;

  std::ofstream to(state_data_file.c_str());
  
  to << std::setprecision(17) << std::scientific;
  to << vertex.size() << "\n";
  for(Vit v = vertex.begin(); v != vertex.end(); ++v)
    to << *v;
}

extern "C" void signal_handler (int sig)
{
  //if(!omp_get_thread_num())
  save_state_data();

  if(sig == SIGUSR1 || sig == SIGUSR2)
    return;

  std::exit(1);
}

double atom_dist_min = 1.5;

std::vector<Atom> state2geom (double r, const Configuration::State& ang) 
{
  D3::Matrix rmatrix;
  quat2mat(ang.orientation(), rmatrix);

  D3::Vector rpos(ang.radius_vector());
  rpos *= r;

  std::vector<Atom> res;
  Atom curr;
  for(int a = 0; a < Structure::fragment(0).size(); ++a) {
    curr = Structure::fragment(0)[a];
    curr /= Phys_const::angstrom;
    res.push_back(curr);
  }

  for(int a = 0; a < Structure::fragment(1).size(); ++a) {
    curr = Structure::fragment(1)[a];
    D3::vprod(Structure::fragment(1)[a], rmatrix, curr);
    curr += rpos;
    for(int b = 0; b < Structure::fragment(0).size(); ++b)
      if(vdistance(curr, Structure::fragment(0)[b]) < atom_dist_min)
	throw Error::Range();
    curr /= Phys_const::angstrom;
    res.push_back(curr);
  }

  return res;
}

int main (int argc, char* argv [])
{
  const char funame [] = "sampling: ";

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  // input
  if (argc < 2) {
    std::cout << funame  << "usage: sampling input_file\n";
    throw Error::Input();
  }

  // base name
  std::string base_name = argv[1];
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    base_name.resize(base_name.size() - 4);

  // random seed initialization
  Random::init();

  // start time
  const std::time_t  start_time = std::time(0);
  const std::clock_t start_cpu  = std::clock();

  // input
  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    throw Error::Input();
  }

  // input parameters
  int    miss_count_max = -1;
  double angle_spacing  = -1.;

  std::vector<double> dist_grid;
  int                 dist_size  = 0;
  double              dist_min   = -1.;
  double              dist_max   = -1.;

  std::ofstream geom_out;

  KeyGroup MainGroup;
 
  Key struc_key("Structure"         );
  Key  dmin_key("DistanceMin[bohr]" );
  Key  dmax_key("DistanceMax[bohr]" );
  Key dsize_key("DistanceSize"      );
  Key dgrid_key("DistanceGrid[bohr]");
  Key  miss_key("MissCountMax"      );
  Key   asp_key("AngularSpacing"    );
  Key   log_key("LogOutput"         );
  Key  stat_key("StateOutput"       );
  Key  geom_key("GeometryOutput"    );
  Key   adm_key("AtomDistanceMin[bohr]");

  std::string token, comment;
  while(from >> token) {
    // log output
    if(log_key == token) {
      if(IO::log.is_open()) {
        std::cerr << funame << token << ": allready initialzed\n";
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
    // angular state output
    else if(stat_key == token) {
      if(state_data_file.size()) {
        std::cerr << funame << token << ": allready initialized\n";
        throw Error::Init();
      }      

      if(!(from >> state_data_file)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);
    }
    // molecular geometry output
    else if(geom_key == token) {
      if(geom_out.is_open()) {
        std::cerr << funame << token << ": allready initialized\n";
        throw Error::Init();
      }      
      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      geom_out.open(stemp.c_str());
      if(!geom_out) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }
    }
    // molecular structure
    else if(struc_key == token) {
      // default log output	
      if(!IO::log.is_open()) {
	stemp = base_name + ".log";
	IO::log.open(stemp.c_str());
	if(!IO::log) {
	  std::cerr << funame << token << ": cannot open " << stemp << " file\n";
	  throw Error::Input();
	}
      }

      std::getline(from, comment);
      Structure::init(from);
    }
    // minimal distance
    else if(dmin_key == token) {
      if(dist_grid.size() || dist_min > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> dist_min)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(dist_min <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // maximal distance
    else if(dmax_key == token) {
      if(dist_grid.size() || dist_max > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> dist_max)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(dist_max <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // distance size
    else if(dsize_key == token) {
      if(dist_grid.size() || dist_size > 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> dist_size)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(dist_size < 2) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // distance grid
    else if(dgrid_key == token) {
      if(dist_grid.size() || dist_size > 0 || dist_max > 0. || dist_min > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      IO::LineInput lin(from);
      while(lin >> dtemp) {

	if(dtemp <= 0.) {
	  std::cerr << funame << token << ": out of range\n";
	  throw Error::Range();
	}

	if(dist_grid.size() && dtemp >= dist_grid.back()) {
	  std::cerr << funame << token << ": should be in descending order\n";
	  throw Error::Range();
	}

	dist_grid.push_back(dtemp);
      }
      
      if(!dist_grid.size()) {
	std::cerr << funame << token << ": no grid\n";
	throw Error::Range();
      }
    }
    // angular spacing
    else if(asp_key == token) {
      if(angle_spacing > 0.) {
	std::cerr << funame << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> angle_spacing)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(angle_spacing <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // miss count max
    else if(miss_key == token) {
      if(miss_count_max > 0) {
	std::cerr << funame << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> miss_count_max)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(miss_count_max <= 0) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // interatomic minimal distance
    else if(adm_key == token) {
      if(!(from >> atom_dist_min)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }

      if(atom_dist_min <= 0) {
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
  from.close();
  from.clear();

  if(miss_count_max < 0) {
    std::cerr << funame << miss_key << ": not initialized\n";
    throw Error::Init();
  }

  if(angle_spacing < 0.) {
    std::cerr << funame << asp_key << ": not initialized\n";
    throw Error::Init();
  }

  if(!dist_grid.size() && (dist_size < 2 || dist_min <= 0. || dist_max <= 0. || dist_max <= dist_min)) {
    std::cerr << funame << "distance grid not initialized properly\n";
    throw Error::Range();
  }

  if(!state_data_file.size())
    state_data_file = base_name + ".stat";

  if(!geom_out.is_open())
    geom_out.open((base_name + ".geom").c_str());

  // distance grid
  if(!dist_grid.size()) {
    dtemp = std::pow(dist_max / dist_min , 1./double(dist_size - 1));

    double dist = dist_max;
    for(int d = 0; d < dist_size; ++d, dist /= dtemp)
      dist_grid.push_back(dist);
  }

  // recover vertex data
  from.open(state_data_file.c_str());
  if(from >> itemp) {
    Configuration::State vtemp;
    for(int v = 0; v < itemp; ++v) {
      if(!(from >> vtemp)) {
	std::cerr << funame << "vertex data file is corrupted\n";
	throw Error::Input();
      }
      vertex.push_back(vtemp);
    }
  }
  from.close();
  from.clear();

  if(Structure::fragment(0).type() != Molecule::NONLINEAR) {
    std::cerr << funame << "first fragment should be nonlinear\n";
    throw Error::Logic();
  }

  // layout and symmetry group
  std::vector<int> layout(1, 3);

  // symmetry group
  std::vector<Symmetry::SpaceGroup> sg;

  SharedPointer<Configuration::GroupBase> symm_group;

  switch(Structure::fragment(1).type()) {
  case Molecule::NONLINEAR:
    layout.push_back(4); // second fragment orientation (quaternion)

    for(int frag = 0; frag < 2; ++frag)
      sg.push_back(*Structure::fragment(frag).symmetry_group);
    symm_group.init(new Configuration::DoubleSpaceGroup(sg));

    break;
  case Molecule::LINEAR:
    layout.push_back(3); // second fragment orientation (vector)
    
    switch(Structure::fragment(1).symmetry_group->size()) {
    case 1:
      // generic non-symmetric linear fragment
      symm_group.init(new Configuration::SpaceGroup(*Structure::fragment(0).symmetry_group, 0));
      break;

    case 2:
      // linear fragment with inversion symmetry
      symm_group.init(new Configuration::SpaceGroup(*Structure::fragment(0).symmetry_group, 1));
      break;
      
    default:
      std::cerr << funame << "linear fragment: wrong case\n";
      throw Error::Logic();
    }
    break;

  case Molecule::MONOATOMIC:
    symm_group.init(new Configuration::SpaceGroup(*Structure::fragment(0).symmetry_group));
    break;
    
  default:
    std::cerr << funame << "second fragment: wrong case\n";
    throw Error::Logic();
  }

  Configuration::State::layout.set(layout);


  IO::Marker funame_marker(funame);

  // signal handing
  sigset_t block_sig;
  sigset_t old_sig;
  sigfillset(&block_sig); // block all signals

  struct sigaction sigact;
  sigact.sa_handler = signal_handler;
  sigact.sa_mask = block_sig;
  sigact.sa_flags = 0;

  sigaction(SIGUSR1, &sigact, 0);
  sigaction(SIGUSR2, &sigact, 0);
  sigaction(SIGINT,  &sigact, 0);
  sigaction(SIGTERM, &sigact, 0);
  sigaction(SIGHUP,  &sigact, 0);
  sigaction(SIGABRT, &sigact, 0);
  sigaction(SIGQUIT, &sigact, 0);
  sigaction(SIGTRAP, &sigact, 0);
  sigaction(SIGALRM, &sigact, 0);
  sigaction(SIGFPE,  &sigact, 0);
  sigaction(SIGILL,  &sigact, 0);
  sigaction(SIGBUS,  &sigact, 0);
  sigaction(SIGSEGV, &sigact, 0);
  sigaction(SIGPIPE, &sigact, 0);

  // configuration states
  Configuration::State  guess;
  std::vector<Configuration::State> guess_orbit(symm_group->size(), guess);

  int miss_count = 0;
  while(miss_count < miss_count_max) {

    // randomly initialize
    Random::orient(guess.radius_vector(), 3);
    Random::orient(guess.orientation(),   4);

    if(!vertex.size()) {
      vertex.push_back(guess);
      continue;
    }
    
    // generate guess orbit
    guess_orbit[0] = guess;

    
    //#pragma omp parallel for default(shared) schedule(static)      
    for(int g = 1; g < symm_group->size(); ++g)
      symm_group->apply(g, guess, guess_orbit[g]);

    // check the distances to other vertices
    bool miss = false;
    for(Vit v = vertex.begin();  v != vertex.end(); ++v) {

      //#pragma omp parallel for default(shared) schedule(static)      
      for(int g = 0; g < symm_group->size(); ++g) 
	if(vdistance(guess_orbit[g], *v) < angle_spacing) {
	  
	  //#pragma omp critical
	  miss = true;
	  break;
	}

      if(miss)
	break;
    }
    
    if(miss) {
      miss_count++;
      continue;
    }

    IO::log << IO::log_offset << "number of vertices = " << std::setw(6) <<  vertex.size() 
	    << ";  miss count = " << std::setw(6) << miss_count 
	    << ";  cpu time [sec] = " << std::setw(6) << std::ceil(double(std::clock() - start_cpu) / CLOCKS_PER_SEC)
	    << ";  elapsed time [sec] = " << std::setw(6) << std::time(0) - start_time
	    << std::endl;

    vertex.push_back(guess);
    miss_count = 0;
  }

  // state output
  sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
  save_state_data();
  sigprocmask(SIG_SETMASK, &old_sig, 0);

  // geometry output
  geom_out << vertex.size() << "   " << Structure::size() << "\n";
  itemp = 0;
  for(Vit v = vertex.begin(); v != vertex.end(); ++v, ++itemp) {
    std::vector<std::vector<Atom> > geom;
    try {
      for(int d = 0; d < dist_grid.size(); ++d)
	geom.push_back(state2geom(dist_grid[d], *v));
    } catch(Error::General) {}
    

    geom_out << itemp << "   " << geom.size() << "\n";
    for(int g = 0; g < geom.size(); ++g) {
      for(std::vector<Atom>::const_iterator at = geom[g].begin(); at != geom[g].end(); ++at)
	geom_out << *at << "\n";

      geom_out << "\n";
    }
  }

  return 0;
}
