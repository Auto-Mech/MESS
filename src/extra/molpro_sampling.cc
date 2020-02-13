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

#include "molpro.hh"
#include "key.hh"
#include "io.hh"
#include "system.hh"

#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <csignal>
#include <cstring>

#include <mpi.h>

#include <fstream>
#include <sstream>

enum {RUN_TAG, STOP_TAG};

std::string scratch_dir;
std::string ener_data_file;
std::vector<Array<double> > ener_data;

int pes_size = 1;

void save_ener_data ()
{
  if(!ener_data.size())
    return;

  std::ofstream ener_out(ener_data_file.c_str());

  ener_out << ener_data.size() << "   " << pes_size << "\n"
	   << std::left << std::scientific;
  for(int v = 0; v < ener_data.size(); ++v) {
    ener_out << std::setw(7) << v << " "
	     << std::setw(3) << ener_data[v].size() 
	     << "\n";

    for(const double* et = ener_data[v].begin(); et != ener_data[v].end(); ++et)
      ener_out << std::setw(15) << *et;
    ener_out << "\n";
  }

  ener_out.close();
}

extern "C" void signal_handler (int sig)
{
  if(sig == SIGUSR1) {
    if(!MPI::COMM_WORLD.Get_rank()) {
      std::cerr << "signal_handler: saving energy data ... ";
      save_ener_data();
      std::cerr << "done" << std::endl;
    }

    return;
  }

  if(!MPI::COMM_WORLD.Get_rank())
    save_ener_data();
  else
    System::call_exe("/bin/rm", "-Rf", scratch_dir.c_str(), (char*) 0);
  
  MPI::Finalize();
  std::exit(1);
}

class ErrLog {
public:
  template <typename C> ErrLog& operator<< (const C&);
};

template <typename C>
ErrLog& ErrLog::operator<< (const C& c) 
{
  if(!MPI::COMM_WORLD.Get_rank())
    std::cerr << c;
  return *this;
}

int make_dir (const std::string& dir, const std::string& funame) {

  std::string token = funame;
  token += dir + ": mkdir: ";

  struct stat wstat;

  if(mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH ))
    switch(errno) {
    case EEXIST:
      if(stat(dir.c_str(), &wstat)) {
	std::cerr << token << strerror(errno) << "\n";
	return 1;
      }

      if(!S_ISDIR(wstat.st_mode)) {
	std::cerr << token << "not a directory\n";
	return 1;
      }
      
      if(!(wstat.st_mode & S_IRUSR) || !(wstat.st_mode & S_IWUSR) || !(wstat.st_mode & S_IXUSR)) {
	std::cerr << token << "wrong permissions\n";
	return 1;
      }

      return 0;

    case ENOENT:
      std::cerr << token << "the parent directory does not exist\n";
      return 1;

    default:
      std::cerr << token << strerror(errno) << "\n";
      return 1;
    }

  return 0;
}

int main (int argc, char* argv[])
{
  MPI::Init(argc, argv);

  const int mpi_size = MPI::COMM_WORLD.Get_size();
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();

  // funame
  std::string  funame = "molpro_sampling: ";
  if(mpi_rank) {
    std::ostringstream node_name;
    node_name << "node" << std::setw(5) <<  mpi_rank << ": ";
    funame += node_name.str();
  }
  else
    funame += "master: ";

  // error log
  ErrLog err_log;

  // usage
  if(argc != 2) {
    err_log << "usage: molpro_sampling input_file\n";
    MPI::Finalize();
    return 1;
  }

  double      dtemp;
  int         itemp;
  bool        btemp;
  std::string stemp;

  std::ifstream from(argv[1]);
  if(!from) {
    err_log << funame << "cannot open " << argv[1] << " file\n";
    MPI::Finalize();
    return 1;
  }

  // input parameters
  std::string work_dir;

  std::ifstream geom_in;

  double ener_threshold = 0.5 * Phys_const::kcal;

  KeyGroup MolproSamplingGroup;

  Key molpro_key("Molpro"                   );
  Key   geom_key("SamplingFile"             );
  Key   ener_key("EnergyOutput"             );
  Key   work_key("MolproWorkDir"            );
  Key    scr_key("MolproScratchDir"         );
  Key  thres_key("EnergyThreshold[kcal/mol]");
  Key    pes_key("PesSize"                  );

  bool input_error = false;
  std::string token, comment;
  while(from >> token) {
    // molpro initialization
    if(molpro_key == token) {
      token += ": ";

      if(Molpro::isinit()) {
	err_log << funame << token << "already initialized\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);

      try {
	Molpro::init(from);
      }
      catch(Error::General) {
	input_error = true;
	break;
      }
    }
    // PES size
    else if(pes_key == token) {
      token += ": ";
      if(!(from >> pes_size)) {
	err_log << funame << token << "corrupted\n";
	input_error = true;
	break;
      }
      if(pes_size < 1) {
	err_log << funame << token << "out of range\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);
    }
    // working directory
    else if(work_key == token) {
      token += ": ";

      if(work_dir.size()) {
	err_log << funame << token << "already initialized\n";
	input_error = true;
	break;
      }

      if(!(from >> work_dir)) {
	err_log << funame << token << "corrupted\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);
    }
    // scratch directory
    else if(scr_key == token) {
      token += ": ";

      if(scratch_dir.size()) {
	err_log << funame << token << "already initialized\n";
	input_error = true;
	break;
      }

      if(!(from >> scratch_dir)) {
	err_log << funame << token << "corrupted\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);
    }
    // geometry configurations
    else if(geom_key == token) {
      token += ": ";

      if(geom_in.is_open()) {
	err_log << funame << token << "already initialized\n";
	input_error = true;
	break;
      }

      if(!(from >> stemp)) {
	err_log << funame << token << "corrupted\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);

      geom_in.open(stemp.c_str());
      if(!geom_in) {
	err_log << funame << token << "cannot open " << stemp << " file\n";
	input_error = true;
	break;
      }
    }
    // energy data input/output
    else if(ener_key == token) {
      token += ": ";

      if(ener_data_file.size()) {
	err_log << funame << token << "already initialized\n";
	input_error = true;
	break;
      }

      if(!(from >> ener_data_file)) {
	err_log << funame << token << "corrupted\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);

      // master
      if(!mpi_rank) {
	std::ifstream ener_in(ener_data_file.c_str());
	if(ener_in.is_open()) {
	  ener_in >> itemp;
	  ener_data.resize(itemp);

	  for(int v = 0; v < ener_data.size(); ++v) {
	    ener_in >> itemp >> itemp;
	    ener_data[v].resize(itemp);
	    for(int d = 0; d < ener_data[v].size(); ++d)
	      ener_in >> ener_data[v][d];
	  }

	  if(!ener_in) {
	    err_log << funame << token << "corrupted\n";
	    input_error = true;
	    break;
	  }

	  ener_in.close();
	}
      }// master
    }
    // energy threshold
    else if(thres_key == token) {
      token += ": ";

      if(!(from >> ener_threshold)) {
	err_log << funame << token << "corrupted\n";
	input_error = true;
	break;
      }
      std::getline(from, comment);

      if(ener_threshold <= 0.) {
	err_log << funame << token << "out of range\n";
	input_error = true;
      }

      ener_threshold *= Phys_const::kcal;
    }
    // unknown keyword
    else if(IO::skip_comment(token, from)) {
      err_log << funame << "unknown keyword " << token << "\n";
      if(!mpi_rank)
	Key::show_all(std::cerr);
      err_log << "\n";
      input_error = true;
      break;
    }
  }

  if(!work_dir.size()) {
    err_log << funame << "no working directory\n";
    input_error = true;
  }

  if(!scratch_dir.size()) {
    err_log << funame << "no scratch directory\n";
    input_error = true;
  }

  if(!ener_data_file.size()) {
    err_log << funame << "no energy data file\n";
    input_error = true;
  }

  if(!geom_in.is_open()) {
    err_log << funame << "no sampling data\n";
    input_error = true;
  }

  // working and scratch directories
  if(!input_error && mpi_rank) {
    // node name
    std::ostringstream node_name;
    node_name << "/node" << mpi_rank;
  
    // set working directory
    input_error |= make_dir(work_dir, funame);

    stemp = work_dir + node_name.str();
    input_error |= make_dir(stemp, funame);

    if(!input_error && chdir(stemp.c_str())) {
      std::cerr << funame << stemp << ": chdir: " << strerror(errno) << "\n";
      input_error = true;
    }

    // set scratch directory
    input_error |= make_dir(scratch_dir, funame);
  
    scratch_dir += node_name.str();
    input_error |= make_dir(scratch_dir, funame);

    if(!input_error)
      Molpro::set_scratch_dir(scratch_dir);
  }

  int vertex_size, atom_size;
  if(!(geom_in >> vertex_size >> atom_size)) {
    err_log << funame << "wrong sampling file format\n";

    input_error = true;
  }

  if(!mpi_rank && ener_data.size() && ener_data.size() > vertex_size) {
    err_log << funame << "sampling data and energies data mismatch\n";
    input_error = true;
  }

  // synchronize the error 
  for(int n = 0; n < mpi_size; ++n) {
    if(mpi_rank == n)
      btemp = input_error;
    MPI::COMM_WORLD.Bcast(&btemp, 1, MPI::BOOL, n);
    input_error |= btemp;
  }

  if(input_error) {
    MPI::Finalize();
    return 0;
  }

  // signal handing
  sigset_t block_sig;
  sigset_t old_sig;
  sigfillset(&block_sig); // block all signals

  struct sigaction sigact;
  sigact.sa_handler = signal_handler;
  sigact.sa_mask = block_sig;
  sigact.sa_flags = 0;

  sigaction(SIGINT,  &sigact, 0);
  sigaction(SIGTERM, &sigact, 0);
  sigaction(SIGUSR1, &sigact, 0);
  sigaction(SIGUSR2, &sigact, 0);
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

  /**********************************************************************************
   ********************************** MASTER ****************************************
   **********************************************************************************/
  if(!mpi_rank) {

    sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
    ener_data.resize(vertex_size);
    sigprocmask(SIG_SETMASK, &old_sig, 0);

    int vertex = 0;
    int node;

    // start work
    for(node = 1; node < mpi_size; ++node) {
      // find new vertex
      while(vertex < vertex_size && ener_data[vertex].size())
	++vertex;

      if(vertex == vertex_size)
	break;

      // send new vertex
      MPI::COMM_WORLD.Send(&vertex, 1, MPI::INT, node, RUN_TAG);
      ++vertex;
    }

    // stop unused nodes
    for(int i = node; i < mpi_size; ++i)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, i, STOP_TAG);

    --node;

    int vertex_data [2];
    MPI::Status stat;
    while(node) {
      // get vertex info
      MPI::COMM_WORLD.Recv(vertex_data, 2, MPI::INT, MPI::ANY_SOURCE, RUN_TAG, stat);

      // get new energies
      Array<double> vtemp(vertex_data[1]);
      MPI::COMM_WORLD.Recv(vtemp, vertex_data[1], MPI::DOUBLE, stat.Get_source(), RUN_TAG);

      // update energy data
      sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
      ener_data[vertex_data[0]] = vtemp;
      sigprocmask(SIG_SETMASK, &old_sig, 0);

      // find new vertex
      while(vertex < vertex_size && ener_data[vertex].size())
	++vertex;

     // send new vertex
      if(vertex < vertex_size) {
	MPI::COMM_WORLD.Send(&vertex, 1, MPI::INT, stat.Get_source(), RUN_TAG);
	++vertex;
      }
      else {
	MPI::COMM_WORLD.Send(0,       0, MPI::INT, stat.Get_source(), STOP_TAG);
 	--node;
      }
    }

    // save energies
    sigprocmask(SIG_SETMASK, &block_sig, &old_sig);
    save_ener_data();
    sigprocmask(SIG_SETMASK, &old_sig, 0);

    std::cerr << funame << "done\n";
  }
  /**********************************************************************************
   *********************************** SLAVE ****************************************
   **********************************************************************************/
  else {
    // open log file
    IO::log.open("molpro.log");

    int vertex;
    MPI::Status stat;

    std::vector<Atom> molec(atom_size);

    int vertex_data [2];
    int dist_size;
    
    while(1) {// work loop
      MPI::COMM_WORLD.Recv(vertex_data, 1, MPI::INT, 0, MPI::ANY_TAG, stat);

      // stop calculation
      if(stat.Get_tag() == STOP_TAG) {
	std::cerr << funame << "cleaning up and exitting\n";
	System::call_exe("/bin/rm", "-Rf", scratch_dir.c_str(), (char*) 0);
	break;
      }
      // run molpro calculation
      else {
	geom_in >> vertex >> dist_size;

	// skip geometries
	while(vertex_data[0] > vertex) {
	  for(int dist = 0; dist < dist_size; ++dist)
	      for(std::vector<Atom>::iterator at = molec.begin(); at != molec.end(); ++at)
		geom_in >> *at;
	  
	  if(!(geom_in >> vertex >> dist_size)) {
	    std::cerr << funame << "geometry input stream is corrupted\n";
	    throw Error::Input();
	  }
	} 

	// remove old wave function file
	Molpro::remove_wfu();
	
	// calculate energies;
	std::vector<double> ener_vec;

	bool molpro_fail = false;

	for(int dist = 0; dist < dist_size; ++dist) {
	  // read geometries
	  for(std::vector<Atom>::iterator at = molec.begin(); at != molec.end(); ++at)
	    geom_in >> *at;
	    
	  // molpro energy
	  if(!molpro_fail) {
	    try {
	      Array<double> ener(pes_size);
	      Molpro::pot(molec, ener);

	      dtemp = ener[0];
	      // energy test
	      if(!dist && (dtemp < -ener_threshold || dtemp > ener_threshold)) {
		molpro_fail = true;
		std::cerr << funame 
			  << "vertex =" << std::setw(5) << vertex_data[0] 
			  << ",  distance =" << std::setw(3) << dist 
			  << ": starting energy," << dtemp / Phys_const::kcal 
			  << " kcal/mol, exceeds the threshold" << std::endl;
	      }
	      else
		for(int e = 0; e < pes_size; ++e)
		  ener_vec.push_back(ener[e]);
	    }
	    catch(Error::Molpro& mess) {
	      molpro_fail = true;
	      std::cerr << funame 
			<< "vertex =" << std::setw(5) << vertex_data[0] 
			<< ",  distance =" << std::setw(3) << dist 
			<< ": " << mess << std::endl;
	    }
	  }
	}

	vertex_data[1] = ener_vec.size();
	Array<double> ener((int)ener_vec.size());
	for(int e = 0; e < ener.size(); ++e)
	  ener[e] = ener_vec[e];

	// send vertex energies
	MPI::COMM_WORLD.Send(vertex_data,  2,           MPI::INT,    0, RUN_TAG);
	MPI::COMM_WORLD.Send(ener,         ener.size(), MPI::DOUBLE, 0, RUN_TAG);
      }// run molpro calculation
    }// work loop
  }// slave

  MPI::Finalize();
  return 0;
}
