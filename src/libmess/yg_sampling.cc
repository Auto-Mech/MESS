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
#include "structure.hh"
#include "dynamic.hh"

#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <list>

#include <mpi.h>

#include <fstream>
#include <sstream>

enum {RUN_TAG, STOP_TAG, FAIL_TAG, HIGH_ENER_TAG, LOW_FREQ_TAG, HIGH_FREQ_TAG, CHECK_TAG, DIR_TAG};

std::string scratch_dir;

int mpi_rank = 0, stop_request = 0;

extern "C" void signal_handler (int sig) { stop_request = 1; }

class ErrLog {
  //
public:
  //
  template <typename C> ErrLog& operator<< (const C&);
};

template <typename C>
ErrLog& ErrLog::operator<< (const C& c) 
{
  if(!mpi_rank)
    //
    std::cerr << c;
  
  return *this;
}

int main (int argc, char* argv[])
{
  double      dtemp;
  int         itemp;
  int         btemp;
  std::string stemp;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &itemp);
   
  const int mpi_size = itemp;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if(mpi_rank) {
    //
    Molpro::mute = 1;

    Structure::mute = 1;
  }
  else {
    //
    Molpro::mute = 0;
    Structure::mute = 0;
  }

  std::cout.precision(3);
  
  std::ostringstream oss;

  // funame
  //
  std::string  funame;

  if(mpi_rank) {
    //
    oss.str("");

    oss << mpi_rank << "-th node: ";

    funame += oss.str();
  }
  else
    //
    funame += "master: ";

  // error log
  //
  ErrLog err_log;

  // usage
  //
  if(argc != 2) {
    //
    err_log << "usage: yg_sampling input_file\n";

    MPI_Finalize();

    return 1;
  }

  std::ifstream from(argv[1]);

  if(!from) {
    //
    err_log << funame << "cannot open " << argv[1] << " file\n";

    MPI_Finalize();

    return 1;
  }

  char hostname [100];
  gethostname(hostname, 99);

  Array<char> host_buff;

  if(!mpi_rank)
    //
    host_buff.resize(100 * mpi_size);

  MPI_Gather(hostname, 100, MPI_CHAR, host_buff, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

  // input parameters
  //
  int fail_count_max = 3;
  
  std::string work_dir;

  double orb_polar_min  = 0., orb_polar_max  =      M_PI;
  double orb_azim_min   = 0., orb_azim_max   = 2. * M_PI;

  double frag_polar_min = 0., frag_polar_max =      M_PI;
  double frag_azim_min  = 0., frag_azim_max  = 2. * M_PI;
  double frag_psi_min   = 0., frag_psi_max   = 2. * M_PI;
  
  int orb_axis = -1, frag_axis = -1;

  double orb_step = -1., frag_step = -1.;

  double dist_min, dist_max, dist_step;

  std::set<double> dist_set;
  
  double ener_min = 1., ener_max = -1.;

  double low_freq_min = -1., high_freq_min = -1., high_freq_max = -1.;

  std::string out_file_name;

  bool is_eref = false;
  double ener_ref;

  double geom_relax_dist = -1.;
  double geom_relax_ener = -1.;

  double ener_incr_max = -1., ener_decr_max = -1.;
  
  KeyGroup YGsamplingGroup;

  Key  struc_key("Structure"                );
  Key molpro_key("Molpro"                   );
  Key   work_key("MolproWorkDir"            );
  Key    scr_key("MolproScratchDir"         );
  Key    out_key("OutputFile"               );
  Key  ostep_key("OrbStep"                  );
  Key    opr_key("OrbPolarRange"            );
  Key    oar_key("OrbAzimuthRange"          );
  Key     ox_key("OrbAxis"                  );
  Key  fstep_key("FragStep"                 );
  Key    fpr_key("FragPolarRange"           );
  Key    far_key("FragAzimuthRange"         );
  Key    fsr_key("FragPsiRange"             );
  Key     fx_key("FragAxis"                 );
  Key drange_key("DistRange[bohr]"          );
  Key erange_key("EnerRange[kcal/mol]"      );
  Key  lfmin_key("LowFreqMin[1/cm]"         );
  Key  hfmin_key("HighFreqMin[1/cm]"        );
  Key  hfmax_key("HighFreqMax[1/cm]"        );
  Key incmax_key("EnerIncrMax[kcal/mol]"    );
  Key decmax_key("EnerDecrMax[kcal/mol]"    );
  Key   eref_key("EnerReference[au]"        );
  Key relaxd_key("GeomRelaxDist[bohr]"      );
  Key relaxe_key("GeomRelaxEner[kcal/mol]"  );
  Key   fail_key("FailCountMax"             );
  
  bool input_error = false;
  std::string token, comment;
 
  while(from >> token) {
    //
    // molecular structure initialization
    //
    if(struc_key == token) {
      //
      token += ": ";

      std::getline(from, comment);

      if(Structure::isinit()) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      
      try {
	//
	Structure::init(from);
      }
      catch(Error::General) {
	//
	input_error = true;
	break;
      }
    }
    // molpro initialization
    //
    else if(molpro_key == token) {
      //
      token += ": ";

      std::getline(from, comment);

      if(Molpro::isinit()) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      
      try {
	Molpro::init(from);
      }
      catch(Error::General) {
	//
	input_error = true;
	break;
      }
    }
    // reference energy
    //
    else if(eref_key == token) {
      //
      token += ": ";

      if(is_eref) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      is_eref = true;
      
      IO::LineInput lin(from);
      
      if(!(lin >> ener_ref)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
    }
    // maximal energy increment
    //
    else if(incmax_key == token) {
      //
      token += ": ";

      if(ener_incr_max > 0.) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      
      IO::LineInput lin(from);
      
      if(!(lin >> ener_incr_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
      if(ener_incr_max <= 0.) {
	//
	err_log << funame << token << "should be positive\n";

	input_error = true;
	break;
      }
      ener_incr_max *= Phys_const::kcal;
    }
    // maximal energy decrement
    //
    else if(decmax_key == token) {
      //
      token += ": ";

      if(ener_decr_max > 0.) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      
      IO::LineInput lin(from);
      
      if(!(lin >> ener_decr_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
      if(ener_decr_max <= 0.) {
	//
	err_log << funame << token << "should be positive\n";

	input_error = true;
	break;
      }
      ener_decr_max *= Phys_const::kcal;
    }
    // geometry relaxation distance range
    //
    else if(relaxd_key == token) {
      //
      token += ": ";

      if(geom_relax_dist > 0.) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      
      IO::LineInput lin(from);
      
      if(!(lin >> geom_relax_dist)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
      if(geom_relax_dist <= 0.) {
	//
	err_log << funame << token << "should be positive\n";

	input_error = true;
	break;
      }
    }
    // geometry relaxation energy range
    //
    else if(relaxe_key == token) {
      //
      token += ": ";

      if(geom_relax_ener > 0.) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }
      
      IO::LineInput lin(from);
      
      if(!(lin >> geom_relax_ener)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
      if(geom_relax_ener <= 0.) {
	//
	err_log << funame << token << "should be positive\n";

	input_error = true;
	break;
      }
      
      geom_relax_ener *= Phys_const::kcal;
    }
    // working directory
    //
    else if(work_key == token) {
      //
      token += ": ";

      if(work_dir.size()) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> work_dir)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
    }
    // scratch directory
    //
    else if(scr_key == token) {
      //
      token += ": ";

      if(scratch_dir.size()) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> scratch_dir)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
    }
    // output file
    //
    else if(out_key == token) {
      //
      token += ": ";

      if(out_file_name.size()) {
	//
	err_log << funame << token << "already initialized\n";
	
	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> out_file_name)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
    }
    // orbital (first fragment) polar angle range
    //
    else if(opr_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> orb_polar_min >> orb_polar_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(orb_polar_min >= orb_polar_max) {
	//
	err_log << funame << token << "first should go minimal value and then the maximal one\n";

	input_error = true;
	break;
      }

      if(orb_polar_min < 0. || orb_polar_max > 180.) {
	//
	err_log << funame << token << "out of range: " << orb_polar_min << ", " << orb_polar_max << ": should be [0-180]\n";

	input_error = true;
	break;
      }

      orb_polar_min *= M_PI / 180.;

      orb_polar_max *= M_PI / 180.;
    }
    // orbital (first fragment) azimuth angle range
    //
    else if(oar_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> orb_azim_min >> orb_azim_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(orb_azim_min >= orb_azim_max) {
	//
	err_log << funame << token << "first should go minimal value and then the maximal one\n";

	input_error = true;
	break;
      }

      if(orb_azim_min < 0. || orb_azim_max > 360.) {
	//
	err_log << funame << token << "out of range: " << orb_azim_min << ", " << orb_azim_max << ": should be [0-360]\n";

	input_error = true;
	break;
      }

      orb_azim_min *= M_PI / 180.;

      orb_azim_max *= M_PI / 180.;
    }
    // orbital (first fragment) (z-)axis
    //
    else if(ox_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(stemp == "x" || stemp == "X" || stemp == "0") {
	//
	orb_axis = 0;
      }
      else if(stemp == "y" || stemp == "Y" || stemp == "1") {
	//
	orb_axis = 1;
      }
      else if(stemp == "z" || stemp == "Z" || stemp == "2") {
	//
	orb_axis = 2;
      }
      else {
	//
	err_log << funame << token << "unknown symbol: " << stemp << ": should be x, y, or z\n";

	input_error = true;
	break;
      }
    }
    // second fragment polar angle range
    //
    else if (fpr_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> frag_polar_min >> frag_polar_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(frag_polar_min >= frag_polar_max) {
	//
	err_log << funame << token << "first should go minimal value and then the maximal one\n";

	input_error = true;
	break;
      }

      if(frag_polar_min < 0. || frag_polar_max > 180.) {
	//
	err_log << funame << token << "out of range: " << frag_polar_min << ", " << frag_polar_max << ": should be [0-180]\n";

	input_error = true;
	break;
      }

      frag_polar_min *= M_PI / 180.;

      frag_polar_max *= M_PI / 180.;
    }
    // second fragment azimuth angle range
    //
    else if(far_key == token) {
      //
      token += ": ";
      
      IO::LineInput lin(from);
      
      if(!(lin >> frag_azim_min >> frag_azim_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(frag_azim_min >= frag_azim_max) {
	//
	err_log << funame << token << "first should go minimal value and then the maximal one\n";

	input_error = true;
	break;
      }

      if(frag_azim_min < 0. || frag_azim_max > 360.) {
	//
	err_log << funame << token << "out of range: " << frag_azim_min << ", " << frag_azim_max << ": should be [0-360]\n";

	input_error = true;
	break;
      }

      frag_azim_min *= M_PI / 180.;

      frag_azim_max *= M_PI / 180.;
    }
    // second fragment psi angle range
    //
    else if(fsr_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> frag_psi_min >> frag_psi_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(frag_psi_min >= frag_psi_max) {
	//
	err_log << funame << token << "first should go minimal value and then the maximal one\n";

	input_error = true;
	break;
      }

      if(frag_psi_min < 0. || frag_psi_max > 360.) {
	//
	err_log << funame << token << "out of range: " << frag_psi_min << ", " << frag_psi_max << ": should be [0-360]\n";

	input_error = true;
	break;
      }

      frag_psi_min *= M_PI / 180.;

      frag_psi_max *= M_PI / 180.;
    }
    // second fragment (z-)axis
    //
    else if(fx_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> stemp)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(stemp == "x" || stemp == "X" || stemp == "0") {
	//
	frag_axis = 0;
      }
      else if(stemp == "y" || stemp == "Y" || stemp == "1") {
	//
	frag_axis = 1;
      }
      else if(stemp == "z" || stemp == "Z" || stemp == "2") {
	//
	frag_axis = 2;
      }
      else {
	//
	err_log << funame << token << "unknown symbol: " << stemp << ": should be x, y, or z\n";

	input_error = true;
	break;
      }
    }
    // orbital (1st fragment) angle step
    //
    else if(ostep_key == token) {
      //
      token += ": ";

      if(orb_step > 0.) {
	//
	err_log << funame << token << "already initialized\n";

	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> orb_step)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(orb_step <= 0.) {
	//
	err_log << funame << token << "out of range: " << orb_step << ": should be positive\n";
	
	input_error = true;
	break;
      }

      orb_step *= M_PI / 180.;
    }
    // 2nd fragment angle step
    //
    else if(fstep_key == token) {
      //
      token += ": ";

      if(frag_step > 0.) {
	//
	err_log << funame << token << "already initialized\n";

	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> frag_step)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(frag_step <= 0.) {
	//
	err_log << funame << token << "out of range: " << frag_step << ": should be positive\n";
	
	input_error = true;
	break;
      }

      frag_step *= M_PI / 180.;
    }
    // distance range (may be several, non-intersecting ones)
    //
    else if(drange_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> dist_min >> dist_max >> dist_step)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(dist_min <= 0. || dist_max <= dist_min || dist_step <= 1.) {
	//
	err_log << funame << token << "out of range: " << dist_min << ", " << dist_max << ", " << dist_step
	  //
		<< ": 0 < first(min) < second(max), third(multiplication step) > 1.\n";
	
	input_error = true;
	break;
      }

      if(dist_set.size() && dist_max >= *dist_set.begin() && dist_min <= *dist_set.rbegin()) {
	//
	err_log << funame << token << "distance intervals intersect\n";
	
	input_error = true;
	break;
      }
      
      dtemp = std::log(dist_max / dist_min);
	
      int n = std::ceil(dtemp / std::log(dist_step));

      dist_step = std::exp(dtemp / (double)n);

      dtemp = dist_min;
	
      for(int i = 0; i < n; ++i, dtemp *= dist_step)
	//
	dist_set.insert(dtemp);
    }
    // energy range
    //
    else if(erange_key == token) {
      //
      token += ": ";

      if(ener_min < 0. || ener_max > 0.) {
	//
	err_log << funame << token << "already initialized\n";

	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> ener_min >> ener_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(ener_min >= 0. || ener_max <= 0.) {
	//
	err_log << funame << token << "out of range: " << ener_min << ", " << ener_max << ": should be negative and positive correspondingly\n";
	
	input_error = true;
	break;
      }

      ener_min *= Phys_const::kcal;
      ener_max *= Phys_const::kcal;
    }
    // minimal conserved mode frequency
    //
    else if(lfmin_key == token) {
      //
      token += ": ";

      if(low_freq_min > 0.) {
	//
	err_log << funame << token << "already initialized\n";

	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> low_freq_min)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(low_freq_min <= 0.) {
	//
	err_log << funame << token << "out of range: " << low_freq_min << ": should be positive\n";
	
	input_error = true;
	break;
      }

      low_freq_min *= Phys_const::incm;
    }
    // minimal conserved mode high frequency
    //
    else if(hfmin_key == token) {
      //
      token += ": ";

      if(high_freq_min > 0.) {
	//
	err_log << funame << token << "already initialized\n";

	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> high_freq_min)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(high_freq_min <= 0.) {
	//
	err_log << funame << token << "out of range: " << high_freq_min << ": should be positive\n";
	
	input_error = true;
	break;
      }

      high_freq_min *= Phys_const::incm;
    }
    // maximal conserved mode frequency
    //
    else if(hfmax_key == token) {
      //
      token += ": ";

      if(high_freq_max > 0.) {
	//
	err_log << funame << token << "already initialized\n";

	input_error = true;
	break;
      }

      IO::LineInput lin(from);
      
      if(!(lin >> high_freq_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }

      if(high_freq_max <= 0.) {
	//
	err_log << funame << token << "out of range: " << high_freq_max << ": should be positive\n";
	
	input_error = true;
	break;
      }

      high_freq_max *= Phys_const::incm;
    }
    // maximal fail count
    //
    else if(fail_key == token) {
      //
      token += ": ";

      IO::LineInput lin(from);
      
      if(!(lin >> fail_count_max)) {
	//
	err_log << funame << token << "corrupted\n";
	
	input_error = true;
	break;
      }
      if(fail_count_max <= 0) {
	//
	err_log << funame << token << "should be positive\n";

	input_error = true;
	break;
      }
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      err_log << funame << "unknown keyword " << token << "\n";
      
      if(!mpi_rank)
	//
	Key::show_all(std::cerr);
      
      err_log << "\n";

      input_error = true;
      break;
    }
  }

  if(!input_error && !work_dir.size()) {
    //
    err_log << funame << "working directory not initialized\n";
    input_error = true;
  }

  if(!input_error && !scratch_dir.size()) {
    //
    err_log << funame << "scratch directory not initialized\n";
    input_error = true;
  }

  if(!input_error && orb_axis < 0) {
    //
    err_log << funame << "orbital (first fragment) (z-)axis not initialized\n";
    input_error = true;
  }
  
  if(!input_error && orb_step < 0.) {
    //
    err_log << funame << "orbital (first fragment) discretization step not initialized\n";
    input_error = true;
  }

  if(!input_error && frag_axis < 0 && Structure::type(1) != Molecule::MONOATOMIC) {
    //
    err_log << funame << "second fragment (z-)axis not initialized\n";
    input_error = true;
  }

  if(!input_error && frag_step < 0 && Structure::type(1) != Molecule::MONOATOMIC) {
    //
    err_log << funame << "second fragment discretization step not initialized\n";
    input_error = true;
  }

  if(!input_error && !dist_set.size()) {
    //
    err_log << funame << "sampling distances not initialized\n";
    input_error = true;
  }
  
  if(!input_error && ener_min > 0.) {
    //
    err_log << funame << "energy range not initialized\n";
    input_error = true;
  }

  if(!input_error && low_freq_min < 0.) {
    //
    err_log << funame << "minimal conserved mode frequency not initialized\n";
    input_error = true;
  }

  // sampling distances initialization
  //
  std::vector<double> dist_data(dist_set.size());

  itemp = 0;
  //
  for(std::set<double>::const_reverse_iterator rit = dist_set.rbegin(); rit != dist_set.rend(); ++rit, ++itemp)
    //
    dist_data[itemp] = *rit;

  // energy data storage
  //
  Array<double> ener_data((int)dist_data.size());

  // angular data initialization
  //
  switch(Structure::type(1)) {
    //
  case Molecule::MONOATOMIC:
    //
    itemp = 0;

    break;
    //
  case Molecule::LINEAR:
    //
    itemp = 3;

    break;
    //
  case Molecule::NONLINEAR:
    //
    itemp = 4;
  }

  Array<double> ang_pos(3 + itemp);

  std::list<Array<double> > ang_list;

  int ix, iy, iz;

  // orbital (1st fragment) polar angle discretization
  //
  const int orb_polar_num  = std::ceil((orb_polar_max - orb_polar_min) / orb_step + 0.5);

  const double orb_polar_step = (orb_polar_max - orb_polar_min) / double(orb_polar_num - 1);
  
  // 2nd fragment polar angle discretization
  //
  const int frag_polar_num = std::ceil((frag_polar_max - frag_polar_min) / frag_step + 0.5);

  const double frag_polar_step = (frag_polar_max - frag_polar_min) / double(frag_polar_num - 1);
  
  // non-linear 2nd fragment psi angle discretization
  //
  const int frag_psi_num   = std::floor((frag_psi_max - frag_psi_min) / frag_step + 0.5);

  const double frag_psi_step = (frag_psi_max - frag_psi_min) / (double)frag_psi_num;

  // orbital (1st fragment) polar angle cycle
  //
  for(int opi = 0; opi < orb_polar_num; ++opi) {
    //
    const double orb_polar_angle = orb_polar_min + (double)opi * orb_polar_step;
    
    const int    orb_azim_num    = std::ceil((orb_azim_max - orb_azim_min) / orb_step * std::sin(orb_polar_angle) + 0.5);

    const double orb_azim_step   = (orb_azim_max - orb_azim_min) / (double)orb_azim_num;

    // orbital (1st fragment) azimuth angle cycle
    //
    for(int oai = 0; oai < orb_azim_num; ++oai) {
      //
      double orb_azim_angle = orb_azim_min + (double)oai * orb_azim_step;

      iz = orb_axis;
      ix = (orb_axis + 1) % 3;
      iy = (orb_axis + 2) % 3;

      ang_pos[iz] = std::cos(orb_polar_angle);
      ang_pos[ix] = std::sin(orb_polar_angle) * std::cos(orb_azim_angle);
      ang_pos[iy] = std::sin(orb_polar_angle) * std::sin(orb_azim_angle);

      // monoatomic 2nd fragment angular state setting
      //
      if(Structure::type(1) == Molecule::MONOATOMIC) {
	//
	ang_list.push_back(ang_pos);
	
	continue;
      }

      // 2nd fragment polar angle cycle
      //
      for(int fpi = 0; fpi < frag_polar_num; ++fpi) {
	//
	const double frag_polar_angle = frag_polar_min + (double)fpi * frag_polar_step;
	
	const int    frag_azim_num = std::ceil((frag_azim_max - frag_azim_min) / frag_step * std::sin(frag_polar_angle) + 0.5);

	const double frag_azim_step = (frag_azim_max - frag_azim_min) / (double)frag_azim_num;

	// 2nd fragment azimuth angle cycle
	//
	for(int fai = 0; fai < frag_azim_num; ++fai) {
	  //
	  const double frag_azim_angle = frag_azim_min + (double)fai * frag_azim_step;

	  iz = frag_axis;
	  ix = (frag_axis + 1) % 3;
	  iy = (frag_axis + 2) % 3;

	  // linear 2nd fragment angular state setting
	  //
	  if(Structure::type(1) == Molecule::LINEAR) {
	    //
	    ang_pos[iz + 3] = std::cos(frag_polar_angle);
	    ang_pos[ix + 3] = std::sin(frag_polar_angle) * std::cos(frag_azim_angle);
	    ang_pos[iy + 3] = std::sin(frag_polar_angle) * std::sin(frag_azim_angle);

	    ang_list.push_back(ang_pos);

	    continue;
	  }

	  for(int fsi = 0; fsi < frag_psi_num; ++fsi) {
	    //
	    const double frag_psi_angle = frag_psi_min + (double)fsi * frag_psi_step;

	    // non-linear second fragment angular state setting
	    //
	    double euler [3];

	    euler[0] = frag_polar_angle;
	    euler[1] = frag_azim_angle;
	    euler[3] = frag_psi_angle;

	    euler2quat(euler, (double*)ang_pos + 3);
	    
	    ang_list.push_back(ang_pos);
	    //
	  }// 2nd fragment psi angle cycle
	  //
	}// 2nd fragment azimuth angle cycle
	//
      }// 2dn fragment polar angle cycle
      //
    }// orbital (1st fragment) azimuth angle cycle
    //
  }// orbital (1st fragment) polar angle cycle
  
  
  std::vector<Array<double> > ang_data(ang_list.size(), Array<double>(ang_pos.size()));

  itemp = 0;
  
  for(std::list<Array<double> >::const_iterator it = ang_list.begin(); it != ang_list.end(); ++it, ++itemp)
    //
    ang_data[itemp] = *it;

  ang_list.clear();

  // working and scratch directories
  //
  if(mpi_rank && !input_error) {
    //
    // node name
    //
    std::ostringstream node_name;

    node_name << "/node" << mpi_rank;
  
    // set working directory
    //
    input_error |= System::make_dir(work_dir, funame);

    stemp = work_dir + node_name.str();

    input_error |= System::make_dir(stemp, funame);

    if(!input_error && chdir(stemp.c_str())) {
      //
      std::cerr << funame << stemp << ": chdir: " << strerror(errno) << "\n";

      input_error = true;
    }

    // set scratch directory
    //
    input_error |= System::make_dir(scratch_dir, funame);
  
    scratch_dir += node_name.str();
    input_error |= System::make_dir(scratch_dir, funame);

    if(!input_error)
      //
      Molpro::set_scratch_dir(scratch_dir);
  }

  // synchronize the error 
  //
  for(int n = 0; n < mpi_size; ++n) {
    //
    if(mpi_rank == n) {
      //
      btemp = input_error;
    }

    MPI_Bcast(&btemp, 1, MPI_INT, n, MPI_COMM_WORLD);

    input_error |= btemp;
  }

  if(input_error) {
    //
    MPI_Finalize();

    return 0;
  }

  // signal handing
  //
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

  MPI_Status stat;

  int ray_index;
  
  /**********************************************************************************
   ********************************** MASTER ****************************************
   **********************************************************************************/
  
  if(!mpi_rank) {
    //
    int count = 0, work = 0;

    int stop_request_reported = 0;
    
    int calculation_completion_reported = 0;

    std::set<int> low_freq_hit, high_freq_hit, high_ener_hit, check_hit, direct_hit, fail_hit;
    
    std::vector<const char*> host(mpi_size);

    for(int i = 0; i < mpi_size; ++i)
      //
      host[i] = (const char*)host_buff + 100 * i;
  
    std::ofstream to(out_file_name.c_str());

    if(!to || stop_request) {
      //
      if(!to)
	//
	std::cerr << funame << "cannot open " << out_file_name << " file for output: exiting ...\n";

      if(stop_request) {
	//
	std::cerr << funame << "stop request has been received: exitting ...\n";

	stop_request_reported = 1;
      }
      
      for(int node = 1; node < mpi_size; ++node)
	//
	MPI_Send(0, 0, MPI_INT, node, STOP_TAG, MPI_COMM_WORLD);
    }
    // sending sampling requests to working nodes
    //
    else {
      //
      to << ang_data.size() << "    " << dist_data.size() << "   " << ang_pos.size() << "\n";

      for(int d = 0; d < dist_data.size(); ++d) {
	//
	if(!(d % 12))
	  //
	  to << "\n";
      
	to << std::setw(10) << dist_data[d];
      }

      to << "\n\n";
    
      // start work
      //
      for(work = 1; work < mpi_size; ++work, ++count) {
	//
	if(count == ang_data.size())
	  //
	  break;
      
	// send new job request
	//
	MPI_Send(&count, 1, MPI_INT, work, RUN_TAG, MPI_COMM_WORLD);
      }

      // stop unused nodes
      //
      for(int i = work; i < mpi_size; ++i)
	//
	MPI_Send(&count, 1, MPI_INT, i, STOP_TAG, MPI_COMM_WORLD);

      /*****************************************************************************
       ************************ SENDING & RECEIVING LOOP ***************************
       *****************************************************************************/
      //
      while(work > 1) {
	//
	// get ray index
	//
	MPI_Recv(&ray_index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

	const int node = stat.MPI_SOURCE;

	const int tag  = stat.MPI_TAG;

	Array<char> mess;

	switch(tag) {
	  //
	case FAIL_TAG:
	  //
	  std::cout <<  node << "-th node (" << host[node] << ") FAILED\n\n";

	  fail_hit.insert(ray_index);
	  
	  //MPI_Send(&count, 1, MPI_INT, node, STOP_TAG, MPI_COMM_WORLD);

	  --work;

	  continue;
	  //
	case HIGH_ENER_TAG:                      // ray sampling: energy went beyond the requested energy range
	  //
	  high_ener_hit.insert(ray_index);
	  
	  continue;
	  //
	case LOW_FREQ_TAG:              	// ray sampling: low frequency has encountered
	  //
	  low_freq_hit.insert(ray_index);

	  continue;
	  //
	case HIGH_FREQ_TAG:              	// ray sampling: high frequency minimum has encountered
	  //
	  high_freq_hit.insert(ray_index);

	  std::cout <<  node << "-th node (" << host[node] << "): " << ray_index << "-th ray sampling: ABSTRACTION\n";

	  MPI_Recv(&itemp, 1, MPI_INT, node, HIGH_FREQ_TAG, MPI_COMM_WORLD, &stat);

	  mess.resize(itemp);
	  
	  MPI_Recv(mess, mess.size(), MPI_CHAR, node, HIGH_FREQ_TAG, MPI_COMM_WORLD, &stat);

	  std::cout << mess;
	  
	  continue;
	  //
	case DIR_TAG:               	        // ray sampling: direct (not-extrapolated)
	  //
	  direct_hit.insert(ray_index);

	  continue;
	  //
	case CHECK_TAG:	                        // ray sampling to check
	  //
	  check_hit.insert(ray_index);
	  
	  std::cout <<  node << "-th node (" << host[node] << "): " << ray_index << "-th ray sampling: CHECK\n";

	  // receive message
	  //
	  MPI_Recv(&itemp, 1, MPI_INT, node, CHECK_TAG, MPI_COMM_WORLD, &stat);

	  mess.resize(itemp);
	  
	  MPI_Recv(mess, mess.size(), MPI_CHAR, node, CHECK_TAG, MPI_COMM_WORLD, &stat);

	  std::cout << mess;

	  // receive energies
	  //
	  MPI_Recv(&itemp, 1, MPI_INT, node, CHECK_TAG, MPI_COMM_WORLD, &stat);

	  ener_data.resize(itemp);
	  
	  MPI_Recv(ener_data, ener_data.size(), MPI_DOUBLE, node, CHECK_TAG, MPI_COMM_WORLD, &stat);
	  
	  MPI_Recv(&itemp, 1, MPI_INT, node, CHECK_TAG, MPI_COMM_WORLD, &stat);

	  std::cout << std::setw(7) << "R, bohr";
	  
	  for(int i = itemp; i < ener_data.size(); ++i)
	    //
	    std::cout << std::setw(7) << dist_data[i];

	  std::cout << "\n" << std::setw(7) << "E, kcal";
	  
	  for(int i = itemp; i < ener_data.size(); ++i) {
	    //
	    dtemp =  ener_data[i] / Phys_const::kcal;

	    std::cout << std::setw(7);
	    
	    if(dtemp > 0.1 || dtemp < -0.1) {
	      //
	      std::cout << dtemp;
	    }
	    else
	      //
	      std::cout << "<0.1";
	  }
	  
	  std::cout << "\n" << std::endl;

	  continue;
	}
	
	// calculation succeeded
	//
	MPI_Recv(&itemp, 1, MPI_INT, node, RUN_TAG, MPI_COMM_WORLD, &stat);

	ener_data.resize(itemp);
	  
	MPI_Recv(ener_data, ener_data.size(), MPI_DOUBLE, node, RUN_TAG, MPI_COMM_WORLD, &stat);
	
	// save new energies
	//
	sigprocmask(SIG_SETMASK, &block_sig, &old_sig);

	to << std::left << std::setw(15) << ray_index <<  std::right;

	for(int i = 0; i < ang_data[ray_index].size(); ++i)
	  //
	  to << std::setw(15) << ang_data[ray_index][i];

	to << "\n" << ener_data.size();

	for(int e = 0; e < ener_data.size(); ++e) {
	  //
	  if(!(e % 8))
	    //
	    to << "\n";
	
	  to << std::setw(15) << ener_data[e] / Phys_const::kcal;
	}

	to << "\n\n";
	
	sigprocmask(SIG_SETMASK, &old_sig, 0);
	
	// stop sampling
	//
	if(stop_request || count == ang_data.size()) {
	  //
	  if(stop_request && !stop_request_reported) {
	    //
	    std::cout << funame << "stop request has been received: rounding up ...\n\n";

	    stop_request_reported = 1;
	  }
	  if(count == ang_data.size() && !calculation_completion_reported) {
	    //
	    std::cout << funame << "calculation request completed: rounding up ...\n\n";

	    calculation_completion_reported = 1;
	  }
	  
	  MPI_Send(&count, 1, MPI_INT, node, STOP_TAG, MPI_COMM_WORLD);

	  --work;
	  
	  continue;
	}
	// restart failed ray samplings
	//
	else if(fail_hit.size()) {
	  //
	  itemp = *fail_hit.begin();
	  
	  MPI_Send(&itemp, 1, MPI_INT, node, RUN_TAG, MPI_COMM_WORLD);

	  fail_hit.erase(itemp);
	}
	// next ray sampling
	//
	else {
	  //
	  // start new job
	  //	  
	  MPI_Send(&count, 1, MPI_INT, node, RUN_TAG, MPI_COMM_WORLD);

	  ++count;
	}
      }
    
      // samplings to check
      //
      if(check_hit.size()) {
	//
	std::cout << "samplings to check: " << check_hit.size() << "\n";
	//
	itemp = 0;
	//
	for(std::set<int>::const_iterator hit = check_hit.begin(); hit != check_hit.end(); ++hit, ++itemp) {
	  //
	  if(itemp && !(itemp % 20))
	    //
	    std::cout << "\n";
	     
	  std::cout << std::setw(6) << *hit;
	}
	std::cout << "\n\n";
      }
      
      // direct (not-extrapolated) samplings
      //
      if(direct_hit.size()) {
	//
	std::cout << "direct samplings: " << direct_hit.size() << "\n";
	//
	itemp = 0;
	//
	for(std::set<int>::const_iterator hit = direct_hit.begin(); hit != direct_hit.end(); ++hit, ++itemp) {
	  //
	  if(itemp && !(itemp % 20))
	    //
	    std::cout << "\n";
	     
	  std::cout << std::setw(6) << *hit;
	}
	std::cout << "\n\n";
      }
      
      // high energy samplings
      //
      if(high_ener_hit.size()) {
	//
	std::cout << "high energy samplings: " << high_ener_hit.size() << "\n";
	//
	itemp = 0;
	//
	for(std::set<int>::const_iterator hit = high_ener_hit.begin(); hit != high_ener_hit.end(); ++hit, ++itemp) {
	  //
	  if(itemp && !(itemp % 20))
	    //
	    std::cout << "\n";
	     
	  std::cout << std::setw(6) << *hit;
	}
	std::cout << "\n\n";
      }
      
      // low frequency samplings
      //
      if(low_freq_hit.size()) {
	//
	std::cout << "low frequency samplings: " << low_freq_hit.size() << "\n";
	//
	itemp = 0;
	//
	for(std::set<int>::const_iterator hit = low_freq_hit.begin(); hit != low_freq_hit.end(); ++hit, ++itemp) {
	  //
	  if(itemp && !(itemp % 20))
	    //
	    std::cout << "\n";
	     
	  std::cout << std::setw(6) << *hit;
	}
	std::cout << "\n\n";
      }
      
      // high frequency minimum samplings
      //
      if(high_freq_hit.size()) {
	//
	std::cout << "high frequency minimum (abstraction) samplings: " << high_freq_hit.size() << "\n";
	//
	itemp = 0;
	//
	for(std::set<int>::const_iterator hit = high_freq_hit.begin(); hit != high_freq_hit.end(); ++hit, ++itemp) {
	  //
	  if(itemp && !(itemp % 20))
	    //
	    std::cout << "\n";
	     
	  std::cout << std::setw(6) << *hit;
	}
	std::cout << "\n\n";
      }
      
      std::cout << funame << "done\n";
    }
  }  
  /**********************************************************************************
   *********************************** SLAVE ****************************************
   **********************************************************************************/
  //
  else {
    //
    int fail_count = 0;
    
    // open log file
    //
    IO::log.open("molpro.log");

    IO::log.precision(3);
    
    // reactive complex
    //
    std::vector<Atom> rc(Structure::size());

    itemp = 0;

    for(int f = 0; f < 2; ++f)
      //
      for(int i = 0; i < Structure::size(f); ++i, ++itemp)
	//
	rc[itemp] = Structure::fragment(f)[i];

    // generalized coordinates
    //
    Dynamic::Coordinates dc;

    for(int i = 0; i < 4; ++i)
      //
      if(!i) {
	//
	dc.write_ang_pos(0)[i] = 1.;
      }
      else
	//
	dc.write_ang_pos(0)[i] = 0.;

    /*************************************************************************
     ****************************** WORK LOOP ********************************
     *************************************************************************/
    //
    while(1) {
      //
      MPI_Recv(&ray_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

      // stop calculation
      //
      if(stat.MPI_TAG == STOP_TAG) {
	//
	// sending high energy and low frequency hits ##
	//
	IO::log << funame << "cleaning up and exitting\n";

	System::call_exe("/bin/rm", "-Rf", scratch_dir.c_str(), (char*) 0);

	break;
      }
   
      oss.str("");
      
      oss << ray_index << "-th ray sampling";

      IO::log << "\n";
      
      IO::Marker ray_marker(oss.str().c_str());
      
      ang_pos = ang_data[ray_index];

      IO::log << "angular orientation:";

      for(int i = 0; i < ang_pos.size(); ++i)
	//
	IO::log << "   " << ang_pos[i];

      IO::log << "\n";
      
      // orient second fragment
      //
      for(int i = 0; i < Structure::pos_size(1); ++i)
	//
	dc.write_ang_pos(1)[i] = ang_pos[i + 3];

      // relative atomic coordinates in the laboratory (1st fragment) reference frame
      //
      std::vector<D3::Vector> rel_pos [2];

      // fragment-associated transitional and conservative modes
      //
      Lapack::Matrix frag_mode [2];
      int            tran_size [2];

      for(int frag = 0; frag < 2; ++frag) {
	//
	dc.set_rel_pos(frag, rel_pos[frag]);

	frag_mode[frag].resize(3 * Structure::size(frag));

	frag_mode[frag] = 0.;
	
	switch(Structure::type(frag)) {
	  //
	case Molecule::MONOATOMIC:
	  //
	  itemp = 0;
	  
	  break;
	  //
	case Molecule::LINEAR:
	  //
	  itemp = 2;

	  break;
	  //
	case Molecule::NONLINEAR:
	  //
	  itemp = 3;
	}
	tran_size[frag] = 3 + itemp;

	// transitional modes setting
	//
	for(int a = 0; a < Structure::size(frag); ++a) {
	  //
	  const double mass_sqrt = std::sqrt(Structure::fragment(frag)[a].mass());

	  // translational modes
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    frag_mode[frag](i + 3 * a, i) = mass_sqrt;

	  if(Structure::type(frag) == Molecule::MONOATOMIC)
	    //
	    continue;
	
	  // find the most colinear axis with the linear fragment ...
	  //
	  int imax = -1;

	  dtemp = -1.;
      
	  if(Structure::type(frag) == Molecule::LINEAR)
	    //
	    for(int i = 0; i < 3; ++i) 
	      //
	      if(dtemp < std::fabs(dc.ang_pos(frag)[i])) {
		//
		dtemp = std::fabs(dc.ang_pos(frag)[i]);

		imax = i;
	      }

	  itemp = 0;

	  // ... and exclude it from consideration
	  //
	  for(int i = 0; i < 3; ++i) {
	    //
	    if(i == imax) {
	      //
	      ++itemp;

	      continue;
	    }
	    
	    frag_mode[frag]((i + 2) % 3 + 3 * a, i - itemp + 3) =  mass_sqrt * rel_pos[frag][a][(i + 1) % 3];
	    frag_mode[frag]((i + 1) % 3 + 3 * a, i - itemp + 3) = -mass_sqrt * rel_pos[frag][a][(i + 2) % 3];
	  }
	}

	// orthogonalize conserved modes to transitional ones
	//
	for(int m = 0; m < frag_mode[frag].size(); ++m) {
	  //
	  if(m >= tran_size[frag]) {
	    //
	    dtemp = 1.;

	    // for conserved mode find initial orth with the smallest projection to already diagonalized subspace
	    //
	    for(int i = 0; i < frag_mode[frag].size(); ++i) {
	      //
	      if(dtemp > vdot(frag_mode[frag].row(i))) {
		//
		dtemp = vdot(frag_mode[frag].row(i));

		itemp = i;
	      }
	    }
	    frag_mode[frag](itemp, m) = 1.;
	  }
	
	  // all modes
	  //
	  for(int j = 0; j < m; ++j)
	    //
	    ::orthogonalize(&frag_mode[frag](0, m), &frag_mode[frag](0, j), frag_mode[frag].size());

	  ::normalize(&frag_mode[frag](0, m), frag_mode[frag].size());
	}
      }

      int con_size [2];
      
      for(int f = 0; f < 2; ++f)
	//
	con_size[f] = frag_mode[f].size() - tran_size[f];

      itemp = con_size[0] + con_size[1];
      
      Lapack::Matrix con_mode(3 * Structure::size(), itemp);

      con_mode = 0.;

      for(int f = 0; f < 2; ++f)
	//
	for(int m = 0; m < con_size[f]; ++m)
	  //
	  for(int i = 0; i < frag_mode[f].size(); ++i) {
	    //
	    dtemp = frag_mode[f](i, m + tran_size[f]);
	    
	    if(!f) {
	      //
	      con_mode(i, m) = dtemp;
	    }
	    else
	      //
	      con_mode(i + frag_mode[0].size(), m + con_size[0]) = dtemp;
	  }
      
      // conservative modes coefficients
      //
      Array<double> con_mod_coef(con_mode.size2());

      con_mod_coef = 0.;

      // atomic coordinates at the start
      //
      for(int a = 0; a < rc.size(); ++a)
	//
	for(int i = 0; i < 3; ++i)
	  //
	  if(a < Structure::size(0)) {
	    //
	    rc[a][i] = rel_pos[0][a][i];
	  }
	  else
	    //
	    rc[a][i] = rel_pos[1][a - Structure::size(0)][i] + dist_data[0] * ang_pos[i];
      
      // remove old wave function file
      //
      Molpro::remove_wfu();

      int ener_size = 0;

      int geom_relax_start = 0, geom_relax = 0;
      
      Lapack::Vector con_freq;

      Lapack::Vector cshift;

      std::ostringstream check_mess;

      int high_freq_report = 0;

      //int high_ener = 0;

      /********************************************************************************
       ******************************** DISTANCE CYCLE ********************************
       ********************************************************************************/
      //
      for(int d = 0; d < dist_data.size(); ++d) {
	//
	const double& dist = dist_data[d];
	
	double& ener_value = ener_data[d];

	/*
	// energy extrapolation
	//
	if(last_d) {
	  //
	  // linear extrapolation
	  //
	  if(last_d < 2) {
	    //
	    ener_value = ener_data[last_d] + ener_grad * (dist - dist_data[last_d]);
	  }
	  // quadratic extrapolation
	  //
	  else {
	    //
	    ener_value = 0.;
	    
	    for(int i = 0; i < 3; ++i) {
	      //
	      dtemp = dist_data[last_d - i];
	      
	      ener_value += ener_data[last_d - i]
		//
		* (dist  - dist_data[last_d - (i + 1) % 3]) * (dist  - dist_data[last_d - (i + 2) % 3])
		//
		/ (dtemp - dist_data[last_d - (i + 1) % 3]) / (dtemp - dist_data[last_d - (i + 2) % 3]);
	    }
	  }

	  IO::log << "distance [bohr] = " << std::setw(4) << dist << ",  energy [kcal/mol] = "
	    //
		  << std::setw(6) << ener_value / Phys_const::kcal << "\n";

	  continue;
	}
	*/
	  
	// cm shift reactive complex atoms coordinates update
	//
	if(d) {
	  //
	  dtemp = dist - dist_data[d - 1];
	  
	  for(int a = Structure::size(0); a < rc.size(); ++a)  
	    //
	    for(int i = 0; i < 3; ++i)
	      //
	      rc[a][i] += dtemp * ang_pos[i];
	}
	
	// molpro caculation
	//
	try {
	  //
	  Array<double> ener(1);

	  // initial wave function guess
	  //
	  if(!d)
	    //
	    Molpro::pot(rc, ener, Molpro::WFU | Molpro::GUESS);

	  // molpro energy for rigid fragments
	  //
	  if(!geom_relax || d < 3) {
	    //
	    Molpro::pot(rc, ener, Molpro::WFU);
	  }
	  // molpro energy with geometry relaxation correction
	  //
	  else {
	    //
	    if(geom_relax++ == 1) {
	      //
	      geom_relax_start = d;
	      
	      IO::log << "\nswitching to geometry relaxation calculation...\n";

	      for(int a = 0; a < rc.size(); ++a)
		//
		IO::log << rc[a] << "\n";

	      IO::log << "\n";
	    }
	    
	    Molpro::pot(rc, ener, Molpro::WFU | Molpro::RELAX);

	    // mass-weighted hessian from molpro output
	    //
	    Lapack::SymmetricMatrix hess = Molpro::hessian  (rc.size());

	    // gradient cartesian components from molpro output
	    //
	    Lapack::Vector grad = Molpro::gradients(rc.size());

	    /*****************************************************************************
	     * converting mass-weighted hessian to atomic units, hess /= Phys_const::amu, *
	     *                                                                            *
	     * seems to be unnecessary: it has already been done in contrast with molpro  *
	     *                                                                            *
	     * assertion                                                                  *
	     ******************************************************************************/
	    
	    // mass-weighted gradients
	    //
	    for(int a = 0; a < rc.size(); ++a) {
	      //
	      dtemp = std::sqrt(rc[a].mass());

	      for(int i = 0; i < 3; ++i)
		//
		grad[i + 3 * a] /= dtemp;
	    }
	  
	    // conserved modes projected gradients
	    //
	    Lapack::Vector con_grad(con_mode.size2());

	    for(int m = 0; m < con_mode.size2(); ++m)
	      //
	      con_grad[m] = ::vdot(&con_mode(0, m), grad, con_mode.size1());

	    // conserved modes projected Hessian
	    //
	    Lapack::SymmetricMatrix con_hess(con_mode.size2());

	    con_hess = 0.;
	
	    for(int m = 0; m < con_mode.size2(); ++m)
	      //
	      for(int n = m; n < con_mode.size2(); ++n)
		//
		for(int i = 0; i < con_mode.size1(); ++i)
		  //
		  for(int j = 0; j < con_mode.size1(); ++j)
		    //
		    con_hess(m, n) += con_mode(i, m) * con_mode(j, n) * hess(i, j);

	    // conserved mode frequencies
	    //
	    con_freq = con_hess.eigenvalues();

	    for(int m = 0; m < con_freq.size(); ++m) {
	      //
	      dtemp = con_freq[m];
	      
	      con_freq[m] = dtemp < 0. ? -std::sqrt(-dtemp) : std::sqrt(dtemp);
	    }

	    /*************************************************************************
	     ************** LOW CONSERVATIVE MODE FREQUENCY EXCEPTION ****************
	     *************************************************************************/
	    //
	    if(con_freq[0] < low_freq_min) {
	      //
	      MPI_Send(&ray_index, 1, MPI_INT, 0, LOW_FREQ_TAG, MPI_COMM_WORLD);

	      // apply, if available, previous conserved modes shift
	      //
	      if(cshift.size())
		//
		for(int a = 0; a < rc.size(); ++a) {
		  //
		  dtemp = std::sqrt(rc[a].mass());
	    
		  for(int i = 0; i < 3; ++i)
		    //
		    rc[a][i] -= cshift[i + 3 * a] / dtemp;
		}
	      
	      // potential energy recalculation
	      //
	      Molpro::pot(rc, ener, Molpro::WFU);
	      
	      ener_value = ener[0] - ener_ref;
	      
	      oss.str("");
	   
	      oss << "low conserved mode frequency\n";

	      dtemp = ener_value - ener_data[d - 1];

	      btemp = 0;
	      
	      if(ener_incr_max > 0. && dtemp > ener_incr_max || ener_decr_max > 0. && -dtemp > ener_decr_max) {
		//
		btemp = 1;
		
		oss << "energy increment/decrement is too large: " << dtemp / Phys_const::kcal << " kcal/mol\n";
	      }
	  
	      oss << "distance [bohr] = " << std::setw(4) << dist << ",  energy [kcal/mol] = " << std::setw(6) << ener_value / Phys_const::kcal;

	      oss << ",  frequencies [1/cm]:";
	
	      for(int m = 0; m < con_freq.size(); ++m)
		//
		oss << "  " << (int)std::floor(con_freq[m] / Phys_const::incm + 0.5);
	      
	      oss << "\n";
	      
	      oss << "geometry [angstrom]:\n";

	      for(int a = 0; a < rc.size(); ++a)
		//
		oss << rc[a] << "\n";

	      if(con_freq.back() < high_freq_min) {
		//
		MPI_Send(&ray_index, 1, MPI_INT, 0, HIGH_FREQ_TAG, MPI_COMM_WORLD);

		itemp = oss.str().size() + 1;
		
		MPI_Send(&itemp, 1, MPI_INT, 0, HIGH_FREQ_TAG, MPI_COMM_WORLD);

		MPI_Send(oss.str().c_str(), itemp, MPI_CHAR, 0, HIGH_FREQ_TAG, MPI_COMM_WORLD);
	      }
		
	      if(btemp) {
		//
		oss << "returning one step back\n";
		
		IO::log << "\n" << oss.str() << "\n";
		
		check_mess << oss.str();
		
		ener_size = d;

		break;
	      }
	      
	      if(cshift.size())
		//
		oss << "applying previous conserved mode shift\n";
	      
	      IO::log << "\n" << oss.str() << "\n";

	      check_mess << oss.str();
		
	      ener_size = d + 1;

	      break;
	      //
	    }// low frequency exception
	  
	    IO::log << "conserved modes frequencies [1/cm]:";
	
	    for(int m = 0; m < con_freq.size(); ++m)
	      //
	      IO::log << "   " << (int)std::floor(con_freq[m] / Phys_const::incm + 0.5);
	      
	    IO::log << "\n";

	    // conserved modes shift
	    //
	    cshift = Lapack::Cholesky(con_hess).invert(con_grad);

	    con_mod_coef -= cshift;
	
	    cshift = con_mode * cshift;

	    // conserved modes shift reactive complex atomic coordinates update
	    //
	    for(int a = 0; a < rc.size(); ++a) {
	      //
	      dtemp = std::sqrt(rc[a].mass());
	    
	      for(int i = 0; i < 3; ++i)
		//
		rc[a][i] -= cshift[i + 3 * a] / dtemp;
	    }

	    // potential energy recalculation
	    //
	    Molpro::pot(rc, ener, Molpro::WFU);
	  
	    // conserved modes current shift output
	    //
	    IO::log << "conserved modes current displacemnts [angstrom]:\n";

	    for(int a = 0; a < rc.size(); ++a) {
	      //
	      Atom atemp = rc[a];
	    
	      dtemp = std::sqrt(atemp.mass());

	      for(int i = 0; i < 3; ++i)
		//
		atemp[i] = cshift[i + 3 * a] / dtemp;

	      IO::log << atemp << "\n";
	    }

	    // overall conserved modes shift output
	    //
	    Lapack::Vector fshift = con_mode * con_mod_coef;
	
	    IO::log << "conserved modes total displacements [angstrom]:\n";

	    for(int a = 0; a < rc.size(); ++a) {
	      //
	      Atom atemp = rc[a];
	    
	      dtemp = std::sqrt(atemp.mass());

	      for(int i = 0; i < 3; ++i)
		//
		atemp[i] = fshift[i + 3 * a] / dtemp;

	      IO::log << atemp << "\n";
	    }

	    // high frequency check
	    //
	    if(high_freq_max > 0. && con_freq.back() > high_freq_max) {
	      //
	      if(!high_freq_report) {
		//
		check_mess << "high conserved mode frequency\n";
		
		high_freq_report = 1;
	      }
	      
	      check_mess << "distance [bohr] = " << dist << ",  energy [kcal/mol] = "<< ener[0] / Phys_const::kcal << ",  frequencies [1/cm]:";
	
	      for(int m = 0; m < con_freq.size(); ++m)
		//
		check_mess << "   " << (int)std::floor(con_freq[m] / Phys_const::incm + 0.5);

	      check_mess << "\n";
	      //
	    }// high frequency check
	    //
	  }// geometry relaxation
	  
	  // energy reference
	  //
	  if(!d && !is_eref) {
	    //
	    ener_ref = ener[0];
	  }
	  ener_value = ener[0] - ener_ref;
	  
	  ener_size = d + 1;

	  if((geom_relax_dist < 0. || dist < geom_relax_dist) && std::fabs(ener_value) > geom_relax_ener && !geom_relax)
	    //
	    geom_relax = 1;
	}
	// molpro failed
	//
	catch(Error::Molpro& mess) {
	  //
	  oss.str("");
	  
	  oss << mess;

	  oss << "geometry [angstrom]:\n";

	  for(int a = 0; a < rc.size(); ++a)
	    //
	    oss << rc[a] << "\n";

	  check_mess << oss.str();

	  IO::log << "\n" << oss.str() << "\n";

	  ener_size = d;

	  if(!d)
	    //
	    ++fail_count;

	  break;
	}

	// energy output
	//
	dtemp = ener_value / Phys_const::kcal;

	IO::log << "distance [bohr] = " << std::setw(4) << dist << ",  energy [kcal/mol] = " << std::setw(6);

	if(dtemp >= -0.1 && dtemp <= 0.1) {
	  //
	  IO::log << "<0.1";
	}
	else
	  //
	  IO::log << dtemp;

	IO::log << std::endl;

	if(!d)
	  //
	  continue;
	
	/******************************************************************************************
	 *************************** LARGE ENERGY INCREMENT EXCEPTION ***************************** 
	 ******************************************************************************************/
	//
	if(ener_incr_max > 0. && ener_value - ener_data[d - 1] > ener_incr_max) {
	  //
	  oss.str("");
	  
	  oss << "large energy increment: molpro seems to converge to the wrong electronic state\n";

	  oss << "distance [bohr] = " << std::setw(4) << dist << ",  energy [kcal/mol] = " << std::setw(6) << ener_value / Phys_const::kcal;

	  if(con_freq.size()) {
	    //
	    oss << ",  frequencies [1/cm]:";
	
	    for(int m = 0; m < con_freq.size(); ++m)
	      //
	      oss << "   " << (int)std::floor(con_freq[m] / Phys_const::incm + 0.5);
	  }
	  oss << "\n";
	  
	  oss << "geometry [angstrom]:\n";

	  for(int a = 0; a < rc.size(); ++a)
	    //
	    oss << rc[a] << "\n";

	  check_mess << oss.str();
	  
	  IO::log << "\n" << oss.str() << "\n";

	  ener_size = d;

	  break;
	}
	
	// check energy decrement
	//
	dtemp = ener_data[d - 1] - ener_value;
	  
	if(ener_decr_max > 0. &&  dtemp > ener_decr_max) {
	  //
	  oss.str("");
	  
	  oss << "large energy decrement\n";

	  oss << "distance [bohr] = " << dist << ", energy [kcal/mol] = " << ener_value / Phys_const::kcal;

	  if(con_freq.size()) {
	    //
	    oss << ", frequencies [1/cm]:";
	
	    for(int m = 0; m < con_freq.size(); ++m)
	      //
	      oss << "   " << (int)std::floor(con_freq[m] / Phys_const::incm + 0.5);
	  }
	  oss << "\n";

	  oss << "geometry [angstrom]:\n";

	  for(int a = 0; a < rc.size(); ++a)
	    //
	    oss << rc[a] << "\n";

	  check_mess << oss.str();

	  IO::log << "\n" << oss.str() << "\n";

	  ener_size = d;
	  
	  break;
	}
	
	/***********************************************************************
	 *********************** ENERGY OUT OF RANGE ***************************
	 ***********************************************************************/
	//
	if(ener_value > ener_max || ener_value < ener_min) {
	  //
	  MPI_Send(&ray_index, 1, MPI_INT, 0, HIGH_ENER_TAG, MPI_COMM_WORLD);

	  /*
	  // proceeding with linear extrapolation
	  //
	  IO::log << "\nenergy out of range: proceeding with linear extrapolation\n\n";

	  ener_size = dist_data.size();

	  dtemp = (ener_value - ener_data[d - 1]) / (dist - dist_data[d - 1]);

	  for(int e = d + 1; e < ener_size; ++e) {
	    //
	    high_ener = 1;
	  
	    ener_data[e] = ener_value + dtemp * (dist_data[e] - dist);
	  }
	  */
	  
	  IO::log << "\nenergy out of range\n\n";

	  ener_size = d;
	  
	  break;
	}//
	//
      }// distance cycle

      // node failed
      //
      if(fail_count >= fail_count_max) {
	//
	MPI_Send(&ray_index,  1, MPI_INT, 0, FAIL_TAG, MPI_COMM_WORLD);

	break;
      }
	
      // check message
      //
      if(check_mess.str().size()) {
	//
	MPI_Send(&ray_index,  1, MPI_INT, 0, CHECK_TAG, MPI_COMM_WORLD);

	itemp = check_mess.str().size() + 1;

	MPI_Send(&itemp, 1, MPI_INT, 0, CHECK_TAG, MPI_COMM_WORLD);

	MPI_Send(check_mess.str().c_str(), itemp, MPI_CHAR, 0, CHECK_TAG, MPI_COMM_WORLD);

	// energy data
	//
	MPI_Send(&ener_size, 1, MPI_INT, 0, CHECK_TAG, MPI_COMM_WORLD);

	MPI_Send(ener_data, ener_size, MPI_DOUBLE, 0, CHECK_TAG, MPI_COMM_WORLD);

	MPI_Send(&geom_relax_start, 1, MPI_INT, 0, CHECK_TAG, MPI_COMM_WORLD);
      }
      
      if(ener_size == dist_data.size())
	//
	MPI_Send(&ray_index, 1,      MPI_INT, 0, DIR_TAG, MPI_COMM_WORLD);

      MPI_Send(&ray_index,   1,      MPI_INT, 0, RUN_TAG, MPI_COMM_WORLD);
      
      MPI_Send(&ener_size,   1,      MPI_INT,    0, RUN_TAG, MPI_COMM_WORLD);
      
      MPI_Send(ener_data, ener_size, MPI_DOUBLE, 0, RUN_TAG, MPI_COMM_WORLD);
      //
    }// work loop
    //
  }// slave

  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();
  return 0;
}
