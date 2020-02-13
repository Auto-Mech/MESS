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

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<list>

#include "opt.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"
#include "random.hh"
#include "lapack.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "zopt: ";

  using namespace Opt;
  
  if(argc < 2) {
    //
    std::cout << "usage: opt input_file\n";

    return 0;
  }

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;

  // base name
  //
  std::string base_name = argv[1];
  
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    //
    base_name.resize(base_name.size() - 4);

  // log stream
  //
  stemp = base_name + ".log";
  
  IO::log.open(stemp.c_str());
  
  if(!IO::log) {
    //
    std::cerr << funame << ": cannot open log file " << stemp << " for writing\n";
    
    throw Error::Input();
  }

  // out stream
  //
  stemp = base_name + ".out";
  
  IO::out.open(stemp.c_str());
  
  if(!IO::out) {
    //
    std::cerr << funame << ": cannot open output file " << stemp << " for writing\n";
    
    throw Error::Input();
  }

  /********************************************************************************************
   ************************************* INPUT PARAMETERS *************************************
   ********************************************************************************************/

  // z-matrix
  //
  Coord::ZMat zmat;

  // fluxional modes and their sampling limits
  //
  ZOpt::mode_t flux_modes;

  // non-fluxional (spectator) modes and their sampling limits
  //
  ZOpt::mode_t spec_modes;

  // non-fluxional modes initial values
  //
  typedef std::map<std::string, double> spec_t;

  spec_t spec_init;

  // fluxional modes samplings number
  //
  int flux_num = 0;

  // non-fluxional (tight) modes samplings number
  //
  int spec_num = 0;

  // temperatures list
  //
  typedef std::list<double> temp_t;
  
  temp_t temperature;

  // zero energy
  //
  double ground_ener = 0.;
  
  /********************************************************************************************
   ******************************************* INPUT ******************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key      zmat_key("Z-matrix"              );
  Key       pot_key("Potential"             );
  Key      temp_key("Temperature[K]"        );
  Key  flux_mod_key("LooseModes"            );
  Key  spec_mod_key("TightModes"            );
  Key  flux_num_key("LooseSamplings"        );
  Key  spec_num_key("TightSamplings"        );
  Key  opt_step_key("OptimizationStep"      );
  Key   opt_num_key("OptimizationIteration" );
  Key  grad_tol_key("GradientTolerance"     );
  Key      diff_key("DifferentiationStep"   );
  Key    ground_key("GroundEnergy[kcal/mol]");
  
  // input stream
  //
  std::ifstream from(argv[1]);
  
  if(!from && base_name == argv[1]) {
    //
    stemp = base_name + ".inp";
    
    from.open(stemp.c_str());
  }

  if(!from) {
    //
    std::cerr << funame << "input file " << argv[1] << " is not found\n";

    return 1;
  }

  std::string token, comment;

  // input cycle
  //
  while(from >> token) {
    //
    if(zmat_key == token) {
      //
      if(zmat.isinit()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }
      
      zmat.init(from);
    }
    // potential
    //
    else if(pot_key == token) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      std::getline(from, comment);

      if(stemp == "Harding") {
	//
	Opt::ZOpt::pot.init(new Opt::HardPot(from));
      }
      // unknown potential type
      //
      else {
	//
	std::cerr << funame << token << ": unknown potential type: " << stemp << "\n";

	throw Error::Init();
      }
    }
    // ground energy
    //
    else if(ground_key == token) {
      //
      if(!(from >> ground_ener)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      ground_ener *= Phys_const::kcal;
      
      std::getline(from, comment);
    }
    // number of fluxional samplings
    //
    else if(flux_num_key == token) {
      //
      if(flux_num) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!(from >> flux_num)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(flux_num <= 0) {
	//
	std::cerr << funame << token << ": out of range: " << flux_num << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
    }
    // number of non-fluxional modes samplings
    //
    else if(spec_num_key == token) {
      //
      if(spec_num) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!(from >> spec_num)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(spec_num <= 0) {
	//
	std::cerr << funame << token << ": out of range: " << spec_num << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
    }
    // numeric differentiation step
    //
    else if(diff_key == token) {
      //
      if(!(from >> ZOpt::diff_step)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(ZOpt::diff_step <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << ZOpt::diff_step << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
    }
    // energy gradient tolerance
    //
    else if(grad_tol_key == token) {
      //
      if(!(from >> ZOpt::grad_tol)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(ZOpt::grad_tol <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << ZOpt::grad_tol << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
    }
    // maximal optimization step
    //
    else if(opt_step_key == token) {
      //
      if(!(from >> ZOpt::max_opt_step)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(ZOpt::max_opt_step <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << ZOpt::max_opt_step << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
    }
    // maximal optimization # of iterations
    //
    else if(opt_num_key == token) {
      //
      if(!(from >> ZOpt::max_opt_count)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(ZOpt::max_opt_count <= 0) {
	//
	std::cerr << funame << token << ": out of range: " << ZOpt::max_opt_count << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
    }
    // temperatures list
    //
    else if(temp_key == token) {
      //
      if(temperature.size()) {
	//
	std::cerr << funame << "temperature list has been already defined\n";
	
	throw Error::Input();
      }
      
      IO::LineInput line_input(from);
      
      while(line_input >> dtemp) {
	//
	if(dtemp <= 0.) {
	  //
	  std::cerr << funame << token << ": " << temperature.size() + 1 << "-th temperature out of range\n";

	  throw Error::Range();
	}
	
	temperature.push_back( dtemp * Phys_const::kelv);
      }
      
      if(!temperature.size()) {
	//
        std::cerr << funame << token << ": corrupted\n";

        throw Error::Input();
      }
    }
    // fluxional modes
    //
    else if(flux_mod_key == token) {
      //
      if(!zmat.isinit()) {
	//
	std::cerr << funame << token << ": z-matrix should be initialized first\n";

	throw Error::Init();
      }

      if(flux_modes.size()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }
      
      std::getline(from, comment);
      
      while (1) {
	//
	IO::LineInput lin(from);
	
	if(!(lin >> stemp)) {
	  //
	  std::cerr << funame << token << ": cannot read the variable name (empty line?)\n";

	  throw Error::Input();
	}

	if(stemp == IO::end_key())
	  //
	  break;
	  
	if(!zmat.isvar(stemp)) {
	  //
	  std::cerr << funame << token << ": unknown z-matrix variable: " << stemp << "\n";

	  throw Error::Init();
	}

	if(flux_modes.find(stemp) != flux_modes.end() || spec_modes.find(stemp) != spec_modes.end()) {
	  //
	  std::cerr << funame << token << ": duplicated z-matrix variable: " << stemp << "\n";

	  throw Error::Init();
	}

	double vmin, vmax;
	
	if(!(lin >> vmin >> vmax))
	  //
	  switch(zmat.type(stemp)) {
	    //
	  case Coord::DISTANCE:
	    //
	    std::cerr << funame << token << ": " << stemp
		      << ": sampling limits should be explicitely defined for distance variables\n";

	    throw Error::Range();

	  case Coord::ANGLE:
	    //
	    IO::log << IO::log_offset << funame << token << ": " <<  stemp << ": WARNING: default sampling limits [0, 180] will be used\n";

	    vmin = 0.;

	    vmax = 180.;

	    break;

	  case Coord::DIHEDRAL:
	    //
	    IO::log << IO::log_offset << funame << token << ": " <<  stemp << ": WARNING: default sampling limits [0, 360] will be used\n";

	    vmin = 0.;

	    vmax = 360.;

	    break;
	  }
	    
	if(vmin < 0. || vmax <= vmin) {
	  //
	  std::cerr << funame << token << ": " << stemp << ": sampling limits out of range: " << vmin << ", " << vmax << "\n";

	  throw Error::Range();
	}

	
	switch(zmat.type(stemp)) {
	  //
	case Coord::DISTANCE:
	  //
	  dtemp = Phys_const::angstrom;

	  break;

	case Coord::ANGLE:
	  //
	  if(vmax > 180.) {
	    //
	    std::cerr << funame << token << ": " << stemp << ": sampling limits out of range: " << vmin << ", " << vmax << "\n";

	    throw Error::Range();
	  }

	  dtemp = M_PI / 180.;

	  break;
	  
	case Coord::DIHEDRAL:
	  //
	  if(vmax > 360.) {
	    //
	    std::cerr << funame << token << ": " << stemp << ": sampling limits out of range: " << vmin << ", " << vmax << "\n";

	    throw Error::Range();
	  }

	  dtemp = M_PI / 180.;

	  break;
	}

	vmin *= dtemp;

	vmax *= dtemp;
	
	flux_modes[stemp] = std::make_pair(vmin, vmax - vmin);
      }
    }
    // spectator modes
    //
    else if(spec_mod_key == token) {
      //
      if(!zmat.isinit()) {
	//
	std::cerr << funame << token << ": z-matrix should be initialized first\n";

	throw Error::Init();
      }

      if(spec_modes.size()) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      std::getline(from, comment);

      while (1) {
	//
	IO::LineInput lin(from);
	
	if(!(lin >> stemp)) {
	  //
	  std::cerr << funame << token << ": cannot read the variable name (empty line?)\n";

	  throw Error::Input();
	}

	if(stemp == IO::end_key())
	  //
	  break;
	  
	if(!zmat.isvar(stemp)) {
	  //
	  std::cerr << funame << token << ": unknown z-matrix variable: " << stemp << "\n";

	  throw Error::Init();
	}

	if(spec_modes.find(stemp) != spec_modes.end() || flux_modes.find(stemp) != flux_modes.end()) {
	  //
	  std::cerr << funame << token << ": duplicated z-matrix variable: " << stemp << "\n";

	  throw Error::Init();
	}

	for(int i = 0; i < 3; ++i) {
	  //
	  if(!(lin >> dtemp)) {
	    //
	    std::cerr << funame << token << ": cannot read " << stemp << " ";

	    switch(i) {
	      //
	    case 0:
	      //
	      std::cerr << "initial value";

	      break;

	    case 1:
	      //
	      std::cerr << "lower limit";

	      break;

	    case 2:
	      //
	      std::cerr << "upper limit";

	      break;
	    }
	    
	    std::cerr << "\n";

	    throw Error::Input();
	  }
	
	  if(dtemp < 0.) {
	    //
	    std::cerr << funame << token << ": " << stemp << ": out of range: " << dtemp << "\n";

	    throw Error::Range();
	  }
	
	  switch(zmat.type(stemp)) {
	    //
	  case Coord::DISTANCE:
	    //
	    dtemp *= Phys_const::angstrom;

	    break;

	  case Coord::ANGLE:
	    //
	    if(dtemp >= 180.) {
	      //
	      std::cerr << funame << token << ": " << stemp << ": out of range: " << dtemp << "\n";

	      throw Error::Range();
	    }

	    dtemp *= M_PI / 180.;

	    break;
	  
	  case Coord::DIHEDRAL:
	    //
	    if(dtemp >= 360.) {
	      //
	      std::cerr << funame << token << ": " << stemp << ": out of range: " << dtemp << "\n";

	      throw Error::Range();
	    }

	    dtemp *= M_PI / 180.;

	    break;
	  }

	  switch(i) {
	    //
	  case 0:
	    //
	    spec_init[stemp] = dtemp;

	    break;

	  case 1:
	    //
	    if(dtemp > spec_init[stemp]) {
	      //
	      std::cerr << funame << token << ": " << stemp << " lower limit, " << dtemp << ", bigger than the initial value: " << spec_init[stemp] << "\n";

	      throw Error::Range();
	    }

	    spec_modes[stemp].first = dtemp;

	    break;
	  
	  case 2:
	    //
	    if(dtemp < spec_init[stemp]) {
	      //
	      std::cerr << funame << token << ": " << stemp << " upper limit, " << dtemp << ", smaller than the initial value: " << spec_init[stemp] << "\n";

	      throw Error::Range();
	    }

	    spec_modes[stemp].second = dtemp;

	    break;
	  }
	}
      }
    }
    // unknown key
    //
    else if(IO::skip_comment(token, from)) {
      //
      std::cerr << funame << "unknown keyword " << token << "\n";
      
      Key::show_all(std::cerr);
      
      std::cerr << "\n";
      
      throw Error::Init();
    }
  }// input cycle

  // checking
  //
  if(!zmat.isinit()) {
    //
    std::cerr << funame << "zmatrix not initialized\n";

    throw Error::Init();
  }
  
  if(!Opt::ZOpt::pot) {
    //
    std::cerr << funame << "potential not initialized\n";

    throw Error::Init();
  }

  if(!flux_modes.size()) {
    //
    std::cerr << funame << "no fluxional modes specified\n";

    throw Error::Init();
  }

  if(!spec_modes.size()) {
    //
    std::cerr << funame << "no spectator modes specified\n";

    throw Error::Init();
  }
  
  // check that all zmatrix variables initialized
  //
  for(int r = 0; r < zmat.size(); ++r)
    //
    for(int v = 0; v < zmat[r].size(); ++v) {
      //
      const std::string& var = zmat[r].var(v);
  
      if(flux_modes.find(var) == flux_modes.end() && spec_modes.find(var) == spec_modes.end()) {
	//
	std::cerr << funame << "z-matrix variable " << var << " not initialized\n";

	throw Error::Init();
      }
    }

  if(!flux_num) {
    //
    std::cerr << funame << "number of samplings not initialized\n";

    throw Error::Init();
  }

  if(!temperature.size()) {
    //
    std::cerr << funame << "temperature list not initialized\n";

    throw Error::Init();
  }
  
  Random::init();

  Coord::ZData zdata(zmat);

  // non-fluxional modes initial values assignment
  //
  for(spec_t::const_iterator sit = spec_init.begin(); sit != spec_init.end(); ++sit)
    //
    zdata[sit->first] = sit->second;

  // statistical weight
  //
  std::vector<double> weight(temperature.size());

  // statistical weight uncertainty
  //
  std::vector<double> werr(temperature.size());
  
  // anharmonic correction average
  //
  std::vector<double> ahc_mean(temperature.size());

  // anharmonic correction root mean square
  //
  std::vector<double> ahc_rms(temperature.size());
  
  // fluxional modes sampling cycle
  //
  int flux_count = 0;
  
  while(flux_count < flux_num) {
    //
    // fluxional modes assignment
    //
    for(ZOpt::mode_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit)
      //
      zdata[fit->first] = fit->second.first + fit->second.second * Random::flat();

    try {
      //
      Opt::ZOpt zopt(zdata, spec_modes);

      IO::log << IO::log_offset << "Sampling " << flux_count + 1 << "\n";
    
      IO::out << "non-fluxional coordinates at minimum:\n";

      for(ZOpt::mode_t::const_iterator sit = spec_modes.begin(); sit != spec_modes.end(); ++sit)
	//
	IO::out << std::setw(10) << sit->first << std::setw(13) << zopt.zmin()[sit->first] / zmat.conv_factor(sit->first) << "\n";

      IO::out << "\n";
    
      IO::out << "anharmonic correction:\n"
	      << std::setw(5)  << "T, K"
	      << std::setw(13) << "Corr"
	      << std::setw(13) << "Err, %"
	      << "\n";
    
      double relerr; // anharmonic correction relative error

      itemp = 0;

      for(temp_t::const_iterator t = temperature.begin(); t != temperature.end(); ++t, ++itemp) {
	//
	double ahc =  zopt.anharmonic_correction(*t, spec_num, &relerr);

	IO::out << std::setw(5) << *t / Phys_const::kelv
		<< std::setw(13) << ahc
		<< std::setw(13) << relerr * 100.
		<< "\n";
      
	dtemp = std::exp(-(zopt.ener_min() - ground_ener) / *t) * zopt.mass_factor() / zopt.fc_factor();

	weight[itemp]   += dtemp;

	werr[itemp]     += dtemp * dtemp;
      
	ahc_mean[itemp] += dtemp * ahc;

	ahc_rms[itemp]  += dtemp * ahc * ahc;
      }


      IO::out << std::endl;

      ++flux_count;
    }
    catch(ZOpt::NoConv) {
      //
      IO::log << IO::log_offset << "Constrained minimum search failed\n";
    }
    //
  }// fluxional modes sampling cycle

  // statistical weight normalizaiton factor
  //
  dtemp = 2. / std::pow(2. * M_PI, double(flux_modes.size() - 1) / 2.);

  for(ZOpt::mode_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit)
    //
    dtemp *= fit->second.second;

  const double nfac = dtemp;
  
  // normalization and output
  //
  IO::out << "statistical weight in harmonic approximation, anharmonic correction average and its root-mean square deviation:\n"
	  << std::setw(5)  << "T, K"
	  << std::setw(13) << "Weight"
	  << std::setw(13) << "Err, %"
	  << std::setw(13) << "Anharm corr"
	  << std::setw(13) << "Anharm rms"
	  << "\n";
  
  itemp = 0;
  
  for(temp_t::const_iterator t = temperature.begin(); t != temperature.end(); ++t, ++itemp) {
    //
    ahc_mean[itemp] /= weight[itemp];

    ahc_rms[itemp] /= weight[itemp];

    ahc_rms[itemp] = std::sqrt(ahc_rms[itemp] - ahc_mean[itemp] * ahc_mean[itemp]);

    weight[itemp] /= (double)flux_num;

    werr[itemp]   /= (double)flux_num;

    werr[itemp] = std::sqrt((werr[itemp] / weight[itemp] / weight[itemp] - 1.) / (double)flux_num);

    weight[itemp] *= nfac *  std::pow(*t, double(2 * spec_modes.size() + flux_modes.size() + 3) / 2.);

    IO::out << std::setw(5)  << *t / Phys_const::kelv
	    << std::setw(13) << weight[itemp]
	    << std::setw(13) << werr[itemp] * 100.
	    << std::setw(13) << ahc_mean[itemp]
	    << std::setw(13) << ahc_rms[itemp]
	    << "\n";
    
  }
  
  return 0;
}
