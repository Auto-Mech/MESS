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

#include "opt.hh"
#include "key.hh"
#include "units.hh"
#include "io.hh"
#include "random.hh"
#include "lapack.hh"

int main (int argc, char* argv [])
{
  const char funame [] = "opt: ";

  if (argc < 2) {
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
   ************************************* INPUT ************************************************
   ********************************************************************************************/

  KeyGroup Main;

  Key zmat_key("Z-matrix"              );
  Key flux_key("FluxionalModes"        );
  Key spec_key("SpectatorModes"        );
  Key  opt_key("OptimizationParameters");
  Key  pot_key("Potential"             );
  Key samp_key("SamplingNumber"        );

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
  
  Coord::ZMat zmat;
  
  typedef std::map<std::string, std::pair<double, double> > flux_t;
  
  flux_t flux_modes;

  typedef std::map<std::string, double> spec_t;

  spec_t spec_modes;

  Opt::CcOpt opt;

  int samp_size = 0;
  
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
    // fluxional modes
    //
    else if(flux_key == token) {
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

	double vmin, vrange;
	
	if(!(lin >> vmin >> vrange))
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

	    vrange = 180.;

	    break;

	  case Coord::DIHEDRAL:
	    //
	    IO::log << IO::log_offset << funame << token << ": " <<  stemp << ": WARNING: default sampling limits [0, 360] will be used\n";

	    vmin = 0.;

	    vrange = 360.;

	    break;
	  }
	    
	if(vmin < 0. || vrange <= 0.) {
	  //
	  std::cerr << funame << token << ": " << stemp << ": sampling limits out of range: " << vmin << ", " << vrange << "\n";

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
	  if(vrange + vmin > 180.) {
	    //
	    std::cerr << funame << token << ": " << stemp << ": sampling limits out of range: " << vmin << ", " << vrange << "\n";

	    throw Error::Range();
	  }

	  dtemp = M_PI / 180.;

	  break;
	  
	case Coord::DIHEDRAL:
	  //
	  if(vrange + vmin > 360.) {
	    //
	    std::cerr << funame << token << ": " << stemp << ": sampling limits out of range: " << vmin << ", " << vrange << "\n";

	    throw Error::Range();
	  }

	  dtemp = M_PI / 180.;

	  break;
	}

	vmin *= dtemp;

	vrange *= dtemp;
	
	flux_modes[stemp] = std::make_pair(vmin, vrange);

	opt.add_constrain(ConstSharedPointer<Coord::CartFun>(new Coord::Internal(zmat.size(), zmat.sign(stemp))));
      }
    }
    // spectator modes
    //
    else if(spec_key == token) {
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

	if(!(lin >> dtemp)) {
	  //
	  std::cerr << funame << token << ": " << stemp << ": cannot read variable value\n";

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

	spec_modes[stemp] = dtemp;
      }
    }
    // optimizer setup (optional)
    //
    else if(opt_key == token) {
      //
      std::getline(from, comment);

      opt.set(from);
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
	opt.set_pot(ConstSharedPointer<Coord::CartFun>(new Opt::HardPot(from)));
      }
      // unknown potential type
      //
      else {
	//
	std::cerr << funame << token << ": unknown potential type: " << stemp << "\n";

	throw Error::Init();
      }
    }
    // number of samplings
    //
    else if(samp_key == token) {
      //
      if(samp_size) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!(from >> samp_size)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(samp_size <= 0) {
	//
	std::cerr << funame << token << ": out of range: " << samp_size << "\n";

	throw Error::Range();
      }
      
      std::getline(from, comment);
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
    //
  }// input cycle

  // checking
  //
  if(!zmat.isinit()) {
    //
    std::cerr << funame << "zmatrix not initialized\n";

    throw Error::Init();
  }
  
  if(!opt.pot()) {
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

  if(!samp_size) {
    //
    std::cerr << funame << "number of samplings not initialized\n";

    throw Error::Init();
  }

  Lapack::SymmetricMatrix hess(zmat.size() * 3);
  
  Random::init();

  Coord::ZData zdata(zmat);

  for(spec_t::const_iterator sit = spec_modes.begin(); sit != spec_modes.end(); ++sit)
    //
    zdata[sit->first] = sit->second;

  double ener_max, ener_min;

  std::vector<double> fm_max(flux_modes.size()), fm_min(flux_modes.size());
  
  // main output
  //
  for(int r = 0; r < zmat.size(); ++r)
    //
    IO::out << zmat[r].atom().name();

  IO::out << "\n";

  // sampling cycle
  //
  for(int samp = 0; samp < samp_size; ++samp) {
    //
    IO::log << IO::log_offset << "Sampling " << samp + 1 << "\n";
    
    // fluxional modes assignment
    //
    for(flux_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit)
      //
      zdata[fit->first] = fit->second.first + fit->second.second * Random::flat();

    Coord::Cartesian x(zdata);

    opt.execute(x);

    if(!samp) {
      //
      zdata.import(x);

      IO::log << IO::log_offset << "Spectator modes:\n";
      
      IO::log << IO::log_offset
	      << std::setw(10) << "Var"
	      << std::setw(13) << "Original"
	      << std::setw(13) << "new"
	      << "\n";

      for(spec_t::const_iterator sit = spec_modes.begin(); sit != spec_modes.end(); ++sit)
	//
	IO::log << IO::log_offset
		<< std::setw(10) << sit->first
		<< std::setw(13) << sit->second
		<< std::setw(13) << zdata[sit->first]
		<< "\n";
    }      
      
    /****************************************** OUTPUT *********************************************/

    // sampling label
    //
    IO::out << "Sampling " << samp + 1 << "\n";

    // energy (kcal/mol)
    //
    IO::out << "Energy\n";

    double ener =  opt.pot()->evaluate(x) / Phys_const::kcal;

    if(!samp || ener > ener_max) {
      //
      ener_max = ener;

      itemp = 0;
      
      for(flux_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit, ++itemp)
	//
	fm_max[itemp] = zdata[fit->first];
    }
    
    if(!samp || ener < ener_min) {
      //
      ener_min = ener;

      itemp = 0;
      
      for(flux_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit, ++itemp)
	//
	fm_min[itemp] = zdata[fit->first];
    }
    
    IO::out << ener << "\n";

    // geometry (angstrom)
    //
    zmat.cm_shift(x);
    
    IO::out << "Geometry\n" << zmat.size() << "\n\n";
    
    for(int a = 0; a < zmat.size(); ++a) {
      //
      IO::out << std::setw(3) << std::left << zmat[a].atom().name() << std::right;

      for(int i = 0; i < 3; ++i)
	//
	IO::out << std::setw(13) << x.atom_pos(a)[i] / Phys_const::angstrom;

      IO::out << "\n";
    }

    // energy gradient
    //
    IO::out << "Gradient";

    for(int i = 0; i < x.size(); ++i) {
      //
      if(!(i % 3))
	//
	IO::out << "\n";

      IO::out << std::setw(13) << opt.pot()->grad(x, i);
    }

    IO::out << "\n";

    // hessian
    //
    for(int i = 0; i < x.size(); ++i)
      //
      for(int j = i; j < x.size(); ++j)
	//
	hess(i, j) = opt.pot()->hess(x, i, j);

    IO::out << "Hessian\n";

    for(int i = 0; i < x.size(); ++i) {
      //
      for(int j = 0; j < x.size(); ++j)
	//
	IO::out << std::setw(13) << hess(i, j);

      IO::out << "\n";
    }
    
    IO::out << std::endl;
  }

  // minimal energy configuration
  //
  IO::log << IO::log_offset << "Minimal energy configuration:\n";

  itemp = 0;
      
  for(flux_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit, ++itemp)
    //
    zdata[fit->first] = fm_min[itemp];

  Coord::Cartesian xmin(zdata);

  opt.use_constrain = 0;

  opt.execute(xmin);

  IO::log << IO::log_offset << "Energy[kcal/mol] = " << std::setw(13) << opt.pot()->evaluate(xmin) / Phys_const::kcal << "\n";

  // cartesian coordinates
  //
  IO::log << IO::log_offset << "Cartesian coordinates[angstrom]:\n";

  for(int a = 0; a < zmat.size(); ++a) {
    //
    IO::log << IO::log_offset << std::setw(3) << zmat[a].atom().name();

    for(int i = 0; i < 3; ++i)
      //
      IO::log << std::setw(13) << xmin.atom_pos(a)[i] / Phys_const::angstrom;

    IO::log << "\n";
  }

  // hessian
  //
  for(int i = 0; i < hess.size(); ++i)
    //
    for(int j = i; j < hess.size(); ++j)
      //
      hess(i, j) = opt.pot()->hess(xmin, i, j) / std::sqrt(zmat[i / 3].atom().mass() * zmat[j / 3].atom().mass());

  // frequencies
  //
  Lapack::Vector eval = hess.eigenvalues();

  IO::log << IO::log_offset << "Frequencies[1/cm]:";
  
  for(int i = 0; i < eval.size(); ++i) {
    //
    dtemp = eval[i];

    dtemp = dtemp >= 0.? std::sqrt(dtemp) : -std::sqrt(-dtemp);
    
    if(!(i % 3))
      //
      IO::log << "\n" << IO::log_offset;
    
    IO::log << std::setw(13) << dtemp / Phys_const::incm;
  }

  IO::log << "\n";

  // center-of-mass shift
  //
  zmat.cm_shift(xmin);
  
  // inertia matrix
  //
  Lapack::SymmetricMatrix imat(3);

  dtemp = 0.;
  
  for(int a = 0; a < zmat.size(); ++a)
    //
    dtemp += zmat[a].atom().mass() * vdot(xmin.atom_pos(a), 3);

  
  imat = dtemp;

  for(int a = 0; a < zmat.size(); ++a)
    //
    for(int i = 0; i < 3; ++i)
      //
      for(int j = i; j < 3; ++j)
	//
	imat(i, j) -= zmat[a].atom().mass() * xmin.atom_pos(a)[i] * xmin.atom_pos(a)[j];

  // rotational constants
  //
  Lapack::Vector imom = imat.eigenvalues();

  IO::log << IO::log_offset << "Rotational constants[1/cm]:";

  for(int i = 0; i < 3; ++i)
    //
    if(imom[i] > 1.)
      //
      IO::log << std::setw(13) << 0.5 / imom[i] / Phys_const::incm;

  IO::log << "\n";
  
  // z-matrix variables
  //
  zdata.import(xmin);

  IO::log << IO::log_offset << "z-matrix variables (angstrom, degrees):\n";

  for(flux_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit, ++itemp) {
    //
    IO::log << IO::log_offset << std::setw(5) << fit->first << " = ";

    dtemp = zdata[fit->first];
    
    switch(zmat.type(fit->first)) {
      //
    case Coord::DISTANCE:
      //
      dtemp /= Phys_const::angstrom;

      break;

    default:
      //
      dtemp *= 180. / M_PI;
    }

    IO::log << std::setw(13) << dtemp << "\n";
  }

  for(spec_t::const_iterator sit = spec_modes.begin(); sit != spec_modes.end(); ++sit) {
    //
    IO::log << IO::log_offset << std::setw(5) << sit->first << " = ";

    dtemp = zdata[sit->first];
    
    switch(zmat.type(sit->first)) {
      //
    case Coord::DISTANCE:
      //
      dtemp /= Phys_const::angstrom;

      break;

    default:
      //
      dtemp *= 180. / M_PI;
    }

    IO::log << std::setw(13) << dtemp << "\n";
  }
  
  IO::log << "\n";
    
  IO::log << "Maximal energy sampling:\n";

  IO::log << "Energy[kcal/mol] = " << std::setw(13) << ener_max << "\n";

  itemp = 0;
      
  for(flux_t::const_iterator fit = flux_modes.begin(); fit != flux_modes.end(); ++fit, ++itemp) {
    //
    IO::log << IO::log_offset << std::setw(5) << fit->first << " = ";

    dtemp = fm_max[itemp];
    
    switch(zmat.type(fit->first)) {
      //
    case Coord::DISTANCE:
      //
      dtemp /= Phys_const::angstrom;
      
      break;

    default:
      //
      dtemp *= 180. / M_PI;
    }

    IO::log << std::setw(13) << dtemp << "\n";
  }
  
  IO::log << "\n";
    
  return 0;
}
