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

#include "harmonic.hh"
#include "key.hh"
#include "harm_pot.hh"
#include "io.hh"
#include "units.hh"

#include <fstream>

int main (int argc, char* argv [])
{
  const char funame [] = "main: ";

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;

  if (argc < 2) {
    std::cout << funame  << "usage: harm_fit input_file\n";
    return 1;
  }

  // base name
  std::string base_name = argv[1];
  if(base_name.size() >= 4 && !base_name.compare(base_name.size() - 4, 4, ".inp", 4))
    base_name.resize(base_name.size() - 4);

  // layout
  std::vector<int> vtemp(2);
  vtemp[0] = 3;
  vtemp[1] = 4;
  Configuration::State::layout.set(vtemp);

  // input
  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    throw Error::Input();
  }

  // input parameters
  int rdim = -1;
  int qdim = -1;

  std::vector<double> dist_grid;
  int                 dist_size  = 0;
  double              dist_min   = -1.;
  double              dist_max   = -1.;

  std::vector<Configuration::State>    vertex;
  std::vector<Array<double> >       ener_data;

  std::string xdir;

  EnergyConverter convert;

  KeyGroup HarmonicExpansionFit;
 
  Key  rdim_key("OrbitalExpansionSize" );
  Key  qdim_key("FragmentExpansionSize");
  Key  dmin_key("DistanceMin[bohr]"    );
  Key  dmax_key("DistanceMax[bohr]"    );
  Key dsize_key("DistanceSize"         );
  Key dgrid_key("DistanceGrid[bohr]"   );
  Key   log_key("LogOutput"            );
  Key  vert_key("VertexData"           );
  Key  ener_key("EnergyData"           );
  Key escal_key("EnergyScale[kcal/mol]");
  Key  xdir_key("ExpansionDirectory"  );

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
    // vertex data
    else if(vert_key == token) {
      if(vertex.size()) {
        std::cerr << funame << token << ": allready initialized\n";
        throw Error::Init();
      }      

      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      std::ifstream vertex_in(stemp.c_str());
      if(!vertex_in) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }

      vertex_in >> itemp;
      vertex.resize(itemp);
      for(int v = 0; v < vertex.size(); ++v)
	vertex_in >> vertex[v];

      if(!vertex_in) {
	std::cerr << funame << token << ": vertex file corrupted\n";
	throw Error::Input();
      }
    }
    // energy data
    else if(ener_key == token) {
      if(ener_data.size()) {
        std::cerr << funame << token << ": allready initialized\n";
        throw Error::Init();
      }      

      if(!(from >> stemp)) {
        std::cerr << funame << token << ": corrupted\n";
        throw Error::Input();
      }
      std::getline(from, comment);

      std::ifstream ener_in(stemp.c_str());
      if(!ener_in) {
        std::cerr << funame << token << ": cannot open " << stemp << " file\n";
        throw Error::Input();
      }

      ener_in >> itemp;
      ener_data.resize(itemp);
      for(int v = 0; v < ener_data.size(); ++v) {
	ener_in >> itemp >> itemp;
	ener_data[v].resize(itemp);
	for(int d = 0; d < ener_data[v].size(); ++d)
	  ener_in >> ener_data[v][d];
      }

      if(!ener_in) {
	std::cerr << funame << token << ": energy data file is corrupted\n";
	throw Error::Input();
      }
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
      std::getline(from, comment);

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
      std::getline(from, comment);

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
      std::getline(from, comment);

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
    // orbital expansion dimension
    else if(rdim_key == token) {
      if(rdim >= 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> rdim)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(rdim < 1) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // fragment expansion dimension
    else if(qdim_key == token) {
      if(qdim >= 0) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> qdim)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(qdim < 1) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    // energy conversion scale
    else if(escal_key == token) {
      if(convert.scale > 0.) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }

      if(!(from >> convert.scale)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(convert.scale <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      
      convert.scale *= Phys_const::kcal;
    }
    // harmonic expansion directory
    else if(xdir_key == token) {
      if(!(from >> xdir)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
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

  // default log output	
  if(!IO::log.is_open()) {
    stemp = base_name + ".log";
    IO::log.open(stemp.c_str());
    if(!IO::log) {
      std::cerr << funame << token << ": cannot open " << stemp << " file\n";
      throw Error::Open();
    }
  }

  // vertex data
  if(!vertex.size()) {
    std::cerr << "no vertex data\n";
    throw Error::Init();
  }

  // energy data
  if(ener_data.size() != vertex.size()) {
    std::cerr << "no proper energy data\n";
    throw Error::Init();
  }

  // distance grid
  if(!dist_grid.size() && (dist_size < 2 || dist_min <= 0. || dist_max <= 0. || dist_max <= dist_min)) {
    std::cerr << funame << "distance grid not initialized properly\n";
    throw Error::Range();
  }

  if(!dist_grid.size()) {
    dtemp = std::pow(dist_max / dist_min , 1./double(dist_size - 1));

    double dist = dist_max;
    for(int d = 0; d < dist_size; ++d, dist /= dtemp)
      dist_grid.push_back(dist);
  }

  // energy data vs distance grid
  for(int v = 0; v < ener_data.size(); ++v)
    if(ener_data[v].size() > dist_grid.size()) {
      std::cerr << funame << "energy data inconsistent with the distance grid\n";
      throw Error::Init();
    }

  // expansion dimensions
  if(rdim < 0 || qdim < 0) {
    std::cerr << funame << "expansion dimensions not initialized properly\n";
    throw Error::Init();
  }

  // harmonic expansion
  ConstSharedPointer<HarmonicExpansion> expansion [2];
  itemp = 0;
  for(int pack = 0; pack < 2; itemp += expansion[pack++]->size()) {
    // harmonic expansion data file
    std::ostringstream xname;
    if(xdir.size())
      xname << xdir << "/";
    xname << "hex_" << rdim - pack << "_" << qdim << ".vec";

    std::ifstream xin(xname.str().c_str());
    if(!xin) {
      std::cerr << funame << "cannot open " << xname.str() << " file\n";
      throw Error::Open();
    }

    expansion[pack].init(new HarmonicExpansion(xin, rdim - pack, qdim));
  }

  const int xsize = itemp;

  // harmonic expansion on the vertex grid
  Lapack::Matrix vx_mat(vertex.size(), xsize);

  itemp = 0;
  for(int  pack = 0; pack < 2; itemp += expansion[pack++]->size()) {

#pragma omp parallel for default(shared) schedule(dynamic)      

    for(int v = 0; v < vertex.size(); ++v) {
      for(int x = 0; x < expansion[pack]->size(); ++x)
	vx_mat(v, x + itemp) = (*expansion[pack])(x, vertex[v]);
    }
  }

  Lapack::Matrix hex_coef(xsize, dist_grid.size());

  // distance cycle
  IO::log << std::setprecision(3);
  for(int dist = 0; dist < dist_grid.size(); ++dist) {

    // right-hand-side (energy) vector
    Lapack::Vector rhs_vec(xsize);
    rhs_vec = 0.;

    // left-hand-side (correlation) matrix
    Lapack::SymmetricMatrix lhs_mat(xsize);
    lhs_mat = 0.;

    // vertex cycle
    for(int v = 0; v < vertex.size(); ++v)
      if(ener_data[v].size() > dist) {
	// energy transformed
	dtemp = convert(ener_data[v][dist]);

	// harmonic expansion cycle	
#pragma omp parallel for default(shared) schedule(dynamic)      

	for(int x = 0; x < xsize; ++x) {
	  rhs_vec[x] +=  dtemp * vx_mat(v, x);

	  for(int y = x; y < xsize; ++y)
	    lhs_mat(x, y) += vx_mat(v, x) * vx_mat(v, y);
	}// harmonic expansion cycle
      }// vertex cycle

    Lapack::Vector vtemp = Lapack::Cholesky(lhs_mat).invert(rhs_vec);
    hex_coef.column(dist) = vtemp;
    Lapack::Vector fit_ener  = vx_mat * vtemp;

    double rms_dev = 0.;
    double max_dev = 0.;
    double ener_max, ener_min;
    int    vmax, vmin;

    // vertex cycle
    itemp = 0;

    std::ostringstream corr_name;
    corr_name << base_name << ".corr." << dist;
    std::ofstream corr_out(corr_name.str().c_str());

    for(int v = 0; v < vertex.size(); ++v)
      if(ener_data[v].size() > dist) {
	// energy transformed
	dtemp = convert(ener_data[v][dist]); 

	corr_out << std::setw(15) << dtemp / Phys_const::kcal 
		 << std::setw(15) << fit_ener[v] / Phys_const::kcal
		 << "\n";

	// energy maximum
	if(!itemp || dtemp > ener_max)
	  ener_max = dtemp;

	// energy minimum
	if(!itemp || dtemp < ener_min) {
	  ener_min = dtemp;
	  vmin = v;
	}

	dtemp -= fit_ener[v];

	// root mean square deviation
	rms_dev += dtemp * dtemp;

	// maximal deviation
	dtemp = dtemp >= 0. ? dtemp : -dtemp;
	if(!itemp || dtemp > max_dev) {
	  max_dev = dtemp;
	  vmax = v;
	}

	++itemp;
      }// vertex cycle

    if(itemp)
      rms_dev = std::sqrt(rms_dev / (double)itemp);

    IO::log << "distance [bohr] = " << dist_grid[dist] << "\n"
	    << "max deviation [kcal/mol] = " 
	    << max_dev  / Phys_const::kcal << "\n"
	    << "max deviation vertex  = " 
	    << vmax << "\n"
	    << "max deviation vertex (scaled) energy [kcal/mol] = " 
	    << convert(ener_data[vmax][dist]) / Phys_const::kcal << "\n"
	    << "max deviation vertex (scaled) energy fit [kcal/mol] = " 
	    << fit_ener[vmax] / Phys_const::kcal << "\n"
	    << "max deviation vertex unscaled energy [kcal/mol] = " 
	    << ener_data[vmax][dist] / Phys_const::kcal << "\n"
	    << "max deviation vertex unscaled energy fit [kcal/mol] = " 
	    << convert(fit_ener[vmax], 1) / Phys_const::kcal << "\n"
	    << "(scaled) rms deviation [kcal/mol] = "
	    << rms_dev  / Phys_const::kcal << "\n"
	    << "(scaled) energy max [kcal/mol] = "  
	    << ener_max / Phys_const::kcal << "\n"
	    << "unscaled energy max [kcal/mol] = "  
	    << convert(ener_max, 1) / Phys_const::kcal << "\n"
	    << "(scaled) energy min [kcal/mol] = " 
	    << ener_min / Phys_const::kcal << "\n"
	    << "(scaled) energy min fit [kcal/mol] = " 
	    << fit_ener[vmin] / Phys_const::kcal << "\n"
	    << "unscaled energy min [kcal/mol] = " 
	    << convert(ener_min, 1) / Phys_const::kcal << "\n"
	    << "unscaled energy min fit [kcal/mol] = " 
	    << convert(fit_ener[vmin], 1) / Phys_const::kcal << "\n"
	    << std::endl;

  }// distance cycle
    
  // output
  std::ofstream fit_out((base_name + ".out").c_str());

  // monom dimensions
  fit_out << rdim << " " << qdim << "\n"; 

  // distance grid
  fit_out << dist_grid.size() << "\n";
  for(int dist = 0; dist < dist_grid.size(); ++dist)
    fit_out << dist_grid[dist] << " ";
  fit_out << "\n";

  // energy conversion scale
  fit_out << convert.scale << "\n";

  fit_out << std::setprecision(10);
  // energy expansion coefficients
  for(int i = 0; i < hex_coef.size1(); ++i) {
    for(int j = 0; j < hex_coef.size2(); ++j)
      fit_out << hex_coef(i, j) << " ";
    fit_out << "\n";
  }

  return 0;
}

