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

#include "harm_pot.hh"
#include "harmonic.hh"
#include "slatec.hh"
#include "io.hh"
#include "structure.hh"

#include<fstream>

EnergyConverter                       convert_energy;
ConstSharedPointer<HarmonicExpansion> harmonic_expansion [2];
std::vector<Slatec::Spline>           expansion_coefficient;

double EnergyConverter::operator() (double ener, int back) const
{
  const char funame [] = "EnergyConverter::backward: ";

  if(scale <= 0.)
    return ener;

  double dtemp = ener / scale;
  dtemp = dtemp * dtemp;

  if(dtemp < 1.e-7)
    return ener;
  
  if(back) {
    if(dtemp > 50.) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }

    dtemp = scale * std::sqrt(std::exp(dtemp) - 1.);
  }
  else
    dtemp = scale * std::sqrt(std::log(1. + dtemp));

    if(ener < 0.)
      return -dtemp;
    
    return dtemp;
}

extern "C" void harm_init_ (const char* file)
{
  const char funame [] = "harm_init: ";

  int         itemp;
  double      dtemp;
  std::string stemp;

  // layout
  std::vector<int> lt(2);
  lt[0] = 3;
  lt[1] = 4;

  Configuration::State::layout.set(lt);

  // input
  std::ifstream from(file);
  if(!from) {
    std::cerr << funame << "cannot open " << file << " file";
    std::exit(1);
  }

  // structure
  from >> stemp;
  if(stemp != "Structure") {
    std::cerr << funame << "structure should go first\n";
    std::exit(1);
  }
  std::getline(from, stemp);
  Structure::init(from);

  // monom dimensions
  int rdim, qdim;
  from >> rdim >> qdim; 

  // distance grid
  from >> itemp;
  Array<double> dist_grid(itemp);

  for(int dist = 0; dist < dist_grid.size(); ++dist)
    from >> dist_grid[dist_grid.size() - dist - 1];


  // energy conversion scale
  from >> convert_energy.scale;

  // harmonic expansion vectors
  itemp = 0;
  for(int pack = 0; pack < 2; itemp += harmonic_expansion[pack++]->size()) {
    std::ostringstream xname;
    xname << "hex_" << rdim - pack << "_" << qdim << ".vec";
    std::ifstream xin(xname.str().c_str());
    if(!xin) {
      std::cerr << funame << "cannot open " << xname.str() << " file\n";
      std::exit(1);
    } 
    harmonic_expansion[pack].init(new HarmonicExpansion(xin, rdim - pack, qdim));
  }

  const int xsize = itemp;
  expansion_coefficient.resize(xsize);

  // energy expansion coefficients
  Array<double> vtemp(dist_grid.size());
  for(int x = 0; x < xsize; ++x) {
    for(int d = 0; d < dist_grid.size(); ++d)
      from >> vtemp[dist_grid.size() - d - 1];
    expansion_coefficient[x].init(dist_grid, vtemp, dist_grid.size());
  }

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    std::exit(1);
  }
}

double harm_fit (double dist, const Configuration::State& stat)
{
  const char funame [] = "harm_fit: ";

  int    itemp;
  double dtemp;

  if(dist < expansion_coefficient[0].arg_min() || dist > expansion_coefficient[0].arg_max()) {
    std::cerr << funame << "distance out of range\n";
    std::exit(1);
  }

  double res = 0.;
  itemp = 0;
  for(int  pack = 0; pack < 2; ++pack)
    for(int x = 0; x < harmonic_expansion[pack]->size(); ++x, ++itemp)
	res += (*harmonic_expansion[pack])(x, stat) * expansion_coefficient[itemp](dist);

  return convert_energy(res, 1);
}

extern "C" void no2_ch3_pot_ (const double& nc_dist, const double& ang, double& res)
{
  const char funame [] = "no2_ch3_pot_: ";

  double dtemp;

  // center-of-mass-to-N vector
  D3::Vector nx = Structure::fragment(0)[0];
  nx *= -1.;
  nx.normalize();

  // center-of-mass-to-O vector
  D3::Vector ny = Structure::fragment(0)[1];
  ny.orthogonalize(nx);
  ny.normalize();
  
  dtemp = ang * M_PI / 180.;
  D3::Vector n2z = std::cos(dtemp) * nx + std::sin(dtemp) * ny;
  D3::Vector n2x = std::cos(dtemp) * ny - std::sin(dtemp) * nx;
  D3::Vector n2y = D3::vprod(n2z, n2x);
  
  Quaternion orient(D3::Matrix(Structure::fragment(1)[1], Structure::fragment(1)[2]).transpose() * D3::Matrix(n2x, n2y));
  
  D3::Vector rvec = Structure::fragment(0)[0] + nc_dist * n2z;
  
  double dist = rvec.normalize();
  //std::cout << "dist = " << dist << "\n";

  Configuration::State stat;
  for(int i = 0; i < 3; ++i)
    stat.radius_vector()[i] = rvec[i];

  for(int i = 0; i < 4; ++i)
    stat.orientation()[i] = orient[i];

  res = harm_fit(dist, stat) / Phys_const::kcal;
}

