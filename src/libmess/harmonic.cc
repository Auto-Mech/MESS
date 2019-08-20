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
#include "io.hh"

/*******************************************************************************************
 ****************** HARMONIC EXPANSION OF THE ANGULAR PART OF THE POTENTIAL *****************
 *******************************************************************************************/

std::ostream& operator<< (std::ostream& to, const HarmonicExpansion& hex)
{
  to << hex._expansion.size() << "\n";
  for(int x = 0; x < hex._expansion.size(); ++x) {
    for(std::map<int, double>::const_iterator it = hex._expansion[x].begin(); it != hex._expansion[x].end(); ++it) 
      to << it->first << " " << it->second << " ";
    to << "\n";
  }
  return to;
}

HarmonicExpansion::HarmonicExpansion (std::istream& from, int rdim, int qdim) 
   
  :  _rmonom(rdim, 3), _qmonom(2 * qdim, 4)
{
  const char funame [] = "HarmonicExpansion::HarmonicExpansion: ";

  
  int         itemp;
  double      dtemp;
  std::string stemp;

  from >> itemp;
  std::getline(from, stemp);

  _expansion.resize(itemp);

  for(int x = 0; x < _expansion.size(); ++x) {
    IO::LineInput lin(from);
    while(lin >> itemp >> dtemp)
      _expansion[x][itemp] = dtemp;
  }

  if(!from) {
    std::cerr << funame << "corrupted\n";
    throw Error::Input();
  }
}

HarmonicExpansion::HarmonicExpansion (const Configuration::DoubleSpaceGroup& symm_group, int rdim, int qdim) 
   
  :  _rmonom(rdim, 3), _qmonom(2 * qdim, 4)
{
  const char funame [] = "HarmonicExpansion::HarmonicExpansion: ";

  static const double tolerance = 1.e-10;

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;

  const int xsize = _rmonom.linear_size() * _qmonom.linear_size();
  IO::log << IO::log_offset << "(" << rdim << ", " << qdim << ")-expansion dimension = " << xsize << std::endl;

  Lapack::Matrix project(xsize);
  project = 0.;

  Lapack::SymmetricMatrix product(xsize);
  product = 0.;

  Lapack::Matrix a(xsize);

  for(int g = 0; g < symm_group.size(); ++g) {

    IO::log << IO::log_offset << std::setw(3) << g + 1 << " out of " << symm_group.size()
	    << " symmetry elements" << std::endl;

    Lapack::Matrix rmatrix = _rmonom(symm_group.rmatrix(g));
    Lapack::Matrix qmatrix = _qmonom(symm_group.qmatrix(g));
      
#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)      

    for(int i = 0; i < xsize; ++i)
      for(int j = 0; j < xsize; ++j) {
	dtemp = rmatrix(i / _qmonom.linear_size(), j / _qmonom.linear_size()) 
	  *     qmatrix(i % _qmonom.linear_size(), j % _qmonom.linear_size());
	a(i, j)        = dtemp;
	project(i, j) += dtemp;
      }

#pragma omp parallel for default(shared) schedule(dynamic)      

    for(int i = 0; i < xsize; ++i)
      for(int j = i; j < xsize; ++j)
	product(i, j) += vdot(a.row(i), a.row(j));
  }

    
  product /= (double)symm_group.size();
  project /= (double)symm_group.size();
  project  = project * product;

  IO::log << IO::log_offset << "diagonalizing symmetry group projector ... ";
  Lapack::Matrix evec;
  Lapack::Vector eval = Lapack::diagonalize(Lapack::SymmetricMatrix(project), product, &evec);
  IO::log << "done" << std::endl;

  // checking eigenvalues
  for(int i = 0; i < xsize; ++i) {
    dtemp = eval[i];
    if(dtemp > - tolerance && dtemp < tolerance)
      continue;

    dtemp -= 1.;
    if(dtemp > - tolerance && dtemp < tolerance)
      continue;

    std::cerr << funame << "not projector\n";
    throw Error::Range();
  }

  itemp = -1;
  for(int i = 0; i < xsize; ++i) {
    if(itemp < 0 && eval[i] > 1. - tolerance) {
      itemp = 0;
      _expansion.resize(xsize - i);
    }
    if(itemp >= 0) {
      for(int j = 0; j < xsize; ++j) {
	dtemp = evec(j, i);
	if(dtemp < -tolerance || dtemp > tolerance)
	  _expansion[itemp][j] = dtemp;
      }
      ++itemp;
    }
  }

  itemp = 0;
  for(int x = 0; x < _expansion.size(); ++x)
    itemp += _expansion[x].size();

  IO::log << IO::log_offset << "symmetric monomials number = "   << _expansion.size()
	  << ";  total number of non-zero terms = " << itemp 
	  << std::endl;
}

double HarmonicExpansion::operator() (int term, const Configuration::State& state) const 
{
  const char funame [] = "HarmonicExpansion::operator(): ";

  double dtemp;

  //check layout
  if(state.size() != 7) {
    std::cerr << funame << "wrong layout\n";
    throw Error::Logic();
  }
  
  if(term < 0 || term >= size()) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  Lapack::Matrix rfactor(_rmonom.rank(), 3);
  for(int i = 0; i < 3; ++i) {
    dtemp = state.radius_vector()[i];
    for(int j = 0; j < _rmonom.rank(); ++j, dtemp *= state.radius_vector()[i])
      rfactor(j, i) = dtemp;
  }

  Lapack::Matrix qfactor(_qmonom.rank(), 4);
  for(int i = 0; i < 4; ++i) {
    dtemp = state.orientation()[i];
    for(int j = 0; j < _qmonom.rank(); ++j, dtemp *= state.orientation()[i])
      qfactor(j, i) = dtemp;
  }

  double res = 0.;

  //#pragma omp parallel for default(shared) private(dtemp) reduction(+: res) schedule(static)

  for(std::map<int, double>::const_iterator it = _expansion[term].begin(); it != _expansion[term].end(); ++it) {
    dtemp = it->second;

    std::vector<int> multi = _rmonom(it->first / _qmonom.linear_size());
    for(int i = 0; i < 3; ++i)
      if(multi[i])
	dtemp *= rfactor(multi[i] - 1, i);

    multi = _qmonom(it->first % _qmonom.linear_size());
    for(int i = 0; i < 4; ++i)
      if(multi[i])
	dtemp *= qfactor(multi[i] - 1, i);

    res += dtemp;
  }

  return res;
}

