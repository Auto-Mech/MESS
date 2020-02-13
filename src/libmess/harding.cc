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

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <cmath>
#include <cstdlib>

#ifdef OPENMP
#include <omp.h>
#endif

#include "harding.hh"
#include "key.hh"
#include "io.hh"
#include "units.hh"

// interatomic distances counted as (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3), ....

Harding harding_global;

void Harding::init (std::istream& from)
{
  Exception::Base funame = "Harding::init: ";

  try {
    
    if(isinit())
      //
      throw funame << "already initialized";

    IO::Marker funame_marker(funame);

    KeyGroup HardingInitGroup;

    Key atom_key("AtomsNumber"     );
    Key term_key("TermOrderMax"    );
    Key poly_key("PolynomOrderMax" );
    Key symm_key("SymmetryGroup"   );
    Key coef_key("CoefficientsFile");
    Key bond_key("BondTerm"        );
    
    int    itemp;
    double dtemp;

    std::vector<std::vector<int> > bond_def;

    std::string line, token, comment;

    std::string     coef_file;

    while(from >> token) {
      //
      // end of input
      //
      if(token == IO::end_key()) {
	//
	std::getline(from, comment);
	//
	break;
      }
      // number of atoms
      //
      else if(token == atom_key) {
	//
	IO::LineInput iss(from);
    
	if(!(iss >> _atom_size))
	  //
	  throw funame << token << ": corrupted";

	if(_atom_size < 2)
	  //
	  throw funame << token << ": out of range: " << _atom_size;

	_dist_size = _atom_size * (_atom_size - 1) / 2;
      }
      // maximal polynomial order
      //
      else if(token == poly_key) {
	//
	IO::LineInput iss(from);

	if(!(iss >> poly_order_max))
	  //
	  throw funame << token << ": corrupted";

	if(poly_order_max <= 0)
	  //
	  throw funame << token << ": out of range: " << poly_order_max;
      }
      // maximal term order
      //
      else if(token == term_key) {
	//
	IO::LineInput iss(from);

	if(!(iss >> term_order_max))
	  //
	  throw funame << token << ": corrupted";

	if(term_order_max <= 0)
	  //
	  throw funame << token << ": out of range: " << term_order_max;
      }
      // bond expansion term
      //
      else if(token == bond_key) {
	//
	IO::LineInput iss(from);

	std::vector<int> bd(3);
	
	for(int i = 0; i < 3; ++i)
	  //
	  if(!(iss >> bd[i]))
	    //
	    throw funame << token << ": corrupted";

	bond_def.push_back(bd);
      }
      // symmetry group
      //
      else if(token == symm_key) {
	//
	IO::LineInput iss(from);
    
	if(!_atom_size)
	  //
	  throw funame << token << ": number of atoms should be set first";

	std::set<Permutation> atom_group;
	
	try {
	  //
	  while(1)
	    //
	    atom_group.insert(Permutation(Permutation::read_orbit(iss), _atom_size));
	} 
	catch(Exception::Eof) {}

	atom_group = permutation_group(atom_group);
	
	for(std::set<Permutation>::const_iterator g = atom_group.begin(); g != atom_group.end(); ++g) {
	  //
	  std::vector<int> gg(_dist_size);
	  
	  for(int i2 = 1; i2 < _atom_size; ++i2)
	    //
	    for(int i1 = 0; i1 < i2; ++i1) {
	      //
	      int j1 = (*g)[i1];
	      
	      int j2 = (*g)[i2];
	      
	      if(j1 > j2)
		//
		std::swap(j1, j2);
	      
	      gg[i2 * (i2 - 1) / 2 + i1] = j2 * (j2 - 1) / 2 + j1;
	    }

	  symm_group.insert(Permutation(gg));
	}
      }
      // expansion coefficients file
      //
      else if(token == coef_key) {
	//
	if(coef_file.size())
	  //
	  throw funame << token << ": already initialized";

	IO::LineInput iss(from);

	if(!(iss >> coef_file))
	  //
	  throw funame << token << ": corrupted";
      }
      // unknown keyword
      //
      else {
	throw funame << "unknown key: " << token << "\n" << Key::show_all();
      }
    }

    if(!_atom_size)
      throw funame << "number of atoms has not been defined";

    if(!poly_order_max)
      throw funame  << "polynomial order maximum has not been initialized";

    if(!term_order_max || term_order_max > poly_order_max)
      term_order_max = poly_order_max;

    if(!symm_group.size())
      symm_group.insert(Permutation(_dist_size));

    // setting single term maps
    std::vector<std::vector<int> > term_multi_map;
    term_multi_map.push_back(std::vector<int>(_dist_size));

    std::vector<int> term_order_map;
    term_order_map.push_back(0);

    int multi_num = _dist_size * term_order_max;

    for(int d = 0; d < _dist_size; ++d) {
      const int term_max = term_multi_map.size();
      term_max_map.push_back(term_max);

      if(d == _dist_size - 1) {
	term_term_map.resize(term_max);
	term_pow_map.resize(term_max);
      }

      for(int term = 0; term < term_max; ++term) {
	itemp = poly_order_max - term_order_map[term];
	itemp = itemp < term_order_max ? itemp : term_order_max;
	const int pow_max = itemp;

	multi_num += pow_max;

	if(d == _dist_size - 1) {
	  term_term_map[term] = term_multi_map.size() - term_max;
	  term_pow_map[term]  = pow_max;
	}

	for(int i = 0; i < pow_max; ++i) {
	  term_order_map.push_back(term_order_map[term] + i + 1);
	  term_multi_map.push_back(term_multi_map[term]);
	  term_multi_map.back()[d] = i + 1;
	}    
      }
    }

    term_val_map_size = term_multi_map.size();

    std::map<std::vector<int>, int> multi_term_map;
    for(int term = 0; term < term_multi_map.size(); ++term)
      multi_term_map[term_multi_map[term]] = term;

    while(multi_term_map.size()) {
      std::vector<int> multi = multi_term_map.begin()->first;
      std::vector<int> term_group;

      for(std::set<Permutation>::const_iterator g = symm_group.begin(); g != symm_group.end(); ++g) {
	std::vector<int> curr = (*g)(multi);
	std::map<std::vector<int>, int>::iterator mp = multi_term_map.find(curr);

	if(mp != multi_term_map.end()) {
	  term_group.push_back(mp->second);
	  multi_term_map.erase(mp);
	}
      }

      if(!term_group.size())
	throw funame << "no terms";

      if(symm_group.size() % term_group.size())
	throw funame << "number of terms in polynomial is not a divider of the symmetry group size";

      poly_map.push_back(term_group);
    }

    // additional bond expansion terms
    //
    for(int b = 0; b < bond_def.size(); ++b) {
      int a0 = bond_def[b][0];
      int a1 = bond_def[b][1];
      if(a0 < 0 || a1 < 0 || a0 >= _atom_size || a1 >= _atom_size || a0 == a1)
	throw funame << "atomic indices in bond term definition out of range";

      if(a0 > a1)
	std::swap(a0, a1);

      int dd = a1 * (a1 - 1) / 2 + a0;
      for(std::map<std::set<int>, int>::const_iterator bit=bond.begin(); bit != bond.end(); ++bit)
	if(bit->first.find(dd) != bit->first.end())
	  throw funame << "bond between " << a0 << "-th and " << a1 << "-th atoms already used in bond expansion";

      std::set<int> bb;
      for(std::set<Permutation>::const_iterator g = symm_group.begin(); g != symm_group.end(); ++g)
	bb.insert((*g)[dd]);

      if(bond_def[b][2] <= term_order_max)
	throw funame << "bond term order should be bigger than term order max";

      bond[bb] = bond_def[b][2];
    }

    int bond_size = 0;
    //
    for(std::map<std::set<int>, int>::const_iterator bit=bond.begin(); bit != bond.end(); ++bit)
      //
      bond_size += bit->second - term_order_max;

    
    _poly_size = poly_map.size() + bond_size;

    // read expansion coefficients
    //
    if(coef_file.size()) {
      //
      std::ifstream cin(coef_file.c_str());
      
      if(!cin)
	throw funame << "cannot open expansion coefficients file " << coef_file;

      coef.resize(_poly_size);

      if(!cin.read(reinterpret_cast<char*>((double*)coef), coef.size() * 8))
	//
	throw funame << coef_file << ": number of coefficients less than expected?";
      
      char c;
      
      if(cin.get(c))
	//
	throw funame << coef_file << ": number of coefficients more than expected?";
    }

    // output
    //
    IO::log << IO::log_offset << "maximal polynomial power        = " << poly_order_max    << "\n"
	    << IO::log_offset << "maximal   monomial power        = " << term_order_max    << "\n"
	    << IO::log_offset << "number of symmetric polynomials = " << poly_map.size()   << "\n"
	    << IO::log_offset << "number of monomial terms        = " << term_val_map_size << "\n";

    if(bond_size)
      //
      IO::log << IO::log_offset << "number of additional bond terms = " << bond_size         << "\n";

    IO::log << IO::log_offset << "number of multiplications       = " << multi_num         << "\n"
	    << IO::log_offset << "number of additions             = " << term_val_map_size - poly_map.size() << "\n"
	    << IO::log_offset << "number of distances             = " << _dist_size         << "\n"
	    << IO::log_offset << "symmetry group size             = " << symm_group.size() << "\n";

#ifdef OPENMP
    itemp = 0;
#pragma omp parallel default(shared) reduction(+: itemp)
    {
#pragma omp master
      {
	IO::log << IO::log_offset << "number of threads               = " << omp_get_num_threads() << "\n";
      }
#pragma omp for schedule(static)
      for(int i = 0; i < 100; ++i)
	itemp += 1;
    }
#endif
  }
  catch(Exception::Base x) {
    x.print();
    std::exit(1);
  }
}

void Harding::update_poly_val (const double* r, double* poly_val) const
{
  const char funame [] = "Harding::update_poly_val: ";

  int    itemp;
  
  double dtemp;

  // updating term value map
  std::vector<double> term_val_map(term_val_map_size);
  term_val_map[0] = 1.;


  for(int d = 0; d < _dist_size; ++d) {
    const int term_max = term_max_map[d];

    std::vector<double> power_term(term_order_max);
    dtemp = r[d];
    for(int i = 0; i < term_order_max; ++i, dtemp *= r[d])
      power_term[i] = dtemp;

    int pow_max, base;

#ifdef OPENMP
    
#pragma omp parallel for default(shared) private(pow_max, base, dtemp) schedule(static)
    
#endif

    for(int term = 0; term < term_max; ++term) {
      // updating terms values for d-th dimension
      pow_max = term_pow_map[term];
      base    = term_term_map[term] + term_max;
      dtemp   = term_val_map[term];

      for(int i = 0; i < pow_max; ++i)
	term_val_map[base + i] = dtemp * power_term[i];
    }
  }

#ifdef OPENMP
#pragma omp parallel for default(shared) private(dtemp) schedule(static)
#endif

  for(int p = 0; p < poly_map.size(); ++p) {
    dtemp = 0.;
    for(std::vector<int>::const_iterator t = poly_map[p].begin(); t != poly_map[p].end(); ++t)
      dtemp += term_val_map[*t];
    poly_val[p] = dtemp;
  }

  // additional bond expansion terms
  int bindex = poly_map.size();
  for(std::map<std::set<int>, int>::const_iterator bit=bond.begin(); bit != bond.end(); ++bit)
    for(int i = term_order_max + 1; i <= bit->second; ++i, ++bindex) {
      dtemp = 0.;
      for(std::set<int>::const_iterator dit = bit->first.begin(); dit != bit->first.end(); ++dit)
	dtemp += std::pow(r[*dit], (double)i);
      poly_val[bindex] = dtemp;
    }
}

double Harding::potential (const double* coord) const
{
  Exception::Base funame = "Harding::potential: ";

  double dtemp;

  try {
    //
    if(!coef.size())
      //
      throw funame << "expansion coefficients have not been initialized";

    Array<double> dist(_dist_size);

    int count = 0;

    for(int a1 = 1; a1 < _atom_size; ++a1) {
      //
      for(int a0 = 0; a0 < a1; ++a0, ++count) {
	//
	double rr = 0.;
	
	for(int i = 0; i < 3; ++i) {
	  //
	  dtemp = coord[a0 * 3 + i] - coord[a1 * 3 + i];
	  
	  rr += dtemp * dtemp;
	}
	
	dist[count] = std::sqrt(rr) / Phys_const::angstrom;
      }
    }

    Array<double> poly_val(_poly_size);
    
    update_poly_val(dist, poly_val);

    double res = 0.;

#ifdef OPENMP
#pragma omp parallel for default(shared) reduction(+: res) schedule(static)
#endif

    for(int i = 0; i < _poly_size; ++i)
      //
      res += poly_val[i] * coef[i];
    
    return res * Phys_const::kcal;
  }
  catch(Exception::Base x) {
    //
    x.print();
    
    std::exit(1);
  }
}

extern "C" void harding_poly_size_ (int& s)
{
  s = harding_global.poly_size();
}

extern "C" void harding_dist_size_ (int& s)
{
  s = harding_global.dist_size();
}

extern "C" void harding_atom_size_ (int& s)
{
  s = harding_global.atom_size();
}

extern "C" void harding_update_ (const double* r, double* poly_val)
{
  harding_global.update_poly_val(r, poly_val);
}

extern "C" void harding_init_ (std::istream& from)
{
  harding_global.init(from);
}

extern "C" void harding_pot_ (const double* coord, double& res) 
{
  res = harding_global.potential(coord);
}
