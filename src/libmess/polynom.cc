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

#include <cmath>
#include <omp.h>

#include "polynom.hh"
#include "multindex.hh"
#include "permutation.hh"

/************************************************************************************************
 ********************************** SYMMETRIC POLYNOMIAL ****************************************
 ************************************************************************************************/

void Polynom::SymPol::init(int dim, int tax, int pax) 
{
  const char funame [] = "Polynom::SymPol::init: ";

  _init = true;

  if(dim < 1 || tax < 1 || pax < 1) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
  
  _dim = dim;
  _term_order_max = tax;
  _poly_order_max = pax;

  int  itemp;
  double dtemp;

  // setting single term maps
  std::vector<std::vector<int> > term_multi_map;
  term_multi_map.push_back(std::vector<int>(dimension()));

  std::vector<int> term_order_map;
  term_order_map.push_back(0);

  for(int d = 0; d < dimension(); ++d) {

    const int term_max = term_multi_map.size();
    _term_max_map.push_back(term_max);

    if(d == dimension() - 1) {
      _term_new_map.resize(term_max);
      _term_pow_map.resize(term_max);
      }

    for(int term = 0; term < term_max; ++term) {

      itemp = poly_order_max() - term_order_map[term];
      itemp = itemp < term_order_max() ? itemp : term_order_max();
      const int pow_max = itemp;

      if(d == dimension() - 1) {
	_term_new_map[term] = term_multi_map.size() - term_max;
	_term_pow_map[term]  = pow_max;
      }

      for(int i = 0; i < pow_max; ++i) {
	term_order_map.push_back(term_order_map[term] + i + 1);
	term_multi_map.push_back(term_multi_map[term]);
	term_multi_map.back()[d] = i + 1;
      }    
    }
  }

  _term_val_map_size = term_order_map.size();

  itemp = 1;
  for(int i = 0; i < dimension(); ++i)
    itemp *= term_order_max() + 1;

  std::vector<int> one_term_map(itemp, -1);
  _one_poly_map.resize(itemp, -1);

  for(int term = 0; term < term_multi_map.size(); ++term)
    one_term_map[multi2one(term_multi_map[term], term_order_max() + 1)] = term;

  // setting polynomial maps
  int poly_ind = 0;

  for(MultiIndex pi(dimension(), term_order_max(), MultiIndex::ORDERED); !pi.end(); ++pi)
    if(pi.cumulative_index() <= poly_order_max()) {
      //  one-dimensional to polynomial map 
      _one_poly_map[multi2one(pi.base(), term_order_max() + 1)] = poly_ind;
      
      std::set<std::vector<int> > mono_pool;
      int perm_group_size = 0;
      for(Permutation perm(dimension()); !perm.end(); ++perm, ++perm_group_size)
	mono_pool.insert(perm(pi.base()));

      // term indices
      std::vector<int> term;
      for(std::set<std::vector<int> >::const_iterator sit = mono_pool.begin(); sit != mono_pool.end(); ++sit) {
	itemp = one_term_map[multi2one(*sit, term_order_max() + 1)];
	if(itemp < 0) {
	  std::cerr << funame << "cannot find the term in the monomial map\n";
	  throw Error::Range();
	}
	term.push_back(itemp);
      }

      if(perm_group_size % term.size()) {
	std::cerr << funame << "number of monomials does not divide the number of permutations\n";
	throw Error::Range();
      }

      // polynomial to single term map
      _poly_term_map.push_back(term);
      // polynomial to symmetry factor map
      _poly_factor.push_back(perm_group_size / term.size());
      
      ++poly_ind;
    }

  std::cout << "initialization of the symmetric polynomial of " << dim << "-th dimension:\n"
	    << "   monomials   number  = " << _term_val_map_size << "\n"
	    << "   polynomials number  = " << poly_ind << "\n"
	    << "done\n";
}

void Polynom::SymPol::update_map (const std::vector<double>& args, std::vector<double>& poly_val_map) const 
{
  const char funame [] = "Polynom::SymPol::update_map: ";

  _isinit();

  int    itemp;
  double dtemp;

  if(dimension() != args.size()) {
    std::cerr << funame << "dimensions mismatch";
    throw Error::Range();
  }

  // updating term value map
  std::vector<double> term_val_map(_term_val_map_size);
  term_val_map[0] = 1.;

  std::vector<std::vector<double> > power_term(dimension(), std::vector<double>(term_order_max()));

#ifndef DEBUG
#pragma omp parallel for default(shared) private(dtemp) schedule(static)
#endif

  for(int d = 0; d < dimension(); ++d) {
    dtemp = args[d];
    for(int i = 0; i < term_order_max(); ++i, dtemp *= args[d])
      power_term[d][i] = dtemp;
  }

  for(int d = 0; d < dimension(); ++d) {

    const int term_max = _term_max_map[d];

#ifndef DEBUG
#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic, 1)
#endif

    for(int term = 0; term < term_max; ++term) {

      // updating terms values for d-th dimension
      const int pow_max = _term_pow_map[term];
      const int base    = _term_new_map[term] + term_max;

      dtemp = term_val_map[term];
      for(int i = 0; i < pow_max; ++i)
	term_val_map[base + i] = dtemp * power_term[d][i];
    }    
  }

  // updating values of polynomials
  poly_val_map.resize(_poly_term_map.size());

#ifndef DEBUG
#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic, 10)
#endif

  for(int p = 0; p < _poly_term_map.size(); ++p) {
    dtemp = 0.;
    for(std::vector<int>::const_iterator vit = _poly_term_map[p].begin(); 
	vit != _poly_term_map[p].end(); ++vit)
      dtemp += term_val_map[*vit];
    dtemp *= (double)_poly_factor[p];
    poly_val_map[p] = dtemp;	
  }
}

double Polynom::SymPol::operator() (const std::multiset<int>& pr, const std::vector<double>& poly_val_map) const 
{
  const char funame [] = "Polynom::SymPol::operator(): ";

  _isinit();

  int itemp;

  if(pr.size() != dimension()) {
    std::cerr << funame << "dimensions mismatch\n";
    throw Error::Range();
  }

  // converting multi-dimensional index to one-dimensional linear index
  itemp = multi2one(pr, term_order_max() + 1);

  if(itemp >= _one_poly_map.size()) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
  
  // converting linear index to actual symmetric polynomial map index 
  itemp = _one_poly_map[itemp];

  if(itemp < 0) {
    std::cerr << funame << "polinomial is not in the map\n";
    throw Error::Range();
  }

  return poly_val_map[itemp];
}

void Polynom::SymPol::_isinit () const 
{
  const char funame [] = "Polynom::SymPol::_isinit: ";

  if(!_init) {
    std::cerr << funame << "not initialized\n";
    throw Error::Range();  
  }
}

