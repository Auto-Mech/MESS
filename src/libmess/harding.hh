#ifndef HARDING_HH
#define HARDING_HH

#include <iostream>
#include <map>
#include <vector>
#include <set>

#include "array.hh"
#include "permutation.hh"

class Harding {
  //
  // input parameters
  //
  int                   term_order_max; // maximal order of individual power terms 
  int                   poly_order_max; // overall expansion polynomial order maximum
  Array<double>         coef;           // expansion coefficients
  std::set<Permutation> symm_group;     // symmetry group of allowed permutations

  int         _atom_size;

  // additional bond expansion terms
  //
  std::map<std::set<int>, int>  bond;

  /* generated parameters */

  // dimensionality of the fitting configurational space
  //
  int _dist_size;

  std::vector<std::vector<int> > poly_map;
  std::vector<int>               term_term_map;
  std::vector<int>               term_pow_map;
  std::vector<int>               term_max_map;

  int term_val_map_size;

  int _poly_size; // number of symmetric polynomials

public:

  Harding () : term_order_max(0), poly_order_max(0), _atom_size(0) {}

  void init (std::istream&);

  bool isinit () const { return _atom_size; }

  Harding (std::istream& from) : term_order_max(0), poly_order_max(0), _atom_size(0) { init(from); }

  void update_poly_val (const double*, double*) const;

  double potential (const double*) const;

  int poly_size () const { return _poly_size; }

  int dist_size () const { return _dist_size; }
  
  int atom_size () const { return _atom_size; }
};

extern "C" void harding_init_        (std::istream&);
extern "C" void harding_update_      (const double*, double*);
extern "C" void harding_pot_         (const double*, double&); 
extern "C" void harding_poly_size_   (int&);
extern "C" void harding_dist_size_   (int&);
extern "C" void harding_atom_size_   (int&);

#endif
