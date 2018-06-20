#ifndef HARDING_HH
#define HARDING_HH

#include <iostream>

extern "C" void harding_init_        (std::istream&);
extern "C" void harding_update_      (const double*, double*);
extern "C" void harding_pot_         (const double*, double&); 
extern "C" void harding_poly_size_   (int&);
extern "C" void harding_dist_size_   (int&);
extern "C" void harding_atom_size_   (int&);

#endif
