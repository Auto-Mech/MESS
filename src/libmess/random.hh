/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef RANDOM_HH
#define RANDOM_HH

#include "error.hh"

namespace Random {

  void   init      ();
  void   init      (int);
  void   send_seed (int);
  double flat      ();             // rand_flat ();
  double norm      ();             // rand_norm ();
  double exp       ();             // rand_exp ();
  void   orient    (double*, int); // rand_orient (double*, int);
  double vol       (double*, int); 

  void spherical_layer (double* vec, int dim, double rmin, double rmax) throw(Error::General);
}


#endif
