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
  void   volume    (double*, int); 

  void spherical_layer (double* vec, int dim, double rmin, double rmax) ;
}


#endif
