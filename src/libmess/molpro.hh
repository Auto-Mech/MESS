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

#ifndef MOLPRO_HH
#define MOLPRO_HH

#include "atom.hh"
#include "lapack.hh"

#include <string>
#include <iostream>
#include <vector>

namespace Molpro {
  //
  extern int mute;
  
  enum {GUESS = 1, WFU = 2, RELAX = 4};
  
  // initialization
  //
  void init (std::istream&);
  bool isinit ();

  // molpro potential energy
  //
  void pot (const std::vector<Atom>&, Array<double>&, int flags = 0);

  bool is_guess ();
  
  void set_scratch_dir (const std::string&);
  
  void remove_wfu ();

  Lapack::SymmetricMatrix hessian (int);

  Lapack::Vector gradients (int);
}

#endif
