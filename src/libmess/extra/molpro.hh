/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef MOLPRO_HH
#define MOLPRO_HH

#include "atom.hh"

#include <string>
#include <iostream>
#include <vector>

namespace Molpro {
  // initialization
  void init (std::istream&) throw(Error::General);
  bool isinit ();

  // molpro potential energy
  void pot (const std::vector<Atom>&, Array<double>&, int flags = 0) throw(Error::General);

  void set_scratch_dir (const std::string& s);
  void remove_wfu      ();

}

#endif
