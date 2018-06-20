/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef MPACK_HH
#define MPACK_HH

#include "lapack.hh"

namespace Mpack {
  //
  Lapack::Vector dd_eigenvalues (Lapack::SymmetricMatrix, Lapack::Matrix* =0) throw(Error::General);

  Lapack::Vector qd_eigenvalues (Lapack::SymmetricMatrix, Lapack::Matrix* =0) throw(Error::General);
}

#endif
