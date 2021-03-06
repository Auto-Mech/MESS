/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2016, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef MPACK_HH
#define MPACK_HH

#include "lapack.hh"

namespace Mpack {
  //
  Lapack::Vector dd_eigenvalues (Lapack::SymmetricMatrix, Lapack::Matrix* =0) ;

  Lapack::Vector qd_eigenvalues (Lapack::SymmetricMatrix, Lapack::Matrix* =0) ;
}

#endif
