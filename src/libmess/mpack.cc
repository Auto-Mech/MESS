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

#include "mpack.hh"

int Mpack::mp_type = Mpack::DOUBLE;

template<>
void Mpack::Rsyev<double> (const char *jobz, const char *uplo, MPACK_INT n, double* a, MPACK_INT lda,
  //
  double* w, double* work, MPACK_INT lwork, MPACK_INT *info)
{
  Lapack::int_t _info;

  dsyev_(*jobz, *uplo, n, a, lda, w, work, lwork, _info);

  *info = _info;
}

template<>
Lapack::Vector Mpack::eigenvalues<double> (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec) 
{
  const char funame [] = "Mpack::eigenvalues<double>: ";

  return mat.eigenvalues(evec);
}
