

#ifndef MPACK_HH
#define MPACK_HH

#include "lapack.hh"

namespace Mpack {
  //
  Lapack::Vector dd_eigenvalues (Lapack::SymmetricMatrix, Lapack::Matrix* =0) throw(Error::General);

  Lapack::Vector qd_eigenvalues (Lapack::SymmetricMatrix, Lapack::Matrix* =0) throw(Error::General);
}

#endif
