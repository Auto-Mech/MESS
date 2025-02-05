#ifndef OFFLOAD_HH
#define OFFLOAD_HH

#include "lapack.hh"

namespace Offload {
  //
  namespace Omp {
    //
    Lapack::Vector eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec = 0);
  }

  namespace Cuda {
    //
    class Init {
      //
    public:
      //
      Init (int =0);

      ~Init();
    };
      
    Lapack::Vector eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec = 0);
  }

  namespace Sycl {
    //
    Lapack::Vector eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec = 0);
  }
}

#endif
