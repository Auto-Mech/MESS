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

  Lapack::Vector eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec = 0);

  inline Lapack::Vector eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec)
  {
#if defined OFFLOAD_CUDA

    return Cuda::eigenvalues(m, evec);

#elif defined OFFLOAD_SYCL

    return Sycl::eigenvalues(m, evec);

#elif defined OFFLOAD_OMP

    return Omp::eigenvalues(m, evec);

#else

    return m.eigenvalues(evec);

#endif
  }
}

#endif
