#include "offload.hh"
#include "io.hh"
#include <cuda_runtime.h>
#include <cusolverDn.h>

#define CUDA_CHECK(err)                                                                            \
    do {                                                                                           \
        cudaError_t err_ = (err);                                                                  \
        if (err_ != cudaSuccess) {                                                                 \
            printf("CUDA error %d at %s:%d\n", err_, __FILE__, __LINE__);                          \
	    fflush(stdout);									   \
            throw std::runtime_error("CUDA error");                                                \
        }                                                                                          \
    } while (0)

#define CUSOLVER_CHECK(err)                                                                        \
    do {                                                                                           \
        cusolverStatus_t err_ = (err);                                                             \
        if (err_ != CUSOLVER_STATUS_SUCCESS) {                                                     \
            printf("cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__);                      \
	    fflush(stdout);									   \
            throw std::runtime_error("cusolver error");                                            \
        }                                                                                          \
    } while (0)

namespace Offload {
  //
  namespace Cuda {
    //
    cusolverDnHandle_t handle = NULL;
    cusolverDnParams_t params = NULL;
  }
}

Offload::Cuda::Init::Init (int gpu)
{
  const char funame [] = "Offload::Cuda::Init::Init: ";
  
  int      gpu_size  = 0;  
  CUDA_CHECK(cudaGetDeviceCount(&gpu_size));
  
  if(gpu_size <= 0) {
    //
    std::cerr << funame << "no GPUs: " << gpu_size << std::endl;

    throw Error::Init();
  }
  
  if(gpu < 0) {
    //
    std::cerr << funame << "gpu id out of range: " << gpu << std::endl;

    throw Error::Range();
  }

  gpu %= gpu_size;
  
  // initialize gpu
  //
  unsigned gpu_flags = 0;
  CUDA_CHECK(cudaGetDeviceFlags(&gpu_flags));

#if CUDART_VERSION >= 12020

  CUDA_CHECK(cudaInitDevice(gpu, gpu_flags, 0));

#endif

  CUDA_CHECK(cudaSetDevice(gpu));

  // initialize cusolver library
  //
  CUSOLVER_CHECK(cusolverDnCreate(&handle));
  CUSOLVER_CHECK(cusolverDnCreateParams(&params));
}

Offload::Cuda::Init::~Init ()
{
  cusolverDnDestroyParams(params);
  cusolverDnDestroy(handle);
  cudaDeviceReset();
}

Lapack::Vector Offload::Cuda::eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec)
{  
  const char funame [] = "Offload::Cuda::eigenvalues: ";

  int itemp;

  if(!m.isinit()) {
    //
    IO::log << IO::log_offset << funame << "not initialized" << std::endl;

    throw Error::Init();
  }

  Lapack::Matrix a = m;
  Lapack::Vector res(m.size());

  int64_t n  = m.size();
  int64_t nn = n * n;

  cudaDataType       type = CUDA_R_64F;
  cusolverEigMode_t  jobz = CUSOLVER_EIG_MODE_NOVECTOR;
  cublasFillMode_t   uplo = CUBLAS_FILL_MODE_LOWER;

  if(evec) {
    //
    jobz = CUSOLVER_EIG_MODE_VECTOR;

    evec->resize(n);
  }
  
  cudaDeviceProp prop;

  cudaGetDeviceProperties(&prop, 0);

  IO::log << IO::log_offset << "gpu name: " << prop.name << std::endl;
  
  cudaGetDevice(&itemp);

  IO::log << IO::log_offset << "gpu id: " << itemp << std::endl;
  
  // initialize matrix on device
  //
  void *matrixOnDevice, *eigenValuesOnDevice; 
  int *infoOnDevice;

  CUDA_CHECK(cudaMalloc(&matrixOnDevice,     nn * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&eigenValuesOnDevice, n * sizeof(double)));
  CUDA_CHECK(cudaMalloc((void **)&infoOnDevice,   sizeof(int)));

  CUDA_CHECK(cudaMemcpy(matrixOnDevice, a, nn * sizeof(double), cudaMemcpyHostToDevice));

  // work space for eigenvalue solver
  //
  void *workSpaceOnDevice = NULL, *workSpaceOnHost = NULL;
  size_t workSizeOnDevice = 0, workSizeOnHost = 0;

  CUSOLVER_CHECK(cusolverDnXsyevd_bufferSize(handle, params, jobz, uplo, n, type, matrixOnDevice, n, type, 
                                             eigenValuesOnDevice, type, &workSizeOnDevice, &workSizeOnHost));

  IO::log << IO::log_offset << "requested work space on device: " << workSizeOnDevice << std::endl;

  IO::log << IO::log_offset << "requested work space on host:   " << workSizeOnHost << std::endl;

  if(workSizeOnDevice)
    //
    CUDA_CHECK(cudaMalloc(&workSpaceOnDevice, workSizeOnDevice));

  if(workSizeOnHost)
    //
    workSpaceOnHost = malloc(workSizeOnHost);

  CUSOLVER_CHECK(cusolverDnXsyevd(handle, params, jobz, uplo, n, type, matrixOnDevice, n, type, eigenValuesOnDevice, 
                                  type, workSpaceOnDevice, workSizeOnDevice, workSpaceOnHost, workSizeOnHost, infoOnDevice));

  int info;
  
  CUDA_CHECK(cudaMemcpy(&info, infoOnDevice, sizeof(int), cudaMemcpyDeviceToHost));

  if(info)
    //
    std::cerr << funame << "dsyevd: failed: see the log file\n";
  
  if(info < 0) {
    //
    IO::log << IO::log_offset << funame << "dsyevd: " << -info << "-th argument has illegal value" << std::endl;

    throw Error::Logic();
  }
  else if(info > 0) {
    //
    IO::log << IO::log_offset << funame << "dsyevd: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero" << std::endl;

    throw Error::Lapack();
  }
  
  CUDA_CHECK(cudaMemcpy(res,  eigenValuesOnDevice, n * sizeof(double), cudaMemcpyDeviceToHost));
  
  if(evec)
    //
    CUDA_CHECK(cudaMemcpy(*evec, matrixOnDevice, nn * sizeof(double), cudaMemcpyDeviceToHost));

  // free memory on device
  //
  CUDA_CHECK(cudaFree(matrixOnDevice));
  CUDA_CHECK(cudaFree(eigenValuesOnDevice));
  CUDA_CHECK(cudaFree(workSpaceOnDevice));
  CUDA_CHECK(cudaFree(infoOnDevice));

  // ... and on host
  //
  if(workSizeOnHost)
    //
    free(workSpaceOnHost);

  return res;
}
