#include "offload.hh"
#include "io.hh"
#include <sycl/sycl.hpp>

#ifdef ONE_MKL
#include <oneapi/mkl.hpp>
#define MATH mkl
#else
#include <oneapi/math.hpp>
#define MATH math
#endif


Lapack::Vector Offload::Sycl::eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec)
{  
  const char funame [] = "Offload::Sycl::eigenvalues: ";

  //using namespace oneapi::MATH;
  
  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  Lapack::Matrix a = m;

  Lapack::Vector res(m.size());

  oneapi::MATH::uplo uplo = oneapi::MATH::uplo::upper;
  
  oneapi::MATH::job  jobz = oneapi::MATH::job::N;

  if(evec) {
    //
    jobz = oneapi::MATH::job::V;
  }
  
  int64_t n = m.size();

  int64_t nn = n * n;

  // compute on device
  //
  {
    sycl::queue queue(sycl::gpu_selector_v);

    IO::log << IO::log_offset << funame << "running on: "
	    << queue.get_device().get_info<sycl::info::device::name>()
	    << std::endl;

    auto lwork = oneapi::MATH::lapack::syevd_scratchpad_size<double>(queue, jobz, uplo, n, n);

    IO::log << IO::log_offset << funame << "requested lwork = " << lwork << std::endl;
 
    auto work          = sycl::malloc_device<double>(lwork, queue);

    auto res_on_device = sycl::malloc_device<double>(n, queue);
  
    auto   a_on_device = sycl::malloc_device<double>(nn, queue);

    queue.copy((const double*)a, a_on_device, nn).wait_and_throw();

    oneapi::MATH::lapack::syevd(queue, jobz, uplo, n, a_on_device, n,
			       //
			       res_on_device, work, lwork).wait_and_throw();

    queue.copy(res_on_device, (double*)res, n);

    if(evec)
      //
      queue.copy(a_on_device, (double*)*evec, nn);

    queue.wait_and_throw();

    sycl::free(a_on_device,   queue);

    sycl::free(res_on_device, queue);

    sycl::free(work,          queue);
  }
  
  return res;
}
