#include "offload.hh"
#include <mkl.h>
#include <mkl_omp_offload.h>

Lapack::Vector Offload::Omp::eigenvalues (Lapack::SymmetricMatrix m, Lapack::Matrix* evec)
{
  const char funame [] = "Offload::Omp::eigenvalues: ";

  typedef MKL_INT int_t;

  if(!m.isinit()) {
    //
    std::cerr << funame << "not initialized\n";
    
    throw Error::Init();
  }

  Lapack::Matrix a;
  
  char jobz = 'N';
  
  if(evec) {
    //
    a = *evec;
    
    jobz = 'V';
  }
  else
    //
    a.resize(m.size());

  a = m;
  
  Lapack::Vector res(m.size());

  int_t lwork, liwork, info = 0, n = m.size();

  //work space query
  //
  lwork = -1;

  liwork = -1;

  Array<double> work(1);

  Array<int_t> iwork(1);

  double* a_pointer = a;

  double* res_pointer = res;

  int_t* info_pointer = &info;

  double*  work_pointer = work;

  int_t* iwork_pointer = iwork;
  
#pragma omp target data map(tofrom : a_pointer[0 : n * n], res_pointer[0 : n], info_pointer[0 : 1], work_pointer[0 : 1], iwork_pointer[0 : 1])
  
   {
     
#pragma omp target variant dispatch use_device_ptr(a_pointer, res_pointer, info_pointer, work_pointer, iwork_pointer)
     
     dsyevd(&jobz, "U", &n , a_pointer, &n, res_pointer, work_pointer, &lwork, iwork_pointer, &liwork, info_pointer);
   }
   
  // matrix diagonalization
  //
  lwork = (int_t)work[0];

  work.resize(lwork);

  liwork = iwork[0];
  
  iwork.resize(liwork);
  
  work_pointer = work;

  iwork_pointer = iwork;
  
#pragma omp target data map(tofrom : a_pointer[0 : n * n], res_pointer[0 : n], info_pointer[0 : 1], work_pointer[0 : lwork], iwork_pointer[0 : liwork])
  
   {
     
#pragma omp target variant dispatch use_device_ptr(a_pointer, res_pointer, info_pointer, work_pointer, iwork_pointer)
     
     dsyevd(&jobz, "U", &n, a, &n, res_pointer, work_pointer, &lwork, iwork_pointer, &liwork, info_pointer);
   }
   
  
  if(info < 0) {
    //
    std::cerr << funame << "dsyevd: " << -info << "-th argument has illegal value\n";

    throw Error::Logic();
  }
  else if(info > 0) {
    //
    std::cerr << funame << "dsyevd: " << info 
       //
              << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero\n";

    throw Error::Lapack();
  }

  return res;
}
