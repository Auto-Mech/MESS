#include "mpack.hh"
#include "io.hh"

#include <mpack/mpack_config.h>
#include <qd/dd_real.h>
#include <mpack/dd_complex.h>
#include <mpack/mlapack_dd.h>
#include <qd/qd_real.h>
#include <mpack/qd_complex.h>
#include <mpack/mlapack_qd.h>

Lapack::Vector Mpack::dd_eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec) throw(Error::General)
{
  const char funame [] = "Mpack::dd_eigenvalues: ";

  if(!mat.isinit()) {
    //
    ErrOut err_out;

    err_out << funame << "not initialized";
  }

  const char* jobz = "N";
  //
  if(evec)
    //
    jobz = "V";

  mpackint lwork, info;

  Array<dd_real> a(mat.size() * mat.size());

  Array<dd_real> w(mat.size());

  for(int i = 0; i < mat.size(); ++i)
    //
    for(int j = i; j < mat.size(); ++j)
      //
      a[i + j * mat.size()] = mat(i, j);

  //work space query
  //
  lwork = -1;

  Array<dd_real> work(1);

  Rsyev(jobz, "U", mat.size(), a, mat.size(), w, work, lwork, &info);

  // matrix diagonalization
  //
  lwork = (int)work[0].x[0];

  work.resize(std::max((mpackint)1, lwork));

  Rsyev(jobz, "U", mat.size(), a, mat.size(), w, work, lwork, &info);

  if(info < 0) {
    //
    ErrOut err_out;

    err_out << funame << "Rsyev: " << -info << "-th argument has an illegal value";
  }
  else if(info > 0) {
    //
    ErrOut err_out;

    err_out << funame << "Rsyev: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero";
  }

  Lapack::Vector res(mat.size());
  
  for(int i = 0; i < mat.size(); ++i)
    //
    res[i] = w[i].x[0];

  if(evec) {
    //
    for(int i = 0; i < mat.size(); ++i) {
      //
      for(int j = 0; j < mat.size(); ++j) {
	//
	(*evec)(i, j) = a[i + j * mat.size()].x[0];
      }
    }
  }

  return res;
}

Lapack::Vector Mpack::qd_eigenvalues (Lapack::SymmetricMatrix mat, Lapack::Matrix* evec) throw(Error::General)
{
  const char funame [] = "Mpack::qd_eigenvalues: ";

  if(!mat.isinit()) {
    //
    ErrOut err_out;

    err_out << funame << "not initialized";
  }

  const char* jobz = "N";
  //
  if(evec)
    //
    jobz = "V";

  mpackint lwork, info;

  Array<qd_real> a(mat.size() * mat.size());

  Array<qd_real> w(mat.size());

  for(int i = 0; i < mat.size(); ++i)
    //
    for(int j = i; j < mat.size(); ++j)
      //
      a[i + j * mat.size()] = mat(i, j);

  //work space query
  //
  lwork = -1;

  Array<qd_real> work(1);

  Rsyev(jobz, "U", mat.size(), a, mat.size(), w, work, lwork, &info);

  // matrix diagonalization
  //
  lwork = (int)work[0].x[0];

  work.resize(std::max((mpackint)1, lwork));

  Rsyev(jobz, "U", mat.size(), a, mat.size(), w, work, lwork, &info);

  if(info < 0) {
    //
    ErrOut err_out;

    err_out << funame << "Rsyev: " << -info << "-th argument has an illegal value";
  }
  else if(info > 0) {
    //
    ErrOut err_out;

    err_out << funame << "Rsyev: " << info << "-th off-diagonal  elements  of intermediate tridiagonal form did not converge to zero";
  }

  Lapack::Vector res(mat.size());
  
  for(int i = 0; i < mat.size(); ++i)
    //
    res[i] = w[i].x[0];

  if(evec) {
    //
    for(int i = 0; i < mat.size(); ++i) {
      //
      for(int j = 0; j < mat.size(); ++j) {
	//
	(*evec)(i, j) = a[i + j * mat.size()].x[0];
      }
    }
  }

  return res;
}

