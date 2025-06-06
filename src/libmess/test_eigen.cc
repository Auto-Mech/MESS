/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#include "lapack.hh"
#include "d3.hh"

#include <cstdlib>
#include <ctime>

int main ()
{
  // random number initialization
  srand48(std::time(0));

  int n = 3;
  Lapack::SymmetricMatrix m(n);

  for(int i = 0; i < n; ++i)
    for(int j = i; j < n; ++j)
      m(i, j) = drand48();

  std::cout << "the initial matrix is:\n\n";
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j)
	std::cout << std::setw(10) << m(i, j);
      std::cout << "\n";
    }
  std::cout << "\n";

  Lapack::Matrix evec(n, n);
  Lapack::Vector eval = m.eigenvalues(&evec);

  std::cout << "eigenvalues: ";
  for(int i = 0; i < n; ++i)
      std::cout << std::setw(11) << eval[i];
  std::cout << "\n\n";


  std::cout << "Transformation determinant by LU-factorization   = " 
	    << Lapack::LU(evec).det() << "\n"
	    << "Transformation determinant through vector product = "
	    <<  D3::volume(D3::Vector(&evec(0, 0)), D3::Vector(&evec(0, 1)),D3::Vector(&evec(0, 2)))
	    << "\n";

  return 0;
}
