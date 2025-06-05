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

#include "random.hh"
#include "linpack.hh"
#include "array.hh"
#include "io.hh"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <omp.h>

namespace Random {
  //
  Array<unsigned short> buffer;
}

/****************************************************************
 ******************* Random Number Generators *******************
 ****************************************************************/

void Random::init ()
{
  const char funame [] = "Random::init: ";

  int itemp = std::time(0);

#ifdef OPENMP

#pragma omp parallel
  {
#pragma omp master
    //
    buffer.resize(3 * omp_get_num_threads());

#pragma omp barrier

    // may be initialization for an individual thread should be more sofisticated
    //
    (int&)buffer[3 * omp_get_thread_num() + 1] = itemp + omp_get_thread_num();

    buffer[3 * omp_get_thread_num()] = 0x330E;
  }
  IO::log << IO::log_offset << funame << buffer.size() / 3 << " independent channels have been initialized\n\n";

#else

  srand48(std::time(0));

#endif
}

void Random::init (int i)
{
  const char funame [] = "Random::init: ";

#ifdef OPENMP

#pragma omp parallel
  {
#pragma omp master
    //
    buffer.resize(3 * omp_get_num_threads());

#pragma omp barrier
    //
    // may be initialization for an individual thread should be more sofisticated
    //
    (int&)buffer[3 * omp_get_thread_num() + 1] = i + omp_get_thread_num();

    buffer[3 * omp_get_thread_num()] = 0x330E;
  }
  IO::log << IO::log_offset << funame << buffer.size() / 3 << " independent channels have been initialized\n\n";

#else

  srand48(i);

#endif
}

double Random::flat ()
{
  const char funame [] = "Random::flat()";

#ifdef OPENMP

  return erand48(3 * omp_get_thread_num() + (unsigned short*)buffer);
  
#else

  return  drand48 ();

#endif

}

// normal distribution with unity standard deviation
//
double Random::norm ()
{
  double dtemp, r, a;

  static int saved = 0;

  static double saved_value;
 
  if(saved) {
    //
    saved = 0;

    return saved_value;
  }

  dtemp = Random::flat();

  if(dtemp <= 0.)
    //
    return 100.;

  if(dtemp >= 1.)
    //
    return 0.;
  
  r = std::sqrt(-2. * std::log(dtemp));

  a = 2. * M_PI * Random::flat();

  saved_value = r * std::sin(a);

  saved = 1;

  return r * std::cos(a);
}

// exponential RNG
//
double Random::exp ()
{
  return -std::log(Random::flat());
}

// randomly orients vector _vec_ of the dimension _dim_ of unity length
//
void Random::orient (double* vec, int dim)
{
  const char funame [] = "Random::orient: ";

  if(dim < 2) {
    //
    std::cerr << funame << "out of range: " << dim << "\n";

    throw Error::Range();
  }

//#pragma omp critical (omp_random_section)
  //
  for(int i = 0; i < dim; ++i)
    //
    vec[i] = norm();

  ::normalize(vec, dim);
}

// random point inside the sphere of unity radius
//
void Random::volume (double* vec, int dim)
{ 
  const char funame [] = "Random::vol: ";

  if(dim < 2) {
    //
    std::cerr << funame << "out of range: " << dim << "\n";

    throw Error::Range();
  }

  orient(vec, dim);

  const double r = std::pow(flat(), 1. / (double)dim);

  for(int i = 0; i < dim; ++i)
    //
    vec[i] *= r;
}

// random point inside the spherical layer 
//
void Random::spherical_layer (double* vec, int dim, double rmin, double rmax) 
{ 
  const char funame [] = "Random::spherical_layer: ";

  if(rmin <= 0. || rmax <= 0. || rmin >= rmax || dim < 2 || vec == 0) {
    //
    std::cerr << funame << "out of range\n";

    throw Error::Range();
  }

  orient(vec, dim);

  rmin = std::pow(rmin, (double)dim);

  rmax = std::pow(rmax, (double)dim);
  
  double r = std::pow((rmax - rmin) * flat() + rmin, 1. / (double)dim);

  for(int i = 0; i < dim; ++i)
    //
    vec[i] *= r;
}

