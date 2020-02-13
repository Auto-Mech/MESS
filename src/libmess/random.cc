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
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

/****************************************************************
 ******************* Random Number Generators *******************
 ****************************************************************/

void Random::init ()
{
  srand48(std::time(0));
}

void Random::init (int i)
{
  srand48(i);
}

double Random::flat ()
{
    return drand48 ();
}

// normal distribution RNG
//
double Random::norm ()
{
  double dtemp;

  dtemp = Random::flat();

  if(dtemp == 0.)
    //
    return 100.;

  if(dtemp == 1.)
    //
    return 0.;
  
  return std::sqrt(-2. * std::log(dtemp)) * std::cos(2. * M_PI * Random::flat());
}

double Random::exp ()
{// RNG for the distribution exp(-u**2/2)u*du
  return std::sqrt(-2.0 * std::log(Random::flat()));
}

void Random::orient (double* vec, int dim)
{ // randomly orients vector _vec_ of the dimension _dim_ of unity length
  double norm, temp;
  do {
    norm = 0.0;
    for (int i = 0; i < dim; ++i)
      {
	temp = 2.0 * Random::flat () - 1.0;
	norm += temp * temp;
	vec[i] = temp;
      }
  } while (norm >= 1.0 || norm < 0.0001);

  norm = std::sqrt(norm);
  for (int i = 0; i < dim; ++i)
    vec[i] /= norm;
}

double Random::vol (double* vec, int dim)
{ // random point inside the sphere of unity radius
  double norm, temp;
  do {
    norm = 0.0;
    for (int i = 0; i < dim; ++i) {
      temp = 2.0 * Random::flat () - 1.0;
      norm += temp * temp;
      vec[i] = temp;
    }
  } while (norm > 1.0);

  return norm;
}

// random point inside the spherical layer 
void Random::spherical_layer (double* vec, int dim, double rmin, double rmax) 
{ 
  const char funame [] = "Random::spherical_layer: ";

  if(rmin <= 0. || rmax <= 0. || rmin >= rmax || dim <= 0 || vec == 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
    
  double norm, dtemp;
  
  const double rmin2 = rmin * rmin;
  const double rmax2 = rmax * rmax;

  double* const end = vec + dim;

  do {
    norm = 0.;
    for (double* it = vec; it != end; ++it) {
      dtemp = (2. * Random::flat () - 1.) * rmax;
      *it = dtemp;
      norm += dtemp * dtemp;
    }
  } while (norm < rmin2 || norm > rmax2);
}

