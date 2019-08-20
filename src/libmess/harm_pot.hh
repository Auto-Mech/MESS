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

#ifndef HARM_POT_HH
#define HARM_POT_HH

#include "configuration.hh"

struct EnergyConverter {
  double scale;
  EnergyConverter () : scale(-1.) {}

  double operator() (double, int =0) const;
};

double harm_fit (double, const Configuration::State&);

extern "C" void   harm_init_ (const char*);
extern "C" void no2_ch3_pot_ (const double&, const double&, double&);

#endif
