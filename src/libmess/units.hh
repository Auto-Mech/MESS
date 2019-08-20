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

#ifndef UNITS_HH
#define UNITS_HH

#include <string>
#include "error.hh"

/************************** Physical Constants **************************/

struct Phys_const
{
  // conversion factors
  static const double  amu; // mass,      amu

  static const double  kelv; // energy,    Kelvin
  static const double  kcal; // energy,    kcal/mol
  static const double kjoul; // energy,    kJ/mol
  static const double  incm; // energy,    1/cm
  static const double    ev; // energy,    eV

  static const double herz; // frequency, Herz

  static const double  angstrom; // distance,  Angstrom
  static const double  cm;       // distance,  cm

  static const double avogadro; // avagadro number
  static const double light_speed; // speed of light

  static const double  bar; // pressure,  bar
  static const double  atm; // pressure,  atmosphere
  static const double  tor; // pressure,  torr

  static double bohr2ang (double x) { return x / angstrom; }
  static double ang2bohr (double x) { return x * angstrom; }

  static double amu2au (double x) { return x*amu; }
  static double au2amu (double x) { return x/amu; }

  static double kelv2hart (double x) { return x*kelv; }
  static double hart2kelv (double x) { return x/kelv; }

  static double kcal2hart (double x) { return x*kcal; }
  static double hart2kcal (double x) { return x/kcal; }

  static double str2fac (const std::string&) ;
};

#endif
