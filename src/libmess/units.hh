/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
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

  static double str2fac (const std::string&) throw(Error::General);
};

#endif
