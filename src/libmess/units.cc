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

#include "units.hh"
#include <iostream>


const double Phys_const::kelv = 3.166829e-6;    // energy,    Kelvin

const double Phys_const::kcal = 1.59362e-3;     // energy,    kcal/mol

const double Phys_const::kjoul = 3.8088e-4;     // energy,    kJ/mol

const double Phys_const::incm = 4.55633e-6;     // energy,    1/cm

const double Phys_const::ev   = 3.67502e-2;     // energy,    eV

const double Phys_const::angstrom = 1.88971616463;         // distance,  angstrom

const double Phys_const::cm   = 1.88971616463e+8;      // distance,  cm

//const double Phys_const::bohr = 0.529177249;    // distance,  au in angstrom

const double Phys_const::herz = 2.41888e-17;    // frequency, Herz

const double Phys_const::amu  = 1822.8885;      // mass,      atomic mass unit

const double Phys_const::bar  = 3.3988e-9;      // pressure,  bar

const double Phys_const::atm  = 3.4438e-9;      // pressure,  atmosphere

const double Phys_const::tor  = 4.5313e-12;     // pressure,  torr

const double Phys_const::avogadro = 6.02252e+23; // avogadro number

const double Phys_const::light_speed = 137.0388; // speed of light in atomic units (inverse fine-structure unit)

double Phys_const::str2fac(const std::string& unit)  
{
  const char funame [] = "Phys_const::str2fac: ";

  if(unit == "Kelvin" || unit == "kelvin" || unit == "kelv")
    //
    return kelv;
  
  if(unit == "kcal" || unit == "kcal/mol")
    //
    return kcal;
  
  if(unit == "incm" || unit == "invcm" || unit == "1/cm")
    //
    return incm;

  if(unit == "kJ/mol" || unit == "kj/mol" || unit == "kjoul/mol" || unit == "kJoul/mol")
    //
    return kjoul;

  if(unit == "eV" || unit == "ev")
    //
    return ev;
  
  if(unit == "Angstrom" || unit == "angstrom")
    //
    return angstrom;

  if(unit == "cm")
    //
    return cm;
  
  if(unit == "Bohr" || unit == "bohr" || unit == "hartree" || unit == "Hartree" || unit == "au" || unit == "nu")
    //
    return 1.;

  if(unit == "amu")
    //
    return amu;
  
  if(unit == "bar")
    //
    return bar;

  if(unit == "tor" || unit == "torr")
    //
    return tor;

  if(unit == "atm")
    //
    return atm;

  if(unit == "herz" || unit == "Herz")
    //
    return herz;
  
  if(unit == "%" || unit == "percent" || unit == "Percent")
    //
    return .01;

  std::cerr << funame << "unknown unit: " << unit << "\n";
  
  throw Error::Form();
}
