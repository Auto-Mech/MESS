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

double Phys_const::str2fac(const std::string& unit) throw(Error::General) 
{
  const char funame [] = "Phys_const::str2fac: ";

  if(unit == "Kelvin" || unit == "kelvin" || unit == "kelv")
    return Phys_const::kelv;
  else if(unit == "kcal" || unit == "kcal/mol")
    return Phys_const::kcal;
  else if(unit == "Angstrom" || unit == "angstrom")
    return Phys_const::angstrom;
  else if(unit == "incm" || unit == "invcm")
    return incm;
  else if(unit == "amu")
    return amu;
  else if(unit == "Bohr" || unit == "bohr" || 
	  unit == "hartree" || unit == "Hartree" || 
	  unit == "au" || unit == "nu") 
    return 1.;
  else if(unit == "%" || unit == "percent" || unit == "Percent")
    return .01;

  std::cerr << funame << "unknown unit: " << unit << "\n";
  throw Error::Form();
}
