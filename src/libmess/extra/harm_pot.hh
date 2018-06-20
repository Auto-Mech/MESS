/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
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
