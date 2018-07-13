

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
