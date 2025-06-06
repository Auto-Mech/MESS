/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2017, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>

#include "io.hh"
#include "slatec.hh"

double x_0, dx =.000001;

double pot (double x)
{
  return x * x * (0.5 - x / x_0 / 3.);
}

double action (double ener)
{
  double dtemp;

  if(ener <= 0.)
    //
    return 2. * M_PI * ener;

  double x = 0.;

  double res = 0.;

  while(x < x_0) {
    //
    dtemp = 2. * (ener - pot(x));

    if(dtemp <= 0.)
      //
      break;

    res += std::sqrt(dtemp);

    x += dx;
  }

  x = -dx;

  while((dtemp = 2. * (ener - pot(x))) > 0.) {
    //
    res += std::sqrt(dtemp);

    x -= dx;
  }

  return 2. * res * dx;
  
}
  
int main (int argc, char* argv [])
{
  const char funame [] = "tunnel: ";

  if (argc != 2) {
    //
    std::cout << "usage: tunnel V^#\n";
    
    return 0;
  }

  int                 itemp;
  double              dtemp;
  std::string         stemp;

  IO::out.open("tunnel.out");

  
  const double ener_max = (double)IO::String(argv[1]);

  x_0 = std::sqrt(ener_max * 6.);

  int emax = 150;

  double ener_step = 1.1;

  Array<double> xtable(emax + 1);
  
  Array<double> ytable(emax + 1);

  xtable[0] = 0.;

  ytable[0] = 0.;
  
  double ener = ener_max;

  for(int e = 1; e < emax; ++e) {
    //
    ener /= ener_step;
    
    dtemp = ener_max - ener;
    
    xtable[e] = dtemp;

    ytable[e] = action(dtemp);
  }

  xtable[emax] = ener_max;

  ytable[emax] = action(ener_max);
  
  // energy-action
  //
  Slatec::Spline ener_act_spl(xtable, ytable, emax + 1);

  ener = 0.;

  for(int e = 0; e < emax; ++e) {
    //
    ytable[e] = ener_act_spl(xtable[e], 1);
  }

  // beta-energy
  //
  Slatec::Spline beta_ener_spl(ytable, xtable, emax);

  int bmax = 100;

  double beta = ytable[1];

  double bstep = (ytable[emax - 1] - beta) / bmax;

  IO::out << ener_max << "\n";
  
  for(int b = 0; b < bmax; ++b, beta += bstep) {
    //
    const double temperature = 1. / beta;

    const double eb = beta_ener_spl(beta);

    const double d2s = std::sqrt(2. * M_PI / ener_act_spl(eb, 2));

    const double pow_ref = ener_act_spl(eb) - eb * beta;
    
    // energy integration
    //
    double tfactor = 1.;

    emax = 100;
    
    ener_step = d2s / emax;

    ener = eb;
    
    while((ener += ener_step) < ener_max) {
      //
      dtemp = ener_act_spl(ener) - ener * beta - pow_ref;

      if(dtemp > 50.)
	//
	break;
      
      tfactor += std::exp(-dtemp);
    }

    ener = eb;
    
    while((ener -= ener_step) > 0.) {
      //
      dtemp = ener_act_spl(ener) - ener * beta - pow_ref;

      if(dtemp > 50.)
	//
	break;
      
      tfactor += std::exp(-dtemp);
    }

    // normalization
    //
    tfactor /= emax;

    IO::out << std::setw(15) << beta / 2. / M_PI
	    << std::setw(15) << 1. - eb / ener_max
	    << std::setw(15) << d2s / ener_max
	    << std::setw(15) << tfactor
	    << std::setw(15) <<  1. / beta / std::sqrt(2. * M_PI * (ener_max - eb))
	    << std::setw(15) << (ener_max - eb) * std::exp(beta)
	    << "\n";
  }

  return 0;
}
