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

#include <fstream>

#include "lr.hh"
#include "read.hh"
#include "slatec.hh"
#include "slatec.h"

#include <cmath>
#include <sstream>

enum calc_t {E_LEVEL, J_LEVEL, M_LEVEL, K_LEVEL};

class ShortDoubleOut : public std::ostringstream {
 public:
  ShortDoubleOut& operator<< (double);
  ShortDoubleOut& operator<< (const   std::string& s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (const          char* s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (const   signed char* s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (const unsigned char* s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (                 int s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (               short s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (                long s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (      unsigned   int s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (      unsigned short s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
  ShortDoubleOut& operator<< (      unsigned  long s) { *static_cast<std::ostringstream*>(this) << s; return *this; }
};

ShortDoubleOut& ShortDoubleOut::operator<< (double d)
{
  std::ostringstream& oss = *this; 
  if( d >= 10. || d <= -10.)
    oss << (int)d; 
  else
    oss << d;
  return *this;
}

class FreeEnergy: public Math::GradientSearch { // S - E/T, E = exp(x)

  const Slatec::Spline& _entropy;
public:
  double beta;  // T^-1

  double operator() (double, int) const ;

  FreeEnergy (const Slatec::Spline& e, double temperature, double tol) 
    : Math::GradientSearch(tol), _entropy(e), beta(1./temperature) {}
};

double FreeEnergy::operator() (double x, int n) const 
{
  return _entropy(x, n) - beta * std::exp(x);
}

/*****************************************************************************************
 ************************************** MAIN *********************************************
 *****************************************************************************************/

int main (int argc, char* argv [])
{
  const char funame [] = "lr: ";

  static const double exp_max = 50.;
  static const double one_eps = 1. + 1.e-10;
 
  int    itemp;
  bool   btemp;
  double dtemp;

  if (argc < 2) {
    std::cout << funame  << "usage: lr input_file\n";
    return 1;
  }

  std::map<std::string, Read> input;
  std::map<std::string, Read>::iterator idit;
 
  double ener_max, ener_min, temp_max, temp_min;
  int coarse_ener_size, medium_ener_size, fine_ener_size, coarse_angle_size, medium_angle_size, temp_size;
  int j_size, m_size, k_size; // number of points for j, m, and k ranges
  double dist_guess; // initial guess
  std::string calc_key;

  input ["Structure"              ] = Read(Structure::init);
  input ["Main"                   ] = Read(LongRange::init);
  input ["Potential"              ] = Read(LongRange::pot);
  input ["EnergyMaximum[kcal/mol]"] = Read(ener_max);
  input ["EnergyMinimum[kcal/mol]"] = Read(ener_min);
  input ["CoarseEnergySize"       ] = Read(coarse_ener_size);        // energy grid for correction factor 
  input ["MediumEnergySize"       ] = Read(medium_ener_size);        // energy grid for e-resolved number of states
  input ["FinestEnergySize"       ] = Read(fine_ener_size);          // energy grid for rate constant
  input ["CoarseAngularSize"      ] = Read(coarse_angle_size, 10);   // angular grid for correction factor
  input ["MediumAngularSize"      ] = Read(medium_angle_size, 20);   // angular grid for e-resolved number of states
  input ["TemperatureMaximum[K]"  ] = Read(temp_max);
  input ["TemperatureMinimum[K]"  ] = Read(temp_min);
  input ["TemperatureSize"        ] = Read(temp_size);
  input ["J-Size"                 ] = Read(j_size, 40);               // number of grid points for J integration
  input ["M-Size"                 ] = Read(m_size, 10);               // number of grid points for M integration
  input ["K-Size"                 ] = Read(k_size, 10);               // number of grid points for K integration
  input ["DistanceGuess[au]"      ] = Read(dist_guess);
  input ["CalculationLevel"       ] = Read(calc_key, "E");

  // actual input
  std::ifstream from(argv[1]);
  if(!from) {
    std::cerr << funame << "input file " << argv[1] << " is not found\n";
    return 1;
  }

  std::string key, comment;
  while(from >> key) {
    idit = input.find(key);
    if (idit == input.end()) {
      std::cerr << funame << "WARNING: did not find the key: " << key
		<< ", taking it as a comment\n";
      getline(from, comment);
    }
    else
      from >> idit->second;
  }
  from.close();
  from.clear();

  // check if all parameters were initialized
  btemp = true;
  for(idit = input.begin(); idit != input.end(); ++idit)
    if(!idit->second.is_init()) {
      std::cerr << funame << idit->first << " is not initialized\n";
      btemp = false;
    }
  if(!btemp)
    return 1;

  // checking for default values
  std::cout << funame << "Default parameters:\n" << std::left;
  for(idit = input.begin(); idit != input.end(); ++idit)
    if(idit->second.is_default())
      std::cout << "   " << std::setw(20) << idit->first << " = " << idit->second << "\n";
  std::cout << "\n" << std::right;

  calc_t clevel;
  if(     calc_key == "E")
    clevel = E_LEVEL;
  else if(calc_key == "J")
    clevel = J_LEVEL;
  else if(calc_key == "M") {
    clevel = M_LEVEL;
    if(!LongRange::symmetric_fragment_size()) {
      std::cerr << funame << "sorry, there are no symmetric fragments\n";
      throw Error::Input();
    }
  }
  else if(calc_key == "K")
    clevel = K_LEVEL;
  else {
    std::cerr << funame << "unknown calculation level: " << calc_key 
	      << ", available levels: E, J, M, or K\n";
    throw Error::Input();
  }

  ener_max *= Phys_const::kcal;
  ener_min *= Phys_const::kcal;
  temp_max *= Phys_const::kelv;
  temp_min *= Phys_const::kelv;

  double ener, ener_step;

  Array<double> x_data, y_data;
  Slatec::Spline corr;

  // coarse grid calculation of the optimization correction factor

  LongRange::ENumberOpt opt_e(dist_guess);

  if(clevel > E_LEVEL) {//correction factor to E-resolved number of states

    LongRange::set_angular_grid(coarse_angle_size);
    x_data.resize(coarse_ener_size);
    y_data.resize(coarse_ener_size);
    
    std::cout << std::setw(20) << "Energy, kcal/mol"
	      << std::setw(20) << "Correction Factor"
	      << std::setw(20) << "TS Distance, bohr"
	      << std::setw(10)  << "J_max/rms";
    if(clevel > J_LEVEL)
      for(int frag = 0; frag < LongRange::symmetric_fragment_size(); ++frag) 
	std::cout << std::setw(10) << "M_max/rms";
    if(clevel > M_LEVEL)
      std::cout << std::setw(10) << "K_max/rms";
    std::cout << std::endl;

    // j-resolved level
    LongRange::JNumberFix fix_j (dist_guess);
    LongRange::JIntegral  fix_ej(&fix_j, ener_max);
    
    LongRange::JNumberOpt opt_j (dist_guess);
    LongRange::JIntegral  opt_ej(&opt_j, ener_max);

    // m-resolved number of states correction factor
    LongRange::MNumberFix fix_m  (dist_guess);
    LongRange::MIntegral  fix_jm (&fix_m);
    LongRange::JIntegral  fix_ejm(&fix_jm, ener_max);

    LongRange::MNumberOpt opt_m  (dist_guess);
    LongRange::MIntegral  opt_jm (&opt_m);
    LongRange::JIntegral  opt_ejm(&opt_jm, ener_max);

    // k-resolved number of states correction factor
    LongRange::KNumberFix fix_k   (dist_guess);
    LongRange::KIntegral  fix_mk  (&fix_k);
    LongRange::MIntegral  fix_jmk (&fix_mk);
    LongRange::JIntegral  fix_ejmk(&fix_jmk, ener_max);

    LongRange::KNumberOpt opt_k   (dist_guess);
    LongRange::KIntegral  opt_mk  (&opt_k);
    LongRange::MIntegral  opt_jmk (&opt_mk);
    LongRange::JIntegral  opt_ejmk(&opt_jmk, ener_max);

    ener      = ener_max;
    ener_step = std::pow(ener_max / ener_min, 1. / double(coarse_ener_size - 1));

    ShortDoubleOut oss;
    oss.precision(1);

    for(int i = 0; i < coarse_ener_size; ++i, ener /= ener_step) {// coarse grid cycle

      int ii = coarse_ener_size - i - 1;
      x_data[ii] = std::log(ener);

      // optimize distance
      opt_e.set_ener(ener);
      opt_e();
      double dist = opt_e.distance();

      // number of states correction factor
      double opt_val, fix_val;
      switch(clevel) {
      case J_LEVEL:

	fix_j.set_dist(dist);
	fix_ej.set_ener(ener);
    
	opt_j.set_dist(dist);
	opt_ej.set_ener(ener);

	opt_val = opt_ej();
	fix_val = fix_ej();

	y_data[ii] = opt_val / fix_val;

	// output
	std::cout << std::setw(20) << ener / Phys_const::kcal
		  << std::setw(20) << y_data[ii]
		  << std::setw(20) << dist;

	oss.str("");
	oss << opt_ej.arg_max << "/" << std::sqrt(opt_ej.arg_var / opt_val);
	std::cout << std::setw(10)  << oss.str();

	// set up new integration steps
	if(opt_ej.arg_max > 0.)
	  LongRange::JIntegral::step = opt_ej.arg_max / (double)j_size;
	opt_ej.arg_max = -1.;
	opt_ej.arg_var = 0.;
	fix_ej.arg_max = -1.;
	fix_ej.arg_var = 0.;
	
	std::cout << std::endl;
	break;

      case M_LEVEL:

	fix_m.set_dist(dist);
	fix_ejm.set_ener(ener);

	opt_m.set_dist(dist);
	opt_ejm.set_ener(ener);

	opt_val = opt_ejm();
	fix_val = fix_ejm();

	y_data[ii] = opt_val / fix_val;

	// output
	std::cout << std::setw(20) << ener / Phys_const::kcal
		  << std::setw(20) << y_data[ii]
		  << std::setw(20) << dist;
	
	oss.str("");
	oss <<  opt_ejm.arg_max << "/" << std::sqrt(opt_ejm.arg_var / opt_val);
	std::cout << std::setw(10)  << oss.str();

	// set up new integration steps
	if(opt_ejm.arg_max > 0.)
	  LongRange::JIntegral::step = opt_ejm.arg_max / (double)j_size;
	opt_ejm.arg_max = -1.;
	opt_ejm.arg_var =  0.;
	fix_ejm.arg_max = -1.;
	fix_ejm.arg_var =  0.;
	
	for(int frag = 0; frag < LongRange::symmetric_fragment_size(); ++frag) {
	  oss.str("");
	  oss << opt_jm.arg_max[frag] << "/" << std::sqrt(opt_jm.arg_var[frag] / opt_val);
	  std::cout << std::setw(10)  << oss.str();

	  // set up new integration steps
	  if(opt_jm.arg_max[frag] > 0.)
	    LongRange::MIntegral::step[frag] = opt_jm.arg_max[frag] / (double)m_size;
	  opt_jm.arg_max[frag] = -1.;
	  opt_jm.arg_var[frag] =  0.;
	  fix_jm.arg_max[frag] = -1.;
	  fix_jm.arg_var[frag] =  0.;
	}

	std::cout << std::endl;
	break;

      case K_LEVEL:
      
	fix_k.set_dist(dist);
	fix_ejmk.set_ener(ener);

	opt_k.set_dist(dist);
	opt_ejmk.set_ener(ener);

	opt_val = opt_ejmk();
	fix_val = fix_ejmk();

	y_data[ii] = opt_val / fix_val;

	// output
	std::cout << std::setw(20) << ener / Phys_const::kcal
		  << std::setw(20) << y_data[ii]
		  << std::setw(20) << dist;

	oss.str("");
	oss << opt_ejmk.arg_max << "/" << std::sqrt(opt_ejmk.arg_var / opt_val);
	std::cout << std::setw(10)  << oss.str();

	// set up new integration steps
	if(opt_ejmk.arg_max > 0.)
	  LongRange::JIntegral::step = opt_ejmk.arg_max / (double)j_size;
	opt_ejmk.arg_max = -1.;
	opt_ejmk.arg_var =  0.;
	fix_ejmk.arg_max = -1.;
	fix_ejmk.arg_var =  0.;

	for(int frag = 0; frag < LongRange::symmetric_fragment_size(); ++frag) {
	  oss.str("");
	  oss <<  opt_jmk.arg_max[frag] << "/" << std::sqrt(opt_jmk.arg_var[frag] / opt_val);
	  std::cout << std::setw(10)  << oss.str();

	  // set up new integration steps
	  if(opt_jmk.arg_max[frag] > 0.)
	    LongRange::MIntegral::step[frag] = opt_jmk.arg_max[frag] / (double)m_size;
	  opt_jmk.arg_max[frag] = -1.;
	  opt_jmk.arg_var[frag] =  0.;
	  fix_jmk.arg_max[frag] = -1.;
	  fix_jmk.arg_var[frag] =  0.;
	}

	oss.str("");
	oss <<  opt_mk.arg_max << "/" <<  std::sqrt(opt_mk.arg_var / opt_val);
	std::cout << std::setw(10)  << oss.str();

	// set up new integration steps
	if(opt_mk.arg_max > 0.)
	  LongRange::KIntegral::step = opt_mk.arg_max / (double)k_size;
	opt_mk.arg_max = -1.;
	opt_mk.arg_var =  0.;
	fix_mk.arg_max = -1.;
	fix_mk.arg_var =  0.;

	std::cout << std::endl;
	break;
      }
    } // coarse grid cycle
    std::cout << std::endl;

    corr.init(x_data, y_data, coarse_ener_size);
  }// !e-level

  /******************************* number of states calculation *********************************/
  LongRange::set_angular_grid(medium_angle_size);

  x_data.resize(medium_ener_size);
  y_data.resize(medium_ener_size);
  
  opt_e.set_dist(dist_guess);

  std::cout << std::setw(20) << "Energy, kcal/mol"
	    << std::setw(20) << "States Number"
	    << std::setw(20) << "TS Distance, bohr"
	      << std::endl;

  ener_max /= one_eps;
  ener_min *= one_eps;

  int ee;
  ener      = ener_max;
  ener_step = std::pow(ener_max / ener_min, 1. / double(medium_ener_size - 1));

  for(int e = 0; e < medium_ener_size; ++e, ener /= ener_step) {

    ee = medium_ener_size - e - 1;
    x_data[ee] = std::log(ener);

    opt_e.set_ener(ener);
    y_data[ee] = opt_e();
    if(clevel != E_LEVEL)
      y_data[ee] *= corr(x_data[ee], 0);

    std::cout << std::setw(20) << ener / Phys_const::kcal
	      << std::setw(20) << y_data[ee]
	      << std::setw(20) << opt_e.distance()
	      << std::endl;

    y_data[ee] = std::log(y_data[ee]);
  }
  std::cout << std::endl;

  Slatec::Spline ne_log(x_data, y_data, medium_ener_size);


  /******************************* rate calculation *********************************/

  Array<double> ener_data(fine_ener_size + 1);
  ener_data[0] = 0.;

  ener_max /= one_eps;
  ener_min *= one_eps;

  ener      = ener_min;
  ener_step = std::pow(ener_max / ener_min, 1. / double(fine_ener_size - 1));

  for(int e = 0; e < fine_ener_size; ++e, ener *= ener_step)
    ener_data[e + 1] = ener;

  Array<double> int_data(fine_ener_size + 1);
  int_data[0] = 0.;

  //normalization factor
  double nfac = std::sqrt(2. * M_PI) / Structure::mass() / std::sqrt(Structure::mass());

  for(int frag = 0; frag < 2; ++frag) {
    switch(LongRange::symmetry_type(frag)) {
    case LongRange::SPHERICAL:

      break;

    case LongRange::LINEAR:

      nfac /= 2. * Structure::fragment(frag).imom(1);
      break;
	
    default: // nonlinear
	
      for(int i = 0; i < 3; ++i)
	nfac /= Structure::fragment(frag).imom_sqrt(i);
      nfac /= 2. * std::sqrt(2. * M_PI);
      break;
    }
  }
  
  std::cout << std::setw(20) << "Temperature, K"
	    << std::setw(20) << "k, 10^-11 cm^3/sec"
	    << std::setw(20) << "E_max, kcal/mol"
	    << std::setw(23) << "high energy cutoff, %"
	    << std::endl;

  double xmax = std::log(ener_max);
  double xguess = xmax;
  FreeEnergy free_ener(ne_log, 1./temp_max, 1.e-5);

  double rate;
  double temperature = temp_max;
  double temp_step   = std::pow(temp_max / temp_min, 1. / double(temp_size - 1));

  for(int t = 0; t < temp_size; ++t, temperature /= temp_step) {// temperature cycle
    ee = 1;
    for(int e = 0; e < fine_ener_size; ++e, ++ee) {
      ener = ener_data[ee];

      dtemp = ener / temperature;
      if(dtemp > exp_max)
	break;

      int_data[ee] = std::exp(-dtemp + ne_log(std::log(ener)));  
    }

    if(ee > 1) {
      davint_(ener_data, int_data, ee, 0., ener_data[ee - 1], rate, itemp);
      if(itemp != 1) {
	std::cerr << "main: davint integration error\n";
	throw Error::Run();
      }
    }
    else
      rate = 0.;

    // cutoff error calculation
    double cut_err = -1;
    free_ener.beta = 1. / temperature;
    try {
      xguess = free_ener.find(xguess);
      dtemp = free_ener(xguess, 0) - free_ener(xmax, 0);
      if(dtemp > 0.) {
	cut_err = erfc(std::sqrt(dtemp)) / 2.;
      }
      else
	std::cerr << funame << "free energy is not at minimum\n";
    }
    catch(Math::Exception) {
      std::cerr << funame << "gradient search for minimal free energy failed\n";
    }
    
    // normalization and output
    rate *= nfac / std::pow(temperature, (double)LongRange::eff_dof() / 2.);

    std::cout << std::setw(20) << temperature / Phys_const::kelv
	      << std::setw(20) << rate * 612. // 10^-11 cm^3/sec
	      << std::setw(20) << std::exp(xguess) / Phys_const::kcal;
    if(cut_err < 1.e-5)
      std::cout << std::setw(23) << "0";
    else
      std::cout << std::setw(23) << cut_err * 100;
    std::cout << std::endl;

    } // temperature cycle

    return 0;
  
  }
