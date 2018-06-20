/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef HARMONIC_HH
#define HARMONIC_HH

#include "monom.hh"
#include "configuration.hh"

/*******************************************************************************************
 ****************** HARMONIC EXPANSION OF THE ANGULAR PART OF THE POTENTIAL *****************
 *******************************************************************************************/

class HarmonicExpansion {
  Monom  _rmonom;
  Monom  _qmonom;
  std::vector<std::map<int, double> > _expansion;

public:
  HarmonicExpansion (const Configuration::DoubleSpaceGroup&, int, int) throw(Error::General);
  HarmonicExpansion (std::istream&,                          int, int) throw(Error::General);
  
  int size () const { return _expansion.size(); }

  double operator () (int, const Configuration::State&) const throw(Error::General);
    
  friend std::ostream& operator<< (std::ostream&, const HarmonicExpansion&); 
};

std::ostream& operator<< (std::ostream&, const HarmonicExpansion&); 

#endif
