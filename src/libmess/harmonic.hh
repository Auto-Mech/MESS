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
  HarmonicExpansion (const Configuration::DoubleSpaceGroup&, int, int) ;
  HarmonicExpansion (std::istream&,                          int, int) ;
  
  int size () const { return _expansion.size(); }

  double operator () (int, const Configuration::State&) const ;
    
  friend std::ostream& operator<< (std::ostream&, const HarmonicExpansion&); 
};

std::ostream& operator<< (std::ostream&, const HarmonicExpansion&); 

#endif
