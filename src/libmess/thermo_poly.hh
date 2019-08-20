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

#ifndef THERMO_POLY_HH
#define THERMO_POLY_HH

class ThermoPoly {
  //
  typedef std::vector<double>                     _pos_t;

  typedef std::map<std::multiset<int>, double> _potex_t ;

  _potex_t _global_potex;

  std::vector<double> _frequency;

  void _read_potex (std::ifstream&, _potex_t&);

  _potex_t _local_potex (const _pos_t&, int) const;

  double   _pot (const _pos_t&) const;

public:

  void init (std::istream&);

  std::map<int, double> weight (double temperature) const;
 
};
