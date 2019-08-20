/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2019, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef OPT_HH
#define OPT_HH

#include "coord.hh"
#include "shared.hh"
#include "harding.hh"

// molecular geometry optimization

namespace Opt {

  using namespace Coord;

  ConstSharedPointer<CartFun> new_cartesian_function (std::istream&);

  /********************* POTENTIALS ************************/
  
  // wrapper for Harding potential
  //
  class HardPot : public CartFun, private Harding {
    //
  public:
    //
    HardPot(std::istream& from) : Harding(from) {}
    
    double evaluate (const Cartesian&) const;
  };
    
  /********************* CONSTRAINS ************************/
  
  // internal coordinate: distance, angle, or dihedral
  //
  class Internal : public CartFun, private std::vector<int> {
    //
  public:
    //
    Internal(std::istream& from);
    
    double evaluate (const Cartesian&) const;
  };
    
  // constrained optimization in cartesian configurational space
  //
  class CcOpt {
    //
    // potential
    //
    ConstSharedPointer<CartFun> _pot;

    // set of constrains
    //
    std::vector<ConstSharedPointer<CartFun> > _constrain;

    int _atom_size;
    
    // energy gradient tolerance
    //
    double _grad_tol;

    // initial time step
    //
    double _time_step;

    // integration time
    //
    double _time_length;
    
  public:
    //
    CcOpt (std::istream&);

    // potential
    //
    const CartFun& pot() const { return *_pot; }
    
    // RHS of differential equation calculator
    //
    void operator () (const Cartesian&, Cartesian&, double) const;

    // find a minimum
    //
    void execute (Cartesian&) const;
  };
}
#endif
