#include "limits.hh"
#include "error.hh"

#include <iostream>

namespace Limits {
  //
  double _exp_pow_max = 200.;

  double exp_pow_max () { return _exp_pow_max; }
}

void Limits::set_exp_pow_max (double v) {
  //
  const char funame [] = "Limits::set_exp_pow_max: ";

  if(v <= 0.) {
    //
    std::cerr << funame << "out of range: " << v << "\n";

    throw Error::Range();
  }

  _exp_pow_max = v;
}
