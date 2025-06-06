#include "mpack_dd.hh"

int main () {

  Mpack_dd::SymmetricMatrix m(10);

  m = 0.;

  Mpack_dd::SymmetricMatrix n = m.copy();

  Mpack_dd::Vector v = m.eigenvalues();

  return 0;
}
