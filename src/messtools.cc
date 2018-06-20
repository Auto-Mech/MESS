#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "libmess/atom.hh"
#include "libmess/units.hh"

namespace py = pybind11;


// Helpers
AtomBase _atom_base(const std::string& l, int iso = 0) {
    AtomBase a;
    if (iso == 0)
        a = AtomBase(l);
    else
        a = AtomBase(l, iso);
    return a;
}


Atom _atom(const std::string& l, const std::vector<double>& xyz, int iso = 0) {
    if (xyz.size() != 3)
        throw std::invalid_argument("Coordinate vector must have 3 elements.");

    Atom a = Atom(l);
    if (iso == 0)
        a = Atom(l);
    else
        a = Atom(l, iso);
    a[0] = xyz[0];
    a[1] = xyz[1];
    a[2] = xyz[2];
    return a;
}


// Functions for python binding
double mass(const std::string& l, int iso = 0) {
    AtomBase a = _atom_base(l, iso);
    return a.mass();
}


PYBIND11_MODULE(messtools, module) {
    module.def("mass", &mass,
               "Get the mass of an isotope",
               py::arg("l"), py::arg("iso")=0);
}
