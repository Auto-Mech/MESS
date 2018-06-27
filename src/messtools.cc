#include <iostream>
#include <vector>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "libmess/atom.hh"
#include "libmess/units.hh"
#include "libmess/model.hh"

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


std::vector<Atom> _atoms(const std::vector<std::string>& ls,
                         const std::vector<std::vector<double>>& xyzs,
                         const std::vector<int>& isos = std::vector<int>()) {
    unsigned int n = ls.size();

    if (xyzs.size() != n || (isos.size() !=n && isos.size() != 0))
        throw std::invalid_argument("Molecule arguments do not match");

    std::vector<Atom> atoms;
    for (unsigned int i = 0; i < n; ++i) {
        Atom a = _atom(ls[i], xyzs[i]);
        atoms.push_back(a);
    }
    return atoms;
}


// Functions for python binding
double mass(const std::string& l, int iso = 0) {
    AtomBase a = _atom_base(l, iso);
    return a.mass();
}


py::array_t<double> partition_function(
        const std::vector<double>& temps,
        const std::vector<std::string>& labels,
        const std::vector<std::vector<double>>& coords,
        const std::vector<double>& freqs,
        const std::vector<std::vector<double>>& anharm_consts,
        const std::vector<double>& rot_consts,
        const std::vector<std::vector<double>>& rovib_coups,
        const std::map<std::string, double>& rot_distorts,
        double sym_fac=1.) {
    const std::vector<Atom> atoms = _atoms(labels, coords);
    double mass = 0.;
    for (auto a : atoms) {
        mass += a.mass();
    }

    std::vector<double> elevels({0.});
    std::vector<int> edegens({1});
    std::string name = "";
    double eground = 0.;
    double emax = -1.;
    Model::RRHO species(atoms, name, Model::NOSTATES, mass, eground, elevels,
            edegens, freqs, emax, anharm_consts, rot_consts, rovib_coups,
            rot_distorts, sym_fac);


    const double vol = Phys_const::cm * Phys_const::cm * Phys_const::cm;
    double m = species.mass();
    std::vector<double> qs;
    for(auto t: temps) {
        double q = std::log(
                species.weight(t) *
                std::pow(m * t / 2. / M_PI, 1.5) *
                vol);
        qs.push_back(q);
    }

    return py::array_t<double>({qs.size()}, qs.data());
}


PYBIND11_MODULE(messtools, module) {
    module.def("mass", &mass,
               "Get the mass of an isotope",
               py::arg("l"), py::arg("iso")=0);
    module.def("partition_function",
               &partition_function,
               py::arg("temps"), py::arg("labels"), py::arg("coords"),
               py::arg("freqs"), py::arg("anharm_consts"),
               py::arg("rot_consts"), py::arg("rovib_coups"),
               py::arg("rot_distorts"), py::arg("sym_fac")=1.);
}
