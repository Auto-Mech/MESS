

#ifndef STRUCTURE_HH
#define STRUCTURE_HH

#include "shared.hh"
#include "atom.hh"
#include "lapack.hh"
#include "io.hh"
#include "symmetry.hh"
#include "configuration.hh"

#include <vector>
#include <iostream>

// Molecule
class Molecule : private std::vector<Atom>, public IO::Read {

  std::string _name;
  int         _type;
  int          _top;
  double      _mass;

  double _iner_mom   [3]; // inertia moments;
  double _imom_sqrt  [3];
 
  // principal axes frame
  D3::Frame _std_frame;

  // multipole moments
  int                 _charge;
  std::vector<double> _dipole;
  std::vector<double> _quadrupole;
  std::vector<double> _polarizability;

  D3::Vector _dipol_vec; // dipole vector
  D3::Matrix  _quad_mat; // quadrupole matrix
  D3::Matrix _polar_mat; // polarizability matrix

public:

  enum {MONOATOMIC, LINEAR, NONLINEAR};
  enum {ASYM, OBLATE, PROLATE, SPHERICAL};

  Molecule () {}

  // reads geometry and possibly other features from the input stream
  // and uniquely orients molecular structure to the principal axes 
  void read (std::istream&) throw(Error::General);

  // output the molecule geometry etc.
  void print (std::ostream&) const throw(Error::General);

  Molecule (std::istream& from) throw(Error::General) { read(from); }

  double mass () const { return _mass; }
  int    type () const { return _type; }
  int    top  () const { return _top;  }
  
  int size () const { return std::vector<Atom>::size(); }
  const Atom& operator[] (int i) const { return std::vector<Atom>::operator[](i); }

  static double tolerance;
  ConstSharedPointer<Symmetry::SpaceGroup> symmetry_group;
  std::set<Permutation> permutation_symmetry_group (int f) const { return ::permutation_symmetry_group(*this, tolerance, f); }

  // inertia moments
  double imom       (int) const throw(Error::General);
  double imom_sqrt  (int) const throw(Error::General);

  // dimensions of angular vector & angular velocity
  int pos_size () const throw(Error::General);
  int vel_size () const throw(Error::General);

  std::vector<int> ref_group;
  D3::Matrix      ref_orient;

  bool is_equal (std::vector<Atom>, D3::Vector&, Quaternion&) const;

  // dynamic variables derivatives
  void set_dvd (const D3::Vector& torque, const double* pos, const double* vel,
			double* pos_drv, double* vel_drv) const throw(Error::General);

  // number of rotational degrees of freedom
  int tm_dof () const throw(Error::General); 

  const std::string& name () const { return _name; }

  // from original frame to standard frame vector transformation
  //void ov2sv (const double*, double*) const;
  //void op2sp (const double*, double*) const;

  int charge () const { return _charge; }
  const std::vector<double>& dipole         () const { return _dipole;         }
  const std::vector<double>& quadrupole     () const { return _quadrupole;     }
  const std::vector<double>& polarizability () const { return _polarizability; }

  const D3::Vector&           dipole_vector () const { return _dipol_vec; }
  const D3::Matrix&       quadrupole_matrix () const { return _quad_mat; }
  const D3::Matrix&   polarizability_matrix () const { return _polar_mat; }
};

inline std::ostream& operator<< (std::ostream& to, const Molecule& m) throw(Error::General)
{
  m.print(to); 
  return to;
}

inline int Molecule::pos_size () const throw(Error::General)
{
  const char funame [] = "Molecule::pos_size: ";

  switch(_type) {
  case MONOATOMIC:
    return 0;
  case LINEAR:
    return 3;
  case NONLINEAR:
    return 4;
  default:
    std::cerr << funame << "unknown molecular type\n";
    throw Error::Logic();
  }  
}

inline int Molecule::vel_size () const throw(Error::General)
{
  const char funame [] = "Molecule::vel_size: ";

  switch(_type) {
  case MONOATOMIC:
    return 0;
  case LINEAR:
    return 3;
  case NONLINEAR:
    return 3;
  default:
    std::cerr << funame << "unknown molecular type\n";
    throw Error::Logic();
  }  
}

inline int Molecule::tm_dof () const throw(Error::General)
{
  const char funame [] = "Molecule::dof: ";

  switch(_type) {
  case MONOATOMIC:
    return 0;
  case LINEAR:
    return 2;
  case NONLINEAR:
    return 3;
  default:
    std::cerr << funame << "unknown molecular type\n";
    throw Error::Logic();
  }  
}

namespace Structure {

  const Molecule& fragment (int frag) throw(Error::General);
  int type (int frag);
  int top  (int frag);

  int size    () throw(Error::General); // total number of atoms
  double mass () throw(Error::General); // reduced mass
  double mass_sqrt () throw(Error::General); // reduced mass square root
 
  // Offset for cm-to-cm vector relative to
  // the beginning of the cooridinates data section (0)
  int orb_pos () throw(Error::General);

  // Offset for orientational coordinates of the fragments
  // relative to the beginning of the cooridinates data section (0)
  int ang_pos (int frag) throw(Error::General);

  // Offset for cm-to-cm velocity relative 
  // to the beginning of the velocities data section (pos_size)
  int orb_vel () throw(Error::General);

  // Offset for angular velicities of the fragments 
  // relative to the beginning of the velocities data (pos_size)
  int ang_vel (int frag) throw(Error::General);

  int pos_size () throw(Error::General); // size of the coordinates data section
  int pos_size (int frag);
  int vel_size () throw(Error::General); // size of the velocities  data section
  int vel_size (int frag);
  int dv_size  () throw(Error::General); // total dynamical variables dimension

  int tm_dof  () throw(Error::General); // number of degrees of freedom for transitional modes
  int tm_dof  (int frag);

  // setup
  bool isinit ();
  void init (std::istream&) throw(Error::General);
}

#endif
