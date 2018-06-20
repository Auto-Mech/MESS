/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include "structure.hh"
#include "units.hh"
#include "key.hh"

#include <sstream>
#include <cmath>
#include <cstdio>

double Molecule::tolerance = 1.e-3;

namespace Structure {

  Molecule _fragment [2];

  int    _size;
  double _mass;
  double _mass_sqrt;

  int _ang_pos [2];
  int _ang_vel [2];

  int _pos_size;
  int _vel_size;
  int  _dv_size;

  int _tm_dof;

  bool _isinit = false;

  // setup
  bool isinit () { return _isinit; }

}

const Molecule& Structure::fragment (int frag) throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _fragment[frag];
}

int Structure::type (int frag) { return fragment(frag).type(); }
int Structure::top  (int frag) { return fragment(frag).top();  }

int Structure::size    () throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _size;
}

double Structure::mass () throw(Error::General)
{ 
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _mass; 
}

double Structure::mass_sqrt () throw(Error::General)
{ 
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _mass_sqrt; 
}

int Structure::orb_pos () throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return 0;
}
int Structure::ang_pos (int frag) throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _ang_pos[frag];
}

int Structure::orb_vel () throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return 0;
}
int Structure::ang_vel (int frag) throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _ang_vel[frag]; 
}

int Structure::pos_size () throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _pos_size;
}
int Structure::vel_size () throw(Error::General)
{ 
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _vel_size;
}

int Structure::pos_size (int frag) 
{ 
  return fragment(frag).pos_size(); 
}

int Structure::vel_size (int frag) 
{ 
  return fragment(frag).vel_size(); 
}

int Structure::dv_size () throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _dv_size;
}

int Structure::tm_dof(int frag) {
  return fragment(frag).tm_dof();
}

int Structure::tm_dof () throw(Error::General)
{
  if(!isinit()) {
    std::cerr << "Structure::*(*): was not initialized\n";
    throw Error::Init();
  }

  return _tm_dof;
}

void Structure::init (std::istream& from) throw(Error::General)
{
  const char funame [] = "Structure::init: ";

  if(_isinit) {
    std::cerr << funame << "has been already initialized\n";
    throw Error::Init();
  }
  _isinit = true;


  int    itemp;
  double dtemp;

  std::string token, line, comment;

  KeyGroup StructureGroup;

  Key frag_key("Fragment" );
  Key  tol_key("Tolerance");

  int frag = 0;
  while(from >> token)
    if(token == IO::end_key())
      break;
    else if(token == frag_key) { // fragments
      if(frag > 1) {
	std::cerr << funame << "too many fragments\n";
	throw Error::Form();
      }
      _fragment[frag].read(from);
      ++frag;
    }
    else if(token == tol_key) {
      if(!(from >> Molecule::tolerance)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      if(Molecule::tolerance <= 0. || Molecule::tolerance >= 1.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
    }
    else {
      std::cerr << funame << "unknown key: " << token << "\nAvailable keys:\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }  
  
  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Form();
  }

  if(frag < 2) {
    std::cerr << funame << "fragments were not initialized\n";
    throw Error::Init();
  }
    
  _size = fragment(0).size() + fragment(1).size();
  _mass = fragment(0).mass() * fragment(1).mass() /
    (fragment(0).mass() + fragment(1).mass());
  _mass_sqrt = std::sqrt(_mass);

  _ang_pos[0] = 3;
  _ang_pos[1] = _ang_pos[0] + fragment(0).pos_size();
  _pos_size   = _ang_pos[1] + fragment(1).pos_size();

  _ang_vel[0] = 3;
  _ang_vel[1] = _ang_vel[0] + fragment(0).vel_size();
  _vel_size   = _ang_vel[1] + fragment(1).vel_size();

  _dv_size = _pos_size + _vel_size;

  _tm_dof = 3 + fragment(0).tm_dof() + fragment(1).tm_dof();

  for(frag = 0; frag < 2; ++frag)
    IO::log << frag << "-th fragment: " << fragment(frag) << "\n";      
}

double Molecule::imom (int i) const throw(Error::General)
{
    const char funame [] = "Molecule::imom: ";

    switch(type()) {
    case MONOATOMIC:
      std::cerr << funame << "should not be used for atomic fragment\n";
      throw Error::Logic();
    case LINEAR:
      return _iner_mom[2];
    case NONLINEAR:
      return _iner_mom[i];
    default:
      std::cerr << funame << "unknown molecular type\n";
      throw Error::Logic();
    }
}

double Molecule::imom_sqrt (int i) const throw(Error::General)
{
    const char funame [] = "Molecule::imom_sqrt: ";

    switch(type()) {
    case MONOATOMIC:
      std::cerr << funame << "should not be used for atomic fragment\n";
      throw Error::Logic();
    case LINEAR:
      return _imom_sqrt[2];
    case NONLINEAR:
      return _imom_sqrt[i];
    default:
      std::cerr << funame << "unknown molecular type\n";
      throw Error::Logic();
    }
}

void Molecule::set_dvd (const D3::Vector& torque, const double* pos, const double* vel,
			double* pos_drv, double* vel_drv) const throw(Error::General)
{
  const char funame [] = "Molecule::set_dvd: ";

  double dtemp;
  D3::Vector vtemp;
  const double* qv;
  double*      qvd;

  switch(type()) {
  case MONOATOMIC:
    return;

  case LINEAR:                    // laboratory frame is assumed !!!
  // angular vector derivative
  D3::vprod(vel, pos, pos_drv);
  // angular velocity derivative
  for(int i = 0; i < 3; ++i)
    vel_drv[i] = torque[i] / imom(2);
  return;

  case NONLINEAR:                 // molecular frame is assumed !!!

    /*
      pos_drv [0] = (- pos[1] * vel[0] - pos[2] * vel[1] - pos[3] * vel[2]) / 2.0;
      pos_drv [1] = (- pos[3] * vel[1] + pos[2] * vel[2] + pos[0] * vel[0]) / 2.0;
      pos_drv [2] = (- pos[1] * vel[2] + pos[3] * vel[0] + pos[0] * vel[1]) / 2.0;
      pos_drv [3] = (- pos[2] * vel[0] + pos[1] * vel[1] + pos[0] * vel[2]) / 2.0;
    */

    qv  = pos     + 1; // vector part of the quaternion
    qvd = pos_drv + 1;

    // dq0 = - qv * omega / 2.
    *pos_drv = - vdot(qv, vel, 3) / 2.;

    // dqv = (qv x omega + q0 * omega) / 2.
    D3::vprod(qv, vel, vtemp); 
    for(int i = 0; i < 3; ++i)
      qvd[i] = (vtemp[i] + *pos * vel[i]) / 2.;

    // rotational frequencies derivatives
    if(top() == SPHERICAL) {
      for(int i = 0; i < 3; ++i)
	vel_drv[i] = torque [i]/imom(i);
    }
    else 
      for(int i = 0; i < 3; ++i) {
	dtemp = (imom((i + 1) % 3) - imom((i + 2) % 3)) / imom(i);               // euler coefficients
	vel_drv[i] = torque [i]/imom(i)  + dtemp * vel[(i + 1) % 3] * vel[(i + 2) % 3]; 
      }
     
    return;

  default:
    std::cerr << funame << "unknown molecular type\n";
    throw Error::Logic();
  }
}

// reads molecular geometry from the input stream and uniquely 
// orients the molecular structure to the principal axes 
void Molecule::read(std::istream& from) throw(Error::General)
{
  const char funame [] = "Molecule::read: ";

  static const double rel_tol = 1.e-5;

  double     dtemp, dmax;
  int        itemp;
  D3::Vector vtemp;
  bool       btemp;

  KeyGroup MoleculeGroup;

  Key geom_key("Geometry");             // molecular geometry
  Key flip_key("FlipAxis");             // 180-degree rotation around flip axis in the standard orientation
  Key unit_key("Angstrom");             // use angstrom as distance units
  Key char_key("Charge[au]");           // charge
  Key  dip_key("DipoleMoment[au]");     // dipole moment
  Key quad_key("QuadrupoleMoment[au]"); // quadrupole moment
  Key  pol_key("Polarizability[au]");   // dipolar polarizability
  
  std::string token, line, comment, unit, stemp;
  int flip_axis = -1;

  std::vector<double> dipol, quad, polar;

  // read name 
  from >> _name;
  std::getline(from, comment);

  bool angstrom = false;
  while(from >> token)
    if(token == IO::end_key())
      break;
    else if(token == geom_key) { // geometry section
      if(!(from >> itemp)) {
	std::cerr << funame << "could not read number of atoms for the " << _name << " fragment\n";
	throw Error::Form();
      }
      std::getline(from, comment);

      resize(itemp);
      for(iterator at = begin(); at != end(); ++at)
	from >> (*at);
    }
    else if(token == unit_key) { //  use angstrom unit
      angstrom = true;
    }
    else if(token == flip_key) { //  which axes to flip
      if(!(from >> stemp)) {
	std::cerr << funame << "corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      switch(stemp[0]) {
      case 'X':
	flip_axis = 0;
	break;
      case 'Y':
	flip_axis = 1;
	break;
      case 'Z':
	flip_axis = 2;
	break;
      default:
	std::cerr << funame << token << ": unknown axis symbol: " << stemp[0] << ": use X, Y, or Z\n";
	throw Error::Range();
      }
    }
    else if(token == char_key) { // charge
      from >> _charge;
      if(!from) {
	std::cerr << funame << token << ": input stream is corrupted\n";
	throw Error::Input();
      }
    }
    else if(token == dip_key) { // dipole moment
      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> dtemp)
	_dipole.push_back(dtemp);
    }
    else if(token == quad_key) { // quadrupole moment
      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> dtemp)
	_quadrupole.push_back(dtemp);
    }
    else if(token == pol_key) { // polarizability
      std::getline(from, line);
      std::istringstream iss(line);
      while(iss >> dtemp)
	_polarizability.push_back(dtemp);
    }
    else {
      std::cerr << funame << "unknown key: " << token << "\nAvailable keys:\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }  
  
  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Form();
  }

  if(!size()) {
    std::cerr << funame << "no atoms found\n";
    throw Error::Logic();
  }

  if(angstrom) { // scale coordinates
    for(iterator at = begin(); at != end(); ++at)
      *at *= Phys_const::angstrom;
  }

  if(_charge && (_dipole.size() || _quadrupole.size() || _polarizability.size())) {
    std::cerr << funame << "ERROR: with charge no other multipole moments shoud be defined\n";
    throw Error::Logic();
  }

  /****************** Atom *******************/

  if(size() == 1) {//atom
    _type = MONOATOMIC;
    _top = SPHERICAL;
    _mass = begin()->mass();
    (D3::Vector&)(*begin()) = 0.;

    // multipole moments
    if(_dipole.size() || _quadrupole.size()) {
      std::cerr << funame << "ERROR: there shoud not be dipole or quadrupole moments for monoatomic fragment\n";
      throw Error::Logic();
    }

    if(_polarizability.size() > 1) {
      std::cerr << funame << "ERROR: polarizability for monoatomic fragment should have one component: isotropic\n";
      throw Error::Logic();
    }
    
    return;
  }// atom

  /****************** Atom end *******************/

  // center of mass position
  _std_frame.origin = 0.;
  _mass = 0.;
  for (const_iterator cat = begin(); cat != end(); ++cat) {
    _mass  += cat->mass();
    _std_frame.origin += cat->mass() * (*cat);
  }
  _std_frame.origin /= _mass;

  for (iterator at = begin(); at != end(); ++at)
    *at -= _std_frame.origin;

  // inertia moment matrix
  Lapack::SymmetricMatrix imm (3);
  imm = 0.;
  for (const_iterator cat = begin(); cat != end(); ++cat)
    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j)
	imm(i, j) -= cat->mass() * (*cat)[i] * (*cat)[j];

  dtemp = 0.;
  for(int i = 0; i < 3; ++i)
    dtemp -= imm(i, i);

  for(int i = 0; i < 3; ++i)
    imm(i, i) += dtemp;

  // principal inertia moments and vectors
  Lapack::Vector tim(3);
  Lapack::Matrix evec(3);
  tim = imm.eigenvalues(&evec);
  
  // the transformation is a rotation
  if(Lapack::LU(evec).det() < 0.)
    evec.column(1) *= -1.;
    
  // standard orientation matrix
  for (int i = 0; i < 3; ++i)
    _std_frame.orient.row(i) = evec.column(i);

  // manual flip
  if(flip_axis >= 0)
    for(int i = 0; i < 3; ++i)
      if(i != flip_axis)
	_std_frame.orient.row(i) *= -1.;
 
  // symmetry group
  std::string symm_group_name;
  std::vector<Symmetry::SpaceElement> symm_group_base(1, Symmetry::unity());

 /****************** Linear molecule *******************/

  if( tim[0] / tim[2] < tolerance) {// linear molecule

    _type = LINEAR;
    _top = PROLATE;

    for (iterator at = begin(); at != end(); ++at) {
      (*at)[0] = vdot(_std_frame.orient.row(0) ,(const double*)(*at));

      for(int i = 1; i < 3; ++i)
	(*at)[i] = 0.;
    }

    dtemp = std::sqrt(tim[2]);
    for(int i = 0; i < 3; ++i) {
      _iner_mom [i] = tim[2];
      _imom_sqrt[i] = dtemp;
    }

    // multipole moments 
    if(_dipole.size() > 1) { 
      std::cerr << funame << "ERROR: dipole moment for linear molecule should have one component\n";
      throw Error::Logic();
    }
    if(_quadrupole.size() > 1) {
      std::cerr << funame << "ERROR: quadrupole moment for linear molecule should have one component\n";
      throw Error::Logic();
    }
    if(_polarizability.size() > 2) {
      std::cerr << funame << "ERROR: polarizability for linear molecule may have two components: isotropic and anisotropic\n";
      throw Error::Logic();
    }
        
    if(_dipole.size()) {
      std::cerr << funame << "WARNING: orientation of the symmetry axis  may have flipped; "
	"check the sign of the dipole\n";
    }

    // symmetry
    itemp = size() / 2;
    for(int i = 0; i < itemp; ++i) {
      const Atom& a = *(begin() + i);
      const Atom& b = *(begin() + size() - 1 - i);
      if(a != b || !are_equal(a[0], -b[0], tolerance)) {
	itemp = 0;
	break;
      }
    }
    
    if(itemp) {
      symm_group_base.push_back(Symmetry::inversion());
      symm_group_name = "I";
      symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, 2));
    }
    else {
      //symm_group_base.push_back(Symmetry::unity());
      symm_group_name = "E";
      symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, 1));
    }
    return;
  }// linear molecule

  /********************* NONLINEAR MOLECULE ***********************/

  _type = NONLINEAR;

  std::set<Permutation> full_perm_group  = ::permutation_symmetry_group(*this, tolerance);

  // find reference atoms
  ref_group.resize(2, -1);
  dmax = 0;
  for(const_iterator a = begin(); a != end(); ++a)
    for(const_iterator b = begin(); b !=     a; ++b) {
      dtemp = vdot(D3::vprod(*a, *b));
      if(dtemp > dmax) {
	dmax = dtemp;
	ref_group[0] = a - begin();
	ref_group[1] = a - begin();
      }
    }
  
  if(std::sqrt(dmax) < tolerance) {
    std::cerr << funame << "cannot find non-colinear atoms\n";
    throw Error::Logic();
  }

  /************************ SPHERICAL TOP *************************/

  if((tim[2] - tim[0]) / tim[2] < tolerance) {// spherical top
    
    _top = SPHERICAL;

    for(int i = 0; i < 3; ++i) {
      _iner_mom  [i] = tim[2];
      _imom_sqrt [i] = std::sqrt(tim[2]);
    }

    // multipole moments
    if(_dipole.size() || _quadrupole.size()) {
      std::cerr << funame << "there shoud not be dipole or quadrupole moments for spherical top\n";
      throw Error::Logic();
    }

    if(_polarizability.size() > 1) {
      std::cerr << funame << "polarizability for spherical top may have one component: isotropic one\n";
      throw Error::Logic();
    }

    ref_orient = D3::Matrix((*this)[ref_group[0]], (*this)[ref_group[1]]);

    // symmetry group    
    itemp = full_perm_group.size();
    if(itemp != 12 && itemp != 24 && itemp != 48 && itemp != 60 && itemp != 120) {
      std::cerr << funame << "unknown symmetry group dimension for spherical top: " << itemp << "\n";
      throw Error::Logic();
    }

    std::vector<int> max_volume_group(3);

    double max_volume = 0.;
    for(int i = 0; i < size(); ++i)
      for(int j = 0; j < i; ++j)
	for(int k = 0; k < j; ++k) { 
	  dtemp = D3::volume((*this)[i], (*this)[j], (*this)[k]);
	  if(dtemp > max_volume) {
	    max_volume_group[0] = i;
	    max_volume_group[1] = j;
	    max_volume_group[2] = k;
	    max_volume = dtemp;
	  }
	  else if(dtemp < -max_volume) {
	    max_volume_group[0] = j;
	    max_volume_group[1] = i;
	    max_volume_group[2] = k;
	    max_volume = -dtemp;
	  }
	}

    // no reorientation is necessary
    _std_frame.orient = 1.;

    // symmetry classes
    std::vector<std::set<Permutation> > symm_class;
    for(std::set<Permutation>::const_iterator git = full_perm_group.begin(); git != full_perm_group.end(); ++git) {
      btemp = false;
      for(std::vector<std::set<Permutation> >::const_iterator sit = symm_class.begin(); sit != symm_class.end(); ++sit)
	if(sit->find(*git) != sit->end()) {
	  btemp = true;
	  break;
	}
      if(btemp)
	continue;

      symm_class.push_back(std::set<Permutation>());
      for(std::set<Permutation>::const_iterator it = full_perm_group.begin(); it != full_perm_group.end(); ++it)
	symm_class.back().insert((*it) * (*git) * it->invert());      
    }

    std::cerr << "there are " << symm_class.size() << " classes:";
    for(int i = 0; i < symm_class.size(); ++i)
      std::cerr << "   " << symm_class[i].size();

    return;
  }// spherical top

  /****************************************************************/

  // find symmetry axis
  int symm_axis = -1;
  if((tim[1] - tim[0])/tim[1] < tolerance)
    symm_axis = 2;
  else if((tim[2] - tim[1])/tim[1] < tolerance)
    symm_axis = 0;

  /************************ SYMMETRIC TOP *************************/

  if(symm_axis >= 0) {// symmetric top

    // molecule type
    if(!symm_axis) 
      _top = PROLATE;
    else
      _top = OBLATE;

    D3::Vector y_axis;
    for (const_iterator at = begin(); at != end(); ++at) {
      D3::vprod(&_std_frame.orient(symm_axis, 0), *at, y_axis);
      dtemp = normalize(y_axis);
      if(at == begin() || dtemp > dmax) {
	dmax = dtemp;
	_std_frame.orient.row(1) = y_axis;
      }
    }

    if(dmax <  tolerance) {
      std::cerr << funame << "could not find an atom off the symmetry axis\n";
      throw Error::Run();
    }

    if(!symm_axis)
      D3::vprod(&_std_frame.orient(0, 0), &_std_frame.orient(1, 0), &_std_frame.orient(2, 0));
    else
      D3::vprod(&_std_frame.orient(1, 0), &_std_frame.orient(2, 0), &_std_frame.orient(0, 0));
      

    _iner_mom[symm_axis] = tim[symm_axis];
    _iner_mom[1] = _iner_mom[2 - symm_axis] = tim[1];

    // multipole moments 
    if(_dipole.size() > 1) { 
      std::cerr << funame << "ERROR: dipole moment for symmetric top should have one component\n";
      throw Error::Logic();
    }
    if(_quadrupole.size() > 1) {
      std::cerr << funame << "ERROR: quadrupole moment for symmetric top should have one component\n";
      throw Error::Logic();
    }
    if(_polarizability.size() > 2) {
      std::cerr << funame 
		<< "ERROR: polarizability for symmetric top may have one or two components: isotropic and anisotropic\n";
      throw Error::Logic();
    }

    if(_dipole.size()) {
      std::cerr << funame << "WARNING: orientation of the symmetry axis may have flipped; "
	"check the sign of the dipole\n";
    }
  }// symmetric top

  /************************ ASYMMETRIC TOP *************************/

  else {// asymmetric top

    _top = ASYM;

    for (int i = 0; i < 3; ++i)
      _iner_mom[i] = tim[i];

    // multipole moments 
    switch(_dipole.size()) {
    case 0:
      break;
    case 3:
      break;
    default:
      std::cerr << funame << "ERROR: dipole moment for asymmetric top should have three components: dx, dy, dz\n";
      throw Error::Logic();
    }

    switch(_quadrupole.size()) {
    case 0:
      break;
    case 5:
      _quadrupole.push_back(- _quadrupole[0] - _quadrupole[3]);
      break;
    default:
      std::cerr << funame 
		<< "ERROR: quadrupole moment for asymetric top should have five components: qxx, qxy, qxz, qyy, qyz\n";
      throw Error::Logic();
    }

    switch(_polarizability.size()) {
    case 0:
      break;
    case 1:
      break;
    case 6:
      break;
    default:
      std::cerr << funame << "ERROR: polarizability for asymetric top may have either six"
	" (pxx, pxy, pxz, pyy, pyz, pzz) or one (isotropic) component\n";
      throw Error::Logic();
    }

    // reoriented multipole moments

    // dipole
    if(_dipole.size() == 3) {
      D3::Vector dv;
      for(int i = 0; i < 3; ++i)
	dv[i] = _dipole[i];

      _dipol_vec = _std_frame.orient * dv;

      for(int i = 0; i < 3; ++i)
	_dipole[i] = _dipol_vec[i];
    }

    // quadrupole
    if(_quadrupole.size()== 6) {
      D3::Matrix mm;
      int ii = 0;
      for(int i = 0; i < 3; ++i)
	for(int j = i; j < 3; ++j) {
	  mm(i, j) = _quadrupole[ii++];
	  if(i != j)
	    mm(j, i) = mm(i, j);
	}
      _quad_mat = _std_frame.orient * mm * _std_frame.orient.transpose();

      ii = 0;
      for(int i = 0; i < 3; ++i)
	for(int j = i; j < 3; ++j) {
	  _quadrupole[ii++] = _quad_mat(i, j);
	}
    }

    // polarizability
    if(_polarizability.size() == 6) {
      D3::Matrix mm;
      int ii = 0;
      for(int i = 0; i < 3; ++i)
	for(int j = i; j < 3; ++j) {
	  mm(i, j) = _polarizability[ii++];
	  if(i != j)
	    mm(j, i) = mm(i, j);
	}
      _polar_mat = _std_frame.orient * mm * _std_frame.orient.transpose();

      ii = 0;
      for(int i = 0; i < 3; ++i)
	for(int j = i; j < 3; ++j) {
	  _polarizability[ii++] = _polar_mat(i, j);
	}
    }
  }// asymmetric top

  // atom coordinates in the molecular frame
  for(iterator at = begin(); at != end(); ++at)
    (D3::Vector&)(*at) = _std_frame.orient * (*at);

  ref_orient = D3::Matrix((*this)[ref_group[0]], (*this)[ref_group[1]]);

  for (int i = 0; i < 3; ++i)
    _imom_sqrt[i] = std::sqrt(_iner_mom[i]);

  /************************************* SYMMETRY GROUP *************************************/ 

  // is molecule plane
  bool is_plain = true;
  for(const_iterator at = begin(); at != end(); ++at) {
    dtemp = (*at)[2];
    dtemp = dtemp < 0. ? -dtemp : dtemp;
    if(dtemp > tolerance) {
      is_plain = false;
      break;
    }
  }

  // symmetric top
  if(symm_axis >= 0) {

    // plane molecule
    if(is_plain) {

      if(symm_axis != 2) {
	std::cerr << funame << "plain symmetric top: symmetry axis should be Z\n";
	throw Error::Logic();
      }

      if(!is_symmetric(*this, Symmetry::reflection(symm_axis), tolerance)) {
	std::cerr << funame << "plain symmetric top: no h-plane symmetry\n";
	throw Error::Logic();
      }

      // add the h-plane symmetry element
      symm_group_base.push_back(Symmetry::reflection(symm_axis));

      // major rotation order
      int rotation_order = 0;
      for(std::set<Permutation>::const_iterator pit = full_perm_group.begin(); pit != full_perm_group.end(); ++pit) {
	itemp = pit->orbit_max().size();
	if(pit == full_perm_group.begin() || itemp > rotation_order)
	  rotation_order = itemp;
      }
      
      if(rotation_order < 3) {
	std::cerr << funame << "plain symmetric top: for symmetric top major rotation order should be >= 3\n";
	throw Error::Logic();
      }
      
      if(full_perm_group.size() % rotation_order) {
	std::cerr << funame << "plain symmetric top: rotation order should divide the permutation group size in whole\n";
	throw Error::Logic();
      }

      std::ostringstream oss;
      oss << rotation_order;
      std::string rotation_order_string = oss.str();

      // major_rotation
      Symmetry::SpaceElement major_rotation = Symmetry::rotation(symm_axis, 2. * M_PI / (double)rotation_order);

      // check that the major rotation is the symmetry element
      Permutation rot_perm = is_symmetric(*this, major_rotation, tolerance);
      if(!rot_perm) {
	std::cerr << funame << "plain symmetric top: " 
		  << rotation_order << "-fold rotation is not a symmetry element\n";
	throw Error::Logic();
      }

      // add major rotation element to the symmetry group base
      symm_group_base.push_back(major_rotation);

      // remove major rotations from symmetry permutations
      Permutation perm(size());
      for(int r = 0; r < rotation_order; ++r, perm = rot_perm * perm) {
	std::set<Permutation>::iterator pit = full_perm_group.find(perm);
	if(pit == full_perm_group.end()) {
	  std::cerr << funame << "plain symmetric top: rotation element not found\n";
	  throw Error::Logic();
	}
	else
	  full_perm_group.erase(pit);
      }

      // no in-plane rotations / reflections
      if(!full_perm_group.size()) {
	// symmetry group is C_nh
	symm_group_name = "C" + rotation_order_string + "h";
	symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, 2 * rotation_order));

	return;
      }

      // there are in-plane rotations / v-reflections
      if(full_perm_group.size() != rotation_order) {
	std::cerr << funame << "plain symmetric top: symmetry group unknown\n";
	throw Error::Logic();
      }

      // find a pair of atoms coupled by v-reflection
      std::vector<int> orbit = full_perm_group.begin()->orbit_max();
      if(orbit.size() != 2) {
	std::cerr << funame << "plain symmetric top: not v-reflection\n";
	throw Error::Logic();
      }
      
      D3::Vector n = (*this)[orbit[1]] - (*this)[orbit[0]];
      n[symm_axis]  =  0.;
      normalize(n);
      
      // check that v-reflection is a symmetry element
      if(*full_perm_group.begin() != is_symmetric(*this, Symmetry::reflection(n), tolerance)) {
	std::cerr << funame << "plain symmetric top: v-reflection is not a symmetry element\n";
	throw Error::Logic();
      }
      
      // add v-reflection to the symmetry elements
      symm_group_base.push_back(Symmetry::reflection(n));

      // symmetry group is D_nh
      symm_group_name = "D" + rotation_order_string + "h";
      symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, 4 * rotation_order));

      return;
    }// plain molecule

    // non-plain molecule

    // find atoms not in plane
    dmax = 0.;
    std::vector<int> frame_group(3);
    for(const_iterator a = begin(); a != end(); ++a)
      for(const_iterator b = begin(); b   !=   a; ++b)
	for(const_iterator c = begin(); c   !=   b; ++c) {
	  dtemp = D3::volume(*a, *b, *c);
	  if(dtemp > dmax) {
	    frame_group[0] = a - begin();
	    frame_group[1] = b - begin();
	    frame_group[2] = c - begin();
	    dmax  = dtemp;
	  }
	  else if(dtemp < -dmax) {
	    frame_group[1] = a - begin();
	    frame_group[0] = b - begin();
	    frame_group[2] = c - begin();
	    dmax  = -dtemp;
	  }
	}

    if(dmax < tolerance) {
      std::cerr << funame << "non-plain symmetric top: cannot find out-of-plane atoms\n";
      throw Error::Logic();
    }

    // extract rotations into separate permutation subgroup 
    std::set<Permutation> rot_perm_group;
    for(std::set<Permutation>::const_iterator pit = full_perm_group.begin(); pit != full_perm_group.end(); ++pit)
      if(D3::volume((*this)[(*pit)[frame_group[0]]],
		    (*this)[(*pit)[frame_group[1]]],
		    (*this)[(*pit)[frame_group[2]]]) > 0.)
	rot_perm_group.insert(*pit);
    
    for(std::set<Permutation>::const_iterator pit = rot_perm_group.begin(); pit != rot_perm_group.end(); ++pit)
      full_perm_group.erase(*pit);

    if(full_perm_group.size() && full_perm_group.size() != rot_perm_group.size()) {
      std::cerr << funame << "non-plain symmetric top: numbers of symmetry elements with and without inversion mismatch\n";
      throw Error::Logic();
    }

    // major rotation order
    int rotation_order = 0;
    for(std::set<Permutation>::const_iterator pit = rot_perm_group.begin(); pit != rot_perm_group.end(); ++pit) {
      itemp = pit->orbit_max().size();
      if(pit == full_perm_group.begin() || itemp > rotation_order)
	rotation_order = itemp;
    }
      
    if(rotation_order < 3) {
      std::cerr << funame << "non-plain symmetric top: major rotation order should be >= 3\n";
      throw Error::Logic();
    }
      
    if(rot_perm_group.size() % rotation_order) {
      std::cerr << funame << "non-plain symmetric top: rotation order should divide the rotational subgroup size in whole\n";
      throw Error::Logic();
    }

    std::ostringstream oss;
    oss << rotation_order;
    std::string rotation_order_string = oss.str();

    // major_rotation
    Symmetry::SpaceElement major_rotation = Symmetry::rotation(symm_axis, 2. * M_PI / (double)rotation_order);
    
    // check that the major rotation is the symmetry element
    Permutation rot_perm = is_symmetric(*this, major_rotation, tolerance);
    if(!rot_perm) {
      std::cerr << funame << "non-plain symmetric top: " << rotation_order << "-fold rotation is not a symmetry element\n";
      throw Error::Logic();
    }

    // add the major rotation symmetry element to the symmetry group basis
    symm_group_base.push_back(major_rotation);

    // remove major rotation elements from the group of rotations
    Permutation perm(size());
    for(int r = 0; r < rotation_order; ++r, perm = rot_perm * perm) {
      std::set<Permutation>::iterator pit = rot_perm_group.find(perm);
      if(pit == rot_perm_group.end()) {
	std::cerr << funame << "non-plain symmetric top: rotation element not found\n";
	throw Error::Logic();
      }
      else
	rot_perm_group.erase(pit);
    }

    // plane symmetry element
    Permutation plane_perm;
    if(full_perm_group.size())
      plane_perm = is_symmetric(*this, Symmetry::reflection(symm_axis), tolerance);

    if(plane_perm) {
      // add the symmetry plane reflection element to the symmetry group basis
      symm_group_base.push_back(Symmetry::reflection(symm_axis));

      // remove major mirrored rotations from reflection permutations
      perm = Permutation(size());
      for(int r = 0; r < rotation_order; ++r, perm = rot_perm * perm) {
	std::set<Permutation>::iterator pit = full_perm_group.find(perm * plane_perm);
	if(pit == full_perm_group.end()) {
	  std::cerr << funame << "non-plain symmetric top: mirrored rotation element not found\n";
	  throw Error::Logic();
	}
	else
	  full_perm_group.erase(pit);
      }
    }

    // only major rotations are available (no in-plane rotations)
    if(!rot_perm_group.size()) {
      
      // either  h-reflection and corresponding mirrored rotations or no reflections at all
      if(!full_perm_group.size()) {
	if(plane_perm) {
	  itemp =  2 * rotation_order;
	  // symmetry group is C_nh
	  symm_group_name = "C" + rotation_order_string + "h";
	}
	else {
	  itemp = rotation_order;
	  // symmetry group is C_n
	  symm_group_name = "C" + rotation_order_string;
	}
	// initialize symmetry group
	symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, itemp));

	return;
      }

      // no plane symmetry but there are v-reflections

      // find a pair of atoms coupled by v-reflection
      std::vector<int> orbit = full_perm_group.begin()->orbit_max();
      if(orbit.size() != 2) {
	std::cerr << funame << "non-plain symmetric top: not v-reflection\n";
	throw Error::Logic();
      }
      
      D3::Vector n = (*this)[orbit[1]] - (*this)[orbit[0]];
      n[symm_axis]  =  0.;
      normalize(n);
      
      // check that v-reflection is a symmetry element
      if(*full_perm_group.begin() != is_symmetric(*this, Symmetry::reflection(n), tolerance)) {
	std::cerr << funame << "non-plain symmetric top: v-reflection is not a symmetry element\n";
	throw Error::Logic();
      }
      
      // add v-reflection to the symmetry elements
      symm_group_base.push_back(Symmetry::reflection(n));
      symm_group_name = "C" + rotation_order_string + "v";
      // symmetry group is C_nv
      symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, 2 * rotation_order));
      
      return;
    }

    // there are in-plane rotations
    if(rot_perm_group.size() != rotation_order) {
      std::cerr << funame << "non-plain symmetric top: symmetry group unknown\n";
      throw Error::Logic();
    }

    // find the axis of one of the in-plane rotations
    D3::Vector n;
    dmax = 0.;
    std::set<std::vector<int> > orbit = rot_perm_group.begin()->orbit();
    for(std::set<std::vector<int> >::const_iterator ot = orbit.begin(); ot != orbit.end(); ++ot) 
      switch(ot->size()) {
      case 1:
	vtemp = (*this)[(*ot)[0]];
	dtemp = vlength(vtemp);
	if(dtemp > dmax) {
	  n    = vtemp;
	  dmax = dtemp;
	}

	break;
      case 2:
	vtemp = (*this)[(*ot)[0]] + (*this)[(*ot)[1]];
	dtemp = vlength(vtemp) / 2.;
	if(dtemp > dmax) {
	  n    = vtemp;
	  dmax = dtemp;
	}

	break;
      default:
	std::cerr << funame << "non-plain symmetric top: not in-plane rotation\n";
	throw Error::Logic();
      }
      
    if(dmax < tolerance) {
      std::cerr << funame << "non-plain symmetric top: in-plane rotation axis not found\n";
      throw Error::Logic();
    }

    n[symm_axis] = 0.;
    normalize(n);

    // check that in-plane rotation is a symmetry element
    if(*rot_perm_group.begin() != is_symmetric(*this, Symmetry::rotation(n), tolerance)) {
      std::cerr << funame << "non-plain symmetric top: in-plane rotation is not a symmetry element\n";
      throw Error::Logic();
    }

    // add in-plane rotation to the symmetry group base
    symm_group_base.push_back(Symmetry::rotation(n));

    // h-reflection symmetry or no elements with inversion at all
    if(plane_perm || !full_perm_group.size()) {
      if(plane_perm) {
	// symmetry group is D_nh
	itemp = 4 * rotation_order;
	symm_group_name = "D" + rotation_order_string + "h";
      }
      else {
	// symmetry group is D_n
	itemp = 2 * rotation_order;
	symm_group_name = "D" + rotation_order_string;
      }
      // initialize symmetry group
      symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, itemp));

      return;
    }

    // no h-reflection symmetry but there are in-plane rotations, v-reflections and mirrored rotations

    // mirrored rotation to half-angle
    Symmetry::SpaceElement mirr_rot = Symmetry::rotation(symm_axis, M_PI / (double)rotation_order) 
      * Symmetry::reflection(symm_axis);

    // check if the mirrored rotation to half-angle is the symmetry element
    if(!is_symmetric(*this, mirr_rot, tolerance)) {
      std::cerr << funame << "non-plain symmetric top: mirrored rotation is not a symmetry element\n";
      throw Error::Logic();
    }
    
    // add mirrored rotation to symmetry group base
    symm_group_base.push_back(mirr_rot);

    // D_nd symmetry group
    symm_group_name = "D" + rotation_order_string + "d";
    symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base, 4 * rotation_order));

    return;
  } // symmetric top

  // asymmetric top
  int ref = 0, rot = 0;
  for(int i = 0; i < 3; ++i) {
    if(is_symmetric(*this, Symmetry::reflection(i), tolerance)) {
      symm_group_base.push_back(Symmetry::reflection(i));
      ++ref;
    }
    if(is_symmetric(*this, Symmetry::rotation(i), tolerance)) {
      symm_group_base.push_back(Symmetry::rotation(i));
      ++rot;
    }
  }

  int inv = 0;
  if(is_symmetric(*this, Symmetry::inversion(), tolerance)) {
    symm_group_base.push_back(Symmetry::inversion());
    ++inv;
  }

  if(ref == 0      && rot == 0 && inv == 0)
    symm_group_name = "C1";
  else if(ref == 1 && rot == 0 && inv == 0)
    symm_group_name = "Cs";
  else if(ref == 0 && rot == 1 && inv == 0)
    symm_group_name = "C2";
  else if(ref == 0 && rot == 0 && inv == 1)
    symm_group_name = "Ci";
  else if(ref == 2 && rot == 1 && inv == 0)
    symm_group_name = "C2v";
  else if(ref == 1 && rot == 1 && inv == 1)
    symm_group_name = "C2h";
  else if(ref == 0 && rot == 3 && inv == 0)
    symm_group_name = "D2";
  else if(ref == 3 && rot == 3 && inv == 1)
    symm_group_name = "D2h";
  else {
    std::cerr << funame << "asymmetric top: symmetry group is unknown\n";
    throw Error::Logic();
  }
 
  symmetry_group.init(new Symmetry::SpaceGroup(symm_group_name, symm_group_base));

}// read

bool Molecule::is_equal (std::vector<Atom> mol, D3::Vector& shift, Quaternion& orient) const
{
  int        itemp;
  double     dtemp;
  D3::Vector vtemp;

  if(mol.size() != size())
    return false;

  // molecular displacement
  shift = 0.;
  for(int a = 0; a < size(); ++a)
    shift += mol[a].mass() * mol[a];
  shift /= mass();

  for(int a = 0; a < size(); ++a)
    mol[a] -= shift;

  D3::Matrix mol_orient(mol[ref_group[0]], mol[ref_group[1]]);
  
  mol_orient = ref_orient.transpose() * mol_orient;
  for(int a = 0; a < size(); ++a) {
    mol[a] = mol_orient * mol[a];

    if(!are_equal(mol[a], (*this)[a], tolerance))
      return false;
  }

  orient = (Quaternion)mol_orient;
  return true;  
}

void Molecule::print (std::ostream& to) const throw(Error::General)
{
  const char funame [] = "Molecule::print: ";

  to << name() << "\n";

  switch(type()) {
  case MONOATOMIC:
    to << "   Monoatomic\n";
    break;
  case LINEAR:
    to << "   Linear\n";
    to << "   Inertia moment = " << imom(0) / Phys_const::amu << " amu*bohr*bohr\n";
    break;
  case NONLINEAR:
    switch(top()) {
    case ASYM:
      to << "   Asymmetric top\n";
      break;
    case SPHERICAL:
      to << "   Spherical top\n";
      break;
    case PROLATE:
      to << "   Prolate symmetric top\n";
      break;
    case OBLATE:
      to << "   Oblate symmetric top\n";
    }
    to << "   Inertia moments = {"; 
    for(int i = 0; i < 3; ++i) {
      if(i)
	to << ", "; 
      to << imom(i) / Phys_const::amu; 
    }
    to << "} amu*bohr*bohr\n";
  }

  if(symmetry_group)
    to << "   Symmetry group = " << symmetry_group->name() << 
      ", group size = "          << symmetry_group->size() << "\n";

  to << "   Mass = " << mass() / Phys_const::amu << " amu\n";
  to << "   Standard orientation geometry:\n";
  for(std::vector<Atom>::const_iterator at = begin(); at != end(); ++at)
    to << "      " << *at << "\n";
  to << "\n";

  if(charge() || dipole().size() || quadrupole().size() || polarizability().size()) {// multipole moments
    to  << "   Multipole moments:\n";

    // Charge
    if(charge())
      to << "      Charge = " << charge() << "\n";

    if(dipole().size()) { // dipole

      to << std::setw(15) << "DX"
	 << std::setw(15) << "DY"
	 << std::setw(15) << "DZ"
	 << "\n";

      switch(dipole().size()) {
      case 1:

	if(type() == LINEAR || type() == NONLINEAR && top() == PROLATE)
	  to << std::setw(15) << dipole()[0]
	     << std::setw(15) << 0.
	     << std::setw(15) << 0.
	     << "\n";
	else if(type() == NONLINEAR && top() == OBLATE)
	  to << std::setw(15) << 0.
	     << std::setw(15) << 0.
	     << std::setw(15) << dipole()[0]
	     << "\n";
	else {
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;

      case 3:

	if(type() == NONLINEAR && top() == ASYM) {
	  for(int i = 0; i < 3; ++i)
	    to << std::setw(15) << dipole()[i];
	  to << "\n";
	}
	else {
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;

      default:

	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
      to << "\n";
    } // dipole

    if(quadrupole().size()) { // quadrupole

      to << std::setw(15) << "QXX"
	 << std::setw(15) << "QXY"
	 << std::setw(15) << "QXZ"
	 << std::setw(15) << "QYY"
	 << std::setw(15) << "QYZ"
	 << std::setw(15) << "QZZ"
	 << "\n";

      switch(quadrupole().size()) {
      case 1: 

	if(type() == LINEAR || type() == NONLINEAR && top() == PROLATE) {
	  for(int i = 0; i < 6; ++i)
	    if(i == 0) 
	      to << std::setw(15) << quadrupole()[0];
	    else if(i == 3 || i == 5)
	      to << std::setw(15) << -quadrupole()[0] / 2.;
	    else 
	      to << std::setw(15) << 0.;
	  to << "\n";
	}
	else if(type() == NONLINEAR && top() == OBLATE) {
	  for(int i = 0; i < 6; ++i)
	    if(i == 5) 
	      to << std::setw(15) << quadrupole()[0];
	    else if(i == 3 || i == 0)
	      to << std::setw(15) << -quadrupole()[0] / 2.;
	    else 
	      to << std::setw(15) << 0.;
	  to << "\n";
	}
	else {
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;

      case 6:

	if(type() == NONLINEAR && top() == ASYM) {
	  for(int i = 0; i < 6; ++i)
	    to << std::setw(15) << quadrupole()[i];
	  to << "\n";
	}
	else {
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;

      default:

	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
      to << "\n";
    }// quadrupole

    if(polarizability().size()) { // polarizability

      to << std::setw(15) << "PXX"
	 << std::setw(15) << "PXY"
	 << std::setw(15) << "PXZ"
	 << std::setw(15) << "PYY"
	 << std::setw(15) << "PYZ"
	 << std::setw(15) << "PZZ"
	 << "\n";

      switch(polarizability().size()) {
      case 1:

	for(int i = 0; i < 6; ++i)
	  if(i == 0 || i == 3 || i == 5) 
	    to << std::setw(15) << polarizability()[0];
	  else
	    to << std::setw(15) << 0.;
	to << "\n";
	break;

      case 2: 

	if(type() == LINEAR || type() == NONLINEAR && top() == PROLATE) {
	  for(int i = 0; i < 6; ++i)
	    if(i == 0) 
	      to << std::setw(15) << polarizability()[0] + polarizability()[1];
	    else if(i == 3 || i == 5)
	      to << std::setw(15) << polarizability()[0] - polarizability()[1] / 2.;
	    else 
	      to << std::setw(15) << 0.;
	  to << "\n";
	}
	else if(type() == NONLINEAR && top() == OBLATE) {
	  for(int i = 0; i < 6; ++i)
	    if(i == 5) 
	      to << std::setw(15) << polarizability()[0] + polarizability()[1];
	    else if(i == 3 || i == 0)
	      to << std::setw(15) << polarizability()[0] - polarizability()[1] / 2.;
	    else 
	      to << std::setw(15) << 0.;
	  to << "\n";
	}
	else {
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;

      case 6:

	if(type() == NONLINEAR && top() == ASYM) {
	  for(int i = 0; i < 6; ++i)
	    to << std::setw(15) << polarizability()[i];
	  to << "\n";
	}
	else {
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;

      default:

	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
      to << "\n";
    }// polarizability

  }// multipole moments
}
