#include "potential.hh"
#include "units.hh"
#include "key.hh"

#include <sstream>

Potential::Wrap default_pot;

void Potential::Wrap::read (std::istream& from) throw(Error::General)
{
  const char funame [] = "Potential::Wrap::read: ";

  KeyGroup PotentialWrap;

  Key  anal_key("Analytic"       );
  //Key  harm_key("Harmonic"       );
  Key    cl_key("ChargeLinear"   );
  Key    cn_key("ChargeNonlinear");
  Key    dd_key("DipoleDipole"   );
  Key multi_key("Multipole"      );

  if(_fun) {
    std::cerr << funame << "potential has been already initialized\n";
    throw Error::Init();
  }

  std::string token, comment;
  while(from >> token) {
    if(anal_key == token) {// analytic
      std::getline(from, comment);
      _fun = ConstSharedPointer<Base>(new Analytic(from));
      return;
    }// analytic
    /*
    else if(harm_key == token) {// Harmonic
      std::getline(from, comment);
      _fun = ConstSharedPointer<Base>(new Harmonic(from));
      return;
    }
    */
    else if(cl_key == token) {// charge-linear molecule
      std::getline(from, comment);
      _fun = ConstSharedPointer<Base>(new ChargeLinear(from));
      return;
    }// charge-linear molecule
    else if(cn_key == token) {// charge-nonlinear molecule
      std::getline(from, comment);
      _fun = ConstSharedPointer<Base>(new ChargeNonlinear(from));
      return;
    }// charge-nonlinear molecule
    else if(dd_key == token) {// dipole-dipole
      std::getline(from, comment);
      _fun = ConstSharedPointer<Base>(new DipoleDipole(from));
      return;
    }// dipole-dipole
    else if(multi_key == token) {// multipole
      std::getline(from, comment);
      _fun = ConstSharedPointer<Base>(new Multipole(from));
      return;
    }// multipole
    else {
      std::cerr << funame << "unknown keyword " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }

  std::cerr << funame << "stream is corrupted\n";
  throw Error::Init();
  
}

/*
Potential::Harmonic::Harmonic (std::istream& from) throw(Error::General)
{
  const char funame [] = "Potential::Harmonic::Harmonic: ";

  int         itemp;
  double      dtemp;
  std::string stemp;

  // layout
  std::vector<int> lt(2);
  lt[0] = 3;
  lt[1] = 4;

  Configuration::State::layout.set(lt);

  // input
  std::ifstream from(file);
  if(!from) {
    std::cerr << funame << "cannot open " << file << " file";
    std::exit(1);
  }

  // monom dimensions
  int rdim, qdim;
  from >> rdim >> qdim; 

  // distance grid
  from >> itemp;
  Array<double> dist_grid(itemp);

  for(int dist = 0; dist < dist_grid.size(); ++dist)
    from >> dist_grid[dist_grid.size() - dist - 1];


  // energy conversion scale
  from >> _convert.scale;

  // harmonic expansion vectors
  itemp = 0;
  for(int pack = 0; pack < 2; itemp += _harmonic_expansion[pack++]->size()) {
    std::ostringstream xname;
    xname << "hex_" << rdim - pack << "_" << qdim << ".vec";
    std::ifstream xin(xname.str().c_str());
    if(!xin) {
      std::cerr << funame << "cannot open " << xname.str() << " file\n";
      std::exit(1);
    } 
    harmonic_expansion[pack].init(new HarmonicExpansion(xin, rdim - pack, qdim));
  }

  const int xsize = itemp;
  expansion_coefficient.resize(xsize);

  // energy expansion coefficients
  Array<double> vtemp(dist_grid.size());
  for(int x = 0; x < xsize; ++x) {
    for(int d = 0; d < dist_grid.size(); ++d)
      from >> vtemp[dist_grid.size() - d - 1];
    expansion_coefficient[x].init(dist_grid, vtemp, dist_grid.size());
  }

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    std::exit(1);
  }
}
*/

Potential::Analytic::Analytic (std::istream& from) throw(Error::General) 
  :  _pot_ener(0), _pot_init(0), _corr_ener(0), _corr_init(0),
     _dist_incr(1.e-4),  _angl_incr(1.e-4)
{    
  const char funame [] = "Potential::Analytic::Analytic: ";

  double dtemp;
  int    itemp;
  
  KeyGroup AnaliticPotential;

  Key  pot_libr_key("Library");
  Key  pot_ener_key("EnergyMethod");
  Key  pot_init_key("InitMethod");
  Key  pot_data_key("InitData");
  Key  pot_rpar_key("ParameterReal");
  Key  pot_ipar_key("ParameterInteger");

  Key corr_libr_key("CorrectionLibrary");
  Key corr_ener_key("CorrectionEnergyMethod");
  Key corr_init_key("CorrectionInitMethod");
  Key corr_data_key("CorrectionInitData");
  Key corr_rpar_key("CorrectionParameterReal");
  Key corr_ipar_key("CorrectionParameterInteger");

  Key dist_key("DistanceIncrement");
  Key angl_key("AngularIncrement");

  std::string token, line, comment, stemp;

  std::string pot_data, corr_data;
  while(from >> token) {// read cycle
    // input end
    if(token == IO::end_key()) {
      std::getline(from, comment);
      break;
    }
    // potential library name
    else if(token == pot_libr_key) {

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      _pot_libr.open(stemp);
    }
    // potential energy method
    else if(token == pot_ener_key) {

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      _pot_ener = (ener_t)_pot_libr.member(stemp);
    }
    // potential initialization method
    else if(token == pot_init_key) {

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      _pot_init = (init_t)_pot_libr.member(stemp);
    }
    // potential initialization data
    else if(token == pot_data_key) {

      if(!(from >> pot_data)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // potential real parameters
    else if(token == pot_rpar_key) {

      IO::LineInput lin(from);
      std::vector<double> rpar;

      while (lin >> dtemp)
	rpar.push_back(dtemp);

      _pot_rpar = rpar;
    }
    // potential integer parameters
    else if(token == pot_ipar_key) {

      IO::LineInput lin(from);
      std::vector<int> ipar;

      while (lin >> itemp)
	ipar.push_back(itemp);

      _pot_ipar = ipar;
    }
    // correction library name
    else if(token == corr_libr_key) {

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      _corr_libr.open(stemp);
    }
    // correction energy method
    else if(token == corr_ener_key) {

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      _corr_ener = (ener_t)_corr_libr.member(stemp);
    }
    // correction initialization method
    else if(token == corr_init_key) {

      if(!(from >> stemp)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      _corr_init = (init_t)_corr_libr.member(stemp);
    }
    // correction initialization data
    else if(token == corr_data_key) {

      if(!(from >> corr_data)) {
	std::cerr << funame << token << ": is corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // correction real parameters
    else if(token == corr_rpar_key) {

      IO::LineInput lin(from);
      std::vector<double> rpar;

      while (lin >> dtemp)
	rpar.push_back(dtemp);

      _corr_rpar = rpar;
    }
    // correction integer parameters
    else if(token == corr_ipar_key) {

      IO::LineInput lin(from);
      std::vector<int> ipar;

      while (lin >> itemp)
	ipar.push_back(itemp);

      _corr_ipar = ipar;
    }
    // distance displacement
    else if(token == dist_key) {

      if(!(from >> _dist_incr)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Init();
      }
      
      if(_dist_incr <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // angular displacement
    else if(token == angl_key) {

      if(!(from >> _angl_incr)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Init();
      }

      if(_angl_incr <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // unknown key
    else {
      std::cerr << funame << "unknown key: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }// read cycle

  // checking input
  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Form();
  }

  if(!_pot_ener) {
    std::cerr << funame << "no energy calculation method was provided\n";
    throw Error::Init();
  }

  if(_pot_init)
    _pot_init(pot_data.c_str());

  if(_corr_init)
    _corr_init(corr_data.c_str());
    
}

void Potential::Analytic::_dc2cart (const Dynamic::Coordinates& dc, Array_2<double>& coord)
{
  int at_shift = Structure::fragment(0).size();
  
  for(int frag = 0; frag < 2; ++frag)
    for(int at = 0; at < Structure::fragment(frag).size(); ++at)
      for(int i = 0; i < 3; ++i)
	if(!frag)
	  coord(i, at)            = dc.rel_pos(frag)[at][i];
	else
	  coord(i, at + at_shift) = dc.rel_pos(frag)[at][i] + dc.orb_pos(i);
}

double Potential::Analytic::_tot_ener (const double* coord) const throw(Error::General)
{
  const char funame [] = "Potential::Analytic::_tot_ener: ";

  int ifail = 0;

  double res = _pot_ener(coord,  _pot_rpar,  _pot_ipar, ifail);
  if(ifail)
    throw Error::Run();

  if(_corr_ener) {
    res += _corr_ener(coord, _corr_rpar, _corr_ipar, ifail);
    if(ifail)
      throw Error::Run();
  }
  

  return res;
}

double Potential::Analytic::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const throw(Error::General)
{
  static const char funame [] = "Potential::Analytic::operator(): ";
  
  int    itemp;
  double dtemp;

  Array_2<double> coord(3, Structure::size());
  _dc2cart(dc, coord);

  if(!torque)
    return _tot_ener(coord);

  D3::Vector& force = *torque;
  torque += 1;

  D3::Vector vtemp;
  D3::Vector atom_force;

  // initialization
  force = 0.;
  for(int frag = 0; frag < 2; ++frag)
    torque[frag] = 0.;
	
  // numerical gradient

  double ener_val = _tot_ener(coord);

  const int sfrag = Structure::fragment(0).size() < Structure::fragment(1).size() ? 0 : 1;         // small fragment
  const int lfrag = 1 - sfrag; // large fragment

  const int at_shift = sfrag ? Structure::fragment(0).size() : 0;

  // force calculation on small fragment
  double incr2 = 2. * _dist_incr;

  for(int i = 0; i < 3; ++i) {
    for(int at = 0; at < Structure::fragment(sfrag).size(); ++at)
      coord(i, at + at_shift) -= _dist_incr;
    force[i] += _tot_ener(coord);

    for(int at = 0; at < Structure::fragment(sfrag).size(); ++at)
      coord(i, at + at_shift) += incr2;
    force[i] -= _tot_ener(coord);

    for(int at = 0; at < Structure::fragment(sfrag).size(); ++at)
      coord(i, at + at_shift) -= _dist_incr;
	
    if(sfrag)
      force[i] /=  incr2;
    else
      force[i] /= -incr2;
  }

  // torque calculation
  if (Structure::fragment(sfrag).type() == Molecule::MONOATOMIC) {// small fragment is an atom
    switch(Structure::fragment(lfrag).type()) {
    case Molecule::LINEAR:
      D3::vprod(force, dc.orb_pos(), torque[lfrag]);
      // enforce orthogonality with the angular vector
      torque[lfrag].orthogonalize(dc.ang_pos(lfrag));
      break;
    case Molecule::NONLINEAR:
      D3::vprod(force, dc.orb_pos(), vtemp);
      dc.lf2mf(lfrag, vtemp, torque[lfrag]);
      break;
    }
    return ener_val;
  }// atom

  // torque on the small fragment
  incr2 = 2. * _angl_incr;
  double cos_val = std::cos(_angl_incr) - 1.;
  double sin_val = std::sin(_angl_incr);

  for(int i = 0; i < 3; ++i) {
    const int i1 = (i + 1) % 3;
    const int i2 = (i + 2) % 3;

    for(int at = 0; at < Structure::fragment(sfrag).size(); ++at) {
      coord(i1, at + at_shift) += cos_val * dc.rel_pos(sfrag)[at][i1] + sin_val * dc.rel_pos(sfrag)[at][i2];
      coord(i2, at + at_shift) += cos_val * dc.rel_pos(sfrag)[at][i2] - sin_val * dc.rel_pos(sfrag)[at][i1];
    }
    torque[sfrag][i] += _tot_ener(coord);

    for(int at = 0; at < Structure::fragment(sfrag).size(); ++at) {
      coord(i1, at + at_shift) -= 2. * sin_val * dc.rel_pos(sfrag)[at][i2];
      coord(i2, at + at_shift) += 2. * sin_val * dc.rel_pos(sfrag)[at][i1];
    }
    torque[sfrag][i] -= _tot_ener(coord);

    if(!sfrag)
      for(int at = 0; at < Structure::fragment(sfrag).size(); ++at) {
	coord(i1, at + at_shift) = dc.rel_pos(sfrag)[at][i1];
	coord(i2, at + at_shift) = dc.rel_pos(sfrag)[at][i2];
      }
    else
      for(int at = 0; at < Structure::fragment(sfrag).size(); ++at) {
	coord(i1, at + at_shift) = dc.rel_pos(sfrag)[at][i1] + dc.orb_pos(i1);
	coord(i2, at + at_shift) = dc.rel_pos(sfrag)[at][i2] + dc.orb_pos(i2);
      }

    torque[sfrag][i] /= incr2;
  }

  // torque on the large fragment
  D3::vprod(force, dc.orb_pos(), vtemp);
  for(int i = 0; i < 3; ++i)
    torque[lfrag][i] = vtemp[i] - torque[sfrag][i];
    
  // convert torques to molecular frame for nonlinear fragments
  for(int frag = 0; frag < 2; ++frag)
    switch(Structure::fragment(frag).type()) {
    case Molecule::LINEAR:
      // enforce orthogonality with the angular vector
      torque[frag].orthogonalize(dc.ang_pos(frag));
      break;
    case Molecule::NONLINEAR:
      dc.lf2mf(frag, torque[frag], vtemp);
      torque[frag] = vtemp;
      break;
    }

  return ener_val;
}

/*************************************************************************
 ******************************** Charge-Linear **************************
 *************************************************************************/

double Potential::ChargeLinear::_pot (double r, double x, _mode_t mode) const throw(Error::General)
{
  static const char funame [] = "Potential::ChargeLinear::_pot: "; 
  static const double c3 = 1./3.;

  const double r2 = r  * r;
  const double r3 = r2 * r;
  const double r4 = r3 * r;
  const double r5 = r4 * r;

  switch(mode) {
  case POT_VALUE: 

    return - _dipole * x / r2 + (0.75 * _quadrupole / r3  - 0.5 * _anisotropic_polarizability / r4)
      * (x * x - c3) - 0.5 * _isotropic_polarizability / r4;

  case DIST_DERIV: // -dV/dR

    return - 2. * _dipole * x / r3 + (2.25 * _quadrupole / r4  - 2. * _anisotropic_polarizability / r5) 
	     * (x * x - c3) - 2. * _isotropic_polarizability / r5;
    
  case ANGLE_DERIV: // -dV/dX, X = cos(theta)

    return _dipole / r2 - (1.5 * _quadrupole / r3  - _anisotropic_polarizability / r4) * x;

  default:
    std::cerr << funame << "wrong case\n";
    throw Error::Logic();
  }
}

double Potential::ChargeLinear::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const throw(Error::General)
{
  D3::Vector nr(dc.orb_pos());
  double dist = nr.normalize();

  D3::Vector nl(dc.ang_pos(1));
  nl.normalize();

  double cos_theta =  vdot(nl, nr);

  // energy calculation
  double ener_val = _pot(dist, cos_theta, POT_VALUE);

  if(!torque)
    return ener_val;
  D3::Vector& force = *torque;
  torque += 1;

  const double vr = _pot(dist, cos_theta, DIST_DERIV);
  const double va = _pot(dist, cos_theta, ANGLE_DERIV);

  // force calculation
  force = nr;
  force *= vr - cos_theta * va / dist;
  D3::Vector vtemp = nl;
  vtemp *= va / dist;
  force += vtemp;

  // torque calculation
  torque[0] = 0.;
  D3::vprod(nl, nr, torque[1]);
  torque[1] *= va;

  return ener_val;
}

Potential::ChargeLinear::ChargeLinear (std::istream& from) throw(Error::General)
  : _charge(1.), _dipole(0.), _quadrupole(0.), _isotropic_polarizability(0.), _anisotropic_polarizability(0.)
{
  static const char funame [] = "Potential::ChargeLinear::ChargeLinear (std::istream&): ";

  if(Structure::fragment(1).type() != Molecule::LINEAR) {
    std::cerr << funame << "second fragment should be linear for this potential\n";
    throw Error::Logic();
  }

  IO::Marker funame_marker(funame);

  KeyGroup ChargeLinearPotential;

  Key  chr_key("Charge[au]");
  Key  dip_key("DipoleMoment[au]");
  Key quad_key("QuadrupoleMoment[au]");
  Key isot_key("IsotropicPolarizability[au]");
  Key anis_key("AnisotropicPolarizability[au]");

  std::string token, line, comment;
  while(from >> token) {// read cycle
    if(token == IO::end_key()) {
      std::getline(from, comment);
      break;
    }
    else if(token == chr_key) {
      from >> _charge;
      if(!from) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    else if(token == dip_key) { // dipole
      from >> _dipole;
      if(!from) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    else if(token == quad_key) { // quadrupole
      from >> _quadrupole;
      if(!from) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    else if(token == isot_key) { // isotropic polarizability
      from >> _isotropic_polarizability;
      if(!from) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    else if(token == anis_key) { // anisotropic polarizability
      from >> _anisotropic_polarizability;
      if(!from) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Input();
    }
  }// read cycle

  IO::log << "Charge-Linear potential:\n"
	    << "  Charge                     = " << _charge                     << "\n"
	    << "  Dipole                     = " << _dipole                     << "\n"
	    << "  Qudrupole                  = " << _quadrupole                 << "\n"
	    << "  Isotropic   Polarizability = " << _isotropic_polarizability   << "\n"
	    << "  Anisotropic Polarizability = " << _anisotropic_polarizability << "\n\n";

  // renormalize multipole moments
  _dipole                     *= _charge;
  _quadrupole                 *= _charge;
  _isotropic_polarizability   *= _charge * _charge;
  _anisotropic_polarizability *= _charge * _charge;
}

/*************************************************************************
 ***************************** Charge-Nonlinear **************************
 *************************************************************************/

double Potential::ChargeNonlinear::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const throw(Error::General)
{
  D3::Vector lfr(dc.orb_pos());
  D3::Vector mfr; // r in molecular frame
  dc.lf2mf(1, lfr, mfr);

  const double r = lfr.vlength();
  const double r2 = 1. / r / r;
  const double r3 = r2 / r;
  const double r5 = r3 * r2;
  const double r6 = r5 / r;
  const double r7 = r6 / r;
  const double r8 = r7 / r;

  D3::Vector mfq = _quadrupole * mfr;
  D3::Vector mfp = _polarizability * mfr;
 
  const double prd = vdot(mfr, _dipole);
  const double prq = vdot(mfr, mfq);
  const double prp = vdot(mfr, mfp);

  double ener_val = -prd * r3 + prq * r5 / 2.  - prp * r6 / 2.;

  if(!torque)
    return ener_val;

  D3::Vector& force = *torque;
  torque += 1;

  torque[0] = 0.;

  D3::Vector vtemp;
  // dipole
  force = r3 * _dipole - (3. * prd * r5) * mfr;

  D3::vprod(_dipole, mfr, vtemp);
  vtemp *= r3;
  torque[1] = vtemp;

  //quadrupole
  force += (2.5 * prq * r7) * mfr - r5 * mfq;

  D3::vprod(mfr, mfq, vtemp);
  vtemp *= r5;
  torque[1] += vtemp;

  // polarizability
  force += r6 * mfp - (3. * prp * r8) * mfr;

  D3::vprod(mfp, mfr, vtemp);
  vtemp *= r6;
  torque[1] += vtemp;

  // convert force to laboratory frame
  dc.mf2lf(1, force, vtemp);
  force = vtemp;

  return ener_val;
}

Potential::ChargeNonlinear::ChargeNonlinear (std::istream& from) throw(Error::General)
  : _charge(1.)
{
  static const char funame [] = "Potential::ChargeNonlinear::ChargeNonlinear (std::istream&): ";

  if(Structure::fragment(1).type() != Molecule::NONLINEAR || Structure::fragment(1).top() == Molecule::SPHERICAL) {
    std::cerr << funame << "second fragment should be nonlinear for and non-spherical this potential\n";
    throw Error::Logic();
  }

  IO::Marker funame_marker(funame);

  double dtemp;
  std::vector<double> vtemp;

  KeyGroup ChargeNonlinearPotential;

  Key   chr_key("Charge[au]");
  Key   dip_key("DipoleMoment[au]");
  Key  quad_key("QuadrupoleMoment[au]");
  Key polar_key("Polarizability[au]");

  for(int i = 0; i < 3; ++i)
    _dipole[i] = 0.;
  for(int i = 0; i < 3; ++i)
    for(int j = i; j < 3; ++j) {
      _quadrupole(i, j) = 0.;
      _quadrupole(j, i) = 0.;
      _polarizability(i, j) = 0.;
      _polarizability(j, i) = 0.;
    }

  std::string token, line, comment;
  while(from >> token) {// read cycle
    if(token == IO::end_key()) {
      std::getline(from, comment);
      break;
    }
    else if(token == chr_key) {
      from >> _charge;
      if(!from) {
	std::cerr << funame << token << ": bad input\n";
	throw Error::Init();
      }
      std::getline(from, comment);
    }
    else if(token == dip_key) { // dipole
      std::getline(from, line);
      std::istringstream iss(line);
      vtemp.clear();
      while(iss >> dtemp)
	vtemp.push_back(dtemp);

      switch(Structure::fragment(1).top()) {
      case Molecule::ASYM:
	if(vtemp.size() == 3) {
	  for(int i = 0; i < 3; ++i)
	    _dipole[i] = vtemp[i];
	}
	else {
	  std::cerr << funame << "dipole moment of non-spherical fragment should have three components\n";
	  throw Error::Input();
	}
	break;
      case Molecule::PROLATE:
	if(vtemp.size() == 1) {
	  _dipole[0] = vtemp[0];
	}
	else if(vtemp.size() == 3) {
	  std::cerr << funame << "WARNING: there should be only one component for the dipole moment "
		    << "along the symmetry axis, ignoring others\n";
	  _dipole[0] = vtemp[0];
	}
	else {
	  std::cerr << funame << "dipole moment should have one component along the symmetry axis\n";
	  throw Error::Input();
	}
	break;
      case Molecule::OBLATE:
	if(vtemp.size() == 1) {
	  _dipole[2] = vtemp[0];
	}
	else if(vtemp.size() == 3) {
	  std::cerr << funame << "WARNING: there should be only one component for the dipole moment "
		    << "along the symmetry axis, ignoring others\n";
	  _dipole[2] = vtemp[2];
	}
	else {
	  std::cerr << funame << "dipole moment should have one component along the symmetry axis\n";
	  throw Error::Input();
	}
      }
    }
    else if(token == quad_key) { // quadrupole
      std::getline(from, line);
      std::istringstream iss(line);
      vtemp.clear();
      while(iss >> dtemp)
	vtemp.push_back(dtemp);
      int k = 0;
      switch(Structure::fragment(1).top()) {
      case Molecule::ASYM:
	if(vtemp.size() == 5) {
	for(int i = 0; i < 2; ++i)
	  for(int j = i; j < 3; ++j) {
	    _quadrupole(i, j) = vtemp[k++];
	    if(i != j)
	      _quadrupole(j, i) = _quadrupole(i, j);
	  }
	}
	else {
	  std::cerr << funame << "quadrupole moment of non-spherical top should have five components: qxx, qxy, qxz, qyy, qyz\n";
	  throw Error::Input();
	}
	break;
      case Molecule::PROLATE:
	if(vtemp.size() == 1) {
	  _quadrupole(0, 0) = vtemp[0];
	  _quadrupole(1, 1) = -_quadrupole(0, 0) / 2.;
	  _quadrupole(2, 2) = -_quadrupole(0, 0) / 2.;
	}
	else {
	  std::cerr << funame << "quadrupole moment of symmetric top should have one component along the symmetry axis\n";
	  throw Error::Input();
	}
	break;
      case Molecule::OBLATE:
	if(vtemp.size() == 1) {
	  _quadrupole(2, 2) = vtemp[0];
	  _quadrupole(0, 0) = -_quadrupole(2, 2) / 2.;
	  _quadrupole(1, 1) = -_quadrupole(2, 2) / 2.;
	}
	else {
	  std::cerr << funame << "quadrupole moment of symmetric top should have one component along the symmetry axis\n";
	  throw Error::Input();
	}
	break;
      default:
	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
    }
    else if(token == polar_key) { //  polarizability
      std::getline(from, line);
      std::istringstream iss(line);
      vtemp.clear();
      while(iss >> dtemp)
	vtemp.push_back(dtemp);
      int k = 0;

      switch(Structure::fragment(1).top()) {
      case Molecule::ASYM:
	if(vtemp.size() == 6) {
	for(int i = 0; i < 3; ++i)
	  for(int j = i; j < 3; ++j) {
	    _polarizability(i, j) = vtemp[k++];
	    if(i != j)
	      _polarizability(j, i) = _polarizability(i, j);
	  }
	}
	else {
	  std::cerr << funame << "polarizability of non-spherical top should have six components: pxx, pxy, pxz, pyy, pyz, pzz\n";
	  throw Error::Input();
	}
	break;
      case Molecule::PROLATE:
	if(vtemp.size() == 3) {
	  for(int i = 0; i < 3; ++i)
	    _polarizability(i, i) = vtemp[i];
	  if(_polarizability(2, 2) != _polarizability(1, 1)) {
	    std::cerr << funame << "p_zz should be equal to p_yy\n";
	    throw Error::Input();
	  }
	}
	else {
	  std::cerr << funame << "polarizability of prolate symmetric top should have three components: p_xx, p_yy, p_zz = p_yy\n";
	  throw Error::Input();
	}
	break;
      case Molecule::OBLATE:
	if(vtemp.size() == 3) {
	  for(int i = 0; i < 3; ++i)
	    _polarizability(i, i) = vtemp[i];
	  if(_polarizability(0, 0) != _polarizability(1, 1)) {
	    std::cerr << funame << "p_xx should be equal to p_yy\n";
	    throw Error::Input();
	  }
	}
	else {
	  std::cerr << funame << "polarizability of oblate symmetric top should have three components: p_xx, p_yy = p_xx, p_zz\n";
	  throw Error::Input();
	}
	break;
      default:
	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
    }
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Init();
    }
  }// read cycle

  // print out
  IO::log << "Charge-Nonlinear potential:\n";
  IO::log << "  Charge = " << _charge << "\n";
  IO::log << std::setw(15) << "DX"
	    << std::setw(15) << "DY"
	    << std::setw(15) << "DZ" 
	    << "\n";
  for(int i = 0; i < 3; ++i)
    IO::log << std::setw(15) << _dipole[i];
  IO::log << "\n\n";
  IO::log << std::setw(15) << "QXX"
	    << std::setw(15) << "QXY"
	    << std::setw(15) << "QXZ"
	    << std::setw(15) << "QYY"
	    << std::setw(15) << "QYZ"
	    << std::setw(15) << "QZZ"
	    << "\n";
  for(int i = 0; i < 3; ++i)
    for(int j = i; j < 3; ++j)
      IO::log << std::setw(15) << _quadrupole(i, j);
  IO::log << "\n\n";
  IO::log << std::setw(15) << "PXX"
	    << std::setw(15) << "PXY"
	    << std::setw(15) << "PXZ"
	    << std::setw(15) << "PYY"
	    << std::setw(15) << "PYZ"
	    << std::setw(15) << "PZZ"
	    << "\n";
  for(int i = 0; i < 3; ++i)
    for(int j = i; j < 3; ++j)
      IO::log << std::setw(15) << _polarizability(i, j);
  IO::log << "\n\n";

  // renormalize multipole moments
  for(int i = 0; i < 3; ++i)
    _dipole[i] *= _charge;

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      _quadrupole    (i, j) *= _charge;
      _polarizability(i, j) *= _charge * _charge;
    }
}

  /*************************************************************************
   ******************************** Dipole-Dipole **************************
   *************************************************************************/

  double Potential::DipoleDipole::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const throw(Error::General)
  {
    static const char funame [] = "Potential::DipoleDipole::operator(): ";
  
    double dtemp;
    D3::Vector vtemp;

    D3::Vector lfr(dc.orb_pos());
    const double r = lfr.vlength();
    const double r2  = 1. / r / r ;
    const double r3  = r2 / r;
    const double r5  = r3 * r2;
    const double r6  = r5 / r;
    const double r7  = r6 / r;
    const double r8  = r7 / r;
    const double r10 = r8 * r2;

    D3::Vector lfd [2]; // dipole moments in the laboratory frame
    for(int frag = 0; frag < 2; ++frag)
      switch(Structure::fragment(frag).type()) {
      case Molecule::LINEAR:
	lfd[frag] = dc.ang_pos(frag);
	lfd[frag].normalize();
	lfd[frag] *= _dipole[frag][0];
	break;
      case Molecule::NONLINEAR:
	dc.mf2lf(frag, _dipole[frag], lfd[frag]);
	break;
      default:
	std::cerr << funame << "you are in trouble. Ha, ha, ha ...\n";
	throw Error::Logic();
      }
  
    int other;     // the other fragment
    double dr [2]; // d * r scalar product
    for(int frag = 0; frag < 2; ++frag)
      dr[frag] = vdot(lfr, lfd[frag]);
   
    double dd = vdot(lfd[0], lfd[1]);
    double ener_val =  dd * r3 - 3. * dr[0] * dr[1] * r5 // dipole - dipole         term
      - _dispersion * r6;                      // dispersion              term

    for(int frag = 0; frag < 2; ++frag) {      // dipole - induced dipole term
      other = 1 - frag;
      ener_val -= 1.5 * _polarizability[frag] * dr[other] * dr[other] * r8;
    }

    if(!torque)
      return ener_val;

    // gradient calculation
    D3::Vector& force = *torque;
    torque += 1; 

    D3::Vector ddv;                      // dipole moments vector product d1 x d2
    D3::vprod(lfd[0], lfd[1], ddv);
    D3::Vector drv [2];                  // d x r vector product
    for(int frag = 0; frag < 2; ++frag)
      D3::vprod(lfd[frag], lfr, drv[frag]);

    // force
    dtemp = -6. * _dispersion * r8               // dispersion              term
      + 3. * dd * r5 - 15. * dr[0] * dr[1] * r7; // dipole-dipole           term
    for(int frag = 0; frag < 2; ++frag) {
      other = 1 - frag;
      dtemp -= 12. * dr[frag] * dr[frag]         // dipole - induced dipole terms 
	* _polarizability[other] * r10;
    }
    force = dtemp * lfr;

    for(int frag = 0; frag < 2; ++frag) {
      other = 1 - frag;
      dtemp = 3. * dr[other] * r5 +                  // dipole - dipole term
	3. * _polarizability[other] * dr[frag] * r8; // dipole - induced dipole term
      force += dtemp * lfd[frag];
    }

    // ... and torque
    for(int frag = 0; frag < 2; ++frag) {
      other = 1 - frag;
      torque[frag] = ddv;                                // dipole-dipole term
      if(frag)
	torque[frag] *= r3;
      else
	torque[frag] *= -r3;
      
      dtemp = 3. * dr[other] * r5                       // dipole - dipole term
	+ 3. * _polarizability[other] * dr[frag] * r8;  // dipole - induced dipole term
      torque[frag] += dtemp * drv[frag];
    }  

    // for nonlinear fragments convert torques to their molecular frames
    for(int frag = 0; frag < 2; ++frag)
      if(Structure::fragment(frag).type() == Molecule::NONLINEAR) {
	dc.lf2mf(frag, torque[frag], vtemp);
	torque[frag] = vtemp;
      }

    return ener_val;
  }

Potential::DipoleDipole::DipoleDipole (std::istream& from) throw(Error::General)
{
  static const char funame [] = "Potential::DipoleDipole::DipoleDipole: ";

  for(int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC || Structure::fragment(frag).top() == Molecule::SPHERICAL) {
      std::cerr << funame << "this potential should not be used for monoatomic or spherical fragments\n";
      throw Error::Logic();
    }

  IO::Marker funame_marker(funame);

  double dtemp;
  std::vector<double> vtemp;

  KeyGroup DipoleDipolePotential;

  Key dip_key("DipoleMoment[au]");
  Key pol_key("Polarizability[au]");
  Key dsp_key("Dispersion[au]");

  // initialization
  for(int frag = 0; frag < 2; ++frag) {
    _dipole        [frag] = 0.;
    _polarizability[frag] = 0.;
  }
  _dispersion = 0.;

  int dip_count = 0, pol_count = 0;
  std::string token, line, comment;
  while(from >> token) {// read cycle
    if(token == IO::end_key()) {
      std::getline(from, comment);
      break;
    }
    else if(token == dip_key) {
      if(dip_count >= 2) {
	std::cerr << funame << token << ": too many dipoles ... exiting\n";
	throw Error::Input();
      }

      std::getline(from, line);
      std::istringstream iss(line);
      vtemp.clear();
      while(iss >> dtemp)
	vtemp.push_back(dtemp);
      
      switch(Structure::fragment(dip_count).type()) {
      case Molecule::LINEAR:
	if(vtemp.size() == 1) {
	  _dipole[dip_count][0] = vtemp[0];
	}
	else if(vtemp.size() == 3) {
	  std::cerr << funame << "WARNING: there should be only one component for the dipole moment of "
		    << dip_count << "-th linear fragment, ignoring others\n";
	  _dipole[dip_count][0] = vtemp[0];
	}
	else {
	  std::cerr << funame << token << ": dipole moment is unreadable\n";
	  throw Error::Input();
	}
	break;
      case Molecule::NONLINEAR:
	switch(Structure::fragment(dip_count).top()) {
	case Molecule::ASYM:
	  if(vtemp.size() == 3) {
	    for(int i = 0; i < 3; ++i)
	      _dipole[dip_count][i] = vtemp[i];
	  }
	  else {
	    std::cerr << funame << "dipole moment should have three components\n";
	    throw Error::Input();
	  }
	  break;
	case Molecule::PROLATE:
	  if(vtemp.size() == 1) {
	    _dipole[dip_count][0] = vtemp[0];
	  }
	  else if(vtemp.size() == 3) {
	    std::cerr << funame << "WARNING: there should be only one component for the dipole moment of "
		      << dip_count << "-th fragment along the symmetry axis, ignoring others\n";
	    _dipole[dip_count][0] = vtemp[0];
	  }
	  else {
	    std::cerr << funame << "dipole moment should have one component along the symmetry axis\n";
	    throw Error::Input();
	  }
	  break;
	case Molecule::OBLATE:
	  if(vtemp.size() == 1) {
	    _dipole[dip_count][2] = vtemp[0];
	  }
	  else if(vtemp.size() == 3) {
	    std::cerr << funame << "WARNING: there should be only one component for the dipole moment of "
		      << dip_count << "-th fragment along the symmetry axis, ignoring others\n";
	    _dipole[dip_count][2] = vtemp[2];
	  }
	  else {
	    std::cerr << funame << "dipole moment should have one component along the symmetry axis\n";
	    throw Error::Input();
	  }
	  break;
	default:
	  std::cerr << funame << "should not be here\n";
	  throw Error::Logic();
	}
	break;
      default:
	std::cerr << funame << token << ": should not be here\n";
	throw Error::Logic();
      }
      dip_count++;
    }
    else if(token == pol_key) {
      if(pol_count >= 2) {
	std::cerr << funame << token << ": too many polarizabilities ... exiting\n";
	throw Error::Input();
      }
      from >> _polarizability[pol_count];
      if(!from) {
	std::cerr << funame << token << ": input stream is corrupted\n";
	throw Error::Input();
      }
      pol_count++;
    }
    else if(token == dsp_key) {
      from >> _dispersion;
      if(!from) {
	std::cerr << funame << token << ": input stream is corrupted\n";
	throw Error::Input();
      }
    }
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Input();
    }
  }// read cycle

  if(dip_count < 2)
    std::cerr << funame << "WARNING: not all dipole have been initialized, assuming zeros\n";

  if(pol_count < 2)
    std::cerr << funame << "WARNING: not all polarizabilities have been initialized, assuming zeros\n";

  IO::log << "dipole-dipole potential:\n";
  IO::log << "  Dispersion = " << _dispersion << "\n\n";
  IO::log << std::setw(10) << "Fragment"
	    << std::setw(15) << "DX"
	    << std::setw(15) << "DY"
	    << std::setw(15) << "DZ" 
	    << std::setw(20) << "Polarizability"
	    << "\n";

  for(int frag = 0; frag < 2; ++frag) {
    IO::log << std::setw(10) << frag;
    for(int i = 0; i < 3; ++i)
      IO::log << std::setw(15) << _dipole[frag][i];
    IO::log << std::setw(20) << _polarizability[frag];
    IO::log << "\n";
  }

  // renormalizing dispersion term with dipole-induced dipole contribution

  double dd; // dipole moment squared
  int other; // the other fragment
  for(int frag = 0; frag < 2; ++frag) {
    other = 1 - frag;
    switch(Structure::fragment(frag).type()) {
    case Molecule::LINEAR:
      dd = _dipole[frag][0] * _dipole[frag][0];
      break;
    case Molecule::NONLINEAR:
      dd = vdot(_dipole[frag]);
      break;
    default:
      std::cerr << funame << ": should not be here\n";
      throw Error::Logic();
    }
    _dispersion += dd * _polarizability[other] / 2.;
  } 

}

/*************************************************************************
 ******************************** Multipole ******************************
 *************************************************************************/

double Potential::Multipole::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const throw(Error::General)
{
  static const char funame [] = "Potential::DipoleDipole::operator(): ";
  
  // initialization
  double ener_val = 0.;
  if(torque)
    for(int i = 0; i < 3; ++i)
      torque[i] = 0.;

  double dtemp;
  D3::Vector vtemp;
  int itemp, frag, other;

  const D3::Vector r(dc.orb_pos());

  const double dist = vlength(dc.orb_pos(), 3);
  const double r2   = 1. / dist / dist ;
  const double r3   = r2 / dist;
  const double r5   = r3 * r2;
  const double r6   = r5 / dist;
  const double r7   = r6 / dist;
  const double r8   = r7 / dist;
  const double r9   = r8 / dist;
  const double r10 = r8 * r2;
  const double r12 = r10 * r2;

  for(frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).charge())
      break;
    
  if(frag < 2) { // charge-multipole interaction
    other = 1 - frag;

    // charge
    const double charge = (double)Structure::fragment(frag).charge();

    D3::Vector d, qr, pr ;
    double rd, rqr, rpr;

    // dipole
    if(Structure::fragment(other).dipole().size()) {

      dc.dipole(other, d);
      d *= charge;
      rd = vdot(r, d, 3);
      if(frag)
	ener_val += rd * r3;
      else
	ener_val -= rd * r3;
    }

    // quadrupole
    if(Structure::fragment(other).quadrupole().size()) {

      dc.quadrupole_vector_product(other, r, qr);
      qr *= charge;
      rqr = vdot(r, qr, 3);
      ener_val += rqr * r5 / 2.; 

    }

    // polarizability
    if(Structure::fragment(other).polarizability().size()) {

      dc.polarizability_vector_product(other, r, pr);
      pr *= charge * charge;
      rpr = vdot(r, pr, 3);
      ener_val -= rpr * r6 / 2.;

    }

    if(!torque)
      return ener_val;

    D3::Vector& force = *torque;
    torque += 1;
      
    // dipole
    if(Structure::fragment(other).dipole().size()) {

      // force
      if(frag)
	force -= r3 * d  - (3. * r5 * rd) * r;
      else
	force += r3 * d  - (3. * r5 * rd) * r;

      // torque
      D3::vprod(d, r, vtemp);
      vtemp *= r3;
      if(frag)
	torque[other] -= vtemp;
      else
	torque[other] += vtemp;
    }

    // quadrupole
    if(Structure::fragment(other).quadrupole().size()) {

      // force
      force += (2.5 * rqr * r7) * r - r5 * qr;

      // torque
      D3::vprod(r, qr, vtemp);
      vtemp *= r5;
      torque[other] += vtemp;

    }

    // polarizability
    if(Structure::fragment(other).polarizability().size()) {

      // force
      force += r6 * pr - (3. * r8 * rpr) * r;

      // torque
      D3::vprod(pr, r, vtemp);
      vtemp *= r6;
      torque[other] += vtemp;
	
    }

    if(Structure::fragment(other).type() == Molecule::NONLINEAR) {

      dc.lf2mf(other, torque[other], vtemp);
      torque[other] = vtemp;

    }

    return ener_val;

  }// charge-multipole potential
      
  D3::Vector   d [2]; // dipole moments in the laboratory frame
  D3::Vector  qr [2]; // q times r;
  D3::Vector  pr [2]; // p times r;
  D3::Vector  qd [2]; // q times d;
  D3::Vector  pd [2]; // p times d;
 
  double rd[2], rqr[2], rpr[2], rqd[2], rpd[2], dpd[2];

  for(frag = 0; frag < 2; ++frag) {
    if(Structure::fragment(frag).dipole().size()) {
      dc.dipole(frag, d[frag]);
      other = 1 - frag;
      rd[frag] = vdot(r, d[frag], 3);

      if(Structure::fragment(other).quadrupole().size()) { // dipole-qudrupole interaction
	 dc.quadrupole_vector_product(other, r,         qr[frag]);
	 dc.quadrupole_vector_product(other, d[frag], qd[frag]);
	 rqr[frag] = vdot(r, qr[frag], 3);
	 rqd[frag] = vdot(r, qd[frag], 3);

	 dtemp = 2.5 * rqr[frag] * rd[frag] * r7 - rqd[frag] * r5;
	 if(frag) 
	   ener_val -= dtemp;
	 else 
	   ener_val += dtemp;
      }// dipole-quadrupole interaction

      if(Structure::fragment(other).polarizability().size()) {// dipole - induced dipole interaction
	dc.polarizability_vector_product(other, r,         pr[frag]);
	dc.polarizability_vector_product(other, d[frag], pd[frag]);
	rpr[frag] = vdot(r, pr[frag], 3);
	rpd[frag] = vdot(r, pd[frag], 3);
	dpd[frag] = vdot(d[frag], pd[frag], 3);

	ener_val -=  dpd[frag] * r6 / 2. + 4.5 * rpr[frag] * rd[frag] * rd[frag] * r10 - 3. * rpd[frag] * rd[frag] * r8;
      }
    }// dipole - induced dipole interaction
  }

  double dd;

  if(Structure::fragment(0).dipole().size() && Structure::fragment(1).dipole().size()) {// dipole-dipole interaction

    dd = vdot(d[0], d[1], 3);
    ener_val += dd * r3 - 3. * rd[0] * rd[1] * r5;

  }// dipole-dipole interaction

  // dispersion
  if(_dispersion != 0.)
    ener_val -= _dispersion * r6;

  if(!torque)
    return ener_val;

  D3::Vector& force = *torque;
  torque += 1;

  D3::Vector rdv[2];
  for(frag = 0; frag < 2; ++frag) {

    if(Structure::fragment(frag).dipole().size()) {
      other = 1 - frag;
      D3::vprod(r, d[frag], rdv[frag]);
      if(Structure::fragment(other).quadrupole().size()) { // dipole-qudrupole interaction

	// force
	vtemp  = r5 * qd[frag]  - (5. *  rd[frag] * r7) * qr[frag] - (2.5 * rqr[frag] * r7) * d[frag] 
	  + (17.5 * rqr[frag] * rd[frag] * r9 - 5. * rqd[frag] * r7) * r;
	if(frag) 
	  force -= vtemp;
	else 
	  force += vtemp;
	
	// torque on dipole
	D3::Vector dqrv, rqrv, rqdv;
	D3::vprod(d[frag], qr[frag], dqrv); 
	D3::vprod(r,         qr[frag], rqrv); 
	D3::vprod(r,         qd[frag], rqdv); 

	vtemp = (2.5 * rqr[frag] * r7) * rdv[frag] + r5 * dqrv;
	if(frag)
	  torque[frag] -= vtemp;
	else
	  torque[frag] += vtemp;

	vtemp = (5. * rd[frag] * r7) * rqrv - r5 * rqdv - r5 * dqrv;
	if(frag)
	  torque[other] -= vtemp;
	else
	  torque[other] += vtemp;

      }// dipole-quadrupole interaction

      if(Structure::fragment(other).polarizability().size()) {// dipole - induced dipole interaction
	force += (9. * rpr[frag] * rd[frag] * r10 - 3. * rpd[frag] * r8) * d[frag]
	  - (3. * rd[frag] * r8) * pd[frag]
	  +(9. * rd[frag] * rd[frag] * r10) * pr[frag]
	  +(24. * rpd[frag] * rd[frag] * r10 - 45. * rpr[frag] * rd[frag] * rd[frag] * r12 - 3. * dpd[frag] * r8) * r;

	D3::Vector dpdv, dprv, rpdv, rprv;
	D3::vprod(d[frag], pd[frag], dpdv);
	D3::vprod(d[frag], pr[frag], dprv);
	D3::vprod(r,         pd[frag], rpdv);
	D3::vprod(r,         pr[frag], rprv);
	
	torque[frag] += r6 * dpdv 
	  + (3. * rpd[frag] * r8 - 9. * rpr[frag] * rd[frag] * r10) * rdv[frag]
	  - (3. * rd[frag] * r8) * dprv; 

	torque[other] += - r6 * dpdv
	  - (9. * rd[frag] * rd[frag] * r10) * rprv
	  + (3. * rd[frag] * r8) * rpdv
	  + (3. * rd[frag] * r8) * dprv;
      }
    }// dipole - induced dipole interaction
  }

  if(Structure::fragment(0).dipole().size() && Structure::fragment(1).dipole().size()) {// dipole-dipole interaction

    D3::Vector ddv;
    D3::vprod(d[0], d[1], ddv);

    force += (3. * dd * r5 - 15. * rd[0] * rd[1] * r7) * r;
    for( frag = 0; frag < 2; ++frag) {
      other = 1 - frag;
      force += (3. * rd[other] * r5) * d[frag];
      torque[frag] -= (3. * rd[other] * r5) * rdv[frag];
      if(frag)
	torque[frag] += r3 * ddv;
      else
	torque[frag] -= r3 * ddv;
    }

  }// dipole-dipole interaction

    
  // dispersion
  if(_dispersion != 0.)
    force -= (6. * r8 * _dispersion) * r;  

  for(frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() == Molecule::NONLINEAR) {
      dc.lf2mf(frag, torque[frag], vtemp);
      torque[frag] = vtemp;
    }
   
  return ener_val;
}

Potential::Multipole::Multipole (std::istream& from) throw(Error::General) : _dispersion(0.)
{
  static const char funame [] = "Potential::Multipole::Multipole: ";

  if(Structure::fragment(0).charge() &&  Structure::fragment(1).charge()) {
    std::cerr << funame << "ERROR: two charged fragments have  not been implemented\n";
    throw Error::Logic();
  }

  IO::Marker funame_marker(funame);

  KeyGroup MultipolePotential;
  Key dsp_key("Dispersion[au]");

  std::string token, line, comment;
  while(from >> token) {// read cycle
    if(token == IO::end_key()) {
      std::getline(from, comment);
      break;
    }
    else if(token == dsp_key) {
      from >> _dispersion;
      if(!from) {
	std::cerr << funame << token << ": input stream is corrupted\n";
	throw Error::Input();
      }
    }
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Input();
    }
  }// read cycle

  if(_dispersion != 0.) {
    IO::log << "Multipole potential:\n";
    IO::log << "   Dispersion = " << _dispersion << "\n\n";
  }
}
