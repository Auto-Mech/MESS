#include "lr.hh"
#include "multindex.hh"
#include "read.hh"
#include "slatec.hh"

#include<map>
#include<vector>
#include<cmath>

namespace LongRange {

  bool _global_init = false;
  bool isinit () { return _global_init; }

  sym_t _sym_type [2]; // fragment symmetry 
  int   _eff_dof;      // effective nummber of degrees of freedom
  int   _orient_dim;   // number of angles the potential effectively depends on

  // symmetric fragments data
  std::vector<int> symmetric_fragment;
  std::vector<double> symmetry_axis_imom;

  // potential
  Potential::Wrap pot;

  // angular integration grid
  int _angl_grid_size;
  int  angular_grid_size () { return _angl_grid_size; }
  void set_angular_grid (int g) { _angl_grid_size = g; }


  // correspondence between <fragment index, angle type> and angle index
  enum ang_t { THETA, PHI, PSI };

  std::map<std::pair<int, ang_t>, int> index_map;
  typedef std::map<std::pair<int, ang_t>, int>::const_iterator IndexMapIterator;
  int  euler_angle_index (int frag, ang_t angle) throw(Error::General);
  void ang2dc      (const std::vector<double>&, Dynamic::Coordinates&) throw(Error::General);

  // number of states densities stuff
  int    EDensity::_dim;
  double EDensity::_nfac;

  int    JDensity::_dim;
  double JDensity::_nfac;

  int    MDensity::_dim;
  double MDensity::_nfac;

  int    KDensity::_dim;
  double KDensity::_nfac;

  // integration steps
  double JIntegral::step;
  double MIntegral::step [2];
  double KIntegral::step;

  // centrifugal barrier stuff
  double _amom_min, _amom_max;
  Slatec::Spline _barrier_log;

  //barrier tolerances
  double bar_ener_tol;
  double bar_dist_tol    = 0.01;      // logarithmic distance tolerance
  double bar_search_step = 0.1;

  class BarrierSearch : public Math::MinimumSearch {
  public:

    double amom;
    double _fun (double x) throw(Math::MinimumSearch::Exception);
    
    BarrierSearch () : Math::MinimumSearch(bar_dist_tol, bar_ener_tol) {}
  };

  // some math
  void increment(double&, double);

  // tolerances for number of states optimizaiton
  double Opt::xtol; // logarithmic tolerance for distance
  double Opt::ytol; // logarithmic tolerance for number of states
  double Opt::step; // initial step
}

void LongRange::increment (double& res, double val) {
  if(res > 0.)
    res += val;
  else
    res  = val;
}

double LongRange::j_shift (double j, double d) throw(Error::General)
{ 
  static const char funame [] = "LongRange::j_shift: ";
  
  if(d <= 0.) {
    std::cerr << funame << "negative or zero distance\n";
    throw Error::Range();
  }

  return j * j / 2. / Structure::mass() / d / d;
}
 
double LongRange::m_shift (const std::vector<double>& m) throw(Error::General)
{
  static const char funame [] = "LongRange::m_shift: ";

  if(m.size() != symmetric_fragment_size()) {
    std::cerr << funame << "wrong number of symmetric fragments\n";
    throw Error::Init();
  }
  double res = 0.;

  for(int i = 0; i < m.size(); ++i)
    res += m[i] * m[i] / symmetry_axis_imom[i];

  return res / 2.;
}

double LongRange::BarrierSearch::_fun (double x) throw(Math::MinimumSearch::Exception)
{
  double d = std::exp(x);
  return - minimum_energy(d) - j_shift(amom, d);
}

double LongRange::centrifugal_barrier (double amom) {
  
  static const char funame [] = "LongRange::centrifugal_barrier: ";

  if(amom > _amom_max) {
    std::cerr << funame << "WARNING: angular momentum = " << amom << " is too big, truncating\n";
    amom = _amom_max;
  }

  if(amom < _amom_min) {
    std::cerr << funame << "WARNING: angular momentum = " << amom << " is too small, neglecting\n";
    return 0.;
  }

  return std::exp(_barrier_log(std::log(amom)));
}

int LongRange::symmetric_fragment_size() throw(Error::General)
{
  static const char funame [] = "LongRange::symmetric_fragment_size: ";

#ifdef DEBUG
  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
#endif

  return symmetric_fragment.size();
}

LongRange::sym_t LongRange::symmetry_type (int frag) throw(Error::General)
{
  static const char funame [] = "LongRange::symmetry_type: ";

#ifdef DEBUG

  if(frag < 0 || frag >= 2) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

#endif

  return _sym_type[frag];

}

int LongRange::eff_dof () throw(Error::General)
{
  static const char funame [] = "LongRange::eff_dof: ";

#ifdef DEBUG
  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }
#endif

  return _eff_dof;
}

int LongRange::orientational_dimension () throw(Error::General)
{
  static const char funame [] = "LongRange::orientational_dimension: ";

#ifdef DEBUG

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

#endif

  return _orient_dim;
}

void LongRange::init (std::istream& from) throw(Error::General)
{
  static const char funame [] = "LongRange::init: ";

  static const double one_eps = 1. + 1.e-10;
  
  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not been initialized\n";
    throw Error::Init();
  }

  _global_init = true;

  int itemp;
  bool btemp;
  double dtemp;

  int amom_size;
  double bar_dist;

  std::map<std::string, Read> input;
  std::map<std::string, Read>::iterator idit;

  input["AngularGridSize"                 ] = Read(_angl_grid_size,    10);
  input["J-IntegralStep[au]"              ] = Read(JIntegral::step,    1.);
  input["M-IntegralStep[au]"              ] = Read(MIntegral::step[0], 1.);
  input["K-IntegralStep[au]"              ] = Read(KIntegral::step,    1.);
  input["AngularMomentumMinimum[au]"      ] = Read(_amom_min,       0.001);
  input["AngularMomentumMaximum[au]"      ] = Read(_amom_max,       1000.);
  input["AngularMomentumSize"             ] = Read(amom_size,        1000);
  input["BarrierEnergyTolerance[kcal/mol]"] = Read(bar_ener_tol,    1.e-4);
  input["BarrierLogDistanceTolerance"     ] = Read(bar_dist_tol,    1.e-2);
  input["BarrierSearchStep"               ] = Read(bar_search_step, 1.e-1);
  input["MaximumBarrierDistanceGuess[au]" ] = Read(bar_dist,           3.);  
  input["OptArgumentTolerance"            ] = Read(Opt::xtol,       1.e-2);
  input["OptFunctionTolerance"            ] = Read(Opt::ytol,       1.e-2);
  input["OptStep"                         ] = Read(Opt::step,          1.);
  
  std::string key, comment;
  while(from >> key) {
    if(key == IO::end_key())
      break;
    idit = input.find(key);
    if (idit == input.end()) {
      std::cerr << funame << "WARNING: did not find the key " << key
		<< ", taking it as a comment\n";
      std::getline(from, comment);
    }
    else
      from >> idit->second;
  }

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Form();
  }

  // check if all parameters were initialized
  btemp = true;
  for(idit = input.begin(); idit != input.end(); ++idit)
    if(!idit->second.is_init()) {
      std::cerr << funame << idit->first << " is not initialized\n";
      btemp = false;
    }
  if(!btemp)
    throw Error::Init();

  // default values
  std::cout << funame << "Default parameters:\n" << std::left;
  for(idit = input.begin(); idit != input.end(); ++idit)
    if(idit->second.is_default())
      std::cout << "   " << std::setw(20) << idit->first << " = " << idit->second << "\n";
  std::cout << std::right << "\n";

  // default M-integral step
  MIntegral::step[1] = MIntegral::step[0];

  // changing units
  bar_ener_tol *= Phys_const::kcal;

  // molecular potential symmetry type
  for(int frag = 0; frag < 2; ++frag) {
    if(pot.type() == Potential::MULTIPOLE && Structure::fragment(frag).charge())
      _sym_type[frag] = SPHERICAL;
    else
      switch(Structure::type(frag)) {
      case Molecule::MONOATOMIC:

	_sym_type[frag] = SPHERICAL;
	break;

      case Molecule::LINEAR:

	_sym_type[frag] = LINEAR;
	break;

      case Molecule::NONLINEAR:

	switch(Structure::top(frag)) {
	case Molecule::SPHERICAL:

	  _sym_type[frag] = SPHERICAL;
	  break;

	case Molecule::ASYM:

	  _sym_type[frag] = GENERIC;
	  break;

	case Molecule::PROLATE:

	  _sym_type[frag] = SYMMETRIC;
	  symmetric_fragment.push_back(frag);
	  symmetry_axis_imom.push_back(Structure::fragment(frag).imom(0));
	  break;

	case Molecule::OBLATE:

	  _sym_type[frag] = SYMMETRIC;
	  symmetric_fragment.push_back(frag);
	  symmetry_axis_imom.push_back(Structure::fragment(frag).imom(2));
	  break;
	}
	break;

      default:

	std::cerr << funame << "wrong molecular type\n";
	throw Error::Logic();
      }
  }
  // orientational dimension and effective number of degrees of freedom
  _orient_dim = -1;
  _eff_dof     = 3;
  for( int frag = 0;  frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:

      break;

    case LINEAR:

      _orient_dim += 2;
      _eff_dof    += 2;
      break;

    case SYMMETRIC:

      _orient_dim += 2;
      _eff_dof    += 3;
      break;

    case GENERIC:

      _orient_dim += 3;
      _eff_dof    += 3;
      break;

    default:

      std::cerr << funame << "wrong symmetry type\n";
      throw Error::Logic();
    }
  if(orientational_dimension() < 0) {
    std::cerr << funame << "the case of two spherically symmetric fragments is not implemented\n";
    throw Error::Init();
  }

  // correspondence between <fragment index and angle type> and angle index
  int ang_index = 0;
  itemp = 0; // number of spherical fragments
  for(int frag = 0; frag < 2; ++frag) {
    switch(symmetry_type(frag)) {
    case SPHERICAL:
      
      ++itemp;
      break;

    case LINEAR:

      index_map[std::pair<int, ang_t>(frag, THETA)] = ang_index++;
      break;

    case SYMMETRIC:
      
      index_map[std::pair<int, ang_t>(frag, THETA)] = ang_index++;
      break;

    case GENERIC:
      
      index_map[std::pair<int, ang_t>(frag, THETA)] = ang_index++;
      index_map[std::pair<int, ang_t>(frag,   PSI)] = ang_index++;
      break;
    }
  }

  if(!itemp)
    index_map[std::pair<int, ang_t>(1, PHI)] = ang_index++;

  //number of states densities initialization
  EDensity::init();
  JDensity::init();
  MDensity::init();
  KDensity::init();

  std::cout << "Effective orientational dimension = " << orientational_dimension() << "\n"
	    << "Number of symmetric fragments     = " << symmetric_fragment_size() << "\n" 
	    << "\n";

  // centrifugal barrier

#ifdef DEBUG  
  //  std::cout << "Minimum energy[kcal/mol]:\n" 
  //<< std::setw(15) << "R, bohr"
  //<< std::setw(15) << "Vmin"
  //<< "\n";
#endif


  BarrierSearch bar;
  Array<double> x_data(amom_size);
  Array<double> y_data(amom_size);
  double xguess = std::log(bar_dist);

  double amom = _amom_max;
  double amom_step = std::pow(_amom_max / _amom_min, 1. / double(amom_size - 1));
  
#ifdef DEBUG
  std::cout << "Centrifugal barrier search:\n"
	    << std::setw(15) << "Ang.Mom.[au]"
	    << std::setw(15) << "Distance[au]"
	    << std::setw(15) << "Energy[kcal]"
	    << std::endl;
#endif

  for(int i = 0; i < amom_size; ++i, amom /= amom_step) {

    int ii     = amom_size - i - 1;
    x_data[ii] = std::log(amom);

    bar.amom = amom;
    dtemp = - bar.find(xguess, bar_search_step, &xguess);

#ifdef DEBUG
    std::cout << std::setw(15) << amom
	      << std::setw(15) << std::exp(xguess)
	      << std::setw(15) << dtemp / Phys_const::kcal
	      << std::endl;
#endif

    if(dtemp <= 0.) {
      std::cerr << funame << "negative barrier energy\n";
      throw Error::Logic();
    }

    y_data[ii] = std::log(dtemp);
  }

  _barrier_log.init(x_data, y_data, amom_size);
  _amom_max /= one_eps;
  _amom_min *= one_eps;  

  // output
  std::cout << "Centrifugal barrier:\n";
  std::cout << std::setw(22) << "Angular Momentum[au]"
	    << std::setw(19) << "Energy[kcal/mol]"
	    << "\n";

  amom = 1.;
  while(amom < _amom_max) {
    std::cout << std::setw(22) << amom
	      << std::setw(19) << centrifugal_barrier(amom) / Phys_const::kcal
	      << "\n";
    amom += 1.;
  }
  std::cout << std::endl;
}

void LongRange::EDensity::init () throw(Error::General)
{

  // two degrees of freedom more
  _dim = 2;
  _nfac = 2. * M_PI;

  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:
      
      break;
      
    case LINEAR:
      
      _dim += 2;
      _nfac *= Structure::fragment(frag).imom(1) / 2. / M_PI;
      break;
      
    case SYMMETRIC:
      
      _dim += 3;
      for(int i = 0; i < 3; ++i)
	_nfac *= Structure::fragment(frag).imom_sqrt(i) / std::sqrt(2. * M_PI);
      _nfac *= 2. * M_PI;
      break;

    case GENERIC:

      _dim += 3;
      for(int i = 0; i < 3; ++i)
	_nfac *= Structure::fragment(frag).imom_sqrt(i) / std::sqrt(2. * M_PI);
      break;
    }

  _nfac *= 2. * Structure::mass() / tgamma((double)_dim / 2. + 1.);
}

void LongRange::JDensity::init () throw(Error::General)
{
  _dim = 0;
  _nfac = 2. * M_PI;

  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:
      
      break;
      
    case LINEAR:
      
      _dim += 2;
      _nfac *= Structure::fragment(frag).imom(1) / 2. / M_PI;
      break;
      
    case SYMMETRIC:
      
      _dim += 3;
      for(int i = 0; i < 3; ++i)
	_nfac *= Structure::fragment(frag).imom_sqrt(i) / std::sqrt(2. * M_PI);
      _nfac *= 2. * M_PI;
      break;

    case GENERIC:

      _dim += 3;
      for(int i = 0; i < 3; ++i)
	_nfac *= Structure::fragment(frag).imom_sqrt(i) / std::sqrt(2. * M_PI);
      break;
    }
  
  _nfac /= tgamma((double)_dim / 2. + 1.);
}

void LongRange::MDensity::init () throw(Error::General)
{
  _dim = 0;
  _nfac = 2. * M_PI;
  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:
      
      break;
      
    case LINEAR:
      
      _dim += 2;
      _nfac *= Structure::fragment(frag).imom(1) / 2. / M_PI;
      break;
      
    case SYMMETRIC:
      
      _dim += 2;
      _nfac *= Structure::fragment(frag).imom(1) / 2. / M_PI;
      break;

    case GENERIC:

      _dim += 3;
      for(int i = 0; i < 3; ++i)
	_nfac *= Structure::fragment(frag).imom_sqrt(i) / std::sqrt(2. * M_PI);
      break;
    }

  _nfac /= tgamma((double)_dim / 2. + 1.);
}

void LongRange::KDensity::init () throw(Error::General)
{

  // one degree of freedom less
  _dim = -1;
  _nfac = std::sqrt(2. * M_PI);

  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:
      
      break;
      
    case LINEAR:
      
      _dim += 2;
      _nfac *= Structure::fragment(frag).imom(1) / 2. / M_PI;
      break;
      
    case SYMMETRIC:
      
      _dim += 2;
      _nfac *= Structure::fragment(frag).imom(1) / 2. / M_PI;
      break;

    case GENERIC:

      _dim += 3;
      for(int i = 0; i < 3; ++i)
	_nfac *= Structure::fragment(frag).imom_sqrt(i) / std::sqrt(2. * M_PI);
      break;
    }

  _nfac /= tgamma((double)_dim / 2. + 1.);
}

double LongRange::StatesNumberDensity::_density (const Dynamic::Coordinates& dc, int dim, double* mep) const
{
  double poten = pot(dc);
  if(mep)
    *mep = poten;

  double e = energy() - poten;

  if(e <= 0.)
    return -1.;
  
  return std::pow(e, (double)dim / 2.);
}

void LongRange::KDensity::set_m_proj (int i, double m) 
{ 
  _m_proj[symmetric_fragment[i]] = m; 
}

LongRange::KDensity::KDensity (double e, const std::vector<double>& m, double k) throw(Error::General)
  : StatesNumberDensity(e), _k_proj(k)
{
  static const char funame [] = "LongRange::KDensity::KDensity: ";

  if(m.size() != symmetric_fragment_size()) {
    std::cerr << funame << "wrong number of symmetric fragments\n";
    throw Error::Init();
  }

  for(int i = 0; i < m.size(); ++i)
    _m_proj[symmetric_fragment[i]] = m[i];
}


double LongRange::KDensity::operator() (const Dynamic::Coordinates& dc, double* mep) const
{
  static const char funame [] = "LongRange::KDensity::operator(): ";

  int itemp;
  double dtemp;

  double poten = pot(dc);
  if(mep)
    *mep = poten;

  double e = energy() - poten;

  double imz = 0.;     // inertia moment projection on the interfragment axis
  double k = _k_proj;

  D3::Vector lfn(dc.orb_pos());
  lfn.normalize();
  D3::Vector mfn;

  for(int frag = 0; frag < 2; ++frag) {// frag cycle
    switch(symmetry_type(frag)) {// type
    case SPHERICAL:

      break;

    case LINEAR:

      mfn = dc.ang_pos(frag);
      mfn.normalize();
      dtemp = vdot(lfn, mfn);
      imz += Structure::fragment(frag).imom(1) * (1. - dtemp * dtemp);
      break;

    default:
	
      dc.lf2mf(frag, lfn, mfn);

      switch(Structure::top(frag)) {// top
      case Molecule::ASYM:

	for(int i = 0; i < 3; ++i)
	  imz += Structure::fragment(frag).imom(i) * mfn[i] * mfn[i];
	break;

      case Molecule::PROLATE:
	  
	imz += Structure::fragment(frag).imom(1) * (1. - mfn[0] * mfn[0]);
	k     -= _m_proj[frag] * mfn[0];
	break;
	
      case Molecule::OBLATE:
	
	imz += Structure::fragment(frag).imom(1) * (1. - mfn[2] * mfn[2]);
	k     -= _m_proj[frag] * mfn[2];
	break;
	
      }// top
    }// type
  }// frag cycle

    
  if(imz <= 0.) {
    std::cerr << funame << "WARNING: inertia moment projection on intermolecular frame is not positive:"
	      << imz << "\n";
    return -1.;
  }

  e -= k * k / 2. / imz;
      
  if(e <= 0.)
    return -1.;
  
  return std::pow(e, (double)_dim/ 2.) / std::sqrt(imz);
}

int LongRange::euler_angle_index (int frag, LongRange::ang_t angle) throw(Error::General)
{
  static const char funame [] = "LongRange::euler_angle_index";

#ifdef DEBUG

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

#endif

  IndexMapIterator it = index_map.find(std::pair<int, ang_t>(frag, angle));

  if(it == index_map.end())
    return -1;
  else
    return it->second;
}

// euler angles to dynamic (angular) coordinates converter
void LongRange::ang2dc (const std::vector<double>& angle, Dynamic::Coordinates& dc) throw(Error::General)
{
  static char funame [] = "LongRange::ang2dc: ";

#ifdef DEBUG

  if(!isinit()) {
    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  if(angle.size() != orientational_dimension()) {

    std::cerr << funame << "wrong angle dimension\n";
    throw Error::Logic();

  }

#endif

  int itemp;
  double polar[2]; // polar angles
  double theta, phi, psi;
  double q1 [4], q2 [4], q3 [4];

  for(int frag = 0; frag < 2; ++frag) {// fragment cycle
    switch(symmetry_type(frag)) {// symmetry type
    case SPHERICAL:

      break;

    case LINEAR:

      itemp = euler_angle_index(frag, THETA);

#ifdef DEBUG
      if(itemp < 0) {
	std::cerr << funame << "wrong angle index: " << itemp << "\n";
	throw Error::Logic();
      }
#endif

      polar[0] = angle[itemp]; // theta angle

      itemp = euler_angle_index(frag, PHI); // phi angle
      if(itemp < 0)
	polar[1] = 0.;
      else
	polar[1] = angle[itemp];
      
      polar2cart(polar, q1);
      dc.write_ang_pos(frag, q1);
      break;
		 
    default:// nonlinear

      itemp = euler_angle_index(frag, THETA);

#ifdef DEBUG
	
      if(itemp < 0) {
	std::cerr << funame << "wrong angle index: " << itemp << "\n";
	throw Error::Logic();
      }

#endif
      theta = angle[itemp];
      
      switch(Structure::top(frag)) {// top
      case Molecule::PROLATE:

	axis_quaternion(1, theta - M_PI / 2., q1); // rotation about Y axis on theta-pi/2 to bring X axis to Z position
	dc.write_ang_pos(frag, q1);

	break;

      case Molecule::OBLATE:

	axis_quaternion(1, theta, q1);
	dc.write_ang_pos(frag, q1);
	break;
	
      case Molecule::ASYM:

	itemp = euler_angle_index(frag, PSI);

#ifdef DEBUG
	
	if(itemp < 0) {
	  std::cerr << funame << "wrong angle index: " << itemp << "\n";
	  throw Error::Logic();
	}

#endif
	psi = angle[itemp];

	axis_quaternion(1, theta, q1); //rotation about Y axis on theta angle
	axis_quaternion(2, psi,   q2); //rotation about Z axis on psi angle
	Quaternion::qprod(q1, q2, q3);
	dc.write_ang_pos(frag, q3);
	
	break;

      default:
	
	std::cerr << funame << "wrong top: " << Structure::fragment(frag).top() << "\n";
	throw Error::Logic();
      }// top

      itemp = euler_angle_index(frag, PHI);
      if(itemp < 0)
	break;
      phi = angle[itemp];

      for(int i = 0; i < 4; ++i)
	q1[i] = dc.ang_pos(frag)[i];

      axis_quaternion(2, phi, q2); // rotation about Z axis on phi angle
      Quaternion::qprod(q2, q1, q3);

      dc.write_ang_pos(frag, q3);
    }//symmetry type
  }// fragment cycle
}

// number of states integrator
double LongRange::StatesNumberDensity::integral (double distance, double* mep) const throw(Error::General)
{
  static const char funame [] = "LongRange::StatesNumberDensity::integral: ";

  const double theta_step = M_PI / double(angular_grid_size() + 1);
  const double   phi_step = 2. * M_PI / double(angular_grid_size());

  int itemp;
  double dtemp;

  int   ang_index;
  ang_t ang_type;
  double weight;
  double rho_val;
  double poten;

  std::vector<double> euler_angle(orientational_dimension());

  Dynamic::Coordinates dc;
  for(int i = 0; i < 2; ++i)
    dc.orb_pos(i) = 0.;
  dc.orb_pos(2) = distance;

  double res = -1.;
  MultiIndexConvert multi_grid(std::vector<int>(orientational_dimension(), angular_grid_size()));
  for(int grid_index = 0; grid_index < multi_grid.size(); ++grid_index) {
 
    // set Euler angles and weight
    weight = 1.;
    std::vector<int> grid_point = multi_grid(grid_index);

    for(IndexMapIterator it = index_map.begin(); it != index_map.end(); ++it) {
      ang_index = it->second;
      ang_type  = it->first.second;
      switch(ang_type) {
      case THETA:
	
	itemp = grid_point[ang_index];
	dtemp = theta_step * double(itemp + 1); 

	euler_angle[ang_index] = dtemp;
	weight *= std::sin(dtemp);
	break;

      default: // PHI and PSI

	euler_angle[ang_index] = phi_step * double(grid_point[ang_index]) ; 
	break;

      }
    }

    // conversion of Euler angles to dynamic coordinates
    ang2dc(euler_angle, dc);

    if(mep) {
      rho_val = (*this)(dc, &poten);
      if(!grid_index || poten < *mep)
	*mep = poten;
    }
    else
      rho_val = (*this)(dc);

    // integration
    if(rho_val > 0.)
      increment(res, rho_val * weight);
  }

  // normalization
  if(res > 0.) {
    for(IndexMapIterator it = index_map.begin(); it != index_map.end(); ++it) {
      ang_type = it->first.second;
      switch(ang_type) {
      case THETA:
	
	res *= theta_step;
	break;

      default:

	res *= phi_step; 
	break;

      }
    }
    res *= this->norm_factor();
    if(dynamic_cast<const EDensity*>(this)) {
      res *= distance * distance;
    }
  }

  return res;
}

/*******************************************************************************************
 *************************************** Minimum Energy Search *****************************
 *******************************************************************************************/

double LongRange::minimum_energy (double distance)
{
  static const char funame [] = "LongRange::minimum_energy: ";

  static const int    theta_size = 3;
  static const double theta_step = M_PI / double(theta_size + 1);
  static const int      phi_size = 8;
  static const double   phi_step = 2. * M_PI / double(phi_size);

  int itemp;
  double dtemp;

  int   ang_index;
  ang_t ang_type;
  double ener, min_ener, max_ener;

  std::vector<double>     euler_angle(orientational_dimension());
  std::vector<double> min_euler_angle(orientational_dimension()); 
  std::vector<double> max_euler_angle(orientational_dimension());

  Dynamic::Coordinates dc;
  for(int i = 0; i < 2; ++i)
    dc.orb_pos(i) = 0.;
  dc.orb_pos(2) = distance;

  std::vector<int> grid_size(orientational_dimension());
  for(IndexMapIterator it = index_map.begin(); it != index_map.end(); ++it) {
    ang_index = it->second;
    ang_type  = it->first.second;
    switch(ang_type) {
    case THETA: // theta
	
      grid_size[ang_index] = theta_size;
      break;

    default:   // PHI and PSI

      grid_size[ang_index] = phi_size;
      break;

    }
  }
  MultiIndexConvert multi_grid(grid_size);
  for(int grid_index = 0; grid_index < multi_grid.size(); ++grid_index) {
 
    // set Euler angles
    std::vector<int> grid_point = multi_grid(grid_index);

    for(IndexMapIterator it = index_map.begin(); it != index_map.end(); ++it) {
      ang_index = it->second;
      ang_type  = it->first.second;
      switch(ang_type) {
      case THETA: // theta
	
	euler_angle[ang_index] = theta_step * double(grid_point[ang_index] + 1);
	break;

      default:    // phi and psi

	euler_angle[ang_index] = phi_step * double(grid_point[ang_index]) ; 
	break;

      }
    }

    // conversion of Euler angles to dynamic coordinates
    ang2dc(euler_angle, dc);
    // energy
    ener = pot(dc);

    // minimal energy
    if(!grid_index) {
      min_ener = ener;
      max_ener = ener;
      min_euler_angle = euler_angle;
      max_euler_angle = euler_angle;
    }
    else if(ener < min_ener) {
      min_ener   = ener;
      min_euler_angle = euler_angle;
    }
    else if(ener > max_ener) {
      max_ener   = ener;
      max_euler_angle = euler_angle;
    }
  }  


  // minimal energy search
  if(max_ener <= min_ener)
    std::cerr << funame << "WARNING: potential is flat, no minimum energy search is performed\n";
  else {
    ang2dc(min_euler_angle, dc);
    min_ener = gradient_search(dc, max_ener - min_ener);
  }

  return min_ener;
  
}

// position derivative vector
extern "C" void LongRange::pos2drv (const double& t, const double* pos0, double* drv0, void*, void* dcp)
{
  Dynamic::Coordinates& dc = *static_cast<Dynamic::Coordinates*>(dcp);
  
  const double* pos = pos0;
  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:

      break;

    default:

      dc.write_ang_pos(frag, pos);
      pos += Structure::pos_size(frag);      

    }

  D3::Vector force [3];
  pot(dc, force);

  D3::Vector* const torque = force + 1;


  pos = pos0;
  double* drv = drv0;
  double q [4];
  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:

      break;

    case LINEAR:

      D3::vprod(torque[frag], pos, drv);

      pos += Structure::pos_size(frag);      
      drv += Structure::pos_size(frag);      
      break;

    default: // non-linear, non-spherical fragment

      q[0] = 0.;
      for(int i = 0; i < 3; ++i)
	q[i + 1] = torque[frag][i] / 2.;
      Quaternion::qprod(pos, q, drv);

      pos += Structure::pos_size(frag);      
      drv += Structure::pos_size(frag);      
    }
}

double LongRange::gradient_search(Dynamic::Coordinates& dc, double evar)
{
  static const char funame [] = "LongRange::minimum_energy_search: ";

  static const double  rel_tol = 1.e-5; // relative integration tolerance
  static const double  abs_tol = 1.e-5; // absolute integration tolerance
  static const double torq_tol = 1.e-2; // torque tolerance

  int itemp;

  itemp = 0;
  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:

      break;

    default:

      itemp += Structure::pos_size(frag);      
    }

  Array<double> var(itemp);

  double* pos = var;
  for(int frag = 0; frag < 2; ++frag)
    switch(symmetry_type(frag)) {
    case SPHERICAL:

      break;

    default:

      for(int i = 0; i < Structure::pos_size(frag); ++i)
	pos[i] = dc.ang_pos(frag, i); 
      pos += Structure::pos_size(frag);      
    }

  Slatec::AdamSolver adams(var.size(), pos2drv, rel_tol, abs_tol);
  Slatec::AdamSolver::Mode mode = Slatec::AdamSolver::RESTART;
  double time = 0.;
  double timeout;
  D3::Vector force[3];
  D3::Vector* const torque = force + 1;
  double ener, torq_len;
  const double step = 3. / evar;
  
  //std::cout << "Minimum Energy Search:\n";
  do {
    timeout = time + step;
    adams.run(time, var, timeout, static_cast<void*>(&dc), mode);
    mode = Slatec::AdamSolver::CONTINUE;

    pos = var;
    for(int frag = 0; frag < 2; ++frag)
      switch(symmetry_type(frag)) {
      case SPHERICAL:

	break;

      default:

	dc.write_ang_pos(frag, pos);
	pos += Structure::pos_size(frag);      
      }
    ener = pot(dc, force);

    torq_len = 0.;
    for(int frag = 0; frag < 2; ++frag)
      switch(symmetry_type(frag)) {
      case SPHERICAL:

	break;

      default:

	torq_len += torque[frag].vlength();

      }

  } while(torq_len > evar * torq_tol);

  return ener;
}

/*******************************************************************************************
 *************************************** Number of States Classes **************************
 *******************************************************************************************/

// integral over total angular momentum
double LongRange::JIntegral::_integral() const
{
  double amom = step;
  double val, res = -1.;

  double old_dist = -1.;

  while(1) {
    _num->set_amom(amom);
    
    val = (*_num)();
    if(old_dist < 0.)
      old_dist = _num->distance();

#ifdef DEBUG
    std::cout << "J  = " << std::setw(7) << amom 
	      << "  J-Number = " << std::setw(13) << val
	      << std::endl;
    std::cout << "\n";
#endif

    if(val > 0.) {
      increment(res, val * amom);
      arg_var += val * amom * amom * amom * 2. * step;
    }
    else
      break;

    amom += step;    
  }

  if(res > 0.)
    res *= 2. * step;

  arg_max = amom > arg_max ? amom : arg_max;
  _num->set_dist(old_dist);
  return res;
}


// integral over m projections, if any
double LongRange::MIntegral::_integral (int index) const
{
  static const char funame [] = "LongRange::MIntegral::_integral: ";

  int itemp;
  double dtemp;

  if(index == symmetric_fragment_size()) {
    return (*_num)(); 
  }

  for(int frag = 0; frag < symmetric_fragment_size(); ++frag)
    if(index == frag) {

      int    proj_size;
      double proj, proj_max;
      double ener_max = ener() - _num->ener_min();

      double old_dist = -1.;

      switch(frag) {
      case 0:

	if(ener_max <= 0.) 
	  return -1.;
	proj_max = std::sqrt(2. * symmetry_axis_imom[0] * ener_max);
	break;

      case 1:
	
	ener_max -= _num->m_proj(0) * _num->m_proj(0) / 2. / symmetry_axis_imom[0];
	if(ener_max <= 0.) 
	  return -1.;
	proj_max = std::sqrt(2. * symmetry_axis_imom[1] * ener_max);
	break;
       
      } 
      
      itemp = (int)std::floor(proj_max / step[frag]);
      proj  = -(double)itemp * step[frag];
      proj_size = 2 * itemp + 1;

      double val, res = -1.;
      for(int i = 0; i < proj_size; ++i, proj += step[frag]) {
	_num->set_m_proj(frag, proj);

	val = _integral(frag + 1);

	if(old_dist < 0.)
	  old_dist = _num->distance();


#ifdef DEBUG
	std::cout << "M" << frag << " = " << std::setw(7) << proj 
		  << "  M-Number = " << std::setw(13) << val
		  << "\n";
#endif

	if(val > 0.) {
	  increment(res, val);
	  arg_var[frag] += val * proj * proj * step[frag];
	}
      }
      if(res > 0.)
	res *= step[frag];
      
      arg_max[frag] = proj_max > arg_max[frag] ? proj_max : arg_max[frag];
      _num->set_dist(old_dist);
      return res;
    }
}

// integral over K projection
double LongRange::KIntegral::_integral() const
{
  double m_proj_max = 0.;
  for(int i = 0; i < symmetric_fragment_size(); ++i) 
    m_proj_max += m_proj(i) > 0. ? m_proj(i) : - m_proj(i);

  double val, res = -1.;
  double weight = 1.;
  bool first = true;
  double proj = 0.;

  double old_dist = -1.;

  while(1) {
    _num->set_k_proj(proj);

    val = (*_num)();
    if(old_dist < 0.)
      old_dist = _num->distance();

#ifdef DEBUG
    std::cout << "K  = "          << std::setw(7) << proj 
	      << "  K-Number = " << std::setw(13) << val
	      << "\n";
#endif

    if(val > 0.) {
      increment(res, weight * val);
      arg_var += weight * val * proj * proj * step;
    }
    else if(_num->k_proj() >= m_proj_max)
	break;
 
    proj += step;
    if(first) {
      first = false;
      weight = 2.;
    }
  }

  if(res > 0.)
    res *= step;

  if(!first)
    arg_max = proj > arg_max ? proj : arg_max;
  _num->set_dist(old_dist);
  return res;
}

double LongRange::Opt::_fun (double x) throw (Math::MinimumSearch::Exception)
{
  _opt_step_count++;

  double dist = std::exp(x);
  set_dist(dist);
  double res = integral();

  if(res <= 0.)
    throw Math::MinimumSearch::Exception(dist);

  return std::log(res);
}


double LongRange::Opt::opt_integral () throw(Error::General)
{
  static const char funame [] = "LongRange::Opt::opt_integral: ";

#ifdef DEBUG
  std::cout << funame << "starting guess: distance = " << distance() << "\n";
#endif

  _opt_step_count = 0;

  double x, y;

  try {
    
    if(distance() > 0.) {
      y = find(std::log(distance()), step, &x);
    }
    else {
      std::cerr << funame << "no initial guess\n";
      throw Error::Logic();
    }

    // finding minimum using initial distance as a guess

#ifdef DEBUG
    std::cout << funame << "# of optimization steps = " << _opt_step_count << "\n";
#endif

    set_dist(std::exp(x));
    return std::exp(y);

  }
  catch(Math::MinimumSearch::Exception) {

#ifdef DEBUG
    std::cout << funame << "# of optimization steps = " << _opt_step_count << "\n";
#endif

    return -1.;
  }  
}

