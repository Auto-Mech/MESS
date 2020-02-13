#include "opt.hh"
#include "io.hh"
#include "key.hh"
#include "random.hh"

#include <boost/numeric/odeint.hpp>

namespace boost {
  //
  namespace numeric {
    //
    namespace odeint {
      //
      template<>
      struct is_resizeable<Coord::Cartesian> {
	//
	typedef boost::true_type type;
	
	static const bool value = type::value;
      };
    }
  }
}

/**********************************************************************************************************
 *************************************** HARDING POTENTIAL WRAPPER ****************************************
 **********************************************************************************************************/

double Opt::HardPot::evaluate (const Cartesian& x) const
{
  const char funame [] = "Opt::HardPot::evaluate: ";

  if(atom_size() != x.atom_size()) {
    //
    std::cerr << funame << "dimensions mismatch: " << atom_size() << " vs. " << x.atom_size() << "\n";

    throw Error::Logic();
  }
  
  return potential(x);
}

/**********************************************************************************************************
 ***************************** CONSTRAINED OPTIMIZATION IN CARTESIAN SPACE ********************************
 **********************************************************************************************************/

int Opt::CcOpt::_grad_count;

void Opt::CcOpt::execute (Coord::Cartesian& x) const
{
  const char funame [] = "Opt::CcOpt::execute: ";

  namespace odeint = boost::numeric::odeint;

  double dtemp;
  int    itemp;
  
  IO::Marker fumarker(funame);

  _grad_count = 0;
  
  Cartesian grad(x.size());
  
  (*this)(x, grad);
 
  dtemp = vlength(grad);

  IO::log << IO::log_offset << "initial gradient projection length = " << dtemp << "\n";
  
  if(dtemp <= _grad_tol) {
    //
    IO::log << IO::log_offset << "gradient projection is small enough: no optimization necessary\n";

    return;
  }

  if(_constrain.size() && use_constrain) {
    //
    IO::log << IO::log_offset << "initial constrain values:";

    for(_con_t::const_iterator cit = _constrain.begin(); cit != _constrain.end(); ++cit) {
      //
      dtemp = (*cit)->evaluate(x);

      switch((*cit)->type()) {
	//
      case Coord::DISTANCE:
	//
	dtemp /= Phys_const::angstrom;

	break;
      
      case Coord::ANGLE:
	//
	dtemp *= 180. / M_PI;

	break;

      case Coord::DIHEDRAL:
	//
	dtemp *= 180. / M_PI;

	break;
      }
    
      IO::log << "   " << dtemp;
    }
  
    IO::log << std::endl;
  }
  
  try {
    //
    double time = 0.;
    
    while(1) {
      //
      IO::log << IO::log_offset
	      << std::setw(13) << "time"
	      << std::setw(13) << "gradient"
	      << std::setw(13) << "target"
	      << "\n";

      double tmax = time + _time_length;
      
      itemp = odeint::integrate(*this, x, time, tmax, _time_step);

      time = tmax;
    
      IO::log << IO::log_offset << "steps # = " << itemp << std::endl;
    }
  } catch(_Fin) {
    //
    (*this)(x, grad);

    dtemp = vlength(grad);
    
    IO::log << IO::log_offset << "final gradient projection length  = " << std::setw(13) << dtemp  << "\n";
    
    IO::log << IO::log_offset << "gradient calls #       = " << std::setw(13) << _grad_count << "\n";

    if(_constrain.size() && use_constrain) {
      //
      IO::log << IO::log_offset << "final constrain values:";

      for(_con_t::const_iterator cit = _constrain.begin(); cit != _constrain.end(); ++cit) {
	//
	dtemp = (*cit)->evaluate(x);

	switch((*cit)->type()) {
	  //
	case Coord::DISTANCE:
	  //
	  dtemp /= Phys_const::angstrom;

	  break;
      
	case Coord::ANGLE:
	  //
	  dtemp *= 180. / M_PI;

	  break;

	case Coord::DIHEDRAL:
	  //
	  dtemp *= 180. / M_PI;

	  break;
	}
    
	IO::log << "   " << dtemp;
      }
  
      IO::log << std::endl;
    }
  }
}

void Opt::CcOpt::operator() (const Cartesian& x, Cartesian& dx, double time) const
{
  const char funame [] = "Opt::CcOpt::operator(): ";

  int    itemp;
  double dtemp;

  ++_grad_count;
  
  if(!x.size()) {
    //
    std::cerr << funame << "zero configuration space vector dimension\n";

    throw Error::Init();
  }

  dx.resize(x.size());
  
  for(int i = 0; i < dx.size(); ++i)
    //
    dx[i] = -pot()->grad(x, i);

  if(use_constrain && _constrain.size()) {
    //
    // constrain subspace
    //
    std::vector<Cartesian> cspace(_constrain.size(), x);

    int cc = 0;

    for(_con_t::const_iterator cit = _constrain.begin(); cit != _constrain.end(); ++cit, ++cc) {
      //
      for(int i = 0; i < x.size(); ++i)
	//
	cspace[cc][i] = (*cit)->grad(x, i);

      dtemp = normalize(cspace[cc]);

      if(dtemp == 0.) {
	//
	std::cerr << funame << "zero constrain vector\n";

	throw Error::Range();
      }
    
      for(int d = 0; d < cc; ++d)
	//
	orthogonalize(cspace[cc], cspace[d]);

      if(cc) {
	//
	dtemp = normalize(cspace[cc]);
    
	if(dtemp < _ci_tol) {
	  //
	  std::cerr << funame << "constrains are linearly dependent: " << dtemp << " vs " << _ci_tol << "\n";

	  throw Error::Run();
	}
      }
    }

    // potential gradient component, orthogonal to constrain subspace
    //
    for(cc = 0; cc < cspace.size(); ++cc)
      //
      orthogonalize(dx, cspace[cc]);
  }
  
  if(time != 0.) {
    //
    dtemp = vlength(dx);

    IO::log << IO::log_offset
	    << std::setw(13) << time
	    << std::setw(13) << dtemp
	    << std::setw(13) << _grad_tol
	    << std::endl;
  
    if(dtemp < _grad_tol)
      //
      throw _Fin();
  }
}

void Opt::CcOpt::set (std::istream& from)
{
  const char funame [] = "Opt::CcOpt::set:";
  
    int    itemp;
    double dtemp;

  IO::Marker funame_marker(funame);

  KeyGroup CcOptGroup;

  Key  tol_key("GradientTolerance[a.u.]"       );
  Key  cit_key("ConstrainIndependenceTolerance");
  Key step_key("IntegrationStep"               );
  Key time_key("IntegrationTime"               );
    
  std::string token, comment, stemp;

  while(from >> token) {
    //
    // end of input
    //
    if(token == IO::end_key()) {
      //
      std::getline(from, comment);
      //
      break;
    }
    // gradient tolerance
    //
    else if(token == tol_key) {
      //
      if(!(from >> _grad_tol)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(_grad_tol <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << _grad_tol << "\n";

	throw Error::Range();
      }
    }
    // constrain (local) independence tolerance
    //
    else if(token == cit_key) {
      //
      if(!(from >> _ci_tol)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(_ci_tol <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << _ci_tol << "\n";

	throw Error::Range();
      }
    }
    // initial time step
    //
    else if(token == step_key) {
      //
      if(!(from >> _time_step)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(_time_step <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << _time_step << "\n";

	throw Error::Range();
      }
    }
    // integration time
    //
    else if(token == time_key) {
      //
      if(!(from >> _time_length)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }

      if(_time_length <= 0.) {
	//
	std::cerr << funame << token << ": out of range: " << _time_length << "\n";

	throw Error::Range();
      }
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      std::cerr << funame << "unknown keyword " << token << "\n";
      
      Key::show_all(std::cerr);
      
      std::cerr << "\n";
      
      throw Error::Init();
    }
  }

  if(!from) {
    //
    std::cerr << funame << "corrupted\n";

    throw Error::Input();
  }
}

/********************************************************************************************************
 *********************** OPTIMIZATION AND IMPORTANCE SAMPLING IN Z-MATRIX COORDINATES********************
 ********************************************************************************************************/

// gradient tolerance
//
double Opt::ZOpt::grad_tol = 1.e-5;

// maximal optimization step
//
double Opt::ZOpt::max_opt_step = 0.1;

// maximal number of optimization iterations
//
int Opt::ZOpt::max_opt_count = 100;

// differentiation step
//
double Opt::ZOpt::diff_step = .001;

// potential
//
ConstSharedPointer<Coord::CartFun> Opt::ZOpt::pot;

double Opt::ZOpt::_con_pot (const std::vector<double>& cpos) const
{
  const char funame [] = "Opt::ZOpt::_con_pot: ";

  int    itemp;
  double dtemp;
  
  if(!pot) {
    //
    std::cerr << funame << "potential not initialized\n";

    throw Error::Init();
  }

  if(cpos.size() != _con_modes.size()) {
    //
    std::cerr << funame << "number of conserved modes mismatch: " << cpos.size() << " vs " << _con_modes.size() << "\n";

    throw Error::Logic();
  }

  ZData zpos = _zmin;

  itemp = 0;
  
  for(mode_t::const_iterator cit = _con_modes.begin(); cit != _con_modes.end(); ++cit, ++itemp)
    //
    zpos[cit->first] = cpos[itemp];

  return pot->evaluate((Cartesian)zpos);
}

Lapack::Vector Opt::ZOpt::_con_grad () const
{
  int    itemp;
  double dtemp;
  
  std::vector<double> cpos(_con_modes.size());

  itemp = 0;
  
  for(mode_t::const_iterator cit = _con_modes.begin(); cit != _con_modes.end(); ++cit, ++itemp)
    //
    cpos[itemp] = _zmin[cit->first];

  Lapack::Vector res(_con_modes.size());

  for(int i  = 0; i < _con_modes.size(); ++i) {
    //
    cpos[i] += diff_step;

    res[i] = _con_pot(cpos);

    cpos[i] -= 2. * diff_step;

    res[i] -= _con_pot(cpos);

    cpos[i] += diff_step;
    
    res[i] /= 2. * diff_step;
  }

  return res;
}

Lapack::SymmetricMatrix Opt::ZOpt::_con_hess () const
{
  int    itemp;
  double dtemp;
  
  std::vector<double> cpos(_con_modes.size());

  itemp = 0;
  
  for(mode_t::const_iterator cit = _con_modes.begin(); cit != _con_modes.end(); ++cit, ++itemp)
    //
    cpos[itemp] = _zmin[cit->first];

  const double e0 = _con_pot(cpos);

  Lapack::SymmetricMatrix res(_con_modes.size());

  for(int i  = 0; i < _con_modes.size(); ++i)
    //
    for(int j  = i; j < _con_modes.size(); ++j)
      //
      if(i != j) {
	//
	cpos[i] += diff_step;

	cpos[j] += diff_step;
	
	res(i, j) = _con_pot(cpos);

	cpos[i] -= 2. * diff_step;

	res(i, j) -= _con_pot(cpos);

	cpos[j] -= 2. * diff_step;

	res(i, j) += _con_pot(cpos);

	cpos[i] += 2. * diff_step;

	res(i, j) -= _con_pot(cpos);

	cpos[i] -= diff_step;
	
	cpos[j] += diff_step;
	
	res(i, j) /= 4. * diff_step * diff_step;
      }
      else {
	//
	cpos[i] += diff_step;
	
	res(i, i) = _con_pot(cpos) - 2. * e0;

	cpos[i] -= 2. * diff_step;

	res(i, i) += _con_pot(cpos);

	cpos[i] += diff_step;
	
	res(i, i) /= diff_step * diff_step;
      }

  return res;
}

Opt::ZOpt::ZOpt (const ZData& zinit, const mode_t& cm) : _zmin(zinit), _con_modes(cm)
{
  const char funame [] = "Opt::ZOpt::ZOpt: ";

  IO::Marker fumark(funame);
  
  double dtemp;
  int    itemp;

  int count = 0;

  IO::log << IO::log_offset << std::setw(5) << "#" << std::setw(13) << "Grad[au]" << std::setw(13) << "Target" << std::endl;
  
  while(count++ < max_opt_count) {
    //
    _fc_eval_sqrt = _con_hess().eigenvalues(&_fc_evec);

    Lapack::Vector dx = _con_grad();
    
    dtemp = vlength(dx);

    IO::log << IO::log_offset << std::setw(5) << count << std::setw(13) << dtemp << std::setw(13) << grad_tol << std::endl;
    
    if(dtemp < grad_tol)
      //
      break;

    dx = dx * _fc_evec;

    for(int i = 0; i < dx.size(); ++i) {
      //
      if(_fc_eval_sqrt[i] <= 0.) {
	//
	IO::log << IO::log_offset << i + 1 << "-th force constant eigenvalue is negative" << std::endl;
	
	if(dx[i] < 0.) {
	  //
	  dx[i] = -max_opt_step;
	}
	else if(dx[i] > 0.) {
	  //
	  dx[i] = max_opt_step;
	}
      }
      else {
	//
	dtemp = dx[i] / _fc_eval_sqrt[i];

	if(dtemp < -max_opt_step) {
	  //
	  dx[i] = -max_opt_step;
	}
	else if(dtemp > max_opt_step) {
	  //
	  dx[i] = max_opt_step;
	}
	else
	  //
	  dx[i] = dtemp;
      }
    }
    
    dx = _fc_evec * dx;

    itemp = 0;

    for(mode_t::const_iterator cit = _con_modes.begin(); cit != _con_modes.end(); ++cit, ++itemp) {
      //
      dtemp = _zmin[cit->first] - dx[itemp];

      if(dtemp < cit->second.first) {
	//
	IO::log << IO::log_offset << cit->first << " proposed value is below lower limit, " << dtemp << ": clip it" << std::endl;
	
	_zmin[cit->first] = cit->second.first;
      }
      else if(dtemp > cit->second.second) {
	//
	IO::log << IO::log_offset << cit->first << " proposed value is above upper limit, " << dtemp << ": clip it" << std::endl;
	
	_zmin[cit->first] = cit->second.second;
      }
      else
	//
	_zmin[cit->first] = dtemp;
    }
  }

  if(count > max_opt_count) {
    //
    IO::log << IO::log_offset << "maximal number of iterations reached: no convergence" << std::endl;

    std::cerr << funame << "maximal number of iterations reached: no convergence\n";

    throw NoConv();
  }

  for(int i = 0; i < _con_modes.size(); ++i) {
    //
    dtemp = _fc_eval_sqrt[i];

    if(dtemp <= 0.) {
      //
      std::cerr << funame << "force constant matrix at minimum is not positively defined\n";

      throw Error::Run();
    }

    _fc_eval_sqrt[i] = std::sqrt(dtemp);
  }

  _mass_factor = std::sqrt(product(_zmin.inertia_moments())) / Lapack::Cholesky(_zmin.mobility_matrix()).det_sqrt();

  _ener_min = pot->evaluate((Cartesian)_zmin);
}

double Opt::ZOpt::anharmonic_correction (double temperature, int count_max, double* rerr, int* fail) const
{
  const char funame [] = "Opt::ZOpt::anharmonic_correction: ";

  static const double exp_pow_max = 100.;

  double dtemp;
  int    itemp;
  bool   btemp;
  
  const double tsqrt = std::sqrt(temperature);

  Lapack::Vector vtemp(_con_modes.size());

  double res = 0., var =  0.;

  int skip = 0;

  // sampling cycle
  //
  for(int count = 0; count < count_max; ++count) {
    //
    double eref = 0.;

    for(int i = 0; i < _con_modes.size(); ++i) {
      //
      dtemp = Random::norm();

      eref += dtemp * dtemp;
    
      vtemp[i] = dtemp / _fc_eval_sqrt[i] * tsqrt;
    }

    eref *= temperature / 2.;

    vtemp = _fc_evec * vtemp;

    ZData ztemp = _zmin;

    itemp = 0;

    btemp = false;
  
    for(mode_t::const_iterator cit = _con_modes.begin(); cit != _con_modes.end(); ++cit, ++itemp) {
      //
      dtemp = ztemp[cit->first] + vtemp[itemp];

      // check that the non-fluxional coordinate is in the allowed window
      //
      if(dtemp < cit->second.first || dtemp > cit->second.second) {
	//
	switch(ztemp.type(cit->first)) {
	  //
	case DISTANCE:
	  //
	  dtemp /= Phys_const::angstrom;

	  break;
	  
	default:
	  //
	  dtemp *= 180. / M_PI;
	}

	IO::log << funame << "WARNING: T = " << temperature / Phys_const::kelv
		<< "K, z-matrix variable " << cit->first << " out of limits: " << dtemp << "\n";

	++skip;
	
	btemp = true;

	break;
      }

      ztemp[cit->first] = dtemp;
    }

    if(btemp)
      //
      continue;

    dtemp = (pot->evaluate((Cartesian)ztemp) - eref - _ener_min) / temperature;

    if(dtemp > exp_pow_max)
      //
      continue;

    if(dtemp < -exp_pow_max) {
      //
      IO::log << funame << "WARNING: anharmonic correction too negative at T = " << temperature / Phys_const::kelv << "K: "
	      << dtemp * temperature / Phys_const::kcal << " kcal/mol: skipping the point" << std::endl;

      ++skip;
      
      continue;
    }

    dtemp = std::exp(-dtemp) / Lapack::Cholesky(ztemp.mobility_matrix()).det_sqrt()
      * std::sqrt(product(ztemp.inertia_moments())) / mass_factor();

    res += dtemp;

    var += dtemp * dtemp;
    //
  }// sampling cycle
  //

  res /= (double)count_max;

  var /= (double)count_max;

  // relative error
  //
  if(rerr)
    //
    *rerr = std::sqrt(var - res * res) / res / std::sqrt((double)count_max);

  // failed points #
  //
  if(fail)
    //
    *fail = skip;
  
  return res;
}
