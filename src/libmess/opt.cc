#include "opt.hh"
#include "io.hh"
#include "key.hh"
#include "linpack.hh"

#include <boost/numeric/odeint.hpp>

ConstSharedPointer<Coord::CartFun> Opt::new_cartesian_function (std::istream& from)
{
  const char funame [] = "Opt::new_potential: ";

  std::string comment, token;

  if(!(from >> token)) {
    //
    std::cerr << funame << "corrupted\n";

    throw Error::Input();
  }

  // skip the rest of the line
  //
  std::getline(from, comment);

  // harding potential
  //
  if(token == "Harding")
    //
    return ConstSharedPointer<CartFun>(new HardPot(from));

  if(token == "Internal")
    //
    return ConstSharedPointer<CartFun>(new Internal(from));

  // unknown potential type
  //
  std::cerr << funame << "unknown keyword (potential type): " << token << "\n";

  Key::show_all(std::cerr);

  std::cerr << "\n";
  
  throw Error::Init();
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
 ************************* INTERNAL COORDINATE (DISTANCE, ANGLE, OR DIHEDRAL) *****************************
 **********************************************************************************************************/

Opt::Internal::Internal (std::istream& from)
{
  const char funame [] = "Opt::Internal::Internal: ";

  double dtemp;
  int    itemp;
  
  IO::LineInput lin(from);

  while(lin >> itemp) {
    //
    if(itemp <= 0) {
      //
      std::cerr << funame << size() + 1 << "-th atomic index out of range: " << itemp << "\n";

      throw Error::Range();
    }

    // Fortran style index to C one
    //
    --itemp;
    
    if(!add_dep(itemp)) {
      //
      std::cerr << funame << "duplicated atomic index: " << itemp + 1 << "\n";

      throw Error::Input();
    }

    push_back(itemp);
  }

  if(size() < 2 || size() > 4) {
    //
    std::cerr << funame << "number of atoms out of range: " << size() << "\n";

    throw Error::Input();
  }
}

double Opt::Internal::evaluate (const Cartesian& x) const
{
  const char funame [] = "Opt::Internal::evaluate: ";

  double dtemp;
  int    itemp;

  std::vector<D3::Vector> atom_pos(size());

  for(int i = 0; i < size(); ++i) {
    //
    itemp = (*this)[i];

    if(itemp >= x.size()) {
      //
      std::cerr << funame << i + 1 << "-th atomic index out of range: " << itemp << " vs. " << x.size() << "\n";

      throw Error::Range();
    }

    atom_pos[i] = x.atom_pos(itemp);
  }

  return internal(atom_pos);
}

/**********************************************************************************************************
 ***************************** CONSTRAINED OPTIMIZATION IN CARTESIAN SPACE ********************************
 **********************************************************************************************************/

void Opt::CcOpt::execute (Coord::Cartesian& x) const
{
  const char funame [] = "Opt::CcOpt::execute: ";

  namespace odeint = boost::numeric::odeint;

  double dtemp;
  int    itemp;
  
  IO::Marker funame_marker(funame);

  Cartesian grad(x.size());
  
  (*this)(x, grad, 0.);
 
  double proj = vlength(grad, grad.size());
    
  IO::log << IO::log_offset
	  << "initial gradient projection length = " << proj << "\n\n"
	  << IO::log_offset
	  << std::setw(10) << "steps #"
	  << std::setw(13) << "grad. proj."
	    << "\n";
  do {
    //
    int steps = odeint::integrate(*this, x, 0., _time_length, _time_step);
    
    (*this)(x, grad, 0.);
    
    proj = vlength(grad, grad.size());
    
    IO::log << IO::log_offset
	    << std::setw(10) << steps
	    << std::setw(13) << proj
	    << "\n";
    
  } while(proj > _grad_tol);
}

void Opt::CcOpt::operator() (const Cartesian& x, Cartesian& dx, double) const
{
  const char funame [] = "Opt::CcOpt::operator(): ";

  // constrain subspace
  //
  std::vector<Cartesian> cspace(_constrain.size(), x);
  
  for(int c = 0; c < _constrain.size(); ++c) {
    //
    for(int i = 0; i < dx.size(); ++i)
      //
      cspace[c][i] = _constrain[c]->grad(x, i);

    for(int d = 0; d < c; ++d)
      //
      ::orthogonalize(cspace[c], cspace[d]);
  }

  // potential gradient component, orthogonal to constrain subspace
  //
  for(int i = 0; i < dx.size(); ++i)
    //
    dx[i] = -_pot->grad(x, i);

  for(int c = 0; c < cspace.size(); ++c)
    //
    ::orthogonalize(dx, cspace[c]);
}

Opt::CcOpt::CcOpt (std::istream& from) : _grad_tol(1.e-5), _time_step(0.1), _time_length(10.)
{
  const char funame [] = "Opt::CcOpt::CcOpt:";
  
    int    itemp;
    double dtemp;

  IO::Marker funame_marker(funame);

  KeyGroup CcOptGroup;

  Key  pot_key("Potential"              );
  Key  con_key("Constrain"              );
  Key  tol_key("GradientTolerance[a.u.]");
  Key step_key("IntegrationStep"        );
  Key time_key("IntegrationTime"        );
    
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
    // potential energy surface
    //
    else if(token == pot_key) {
      //
      if(_pot) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      _pot = new_cartesian_function(from);
    }
    // constrain
    //
    else if(token == con_key) {
      //
      _constrain.push_back(new_cartesian_function(from));
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
  
  if(!_pot) {
    //
    std::cerr << funame << "potential not initialized\n";

    throw Error::Init();
  }
}
