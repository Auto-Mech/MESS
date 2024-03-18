#include "model_2d.hh"
#include "key.hh"

/*********************************** 2D RIGID ROTOR **************************************/

double Model::RigidRotor2::eps = 1.e-3;

Model::RigidRotor2::RigidRotor2 (IO::KeyBufferStream& from, const std::vector<Atom>& atom, int m)
  : Core2(m), _symmetry_factor(1.)
{
  const char funame [] = "Model::RigidRotor2::RigidRotor2: ";

  IO::Marker funame_marker(funame);

  int         itemp;
  double      dtemp;
  std::string stemp;
  
  KeyGroup RigidRotor2Model;

  Key      symm_key("SymmetryFactor");

  std::string token, comment;
  
  while(from >> token) {
    //
    // input end
    //
    if(IO::end_key() == token) {
      //
      std::getline(from, comment);
      
      break;
    }
    // symmetry factor (number of symmetry operations)
    //
    else if(symm_key == token) {
      //
      if(!(from >> _symmetry_factor)) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": corrupted";
      }
      
      std::getline(from, comment);

      if(_symmetry_factor <= 0.) {
	//
	ErrOut err_out;
	
	err_out << funame << token << ": should be positive";
      }
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      ErrOut err_out;
      
      err_out << funame << "unknown keyword " << token << "\n" << Key::show_all() << "\n";
    }
  }
  
  if(!from) {
    //
    ErrOut err_out;
    
    err_out << funame << "input corrupted";
  }

  if(!atom.size()) {
    //
    std::cerr << funame << "geometry is not available\n";
      
    throw Error::Init();
  }
  
  Lapack::Vector imom = inertia_moment_matrix(atom).eigenvalues();

  // molecule type:

  // linear rotor
  //
  if(imom[0] / imom[1] < eps) {
    //
    _type = LINEAR_ROTOR;
    
    _rc.push_back(0.5 / imom[2]);
  }
  // spherical top
  //
  else if(imom[0] / imom[2] > 1 - eps) {
    //
    _type = SPHERICAL_TOP;
    
    _rc.push_back(0.5 / imom[2]);
  }
  // oblate symmetric top, B_z < B_x,y
  //
  else if(imom[1] < std::sqrt(imom[0] * imom[2])) {
    //
    _type = OBLATE_TOP; 

    dtemp = std::sqrt(imom[0] * imom[1]);

    _rc.push_back(0.5 / imom[2]); // B_z

    _rc.push_back(0.5 / dtemp);  // B_x,y
  }
  // prolate symmetric top, B_z > B_x,y
  //
  else {
    //
    _type = PROLATE_TOP;

    dtemp = std::sqrt(imom[1] * imom[2]);

    _rc.push_back(0.5 / dtemp);   // B_x,y

    _rc.push_back(0.5 / imom[0]); // B_z
  }
}

double Model::RigidRotor2::ground (double j) const
{
  const char funame [] = " Model::RigidRotor2::ground: ";

  return _rc[0] * j * j;
}

void Model::RigidRotor2::states (Array<double>& stat, double de, double j) const
{
  const char funame [] = " Model::RigidRotor2::states: ";

  double dtemp;
  int    itemp;

  if(mode() != DENSITY && mode() != NUMBER) {
    //
    std::cerr << funame << "wrong mode\n";

    throw Error::Logic();
  }
  
  if(!stat.size()) {
    //
    std::cerr << funame << "states array not initialized\n";

    throw Error::Init();
  }
  
  stat = 0.;

  int n;
  
  switch(_type) {
    //
  case LINEAR_ROTOR:
    //
    switch(mode()) {
      //
    case DENSITY:
      //
      stat[0] = 2.* j / _symmetry_factor / de;

      return;

    case NUMBER:
      //
      stat = 2. * j / _symmetry_factor;

      return;
    }

  case SPHERICAL_TOP:
    //
    switch(mode()) {
      //
    case DENSITY:
      //
      stat[0] = 4.* j * j / _symmetry_factor / de;

      return;

    case NUMBER:
      //
      stat = 4.* j * j / _symmetry_factor;
      
      return;
    }

  case PROLATE_TOP:
    //
    n = std::ceil((_rc[1] - _rc[0]) * j * j / de);

    switch(mode()) {
      //
    case DENSITY:
      //
      dtemp = 4. * j * j / _symmetry_factor / std::sqrt((double)n) / de;
    
      for(int i = 0; i < n; ++i) {
	//
	if(i == stat.size())
	  //
	  break;

	stat[i] = (std::sqrt((double)(i + 1)) - std::sqrt((double)i)) * dtemp;
      }

      return;

    case NUMBER:
      //
      dtemp = 4. * j * j / _symmetry_factor / std::sqrt((double)n);
    
      for(int i = 0; i < n; ++i) {
	//
	if(i == stat.size())
	  //
	  break;

	stat[i] = std::sqrt((double)(i + 1)) * dtemp;
      }

      dtemp = 4. * j * j / _symmetry_factor;

      for(int i = n; i < stat.size(); ++i)
	//
	stat[i] = dtemp;

      return;
    }
    
  case OBLATE_TOP:
    //
    n = std::ceil((_rc[1] - _rc[0]) * j * j / de);

    switch(mode()) {
      //
    case DENSITY:
      //
      dtemp = 4. * j * j / _symmetry_factor / std::sqrt((double)n) / de;
    
      for(int i = 0; i < n; ++i) {
	//
	if(i == stat.size())
	  //
	  break;

	stat[i] = (std::sqrt((double)(n - i)) - std::sqrt((double)(n - i - 1))) * dtemp;
      }

      return;

    case NUMBER:
      //
      dtemp = 4. * j * j / _symmetry_factor / std::sqrt((double)n);
    
      for(int i = 0; i < n; ++i) {
	//
	if(i == stat.size())
	  //
	  break;

	stat[i] = (std::sqrt((double)n) - std::sqrt((double)(n - i - 1))) * dtemp;
      }

      dtemp = 4. * j * j / _symmetry_factor;

      for(int i = n; i < stat.size(); ++i)
	//
	stat[i] = dtemp;

      return;
    }

  default:
    //
    std::cerr << funame << "wrong type\n";

    throw Error::Logic();
  }
}

void Model::RigidRotor2::convolute (Array<double>& stat, double de, double j) const
{
  const char funame [] = " Model::RigidRotor2::states: ";

  double dtemp;
  int    itemp;

  if(mode() != DENSITY && mode() != NUMBER) {
    //
    std::cerr << funame << "wrong mode\n";

    throw Error::Logic();
  }
  
  if(!stat.size()) {
    //
    std::cerr << funame << "states array not initialized\n";

    throw Error::Init();
  }
  
  int n;

  // rotational density of states
  //
  Array<double> rdos;
  
  switch(_type) {
    //
  case LINEAR_ROTOR:
    //
    stat *= 2. * j / _symmetry_factor;

    return;
 
  case SPHERICAL_TOP:
    //
    stat *= 4.* j * j / _symmetry_factor;
      
    return;

  case PROLATE_TOP:
    //
    n = std::ceil((_rc[1] - _rc[0]) * j * j / de);

    rdos.resize(n);

    dtemp = 4. * j * j / _symmetry_factor / std::sqrt((double)n);

#pragma omp parallel for default(shared) schedule(static)
    //
    for(int i = 0; i < n; ++i)
	//
	rdos[i] = (std::sqrt((double)(i + 1)) - std::sqrt((double)i)) * dtemp;

    break;

  case OBLATE_TOP:
    //
    n = std::ceil((_rc[1] - _rc[0]) * j * j / de);

    rdos.resize(n);
    
    dtemp = 4. * j * j / _symmetry_factor / std::sqrt((double)n);

#pragma omp parallel for default(shared) schedule(static)
    //
    for(int i = 0; i < n; ++i)
      //
      rdos[i] = (std::sqrt((double)(n - i)) - std::sqrt((double)(n - i - 1))) * dtemp;

    break;
    
  default:
    //
    std::cerr << funame << "wrong type\n";

    throw Error::Logic();
  }

  Array<double> old_stat(stat);

  stat = 0.;

#pragma omp parallel for default(shared) schedule(static)
  //
  for(int e = 0; e < stat.size(); ++e)
    //
    for(int i = 0; i < rdos.size(); ++i) {
      //
      if(i > e)
	//
	break;
      
      stat[e] += old_stat[e - i] * rdos[i];
    }
}

double Model::RigidRotor2::weight (double t) const
{
  const char funame [] = "Model::RigidRotor2::weight: ";

  switch(_type) {
    //
  case LINEAR_ROTOR:
    //
    return t / _rc[0] / _symmetry_factor;

  case SPHERICAL_TOP:
    //
    return std::sqrt(M_PI * t / _rc[0]) * t / _rc[0] / _symmetry_factor;

  case OBLATE_TOP:
    //
    return std::sqrt(M_PI * t / _rc[0]) * t / _rc[1] / _symmetry_factor;

  case PROLATE_TOP:
    //
    return std::sqrt(M_PI * t / _rc[1]) * t / _rc[0] / _symmetry_factor;

  default:
    //
    std::cerr << funame << "wrong case\n";

    throw Error::Logic();
  }
}
