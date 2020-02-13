#include "coord.hh"
#include "io.hh"

/******************************************************************************************
 ********** INTERNAL COORDINATE: INTERATOMIC DISTANCE, PLANE, OR DIHEDRAL ANGLES **********
 ******************************************************************************************/

double Coord::internal (const std::vector<D3::Vector>& atom_pos)
{
  const char funame [] = "Coord::internal: ";

  double dtemp;
  int    itemp;

  D3::Vector n, v1, v2;
  
  switch(atom_pos.size()) {
    //
    // interatomic distance
    //
  case 2:
    //
    return (atom_pos[1] - atom_pos[0]).vlength();

    // plane angle
    //
  case 3:
    //
    v1 = atom_pos[0] - atom_pos[1];
    
    v2 = atom_pos[2] - atom_pos[1];
    
    dtemp = vdot(v1, v2) / v1.vlength() / v2.vlength();
    
    if(dtemp < -1.)
      //
      dtemp = -1.;
    
    if(dtemp > 1.)
      //
      dtemp = 1.;

    return std::acos(dtemp);

    // dihedral angle
    //
  case 4:
    //
    n  = atom_pos[2] - atom_pos[1];
    
    v1 = atom_pos[0] - atom_pos[1];
    
    v2 = atom_pos[3] - atom_pos[2];

    v1.orthogonalize(n);
    
    v2.orthogonalize(n);
    
    dtemp = v1.vlength() * v2.vlength();
    
    if(dtemp < 1.e-8)
      //
      return 0.;

    dtemp = vdot(v1, v2) / dtemp;
    
    if(dtemp < -1.)
      //
      dtemp = -1.;
  
    if(dtemp > 1.)
      //
      dtemp = 1.;

    dtemp = std::acos(dtemp);

    if(volume(n, v1, v2) > 0.) {
      //
      return dtemp;
    }
    else
      //
      return 2. * M_PI - dtemp;

    // wrong number of atoms
    //
  default:

    std::cerr << funame << "wrong number of atoms: " << atom_pos.size() << "\n";

    throw Error::Range();
  }
}

/*************************************************************************************
 ****************************** CARTESIAN COORDINATES ********************************
 *************************************************************************************/

void Coord::Cartesian::_assert (int s) const
{
  if(s % 3) {
    //
    std::cerr << "Coord::Cartesian:_assert: wrong dimension: " << s <<"\n";

    throw Error::Range();
  }
}

/*********************************************************************************************************
 ***************************************** Z-MATRIX DEFINITION ******************************************
 *********************************************************************************************************/

void Coord::ZMat::init (std::istream& from)
{
  const char funame [] = "Coord::ZMat::init: ";

  if(isinit()) {
    //
    std::cerr << funame << "already initialized\n";

    throw Error::Init();
  }
  
  double dtemp;
  
  int    itemp;

  std::string stemp, comment, token;

  IO::String sym;

  std::map<std::string, int> sym_to_num;
  
  if(!(from >> itemp)) {
    //
    std::cerr << funame << "cannot read number of atoms\n";

    throw Error::Input();
  }

  std::getline(from, comment);
  
  if(itemp <= 0) {
    //
    std::cerr << funame << "number of atoms out of range: " << itemp << "\n";

    throw Error::Range();
  }

  resize(itemp);

  for(iterator rit = begin(); rit != end(); ++rit) {
    //
    const int curr = rit - begin();

    itemp = curr < 3 ? curr : 3;

    rit->resize(itemp);

    IO::LineInput lin(from);

    // atom description
    //
    if(!(lin >> token)) {
      //
      std::cerr << funame << "cannot read atom description\n";

      throw Error::Input();
    }

    if(!std::isalpha(token[0])) {
      //
      std::cerr << funame << curr + 1 << "-th atom description: should start with the atom name: " << token << "\n";

      throw Error::Init();
    }

    std::set<int> rpool;
    
    for(int v = 0; v < rit->size(); ++v) {
      //
      // read reference label
      //
      if(!(lin >> sym)) {
	//
	std::cerr << funame << "cannot read " << curr + 1 << "-th atom " << v + 1 << "-th reference label\n";

	throw Error::Input();
      }

      // reference is line number
      //
      if(std::isdigit(sym[0])) {
	 //
	rit->ref(v) = (int)sym - 1;
      }
      // reference is a symbol
      //
      else if(sym_to_num.find(sym) != sym_to_num.end()) {
	//
	rit->ref(v) = sym_to_num[sym];
      }
      else {
	//
	std::cerr << funame << "cannot find " << curr + 1 << "-th atom " << v + 1 << "-th reference label: " << sym << "\n";

	throw Error::Init();
      }

      if(rit->ref(v) < 0 || rit->ref(v) >= curr) {
	//
	std::cerr << funame << curr + 1 << "-th atom " << v + 1 << "-th reference out of range: " << rit->ref(v) + 1 << "\n";

	throw Error::Range();
      }

      if(!rpool.insert(rit->ref(v)).second) {
	//
	std::cerr << funame << curr + 1 << "-th atom " << v + 1 << "-th reference, " << rit->ref(v) + 1 << ": duplicated\n";

	throw Error::Init();
      }
	
      // read variable name
      //
      if(!(lin >> rit->var(v))) {
	//
	std::cerr << funame << "cannot read " << curr + 1 << "-th atom " << v + 1 << "-th variable name\n";

	throw Error::Input();
      }

      if(!std::isalpha(rit->var(v)[0])) {
	//
	std::cerr << funame << curr + 1 << "-th atom " << v + 1 << "-th variable name should start with letter: " << rit->var(v) << "\n";

	throw Error::Init();
      }

      _vmap_t::const_iterator vit = _var_to_ind.find(rit->var(v));
      
      if(vit != _var_to_ind.end()) {
	//
	std::cerr << funame << curr + 1 << "-th atom " << v + 1 << "-th variable name, " << rit->var(v)
		  << ", duplicated in "<< vit->second.second + 1 << " record, " << vit->second.first + 1 << " variable\n";

	throw Error::Init();
      }

      _var_to_ind[rit->var(v)] = std::make_pair(v, curr);
    }
    
    // check if there is an isotope modifier in the atom description
    //
    std::string::size_type ipos = token.find('.');

    int isot = -1;
    
    if(ipos < token.size()) {
      //
      ++ipos;
      
      if(ipos == token.size()) {
	//
	std::cerr << funame << curr + 1 << "-th atom description: dot position out of range: " << token << "\n";

	throw Error::Range();
      }
      
      isot = (int)IO::String(token.substr(ipos));

      if(isot < 0) {
	//
	std::cerr << funame << curr + 1 << "-th atom description: isotope index out of range: " << token << "\n";

	throw Error::Range();
      }

      token = token.substr(0, ipos - 1);
    }

    for(int c = 1; c < token.size(); ++c) {
      //
      if(std::isdigit(token[c])) {
	//
	if(sym_to_num.find(token) != sym_to_num.end()) {
	  //
	  std::cerr << funame << curr + 1 << "-th atom description: duplicated label: " << token << "\n";

	  throw Error::Init();
	}

	sym_to_num[token] = curr;

	token = token.substr(0, c);

	break;
      }
    }

    if(isot < 0) {
      //
      rit->atom().set(token);
    }
    else
      //
      rit->atom().set(token, isot);
    //
  }// cycle 
  //
}//

std::vector<int> Coord::ZMat::sign (const std::string& var) const
{
  const int vmax = _pos(var).first + 1;
 
  const int a = _pos(var).second;

  std::vector<int> res(vmax + 1, a);

  for(int v = 0; v < vmax; ++v)
    //
    res[v + 1] = at(a).ref(v);

  return res;
}

void Coord::ZMat::cm_shift (Cartesian& x) const
{
  const char funame [] = "Coord::ZMat::cm_shift: ";

  double dtemp;
  int    itemp;
  
  if(x.atom_size() != size()) {
    //
    std::cerr << funame << "sizes mismatch: " << x.atom_size() << " vs. " << size() << "\n";

    throw Error::Range();
  }
  
  for(int i = 0; i < 3; ++i) {
    //
    double mass = 0., shift = 0.;

    for(int a = 0; a < size(); ++a) {
      //
      dtemp = (*this)[a].atom().mass();

      mass += dtemp;

      shift += x.atom_pos(a)[i] * dtemp;
    }

    shift /= mass;
    
    for(int a = 0; a < size(); ++a)
      //
      x.atom_pos(a)[i] -= shift;
  }
}

/*********************************************************************************************************
 ***************************************** Z-MATRIX COORDINATES ******************************************
 *********************************************************************************************************/

Lapack::Vector Coord::ZData::inertia_moments () const
{
  const char funame [] = "Coord::ZData::inertia_moments: ";

  double dtemp;
  int    itemp;
  
  Cartesian x = *this;

  _base.cm_shift(x);

  // inertia matrix
  //
  Lapack::SymmetricMatrix res(3);

  dtemp = 0.;
  
  for(int a = 0; a < _base.size(); ++a)
    //
    dtemp += _base[a].atom().mass() * vdot(x.atom_pos(a), 3);

  
  res = dtemp;

  for(int a = 0; a < _base.size(); ++a)
    //
    for(int i = 0; i < 3; ++i)
      //
      for(int j = i; j < 3; ++j)
	//
	res(i, j) -= _base[a].atom().mass() * x.atom_pos(a)[i] * x.atom_pos(a)[j];

  return res.eigenvalues();
}

Lapack::SymmetricMatrix Coord::ZData::mobility_matrix () const
{
  const char funame [] = "Coord::ZData::mobility_matrix: ";

  double dtemp;

  int    itemp;
  
  if(atom_size() < 3) {
    //
    std::cerr << funame << "number of atoms out of range: " << atom_size() << "\n";

    throw Error::Range();
  }
  
  std::vector<Internal> internal;
  
  for(int a = 1; a < atom_size(); ++a) {
    //
    std::vector<int> sign(1, a);

    for(int i = 0; i < _base[a].size(); ++i) {
      //
      sign.push_back(_base[a].ref(i));

      internal.push_back(Internal(atom_size(), sign));
    }
  }
      
  Cartesian x = *this;

  // gradient matrix
  //
  Lapack::Matrix gm(internal.size(), x.size());

  for(int i = 0; i < internal.size(); ++i)
    //
    for(int j = 0; j < x.size(); ++j)
      //
      gm(i, j) = internal[i].grad(x, j);

  // generalized mobility matrix
  //
  Lapack::SymmetricMatrix res(internal.size());

  res = 0.;

  for(int i = 0; i < internal.size(); ++i)
    //
    for(int j = i; j < internal.size(); ++j)
      //
      for(int a = 0; a < atom_size(); ++a) {
	//
	dtemp = 0.;
	
	for(int c = 0; c < 3; ++c)
	  //
	  dtemp += gm(i, a * 3 + c) * gm(j, a * 3 + c);

	dtemp /= _base[a].atom().mass();

	res(i, j) += dtemp;
      }
  
  return res;
}

void Coord::ZData::_assert(int t, int a) const
{
  const char funame [] = "Coord::ZData::_assert: ";

  if(t < 0 || t > 2 || a <= 0 || a >= atom_size() || t >= a) {
    //
    std::cerr << funame << "indices out of range: " << t << ", " << a << "\n";

    throw Error::Range();
  }
}

void Coord::ZData::import (const Cartesian& cart)
{
  const char funame [] = "Coord::ZData::import: ";

  double dtemp;
  int    itemp;

  if(cart.atom_size() != atom_size()) {
    //
    std::cerr << funame << "sizes mismatch: " << cart.atom_size() << ", " << atom_size() << "\n";

    throw Error::Range();
  }

  for(int a = 1; a < atom_size(); ++a) {
    //
    std::vector<D3::Vector> atom_pos(1);
      
    atom_pos[0] = cart.atom_pos(a);

    int tmax =  a < 3 ? a : 3;
    
    for(int t = 0; t < tmax; ++t) {

      atom_pos.push_back((D3::Vector)cart.atom_pos(_base[a].ref(t)));

      (*this)(t, a) = internal(atom_pos);
    }
  }
}

Coord::ZData::operator Cartesian () const
{
  const char funame [] = "Coord::ZData::operator Cartesian: ";

  double dtemp;
  int    itemp;

  Cartesian res(atom_size() * 3);

  D3::Vector nx, ny, nz;
  
  double cx, cy, cz;

  for(int a = 0; a < atom_size(); ++a) {
    //
    switch(a) {
      //
    case 0:
      //
      for(int i = 0; i < 3; ++i)
	//
	res.atom_pos(a)[i] = 0.;

      break;

    case 1:
      //
      for(int i = 0; i < 3; ++i)
	//
	if(i) {
	  //
	  res.atom_pos(a)[i] = 0.;
	}
	else
	  //
	  res.atom_pos(a)[i] = (*this)(DISTANCE, a);
	  
      break;

    case 2:
      //
      dtemp = (*this)(DISTANCE, a) * std::cos((*this)(ANGLE, a));
      
      switch(_base[a].ref(0)) {
	//
      case 0:
	//
	res.atom_pos(a)[0] = dtemp;
	  
	break;

      case 1:
	//
	res.atom_pos(a)[0] = res.atom_pos(1)[0] - dtemp;

	break;

      default:
	//
	std::cerr << funame << "wrong reference atom: " << _base[a].ref(0) << "\n";
	
	throw Error::Logic();
      }
      
      res.atom_pos(a)[1] = (*this)(DISTANCE, a) * std::sin((*this)(ANGLE, a));

      res.atom_pos(a)[2] = 0.;

      break;

    default:
      // nx
      //
      nx  = (const double*)res.atom_pos(_base[a].ref(1));
      nx -= (const double*)res.atom_pos(_base[a].ref(0));
      
      nx.normalize();

      // ny
      //
      ny  = (const double*)res.atom_pos(_base[a].ref(2));
      ny -= (const double*)res.atom_pos(_base[a].ref(1));
      
      ny.orthogonalize(nx);

      ny.normalize();

      // nz
      //
      D3::vprod(nx, ny, nz);

      cx =  (*this)(DISTANCE, a) * std::cos((*this)(ANGLE, a));

      dtemp = std::sin((*this)(ANGLE, a));

      cy =  (*this)(DISTANCE, a) * dtemp * std::cos((*this)(DIHEDRAL, a));
      
      cz = -(*this)(DISTANCE, a) * dtemp * std::sin((*this)(DIHEDRAL, a));
      
      for(int i = 0; i < 3; ++i)
	//
	res.atom_pos(a)[i] = res.atom_pos(_base[a].ref(0))[i] + cx * nx[i] + cy * ny[i] + cz * nz[i];
      
      break;
    }
  }
  
  return res;
}

/***************************************************************************************************
 ***************************** ABSTRACT FUNCTION IN CARTESIAN SPACE ********************************
 ***************************************************************************************************/

// numerical derivative step
//
double Coord::CartFun::step = 1.e-2;

// gradient length
//
double Coord::CartFun::grad_length (const Cartesian& x) const
{
  double dtemp;
  
  int    itemp;
  
  double res = 0.;

  for(int i = 0; i < x.size(); ++i) {
    //
    dtemp = grad(x, i);

    res += dtemp * dtemp;
  }

  return std::sqrt(res);
}

// gradient component
//
double Coord::CartFun::grad (const Cartesian& x, int i) const
{
  _assert(i, x.size());
  
  if(!_depend(i))
    //
    return 0.;

  Cartesian y = x;

  y[i] += step;

  double res = evaluate(y);

  y[i] -= step + step;

  res -= evaluate(y);

  _adjust_dihedral(res);

  return res / 2. / step;
}

// hessian component
//
double Coord::CartFun::hess (const Cartesian& x, int i, int j) const
{
  _assert(i, x.size());

  _assert(j, x.size());
  
  if(!_depend(i) || !_depend(j))
    //
    return 0.;

  Cartesian y = x;

  double two = step + step;

  double res;
  
  if(i != j) {
    //
    y[i] += step;

    y[j] += step;
    
    res = evaluate(y);

    y[i] -= two;
    
    res -= evaluate(y);
    
    y[j] -= two;
    
    res += evaluate(y);

    y[i] += two;

    res -= evaluate(y);

    
    _adjust_dihedral(res);

    return res / two / two;
  }
  
  res = -2. * evaluate(y);
  
  y[i] += step;

  res += evaluate(y);

  y[i] -= two;

  res += evaluate(y);

  _adjust_dihedral(res);

  return res / step / step;
}

void Coord::CartFun::_adjust_dihedral (double& dv) const
{
  const char funame [] = "Coord::CartFun::_adjust_dihedral: ";

  double dtemp;
  
  if(type() != DIHEDRAL)
    //
    return;

  dtemp = std::nearbyint(dv / M_PI);

  if(dtemp == 0.)
    //
    return;

  dv -= M_PI * dtemp;
  
  IO::log << funame << "dihedral jump: " << std::setw(3) << dtemp << " dihedral value: " << dv << "\n";

}

/**********************************************************************************************************
 ************************* INTERNAL COORDINATE (DISTANCE, ANGLE, OR DIHEDRAL) *****************************
 **********************************************************************************************************/
void Coord::Internal::_assert ()
{
  const char funame [] = "Coord::Internal::_assert: ";

  if(size() < 2 || size() > 4) {
    //
    std::cerr << funame << "number of atomic indices out of range: " << size() << "\n";

    throw Error::Init();
  }

  _type = size() - 2;

  if(atom_size() < 2) {
    //
    std::cerr << funame << "atoms number out of range: " << atom_size() << "\n";

    throw Error::Range();
  }

  for(std::vector<int>::const_iterator cit = begin(); cit != end(); ++cit) {
    //
    if(!add_dep(*cit)) {
      //
      std::cerr << funame << "duplicated atomic index: " << *cit + 1 << "\n";

      throw Error::Init();
    }

    if(*cit < 0 || *cit >= atom_size()) {
      //
      std::cerr << funame << cit - begin() + 1 << "-th atomic index out of range: " << *cit << "\n";

      throw Error::Range();
    }
  }
}

Coord::Internal::Internal (int as, std::istream& from) : _atom_size(as)
{
  const char funame [] = "Coord::Internal::Internal: ";

  int    itemp;
  
  IO::LineInput lin(from);

  while(lin >> itemp)
    //
    push_back(--itemp);

  _assert();
}

double Coord::Internal::evaluate (const Cartesian& x) const
{
  const char funame [] = "Coord::Internal::evaluate: ";

  double dtemp;
  int    itemp;

  if(x.atom_size() != atom_size()) {
    //
    std::cerr << funame << "wrong number of atoms: " << x.atom_size() << "\n";

    throw Error::Range();
  }
    
  std::vector<D3::Vector> atom_pos(size());
  
  for(int i = 0; i < size(); ++i)
    //
    atom_pos[i] = x.atom_pos((*this)[i]);

  return internal(atom_pos);
}

