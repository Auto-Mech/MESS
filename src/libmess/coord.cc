#include "coord.hh"

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
 ***************************************** Z-MATRIX COORDINATES ******************************************
 *********************************************************************************************************/

void Coord::ZMat::_assert(int t, int a) const
{
  const char funame [] = "Coord::ZMat::_assert: ";

  if(t < 0 || t > 2 || a <= 0 || a >= _atom_size || t >= a) {
    //
    std::cerr << funame << "indices out of range: " << t << ", " << a << "\n";

    throw Error::Range();
  }
}

/*******************************************************************************************************
 ***************************************** Z-MATRIX SIGNATURE ******************************************
 *******************************************************************************************************/

void Coord::ZBase::_assert() const
{
  const char funame [] = "Coord::ZBase::_assert: ";

  int itemp;

  if(size() % 3) {
    //
    std::cerr << funame << "size should be multiple of 3\n";

    throw Error::Init();
  }

  for(int a = 1; a < size() / 3; ++a) {
    //
    const int tmax = a < 3 ? a : 3;

    std::set<int> pool;
    
    for(int t = 0; t < tmax; ++t) {
      //
      itemp = (*this)[3 * a + t];
      
      if(itemp < 0 || itemp >= a) {
	//
	std::cerr << funame << t + 1 << "-th reference atom index for " << a + 1 << "-th atom out of range: " << itemp + 1 << "\n";

	throw Error::Range();
      }

      if(!pool.insert(itemp).second) {
	//
	std::cerr << funame << t + 1 << "-th reference atom index for " << a + 1 << "-th atom is already in use (duplicated): " << itemp + 1 << "\n";

	throw Error::Init();
      }
    }
  }
}

void Coord::ZBase::_assert(int t, int a) const
{
  const char funame [] = "Coord::ZBase::_assert: ";

  if(t < 0 || t > 2 || a <= 0 || a >= atom_size() || t >= a) {
    //
    std::cerr << funame << "indices out of range: " << t << ", " << a << "\n";

    throw Error::Range();
  }
}

void Coord::ZBase::init(std::istream& from)
{
  const char funame [] = "Coord::ZBase::init: ";

  int itemp;

  if(!(from >> itemp)) {
    //
    std::cerr << funame << "cannot read atoms #\n";

    throw Error::Input();
  }

  if(itemp <= 0) {
    //
    std::cerr << funame << "atoms # out of range: " << itemp << "\n";

    throw Error::Range();
  }

  resize(3 * itemp);

  const int amax = itemp;

  for(int a = 0; a < amax; ++a) {
    //
    IO::LineInput lin(from);

    if(!(lin >> itemp)) {
      //
      std::cerr << funame << "cannot read " << a + 1 << "-th atom index\n";

      throw Error::Input();
    }

    if(itemp != a + 1) {
      //
      std::cerr << funame << "wrong atom index: " << itemp << ": should be " << a + 1 << "\n";

      throw Error::Init();
    }
    
    const int tmax = a < 3 ? a : 3;

    std::set<int> tpool;
    
    for(int t = 0; t < tmax; ++t) {
      //
      if(!(lin >> itemp)) {
	//
	std::cerr << funame << "cannot read " << t + 1 << "-th reference index for " << a + 1 << "-th atom\n";

	throw Error::Input();
      }

      if(itemp < 1 || itemp > a) {
	//
	std::cerr << funame << t + 1 << "-th reference index for " << a + 1 << "-th atom out of range: " << itemp << "\n";

	throw Error::Input();
      }

      if(!tpool.insert(itemp).second) {
	//
	//
	std::cerr << funame << t + 1 << "-th reference index for " << a + 1 << "-th atom is already in use (duplicated): " << itemp << "\n";

	throw Error::Init();
      }

      (*this)[3 * a + t] = itemp - 1;
    }
  }
}

Coord::ZMat Coord::ZBase::operator() (const Cartesian& cart) const
{
  const char funame [] = "Coord::ZBase::operator(): ";

  double dtemp;
  int    itemp;

  if(cart.atom_size() != atom_size()) {
    //
    std::cerr << funame << "sizes mismatch: " << cart.atom_size() << ", " << atom_size() << "\n";

    throw Error::Range();
  }

  ZMat res(atom_size());

  for(int a = 1; a < atom_size(); ++a) {
    //
    std::vector<D3::Vector> atom_pos(1);
      
    atom_pos[0] = cart.atom_pos(a);

    int tmax =  a < 3 ? a : 3;
    
    for(int t = 0; t < tmax; ++t) {

      atom_pos.push_back((D3::Vector)cart.atom_pos((*this)(t, a)));

      res(t, a) = internal(atom_pos);
    }
  }

  return res;
}

Coord::Cartesian Coord::ZBase::operator() (const ZMat& zmat) const
{
  const char funame [] = "Coord::ZBase::operator(): ";

  double dtemp;
  int    itemp;

  if(zmat.atom_size() != atom_size()) {
    //
    std::cerr << funame << "sizes mismatch: " << zmat.atom_size() << ", " << atom_size() << "\n";

    throw Error::Range();
  }

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
	  res.atom_pos(a)[i] = zmat(0, a);
	  
      break;

    case 2:
      //
      dtemp = zmat(0, a) * std::cos(zmat(1, a));
      
      switch((*this)(0, a)) {
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
	std::cerr << funame << "wrong reference atom: " << (*this)(0, a) << "\n";
	
	throw Error::Logic();
      }
      
      res.atom_pos(a)[1] = zmat(0, a) * std::sin(zmat(1, a));

      res.atom_pos(a)[2] = 0.;

      break;

    default:
      // nx
      //
      nx  = (const double*)res.atom_pos((*this)(1, a));
      nx -= (const double*)res.atom_pos((*this)(0, a));
      
      nx.normalize();

      // ny
      //
      ny  = (const double*)res.atom_pos((*this)(2, a));
      ny -= (const double*)res.atom_pos((*this)(1, a));
      
      ny.orthogonalize(nx);

      ny.normalize();

      // nz
      //
      D3::vprod(nx, ny, nz);

      cx = zmat(0, a) * std::cos(zmat(1, a));

      dtemp = std::sin(zmat(1, a));

      cy =  zmat(0, a) * dtemp * std::cos(zmat(2, a));
      
      cz = -zmat(0, a) * dtemp * std::sin(zmat(2, a));
      
      for(int i = 0; i < 3; ++i)
	//
	res.atom_pos(a)[i] = res.atom_pos((*this)(0, a))[i] + cx * nx[i] + cy * ny[i] + cz * nz[i];
      
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
  
    return res / two / two;
  }
  
  res = -2. * evaluate(y);
  
  y[i] += step;

  res += evaluate(y);

  y[i] -= two;

  res += evaluate(y);

  return res / step / step;
}
