/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#include "dynamic.hh"
#include "units.hh"

#include <cmath>

// maximal allowed angular vector length deviation from unity
//
double Dynamic::Coordinates::ang_len_tol = 2.;

// configuration region which is not sampled 
//
Dynamic::CCP Dynamic::exclude_region;

// Calculates the dynamic variables derivatives;  the torque is assumed 
// to be in laboratory frame for linear fragments and in molecular frame 
// for nonlinear ones

void Dynamic::set_dvd (const D3::Vector* torque, const double* dv, double* dvd)
{
  const char funame [] = "Dynamic::set_dvd: ";

  const D3::Vector& force = *torque;
  torque += 1;

  const double* pos = dv;
  const double* vel = dv + Structure::pos_size();

  double* pos_drv = dvd;
  double* vel_drv = dvd + Structure::pos_size();

  const double* orb_pos = pos + Structure::orb_pos();
  const double* orb_vel = vel + Structure::orb_vel();

  double* orb_pos_drv = pos_drv + Structure::orb_pos();
  double* orb_vel_drv = vel_drv + Structure::orb_vel();

  const double*  ang_pos [2];
  const double*  ang_vel [2];

  double* ang_pos_drv [2];
  double* ang_vel_drv [2];

  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::type(frag) != Molecule::MONOATOMIC) {
      //
      ang_pos [frag] = pos + Structure::ang_pos(frag);
      ang_vel [frag] = vel + Structure::ang_vel(frag);

      ang_pos_drv [frag] = pos_drv + Structure::ang_pos(frag);
      ang_vel_drv [frag] = vel_drv + Structure::ang_vel(frag);
    }
    else {
      //
      ang_pos[frag] = ang_vel[frag] = ang_pos_drv[frag] = ang_vel_drv[frag] = 0;
    }

  // orbital motion
  //
  for(int i = 0; i < 3; ++i) {
    //
    orb_pos_drv[i] = orb_vel[i];
    orb_vel_drv[i] = force[i] / Structure::mass();
  }

  // Fragments rotations
  // for linear fragments lab frame torque is assumed
  // for nonlinear fragments molecular frame torque is assumed

  for(int frag = 0; frag < 2; ++frag)
    //
    Structure::fragment(frag).set_dvd(torque[frag], ang_pos[frag], ang_vel[frag], ang_pos_drv[frag], ang_vel_drv[frag]);
} 

/*********************************************************************************************
 ****************************** DYNAMIC COORDINATES METHODS **********************************
 *********************************************************************************************/

double Dynamic::Coordinates::min_atom_dist = 1.5;

void Dynamic::Coordinates::put (double* dv) const
{
  for(int i = 0; i < size(); ++i)
    //
    dv[i] = (*this)[i];
}

void Dynamic::Coordinates::get (const double* dv, double* len) 
{
  Array<double>::operator=(dv);

  for(int f = 0; f < 2; ++f)
    //
    if(len) {
      //
      len[f] = normalize(f);
    }
    else
      //
      normalize(f);

  for(int f = 0; f < 2; ++f)
    //
    if(Structure::type(f) == Molecule::NONLINEAR)
      //
      _need_mfo_update[f] = true;
}

double Dynamic::Coordinates::get_ang_pos (int f, const double* pos) 
{
  const char funame [] = "Dynamic::Coordinates::get_ang_pos: ";

  if(Structure::type(f) == Molecule::MONOATOMIC) {
    //
    std::cerr << funame << "should not be called for monoatomic molecule\n";

    throw Error::Logic();
  }

  for(int i = 0; i < Structure::pos_size(f); ++i)
    //
    write_ang_pos(f)[i] = pos[i];

  if(Structure::type(f) == Molecule::NONLINEAR)
    //
    _need_mfo_update[f] = true;

  return normalize(f);
}

// normalize fragment's angular vector
//
double Dynamic::Coordinates::normalize (int frag) 
{
  const char funame [] = "Dynamic::Coordinates::normalize: ";

  if(Structure::type(frag) == Molecule::MONOATOMIC) {
    //
    return -1.;
  }

  return ::normalize(write_ang_pos(frag), Structure::pos_size(frag));
}

void Dynamic::Coordinates::mf2lf (int f, const double* mf, double* lf) const 
{
  const char funame [] = "Dynamic::Coordinates::mf2lf: ";

  _update_mfo(f);

  switch(Structure::type(f)) {
    //
  case Molecule::MONOATOMIC:
    //
    for(int i = 0; i < 3; ++i) {
      //
      lf[i] = 0.;

      if(mf[i] != 0.)
	//
	std::cerr << funame << "WARNING: non-zero molecular frame vector for a monoatomic fragment: " << mf[i] << "\n";
    }

    break;
    //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      lf[i] = mf[0] * ang_pos(f)[i];

    break;
    //
  case Molecule::NONLINEAR:
    //
    D3::vprod(mf, _mfo[f], lf);
    
    break;
    //
  default:
    //
    std::cerr << funame << "wrong case\n";

    throw Error::Logic();
  }
}

void Dynamic::Coordinates::lf2mf (int f, const double* lf, double* mf) const 
{
  const char funame [] = "Dynamic::Coordinates::lf2mf: ";

  double dtemp;

  _update_mfo(f);

  switch(Structure::type(f)) {
    //
  case Molecule::MONOATOMIC:
    //
    for(int i = 0; i < 3; ++i) {
      //
      mf[i] = 0.;

      if(lf[i] != 0.)
	//
	std::cerr << funame << "WARNING: non-zero laboratory frame vector for a monoatomic fragment: " << lf[i] << "\n";
    }

    break;
    //
  case Molecule::LINEAR:
    //
    mf[0] = vdot(lf, ang_pos(f), 3);
    
    mf[1] = mf[2] = 0.;

#ifdef DEBUG
    
    dtemp = 1. - mf[0] * mf[0] / ::vdot(lf, 3) / ::vdot(ang_pos(f), 3);

    if(dtemp > 1.e-7)
      //
      std::cerr << funame << "WARNING: lab frame vector and the orientational orth are not colinear\n";

#endif
    
    break;
    //
  case Molecule::NONLINEAR:
    //
    D3::vprod(_mfo[f], lf, mf);

    break;
  }
}

void Dynamic::Coordinates::set_rel_pos(int f, std::vector<D3::Vector>& rel_pos) const
{
  const char funame [] = "Dynamic::Coordinates::set_rel_pos: ";

  rel_pos.resize(Structure::fragment(f).size());
  
  for(int a = 0; a < Structure::fragment(f).size(); ++a)
    //
    mf2lf(f, Structure::fragment(f)[a], rel_pos[a]);
}

bool Dynamic::Coordinates::are_atoms_close () const
{
  double dtemp, dist;

  std::vector<D3::Vector> rel_pos [2];

  for(int f = 0; f < 2; ++f)
    //
    set_rel_pos(f, rel_pos[f]);
  
  for(int a0 = 0; a0 < Structure::fragment(0).size(); ++a0)
    //
    for(int a1 = 0; a1 < Structure::fragment(1).size(); ++a1) {
      //
      dist = 0.;

      for(int i = 0; i < 3; ++i) {
	//
	dtemp = rel_pos[1][a1][i] + orb_pos()[i] - rel_pos[0][a0][i];

	dist += dtemp * dtemp;
      }

      dist = std::sqrt(dist);

      if(dist < min_atom_dist)
	//
	return true;
    }

  return false;
}

double Dynamic::Coordinates::_imm (int f, int i, int j) const 
{
  const char funame [] = "Dynamic::Coordinates::_imm: ";

  double dtemp;

  switch(Structure::type(f)) {
    //
  case Molecule::MONOATOMIC:
    //
    return 0.;

  case Molecule::LINEAR:
    //
    if(i != j) {
      //
      return -Structure::fragment(f).imom(0) * ang_pos(f)[i] * ang_pos(f)[j];
    }
    else
      //
      return  Structure::fragment(f).imom(0) * (1. - ang_pos(f)[i] * ang_pos(f)[j]);

  case Molecule::NONLINEAR:
    //
    dtemp = 0.;

    for(int k = 0; k < 3; ++k) {
      //
      dtemp += Structure::fragment(f).imom(k) * _mfo[f](k, i) * _mfo[f](k, j);
    }

    return dtemp;
    //
  default:
    //
    std::cerr << funame << "wrong case\n";

    throw Error::Logic();
  }  
}

Lapack::SymmetricMatrix Dynamic::Coordinates::imm () const
{
  for(int f = 0; f < 2; ++f)
    //
    _update_mfo(f);

  double r2 = vdot(orb_pos(), 3);

  Lapack::SymmetricMatrix res(3);
  
  for(int i = 0; i < 3; ++i)
    //
    for(int j = i; j < 3; ++j) {
      //
      res(i, j) = -Structure::mass() * orb_pos()[i] * orb_pos()[j];

      if(i == j)
	//
	res(i, j) += Structure::mass() * r2;
      
      for(int f = 0; f < 2; ++f)
	//
	if(Structure::type(f) != Molecule::MONOATOMIC)
	  //
	  res(i, j) += _imm(f, i, j);
    }

#ifdef DEBUG
  
  Lapack::Vector ev = res.eigenvalues();
  
#endif
   
  return res;
}

void Dynamic::Coordinates::print_geom (std::ostream& to, const std::string& indent) const
{
  double     dtemp;
  int        itemp;
  D3::Vector vtemp;
  
  //int old_prec = to.precision(6);

  const int prec = to.precision();

  std::vector<D3::Vector> rel_pos [2];

  to << indent << "Geometry, Angstrom:\n";
  
  for(int f = 0; f < 2; ++f) {
    //
    set_rel_pos(f, rel_pos[f]);

    if(!f) {
      //
      dtemp = -1.;
    }
    else
      //
      dtemp = 1.;
  
    for(int at = 0; at < rel_pos[f].size(); ++at) {
      //
      to << indent << "   " << std::setw(3) << Structure::fragment(f)[at].name();
      
      for(int i = 0; i < 3; ++i)
	//
	to << std::setw(prec + 7) << (rel_pos[f][at][i] + dtemp * Structure::mass() / Structure::fragment(f).mass() * orb_pos()[i]) 
	//
	/ Phys_const::angstrom;
      
      to << "\n";
    }
  }
  
  for(int f = 0; f < 2; ++f)
    //
    if(Structure::type(f) != Molecule::MONOATOMIC) {
      //
      lf2mf(f, orb_pos(), (double*)vtemp);
      
      to << indent << "cm position of the " << 1 - f << "-th fragment in the framework of\n"
	//
	 << indent << "the " << f << "-th fragment, Bohr:";

      for(int i = 0; i < 3; ++i) {
	//
	dtemp = vtemp[i];
	
	if(f)
	  //
	  dtemp = -dtemp;

	to << std::setw(prec + 7) << dtemp;
      }
      
      to << "\n";
    }

  double dist;
  
  to << indent << "Distance matrix, Angstrom:\n";
  
  to << indent << "   " << std::setw(3) << "0\\1";
  
  for(int at1 = 0; at1 < Structure::fragment(1).size(); ++at1)
    //
    to << std::setw(prec + 7) << Structure::fragment(1)[at1].name();
  
  to << "\n";

  for(int at0 = 0; at0 < rel_pos[0].size(); ++at0) {
    //
    to << indent << "   " << std::setw(3) << Structure::fragment(0)[at0].name();
    
    for(int at1 = 0; at1 < rel_pos[1].size(); ++at1) {
      //
      dist = 0.;
      
      for(int i = 0; i < 3; ++i) {
	//
	dtemp = rel_pos[1][at1][i] + orb_pos()[i] - rel_pos[0][at0][i];
	
	dist += dtemp * dtemp;
      }
      
      to << std::setw(prec + 7) << std::sqrt(dist) / Phys_const::angstrom;
    }
    
    to << "\n";
  }

//  to.precision(old_prec);
}

void Dynamic::Coordinates::dipole (int frag, double* d) const 
{
  static const char funame [] = "Dynamic::Coordinates::dipole";

  D3::Vector lfd;
  
  switch(Structure::fragment(frag).dipole().size()) {
    //
  case 3:
    //
    mf2lf(frag, Structure::fragment(frag).dipole_vector(), d);
    
    return;
    //
  case 1:
    //
    switch(Structure::fragment(frag).type()) {
      //
    case Molecule::LINEAR:
      //
      lfd = ang_pos(frag);
      lfd.normalize();
      lfd *= Structure::fragment(frag).dipole()[0];
      
      for(int i = 0; i < 3; ++i)
	//
	d[i] = lfd[i];
      
      return;
      //
    case Molecule::NONLINEAR:
      //
      switch(Structure::fragment(frag).top()) {
	//
      case Molecule::PROLATE:
	//
	lfd[0] = Structure::fragment(frag).dipole()[0];
	
	mf2lf(frag, lfd, d);
	
	return;
	//
      case Molecule::OBLATE:
	//
	lfd[2] = Structure::fragment(frag).dipole()[0];
	
	mf2lf(frag, lfd, d);
	
	return;
	//
      default:
	//
	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
      
      return;
      //
    default:
      //
      std::cerr << funame << "should not be here\n";
      throw Error::Logic();
    }
    
    return;
    //
  default:
    //
    std::cerr << funame << "should not be here\n";
    throw Error::Logic();
  }
}

void Dynamic::Coordinates::quadrupole_vector_product (int frag, const double* v, double* q) const 
{
  static const char funame [] = "Dynamic::Coordinates::quadrupole_vector_product";

  D3::Vector mfv, mfq;
  switch(Structure::fragment(frag).quadrupole().size()) {
    //
  case 6:
    //
    lf2mf(frag, v, mfv);
    
    mfq = Structure::fragment(frag).quadrupole_matrix() * mfv;
    
    mf2lf(frag, mfq, q);
    
    return;
    //
  case 1:
    //
    switch(Structure::fragment(frag).type()) {
      //
    case Molecule::LINEAR:
      //
      mfq  = ang_pos(frag);
      mfq *= 3. * vdot(mfq, v) / vdot(mfq);
      mfq -= v;
      mfq *= Structure::fragment(frag).quadrupole()[0] / 2.;
      
      for(int i = 0; i < 3; ++i)
	//
	q[i] = mfq[i];
      
      return;
      //
    case Molecule::NONLINEAR:
      //
      switch(Structure::fragment(frag).top()) {
	//
      case Molecule::PROLATE:
	//
	lf2mf(frag, v, mfq);
	
	mfq[0] *=  Structure::fragment(frag).quadrupole()[0];
	mfq[1] *= -Structure::fragment(frag).quadrupole()[0] / 2.;
	mfq[2] *= -Structure::fragment(frag).quadrupole()[0] / 2.;
	
	mf2lf(frag, mfq, q);
	
	return;
	//
      case Molecule::OBLATE:
	//
	lf2mf(frag, v, mfq);
	
	mfq[0] *= -Structure::fragment(frag).quadrupole()[0] / 2.;
	mfq[1] *= -Structure::fragment(frag).quadrupole()[0] / 2.;
	mfq[2] *=  Structure::fragment(frag).quadrupole()[0];
	
	mf2lf(frag, mfq, q);
	
	return;
	//
      default:
	//
	std::cerr << funame << "should not be here\n";
	throw Error::Logic();
      }
      
      return;
      //
    default:
      //
      std::cerr << funame << "should not be here\n";
      throw Error::Logic();
    }
    return;

  default:
    //
    std::cerr << funame << "should not be here\n";
    throw Error::Logic();
  }

}

void Dynamic::Coordinates::polarizability_vector_product (int frag, const double* v, double* p) const 
{
  static const char funame [] = "Dynamic::Coordinates::polarizability_vector_product";

  D3::Vector mfv, mfp;
  
  switch(Structure::fragment(frag).polarizability().size()) {
    //
  case 6:
    //
    lf2mf(frag, v, mfv);
    
    mfp = Structure::fragment(frag).polarizability_matrix() * mfv;
    
    mf2lf(frag, mfp, p);
    
    return;
    //
  case 1:
    //
    for(int i = 0; i < 3; ++i)
      //
      p[i] = Structure::fragment(frag).polarizability()[0] * v[i];
    
    return;
    //
  case 2:

    switch(Structure::fragment(frag).type()) {
      //
    case Molecule::LINEAR:
      
      mfp  = ang_pos(frag);
      mfp *= 3. * vdot(mfp, v) / vdot(mfp);
      mfp -= v;
      mfp *= Structure::fragment(frag).polarizability()[1] / 2.;
      
      for(int i = 0; i < 3; ++i)
	//
	p[i] = mfp[i] + Structure::fragment(frag).polarizability()[0] * v[i];
      
      return;
      //
    case Molecule::NONLINEAR:
      //
      switch(Structure::fragment(frag).top()) {
	//
      case Molecule::PROLATE:
	//
	lf2mf(frag, v, mfp);
	
	mfp[0] *= Structure::fragment(frag).polarizability()[0] + Structure::fragment(frag).polarizability()[1];
	mfp[1] *= Structure::fragment(frag).polarizability()[0] - Structure::fragment(frag).polarizability()[1] / 2.;
	mfp[2] *= Structure::fragment(frag).polarizability()[0] - Structure::fragment(frag).polarizability()[1] / 2.;
	
	mf2lf(frag, mfp, p);
	
	return;
	//
      case Molecule::OBLATE:
	//
	lf2mf(frag, v, mfp);
	
	mfp[2] *= Structure::fragment(frag).polarizability()[0] + Structure::fragment(frag).polarizability()[1];
	mfp[1] *= Structure::fragment(frag).polarizability()[0] - Structure::fragment(frag).polarizability()[1] / 2.;
	mfp[0] *= Structure::fragment(frag).polarizability()[0] - Structure::fragment(frag).polarizability()[1] / 2.;
	
	mf2lf(frag, mfp, p);
	
	return;
	//
      default:
	//
	std::cerr << funame << "should not be here\n";
	
	throw Error::Logic();
      }
      
      return;
      //
    default:
      //
      std::cerr << funame << "should not be here\n";
      
      throw Error::Logic();
    }
    
    return;
    //
  default:
    //
    std::cerr << funame << "should not be here\n";
    
    throw Error::Logic();
  }
}

/************************************************************************
 ************************ CONFIGURATIONAL OBSERVABLES *******************
 ************************************************************************/

D3::Vector Dynamic::ZmatObserv::_vector (int a0, int a1, const Coordinates& dc, const std::vector<D3::Vector> rel_pos [])
{
  D3::Vector res = rel_pos[a1 % 2][a1 / 2] - rel_pos[a0 % 2][a0 / 2];

  if(a1 % 2 && !(a0 % 2)) {
    //
    res += dc.orb_pos();
  }
  else if(a0 % 2 && !(a1 % 2))
    //
    res -= dc.orb_pos();

  return res;
}

Dynamic::ZmatObserv::ZmatObserv (const std::vector<int>& ai) : _atom(ai)
{
  const char funame [] = "Dynamic::ZmatObserv::ZmatObserv: ";

  int itemp;

  if(_atom.size() > 4 || _atom.size() < 2) {
    //
    std::cerr << funame << "number of atoms out of range: " << _atom.size()
      //
	      << ": possible values: 2 (distance), 3 (angle), or 4 (dihedral angle)\n";
    
    throw Error::Range();
  }

  std::set<int> frag_pool;
  std::set<int> atom_pool;

  for(int i = 0; i < _atom.size(); ++i) {
    //
    itemp = _atom[i];
    
    if(itemp < 0 || itemp / 2 >= Structure::fragment(itemp % 2).size()) {
      //
      std::cerr << funame << i + 1 << "-th atomic index out of range: " << itemp << "\n";
      
      throw Error::Range();
    }

    if(!atom_pool.insert(itemp).second) {
      //
      std::cerr << funame << "there are identical atoms: " << itemp << "\n";
      
      throw Error::Input();
    }

    frag_pool.insert(itemp % 2);
  }

  if(frag_pool.size() == 1) {
    //
    itemp = *frag_pool.begin();
    
    std::cerr << funame << "all atoms are associated with the ";
    
    if(itemp == 0) {
      //
      std::cerr << "first";
    }
    else
      //
      std::cerr << "second";
    
    std::cerr << " fragment\n";

    throw Error::Range();
  }
}

double Dynamic::ZmatObserv::evaluate (const Coordinates& dc) const
{
  double dtemp;
  int    itemp;

  std::vector<D3::Vector> rel_pos [2];

  for(int f = 0; f < 2; ++f)
    //
    dc.set_rel_pos(f, rel_pos[f]);
  
  
  // distance
  //
  if(_atom.size() == 2) {
    //
    return _vector(_atom[0], _atom[1], dc, rel_pos).vlength();
  }
  // angle
  //
  else if(_atom.size() == 3) {
    //
    const double r0 = _vector(_atom[1], _atom[2], dc, rel_pos).vlength();
    const double r1 = _vector(_atom[2], _atom[0], dc, rel_pos).vlength();
    const double r2 = _vector(_atom[0], _atom[1], dc, rel_pos).vlength();

    return std::acos((r0 * r0 + r2 * r2 - r1 * r1) / 2. / r0 / r2);
  }
  // dihedral angle
  //
  else {
    //
    D3::Vector a1 = _vector(_atom[1], _atom[0], dc, rel_pos);
    D3::Vector b  = _vector(_atom[1], _atom[2], dc, rel_pos);
    D3::Vector a2 = _vector(_atom[2], _atom[3], dc, rel_pos);

    b.normalize();

    a1.orthogonalize(b);
    a2.orthogonalize(b);

    a1.normalize();
    a2.normalize();

    dtemp = std::acos(::vdot(a1, a2));
    
    if(D3::volume(b, a1, a2) < 0.)
      //
      dtemp = -dtemp;

    return dtemp;
  }
}

Dynamic::ObservCondition::ObservCondition (ConstSharedPointer<ObservBase> p, double limit, int m)
  //
  : _observ(p), _mode(m), _range(std::make_pair(limit, limit))
{
  const char funame [] = "Dynamic::ObservCondition::ObservCondition: ";

  if(!_observ) {
    //
    std::cerr << funame << "dynamic observable not initialized\n";
    
    throw Error::Init();
  }

  if(_mode != UPPER_BOUND && _mode != LOWER_BOUND) {
    //
    std::cerr << funame << "wrong mode\n";
    
    throw Error::Range();
  }
}

Dynamic::ObservCondition::ObservCondition (ConstSharedPointer<ObservBase> p, double lower, double upper)
  //
  : _observ(p), _mode(RANGE_BOUND), _range(std::make_pair(lower, upper))
{
  const char funame [] = "Dynamic::ObservCondition::ObservCondition: ";

  if(!_observ) {
    //
    std::cerr << funame << "dynamic observable not initialized\n";
    
    throw Error::Init();
  }

  if(lower >= upper) {
    //
    std::cerr << funame << "lower bound is larger than upper bound\n";
    
    throw Error::Range();
  }
} 

bool Dynamic::ObservCondition::test(const Coordinates& dc) const
{
  const char funame [] = "Dynamic::ObservCondition::test: ";

  double param_value = _observ->evaluate(dc);

  bool res;
  
  switch(_mode) {
    //
  case LOWER_BOUND:
    //
    if(param_value > _range.first)
      //
      return true;
    
    return false;
    //
  case UPPER_BOUND:
    //
    if(param_value < _range.second)
      //
      return true;
    
    return false;
    //
  case RANGE_BOUND:
    //
    if(param_value < _range.second && param_value > _range.first)
      //
      return true;
    
    return false;
    //
  default:
    //
    std::cerr << funame << "unknown mode\n";
    
    throw Error::Logic();
  }
}

/*********************************************************************
 ******************** DYNAMIC MOMENTA METHODS ************************
 **********************************************************************/

void Dynamic::Momenta::_init ()
{
  _orb_vel = begin() + Structure::orb_vel();
  
  for(int frag = 0; frag < 2; ++frag)
    //
    _ang_vel[frag] = begin() + Structure::ang_vel(frag);
}

void Dynamic::Momenta::put (double* dv) const
{
  for(int i = 0; i < size(); ++i)
    //
    dv[i] = (*this)[i];
}

double Dynamic::Momenta::orbital_kinetic_energy () const
{
    return Structure::mass() * vdot(orb_vel(), 3) / 2.;
}

double Dynamic::Momenta::total_kinetic_energy () const
{
    double res = orbital_kinetic_energy();

    for(int frag = 0; frag < 2; ++frag)
 	// 
	res += fragment_rotational_energy(frag);

    return res;
}

double Dynamic::Momenta::fragment_rotational_energy (int frag) const
{
  double dtemp;

  if(Structure::type(frag) == Molecule::MONOATOMIC)
    //
    return 0.;
    
  double res = 0.;

  for(int i = 0; i < 3; ++i) {
    //
    dtemp = ang_vel(frag)[i];

    res += Structure::fragment(frag).imom(i) * dtemp * dtemp;
  }

  return res / 2.;
}

/******************************************************************
 ******************* DYNAMIC VARIABLES METHODS ********************
 ******************************************************************/

double Dynamic::Vars::radial_kinetic_energy () const
{
  double vr = vdot(orb_pos(), orb_vel(), 3);
  
  return Structure::mass() * vr * vr  / vdot(orb_pos(), 3) / 2.;
}

void Dynamic::Vars::orbital_angular_momentum (double* am) const
{
    D3::vprod(orb_pos(), orb_vel(), am);
    
    for(int i = 0; i < 3; ++i)
      //
      am[i] *= Structure::mass();
}

void Dynamic::Vars::fragment_angular_momentum (int frag, double* am) const
{
  D3::Vector vtemp;

  switch(Structure::fragment(frag).type()) {
    //
  case Molecule::MONOATOMIC:
    //
    for(int i = 0; i < 3; ++i)
      //
      am[i] = 0.;
      
    return;
    //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      am[i] = Structure::fragment(frag).imom(0) * ang_vel(frag)[i];
      
    orthogonalize(am, ang_pos(frag), 3);
      
    return;
    //
  case Molecule::NONLINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      vtemp[i] = Structure::fragment(frag).imom(i) * ang_vel(frag)[i];
      
    mf2lf(frag, vtemp, am);
      
    return;
  }
}

void Dynamic::Vars::total_angular_momentum (double* am) const
{
  D3::Vector vtemp;

  orbital_angular_momentum(am);

  for(int frag = 0; frag < 2; ++frag) {
    //
    fragment_angular_momentum(frag, vtemp);

    for(int i = 0; i < 3; ++i)
      //
      am[i] += vtemp[i];
  }
}

double Dynamic::Vars::total_angular_momentum_length () const
{
  D3::Vector am;

  total_angular_momentum(am);

  return am.vlength();
}

double Dynamic::Vars::angular_momentum_k_projection (int frag) const
{
  D3::Vector am;

  fragment_angular_momentum(frag, am);

  return vdot(am, orb_pos()) / ::vlength(orb_pos(), 3);
}

double Dynamic::Vars::angular_momentum_m_projection (int frag) const 
{
  static const char funame [] = "Dynamic::Vars::angular_momentum_m_projection: ";

  if(Structure::type(frag) == Molecule::NONLINEAR)
    //
    switch(Structure::top(frag)) {
      //
    case Molecule::PROLATE:
      //
      return  Structure::fragment(frag).imom(0) * ang_vel(frag)[0];
      //
    case Molecule::OBLATE:
      //
      return  Structure::fragment(frag).imom(2) * ang_vel(frag)[2];
      //
    default:
      //
      break;
    }

  std::cerr << funame << "should be used for symmetric top only\n";

  throw Error::Logic();
}
