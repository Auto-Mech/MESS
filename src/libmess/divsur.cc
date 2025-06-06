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

#include "divsur.hh"
#include "random.hh"
#include "units.hh"
#include "key.hh"

#include <cmath>
#include <sstream>

/********************************************************************************
 *                                    Input/Output                              *
 ********************************************************************************/

void DivSur::FragData::read(int frag, std::istream& from) 
{
  const char funame [] = "DivSur::FragData::read: ";

  if(!Structure::isinit()) {
    //
    std::cerr << funame << "structure has not yet been initialized\n";
    
    throw Error::Init();
  }

  KeyGroup FragmentDataGroup;
  
  Key plane_key ("Plane");
  Key pivot_key ("Pivot");
  Key  cone_key ("Cone" );
  
  std::string token, comment, line, name;
  D3::Vector pos;
  
  double dtemp;
  int    itemp;
  
  clear();

  // fragment name should coincide with the molecular name
  //
  if(!(from >> name)) {
    //
    std::cerr << funame << frag <<"-th fragment: cannot read fragment's name\n";

    throw Error::Input();
  }

  std::getline(from, comment);

  if(name != Structure::fragment(frag).name()) {
    //
    std::cerr << funame << frag <<"-th fragment data section name, " << name
      //
	      << ", differs from the fragment's name, " << Structure::fragment(frag).name() <<"\n";

    throw Error::Input();
  }
    
  // add center of mass as a pivot point
  //
  pos = 0.;
  
  (*this)[name] = ConstSharedPointer<GeomObject>(new GeomPivot(pos));

  // read cycle
  //
  while(from >> token) {
    //
    if(token == IO::end_key()) {
      //
      return;
    }

    if(Structure::type(frag) == Molecule::MONOATOMIC) {
      //
      std::cerr << funame << "no explicit geometric objects for monoatomic fragments\n";
      
      throw Error::Form();
    }

    IO::LineInput lin(from);
    
    if(!(lin >> name)) {
      //
      std::cerr << funame << frag << "-th fragment: " << token << ": corrupted\n";
      
      throw Error::Form();
    }
    
    if(find(name) != end()) {
      //
      std::cerr << funame << frag << "-th fragment: " << name << " name already has been used\n";
      
      throw Error::Form();
    }

    std::vector<double> val;
    
    while(lin >> dtemp)
      //
      val.push_back(dtemp);

    // plane
    //
    if(token == plane_key) {
      //
      switch(Structure::type(frag)) {
	//
      case Molecule::LINEAR:
	//
	pos = 0.; 
	pos[0] = 1.;
	
	if(val.size() != 1) {
	  //
	  std::cerr << funame << frag << "-th fragment: "  << name << " plane: wrong number of arguments ," << val.size()
	    //
		    << ": should be 1: plane displacement along molecular axis\n";
	  
	  throw Error::Form();
	}
	
	(*this)[name] = ConstSharedPointer<GeomObject>(new GeomPlane(pos, val[0]));
	
	break;
	//
      case Molecule::NONLINEAR:
	//
	if(val.size() != 4) {
	  //
	  std::cerr << funame << frag << "-th fragment: " << name << " plane: wrong number of arguments, " << val.size()
	    //
		    << ": should be 4: plane normal and displacement\n";
	  
	  throw Error::Form();
	}
	
	for(int i = 0; i < 3; ++i)
	  //
	  pos[i] = val[i];
	
	(*this)[name] = ConstSharedPointer<GeomObject>(new GeomPlane(pos, val[3]));
      }
    }
    // pivot point
    //
    else if(token == pivot_key) {
      //
      switch(Structure::type(frag)) {
	//
      case Molecule::LINEAR:
	//
	if(val.size() != 1) {
	  //
	  std::cerr << funame << frag << "-th fragment: " << name << " pivot point: wrong number of arguments, " << val.size()
	    //
		    << ": should be 1: pivot point position along molecular frame\n";
	  
	  throw Error::Form();
	}
	
	pos = 0.; 
	pos[0] = val[0];
	
	break;
	//
      case Molecule::NONLINEAR:
	//
	if(val.size() == 1) {
	  //
	  // atomic index
	  //
	  itemp = std::floor(val[0] + 0.5);

	  if(itemp < 0 || itemp >= Structure::fragment(frag).size()) {
	    //
	    std::cerr << funame << frag << "-th fragment: " << name << " pivot point: atomic index out of range: " << itemp << "\n";

	    throw Error::Range();
	  }
	  
	  pos = Structure::fragment(frag)[itemp];
	}
	else if(val.size() == 3) {
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    pos[i] = val[i];
	}
	else {
	  //
	  std::cerr << funame << frag << "-th fragment: " << name << " pivot point: wrong number of arguments, " << val.size()
	    //
		    << ": should be either 1 (atom index) or 3 (vector)\n";
	  
	  throw Error::Form();
	}
      }
      
      // initialize the pivot point
      //
      (*this)[name] = ConstSharedPointer<GeomObject>(new GeomPivot(pos));
    }
    // cone
    //
    else if(token == cone_key) {
      //
      // initialize cone object
      //
      GeomCone* gcp = new GeomCone();

      // input order: head, axis, angle, radius
      //
      switch(Structure::type(frag)) {
	//
      case Molecule::LINEAR:
	//
	pos = 0.;
	
	if(val.size() == 4) {
	  //
	  // cone head position along molecular axis
	  //
	  pos[0] = val[0];

	  gcp->head = pos;

	  // cone axis direction along or opposite to molecular axis
	  //
	  if(val[1] < 0.) {
	    //
	    pos[0] = -1.;
	  }
	  else
	    //
	    pos[0] = 1.;

	  gcp->axis = pos;

	  // cone angle
	  //
	  gcp->proj = std::cos(val[2] * M_PI / 180.);

	  // sampling radius
	  //
	  gcp->radius = val[3];
	}
	else {
	  //
	  std::cerr << funame << frag << "-th fragment: " << name << " cone: wrong number of arguments, " << val.size()
	    //
		    << ": should be 4: head projection, axis direction, cone angle, and sampling radius\n";
	  
	  throw Error::Form();
	}
	
	break;
	//
      case Molecule::NONLINEAR:
	//
	if(val.size() == 6) {
	  //
	  // atomic index
	  //
	  itemp = std::floor(val[0] + 0.5);

	  if(itemp < 0 || itemp >= Structure::fragment(frag).size()) {
	    //
	    std::cerr << funame << frag << "-th fragment: " << name << " cone: atomic index out of range: " << itemp << "\n";

	    throw Error::Range();
	  }

	  // cone head
	  //
	  gcp->head = Structure::fragment(frag)[itemp];

	  // cone axis
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    gcp->axis[i] = val[i + 1];

	  gcp->axis.normalize();

	  // cone angle
	  //
	  gcp->proj = std::cos(val[4] * M_PI / 180.);
	  
	  // sampling radius
	  //
	  gcp->radius = val[5];
	}
	else if(val.size() == 8) {
	  //
	  // cone head
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    gcp->head[i] = val[i];

	  // cone axis
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    gcp->axis[i] = val[i + 3];

	  gcp->axis.normalize();

	  // cone angle
	  //
	  gcp->proj = std::cos(val[6] * M_PI / 180.);
	  
	  // sampling radius
	  //
	  gcp->radius = val[7];
	}
	else {
	  //
	  std::cerr << funame << frag << "-th fragment: " << name << " cone: wrong number of arguments, " << val.size() << "\n";
	  
	  throw Error::Form();
	}
      }
      
      (*this)[name] = ConstSharedPointer<GeomObject>(gcp);
    }
    // unknown keyword
    //
    else {
      //
      std::cerr << funame << "unknown geometric object: " << token << "\n";
      
      Key::show_all(std::cerr);
      std::cerr << "\n";
      
      throw Error::Form();
    }//
    //
  }// read cycle

  if(!from) {
    //
    std::cerr << funame << "input stream is corrupted\n";
    
    throw Error::Form();
  }
}

/*********************************************************************************
 **************************** PRIMITIVE SURFACE METHODS **************************
 *********************************************************************************/

double DivSur::Primitive::distance (const Dynamic::Coordinates& dc) const 
{
  const char funame [] = "DivSur::Primitive::distance: ";

  std::cerr << funame << "this is a fake function, should not be here\n";

  throw Error::Logic();
}

double DivSur::Primitive::_set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_norm) const
{
  const char funame [] = "DivSur::Primitive::_set_rc: ";

  std::cerr << funame << "this is a fake function, should not be here\n";

  throw Error::Logic();
}
 
double DivSur::Primitive::weight (const Dynamic::Coordinates& dc) const
{
  Dynamic::Momenta rc_norm;
  
  return _set_rc(dc, rc_norm);
}

// generalized velocities for the given temperature
//
double DivSur::Primitive::t_velocity (double t, Dynamic::Vars& dv) const
{
  Dynamic::Momenta rc_norm;

  double rc_len =  _set_rc(dv, rc_norm);
  
  _rc2tv(rc_norm, dv, t, dv);

  return rc_len;
}

// generalized velocities for the given energy
//
double DivSur::Primitive::e_velocity (double ener, Dynamic::Vars& dv) const
{
  Dynamic::Momenta rc_norm;
  
  double rc_len =  _set_rc(dv, rc_norm);
  double w      =  _rc2ev(rc_norm, dv, ener, dv);
  
  if(w <= 0.)
    //
    return -1.;
  
  return rc_len  * w;
}

// generalized velocities for the given energy and angular momentum
//
double DivSur::Primitive::j_velocity (double ener, double amom, Dynamic::Vars& dv) const
{
  Dynamic::Momenta rc_norm;
  
  double rc_len =  _set_rc(dv, rc_norm);
  double w      =  _rc2jv(rc_norm, dv, ener, amom, dv);
  
  if(w <= 0.)
    //
    return -1.;
  
  return rc_len  * w;
} 

// generate random generalized velocities from the reaction coordinate vector for the given temperature
//
void DivSur::Primitive::_rc2tv(const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, double t, Dynamic::Momenta& dm)
{
  const char funame [] = "DivSur::Primitive::_rc2tv: ";

  double dtemp;
    
  //#pragma omp critical (omp_random_section)
  //
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    //
    dm[i] = Random::norm();

  // for linear fragments remove the component colinear to the molecular axis
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::fragment(frag).type() == Molecule::LINEAR)
      //
      orthogonalize(dm.ang_vel(frag), dc.ang_pos(frag), 3);

  // adjust reactive component
  //
  orthogonalize(dm.begin(), rc_vec.begin(), Dynamic::Momenta::size());
  
  dtemp = std::sqrt(2. * Random::exp());

  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    //
    dm[i] += dtemp * rc_vec[i];

  // scale to the temperature
  //
  dm *= std::sqrt(t);

  /**************************** CONVERTING TO REGULAR VELOCITIES ***************************/

  // orbital motion
  //
  for (int i = 0; i < 3; ++i)
    //
    dm.orb_vel()[i] /= Structure::mass_sqrt();

  // fragments rotational motion
  //
  for(int f = 0; f < 2; ++f)
    //
    if(Structure::type(f) != Molecule::MONOATOMIC)
      //
      for (int i = 0; i < 3; ++i)
	//
	dm.ang_vel(f)[i] /= Structure::fragment(f).imom_sqrt(i);
}

// generate random generalized velocities from the reaction coordinate vector for the given energy
//
double DivSur::Primitive::_rc2ev(const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, double ener, Dynamic::Momenta& dm)
{
  const char funame [] = "DivSur::Primitive::_rc2ev: ";

  if(ener <= 0.)
    //
    return -1.;

  double dtemp;

  Random::orient(dm, Dynamic::Momenta::size());
  
  // for linear fragments remove the component colinear to the molecular axis
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::fragment(frag).type() == Molecule::LINEAR)
      //
      orthogonalize(dm.ang_vel(frag), dc.ang_pos(frag), 3);

  // remove reaction coordinate
  //
  orthogonalize(dm.begin(), rc_vec.begin(), Dynamic::Momenta::size());

  double nr_ener = ener * std::pow(Random::flat(), 2. / double(Structure::tm_dof() - 1)); // non-reactive modes
  double rc_ener = ener - nr_ener;         // energy of non-reactive modes
  double nr_len = std::sqrt(2. * nr_ener); // amplitude of non-reactive modes
  double rc_len = std::sqrt(2. * rc_ener); // reaction coordinate momentum

  // scale the amplitude of internal motion
  //
  dm *= nr_len / vlength(dm.begin(), Dynamic::Momenta::size());

  // add reaction coordinate motion
  //
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    //
    dm[i] += rc_len * rc_vec[i];

  // converting to regular velocities

  // orbital motion
  //
  for (int i = 0; i < 3; ++i)
    //
    dm.orb_vel()[i] /= Structure::mass_sqrt();

  // fragments rotational motion
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::type(frag) != Molecule::MONOATOMIC)
      //
      for (int i = 0; i < 3; ++i)
	//
	dm.ang_vel(frag)[i] /= Structure::fragment(frag).imom_sqrt(i);

  // statistical weight
  //
  return std::pow(ener, double(Structure::tm_dof() - 1) / 2.);
}

// generate random generalized velocities from the reaction coordinate vector for the given energy and angular momentum
//
double DivSur::Primitive::_rc2jv(const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, double ener, double amom_value, Dynamic::Momenta& dm)
{
  int        itemp;
  double     dtemp;
  D3::Vector vtemp;

  // inertia moments matrix
  //
  Lapack::Cholesky imm(dc.imm());

  // total angular momentum vector
  //
  Lapack::Vector ang_mom(3);
  ang_mom = 0.;
  ang_mom[2] = amom_value;

  // total angular velocity
  //
  Lapack::Vector ang_vel = imm.invert(ang_mom);
  
  // rotational kinetic energy shift
  //
  ener -= ang_mom * ang_vel / 2.;

  if(ener <= 0.)
    //
    return -1.;

  // random generalized velocity
  //
  Random::orient(dm, Dynamic::Momenta::size());

  // for linear fragments remove the component colinear to the molecular axis
  //
  for (int frag = 0; frag < 2; ++frag)
    //
    if(Structure::fragment(frag).type() == Molecule::LINEAR)
      //
      orthogonalize(dm.ang_vel(frag), dc.ang_pos(frag), 3);

  // angular momentum associated with the chosen generalized random velocity
  //
  Lapack::Vector ass_ang_mom(3);

  D3::Vector frag_ang_mom;
  
  // orbital part
  //
  D3::vprod(dc.orb_pos(), dm.orb_vel(), ass_ang_mom);
  
  ass_ang_mom *= Structure::mass_sqrt();
  
  // fragments rotational part
  //
  for(int frag = 0; frag < 2; ++frag) {
    //
    switch (Structure::type(frag)) {
      //
    case Molecule::LINEAR:
      //
      for (int i = 0; i < 3; ++i)
	//
	ass_ang_mom[i] += Structure::fragment(frag).imom_sqrt(0) * dm.ang_vel(frag)[i];
      
      break;
      //
    case Molecule::NONLINEAR:
      //
      for (int i = 0; i < 3; ++i)
	//
	vtemp[i] = Structure::fragment(frag).imom_sqrt(i) * dm.ang_vel(frag)[i];
      
      dc.mf2lf(frag, vtemp, frag_ang_mom);
      
      for (int i = 0; i < 3; ++i)
	//
	ass_ang_mom[i] += frag_ang_mom[i];
    }
  }

  // angular velocity associated with the chosen generalized random velocity
  //
  Lapack::Vector ass_ang_vel = imm.invert(ass_ang_mom);
  
  /************** REMOVE OVERALL ROTATION FROM THE GENERALIZED RANDOM VELOCITY **************/
  
  // orbital part
  //
  D3::vprod(ass_ang_vel, dc.orb_pos(), vtemp);
  
  for(int i = 0; i < 3; ++i)
    //
    dm.orb_vel()[i] -= Structure::mass_sqrt() * vtemp[i];

  // fragments rotational part
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    switch (Structure::type(frag)) {
      //
    case Molecule::MONOATOMIC:
      
      break;
      //
    case Molecule::LINEAR:
      //
      for (int i = 0; i < 3; ++i)
	//
	vtemp[i] = ass_ang_vel[i];
      
      orthogonalize(vtemp, dc.ang_pos(frag), 3);
      
      for (int i = 0; i < 3; ++i)
	//
	dm.ang_vel(frag)[i] -= Structure::fragment(frag).imom_sqrt(0) * vtemp[i];
      
      break;
      //
    case Molecule::NONLINEAR:
      //
      dc.lf2mf(frag, ass_ang_vel, vtemp);
      
      for (int i = 0; i < 3; ++i)
	//
	dm.ang_vel(frag)[i] -= Structure::fragment(frag).imom_sqrt(i) * vtemp[i];
      
      break;
    }

  // remove reactive component
  //
  orthogonalize(dm.begin(), rc_vec.begin(), Dynamic::Momenta::size());

  // energy of internal rotations
  //
  double in_ener = ener * std::pow(Random::flat(), 2. / double(Structure::tm_dof() - 4));
  double rc_ener = ener - in_ener;         // reaction coordinate energy
  double rc_len = std::sqrt(2. * rc_ener); // reaction coordinate momentum
  double in_len = std::sqrt(2. * in_ener); // momentum amplitude of internal rotations

  // scale the amplitude of internal motion
  //
  dm *= in_len / vlength(dm.begin(), Dynamic::Momenta::size());

  // add reaction coordinate motion
  //
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    //
    dm[i] += rc_len * rc_vec[i];

  /************************ ADD OVERALL ROTATION ***********************/
  
  // orbital part
  //
  D3::vprod(ang_vel, dc.orb_pos(), vtemp);
  
  for(int i = 0; i < 3; ++i)
    //
    dm.orb_vel()[i] += Structure::mass_sqrt() * vtemp[i];

  // fragments rotational part
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    switch (Structure::type(frag)) {
      //
    case Molecule::MONOATOMIC:
      
      break;
      //
    case Molecule::LINEAR:
      //
      for (int i = 0; i < 3; ++i)
	//
	vtemp[i] = ang_vel[i];
      
      vtemp.orthogonalize(dc.ang_pos(frag));
      
      for (int i = 0; i < 3; ++i)
	//
	dm.ang_vel(frag)[i] += Structure::fragment(frag).imom_sqrt(0) * vtemp[i];
      
      break;
      //
    case Molecule::NONLINEAR:
      //
      dc.lf2mf(frag, ang_vel, vtemp);
      
      for (int i = 0; i < 3; ++i)
	//
	dm.ang_vel(frag)[i] += Structure::fragment(frag).imom_sqrt(i) * vtemp[i];
    }
  
  /******************** CONVERT TO REGULAR VELOCITIES ********************/
  
  // orbital motion
  //
  for (int i = 0; i < 3; ++i)
    //
    dm.orb_vel()[i] /= Structure::mass_sqrt();

  // fragments rotational motion
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::type(frag) != Molecule::MONOATOMIC)
      //
      for(int i = 0; i < 3; ++i)
	//
	dm.ang_vel(frag)[i] /= Structure::fragment(frag).imom_sqrt(i);
  
  return std::pow(ener, double(Structure::tm_dof() - 4) / 2.) / imm.det_sqrt();
}

/*****************************************************************************
 **************************** SPHERICAL FACET ********************************
 *****************************************************************************/

void DivSur::Sphere::print (std::ostream& to) const
{
  to << "sphere = {";
  
  for(int frag = 0; frag < 2; ++frag)
    //
    to << frag << "-th pivot = " << _pivot[frag] << ", ";

  to << "radius = " << _dist << "}";
}

void DivSur::Sphere::random_orient (Dynamic::Coordinates& dc) const
{
  double dtemp;
  int    itemp;
  double vtemp [3];

  // randomly orient molecules
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::fragment(frag).type() != Molecule::MONOATOMIC)
      //
      Random::orient(dc.write_ang_pos(frag), Structure::fragment(frag).pos_size());

  // randomly orient a vector between the pivot points
  //
  Random::orient(dc.orb_pos(), 3);
  
  for (int i = 0; i < 3; ++i)
    //
    dc.orb_pos()[i] *= _dist;

  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::type(frag) != Molecule::MONOATOMIC) {
      //
      dc.mf2lf(frag, _pivot[frag], vtemp);
      
      if (!frag) {
	//
	for(int i = 0; i < 3; ++i)
	  //
	  dc.orb_pos()[i] += vtemp[i];
      }
      else
	//
	for (int i = 0; i < 3; ++i)
	  //
	  dc.orb_pos()[i] -= vtemp[i];
    }
}

// pp2pp vector, 1->2
//
D3::Vector DivSur::Sphere::_lf_pp_12 (const Dynamic::Coordinates& dc) const
{
  double     dtemp;
  int        itemp;
  D3::Vector vtemp;

  D3::Vector res(dc.orb_pos());

  for(int frag = 0; frag < 2; ++frag)
    //
    if(Structure::fragment(frag).type() != Molecule::MONOATOMIC) {
      //
      dc.mf2lf(frag, _pivot[frag], vtemp);
      
      if (!frag) {
	//
	res -= vtemp;
      }
      else
	//
	res += vtemp;
    }

  return res;
}

// reaction coordinate generalized vector
//
double DivSur::Sphere::_set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_vec) const
{
  double     dtemp;
  int        itemp;
  D3::Vector vtemp;

  // pp2pp vector, 1->2
  //
  D3::Vector lf_pp_12 = _lf_pp_12(dc);

  // orbital part of the reaction coordinate vector
  //
  for (int i = 0; i < 3; ++i)
    //
    rc_vec.orb_vel()[i] = -lf_pp_12[i] / Structure::mass_sqrt();

  // angular part of the reaction coordinate vector
  //
  for(int frag = 0; frag < 2; ++frag)
    //
    switch (Structure::fragment(frag).type()) {
      //
    case Molecule::LINEAR:
      //
      if(frag) {
	//
	D3::vprod(lf_pp_12, dc.ang_pos(frag), rc_vec.ang_vel(frag));
      }
      else
	//
	D3::vprod(dc.ang_pos(frag), lf_pp_12, rc_vec.ang_vel(frag));

      dtemp = _pivot[frag][0] / Structure::fragment(frag).imom_sqrt(0);
      
      for (int i = 0; i < 3; ++i)
	//
	rc_vec.ang_vel(frag)[i] *= dtemp;
      
      break;
      //
    case Molecule::NONLINEAR:
      //
      dc.lf2mf(frag, lf_pp_12, vtemp);
      
      if(frag) {
	//
	D3::vprod(vtemp, _pivot[frag], rc_vec.ang_vel(frag));
      }
      else
	//
	D3::vprod(_pivot[frag], vtemp, rc_vec.ang_vel(frag));

      for (int i = 0; i < 3; ++i)
	//
	rc_vec.ang_vel(frag)[i] /= Structure::fragment(frag).imom_sqrt(i);
    }

  // normalization factor
  //
  return  _dist * normalize(rc_vec.begin(), Dynamic::Momenta::size());
}

/********************************************************************************
 ****************************** CONICAL FACET METHODS ***************************
 ********************************************************************************/

void DivSur::Cone::print (std::ostream& to) const
{
  to << "cone   = {";
    
  for(int frag = 0; frag < 2; ++frag) {
    //
    if(frag)
      //
      to << ", ";
    
    to << frag << "-th fragment";
    
    if(frag == _cone_frag) {
      //
      to << " head = "   << _head
	 << " axis = "   << _axis
	 << " proj = "   << _proj
	 << " radius = " << _radius;
    }
    else
      //
      to << " pivot = " << _pivot;
  }
  
  to << "}";
}

void DivSur::Cone::random_orient (Dynamic::Coordinates& dc) const
{
  const char funame [] = "DivSur::Cone::random_orient: ";

  double dtemp;
  int    itemp;
  double vtemp [3];

  const int  point_frag = 1 - _cone_frag;

  // randomly orient molecules
  //
  for (int frag = 0; frag < 2; ++frag)
    //
    if (Structure::type(frag) != Molecule::MONOATOMIC)
      //
      Random::orient(dc.write_ang_pos(frag), Structure::fragment(frag).pos_size());

  // random point on the cone surface within a certain distance from its head
  //
  double v[2], z;

  dtemp = Random::flat();

  v[0] = std::cos(2. * M_PI * dtemp);
  
  v[1] = std::sin(2. * M_PI * dtemp);

  dtemp = std::pow(Random::flat(), 1. / 3.) * _radius;

  double con_sin = std::sqrt(1. - _proj * _proj);
  
  for(int i = 0; i < 2; ++i)
    //
    v[i] *= dtemp * con_sin;

  z = dtemp * _proj * _radius;
  
  
  D3::Vector n[2];

  switch(Structure::type(_cone_frag)) {
    //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      dc.orb_pos()[i] = dc.ang_pos(_cone_frag)[i];
    
    // vectors orthogonal to the molecular axis in the lab frame
    //
    find_orth(dc.orb_pos(), n);

    // random point on the cone in the lab frame
    //
    dtemp = _head[0] + z * _axis[0];
    
    for(int i = 0; i < 3; ++i) {
      //
      dc.orb_pos()[i] *= dtemp;

      for(int j = 0; j < 2; ++j)
	//
	dc.orb_pos()[i] +=  v[j] * n[j][i]; 
    }
    
    break;
    //
  case Molecule::NONLINEAR:
    //
    // vectors orthogonal to the cone axis in the molecular frame
    //
    find_orth(_axis, n);

    // random point on the cone in the molecular frame
    //
    for(int i = 0; i < 3; ++i) {
      //
      vtemp[i] = _head[i] + z * _axis[i];

      for(int j = 0; j < 2; ++j)
	//
	vtemp[i] += v[j] * n[j][i];
    }
    
    // convert to the lab frame
    //
    dc.mf2lf(_cone_frag, vtemp, dc.orb_pos());
  }

  // shift by the pivot point vector
  //
  if(Structure::type(point_frag) != Molecule::MONOATOMIC) {
    //
    dc.mf2lf(point_frag, _pivot, vtemp);
    
    for(int i = 0; i < 3; ++i)
      //
      dc.orb_pos()[i] -= vtemp[i];
  }

  if(_cone_frag)
    //
    for (int i = 0; i < 3; ++i)
      //
      dc.orb_pos()[i] = -dc.orb_pos()[i];
}

double DivSur::Cone::distance (const Dynamic::Coordinates& dc) const
{
  const char funame [] = "DivSur::Cone::distance: ";

  D3::Vector vtemp;

  const int point_frag = 1 - _cone_frag;

  // plane in the lab frame
  //
  D3::Vector lf_head, lf_axis;

  dc.mf2lf(_cone_frag, _head, lf_head);

  dc.mf2lf(_cone_frag, _axis, lf_axis);

  // pivot point position relative to the cone head in the labooratory frame
  //
  D3::Vector ph_vec;
  
  dc.mf2lf(point_frag, _pivot, ph_vec);

  if(_cone_frag) {
    //
    ph_vec -= dc.orb_pos();
  }
  else
    //
    ph_vec += dc.orb_pos();

  ph_vec -= lf_head;
  
  ph_vec.normalize();

  return vdot(ph_vec, lf_axis) - _proj;
}

double DivSur::Cone::_set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_vec) const
{
  const char funame [] = "DivSur::Cone::_set_rc: ";

  double     dtemp;
  int        itemp;
  D3::Vector vtemp;

  const int point_frag = 1 - _cone_frag;

  // cone head and axis in the laboratory frame
  //
  D3::Vector lf_head, lf_axis;

  dc.mf2lf(_cone_frag, _head, lf_head);

  dc.mf2lf(_cone_frag, _axis, lf_axis);

  // laboratory frame pivot point
  //
  D3::Vector lf_pivot;

  dc.mf2lf(point_frag, _pivot, lf_pivot);

  // p-h vector between pivot point and cone head
  //
  D3::Vector n = lf_pivot;

  if(!_cone_frag) {
    //
    n += dc.orb_pos();
  }
  else
    //
    n -= dc.orb_pos();

  n -= lf_head;

  double ph_len = n.normalize();

  // cosine between the cone axis and p-h vector (reaction coordinate)
  //
  double rc_val = ::vdot(n, lf_axis, 3);

  // b = (a - n * (n*a)) / ph_len
  //
  D3::Vector b = lf_axis;

  b -= rc_val * n;

  b /= ph_len;

  // orbital part of the reaction coordinate vector
  //
  if(!_cone_frag) {
    //
    for(int i = 0; i < 3; ++i)
      //
      rc_vec.orb_vel()[i] = b[i] / Structure::mass_sqrt();
  }
  else
    //
    for(int i = 0; i < 3; ++i)
      //
      rc_vec.orb_vel()[i] = -b[i] / Structure::mass_sqrt();

  // angular part of the reaction coordinate vector for the pivot point fragment
  //
  D3::vprod(lf_pivot, b, vtemp);

  switch(Structure::type(point_frag)) {
    //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      rc_vec.ang_vel(point_frag)[i] = vtemp[i];

    break;
    //
  case Molecule::NONLINEAR:
    //
    dc.lf2mf(point_frag, vtemp, rc_vec.ang_vel(point_frag));
  }

  // angular part of the reaction coordinate vector for the cone fragment
  //
  D3::Vector cone_rc;

  D3::vprod(b, lf_head, cone_rc);

  D3::vprod(lf_axis, n, vtemp);

  cone_rc += vtemp;
  
  switch(Structure::type(_cone_frag)) {
    //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      rc_vec.ang_vel(_cone_frag)[i] = cone_rc[i];

    break;
    //
  case Molecule::NONLINEAR:
    //
    dc.lf2mf(_cone_frag, cone_rc, rc_vec.ang_vel(_cone_frag));
  }

  // normalize angular part of the reaction coordinate vector
  //
  for(int f = 0; f < 2; ++f)
    //
    if(Structure::type(f) != Molecule::MONOATOMIC)
      //
      for (int i = 0; i < 3; ++i)
	//
	rc_vec.ang_vel(f)[i] /= Structure::fragment(f).imom_sqrt(i);
      
  // normalize reaction coordinate vector
  //
  return normalize(rc_vec.begin(), Dynamic::Momenta::size()) * _radius * _radius * _radius / 6.;
}

/******************************************************************************
 ****************************** PLANE FACET METHODS ***************************
 ******************************************************************************/

void DivSur::Plane::print (std::ostream& to) const
{
  to << "plane  = {";
    
  for(int frag = 0; frag < 2; ++frag) {
    //
    if(frag)
      //
      to << ", ";
    
    to << frag << "-th fragment ";
    
    if(frag == _plane_frag) {
      //
      to << "plane = " << _plane;
    }
    else
      //
      to << "pivot = " << _pivot;
  }
  
  to << ", sampling radius = " << _radius;

  to << "}";
}

void DivSur::Plane::random_orient (Dynamic::Coordinates& dc) const
{
  const char funame [] = "DivSur::Plane::random_orient: ";

  if(_radius < 0.) {
    //
    std::cerr << funame << "sampling circle was not set up\n";
    
    throw Error::Init();
  }

  double dtemp;
  int    itemp;
  double vtemp [3];

  const int  point_frag = 1 - _plane_frag;
  
  // randomly orient molecules
  //
  for (int frag = 0; frag < 2; ++frag)
    //
    if (Structure::type(frag) != Molecule::MONOATOMIC)
      //
      Random::orient(dc.write_ang_pos(frag), Structure::fragment(frag).pos_size());

  // random point on the circle
  //
  double v[2];

  dtemp = Random::flat();

  v[0] = std::cos(2. * M_PI * dtemp);
  
  v[1] = std::sin(2. * M_PI * dtemp);

  dtemp = std::sqrt(Random::flat()) * _radius;

  for(int i = 0; i < 2; ++i)
    //
    v[i] *= dtemp;
  
  D3::Vector n[2];

  switch(Structure::type(_plane_frag)) {
    //
  case Molecule::LINEAR:
    //
    // plane normal in the lab frame
    //
    for(int i = 0; i < 3; ++i)
      //
      dc.orb_pos()[i] = dc.ang_pos(_plane_frag)[i];
    
    // vectors orthogonal to the plane normal in the lab frame
    //
    find_orth(dc.orb_pos(), n);
    
    // random point on the plane within the given circle in the lab frame
    //
    for (int i = 0; i < 3; ++i) {
      //
      dc.orb_pos()[i] *= _plane.dist();
      
      for(int j = 0; j < 2; ++j)
	//
	dc.orb_pos()[i] += v[j] * n[j][i]; 
    }
    
    break;
    //
  case Molecule::NONLINEAR:
    //
    // random point on the plane within the given circle in the molecular frame
    //
    for(int i = 0; i < 3; ++i) {
      //
      vtemp[i] = _plane.dist() * _plane.normal()[i];
      
      for(int j = 0; j < 2; ++j)
	//
	vtemp[i] += v[j] * _plane.orth(j)[i];
    }
    
    // transform to the lab frame
    //
    dc.mf2lf(_plane_frag, vtemp, dc.orb_pos());
  }

  // shift by the pivot point vector
  //
  if(Structure::type(point_frag) != Molecule::MONOATOMIC) {
    //
    dc.mf2lf(point_frag, _pivot, vtemp);
    
    for(int i = 0; i < 3; ++i)
      //
      dc.orb_pos()[i] -= vtemp[i];
  }

  if(_plane_frag)
    //
    for (int i = 0; i < 3; ++i)
      //
      dc.orb_pos()[i] = -dc.orb_pos()[i];
}

double DivSur::Plane::distance (const Dynamic::Coordinates& dc) const
{
  const char funame [] = "DivSur::Plane::distance: ";

  D3::Vector vtemp;

  const int  point_frag = 1 - _plane_frag;
  
  // plane normal in the lab frame
  //
  D3::Vector lf_normal;

  dc.mf2lf(_plane_frag, _plane.normal(), lf_normal);

  // pivot point relative to the plane fragment origin
  //
  D3::Vector pp_vec;

  dc.mf2lf(point_frag, _pivot, pp_vec);
    
  if(_plane_frag) {
    //
    pp_vec -= dc.orb_pos();
  }
  else
    //
    pp_vec += dc.orb_pos();
  
  return vdot(pp_vec, lf_normal) - _plane.dist();
}

double DivSur::Plane::_set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_vec) const
{
  const char funame [] = "DivSur::Plane::_set_rc: ";

  double     dtemp;
  int        itemp;
  D3::Vector vtemp;

  const int point_frag = 1 - _plane_frag;
  
  // plane normal in the lab frame
  //
  D3::Vector lf_normal;

  dc.mf2lf(_plane_frag, _plane.normal(), lf_normal);

  // pivot point in the lab frame
  //
  D3::Vector lf_pivot;
  
  dc.mf2lf(point_frag, _pivot, lf_pivot);

  // pivot point relative to the plane fragment origin
  //
  D3::Vector pp_vec = lf_pivot;
  
  if(!_plane_frag) {
    //
    pp_vec += dc.orb_pos();
  }
  else
    //
    pp_vec -= dc.orb_pos();

  /*************************************************************
   ********** REACTION COORDINATE GENERALIZED VECTOR ***********
   *************************************************************/

  // orbital part of the reaction coordinate vector
  //
  if(!_plane_frag) {
    //
    for (int i = 0; i < 3; ++i)
      //
      rc_vec.orb_vel()[i] = lf_normal[i] / Structure::mass_sqrt();
  }
  else
    //
    for (int i = 0; i < 3; ++i)
      //
      rc_vec.orb_vel()[i] = -lf_normal[i] / Structure::mass_sqrt();
  
  // angular part of the reaction coordinate vector for the pivot fragment
  //
  D3::vprod(lf_pivot, lf_normal, vtemp);
	
  switch (Structure::type(point_frag)) {
      //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      rc_vec.ang_vel(point_frag)[i] = vtemp[i];
	
    break;
    //
  case Molecule::NONLINEAR:
    //
    dc.lf2mf(point_frag, vtemp, rc_vec.ang_vel(point_frag));
  }

  // angular part of the reaction coordinate vector for the plane fragment
  //
  D3::vprod(lf_normal, pp_vec, vtemp);
      
  switch (Structure::type(_plane_frag)) {
    //
  case Molecule::LINEAR:
    //
    for(int i = 0; i < 3; ++i)
      //
      rc_vec.ang_vel(_plane_frag)[i] = vtemp[i];
      
    break;
    //
  case Molecule::NONLINEAR:
    //
    dc.lf2mf(_plane_frag, vtemp, rc_vec.ang_vel(_plane_frag));
  }

  // normalize angular part of the reaction coordinate vector
  //
  for(int f = 0; f < 2; ++f)
    //
    if(Structure::type(f) != Molecule::MONOATOMIC)
      //
      for (int i = 0; i < 3; ++i)
	//
	rc_vec.ang_vel(f)[i] /= Structure::fragment(f).imom_sqrt(i);
      
  // normalize reaction coordinate vector
  //
  return normalize(rc_vec.begin(), Dynamic::Momenta::size()) * _radius * _radius / 4.;
}

/*****************************************************************************
 ****************** MULTIPLE SPECIES DIVIDING SURFACE ************************
 *****************************************************************************/

void DivSur::MultiSur::read (std::istream& from) 
{
  const char funame [] = "DivSur::MultiSur::read: ";

  double dtemp;
  int    itemp;

  if(!Structure::isinit()) {
    //
    std::cerr << funame << "molecular structure has not been initialized\n";
    
    throw Error::Init();
  }

  if(_isinit) {
    //
    std::cerr << funame << "already initialized\n";
    
    throw Error::Init();
  }

  _isinit = true;

  std::string spec_key = "Species";
  std::string prim_key = "Primitives";

  std::map<std::string, int> fmap;

  std::string token, comment, name;
  
  /***************************************************************
   *************************** FRAGMENTS *************************
   ***************************************************************/
  
  FragData fragment[2];
  
  for(int f = 0; f < 2; ++f)
    //
    fragment[f].read(f, from);

  /**************************************************************
   *********************** PRIMITIVES ***************************
   **************************************************************/

  if(!(from >> token)) {
    //
    std::cerr << funame << "corrupted: should be " << prim_key << "\n";

    throw Error::Input();
  }

  if(token != prim_key) {
    //
    std::cerr << funame << "unknown key: " << token << ": should be " << prim_key << "\n";

    throw Error::Input();
  }

  std::getline(from, comment);

  KeyGroup PrimGroup;
  Key  plane_key ("Plane");
  Key sphere_key ("Sphere");
  Key   cone_key ("Cone");
  
  // read cycle
  //
  while(from >> token) {
    //
    if(token == IO::end_key())
      //
      break;

    IO::LineInput lin(from);
    
    // read primitive surface name
    //
    if(!(lin >> name)) {
      //
      std::cerr << funame << "cannot read primitive's name\n";

      throw Error::Input();
    }
    
    if(fmap.find(name) != fmap.end()) {
      //
      std::cerr << funame << name << " primitive's name already has been used\n";
      
      throw Error::Form();
    }
    
    fmap[name] = _prim.size();

    ConstSharedPointer<GeomObject> gop[2];
    
    for(int f = 0; f < 2; ++f) {
      //
      if(!(lin >> name)) {
	//
	std::cerr << funame << _prim.size() << "-th primitive: "<< f <<"-th fragment: cannot read geom. object name\n";

	throw Error::Input();
      }

      if(fragment[f].find(name) == fragment[f].end()) {
	//
	std::cerr << funame << f << "-th fragment: cannot find " << name << " geom. object\n";
	
	throw Error::Init();
      }
      
      gop[f] = fragment[f][name];
    }
    
    // plane
    //
    if(token == plane_key) {
      //
      _prim.push_back(SharedPointer<Primitive>(new Plane(gop)));
    }
    // cone
    //
    else if(token == cone_key) {
      //
      _prim.push_back(SharedPointer<Primitive>(new Cone(gop)));
    }
    // sphere
    //
    else if(token == sphere_key) {
      //
      // read radius
      //
      if(!(lin >> dtemp)) {
	//
	std::cerr << funame << _prim.size() << "-th primitive: cannot read radius\n";

	throw Error::Input();
      }

      if(dtemp <= 0.) {
	//
	std::cerr << funame << _prim.size() << "-th primitive: radius out of range: " << dtemp << "\n";

	throw Error::Range();
      }
      
      _prim.push_back(SharedPointer<Primitive>(new Sphere(gop, dtemp)));
    }
    // unknown primitive
    //
    else {
      //
      std::cerr << funame << "unknown primitive: " << token << "\n";
      
      Key::show_all(std::cerr);
      std::cerr << "\n";
      
      throw Error::Form();
    }
    //
  }// read cycle

  if(!from) {
    //
    std::cerr << funame << "input stream is corrupted\n";
    
    throw Error::Form();
  }

  // sampling radius for the planes setup
  //
  double rad = -1.;
  
  for(std::vector<SharedPointer<Primitive> >::iterator sit = _prim.begin(); sit != _prim.end(); ++sit) {
    //
    const Sphere* p = dynamic_cast<Sphere*>((Primitive*)(*sit));
    
    if(p) {
      //
      dtemp = p->max_cm_dist();
      
      if(dtemp > rad)
	//
	rad = dtemp;
    }
  }

  if(rad < 0.) {
    //
    std::cerr << funame << "no spherical surfaces found\n";
    
    throw Error::Init();
  }

  for(std::vector<SharedPointer<Primitive> >::iterator sit = _prim.begin(); sit != _prim.end(); ++sit) {
    //
    Plane* p = dynamic_cast<Plane*>((Primitive*)(*sit));
    
    if(p)
      //
      p->set_circle(rad);
  }
  
  /**********************************************************************
   ************************ SPECIES DEFINITIONS *************************
   **********************************************************************/

  if(!(from >> token)) {
    //
    std::cerr << funame << "corrupted: should be " << spec_key << "\n";

    throw Error::Input();
  }

  if(token != spec_key) {
    //
    std::cerr << funame << "unknown key: " << token << ": should be " << spec_key << "\n";

    throw Error::Input();
  }
  
  if(!_prim.size()) {
    //
    std::cerr << funame << token << ": primitives should be defined first\n";

    throw Error::Init();
  }

  if(_species.size()) {
    //
    std::cerr << funame << token << ": already defined\n";

    throw Error::Init();
  }
      
  if(!(from >> itemp)) {
    //
    std::cerr << funame << token << ": corrupted\n";

    throw Error::Input();
  }
  std::getline(from, comment);

  if(itemp < 1) {
    //
    std::cerr << funame << token << ": out of range\n";

    throw Error::Range();
  }
  _species.resize(itemp);

  for(int s = 0; s < _species.size(); ++s) {
    //
    IO::LineInput lin(from);

    _species[s] = Logical::read_expr(lin);
    _species[s]->init(fmap);
  }

  // end key
  //
  if(!(from >> token)) {
    //
    std::cerr << funame << "corrupted: should be " << IO::end_key() << "\n";

    throw Error::Input();
  }

  if(token != IO::end_key()) {
    //
    std::cerr << funame << "unknown key: " << token << ": should be " << IO::end_key() << "\n";

    throw Error::Input();
  }

  std::getline(from, comment);
    
  // checking
  //
  if(!from) {
    //
    std::cerr << funame << "stream is corrupted\n";
    
    throw Error::Input();
  }
    
  if(!_species.size()) {
    //
    std::cerr << funame << "species not defined\n";
    
    throw Error::Init();
  }

  // Checking species for logical intersections;
  //
  std::vector<bool> pset(_prim.size(), false);
  std::vector<bool> spec_res(_species.size());
  int pi, spec_num;

  while(1) {
    //
    spec_num = 0;
    
    for(int si = 0; si < _species.size(); ++si)
      //
      if(spec_res[si] = _species[si]->evaluate(pset))
	//
	spec_num++;
    
    if(spec_num > 1) {
      //
      std::cerr << funame << ": WARNING: the primitives logical configuration (";
      
      for(pi = 0; pi < _prim.size(); ++pi)
	std::cerr << pset[pi];
      std::cerr << ") is compatible with multiple species (";
      for(int si = 0; si < _species.size(); ++si)
	std::cerr << spec_res[si];
      std::cerr << ")\n";      
    }

    for(pi = 0;  pi < _prim.size(); ++pi) {
      pset[pi] = !pset[pi];
      if(pset[pi])
	break;
    }
    
    if(pi == _prim.size())
      break;
  }
}

void DivSur::MultiSur::print (std::ostream& to, const std::string& prefix) const
{
  to << prefix << "primitives:\n";

  for(int p = 0; p < _prim.size(); ++p) {
    //
    to << prefix;
    _prim[p]->print(to);
    
    to << "\n";
  }
}

int DivSur::MultiSur::classify  (const Dynamic::Coordinates& dc) const 
{
  const char funame [] = "DivSur::MultiSur::classify: ";

  std::vector<bool> prim_test(primitive_size());
  
  for(int i = 0; i < _prim.size(); ++i)
    //
    prim_test[i] = _prim[i]->test(dc);

  std::vector<int> res;
  
  for(int s = 0; s < _species.size(); ++s)
    //
    if(_species[s]->evaluate(prim_test))
      //
      res.push_back(s);

  // bimolecular species
  //
  if(!res.size())
    //
    return _species.size();

  if(res.size() == 1)
    //
    return *res.begin();

  // logical error
  //
  std::ostringstream oss;
  
  oss << funame << "configuration belongs to multiple species:";
  
  for(int i = 0; i < res.size(); ++i)
    //
    oss << " " << res[i];
  
  oss << "\n";

  IO::log << oss.str() << std::flush;

  std::cerr << oss.str();
  
  throw Error::Logic();
}

DivSur::SmpRes DivSur::MultiSur::surface_test (int sur, const Dynamic::Coordinates& dc) const 
{
  const char funame [] = "DivSur::MultiSur::surface_test: ";

  if(!Structure::isinit()) {
    //
    std::cerr << funame << "structure has not yet been initialized\n";
    
    throw Error::Init();
  }

  SmpRes res;
    
  // configuration belongs to the excluded region
  //
  if(Dynamic::exclude_region && Dynamic::exclude_region->test(dc)) {
    //
    res.stat = EXCLUDE;
    
    return res;
  }

  // atoms are too close
  //
  if(dc.are_atoms_close()) {
    //
    res.stat = CLOSE;
    
    return res;
  }

  // to which actual facet does the configuration belong
  //
  std::vector<bool> prim_test(primitive_size());
  
  for(int i = 0; i < primitive_size(); ++i)
    //
    if(i != sur)
      //
      prim_test[i] = _prim[i]->test(dc);

  // species tests
  //
  std::vector<bool> inner_test(_species.size()), outer_test(_species.size());
  for(int s = 0; s < _species.size(); ++s) {
    //
    prim_test[sur] = true;
    inner_test[s] = _species[s]->evaluate(prim_test);

    prim_test[sur] = false;
    outer_test[s] = _species[s]->evaluate(prim_test);
  }

  std::set<int> inner, from, to;
  for(int s = 0; s < _species.size(); ++s)
    //
    if(inner_test[s] && outer_test[s]) {
      //
      inner.insert(s);
    }
    else if(!inner_test[s] && outer_test[s]) {
      //
      from.insert(s);
    }
    else if(inner_test[s] && !outer_test[s])
      //
      to.insert(s);

  // analyze the results

  // inner region
  //
  if(inner.size() <= 1 && !to.size() && !from.size()) {
    //
    res.stat = INNER;
    
    return res;
  }

  // from bimolecular products to a bound species
  //
  if(!inner.size()  && to.size() == 1 && !from.size()) {
    //
    res.face.first  = _species.size();
    res.face.second =     *to.begin();

    res.stat = FACET;

    return res;
  }
  
  // from a bound species to bimolecular products
  //
  if(!inner.size()  && from.size() == 1 && !to.size()) {
    //
    res.face.first   =   *from.begin();
    res.face.second  = _species.size();

    res.stat = FACET;

    return res;
  }
  
  // from a bound species to another bound species
  //
  if(!inner.size()  && from.size() == 1 && to.size() == 1) {
    //
    res.face.first   = *from.begin();
    res.face.second  =   *to.begin();

    res.stat = FACET;

    return res;
  }
  
  // logic failure
  //
  std::ostringstream oss;
  
  oss << funame << sur << "-th surface: sampling error: configuration:";

  if(inner.size()) {
    //
    oss << " inner=(";
    
    for(std::set<int>::iterator it = inner.begin(); it != inner.end(); ++it) {
      //
      if(it != inner.begin())
	//
	oss << ", ";
      
      oss << *it;
    }
    
    oss << ")";
  }
  
  if(from.size()) {
    //
    oss << " from=(";
    
    for(std::set<int>::iterator it = from.begin(); it != from.end(); ++it) {
      //
      if(it != from.begin())
	//
	oss << ", ";
      
      oss << *it;
    }
    oss << ")";
  }

  if(to.size()) {
    //
    oss << " to=(";
    
    for(std::set<int>::iterator it = to.begin(); it != to.end(); ++it) {
      //
      if(it != to.begin())
	//
	oss << ", ";
      
      oss << *it;
    }
    
    oss << ")";
  }

  oss << "\n";

  std::cerr << oss.str();

  IO::log << oss.str();
  
  res.stat = FAIL;
  
  return res;
}

