#include "divsur.hh"
#include "random.hh"
#include "units.hh"
#include "key.hh"

#include <cmath>
#include <sstream>

/********************************************************************************
 *                                    Input/Output                              *
 ********************************************************************************/

void DivSur::FragData::read(int frag, std::istream& from) throw(Error::General)
{
  const char funame [] = "DivSur::FragData::read: ";

  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not yet been initialized\n";
    throw Error::Init();
  }

  KeyGroup FragmentData;
  Key plane_key ("Plane");
  Key pivot_key ("Pivot");

  std::string token, comment, line, obj_name;
  D3::Vector pos;
  double dtemp;

  data.clear();

  // fragment name
  from >> name;
  std::getline(from, comment);

  // add center of mass as a pivot point
  pos = 0.;
  data[name] = ConstSharedPointer<GeomObject>(new GeomPivot(pos));

  while(from >> token) {// read cycle
    if(token == IO::end_key()) {      
      std::getline(from, comment);
      return;
    }

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "no explicit geometrical objects for monoatomic fragments\n";
      throw Error::Form();
    }

    from >> obj_name;
    std::getline(from, line);
    if(!from) {
      std::cerr << funame << frag << "-th fragment: flux is corrupted\n";
      throw Error::Form();
    }
    if(data.find(obj_name) != data.end()) {
      std::cerr << funame << frag << "-th fragment: the " << obj_name 
		<< " geometrical object already exists\n";
      throw Error::Form();
    }

    std::istringstream iss(line);
    std::vector<double> val;
    while(iss >> dtemp)
      val.push_back(dtemp);

    if(token == plane_key) {// plane
      switch(Structure::fragment(frag).type()) {
      case Molecule::LINEAR:
	pos = 0.; 
	pos[0] = 1.;
	if(val.size() != 1) {
	  std::cerr << funame << frag << "-th fragment: wrong number of arguments (" 
		    << val.size() <<") for the " << obj_name << " plane\n";
	  throw Error::Form();
	}
	data[obj_name] = ConstSharedPointer<GeomObject>(new GeomPlane(pos, val[0]));
	break;

      case Molecule::NONLINEAR:
	if(val.size() != 4) {
	  std::cerr << funame << frag << "-th fragment: wrong number of arguments ("
		    << val.size() << ") for the " << obj_name << " plane\n";
	  throw Error::Form();
	}
	for(int i = 0; i < 3; ++i)
	  pos[i] = val[i];
	data[obj_name] = ConstSharedPointer<GeomObject>(new GeomPlane(pos, val[3]));
	break;
      }

    }// plane

    else if(token == pivot_key) {// pivot point
      switch(Structure::fragment(frag).type()) {
      case Molecule::LINEAR:
	if(val.size() != 1 && val.size() != 3) {
	  std::cerr << funame << frag << "-th fragment: wrong number of arguments ("
		    << val.size() << ") for the " << obj_name << " pivot point\n";
	  throw Error::Form();
	}
	pos = 0.; 
	pos[0] = val[0];
	break;

      case Molecule::NONLINEAR:
	if(val.size() != 3) {
	  std::cerr << funame << frag << "-th fragment: wrong number of arguments ("
		    << val.size() << ") for the " << obj_name << " pivot point\n";
	  throw Error::Form();
	}
	for(int i = 0; i < 3; ++i)
	  pos[i] = val[i];
	break;
      }
      data[obj_name] = ConstSharedPointer<GeomObject>(new GeomPivot(pos));
    }// pivot point

    else {// unknown keyword
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }// unknown keyword
  }// read cycle

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Form();
  }
}

void DivSur::PrimSet::read (std::istream& from) throw(Error::General)
{
  const char funame [] = "DivSur::PrimSet::read: ";

  if(!Structure::isinit()) {
    std::cerr << funame << "molecular structure has not been initialized\n";
    throw Error::Init();
  }

  double dtemp;

  // read fragments data
  FragData fragment [2];
  for(int frag = 0; frag < 2; ++frag)
    fragment[frag].read(frag, from);

  KeyGroup PrimSet;
  Key   prim_key ("Primitives");
  Key  plane_key ("Plane");
  Key sphere_key ("Sphere");
  
  // read primitive surfaces
  std::string token;
  from >> token;
  if(token != prim_key) {
    std::cerr << funame << "unknown keyword: " << token << "; next should go "
	      << prim_key << " key\n";
    throw Error::Form();
  }

  while(from >> token) {// read cycle
    if(token == IO::end_key())
      break;

    // read primitive surface name
    std::string prim_name;
    from >> prim_name; 
    if(_fmap.find(prim_name) != _fmap.end()) {
      std::cerr << funame << prim_name << " primitive surface already exists\n";
      throw Error::Form();
    }

    // ... geometrical objects names
    std::string obj_name [2];
    for(int frag = 0; frag < 2; ++frag)
      from >> obj_name[frag];
    
    // find the objects in the fragments data
    std::map<std::string, ConstSharedPointer<GeomObject> >::const_iterator vit;
    ConstSharedPointer<GeomObject> gop [2];
    for(int frag = 0; frag < 2; ++frag) {
      vit = fragment[frag].data.find(obj_name[frag]);
      if(vit == fragment[frag].data.end()) {
	std::cerr << funame << prim_name << " primitive surface: could not find the " 
		  << obj_name[frag] << " geometric object in the " 
		  << frag << "-th fragment (" << fragment[frag].name 
		  << ") data section\n";
	throw Error::Form();
      }
      gop[frag] = vit->second;
    }
    
    if(token == plane_key) {// plane
      _fmap[prim_name] = _sur.size();
      _sur.push_back(SharedPointer<Primitive>(new Plane(gop)));
    }// plane
    else if(token == sphere_key) {// sphere
      from >> dtemp; // read radius
      _fmap[prim_name] = _sur.size();
      _sur.push_back(SharedPointer<Primitive>(new Sphere(gop, dtemp)));
    }// sphere
    else {// unknown primitive
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }// unknown primitive
  }// read cycle

  if(!from) {
    std::cerr << funame << "input stream is corrupted\n";
    throw Error::Form();
  }

  // setup the sampling radius for the planes
  std::vector<SharedPointer<Primitive> >::const_iterator sit;
  const Sphere* sph;
  double rad = -1.;
  for(sit = _sur.begin(); sit != _sur.end(); ++sit) {
    sph = dynamic_cast<Sphere*>((Primitive*)(*sit));
    if(sph) {
      dtemp = sph->max_cm_dist();
      if(dtemp > rad)
	rad = dtemp;
    }
  }

  if(rad < 0.) {
    std::cerr << funame << "no spherical surfaces found\n";
    throw Error::Init();
  }

  Plane*  pln;
  for(sit = _sur.begin(); sit != _sur.end(); ++sit) {
    pln = dynamic_cast<Plane*>((Primitive*)(*sit));
    if(pln)
      pln->set_circle(rad);
  }
}

void DivSur::PrimSet::print (std::ostream& to, const std::string& prefix) const
{
  std::vector<SharedPointer<Primitive> >::const_iterator sit;
  for(sit = _sur.begin(); sit != _sur.end(); ++sit) {
    to << prefix;
    (*sit)->print(to);
    to << "\n";
  }
}

/*********************************************************************
 **************************** Primitive Methods **************************
 *********************************************************************/

double DivSur::Primitive::distance (const Dynamic::Coordinates& dv) const 
{
  const char funame [] = "DivSur::Primitive::distance (const Dynamic::Coordinates&): ";
  std::cerr << funame << "this is a fake function, should not be here\n";
  throw Error::Logic();
}

double DivSur::Primitive::_set_rc (const Dynamic::Coordinates& dv, Dynamic::Momenta& rc_norm) const
{
  const char funame [] = "DivSur::Primitive::_set_rc (const Dynamic::Coordinates&, Array<double>&): ";
  std::cerr << funame << "this is a fake function, should not be here\n";
  throw Error::Logic();
}
 
double DivSur::Primitive::weight (const Dynamic::Coordinates& dv) const
{
  Dynamic::Momenta rc_norm;
  return _set_rc(dv, rc_norm);
}

// generalized velocities for the given temperature
double DivSur::Primitive::t_velocity (double temper, Dynamic::Vars& dv) const
{
  Dynamic::Momenta rc_norm;

  double rc_len =  _set_rc(dv, rc_norm);
  _rc2tv(rc_norm, dv, temper, dv);

  return rc_len;
}

// generalized velocities for the given energy
double DivSur::Primitive::e_velocity (double ener, Dynamic::Vars& dv) const
{
  Dynamic::Momenta rc_norm;
  double rc_len =  _set_rc(dv, rc_norm);
  double w =  _rc2ev(rc_norm, dv, ener, dv);
  if(w <= 0.)
    return -1.;
  return rc_len  * w;
}

// generalized velocities for the given energy and angular momentum 
double DivSur::Primitive::j_velocity (double ener, double amom, Dynamic::Vars& dv) const
{
  Dynamic::Momenta rc_norm;
  double rc_len =  _set_rc(dv, rc_norm);
  double w =  _rc2jv(rc_norm, dv, ener, amom, dv);
  if(w <= 0.)
    return -1.;
  return rc_len  * w;
} 

// generate random generalized velocities from the reaction coordinate vector
// for the given temperature
void DivSur::Primitive::_rc2tv(const Dynamic::Momenta& rc_norm, const Dynamic::Coordinates& dc, 
			       double temper, Dynamic::Momenta& dm)
{
  const char funame [] = "DivSur::Primitive::_rc2tv: ";

  double dtemp;
    
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    dm[i] = Random::norm();

  // for linear fragments remove the component colinear to the molecular axis
  for(int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() == Molecule::LINEAR)
      orthogonalize(dm.ang_vel(frag), dc.ang_pos(frag), 3);

  // adjust reactive component
  orthogonalize(dm.begin(), rc_norm.begin(), Dynamic::Momenta::size());
  dtemp = Random::exp();
  //dtemp = std::sqrt(2. * Random::exp());
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    dm[i] += dtemp * rc_norm[i];

  // scale to the temperature
  dm *= std::sqrt(temper);

  // converting to regular velocities

  // orbital motion
  for (int i = 0; i < 3; ++i)
    dm.orb_vel(i) /= Structure::mass_sqrt();

  // fragments rotational motion
  for(int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;
    case Molecule::LINEAR:
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(0);
      break;

    case Molecule::NONLINEAR:
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(i);
      break;
    }
}

// generate random generalized velocities from the reaction coordinate vector
// for the given energy
double DivSur::Primitive::_rc2ev(const Dynamic::Momenta& rc_norm, const Dynamic::Coordinates& dc, 
			       double ener, Dynamic::Momenta& dm)
{
  const char funame [] = "DivSur::Primitive::_rc2ev: ";

  static const double max_exp = 300.; // maximum exponent in the statistical factor

  if(ener <= 0.)
    return -1.;

  double dtemp;
    
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    dm[i] = Random::norm();

  // for linear fragments remove the component colinear to the molecular axis
  for(int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() == Molecule::LINEAR)
      orthogonalize(dm.ang_vel(frag), dc.ang_pos(frag), 3);

  // remove reaction coordinate mode
  orthogonalize(dm.begin(), rc_norm.begin(), Dynamic::Momenta::size());


  double nr_ener = ener * std::pow(Random::flat(), 2. / double(Structure::tm_dof() - 1)); // non-reactive modes
  double rc_ener = ener - nr_ener;         // energy of non-reactive modes
  double nr_len = std::sqrt(2. * nr_ener); // amplitude of non-reactive modes
  double rc_len = std::sqrt(2. * rc_ener); // reaction coordinate momentum

  // scale the amplitude of internal motion
  dm *= nr_len / vlength(dm.begin(), Dynamic::Momenta::size());

  // add reaction coordinate motion
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    dm[i] += rc_len * rc_norm[i];

  // converting to regular velocities

  // orbital motion
  for (int i = 0; i < 3; ++i)
    dm.orb_vel(i) /= Structure::mass_sqrt();

  // fragments rotational motion
  for(int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;
    case Molecule::LINEAR:
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(0);
      break;

    case Molecule::NONLINEAR:
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(i);
      break;
    }

  // statistical weight
  return std::pow(ener, double(Structure::tm_dof() - 1) / 2.);
}

// generate random generalized velocities from the reaction coordinate vector
// for the given energy and angular momentum 
double DivSur::Primitive::_rc2jv(const Dynamic::Momenta& rc_norm, const Dynamic::Coordinates& dc, 
				  double ener, double amom, Dynamic::Momenta& dm)
{
  double dtemp;
  D3::Vector vtemp, vtemp1;

  // inertia moments matrix
  Lapack::Cholesky imm(dc.imm());

  // total angular momentum vector
  Lapack::Vector ang_mom(3);
  ang_mom = 0.;
  ang_mom[2] = amom;

  // total angular velocity
  Lapack::Vector ang_vel = imm.invert(ang_mom);
  
  // rotational kinetic energy shift
  ener -= ang_mom * ang_vel / 2.;

  if(ener <= 0.)
    return -1.;

  // random generalized velocity
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    dm[i] = Random::norm();

  // for linear fragments remove the component colinear to the molecular axis
  for (int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() == Molecule::LINEAR)
      orthogonalize(dm.ang_vel(frag), dc.ang_pos(frag), 3);

  // calculate angular momentum associated with the chosen generalized random velocity
  Lapack::Vector ass_ang_mom(3);
  // orbital part
  D3::vprod(dc.orb_pos(), dm.orb_vel(), ass_ang_mom);
  ass_ang_mom *= Structure::mass_sqrt();
  // fragments rotational part
  for(int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;
    case Molecule::LINEAR:
      for (int i = 0; i < 3; ++i)
	ass_ang_mom[i] += Structure::fragment(frag).imom_sqrt(0) * dm.ang_vel(frag, i);
      break;

    case Molecule::NONLINEAR:
      for (int i = 0; i < 3; ++i)
	vtemp[i] = Structure::fragment(frag).imom_sqrt(i) * dm.ang_vel(frag, i);
      dc.mf2lf(frag, vtemp, vtemp1);
      for (int i = 0; i < 3; ++i)
	ass_ang_mom[i] += vtemp1[i];
      break;
    }

  // angular velocity associated with the chosen generalized random velocity
  Lapack::Vector ass_ang_vel = imm.invert(ass_ang_mom);
  
  // remove overall rotation from the generalized random velocity
  // orbital part
  D3::vprod(ass_ang_vel, dc.orb_pos(), vtemp);
  for(int i = 0; i < 3; ++i)
    dm.orb_vel(i) -= Structure::mass_sqrt() * vtemp[i];

  // fragments rotational part
  for(int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;

    case Molecule::LINEAR:
      for (int i = 0; i < 3; ++i)
	vtemp[i] = ass_ang_vel[i];
      orthogonalize(vtemp, dc.ang_pos(frag), 3);
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) -= Structure::fragment(frag).imom_sqrt(0) * vtemp[i];
      break;

    case Molecule::NONLINEAR:
      dc.lf2mf(frag, ass_ang_vel, vtemp);
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) -= Structure::fragment(frag).imom_sqrt(i) * vtemp[i];
      break;
    }

  // remove reactive component
  orthogonalize(dm.begin(), rc_norm.begin(), Dynamic::Momenta::size());

  // energy of internal rotations
  double in_ener = ener * std::pow(Random::flat(), 2. / double(Structure::tm_dof() - 4));
  double rc_ener = ener - in_ener;         // reaction coordinate energy
  double rc_len = std::sqrt(2. * rc_ener); // reaction coordinate momentum
  double in_len = std::sqrt(2. * in_ener); // momentum amplitude of internal rotations

  // scale the amplitude of internal motion
  dm *= in_len / vlength(dm.begin(), Dynamic::Momenta::size());

  // add reaction coordinate motion
  for (int i = 0; i < Dynamic::Momenta::size(); ++i)
    dm[i] += rc_len * rc_norm[i];

  // add overall rotation
  // orbital part
  D3::vprod(ang_vel, dc.orb_pos(), vtemp);
  for(int i = 0; i < 3; ++i)
    dm.orb_vel(i) += Structure::mass_sqrt() * vtemp[i];

  // fragments rotational part
  for(int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;

    case Molecule::LINEAR:
      for (int i = 0; i < 3; ++i)
	vtemp[i] = ang_vel[i];
      vtemp.orthogonalize(dc.ang_pos(frag));
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) += Structure::fragment(frag).imom_sqrt(0) * vtemp[i];
      break;

    case Molecule::NONLINEAR:
      dc.lf2mf(frag, ang_vel, vtemp);
      for (int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) += Structure::fragment(frag).imom_sqrt(i) * vtemp[i];
      break;
    }

  // convert to regular velocities
  // orbital motion
  for (int i = 0; i < 3; ++i)
    dm.orb_vel(i) /= Structure::mass_sqrt();

  // fragments rotational motion
  for(int frag = 0; frag < 2; ++frag)
    switch(Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;
    case Molecule::LINEAR:
      for(int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(0);
      break;

    case Molecule::NONLINEAR:
      for(int i = 0; i < 3; ++i)
	dm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(i);
      break;
    }

  //return std::pow(ener, double(Structure::tm_dof() - 4) / 2.) / std::sqrt(imm.det());
  return std::pow(ener, double(Structure::tm_dof() - 4) / 2.) / imm.det_sqrt();
}

/***********************************************************************
 ************************* Spherical facet methods *********************
 ***********************************************************************/

void DivSur::Sphere::print (std::ostream& to) const
{
  to << "Sphere = {";
  for(int frag = 0; frag < 2; ++frag) {
    to << frag << "-th pivot = " << _pivot[frag] << ", ";
  }
  to << "radius = " << _dist << "}";
}

void DivSur::Sphere::random_orient (Dynamic::Coordinates& dv) const
{
  // temporary variables
  double dtemp;
  int itemp;
  double vtemp [3];

  // randomly orient molecules
  for(int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() != Molecule::MONOATOMIC)
      Random::orient(dv.ang_pos(frag), Structure::fragment(frag).pos_size());

  // randomly orient a vector between the pivot points
  Random::orient(dv.orb_pos(), 3);
  for (int i = 0; i < 3; ++i)
    dv.orb_pos(i) *= _dist;

  for(int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() != Molecule::MONOATOMIC) {
      dv.mf2lf(frag, _pivot[frag], vtemp);
      if (frag == 0)
	for (int i = 0; i < 3; ++i)
	  dv.orb_pos(i) += vtemp[i];
      else
	for (int i = 0; i < 3; ++i)
	  dv.orb_pos(i) -= vtemp[i];
    }
}

D3::Vector DivSur::Sphere::_lf_pp_12 (const Dynamic::Coordinates& dv) const
{
  // temporary variables
  double dtemp;
  int itemp;
  D3::Vector vtemp;

  // pp2pp vector, 1->2            
  D3::Vector lf_pp_12; 

    lf_pp_12 = dv.orb_pos();

  for(int frag = 0; frag < 2; ++frag)
    if(Structure::fragment(frag).type() != Molecule::MONOATOMIC) {
      dv.mf2lf(frag, _pivot[frag], vtemp);
      if (frag == 0)
	lf_pp_12 -= vtemp;
      else
	lf_pp_12 += vtemp;
    }

  return lf_pp_12;
}

double DivSur::Sphere::_set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_norm) const
{
  // temporary variables
  double dtemp;
  int itemp;
  D3::Vector vtemp;

  // pp2pp vector, 1->2            
  D3::Vector lf_pp_12 = _lf_pp_12(dc);

  /********* reaction coordinate generalized vector ***********/

  // orbital motion degrees of freedom
  for (int i = 0; i < 3; ++i)
    rc_norm.orb_vel(i) = - lf_pp_12[i] / Structure::mass_sqrt();

  // fragments rotational degrees of freedom
  for (int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;

    case Molecule::LINEAR:
      if (frag == 1)
	D3::vprod(lf_pp_12, dc.ang_pos(frag), rc_norm.ang_vel(frag));
      else
	D3::vprod(dc.ang_pos(frag), lf_pp_12, rc_norm.ang_vel(frag));

      for (int i = 0; i < 3; ++i)
	rc_norm.ang_vel(frag, i) *= _pivot[frag][0] / Structure::fragment(frag).imom_sqrt(0);
      break;

    case Molecule::NONLINEAR:
      dc.lf2mf(frag, lf_pp_12, vtemp);
      if (frag == 1)
	D3::vprod(vtemp, _pivot[frag], rc_norm.ang_vel(frag));
      else
	D3::vprod(_pivot[frag], vtemp, rc_norm.ang_vel(frag));

      for (int i = 0; i < 3; ++i)
	rc_norm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(i);
      break;
    }

  // normalize reaction coordinate vector
  return  _dist * normalize(rc_norm.begin(), Dynamic::Momenta::size());
}

/******************************************************************
 ****************         Plane facet methods      ****************
 ******************************************************************/

void DivSur::Plane::print (std::ostream& to) const
{
    to << "Plane = {";
    for(int frag = 0; frag < 2; ++frag) {
	to << frag << "-th fragment ";
	if(frag == _plane_frag)
	    to << "plane = " << _plane;
	else
	    to << "pivot = " << _pivot;
	if(!frag)
	    to << ", ";
    }
    to << "}";
}

void DivSur::Plane::random_orient (Dynamic::Coordinates& dv) const
{
  const char funame [] = "DivSur::Plane::random_orient: ";

  if(_radius < 0.) {
    std::cerr << funame << "sampling circle was not set up\n";
    throw Error::Init();
  }

  // temporary variables
  double dtemp, res;
  int itemp;
  double vtemp [3];

  // randomly orient molecules
  for (int frag = 0; frag < 2; ++frag)
    if (Structure::fragment(frag).type() != Molecule::MONOATOMIC)
      Random::orient(dv.ang_pos(frag), Structure::fragment(frag).pos_size());

  // random point on the unit circle
  double v [2];
  do {
    res = 0.;
    for(int i = 0; i < 2; ++i) {
      dtemp = 2. * Random::flat() - 1.;
      res += dtemp * dtemp;
      v[i] = dtemp;
    }
  } while (res > 1.);

  D3::Vector n[2];

  switch (Structure::fragment(_plane_frag).type()) {
  case Molecule::MONOATOMIC:
    std::cerr << funame << "monoatomic molecule cannot have plane elements\n";
    throw Error::Logic();
  case Molecule::LINEAR:
    // plane normal in the lab frame
    for (int i = 0; i < 3; ++i)
      dv.orb_pos(i) = dv.ang_pos(_plane_frag, i);
    // vectors orthogonal to the plane normal in the lab frame
    find_orth(dv.orb_pos(), n);
    // random point on the plane within the given circle in the lab frame
    for (int i = 0; i < 3; ++i) {
      dv.orb_pos(i) *= _plane.dist();
      for(int j = 0; j < 2; ++j)
	dv.orb_pos(i) +=  _radius * v[j] * n[j][i]; 
    }
    break;
  case Molecule::NONLINEAR:
    // random point on the plane within the given circle in the molecular frame
    for (int i = 0; i < 3; ++i) {
      vtemp[i] = _plane.dist() * _plane.normal()[i];
      for(int j = 0; j < 2; ++j)
	vtemp[i] += _radius * v[j] * _plane.orth(j)[i];
    }
    // transform to the lab frame
    dv.mf2lf(_plane_frag, vtemp, dv.orb_pos());
    break;
  }

  const int  _point_frag = 1 - _plane_frag;
  if(Structure::fragment(_point_frag).type() != Molecule::MONOATOMIC) {
    dv.mf2lf(_point_frag, _pivot, vtemp);
    for (int i = 0; i < 3; ++i)
      dv.orb_pos(i) -= vtemp[i];
  }

  if (_plane_frag)
    for (int i = 0; i < 3; ++i)
      dv.orb_pos(i) = - dv.orb_pos(i);
}

double DivSur::Plane::distance (const Dynamic::Coordinates& dv) const
{
  const char funame [] = "DivSur::Plane::distance: ";

  D3::Vector vtemp;

  // plane in the lab frame
  D3::Vector lf_plane_normal;

  switch (Structure::fragment(_plane_frag).type()) {
  case Molecule::MONOATOMIC:
    std::cerr << funame << "monoatomic molecule cannot have plane elements\n";
    throw Error::Logic();
  case Molecule::LINEAR:
    lf_plane_normal = dv.ang_pos(_plane_frag);
    break;
  case Molecule::NONLINEAR:
    dv.mf2lf(_plane_frag, _plane.normal(), lf_plane_normal);
    break;
  }

  D3::Vector lf_pivot; // pivot point in the lab frame relatively to the plane fragment
  if (!_plane_frag)
    for (int i = 0; i < 3; ++i)
      lf_pivot[i] =   dv.orb_pos(i);
  else
    for (int i = 0; i < 3; ++i)
      lf_pivot[i] = - dv.orb_pos(i);

  int  _point_frag = 1 - _plane_frag;
  if(Structure::fragment(_point_frag).type() != Molecule::MONOATOMIC) {
    dv.mf2lf(_point_frag, _pivot, vtemp);
    lf_pivot += vtemp;
  }

#ifdef DEBUG
  dv.lf2mf(_plane_frag, lf_pivot, vtemp);
#endif

  return vdot(lf_pivot, lf_plane_normal) - _plane.dist();
}

double DivSur::Plane::_set_rc (const Dynamic::Coordinates& dv, Dynamic::Momenta& rc_norm) const
{
  const char funame [] = "DivSur::Plane::_set_rc: ";

  // temporary variables
  double dtemp;
  int itemp;
  D3::Vector vtemp;

  // plane normal in the lab frame
  D3::Vector lf_plane_normal;

  switch (Structure::fragment(_plane_frag).type()) {
  case Molecule::MONOATOMIC:
    std::cerr << funame << "monoatomic molecule cannot have plane elements\n";
    throw Error::Logic();
  case Molecule::LINEAR:
    lf_plane_normal = dv.ang_pos(_plane_frag);
    break;
  case Molecule::NONLINEAR:
    dv.mf2lf(_plane_frag, _plane.normal(), lf_plane_normal);
    break;
  }

  D3::Vector lf_pivot; // pivot point in the lab frame relatively to the plane fragment, r + p
  if (!_plane_frag)
    for (int i = 0; i < 3; ++i)
      lf_pivot[i] =   dv.orb_pos(i);
  else
    for (int i = 0; i < 3; ++i)
      lf_pivot[i] = - dv.orb_pos(i);

  int  _point_frag = 1 - _plane_frag;
  if(Structure::fragment(_point_frag).type() != Molecule::MONOATOMIC) {
    dv.mf2lf(_point_frag, _pivot, vtemp);
    lf_pivot += vtemp;
  }

  /********* reaction coordinate generalized vector ***********/

  // orbital motion degrees of freedom
  for (int i = 0; i < 3; ++i)
    rc_norm.orb_vel(i) = lf_plane_normal[i] / Structure::mass_sqrt();

  // fragments rotational degrees of freedom
  for (int frag = 0; frag < 2; ++frag)
    switch (Structure::fragment(frag).type()) {
    case Molecule::MONOATOMIC:
      break;
    case Molecule::LINEAR:
      if (frag == _plane_frag)
	D3::vprod(lf_plane_normal, lf_pivot, rc_norm.ang_vel(frag));
      else {
	D3::vprod(dv.ang_pos(frag), lf_plane_normal, rc_norm.ang_vel(frag));
	for (int i = 0; i < 3; ++i)
	  rc_norm.ang_vel(frag, i) *= _pivot[0];
      }
      for (int i = 0; i < 3; ++i)
	rc_norm.ang_vel(frag, i) /=  Structure::fragment(frag).imom_sqrt(0);
      break;
    case Molecule::NONLINEAR:
      if (frag == _plane_frag) {
	dv.lf2mf(frag, lf_pivot, vtemp);
	D3::vprod(_plane.normal(), vtemp, rc_norm.ang_vel(frag));
      }
      else {
	dv.lf2mf(frag, lf_plane_normal, vtemp);
	D3::vprod(_pivot, vtemp, rc_norm.ang_vel(frag));
      }
      for (int i = 0; i < 3; ++i)
	rc_norm.ang_vel(frag, i) /= Structure::fragment(frag).imom_sqrt(i);
      break;
    }

  // normalize reaction coordinate vector
  return normalize(rc_norm.begin(), Dynamic::Momenta::size()) * _radius * _radius / 4.;
}

/***********************************************************************
 ******************** Basic Dividing Surface Methods *******************
 ***********************************************************************/

void DivSur::OneSur::read (std::istream& from) throw(Error::General)
{
  const char funame [] = "DivSur::OneSur::read: ";

  if(_isinit) {
    std::cerr << funame << "already initialized\n";
    throw Error::Init();
  }

  _isinit = true;

  if(!Structure::isinit()) {
    std::cerr << funame << "structure not initialized yet\n";
    throw Error::Init();
  }
  
  KeyGroup BasicSurface;
  Key  sur_key ("Species");

  // read primitives
  PrimSet::read(from);

  std::string token;
  while(from >> token) {
    // input end
    if(token == IO::end_key()) {
      break;
    }
    // species (region) definition
    else if(token == sur_key) {
      if((bool)_prod) {
	std::cerr << funame << token << ": already initialized\n";
	throw Error::Init();
      }
      IO::LineInput lin(from);
      _prod = Logical::read_expr(lin);
      _prod->init(_fmap);
    }
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }
  }

  if(!from) {
    std::cerr << funame << "input corrupted\n";
    throw Error::Input();
  }

  if(!_prod) {
    std::cerr << funame << "species not defined\n";
    throw Error::Init();
  }
}

void DivSur::OneSur::print (std::ostream& to, const std::string& prefix) const
{
  to << prefix << "Primitive surfaces:\n";
  PrimSet::print(to, prefix + "   ");
  to << "\n";
}

// test if the configuration is inside the surface
bool DivSur::OneSur::test (const Dynamic::Coordinates& dc) const
  throw(Error::General)
{
  const char funame [] = "DivSur::OneSur::test: ";

  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not yet been initialized\n";
    throw Error::Init();
  }

    std::vector<bool> facet_test(size());
    for(int i = 0; i < size(); ++i)
      facet_test[i] = PrimSet::test(i, dc);
    return _prod->evaluate(facet_test);
}

DivSur::OneSur::SmpRes DivSur::OneSur::facet_test (int face, const Dynamic::Coordinates& dc) const
  throw(Error::General)
{
  const char funame [] = "DivSur::OneSur::facet_test: ";

  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not yet been initialized\n";
    throw Error::Init();
  }

  SmpRes res;

  if(dc.are_atoms_close()) {
    res.stat = SmpRes::CLOSE;
    return res;
  }

  if(Dynamic::exclude_region && Dynamic::exclude_region->test(dc)) {
    res.stat = SmpRes::EXCLUDE;
    return res;
  }

  std::vector<bool> facet_test(size());
  for(int i = 0; i < size(); ++i)
    if(i != face)
      facet_test[i] = PrimSet::test(i, dc);

  // Boundary test:
  // configuration is in the inner part of the configurational 
  // space in relation to the facet
  facet_test[face] = true;
  res.inner = _prod->evaluate(facet_test);
  // configuration is in the outer part of the configurational 
  // space in relation to the facet
  facet_test[face] = false;
  bool outer = _prod->evaluate(facet_test);

  if(res.inner == outer)
    res.stat = SmpRes::INNER;
  else
    res.stat = SmpRes::FACET;

  return res;
}

// facet to which the configuration is closest
int DivSur::OneSur::classify (const Dynamic::Coordinates& dc) const
  throw(Error::General)
{
  const char funame [] = "DivSur::OneSur::classify: ";

  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not yet been initialized\n";
    throw Error::Init();
  }

  double dist_min, dtemp;
  int    face_min = -1;
  for(int face = 0; face < size(); ++face) {
    dtemp = distance(face, dc);
    dtemp = dtemp >= 0. ? dtemp : -dtemp;
    if(!face || dtemp < dist_min) {
      dist_min = dtemp;
      face_min = face;
    } 
  }
  return face_min;
}

/****************************************************************************
 *                Multiple species dividing surface: MultiSur               *
 ****************************************************************************/

void DivSur::MultiSur::read (std::istream& from) throw(Error::General)
{
  const char funame [] = "DivSur::MultiSur::read: ";

  if(_isinit) {
    std::cerr << funame << "already initialized\n";
    throw Error::Init();
  }

  _isinit = true;

  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not yet been initialized\n";
    throw Error::Init();
  }

  int itemp;
  std::string comment, line;

  DivSur::PrimSet::read(from);

  KeyGroup MultipleSpeciesSurface;
  Key spec_key ("Species");

  std::string token;
  while(from >> token) {// read cycle
    // end input
    if(token == IO::end_key()) {
      break;
    }
    // species
    else if(token == spec_key) {
      if(_species.size()) {
	std::cerr << funame << token << ": already defined\n";
	throw Error::Init();
      }
      
      if(!(from >> itemp)) {
	std::cerr << funame << token << ": corrupted\n";
	throw Error::Input();
      }
      std::getline(from, comment);

      if(itemp < 1) {
	std::cerr << funame << token << ": out of range\n";
	throw Error::Range();
      }
      _species.resize(itemp);

      for(int s = 0; s < _species.size(); ++s) {
	IO::LineInput lin(from);
	_species[s] = Logical::read_expr(lin);
	_species[s]->init(_fmap);
      }
    }// species
    else {
      std::cerr << funame << "unknown keyword: " << token << "\n";
      Key::show_all(std::cerr);
      std::cerr << "\n";
      throw Error::Form();
    }
  }// read cycle

  // checking
  if(!from) {
    std::cerr << funame << "stream is corrupted\n";
    throw Error::Input();
  }
    
  if(!_species.size()) {
    std::cerr << funame << "species not defined\n";
    throw Error::Init();
  }

  //
  // check if the species do not intersect
  //
  std::vector<bool> pset(primitive_size(), false);
  std::vector<bool> spec_res(_species.size());
  int pi;
  int spec_num;

  // Checking species for logical intersections;
  while(1) {//primitives configuration cycle
    spec_num = 0;
    for(int si = 0; si < _species.size(); ++si)
      if(spec_res[si] = _species[si]->evaluate(pset))
	spec_num++;
    if(spec_num > 1) {
      std::cerr << funame << ": WARNING: the primitives logical configuration (";
      for(pi = 0; pi < primitive_size(); ++pi)
	std::cerr << pset[pi];
      std::cerr << ") is compatible with multiple species (";
      for(int si = 0; si < _species.size(); ++si)
	std::cerr << spec_res[si];
      std::cerr << ")\n";      
    }

    for(pi = 0;  pi < primitive_size(); ++pi) {
      pset[pi] = !pset[pi];
      if(pset[pi])
	break;
    }
    if(pi == primitive_size())
      break;
  }//primitives configuration cycle
}

void DivSur::MultiSur::print (std::ostream& to, const std::string& prefix) const
{
  to << prefix << "Primitive surfaces:\n";
  PrimSet::print(to, prefix + "   ");
  to << "\n";
}

int DivSur::MultiSur::classify  (const Dynamic::Coordinates& dc) const 
  throw(Error::General)
{
  const char funame [] = "DivSur::MultiSur::classify: ";

  std::vector<bool> surface_test(primitive_size());
  for(int i = 0; i < primitive_size(); ++i)
    surface_test[i] = DivSur::PrimSet::test(i, dc);

  std::vector<int> res;
  for(int s = 0; s < _species.size(); ++s)
    if(_species[s]->evaluate(surface_test))
      res.push_back(s);
  
  if(!res.size())
    return _species.size();

  if(res.size() == 1)
    return *res.begin();

  std::cerr << funame << "configuration belongs to more than one species: ";
  for(int i = 0; i < res.size(); ++i)
    std::cerr << res[i] << " ";
  std::cerr << "\n";
  throw Error::Logic();
}

bool DivSur::MultiSur::species_test (int spec, const Dynamic::Coordinates& dc) const
  throw(Error::General)
{
  int res = classify(dc);
  if(res == spec)
    return true;
  return false;   
}

DivSur::MultiSur::SmpRes DivSur::MultiSur::facet_test (int sur, const Dynamic::Coordinates& dc) const 
  throw(Error::General)
{
  const char funame [] = "DivSur::MultiSur::facet_test: ";

  if(!Structure::isinit()) {
    std::cerr << funame << "structure has not yet been initialized\n";
    throw Error::Init();
  }

  SmpRes res;
    
  // check if the configuration belongs to the excluded region
  if(Dynamic::exclude_region && Dynamic::exclude_region->test(dc)) {
    res.stat = SmpRes::EXCLUDE;
    return res;
  }

  // check if the atoms are too close
  if(dc.are_atoms_close()) {
    res.stat = SmpRes::CLOSE;
    return res;
  }

  // to which actual facet does the configuration belong 
  std::vector<bool> surface_test(primitive_size());
  for(int i = 0; i < primitive_size(); ++i)
    if(i != sur)
      surface_test[i] = DivSur::PrimSet::test(i, dc);

  // species tests
  std::vector<bool> negative_test(_species.size());
  surface_test[sur] = false;
  for(int s = 0; s < _species.size(); ++s)
    negative_test[s] = _species[s]->evaluate(surface_test);

  std::vector<bool> positive_test(_species.size());
  surface_test[sur] = true;
  for(int s = 0; s < _species.size(); ++s)
    positive_test[s] = _species[s]->evaluate(surface_test);

  // analyze the results
  std::set<int> inner, from, to;
  for(int s = 0; s < _species.size(); ++s)
    if(positive_test[s] && negative_test[s])
      inner.insert(s);
    else if(!positive_test[s] && negative_test[s])
      from.insert(s);
    else if(positive_test[s] && !negative_test[s])
      to.insert(s);

  // inner region of some species
  if((inner.size() == 1 || !inner.size()) && !to.size() && !from.size()) {
    res.stat = SmpRes::INNER;
    return res;
  }

  if(!inner.size()  && to.size() == 1 && !from.size()) {
    // from bimolecular products to a bound species
    res.face.first  = _species.size();
    res.face.second =     *to.begin();
  }
  else if(!inner.size()  && from.size() == 1 && !to.size()) {
    // from a bound species to bimolecular products
    res.face.first   =   *from.begin();
    res.face.second  = _species.size();
  }
  else if(!inner.size()  && from.size() == 1 && to.size() == 1) {
    // from a bound species to another bound species
    res.face.first   = *from.begin();
    res.face.second  =   *to.begin();
  }
  else {
    // logic failure
    std::cerr << funame << "Sampling error for "
	      << sur << "-th surface; dumping the results:\n";

    std::cerr << funame 
	      << "   Configuration is inside of the following species: ";
    for(std::set<int>::iterator it = inner.begin(); it != inner.end(); ++it) {
      if(it != inner.begin())
	std::cerr << ", ";
      std::cerr << *it;
    }
    std::cerr << "\n";

    std::cerr << funame 
	      << "   Configuration is at the interface of the following species: ";

    std::cerr << funame << "      From negative side: ";
    for(std::set<int>::iterator it = from.begin(); it != from.end(); ++it) {
      if(it != from.begin())
	std::cerr << ", ";
      std::cerr << *it;
    }
    std::cerr << "\n";
	
    std::cerr << funame << "      From positive side: ";
    for(std::set<int>::iterator it = to.begin(); it != to.end(); ++it) {
      if(it != to.begin())
	std::cerr << ", ";
      std::cerr << *it;
    }
    std::cerr << "\n";

    std::cerr << "\n";

    res.stat = SmpRes::FAIL;
    return res;
  }

  res.stat = SmpRes::FACET;

  return res;
}

