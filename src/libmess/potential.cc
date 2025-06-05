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

#include "potential.hh"
#include "units.hh"
#include "key.hh"
#include "limits.hh"

#include <sstream>

/*********************************************************************
 ************************ POTENTIAL WRAPPER **************************
 *********************************************************************/

Potential::Wrap default_pot;

void Potential::Wrap::read (std::istream& from) 
{
  const char funame [] = "Potential::Wrap::read: ";

  KeyGroup PotentialWrap;

  Key    tn_key("TorqueNumeric"  );
  Key    cl_key("ChargeLinear"   );
  Key    cn_key("ChargeNonlinear");
  Key    dd_key("DipoleDipole"   );
  Key multi_key("Multipole"      );

  if(_fun) {
    //
    std::cerr << funame << "potential has been already initialized\n";

    throw Error::Init();
  }

  std::string token, comment;

  while(from >> token) {
    //
    // torque-numeric
    //
    if(tn_key == token) {
      //
      std::getline(from, comment);

      _fun = ConstSharedPointer<Base>(new TorqueNumeric(from));
      
      return;
    }
    // charge-linear molecule
    //
    else if(cl_key == token) {
      //
      std::getline(from, comment);

      _fun = ConstSharedPointer<Base>(new ChargeLinear(from));
      return;
    }
    // charge-nonlinear molecule
    //
    else if(cn_key == token) {
      //
      std::getline(from, comment);

      _fun = ConstSharedPointer<Base>(new ChargeNonlinear(from));
      return;
    }
    // dipole-dipole
    //
    else if(dd_key == token) {
      //
      std::getline(from, comment);

      _fun = ConstSharedPointer<Base>(new DipoleDipole(from));
      return;
    }
    // multipole
    //
    else if(multi_key == token) {
      //
      std::getline(from, comment);

      _fun = ConstSharedPointer<Base>(new Multipole(from));
      return;
    }
    // unknown keyword
    //
    else {
      //
      std::cerr << funame << "unknown keyword " << token << "\n";

      Key::show_all(std::cerr);
      std::cerr << "\n";

      throw Error::Init();
    }
  }

  std::cerr << funame << "stream is corrupted\n";
  throw Error::Init();
  
}

Potential::TorqueNumeric::TorqueNumeric (std::istream& from) : _dist_incr(1.e-4),  _ang_incr(1.e-4)
{
  const char funame [] = "Potential::TorqueNumeric::TorqueNumeric: ";

  double dtemp;
  int    itemp;

  KeyGroup TorqueNumericGroup;

  Key dist_key("DistanceIncrement"  );
  Key angl_key("AngularIncrement"   );
  Key bare_key("BareEnergyMethod"   );
  
  std::string token, line, comment, stemp;
  
  // read cycle
  //
  while(from >> token) {
    //
    // input end
    //
    if(token == IO::end_key()) {
      //
      std::getline(from, comment);
      
      break;
    }
    // bare energy potential
    //
    else if(token == bare_key) {
      //
      if(_bare_ener) {
	//
	std::cerr << funame << token << ": already initialized\n";

	throw Error::Init();
      }

      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      if(stemp == "sjk") {
	//
	_bare_ener = ConstSharedPointer<BareEnergy>(new SJK(from));
      }
      // unknown bare energy potential type
      //
      else {
	//
	std::cerr << funame << token << ": unknown bare energy potential type: " << stemp << ": available types: sjk\n";

	throw Error::Init();
      }
    }
    // distance displacement
    //
    else if(token == dist_key) {
      //
      if(!(from >> _dist_incr)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      
      if(_dist_incr <= 0.) {
	//
	std::cerr << funame << token << ": out of range\n";
	
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // angular displacement
    //
    else if(token == angl_key) {
      //
      if(!(from >> _ang_incr)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Init();
      }

      if(_ang_incr <= 0.) {
	//
	std::cerr << funame << token << ": out of range\n";
	
	throw Error::Range();
      }

      std::getline(from, comment);
    }
    // unknown key
    //
    else {
      //
      std::cerr << funame << "unknown key: " << token << "\n";
      
      Key::show_all(std::cerr);
      std::cerr << "\n";
      
      throw Error::Init();
    }
  }// read cycle
  
  // checking input
  //
  if(!from) {
    //
    std::cerr << funame << "input stream is corrupted\n";
    
    throw Error::Form();
  }

  if(!_bare_ener) {
    //
    std::cerr << funame << "no bare energy calculation method was provided\n";
    
    throw Error::Init();
  }
}

void Potential::TorqueNumeric::_set_torque (const Dynamic::Coordinates& _dc, D3::Vector* torque) const 
{
  const char funame [] = "Potential::TorqueNumeric::_set_torque: ";

  D3::Vector& force = *torque;

  torque += 1;

  Dynamic::Coordinates dc = _dc;

  double ener [2];
  
  for(int i = 0; i < 3; ++i) {
    //
    for(int j = 0; j < 2; ++j) {
      //
      dc.orb_pos()[i] = _dc.orb_pos()[i] + ((double)j - 0.5) * _dist_incr;

      ener[j] = _bare_ener->evaluate(dc);
    }

    force[i] = (ener[0] - ener[1]) / _dist_incr;
    
    dc.orb_pos()[i] = _dc.orb_pos()[i];
  }

  double q[4];

  double* v = q + 1;
  
  const double* n;
  
  for(int frag = 0; frag < 2; ++frag) {
    //
    switch(Structure::type(frag)) {
      //
    case Molecule::MONOATOMIC:
      //
      torque[frag] = 0.;

      break;
      //
    case Molecule::LINEAR:
      //
      for(int i = 0; i < 3; ++i) {
	//
	v[i] = 0.;

	// right sign?
	//
	v[(i + 1) % 3] = -_dc.ang_pos(frag)[(i + 2) % 3];
	v[(i + 2) % 3] =  _dc.ang_pos(frag)[(i + 1) % 3];

	for(int j = 0; j < 2; ++j) {
	  //
	  for(int k = 0; k < 3; ++k)
	    //
	    dc.write_ang_pos(frag)[k] = _dc.ang_pos(frag)[k] + ((double)j - 0.5) * v[k] * _ang_incr;

	  dc.normalize(frag);
	  
	  ener[j] = _bare_ener->evaluate(dc);
	}

	torque[frag][i] = (ener[0] - ener[1]) / _ang_incr;

	for(int k = 0; k < 3; ++k)
	  //
	  dc.write_ang_pos(frag)[k] = _dc.ang_pos(frag)[k];
      }
	
      break;
      //
    case Molecule::NONLINEAR:
      //
      n = _dc.ang_pos(frag) + 1;
      
      for(int i = 0; i < 3; ++i) {
	//
	q[0] = -n[i];

	v[i] = _dc.ang_pos(frag)[0];

	// right sign ?
	//
	v[(i + 1) % 3] =  n[(i + 2) % 3];
	v[(i + 2) % 3] = -n[(i + 1) % 3];

	for(int j = 0; j < 2; ++j) {
	  //
	  for(int k = 0; k < 4; ++k)
	    //
	    dc.write_ang_pos(frag)[k] = _dc.ang_pos(frag)[k] + ((double)j - 0.5) * q[k] * _ang_incr / 2.;

	  dc.normalize(frag);
	  
	  ener[j] = _bare_ener->evaluate(dc);
	}

	torque[frag][i] = (ener[0] - ener[1]) / _ang_incr;

	for(int k = 0; k < 4; ++k)
	  //
	  dc.write_ang_pos(frag)[k] = _dc.ang_pos(frag)[k];
      }
    }
  }

#ifdef DEBUG
  /*  
  // test that total torque is zero
  // 
  D3::Vector tot_torque, orb_torque, lab_torque[2];

  for(int f = 0; f < 2; ++f)
    //
    if(Structure::type(f) == Molecule::NONLINEAR) {
      //
      _dc.mf2lf(f, (const double*)torque[f], (double*)lab_torque[f]);
    }
    else
      //
      lab_torque[f] = torque[f];
 
  D3::vprod(_dc.orb_pos(), (const double*)force, (double*)orb_torque);

  tot_torque = orb_torque + lab_torque[0] + lab_torque[1];

  std::cerr << " 0-th torque = " << std::setw(12) << lab_torque[0].vlength()
    //
	    << " 1-st torque = " << std::setw(12) << lab_torque[1].vlength()
    //
	    << " orb. torque = " << std::setw(12) << orb_torque.vlength()
    //
	    << " tot. torque = " << std::setw(12) << tot_torque.vlength()
    //
	    << "\n";
  */
#endif
}

double Potential::TorqueNumeric::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const 
{
  static const char funame [] = "Potential::TorqueNumeric::operator(): ";
 
  _increment_count();
 
  if(torque)
    //
    _set_torque(dc, torque);

  return _bare_ener->evaluate(dc);
}

/*******************************************************************************************************
 ********************************** SJK-STYLE BARE ENERGY POTENTIAL ************************************
 *******************************************************************************************************/

Potential::SJK::SJK (std::istream& from) : _pot_form(SJK_FORM), _pot_ener(0), _pot_init(0), _corr_ener(0), _corr_init(0)
{    
  const char funame [] = "Potential::SJK::SJK: ";

  double dtemp;
  int    itemp;
  
  KeyGroup SjkGroup;

  Key  pot_form_key("Format");
  
  Key  pot_libr_key("Library");
  Key  pot_ener_key("Symbol");
  Key  pot_init_key("InitMethod");
  Key  pot_data_key("InitData");
  Key  pot_rpar_key("ParameterReal");
  Key  pot_ipar_key("ParameterInteger");

  Key corr_libr_key("CorrectionLibrary");
  Key corr_ener_key("CorrectionSymbol");
  Key corr_init_key("CorrectionInitMethod");
  Key corr_data_key("CorrectionInitData");
  Key corr_rpar_key("CorrectionParameterReal");
  Key corr_ipar_key("CorrectionParameterInteger");

  std::string token, line, comment, stemp;

  std::string pot_data, corr_data;

  // read cycle
  //
  while(from >> token) {
    //
    // input end
    //
    if(token == IO::end_key()) {
      //
      std::getline(from, comment);
      
      break;
    }
    // potential format
    //
    else if(token == pot_form_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";

	throw Error::Input();
      }
      std::getline(from, comment);

      if(stemp == "sjk") {
	//
	_pot_form = SJK_FORM;
      }
      // unknown format
      //
      else {
	//
	std::cerr << funame << "unknown format: " << stemp << ": available formats: sjk\n";
	
	throw Error::Range();
      }
    }
    // potential library name
    //
    else if(token == pot_libr_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      _pot_libr.open(stemp);
    }
    // potential energy method
    //
    else if(token == pot_ener_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      _pot_ener = _pot_libr.member(stemp);
    }
    // potential initialization method
    //
    else if(token == pot_init_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      _pot_init = (init_t)_pot_libr.member(stemp);
    }
    // potential initialization data
    //
    else if(token == pot_data_key) {
      //
      if(!(from >> pot_data)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // potential real parameters
    //
    else if(token == pot_rpar_key) {
      //
      IO::LineInput lin(from);
      std::vector<double> rpar;

      while(lin >> dtemp)
	//
	rpar.push_back(dtemp);

      if(!rpar.size()) {
	//
	std::cerr << funame << token << ": empty\n";

	throw Error::Init();
      }

      _pot_rpar = rpar;
    }
    // potential integer parameters
    //
    else if(token == pot_ipar_key) {
      //
      IO::LineInput lin(from);
      std::vector<int> ipar;

      while(lin >> itemp)
	//
	ipar.push_back(itemp);

      if(!ipar.size()) {
	//
	std::cerr << funame << token << ": empty\n";

	throw Error::Init();
      }

      _pot_ipar = ipar;
    }
    // correction library name
    //
    else if(token == corr_libr_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      _corr_libr.open(stemp);
    }
    // correction energy method
    //
    else if(token == corr_ener_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      _corr_ener = _corr_libr.member(stemp);
    }
    // correction initialization method
    //
    else if(token == corr_init_key) {
      //
      if(!(from >> stemp)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);

      _corr_init = (init_t)_corr_libr.member(stemp);
    }
    // correction initialization data
    //
    else if(token == corr_data_key) {
      //
      if(!(from >> corr_data)) {
	//
	std::cerr << funame << token << ": is corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // correction real parameters
    //
    else if(token == corr_rpar_key) {
      //
      IO::LineInput lin(from);
      std::vector<double> rpar;

      while (lin >> dtemp)
	//
	rpar.push_back(dtemp);

      if(!rpar.size()) {
	//
	std::cerr << funame << token << ": empty\n";

	throw Error::Init();
      }

      _corr_rpar = rpar;
    }
    // correction integer parameters
    //
    else if(token == corr_ipar_key) {
      //
      IO::LineInput lin(from);
      std::vector<int> ipar;

      while (lin >> itemp)
	//
	ipar.push_back(itemp);

      if(!ipar.size()) {
	//
	std::cerr << funame << token << ": empty\n";

	throw Error::Init();
      }

      _corr_ipar = ipar;
    }
    // unknown key
    //
    else {
      //
      std::cerr << funame << "unknown key: " << token << "\n";
      
      Key::show_all(std::cerr);
      std::cerr << "\n";
      
      throw Error::Init();
    }
  }// read cycle

  // checking input
  //
  if(!from) {
    //
    std::cerr << funame << "input stream is corrupted\n";
    
    throw Error::Form();
  }

  if(!_pot_ener) {
    //
    std::cerr << funame << "no energy calculation method was provided\n";
    
    throw Error::Init();
  }

  if(_pot_init)
    //
    _pot_init(pot_data.c_str());

  if(_corr_init)
    //
    _corr_init(corr_data.c_str());
}

double Potential::SJK::evaluate (const Dynamic::Coordinates& dc) const 
{
  const char funame [] = "Potential::SJK::evaluate: ";

  std::vector<D3::Vector> rel_pos [2];
  
  for(int f = 0; f < 2; ++f)
    //
    dc.set_rel_pos(f, rel_pos[f]);
  
  double res = 0.;

  if(_format() == SJK_FORM) { 
    //
    Array_3<double> r(2, 100, 3);
    
    for(int f = 0; f < 2; ++f)
      //
      for(int a = 0; a < Structure::fragment(f).size(); ++a)
	//
        for(int i = 0; i < 3; ++i)
	  //
	  if(!f) {
	    //
	    r(f, a, i) = rel_pos[f][a][i];
	  }
	  else
	    //
	    r(f, a, i) = rel_pos[f][a][i] + dc.orb_pos()[i];
    
    res = ((sjk_t)_pot_ener)(r,  _pot_rpar,  _pot_ipar);

    if(_corr_ener)
      //
      res += ((sjk_t)_corr_ener)(r, _corr_rpar, _corr_ipar);
  }

  return res;
}

/********************************************************************************************
 ************************************** YG POTENTIAL ****************************************
 ********************************************************************************************/

void Potential::YG::Base::read_coefficients (const std::string& file)
{
  const char funame [] = "Potential::YG::Base::read_coefficients: ";

  double dtemp;
  int    itemp;
  
  assert();

  std::ifstream from(file.c_str());

  if(!from) {
    //
    std::cerr << funame << "cannot open " << file << " file\n";

    throw Error::Input();
  }

  std::string token, comment;
  
  int basis_size, dist_size, index_size;

  IO::LineInput lin(from);

  if(!(lin >> basis_size >> dist_size >> index_size)) {
    //
    std::cerr << funame << "cannot read # of basis functions OR # of distances OR basis function index size\n";

    throw Error::Input();
  }

  if(basis_size < 2 || dist_size < 2 || index_size < 1) {
    //
    std::cerr << funame << "# of basis functions OR # of distances OR basis function index size out of range: "
      //
	      << basis_size << ", " << dist_size << ", " << index_size << "\n";

    throw Error::Range();
  }

  Array<double> dist(dist_size);
  Array<double> coef(dist_size);
  
  std::vector<int>   index(index_size);
  
  // read distances in bohr
  //
  std::getline(from, comment);

  for(int i = 0; i < dist_size; ++i) {
    //
    if(!(from >> dtemp)) {
      //
      std::cerr << funame << "cannot read " << i << "-th distance\n";

      throw Error::Input();
    }

    if(dtemp <= 0.) {
      //
      std::cerr << funame << "negative distance: " << dtemp << ", " << i << "\n";

      throw Error::Range();
    }
    else if(i && dtemp <= dist[i - 1]) {
      //
      std::cerr << funame << "distances are not in increasing order: " << i << ", " << dtemp << "\n";

      throw Error::Logic();
    }
  }

  std::getline(from, comment);

  for(int b = 0; b < basis_size; ++b) {
    //
    // empty line
    //
    std::getline(from, comment);
    
    // read basis function composite index
    //
    lin.read_line(from);

    for(int i = 0; i < index.size(); ++i) {
      //
      if(!(lin >> itemp)) {
	//
	std::cerr << funame << b << "-th basis function: cannot read " << i << "-th index value\n";

	throw Error::Input();
      }
  
      if(itemp < 0) {
	//
	std::cerr << funame << b << "-th basis function: " << i << "-th index value out of range: "  << itemp << "\n";

	throw Error::Input();
      }

      index[i] = itemp;
    }

    if(_coef.find(index) != _coef.end()) {
      //
      std::cerr << funame << b << "-th basis function: index already has been used\n";

      throw Error::Logic();
    }
    
    // read basis function coefficient distance dependence
    //
    for(int i = 0; i < dist_size; ++i)
      //
      if(!(from >> coef[i])) {
	//
	std::cerr << funame << b << "-th basis function: cannot read " << i << "-th coefficient\n";

	throw Error::Input();
      }

    std::getline(from, comment);
    
    _coef[index].init(dist, coef, dist_size);
  }

  _dist_min = dist.front();
  _dist_max = dist.back();
}

double Potential::YG::Base::evaluate (const Dynamic::Coordinates& dc) const
{
  const char funame [] = "Potential::YG::Base::evaluate: ";

  assert();

  double dtemp;
  int    itemp;
      
  Array<double> vtemp(4);
      
  if(!size()) {
    //
    std::cerr << funame << "basis function coefficients have not been initialized\n";
      
    throw Error::Init();
  }
      
  Array<double> state(3 + Structure::pos_size(1));

  // orbital part of the state vector
  //
  dc.lf2mf(0, dc.orb_pos(), (double*)state);
      
  const double dist = ::normalize((double*)state, 3);

  if(dist < _dist_min || dist > _dist_max) {
    //
    std::cerr << funame << "interfragment distance is out of range: " << dist << ": [" << _dist_min << ", " << _dist_max << "]\n";

    throw Error::Range();
  }

  // angular part of the state vector
  //
  switch(Structure::type(1)) {
    //
  case Molecule::LINEAR:
    //
    dc.lf2mf(0, dc.ang_pos(1), (double*)state +3);

    break;
    //
  case Molecule::NONLINEAR:
    //
    // inverse first fragment orientation
    //
    for(int i = 0; i < 4; ++i)
      //
      if(!i) {
	//
	vtemp[i] = dc.ang_pos(0)[i];
      }
      else
	//
	vtemp[i] = -dc.ang_pos(0)[i];

    // 2nd fragment orientation in the reference frame of the first fragment: q_0^{-1}q_1
    //
    Quaternion::qprod(vtemp, dc.ang_pos(1), (double*)state + 3);
  }
      
  double res = 0.;
      
  for(coef_t::const_iterator cit = _coef.begin(); cit != _coef.end(); ++cit)
    //
    res += cit->second(dist) * basis_function(cit->first, state);

  return res;
}

// axis index to name
//
std::string i2n (int i)
{
  const char funame [] = "i2n: ";

  switch(i) {
    //
  case 0:
    //
    return "X";
    //
  case 1:
    //
    return "Y";
    //
  case 2:
    //
    return "Z";
    //
  default:
    //
    std::cerr << funame << "axis index out of range: " << i << "\n";

    throw Error::Range();
  }
}

ConstSharedPointer<Potential::YG::Base> Potential::YG::pot;
    
void Potential::YG::init(std::istream& from)
{
  const char funame [] = "Potential::YG::init: ";

  int         itemp;
  double      dtemp;
  std::string stemp;

  std::string token, comment, name;
  
  IO::LineInput lin(from);

  if(!(lin >> name)) {
    //
    std::cerr << funame << "cannot read potential name\n";

    throw Error::Input();
  }

  if(name == "Cs+Atom" || name == "C1+Atom") {
    //
    int axis = -1;
    
    Key axis_key("Axis");

    while(from >> token) {
      //
      if(token == IO::end_key()) {
	//
	break;
      }
      // reflection symmetry axis
      //
      else if(token == axis_key) {
	//
	lin.read_line(from);

	if(!(lin >> stemp)) {
	  //
	  std::cerr << funame << token << ": corrupted\n";

	  throw Error::Input();
	}

	if(stemp == "X" || stemp == "x") {
	  //
	  axis = 0;
	}
	else if(stemp == "Y" || stemp == "y") {
	  //
	  axis = 1;
	}
	else if(stemp == "Z" || stemp == "z") {
	  //
	  axis = 2;
	}
	else {
	  //
	  std::cerr << funame << token << ": unknown axis symbol: " << stemp << ": available symbols: X, Y, Z\n";

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

    if(axis < 0) {
      //
      std::cerr << funame << "axis not initialized\n";

      throw Error::Init();
    }

    if(name == "Cs+Atom") {
      //
      pot.init(new Cs_Atom(axis));
    }

    if(name == "C1+Atom") {
      //
      pot.init(new C1_Atom(axis));
    }
  }  
  // unknown potential name
  //
  else {
    //
    std::cerr << funame << "unknown potential name: " << stemp << "\n";

    throw Error::Init();
  }
}

/*******************************************************************************
 ********************************* Cs + ATOM ***********************************
 *******************************************************************************/
//
Potential::YG::Cs_Atom::Cs_Atom (int z) : _x((z + 1) % 3) 
{
  const char funame [] = "Potential::YG::Cs_Atom:Cs_Atom: ";

  double      dtemp;
  int         itemp;
  std::string stemp;

  if(z < 0 || z > 2) {
    //
    std::cerr << funame << "symmetry axis index out of range: " << z << "\n";

    throw Error::Range();
  }
  
  if(!Structure::fragment(0).symmetry_group) {
    //
    std::cerr << funame << "the anchor fragment symmetry group is not initialized\n";

    throw Error::Range();
  }

  stemp = Structure::fragment(0).symmetry_group->name();
  
  if(stemp != "Cs") {
    //
    std::cerr << funame << "the anchor fragment symmetry is not Cs: " << stemp << "\n";

    throw Error::Range();
  }

  if(!Structure::fragment(0).is_symmetric(Symmetry::reflection(z))) {
    //
    std::cerr << funame << "the anchor fragment reflection symmetry axis is not " << i2n(z) << "\n";

    throw Error::Range();
  }

  if(Structure::type(1) != Molecule::MONOATOMIC) {
    //
    std::cerr << funame << "the mobile fragment is not monoatomic\n";

    throw Error::Range();
  }
}

// Cs + Atom basis function
//
double Potential::YG::Cs_Atom::basis_function (const std::vector<int>& index, const Array<double>& state) const
{
  const char funame [] = "Potential::YG::Cs_Atom:basis_function: ";

  double dtemp;
  int    itemp;

  if(index.size() != 2 || index[0] < 0 || index[1] < 0) {
    //
    std::cerr << funame << "index out of range: " << index.size();

    for(int i = 0; i < index.size(); ++i)
      //
      std::cerr << ", " << index[i];

    std::cerr << "\n";

    throw Error::Range();
  }
  
  if(state.size() != 3) {
    //
    std::cerr << funame << "state vector dimension out of range: " << state.size() << "\n";

    throw Error::Logic();
  }

  double res = 1.;

  double lval = 0.;

  itemp = 1;
  
  for(int i = 0; i < 2; ++i) 
    //
    if(!index[i]) {
      //
      itemp = 1 - i;

      if(!index[itemp])
	//
	return 1.;
      
      dtemp = state[(_x + itemp) % 3];

      if(dtemp == 0.)
	//
	return 0.;
      
      lval = std::log(std::fabs(dtemp)) * (double)index[itemp];

      if(lval < -Limits::exp_pow_max())
	//
	return 0.;
      
      res = std::exp(lval);

      if(index[itemp] % 2 && dtemp < 0.) {
	//
	return -res;
      }
      else
	//
	return res;
    }

  lval = 0.;

  itemp = 1;
  
  for(int i = 0; i < 2; ++i) {
    //
    dtemp = state[(_x + i) % 3];

    if(dtemp == 0.)
      //
      return 0.;
    
    lval += std::log(dtemp * dtemp * double(index[0] + index[1]) / (double)index[i]) / 2. * (double)index[i];

    if(index[i] % 2 && dtemp < 0.)
      //
      itemp = -itemp;
  }
  
  if(lval < -Limits::exp_pow_max())
    //
    return 0.;
      
  return std::exp(lval) * double(itemp);
  //
}// Cs_Atom basis function

/*******************************************************************************
 ********************************* C1 + ATOM ***********************************
 *******************************************************************************/
//
Potential::YG::C1_Atom::C1_Atom (int z) : _x((z + 1) % 3) 
{
  const char funame [] = "Potential::YG::Cs_Atom:Cs_Atom: ";

  double      dtemp;
  int         itemp;
  std::string stemp;

  if(z < 0 || z > 2) {
    //
    std::cerr << funame << "symmetry axis index out of range: " << z << "\n";

    throw Error::Range();
  }
  
  if(!Structure::fragment(0).symmetry_group) {
    //
    std::cerr << funame << "the anchor fragment symmetry group is not initialized\n";

    throw Error::Range();
  }

  stemp = Structure::fragment(0).symmetry_group->name();
  
  if(stemp != "C1")
    //
    std::cerr << funame << "WARNING: the anchor fragment symmetry is not C1: " << stemp << "\n";

  if(Structure::type(1) != Molecule::MONOATOMIC) {
    //
    std::cerr << funame << "the mobile fragment is not monoatomic\n";

    throw Error::Range();
  }
}

double Potential::YG::C1_Atom::basis_function (const std::vector<int>& index, const Array<double>& state) const
{
  const char funame [] = "Potential::YG::Cs_Atom:basis_function: ";

  double dtemp;
  int    itemp;

  if(index.size() != 3 || index[0] < 0 || index[1] < 0 || index[2] < 0 || index[2] > 1) {
    //
    std::cerr << funame << "index out of range: " << index.size();

    for(int i = 0; i < index.size(); ++i)
      //
      std::cerr << ", " << index[i];

    std::cerr << "\n";

    throw Error::Range();
  }
  
  if(state.size() != 3) {
    //
    std::cerr << funame << "state vector dimension out of range: " << state.size() << "\n";

    throw Error::Logic();
  }

  double res = 1.;

  if(index[2])
    //
    res = state[(_x + 2) % 3];

  double lval;
  
  for(int i = 0; i < 2; ++i) 
    //
    if(!index[i]) {
      //
      itemp = 1 - i;

      if(!index[itemp])
	//
	return res;
      
      dtemp = state[(_x + itemp) % 3];

      if(dtemp == 0.)
	//
	return 0.;
      
      lval = std::log(std::fabs(dtemp)) * (double)index[itemp];
      
      if(lval < -Limits::exp_pow_max())
	//
	return 0.;
      
      res *= std::exp(lval);

      if(index[itemp] % 2 && dtemp < 0.) {
	//
	return -res;
      }
      else
	//
	return res;
    }

  itemp = 1;

  lval = 0.;
  
  for(int i = 0; i < 2; ++i) {
    //
    dtemp = state[(_x + i) % 3];

    if(dtemp == 0.)
      //
      return 0.;
    
    lval += std::log(dtemp * dtemp * double(index[0] + index[1]) / (double)index[i]) / 2. * (double)index[i];

    if(index[i] % 2 && dtemp < 0.)
      //
      itemp = -itemp;
  }
  
  if(lval < -Limits::exp_pow_max())
    //
    return 0.;
      
  res *= std::exp(lval);

  switch(itemp) {
    //
  case 1:
    //
    return res;
    //
  case -1:
    //
    return -res;
    //
  default:
    //
    std::cerr << funame << "wrong case: " << itemp << "\n";

    throw Error::Range();
  }//
  //
}// C1_Atom basis function

/*************************************************************************
 ******************************** Charge-Linear **************************
 *************************************************************************/

double Potential::ChargeLinear::_pot (double r, double x, _mode_t mode) const 
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

double Potential::ChargeLinear::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const 
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

Potential::ChargeLinear::ChargeLinear (std::istream& from)
  //
  : _charge(1.), _dipole(0.), _quadrupole(0.), _isotropic_polarizability(0.), _anisotropic_polarizability(0.)
{
  static const char funame [] = "Potential::ChargeLinear::ChargeLinear (std::istream&): ";

  if(Structure::fragment(1).type() != Molecule::LINEAR) {
    //
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

  // read cycle
  //
  while(from >> token) {
    //
    if(token == IO::end_key()) {
      //
      std::getline(from, comment);
      
      break;
    }
    // charge
    //
    else if(token == chr_key) {
      //
      if(!(from >> _charge)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
       std::getline(from, comment);
    }
    // dipole
    //
    else if(token == dip_key) {
      //
      if(!(from >> _dipole)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // quadrupole
    //
    else if(token == quad_key) {
      //
      if(!(from >> _quadrupole)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // isotropic polarizability
    //
    else if(token == isot_key) {
      //
      if(!(from >> _isotropic_polarizability)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // anisotropic polarizability
    //
    else if(token == anis_key) {
      //
      if(!(from >> _anisotropic_polarizability)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Input();
      }
      std::getline(from, comment);
    }
    // unknown key
    //
    else {
      //
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
  //
  _dipole                     *= _charge;
  _quadrupole                 *= _charge;
  _isotropic_polarizability   *= _charge * _charge;
  _anisotropic_polarizability *= _charge * _charge;
}

/*************************************************************************
 ***************************** Charge-Nonlinear **************************
 *************************************************************************/

double Potential::ChargeNonlinear::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const 
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
    //
    return ener_val;

  D3::Vector& force = *torque;
  torque += 1;

  torque[0] = 0.;

  D3::Vector vtemp;
  
  // dipole
  //
  force = r3 * _dipole - (3. * prd * r5) * mfr;

  D3::vprod(_dipole, mfr, vtemp);
  vtemp *= r3;
  torque[1] = vtemp;

  //quadrupole
  //
  force += (2.5 * prq * r7) * mfr - r5 * mfq;

  D3::vprod(mfr, mfq, vtemp);
  vtemp *= r5;
  torque[1] += vtemp;

  // polarizability
  //
  force += r6 * mfp - (3. * prp * r8) * mfr;

  D3::vprod(mfp, mfr, vtemp);
  vtemp *= r6;
  torque[1] += vtemp;

  // convert force to laboratory frame
  //
  dc.mf2lf(1, force, vtemp);
  force = vtemp;

  return ener_val;
}

Potential::ChargeNonlinear::ChargeNonlinear (std::istream& from) : _charge(1.)
{
  static const char funame [] = "Potential::ChargeNonlinear::ChargeNonlinear: ";

  if(Structure::fragment(1).type() != Molecule::NONLINEAR || Structure::fragment(1).top() == Molecule::SPHERICAL) {
    //
    std::cerr << funame << "second fragment should be nonlinear for and non-spherical this potential\n";
    
    throw Error::Logic();
  }

  IO::Marker funame_marker(funame);

  double              dtemp;
  std::vector<double> vtemp;

  KeyGroup ChargeNonlinearPotential;

  Key   chr_key("Charge[au]"          );
  Key   dip_key("DipoleMoment[au]"    );
  Key  quad_key("QuadrupoleMoment[au]");
  Key polar_key("Polarizability[au]"  );

  for(int i = 0; i < 3; ++i)
    //
    _dipole[i] = 0.;
  
  for(int i = 0; i < 3; ++i)
    //
    for(int j = i; j < 3; ++j) {
      //
      _quadrupole(i, j) = 0.;
      _quadrupole(j, i) = 0.;
      
      _polarizability(i, j) = 0.;
      _polarizability(j, i) = 0.;
    }

  std::string token, line, comment;

  // read cycle
  //
  while(from >> token) {
    //
    if(token == IO::end_key()) {
      //
      std::getline(from, comment);
      
      break;
    }
    // charge
    //
    else if(token == chr_key) {
      //
      if(!(from >> _charge)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	throw Error::Init();
      }
      std::getline(from, comment);
    }
    // dipole
    //
    else if(token == dip_key) {
      //
      std::getline(from, line);
      
      std::istringstream iss(line);
      
      vtemp.clear();
      
      while(iss >> dtemp)
	//
	vtemp.push_back(dtemp);

      switch(Structure::top(1)) {
	//
      case Molecule::ASYM:
	//
	if(vtemp.size() == 3) {
	  //
	  for(int i = 0; i < 3; ++i)
	    //
	    _dipole[i] = vtemp[i];
	}
	else {
	  //
	  std::cerr << funame << "dipole moment of non-spherical fragment should have three components\n";
	  
	  throw Error::Input();
	}
	
	break;
	//
      case Molecule::PROLATE:
	//
	if(vtemp.size() != 1) {
	  //
	  std::cerr << funame << "there should be only one component for the dipole moment along the symmetry axis\n";

	  throw Error::Init();
	}
	_dipole[0] = vtemp[0];

	break;
	//
      case Molecule::OBLATE:
	//
	if(vtemp.size() != 1) {
	  //
	  std::cerr << funame << "there should be only one component for the dipole moment along the symmetry axis\n";

	  throw Error::Init();
	}
	_dipole[2] = vtemp[0];
      }
    }
    // quadrupole
    //
    else if(token == quad_key) {
      //
      std::getline(from, line);
      
      std::istringstream iss(line);
      
      vtemp.clear();
      while(iss >> dtemp)
	//
	vtemp.push_back(dtemp);
      
      int k = 0;
      switch(Structure::top(1)) {
	//
      case Molecule::ASYM:
	//
	if(vtemp.size() == 5) {
	  //
	for(int i = 0; i < 2; ++i)
	  //
	  for(int j = i; j < 3; ++j) {
	    //
	    _quadrupole(i, j) = vtemp[k++];
	    
	    if(i != j)
	      //
	      _quadrupole(j, i) = _quadrupole(i, j);
	  }
	}
	else {
	  //
	  std::cerr << funame << "quadrupole moment of non-spherical top should have five components: qxx, qxy, qxz, qyy, qyz\n";
	  
	  throw Error::Input();
	}
	
	break;
	//
      case Molecule::PROLATE:
	//
	if(vtemp.size() != 1) {
	  //
	  std::cerr << funame << "quadrupole moment of symmetric top should have one component along the symmetry axis\n";
	  
	  throw Error::Input();
	}
	_quadrupole(0, 0) = vtemp[0];
	_quadrupole(1, 1) = -_quadrupole(0, 0) / 2.;
	_quadrupole(2, 2) = -_quadrupole(0, 0) / 2.;

	break;
	//
      case Molecule::OBLATE:
	//
	if(vtemp.size() != 1) {
	  //
	  std::cerr << funame << "quadrupole moment of symmetric top should have one component along the symmetry axis\n";
	  throw Error::Input();
	}
	_quadrupole(2, 2) = vtemp[0];
	_quadrupole(0, 0) = -_quadrupole(2, 2) / 2.;
	_quadrupole(1, 1) = -_quadrupole(2, 2) / 2.;
	
	break;
	//
      default:
	//
	std::cerr << funame << token << "should not be here\n";
	throw Error::Logic();
      }
    }
    //  polarizability
    //
    else if(token == polar_key) {
      //
      std::getline(from, line);
      
      std::istringstream iss(line);
      vtemp.clear();
      
      while(iss >> dtemp)
	//
	vtemp.push_back(dtemp);
      
      int k = 0;

      switch(Structure::top(1)) {
	//
      case Molecule::ASYM:
	//
	if(vtemp.size() != 6) {
	  //
	  std::cerr << funame << "polarizability of non-spherical top should have six components: pxx, pxy, pxz, pyy, pyz, pzz\n";
	  
	  throw Error::Input();
	}
	for(int i = 0; i < 3; ++i)
	  //
	  for(int j = i; j < 3; ++j) {
	    //
	    _polarizability(i, j) = vtemp[k++];
	    
	    if(i != j)
	      //
	      _polarizability(j, i) = _polarizability(i, j);
	  }

	break;
	//
      case Molecule::PROLATE:
	//
	if(vtemp.size() != 2) {
	  //
	  std::cerr << funame << "polarizability of prolate symmetric top should have two components: p_xx, p_yy\n";
	  throw Error::Input();
	}
	for(int i = 0; i < 2; ++i)
	  //
	  _polarizability(i, i) = vtemp[i];

	_polarizability(2, 2) = _polarizability(1, 1);
	
	break;
	//
      case Molecule::OBLATE:
	//
	if(vtemp.size() != 2) {
	  //
	  std::cerr << funame << "polarizability of oblate symmetric top should have two components: p_yy, p_zz\n";
	  throw Error::Input();
	}
	
	for(int i = 1; i < 3; ++i)
	  //
	  _polarizability(i, i) = vtemp[i - 1];
	
	_polarizability(0, 0) = _polarizability(1, 1);
	
	break;
	//
      default:
	//
	std::cerr << funame << token << ": should not be here\n";
	throw Error::Logic();
      }
    }
    // unknown key
    //
    else {
      //
      std::cerr << funame << "unknown keyword: " << token << "\n";
      
      Key::show_all(std::cerr);
      std::cerr << "\n";
      
      throw Error::Init();
    }
  }// read cycle

  // print out
  //
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
  //
  for(int i = 0; i < 3; ++i)
    //
    _dipole[i] *= _charge;

  for(int i = 0; i < 3; ++i)
    //
    for(int j = 0; j < 3; ++j) {
      //
      _quadrupole    (i, j) *= _charge;
      
      _polarizability(i, j) *= _charge * _charge;
    }
}

  /*************************************************************************
   ******************************** Dipole-Dipole **************************
   *************************************************************************/

double Potential::DipoleDipole::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const 
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
    //
    for(int frag = 0; frag < 2; ++frag)
      if(Structure::fragment(frag).type() == Molecule::NONLINEAR) {
	dc.lf2mf(frag, torque[frag], vtemp);
	torque[frag] = vtemp;
      }

    return ener_val;
  }

Potential::DipoleDipole::DipoleDipole (std::istream& from) 
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

double Potential::Multipole::operator() (const Dynamic::Coordinates& dc, D3::Vector* torque) const 
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

Potential::Multipole::Multipole (std::istream& from)  : _dispersion(0.)
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
