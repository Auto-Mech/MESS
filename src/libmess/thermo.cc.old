#include <list>

#include "thermo.hh"
#include "graph_omp.hh"
#include "key.hh"

// maximal number of importance samplings
//
int     Thermo::Species::count_max        = 1000;

// number of importance samplings to thermolize
//
int     Thermo::Species::count_min        = 10;

// numerical derivative step (in bohr)
//
double  Thermo::Species::numd_step        = 0.01;

// importance sampling step (in bohr)
//
double  Thermo::Species::samp_step        = 0.1;

// low  frequency threshold in qfactor
//
double  Thermo::Species::low_freq_thresh  = 0.1;

// high frequency threshold in qfactor
//
double  Thermo::Species::high_freq_thresh =  5.;

// deep tunneling threshold
//
double  Thermo::Species::deep_tun_thresh =  -0.8;


/********************************************************************************************
 ******************************* GROWING SYMMETRIC TENSOR INDEX ***************************** 
 ********************************************************************************************/

void Thermo::SymIndex::operator++ () 
{
  int itemp;

  for(int i = 0; i < size(); ++i) {
    //
    itemp = ++(*this)[i];
    //
    if(itemp < _range) {
      //
      for(int j = 0; j < i; ++j) {
	//
	(*this)[j] = itemp;
      }
      //
      return;
    }
  }

  itemp = size() + 1;
  //
  clear();
  //
  resize(itemp, 0);
}

Thermo::SymIndex::operator std::multiset<int> () const 
{
  std::multiset<int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    res.insert(*it);
  }

  return res;
}

Thermo::SymIndex::operator std::map<int, int> () const 
{
  std::map<int, int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    ++res[*it];
  }

  return res;
}

/*********************************************************************************************
 ******************************* FIXED SIZE GENERIC TENSOR INDEX ***************************** 
 *********************************************************************************************/

void Thermo::GenIndex::operator++ () 
{
  const char funame [] = "Thermo::GenIndex::operator++: ";

  if(_fin) {
    ErrOut err_out;
    err_out << funame << "already finished";
  }
  
  int itemp;

  for(int i = 0; i < size(); ++i) {
    //
    itemp = ++(*this)[i];
    //
    if(itemp < _range) {
      //
      for(int j = 0; j < i; ++j) {
	//
	(*this)[j] = 0;
      }
      //
      return;
    }
  }
  
  _fin = true;
}

Thermo::GenIndex::operator std::multiset<int> () const 
{
  std::multiset<int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    res.insert(*it);
  }

  return res;
}

Thermo::GenIndex::operator std::map<int, int> () const 
{
  std::map<int, int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    ++res[*it];
  }

  return res;
}

/***********************************************************************************************
***************************** NUMERICAL DIFFERENTIATION MAPPER *********************************
************************************************************************************************/

// converts derivative signature into the displaced configurations mapping
//
Thermo::NumDer::cmap_t Thermo::NumDer::operator() (const der_t& der) const
{
  const char funame [] = "Thermo::NumDer::operator(): ";
  
  _isinit();

  cmap_t cmap;

  if(!der.size()) {
    //
    cmap[std::vector<int>(_size)] = 1;
    //
    return cmap;
  }

  if(der.begin()->first < 0 || der.begin()->first >= _size) {
    //
    ErrOut err_out;
    //
    err_out << funame << "wrong derivative index: " << der.begin()->first;
  }
    
  if(der.begin()->second <= 0) {
    //
    ErrOut err_out;
    //
    err_out << funame << "wrong derivative order: " << der.begin()->second;
  }

  der_t new_der = der;

  if(der.begin()->second <= 2) {
    //
    new_der.erase(new_der.begin());
  }
  else {
    //
    new_der.begin()->second -= 2;
  }

  cmap_t new_cmap = (*this)(new_der);

  for(cmap_t::const_iterator cmit = new_cmap.begin(); cmit != new_cmap.end(); ++cmit) {
    //
    std::vector<int> conf = cmit->first;
	
    // first derivative
    //
    if(der.begin()->second == 1) {
      //
      conf[der.begin()->first] += 1;
      //
      cmap[conf] += cmit->second;

      conf[der.begin()->first] -= 2;
      //
      cmap[conf] -= cmit->second;
    }
    // second derivative
    //
    else {
      //
      cmap[conf] -= 2 * cmit->second;

      conf[der.begin()->first] += 1;
      //
      cmap[conf] += cmit->second;
	
      conf[der.begin()->first] -= 2;
      //
      cmap[conf] += cmit->second;
    }
  }

  // clean configuration map
  //
  cmap_t res;
  //
  for(cmap_t::const_iterator cmit = cmap.begin(); cmit != cmap.end(); ++cmit) {
    //
    if(cmit->second) {
      //
      res.insert(*cmit);
    }
  }
  return res;
}

/************************************************************************************
 ********************************* POTENTIAL ENERGY SURFACE *************************
 ************************************************************************************/

double Thermo::Harding::operator() (const  pos_t& pos) const 
{
  const char funame [] = "Thermo::Harding::operator(): ";

  if(size() * 3 != pos.size()) {
    ErrOut err_out;
    err_out << funame << "number of atoms and number of coordinates mismatch: " << size() << ", " << pos.size();
  }

  // coordinates in angstrom
  //
  Array<double> apos((int)pos.size());

  for(int i = 0; i < pos.size(); ++i)
    //
    apos[i] = pos[i] / Phys_const::angstrom;
    
  double res;
  //
  harding_pot_(apos, res);
  //
  return res * Phys_const::kcal; 
}

Thermo::PotWrap Thermo::new_pot(std::istream& from)
{
  const char funame [] = "Thermo::new_pot: ";

  KeyGroup NewPotGroup;

  Key hard_key("Harding");
  
  std::string token, comment;

  if(!(from >> token)) {
    ErrOut err_out;
    err_out << funame << "corrupted";
  }

  std::getline(from, comment);

  // unknown keyword
  //
  if(hard_key != token) {
    ErrOut err_out;
    err_out << funame << "unknown keyword: " << token << "\n"
	    << Key::show_all();
  }

  // harding potential
  //
  return PotWrap(new Harding(from));
}

/***********************************************************************************************************
 ******************************* CHEMICAL SPECIES PARTITION FUNCTION CORRECTION  *************************** 
 ***********************************************************************************************************/

void Thermo::Species::init (const std::vector<Atom>& at, PotWrap pw, ConWrap cw)
{
  const char funame [] = "Thermo::Species::init: ";

  IO::Marker funame_marker(funame);

  double dtemp;
  int    itemp;
  bool   btemp;
  
  if(at.size() < 2) {
    ErrOut err_out;
    err_out << funame << "not enough atoms: " << at.size();
  }
  //
  _atom = at;
  
  if(!pw) {
    ErrOut err_out;
    err_out << funame << "potential not initialized";
  }
  //
  _pot = pw;

  if(_atom.size() != _pot.size()) {
    ErrOut err_out;
    err_out << funame << "number of atoms in the potential description and in the species mismatch: "
	    << _pot.size() << ", " << _atom.size();
  }

  // numerical derivative configurations signatures 
  //
  for(SymIndex sin(_atom.size() * 3); sin.size() <= Graph::potex_max; ++sin) {
    //
    NumDer::cmap_t cmap = NumDer(_atom.size() * 3)(sin);

    for(NumDer::cmap_t::const_iterator cmit = cmap.begin(); cmit != cmap.end(); ++cmit) {
      //
      _conf_pool.insert(cmit->first);
    }
  }

  // individual mass square root
  //
  std::vector<double> imass_sqrt(_atom.size());

  // cumulative mass square root
  //
  std::vector<double> cmass_sqrt(_atom.size());

  dtemp = 0.;
  //
  for(int i = 0; i < _atom.size(); ++i) {
    //
    dtemp += _atom[i].mass();

    imass_sqrt[i] = std::sqrt(_atom[i].mass());

    cmass_sqrt[i] = std::sqrt(dtemp);
  }

  // excluded cm motion basis
  //
  _cm_orth.resize(_atom.size(), _atom.size() - 1);

  _cm_orth = 0.;
  //
  for(int k = 0; k < _cm_orth.size2(); ++k) {
    //
    dtemp = imass_sqrt[k + 1] / cmass_sqrt[k + 1] / cmass_sqrt[k];

    for(int i = 0; i <= k; ++i) {
      //      
      _cm_orth(i, k) = dtemp;
    }
   
    _cm_orth(k + 1, k) = - cmass_sqrt[k] / cmass_sqrt[k + 1] / imass_sqrt[k + 1];
  }

  // starting configuration
  //
  _start_pos.resize(_atom.size() * 3);

  for(int a = 0; a < _atom.size(); ++a) {
    //
    for(int i = 0; i < 3; ++i) {
      //
      _start_pos[a * 3 + i] = _atom[a][i];
    }
  }

  // reference frequencies
  //
  std::map<std::vector<int>, double> conf_value;

  conf_value[std::vector<int>(_atom.size() * 3)] = _pot(_start_pos);
  
  _potex_t potex = _set_potex(_start_pos, conf_value, 2, _ref_freq);

  IO::log << "\n" <<  IO::log_offset << "reference frequencies, 1/cm:";
  //
  for(int i = 0; i < _ref_freq.size(); ++i) {
    //
    IO::log << "   " << _ref_freq[i] / Phys_const::incm;
  }

  IO::log << "\n\n";

  IO::log << IO::log_offset << "deviations from minimum / sadle point, 1/cm:";

  for(int i = 0; i < _ref_freq.size(); ++i) {
    //
    std::multiset<int> sin;

    sin.insert(i);

    _potex_t::const_iterator pit = potex.find(sin);

    if(pit == potex.end()) {
      //
      ErrOut err_out;
      
      err_out << funame << i << "-th gradient not in the potential expansion";
    }

    if(_ref_freq[i] > Phys_const::incm || _ref_freq[i] < -Phys_const::incm) {
      //
      dtemp = pit->second * pit->second / _ref_freq[i] / _ref_freq[i] / 2.;

      if(_ref_freq[i] > 0.)
	//
	dtemp = -dtemp;

      IO::log << "   " << dtemp / Phys_const::incm; 
    }
    else
      //
      IO::log << "   ***"; 
  } 

  IO::log << "\n\n";
}

std::vector<double> Thermo::Species::_frequencies (const pos_t& init_pos) const
{
  const char funame [] = "Thermo::Species::_frequencies: ";

  double dtemp;
  int    itemp;
  bool   btemp;

  Lapack::SymmetricMatrix force_const(_atom.size() * 3);
  //
  for(int i = 0; i < force_const.size(); ++i) {
    //
    for(int j = i; j < force_const.size(); ++j) {
      //
      pos_t pos = init_pos;

      if(i != j) {
	//
	pos[i] += numd_step;
	//
	pos[j] += numd_step;
	//
	dtemp   = _pot(pos);

	pos[i] -= 2. * numd_step;
	//
	dtemp  -= _pot(pos);

	pos[j] -= 2. * numd_step;
	//
	dtemp  += _pot(pos);

	pos[i] += 2. * numd_step;
	//
	dtemp  -= _pot(pos);

	dtemp /= 4. * numd_step * numd_step;
	
	force_const(i, j) = dtemp;
      }
      else {
	//
	dtemp = -2. * _pot(pos);

	pos[i] += numd_step;
	//
	dtemp  += _pot(pos);

	pos[i] -= 2. * numd_step;
	//
	dtemp  += _pot(pos);

	dtemp /= numd_step * numd_step;

	force_const(i, i) = dtemp;
      }
    }
  }

  // mass waited force constant
  //
  Lapack::SymmetricMatrix mwfc(_atom.size() * 3 - 3);
  //
  for(int i1 = 0; i1 < mwfc.size(); ++i1) {
    //
    int ai1 = i1 / 3;
    
    int ci  = i1 % 3;
    
    for(int j1 = i1; j1 < mwfc.size(); ++j1) {
      //
      int aj1 = j1 / 3;

      int cj  = j1 % 3;
      
      dtemp = 0.;
      //
      for(int ai2 = 0; ai2 < _atom.size(); ++ai2) {
	//
	for(int aj2 = 0; aj2 < _atom.size(); ++aj2)
	  //
	  dtemp += force_const(3 * ai2 + ci, 3 * aj2 + cj) * _cm_orth(ai2, ai1) * _cm_orth(aj2, aj1);
      }

      mwfc(i1, j1) = dtemp;
    }
  }
  
  Lapack::Vector eval = mwfc.eigenvalues();

  std::vector<double> freq(eval.size());
  //
  for(int i = 0; i < eval.size(); ++i) {
    //
    dtemp = eval[i];
   
    if(dtemp >= 0.) {
      //
      dtemp = std::sqrt(dtemp);
    }
    else {
      //
      dtemp = -std::sqrt(-dtemp);
    }

    freq[i] = dtemp;
  }

  return freq;
}

void Thermo::Species::_print_geom (const pos_t& pos) const
{
  const char funame [] = "Thermo::Species::_print_geom: ";

  if(pos.size() != _atom.size() * 3) {
    //
    ErrOut err_out;

    err_out << funame << "configuration dimension and number of atom mismatch: " << pos.size() << ", " << _atom.size();
  }
  
  IO::log << IO::log_offset << "geometry, angstrom:\n";
  //
  for(int a = 0; a < _atom.size(); ++a) {
    //
    IO::log << IO::log_offset << std::setw(2) << _atom[a].name();

    for(int i = 0; i < 3; ++i)
      //
      IO::log << std::setw(15) << pos[a * 3 + i] / Phys_const::angstrom;

    IO::log << "\n";
  }

  IO::log << "\n";
}

// potential exapansion in normal modes coordinates with center-of-mass motion excluded, around centroid position
//
Thermo::Species::_potex_t Thermo::Species::_set_potex (const pos_t&                              cent_pos, 
						       const std::map<std::vector<int>, double>& conf_value, 
						       int                                       potex_max, 
						       std::vector<double>&                      freq
						       ) const
{
  const char funame [] = "Thermo::Species::_set_potex: ";

  double dtemp;
  int    itemp;

  if(cent_pos.size() != _atom.size() * 3) {
    //
    ErrOut err_out;

    err_out << funame << "centroid position dimension and number of atoms mismatch: " << cent_pos.size() << ", " << _atom.size();
  }

  if(potex_max < 2) {
    //
    ErrOut err_out;
    
    err_out << funame << "expansion rank should be no less than two";
  }

  // numerical potential derivatives
  //
  _potex_t pot_der;

  for(SymIndex sin(_atom.size() * 3); sin.size() <= potex_max; ++sin) {
    //
    // derivative signature
    //
    std::map<int, int> dmap = sin;

    // configuration signatures map for numeric differentiation
    //
    NumDer::cmap_t cmap = NumDer(_atom.size() * 3)(dmap);

    dtemp = 0.;
    //
    for(NumDer::cmap_t::const_iterator cmit = cmap.begin(); cmit != cmap.end(); ++cmit) {
      //
      std::map<std::vector<int>, double>::const_iterator cvit = conf_value.find(cmit->first);

      if(cvit == conf_value.end()) {
	//
	pos_t pos = cent_pos;

	if(pos.size() != cmit->first.size()) {
	  //
	  ErrOut err_out;
	 
	  err_out << funame << "number of coodinates and configuration signature dimensions mismatch: " << pos.size() << ", " << cmit->first.size();
	}

	for(int i = 0; i < pos.size(); ++i) {
	  //
	  pos[i] += (double)cmit->first[i] * numd_step;
	}
    
	dtemp += (double)cmit->second * _pot(pos);
      }
      else
	//
	dtemp += (double)cmit->second * cvit->second;
    }

    // correction for odd derivatives
    //
    for(std::map<int, int>::const_iterator dmit = dmap.begin(); dmit != dmap.end(); ++dmit)
      //
      if(dmit->second % 2)
	//
	dtemp /= 2.;

    if(sin.size())
      //
      dtemp /= std::pow(numd_step, (double)sin.size());

    pot_der[sin] = dtemp;
  }

  // mass waited force constant
  //
  Lapack::SymmetricMatrix mwfc(_atom.size() * 3 - 3);

  for(int i1 = 0; i1 < mwfc.size(); ++i1) {
    //
    int ai1 = i1 / 3;
    
    int ci  = i1 % 3;
    
    for(int j1 = i1; j1 < mwfc.size(); ++j1) {
      //
      int aj1 = j1 / 3;

      int cj  = j1 % 3;
      
      dtemp = 0.;
      //
      for(int ai2 = 0; ai2 < _atom.size(); ++ai2) {
	//
	for(int aj2 = 0; aj2 < _atom.size(); ++aj2) {
	  //
	  std::multiset<int> gin;
	  
	  gin.insert(3 * ai2 + ci);
	  
	  gin.insert(3 * aj2 + cj);

	  _potex_t::const_iterator pit = pot_der.find(gin);
	  //
	  if(pit == pot_der.end()) {
	    //
	    ErrOut err_out;

	    err_out << funame << "term is not in the pot_der expansion";
	  }
	  
	  dtemp += pit->second * _cm_orth(ai2, ai1) * _cm_orth(aj2, aj1);
	}
      }
      
      mwfc(i1, j1) = dtemp;
    }
  }

  Lapack::Matrix evec;
  
  Lapack::Vector eval = mwfc.eigenvalues(&evec);

  freq.resize(eval.size());

  for(int i = 0; i < eval.size(); ++i) {
    //
    dtemp = eval[i];
   
    if(dtemp >= 0.) {
      //
      freq[i] = std::sqrt(dtemp);
    }
    else {
      //
      freq[i] = -std::sqrt(-dtemp);
    }
  }

  // transformation to normal mode coordinates
  //
  Lapack::Matrix tran(3 * _atom.size(), 3 * _atom.size() - 3);
  
  for(int i = 0; i < tran.size1() ; ++i) {
    //
    int ai = i / 3;

    int ci = i % 3;

    for(int j = 0; j < tran.size2(); ++j) {
      //
      dtemp = 0.;

      for(int ak = 0; ak < _cm_orth.size2(); ++ak) {
	//
	int k = 3 * ak + ci;
      
	dtemp += _cm_orth(ai, ak) * evec(k, j);
      }
    
      tran(i, j) = dtemp;
    }
  }
      
  // potential expansion in normal mode coordinates
  //
  _potex_t potex;

  for(SymIndex sin(_atom.size() * 3 - 3); sin.size() <= potex_max; ++sin) {
    //
    potex[sin] = 0.;

    _potex_t::iterator pexit = potex.find(sin);
    
    for(GenIndex gin(sin.size(), _atom.size() * 3); !gin.fin(); ++gin) {
      //
      _potex_t::const_iterator pit = pot_der.find(gin);

      if(pit == pot_der.end()) {
	//
	ErrOut err_out;

	err_out << funame << "potential expansion term signature not found in pot_der map";
      }

      dtemp = pit->second;
      //
      for(int i = 0; i < sin.size(); ++i) {
	//
	dtemp *= tran(gin(i), sin(i));
      }

      pexit->second += dtemp;
      //
      //
    } // old potential expansion cycle
    //
    //
  } //  new potential expansion cycle

  return potex;
}

Thermo::Species::_corr_t Thermo::Species::_correction (const pos_t& cent_pos, double temperature, double& qfac) const
{
  const char funame [] = "Thermo::Species::_correction: ";

  IO::Marker funame_marker(funame);
  
  double dtemp;
  int    itemp;

  if(cent_pos.size() != _atom.size() * 3) {
    //
    ErrOut err_out;

    err_out << funame << "centroid position dimension and number of atoms mismatch: " << cent_pos.size() << ", " << _atom.size();
  }

  if(temperature <= 0.) {
    //
    ErrOut err_out;

    err_out << funame << "negative temperature[K]: " << temperature / Phys_const::kelv;
  }
    
  // displaced configurations potential values for numerical derivatives
  //
  std::map<std::vector<int>, double> conf_value;
  //
  for(std::set<std::vector<int> >::const_iterator cit = _conf_pool.begin(); cit != _conf_pool.end(); ++cit) {
    //
    pos_t pos = cent_pos;

    for(int i = 0; i < pos.size(); ++i) {
      //
      pos[i] += (double)(*cit)[i] * numd_step;
    }
    
    conf_value[*cit] = _pot(pos);
  }

  // local frequencies
  //
  std::vector<double> freq;

  // potential expansion in normal mode coordinates
  //
  _potex_t potex = _set_potex(cent_pos, conf_value, Graph::potex_max, freq);

  for(int i = 0; i < freq.size(); ++i) {
    //
    if(freq[i] / temperature / 2. / M_PI < deep_tun_thresh) {
      //
      _deep_t x;
	
      x << funame
	<< "deep tunneling regime: frequency[1/cm] = "
	<< dtemp / Phys_const::incm
	<< "   temperature[K] = "
	<< temperature / Phys_const::kelv
	<< "\n";
      
      IO::log << IO::log_offset << x;
      
      throw x;
    }
  }

  IO::log << "\n";

  IO::log << IO::log_offset << "temperature, K = " << temperature / Phys_const::kelv << "\n\n";

  _print_geom(cent_pos);
  
  IO::log << IO::log_offset << "frequencies, 1/cm:";
  //
  for(int i = 0; i < freq.size(); ++i)
    //
    IO::log << "   " << freq[i] / Phys_const::incm;

  IO::log << "\n\n";

  // local harmonic approximation correction
  //
  qfac = _qfactor(freq, temperature);

  // anharmonic correction
  //
  std::map<int, double> res = Graph::Expansion(freq, potex).centroid_correction(temperature);

  IO::log << "\n";

  IO::log << IO::log_offset << "qfactor = " << qfac << "\n\n";

  IO::log << IO::log_offset << "anharmonic correction:\n";

  IO::log << IO::log_offset << std::setw(2) << "BO"<< std::setw(15) << "Value" << "\n";
  
  for(std::map<int, double>::const_iterator git = res.begin(); git != res.end(); ++git) {
    //
    IO::log << IO::log_offset << std::setw(2) << git->first << std::setw(15) << git->second << "\n";
  }

  IO::log << "\n";

  return res;
}

double Thermo::Species::_qfactor (const std::vector<double>& freq, double temperature) const
{
  const char funame [] = "Thermo::Species::_qfactor: ";

  double dtemp;
  int    itemp;

  if(freq.size() != _ref_freq.size()) {
    //
    ErrOut err_out;

    err_out << funame << "wrong number of frequencies: " << freq.size();
  }

  if(temperature <= 0.) {
    //
    ErrOut err_out;

    err_out << funame << "negative temperature[K]: " << temperature / Phys_const::kelv;
  }

  double res = 1.;

  double pow = 0.;

  for(int v = 0; v < 2; ++v) {
    //
    double fac   = 1.;

    double shift = 0.;

    for(int i = 0; i < freq.size(); ++i) {
      //
      if(!v) {
	// 
	dtemp = freq[i] / temperature / 2.;
      }
      else {
	//
	dtemp = _ref_freq[i] / temperature / 2.;
      }

      if(dtemp <= -3.) {
	ErrOut err_out;
	err_out << funame << "deep tunneling regime: frequency[1/cm] =  " << freq[i] / Phys_const::incm 
	      << "   temperature[K] = " << temperature / Phys_const::kelv;
      }

      if(dtemp < -0.5 * low_freq_thresh) {
	//
	fac *= dtemp / std::sin(dtemp);
      }
      else if( dtemp < 0.5 * low_freq_thresh) {
	;
      }
      else if(dtemp < 0.5 * high_freq_thresh) {
	//
	fac *= dtemp / std::sinh(dtemp);
      }
      else {
	//
	fac *= 2. * dtemp;
	//
	shift -= dtemp;
      }
    }

    if(!v) {
      //
      res *= fac;
      //
      pow += shift;
    }
    else {
      //
      res /= fac;
      //
      pow -= shift;
    }
  }

  res *= std::exp(pow);
    
  return res;
}

Thermo::Species::_corr_t Thermo::Species::weight (double temperature, double& harm_corr) const
{
  const char funame [] = "Thermo::Species::weight: ";

  static const double ener_min = -0.2 * Phys_const::kcal;

  IO::Marker funame_marker(funame);
  
  double dtemp;
  int    itemp;
  bool   btemp;

  if(temperature <= 0.) {
    //
    ErrOut err_out;

    err_out << funame << "negative temperature[K]: " << temperature / Phys_const::kelv;
  }

  IO::log << "\n" << IO::log_offset << "temperature, K = " << temperature / Phys_const::kelv << "\n\n";
  
  std::list<std::vector<double> >   corr_pool;
  
  _corr_t tot_corr;

  harm_corr = 0.;

  pos_t pos   = _start_pos;

  double ener = _pot(pos);
  
  // counts
  //
  int count = 0;

  int miss_count = 0;

  int fail_count = 0;
  
  // sampling loop
  //
  while(count < count_max) {
    //
    pos_t new_pos = pos;
    //
    for(int i = 0; i < pos.size(); ++i) {
      //
      new_pos[i] += samp_step * (drand48() - 0.5);
    }

    double new_ener = _pot(new_pos);
    //
    if(new_ener < ener_min || new_ener > ener && std::exp((ener - new_ener) / temperature) < drand48()) {
      //
      ++miss_count;
      //
      continue;
    }

    IO::log << IO::log_offset
	    << "count # = "   << count
	    << "   miss # = " << miss_count 
	    << "   fail # = " << fail_count
	    << "\n" << std::endl;

    miss_count = 0;

    pos =   new_pos;
    //
    ener = new_ener;

    if(count++ < count_min) {
      //
      continue;
    }

    try {
      //
      double qfac;

      _corr_t anharm_corr = _correction(pos, temperature, qfac);

      harm_corr += qfac;

      std::vector<double> corr(anharm_corr.size() + 1);
      //
      corr[0] = qfac;
      
      itemp = 1;
      //
      for(_corr_t::const_iterator cit = anharm_corr.begin(); cit != anharm_corr.end(); ++cit, ++itemp) {
	//
	tot_corr[cit->first]  += qfac * cit->second;
	
	corr[itemp] = cit->second;
      }

      corr_pool.push_back(corr);
    }
    catch(_deep_t x) {
      //
      ++fail_count;

      std::cerr << x;
    }

    IO::log << "\n";
    //
    //
  } // sampling loop
  
  if(count <= count_min + fail_count) {
    //
    ErrOut err_out;

    err_out << funame << "not enough samplings: " << count;
  }

  harm_corr /= (double)(count - count_min - fail_count);

  for(_corr_t::iterator cit = tot_corr.begin(); cit != tot_corr.end(); ++cit) {
    //
    cit->second  /= (double)(count - count_min - fail_count);
  }

  IO::log << IO::log_offset << "corrections:\n";

  IO::log << IO::log_offset << std::setw(12) << "QFactor";
  //
  for(_corr_t::const_iterator cit = tot_corr.begin(); cit != tot_corr.end(); ++cit) {
    //
    IO::log << std::setw(10) << "BO = " << std::setw(2) << cit->first;
  }
  //
  IO::log << "\n";

  for(std::list<std::vector<double> >::const_iterator pit = corr_pool.begin(); pit != corr_pool.end(); ++pit) {
    //
    IO::log << IO::log_offset;

    if(pit->size() - 1 != tot_corr.size()) {
      //
      std::cerr << funame << "WARNING: number of anharmonic corrections and total corrections number mismatch: " 
		<< pit->size() - 1 << ", " << tot_corr.size() << "\n";
    }
    
    for(std::vector<double>::const_iterator cit = pit->begin(); cit != pit->end(); ++cit) {
      //
      IO::log << std::setw(12) << *cit;
    }
    
    IO::log << "\n";
  }

  IO::log << std::endl;

  return tot_corr;
}
