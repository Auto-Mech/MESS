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
double  Thermo::QFactor::low_freq_thresh  = 1.e-3;

// high frequency threshold in qfactor
//
double  Thermo::QFactor::high_freq_thresh =  5.;

// deep tunneling threshold
//
double  Thermo::QFactor::deep_tun_thresh =  -0.8;


/************************************************************************************
 ********************************* POTENTIAL ENERGY SURFACE *************************
 ************************************************************************************/

double Thermo::Harding::operator() (const  pos_t& pos) const 
{
  const char funame [] = "Thermo::Harding::operator(): ";

  if(size() != pos.size()) {
    //
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

void Thermo::CSpec::init (const std::vector<Atom>& mol, ConstSharedPointer<Coord::CartFun> pp)
{
  const char funame [] = "Thermo::CSpec::init: ";

  IO::Marker funame_marker(funame);

  double dtemp;
  int    itemp;
  bool   btemp;
  
  if(mol.size() < 2) {
    //
    ErrOut err_out;

    err_out << funame << "not enough atoms: " << mol.size();
  }
  
  
  if(!pp) {
    //
    ErrOut err_out;

    err_out << funame << "potential not initialized";
  }
  
  _pot = pp;

  if(mol.size() * 3 != _atom_size()) {
    //
    ErrOut err_out;

    err_out << funame << "number of atoms in the potential definition and in the molecular specification mismatch: "
	    << _atom_size() << ", " << mol.size();
  }

  // configuration grid points for numerical derivatives
  //
  for(SymIndex sin(_atom_size() * 3); sin.size() <= Graph::potex_max; ++sin) {
    //
    NumDer::cmap_t cmap = NumDer(_atom.size() * 3).convert(sin);

    for(NumDer::cmap_t::const_iterator cmit = cmap.begin(); cmit != cmap.end(); ++cmit) {
      //
      _conf_pool.insert(cmit->first);
    }
  }

  // individual mass square root
  //
  std::vector<double> imass_sqrt(_atom_size());

  // cumulative mass square root
  //
  std::vector<double> cmass_sqrt(_atom_size());

  dtemp = 0.;
  //
  for(int i = 0; i < mol.size(); ++i) {
    //
    dtemp += mol[i].mass();

    imass_sqrt[i] = std::sqrt(mol[i].mass());

    cmass_sqrt[i] = std::sqrt(dtemp);
  }

  // excluded cm motion basis
  //
  _cm_orth.resize(_atom_size(), _atom_size() - 1);

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

  set_potex(_start_pos);

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

    Graph::potex_t::const_iterator pit = potex.find(sin);

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
void Thermo::CSpec::set_potex (const pos_t& pos) const
{
  const char funame [] = "Thermo::CSpec::set_potex: ";

  double dtemp;
  int    itemp;

  _potex.clear();

  _freq.clear();
  
  // function-on-the-grid data for numerical derivatives
  //
  NumDer::fdata_t fdata;

  _set_fdata(pos, fdata);
  
  const int cart_size = _atom.size() * 3;
  
  // numerical potential derivatives
  //
  Graph::potex_t pot_der;

  for(SymIndex sin(cart_size); sin.size() <= Graph::potex_max; ++sin)
    //
    pot_der[sin] = NumDer(cart_size)(sin, fdata, numd_step);
  

  // mass-waited hessian without cm motion
  //
  Lapack::SymmetricMatrix mwfc(cart_size - 3);

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

	  Graph::potex_t::const_iterator pit = pot_der.find(gin);
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
  Lapack::Matrix tran(cart_size, cart_size - 3);
  
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
  for(SymIndex sin(cart_size - 3); sin.size() <= Graph::potex_max; ++sin) {
    //
    potex[sin] = 0.;

    Graph::potex_t::iterator pexit = potex.find(sin);
    
    for(GenIndex gin(sin.size(), cart_size); !gin.fin(); ++gin) {
      //
      Graph::potex_t::const_iterator pit = pot_der.find(gin);

      if(pit == pot_der.end()) {
	//
	ErrOut err_out;

	err_out << funame << "potential expansion term signature not found";
      }

      dtemp = pit->second;

      for(int i = 0; i < sin.size(); ++i)
	//
	dtemp *= tran(gin(i), sin(i));

      pexit->second += dtemp;
      
    }// potential expansion in cartesian coordinates
    //
  }// potential expansion in normal modes coordinates
}

// potential-on-the-grid data for numerical derivative calculations
//
void Thermo::CSpec::_set_fdata(const pos_t& pos, NumDer::fdata_t& fdata) const
{
  fdata.clear();
  
  Coord::Cartesian temp_pos(_atom_size * 3);
    
  for(std::set<std::vector<int> >::const_iterator cit = _conf_pool.begin(); cit != _conf_pool.end(); ++cit) {
    //
    temp_pos  = pos;

    for(int i = 0; i < pos.size(); ++i) {
      //
      temp_pos[i] += (double)(*cit)[i] * numd_step;
    }
    
    fdata[*cit] = _pot->evaluate(temp_pos);
  }
}

double Thermo::SpecBase::anharmonic_correction (double temperature) const
{
  const char funame [] = "Thermo::SpecBase::anharmonic_correction: ";

  IO::Marker funame_marker(funame);
  
  double dtemp;
  int    itemp;

  if(temperature <= 0.) {
    //
    ErrOut err_out;

    err_out << funame << "negative temperature[K]: " << temperature / Phys_const::kelv;
  }
    
  // anharmonic correction
  //
  std::map<int, double> res = Graph::Expansion(_freq, _potex).centroid_correction(temperature);

  IO::log << IO::log_offset << "anharmonic correction:\n";

  IO::log << IO::log_offset << std::setw(2) << "BO"<< std::setw(15) << "Value" << "\n";
  
  for(std::map<int, double>::const_iterator git = res.begin(); git != res.end(); ++git) {
    //
    IO::log << IO::log_offset << std::setw(2) << git->first << std::setw(15) << git->second << "\n";
  }

  IO::log << "\n";

  return res.rbegin()->second;
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

  IO::log << IO::log_offset << "Temperature = " << temperature / Phys_const::kelv << "K\n\n";
  
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
    
    _print_geom(pos);
  
    ener = new_ener;

    if(count++ < count_min)
      //
      continue;

    // frequencies & potential expansion in normal mode coordinates
    //
    std::vector<double> freq;

    Graph::potex_t potex;
    
    _set_potex(pos, freq, potex);

    IO::log << IO::log_offset << "frequencies, 1/cm:";

    for(int i = 0; i < freq.size(); ++i)
      //
      IO::log << "   " << freq[i] / Phys_const::incm;

    IO::log << "\n\n";

    if(freq[0] / temperature / 2. / M_PI < deep_tun_thresh) {
      //
      IO::log << IO::log_offset 
	      << "Oops: deep tunneling regime: temperature = "
	      << temperature / Phys_const::kelv
	      << "K\n";

      ++fail_count;
      
      continue;
    }

    double qfac;

    _corr_t anharm_corr = _correction(freq, potex, temperature, qfac);

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
    //
  }// sampling loop
  
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

/****************************************************************************************************************
 ********************************* LOCAL HARMONIC APPROXIMATION CORRECTION (QFACTOR) ****************************
 ****************************************************************************************************************/

// low  frequency threshold in qfactor
//
double  Thermo::QFactor::low_freq_thresh  = 1.e-3;

// high frequency threshold in qfactor
//
double  Thermo::QFactor::high_freq_thresh =  5.;

// deep tunneling threshold
//
double  Thermo::QFactor::deep_tunnel_thresh =  0.8;

Thermo::QFactor::QFactor (const std::vector<double>& freq, double temperature)
{
  const char funame [] = "Thermo::QFactor::QFactor: ";

  double dtemp;
  int    itemp;

  if(temperature <= 0.) {
    //
    std::cerr << funame << "temperature out of range: " << temperature / Phys_const::kelv;

    throw Error::Range();
  }
  
  _shift = 0.;

  _factor = 1.;

  for(int f = 0; f < freq.size(); ++f) {
    //
    dtemp = freq[f] / temperature / 2.;

    if(dtemp < -deep_tunnel_thresh * M_PI) {

      DeepTunnel x;
        
      x << funame
        << "deep tunneling regime: frequency[1/cm] = "
        << freq[f] / Phys_const::incm
        << "   temperature[K] = "
        << temperature / Phys_const::kelv
        << "\n";
      
      IO::log << IO::log_offset << x;
      
      throw x;
    }

    if(dtemp < -low_freq_thresh) {
      //
      _factor *= dtemp / std::sin(dtemp);
    }
    else if(dtemp > high_freq_thresh) {
      //
      _factor *= 2. * dtemp;

      _shift += freq[f] / 2.;
    }
    else if(dtemp > low_freq_thresh) {
      //
      _factor *= dtemp / std::sinh(dtemp);
    }
  }
}

