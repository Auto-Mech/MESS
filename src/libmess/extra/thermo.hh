/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef THERMO_HH
#define THERMO_HH

#include<vector>
#include<map>
#include<set>

#include "lapack.hh"
#include "io.hh"
#include "harding.hh"
#include "atom.hh"
#include "key.hh"

namespace Thermo {

  // symmetric tensor index
  //
  class SymIndex : private std::vector<int> {
    //
    int _range;
    
    void _isinit () const;

    SymIndex ();

  public:

    explicit SymIndex(int r) : _range(r) { _isinit(); }

    void operator++ ();
    void operator++ (int) { operator++(); }

    int range () const { return _range; }

    int size () const { return std::vector<int>::size(); }

    int operator() (int i) const { return (*this)[i]; }

    operator std::multiset<int> () const;
    //
    operator std::map<int, int> () const;
  };

  inline void SymIndex::_isinit () const
  {
    const char funame [] = "Thermo::SymIndex::_isinit: ";
    
    if(_range <= 0) {
      ErrOut err_out;
      err_out << funame << "not initialized properly";
    }
  }

  // generic fixed size multi-index 
  //
  class GenIndex : private std::vector<int> {
    //
    int _range;

    bool _fin;
    
    void _isinit () const;

    GenIndex ();

  public:

    GenIndex(int s, int r) : std::vector<int>(s), _range(r), _fin(false) { _isinit(); }

    void operator++ ();
    void operator++ (int) { operator++(); }

    int range () const { return _range; }

    int size () const { return std::vector<int>::size(); }

    bool fin () const { return _fin; }

    int operator() (int i) const { return (*this)[i]; }

    operator std::multiset<int> () const;
    //
    operator std::map<int, int> () const;
  };

  inline void GenIndex::_isinit () const
  {
    const char funame [] = "Thermo::GenIndex::_isinit: ";
    
    if(_range <= 0) {
      ErrOut err_out;
      err_out << funame << "not initialized properly";
    }
  }

  // numerical derivative generator
  //
  class NumDer {
    //
    // configurational space dimensionality
    //
    int _size;
    
    void _isinit ()      const;

  public:
    //
    // derivative order
    //
    typedef std::map<int, int>               der_t;

    // configuration map
    //
    typedef std::map<std::vector<int>, int>  cmap_t;

    NumDer () : _size(0) {}
    
    explicit NumDer (int s) : _size(s) { _isinit(); }

    void resize (int s) { _size = s; _isinit(); }

    int size () const { return _size; }

    // converts derivative signature to configuration signatures map
    //
    cmap_t operator() (const der_t& der) const;
  };

  inline void NumDer::_isinit () const
  {
    const char funame [] = "Thermo::NumDer::_isinit: ";

    if(_size <= 0) {
      ErrOut err_out;
      err_out << funame << "not initialized properly";
    }
  }

  /************************************************************************************
   ********************************* POTENTIAL ENERGY SURFACE *************************
   ************************************************************************************/

  typedef std::vector<double> pos_t;
  
  class PotBase {
    
  public:
    //
    // potential energy value
    //
    virtual double operator() (const  pos_t& pos) const = 0;
    
    // number of atoms
    //
    virtual int size () const = 0;
  };

  class Harding : public PotBase {
    //
    Harding ();

  public:
    //
    explicit Harding (std::istream& from) { harding_init_(from); }

    double operator() (const  pos_t& pos) const;
    
    int size () const { int res; harding_atom_size_(res); return res; }
  };

  class PotWrap {
    //
    ConstSharedPointer<PotBase> _pot;

    void _isinit () const;
    
  public:
    //
    void init (const PotBase* p);
    
    PotWrap () {}
    PotWrap (const PotBase* p) : _pot(p) {}

    operator bool () { return _pot; }

    // potential energy value
    //
    double operator() (const pos_t& pos) const { _isinit(); return _pot->operator()(pos); }

    // number of atoms
    //
    int size () const { _isinit(); return _pot->size(); }
  };
  
  inline void PotWrap::init (const PotBase* p)
  {
    const char funame [] = "Thermo::PotWrap::init: ";

    if(!p) {
      //
      ErrOut err_out;

      err_out << funame << "zero pointer";
    }
    
    if(_pot) {
      //
      ErrOut err_out;

      err_out << funame << "already initialized";
    }

    _pot.init(p);
  }
  
  inline void PotWrap::_isinit () const
  {
    const char funame [] = "Thermo::PotWrap::_isinit: ";
    
    if(!_pot) {
      //
      ErrOut err_out;

      err_out << funame << "not initialized";
    }
  }

  PotWrap new_pot (std::istream&);

  /************************************************************************************
   ********************************* GEOMETRICAL CONSTRAIN ****************************
   ************************************************************************************/

  class ConBase {
  public:
    virtual bool operator () (const pos_t& pos) const = 0;
  };
  
  class ConWrap {
    //
    ConstSharedPointer<ConBase> _con;

    void _isinit () const;
    
  public:
    void init (const ConBase* p);
    
    ConWrap () {}
    ConWrap (const ConBase* p) : _con(p) {}

    operator bool () { return _con; }
    
    bool operator () (const pos_t& pos) const { _isinit(); return _con->operator()(pos); }
  };

  
  inline void ConWrap::init (const ConBase* p)
  {
    const char funame [] = "Thermo::ConWrap::init: ";

    if(!p) {
      //
      ErrOut err_out;

      err_out << funame << "zero pointer";
    }
    
    if(_con) {
      //
      ErrOut err_out;

      err_out << funame << "already initialized";
    }

    _con.init(p);
  }
  
  inline void ConWrap::_isinit () const
  {
    const char funame [] = "Thermo::ConWrap::_isinit: ";
    
    if(!_con) {
      //
      ErrOut err_out;

      err_out << funame << "not initialized";
    }
  }
    
  class Species {
    //
    // deep tunneling exception
    //
    class _deep_t : public Exception::Base {
    
    public:
      //
      template <typename T>
      
      _deep_t& operator<< (const T& t) { (Exception::Base&)(*this) << t; return *this; } 
    };

    // potential
    //
    PotWrap _pot;

    // geometrical constrain
    //
    ConWrap _con;

    // atom specification and reference geometry
    //
    std::vector<Atom> _atom;

    // reference frequencies
    //
    std::vector<double> _ref_freq;

    // configuration signatures set for numerical derivatives
    //
    std::set<std::vector<int> > _conf_pool;

    // basis orthogonal to cm motion
    //
    Lapack::Matrix _cm_orth;

    // importance sampling starting configuration
    //
    pos_t _start_pos;

    typedef std::map<std::multiset<int>, double> _potex_t;

    // potential expansion in normal mode coordinates with cm motion excluded 
    //
    _potex_t _set_potex (const pos_t&                              cent_pos, 
			 const std::map<std::vector<int>, double>& conf_value, 
			 int                                       potex_max, 
			 std::vector<double>&                      freq
			 ) const;

    // frequencies only
    //
    std::vector<double> _frequencies (const pos_t& pos) const;
    
    // quantum factor omega * beta / 2 / sinh(omega * beta / 2) referenced to minimal energy configuration
    //
    double _qfactor (const std::vector<double>& freq, double temperature) const;

    typedef std::map<int, double> _corr_t;

    // harmonic and anharmonic corrections
    //
    _corr_t _correction (const pos_t& pos, double temperature, double& qfac) const;

    void _print_geom (const pos_t&) const;

  public:
    //
    void init (const std::vector<Atom>&, PotWrap, ConWrap = ConWrap());

    Species   (const std::vector<Atom>&, PotWrap, ConWrap = ConWrap());
    
    // partition function harmonic and total corrections
    //
    _corr_t weight (double temperature, double& harm_corr) const;

    // maximal number of importance samplings
    //
    static int count_max;

    // number of importance samplings to thermolize
    //
    static int count_min;

    // numerical derivative step (bohr)
    //
    static double numd_step;

    // importance sampling step (bohr)
    //
    static double samp_step;

    // low  frequency threshold (for qfactor)
    //
    static double  low_freq_thresh;

    // high frequency threshold (for qfactor)
    //
    static double high_freq_thresh;

    // deep tunneling threshold
    //
    static double deep_tun_thresh;

    //
    //
  };// Species 

  inline Species::Species (const std::vector<Atom>& a, PotWrap p, ConWrap c) { init(a, p, c); }  
    
}
#endif
