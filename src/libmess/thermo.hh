/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2017, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
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
#include "multindex.hh"
#include "coord.hh"
#include "graph_common.hh"

namespace Thermo {

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
    
    // number of coordinates
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
    
    int size () const { int res; harding_atom_size_(res); return res * 3; }
  };

  // potential in cartesian coordinates
  //
  class CartPot: public PotBase {
    //
    ConstSharedPointer<Coord::CartFun> _pot;

  public:
    //
    virtual double operator() (const  pos_t& pos) const;
    
    virtual int size () const;

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

  /*****************************************************************************************
   ***************************** QUANTUM ANHARMONIC CORRECTION *****************************
   *****************************************************************************************/
  
  // deep tunneling exception
  //
  class DeepTunnel : public Exception::Base {
    
  public:
    //
    DeepTunnel (const std::string& s) : Exception::Base(s) {}
    
    DeepTunnel (const char* s) : Exception::Base(s) {}
    
    template <typename T>
    //
    DeepTunnel& operator<< (const T& t) { Exception::Base::operator<<(t); return *this; } 
  };

  // base class for anharmonic correction calculation
  //
  class SpecBase {
    //
  protected:
    //
    Graph::potex_t      _potex;
    
    std::vector<double> _freq;

  public:
    //
    virtual void set_potex (const pos_t& pos) = 0;
    
    double anharmonic_correction (double temperature) const;

    const std::vector<double>& frequencies () const { return _freq; }

    bool isinit () const { return _freq.size(); }
  };

  // initial potential expansion in cartesian coordinates with center-of-mass motion excluded
  //
  class CSpec: public SpecBase {
    //
    ConstSharedPointer<Coord::CartFun> _pot;

    int _atom_size () const { return _pot->atom_size(); }

    void _isinit () const;
    
    // configuration signatures set for numerical derivatives
    //
    std::set<std::vector<int> > _conf_pool;

    // basis orthogonal to cm motion
    //
    Lapack::Matrix _cm_orth;

    // potential-on-the-grid for numerical derivatives
    //
    void _set_fdata(const pos_t& pos, NumDer::fdata_t& fdata) const;
    
    void _print_geom (const pos_t&) const;

  public:
    //
    void init (const std::vector<Atom>&, ConstSharedPointer<Coord::CartFun>);
    
    CSpec (const std::vector<Atom>& mol, ConstSharedPointer<Coord::CartFun> p) { init(mol, p); }
    
    void set_potex (const pos_t& pos);
    
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

    //
    //
  };// CSpec

  inline void CSpec::_isinit () const
  {
    if(!_pot) {
      //
      std::cerr << "Thermo::CSpec: not initialized\n";

      throw Error::Init();
    }
  }

  // local harmonic approximation
  //
  class QFactor {
    //
    // energy shift
    //
    double _shift;

    // quantum correction factor
    //
    double _factor;

  public:

    QFactor (const std::vector<double>&, double);

    double shift () const { return _shift; }

    double factor () const { return _factor; }
    
    // low  frequency threshold 
    //
    static double  low_freq_thresh;

    // high frequency threshold 
    //
    static double high_freq_thresh;

    // deep tunneling threshold
    //
    static double deep_tunnel_thresh;
  };

}
#endif
