/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2025, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#ifndef CROSSRATE_HH
#define CROSSRATE_HH

#include <cmath>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <fstream>

#include "divsur.hh"
#include "logical.hh"
#include "io.hh"
#include "potential.hh"
#include "trajectory.hh"

/*************************************************************
 *      namespace for classical trajectory calculations      *
 *      of the cross rates  between multiple bound and       *
 *      a single bimolecular products                        *
 ************************************************************/

namespace CrossRate {
  //
  // calculation type
  //
  enum job_t {
    TEST_JOB, // test configuration
    STAT_JOB, // statistical calculation
    DYN_JOB   // dynamical calculation
  };

  job_t job ();

  // calculation mode
  //
  enum mode_t { T_MODE, // canonical
		E_MODE, // microcanonical
		J_MODE  // E,J-resolved
  };

  // output
  //
  int rad_en_out ();

  mode_t mode ();

  // temperature for thermal rate calculation
  //
  double temperature ();

  // energy and angular momentum for E and E,J-resolved calculations
  //
  double energy_value ();
  double amom_value ();

  // reactive potential energy
  //
  double   reactive_energy ();
  void set_reactive_energy (double);
  
  // Normalization factor for the flux
  //
  double norm_factor ();

  // reactant species;
  //
  int reactant ();

  int seed ();

  void init (std::istream&);
  bool isinit ();
  
  void print_progress (int&, int);

  // test if the configuration is in a given species region
  //
  class SpecCondition : public Dynamic::Condition
  {
    const DivSur::MultiSur& _ms;
    int _spec;

  public:
    //
    SpecCondition (const DivSur::MultiSur& ms, int s) : _ms(ms), _spec(s) {}
    
    bool test (const Dynamic::Coordinates& dc) const { return _ms.test(_spec, dc); }
  };

  // results of the trajectory propagation
  //
  struct DynRes : public Trajectory::Propagator {
    //
    enum {INIT, DIRECT, RECROSS, POT_FAIL, RUN_FAIL, PASS};

    int stat; // status

    DynRes (const Potential::Wrap& pot, const Dynamic::Vars& dv, int dir) : Trajectory::Propagator(pot, dv, dir), stat(INIT) {}
  };

  // surface sampling with the results of the trajectory propagation, backward and forward
  //
  class DynSmp : public Dynamic::Vars {
    //
    double _ranval;    // random value between 0 and 1
    double _weight;    // statitstical weight
    double _energy;    // configuration potential energy

  public:
    //
    enum {INIT_ERR = 1,LOW_ENER_ERR = 2};

    class Exception {};

    class LowEnerErr : public Exception {};
    class    InitErr : public Exception {};

    DynSmp () {}

    void init (const Potential::Wrap&, const DivSur::MultiSur&, int, const Dynamic::Coordinates&);

    DynSmp (const Potential::Wrap& pot, const DivSur::MultiSur& ms, int prim, const Dynamic::Coordinates& dc) { init(pot, ms, prim, dc); }

    // statistical methods
    //
    double ranval           () const { return _ranval; }
    double weight           () const { return _weight; }
    double potential_energy () const { return _energy; }
    double     total_energy () const { return potential_energy() + total_kinetic_energy(); }

    // dynamical propagators
    //
    SharedPointer<DynRes> forw;
    SharedPointer<DynRes> back;

    ConstSharedPointer<DynRes> dyn_res (int dir) const;
    
    void run_traj (const DivSur::MultiSur&, const DivSur::face_t&, const Dynamic::CCP&);

    bool is_run      ()    const;
    bool is_pot_fail ()    const;
    bool is_run_fail ()    const;
    bool is_exclude  (int) const;

    double ener_deviation (int dir) const;

    double amom_deviation (int dir) const;
    
    void print_traj_results () const;
  };
  
  inline ConstSharedPointer<DynRes> DynSmp::dyn_res (int dir) const
  {
    const char funame [] = "CrossRate::DynSmp::dyn_res: ";

    switch(dir) {
      //
    case BACKWARD:
      //
      return back;

    case FORWARD:
      //
      return forw;

    default:
      //
      std::cerr << funame << "wrong case: " << dir << "\n";

      throw Error::Logic();
    }
  }
    
  // sampling result
  //
  struct SmpRes : public DynSmp, public DivSur::SmpRes {
    //
    enum {SUCCESS, FAKE, FAIL};

    int smp_stat;

    int error;

    SmpRes () : error(0) {}
  };

  inline double DynSmp::ener_deviation (int dir) const
  {
    const char funame [] = "CrossRate::DynSmp::ener_deviation: ";

    double res = dyn_res(dir)->final_energy() - total_energy();

    return res < 0. ? -res : res;
  }

  inline double DynSmp::amom_deviation (int dir) const
  {
    const char funame [] = "CrossRate::DynSmp::amom_deviation: ";

    return (dyn_res(dir)->total_angular_momentum() - total_angular_momentum()).vlength();
  }

  inline bool DynSmp::is_run () const
  { 
    if(forw->stat == DynRes::INIT && back->stat == DynRes::INIT) {
      //
      return false;
    }
    
    return true;
  }

  inline bool DynSmp::is_pot_fail () const
  { 
    if(forw->stat == DynRes::POT_FAIL ||  back->stat == DynRes::POT_FAIL) {
      //
      return true;
    }
    
    return false;
  }

  inline bool DynSmp::is_run_fail () const
  { 
    if(forw->stat == DynRes::RUN_FAIL ||  back->stat == DynRes::RUN_FAIL) {
      //
      return true;
    }
    
    return false;
  }

  inline bool DynSmp::is_exclude (int ward) const
  {
    const char funame [] = "CrossRate::DynSmp::is_exclude: ";

    ConstSharedPointer<DynRes> dr = dyn_res(ward);

    if(is_run_fail() || is_pot_fail() || dr->stat == DynRes::INIT
       //
       || !Dynamic::exclude_region || !Dynamic::exclude_region->test(*dr))
      //
      return false;

    return true;
  }

  /******************************************************************
   * Importance samplings array for a  facet; The faset is defined  *
   * as a surface part which separates two different species (the   *
   * direction is important)                                        *
   ******************************************************************/

  class FacetArray : public std::list<DynSmp>
  {
    long _flux_num; // # of potential energy samplings
    long _fail_num; // # of potential energy samplings failed
    long _fake_num; // # of samplings for which potential energy calculation was skipped

    double _flux; // cumulative flux value 
    double _fvar; // cumulative flax variation

    double               _min_ener;
    Dynamic::Coordinates _min_geom;

    double               _max_weight;
    
  public:
    //
    FacetArray () : _flux_num(0), _fail_num(0), _fake_num(0), _flux(0.), _fvar(0.), _min_ener(-100.) {}

    void merge (const FacetArray&);

    long flux_num () const { return _flux_num; }
    long fail_num () const { return _fail_num; }
    long fake_num () const { return _fake_num; }
    long face_num () const { return _flux_num + _fail_num + _fake_num; }
    long samp_num () const;

    double                      min_ener () const;
    const Dynamic::Coordinates& min_geom () const { return _min_geom; }

    void add_fail () { ++_fail_num; }
    void add_fake () { ++_fake_num; }
    void add_smp (const DynSmp&);

    double flux_value   () const;                                      // flux value
    double flux_var     () const;                                      // flux variance
    double flux_rel_var () const;                                      // flux relative variance
    double flux_rel_err () const { return std::sqrt(flux_rel_var()); } // flux relative error

    void run_traj (const DivSur::MultiSur&, const DivSur::face_t&, Dynamic::CCP);

    int    init_traj_num ()				const;
    int potfail_traj_num ()				const;
    int runfail_traj_num ()				const;
    int exclude_traj_num ()				const;
    int     run_traj_num ()				const;
    int  direct_traj_num ()				const;
    int recross_traj_num (int)				const;
    int  direct_traj_num (int, int* =0)			const;
    int    reac_traj_num (int, int, int* =0, int* =0)	const;

    double dyn_fac (int ward, int spec) const; // recrossing factor
    double dyn_dev (int ward, int spec) const; // recrossing factor standard deviation
    double dyn_var (int ward, int spec) const; // recrossing factor variance
  };

  inline double FacetArray::min_ener () const 
  {
    if(!samp_num()) {
      //
      std::cout << "FacetArray::min_ener: minimal energy has not been initialized yet => STOP\n";
      
      throw Error::Logic();
    }
    return _min_ener; 
  }

  /*****************************************************************
   * Set of facet samplings arrays for a surface; Faset is defined *
   * as the surface region which separates two different species   *
   *****************************************************************/

  class SurArray : public std::map<DivSur::face_t, FacetArray>
  {
    long _fail;       // # of logical failures
    long _inner;      // # of species inner regions hits
    long _exclude;    // # of excluded region hits
    long _close;      // # of close atoms hits

  public:
    //
    SurArray () : _fail(0), _inner(0), _exclude(0), _close(0) {} 

    int merge (const SurArray&);

    long  fail_num () const { return _fail; }
    long inner_num () const { return _inner; }
    long  excl_num () const { return _exclude; }
    long close_num () const { return _close; }

    void add_fail  () { ++_fail;    }
    void add_inner () { ++_inner;   }
    void add_excl  () { ++_exclude; }
    void add_close () { ++_close;   }

    long tot_smp_num () const;

    // Facet volume fraction, its standard deviation, and relative error
    //
    double vol_frac    (const_iterator sit) const;
    double vol_rmsd    (const_iterator sit) const;
    double vol_rel_var (const_iterator sit) const;

    double face_flux (const_iterator sit) const { return vol_frac(sit) * sit->second.flux_value(); }
  };

  inline long SurArray::tot_smp_num() const
  {
    long res = _fail + _inner + _exclude + _close;

    for(const_iterator cit = begin(); cit != end(); ++cit)
      //
      res += cit->second.face_num();
    
    return res;
  }

  inline double SurArray::vol_frac (const_iterator sit) const
  {
    long itemp = sit->second.face_num();
    
    if(!itemp)
      //
      return 0.;
    
    return (double)(itemp) / (double)(tot_smp_num());
  }

  inline double SurArray::vol_rmsd (const_iterator sit) const
  {
    double dtemp;
    
    long   itemp = sit->second.face_num();
    
    if(!itemp)
      //
      return 0.;
    
    dtemp = (double)(itemp) / (double)(tot_smp_num());
    
    return std::sqrt(dtemp * (1. - dtemp));
  }

  inline double SurArray::vol_rel_var (const_iterator sit) const
  {
    double dtemp;    
    long   itemp = sit->second.face_num();
    
    if(!itemp)
      //
      return 0.;
    
    return  1. / (double)(itemp)- 1. / (double)(tot_smp_num());    
  }

  // samplings data & trajectory results for the multiple species dividing surface
  //
  class MultiArray : public std::vector<SurArray>
  {
    const DivSur::MultiSur& _ms;
    
    Potential::Wrap         _pot;

    Dynamic::CCP            _stop;
    
    long _stat_pot_hit, _traj_pot_hit;

    int _sample (iterator, const std::set<DivSur::face_t>&);

    bool _work ();

    void _get_stat_flux(std::vector<double>&, 
			std::vector<double>&,
			std::vector<double>&) const;

    void _get_reac_flux (std::map<DivSur::face_t, double>&, 
			 std::map<DivSur::face_t, double>&,
			 std::map<DivSur::face_t, double>&) const;

    // print samplings results
    //
    void _print_sampling_results () const;

    // print trajectories results
    //
    void _print_traj_results () const;
    void _analize () const;

  public:
    //
    static int  min_sur_size; // minimum number of surface samplings to find facets
    static int  min_imp_size; // minimum importance samplings array size for the facet
    static int  max_imp_size; // maximum number of importance samplings for the facet; 
    static int  min_pot_size; // minimum number of facet samplings before the accuracy can be estimated
    static int  smp_set_size; // number of samplings performed at once

    static double face_rel_tol;  // uniform relative tolerance for a facet statistical flux
    static double spec_rel_tol;  // relative tolerance for a species statistical flux
    static double reac_rel_tol;  // relative tolerance for a species reactive flux
    static double tran_rel_tol;  // relative tolerance for the reactive transition flux

    static double ener_dev_max;  // trajectory energy deviation threshold
    static double amom_dev_max;  // trajectory angular momentum deviation threshold
    
    static std::vector<int> reactive_transition; // reactive transition for tolerance evaluation

    MultiArray(const DivSur::MultiSur&, Potential::Wrap, Dynamic::CCP);

    // run trajectories which have not been run
    //
    void run_traj ();
   //
  };
}// CrossRate

#endif

