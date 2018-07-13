

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

    // calculation type
    enum job_t {
	TEST_JOB, // test configuration
	STAT_JOB, // statistical calculation
	DYN_JOB   // dynamical calculation
    };

    job_t job ();

  // calculation mode
  enum mode_t { T_MODE, // canonical
		E_MODE, // microcanonical
		J_MODE  // E,J-resolved
  };

  // output
  int rad_en_out ();

  mode_t mode ();

  // temperature for thermal rate calculation
  double temperature ();

  // energy and angular momentum for E and E,J-resolved calculations
  double energy_value ();
  double amom_value ();

  // reactive potential energy
  double   reactive_energy ();
  void set_reactive_energy (double);
  
  // Normalization factor for the flux
  double norm_factor ();

  // reactant species;
  int reactant ();

  int seed ();

  //bool isthermal (); // which calculation

  void init (std::istream&) throw(Error::General);
  bool isinit ();
  
  void print_progress (int&, int);

  extern std::ofstream xout;

  // test if the configuration is in a given species region
  class SpecCondition : public Dynamic::Condition
  {
    const DivSur::MultiSur& _ms;
    int _spec;

  public:
    SpecCondition (const DivSur::MultiSur& ms, int s) : _ms(ms), _spec(s) {}
    bool test (const Dynamic::Coordinates& dc) const { return _ms.species_test(_spec, dc); }
  };

  // results of the trajectory propagation
  struct DynRes : public Trajectory::Propagator {
    enum Stat {INIT, DIRECT, RECROSS, POT_FAIL, RUN_FAIL, EXCLUDE, PASS};

    Stat stat; // status

    DynRes (Potential::Wrap pot, const Dynamic::Vars& dv, int dir) : Trajectory::Propagator(pot, dv, dir), stat(INIT) {}
  };

  // surface sampling with the results of the 
  // trajectory propagation, backward and forward
  class DynSmp : public Dynamic::Vars {
    double _ranval;      // random value between 0 and 1
    double _weight;    // statitstical weight
    double _energy;    // configuration potential energy

  public:
    DynSmp (Potential::Wrap, const DivSur::MultiSur&, int, const Dynamic::Coordinates&) throw(Error::General);

    // statistical methods
    double ranval           () const { return _ranval; }
    double weight           () const { return _weight; }
    double potential_energy () const { return _energy; }
    double     total_energy () const { return potential_energy() + total_kinetic_energy(); }

    // dynamical methods
    SharedPointer<DynRes> forw;
    SharedPointer<DynRes> back;

    void run_traj (const DivSur::MultiSur&, const DivSur::face_t&, Dynamic::CCP);

    bool is_run      () const;
    bool is_pot_fail () const;
    bool is_run_fail () const;
    bool is_exclude  () const;
  };
  
  inline bool DynSmp::is_run () const
  { 
    if(forw->stat == DynRes::INIT && back->stat == DynRes::INIT)
      return false;
    else
      return true;
  }

  inline bool DynSmp::is_pot_fail () const
  { 
    if(forw->stat == DynRes::POT_FAIL ||  back->stat == DynRes::POT_FAIL)
      return true;
    else
      return false;
  }

  inline bool DynSmp::is_run_fail () const
  { 
    if(forw->stat == DynRes::RUN_FAIL ||  back->stat == DynRes::RUN_FAIL)
      return true;
    else
      return false;
  }

  inline bool DynSmp::is_exclude () const
  { 
    if(forw->stat == DynRes::EXCLUDE || back->stat == DynRes::EXCLUDE)
      return true;
    else
      return false;
  }

  /******************************************************************
   * Importance samplings array for a  facet; The faset is defined  *
   * as a surface part which separates two different species (the   *
   * direction is important)                                        *
   ******************************************************************/

  class FacetArray : private std::list<DynSmp>
  {
    int _flux_num; // # of potential energy samplings
    int _fail_num; // # of potential energy samplings failed
    int _fake_num; // # of samplings for which potential energy calculation was skipped

    double _max_weight;
    double _min_ener;
    Dynamic::Coordinates _min_geom;

    double _flux; // cumulative flux value 
    double _fvar; // cumulative flax variation

  public:
    FacetArray () : _flux_num(0), _fail_num(0), _fake_num(0), _flux(0.), _fvar(0.), _min_ener(-100.) {}

    int flux_num () const { return _flux_num; }
    int fail_num () const { return _fail_num; }
    int fake_num () const { return _fake_num; }
    int face_num () const { return _flux_num + _fail_num + _fake_num; }
    int samp_num () const;

    double min_ener () const throw(Error::General);
    const Dynamic::Coordinates& min_geom () const { return _min_geom; }

    void add_fake () { ++_fake_num; }
    bool add_smp  (Potential::Wrap, const DivSur::MultiSur&, int, const Dynamic::Coordinates&);

    // flux value
    double flux_val () const;
    // flux variance
    double flux_var    () const;
    // flux relative variance
    double flux_rel_var () const;
    // flux relative error
    double flux_rel_err () const { return std::sqrt(flux_rel_var()); }


    // importance samplings size
    int size () const { return std::list<DynSmp>::size(); }

    typedef std::list<DynSmp>::iterator iterator;
    typedef std::list<DynSmp>::const_iterator const_iterator;

    iterator begin () { return std::list<DynSmp>::begin(); }
    iterator end   () { return std::list<DynSmp>::end(); }

    const_iterator begin () const { return std::list<DynSmp>::begin(); }
    const_iterator end   () const { return std::list<DynSmp>::end(); }

    void run_traj (const DivSur::MultiSur&, const DivSur::face_t&, Dynamic::CCP);

    int    init_traj_num ()                   const;
    int potfail_traj_num ()                   const;
    int runfail_traj_num ()                   const;
    int exclude_traj_num ()                   const;
    int     run_traj_num ()                   const;
    int recross_traj_num (int ward)           const;
    int    reac_traj_num (int ward, int spec) const;

    double dyn_fac (int ward, int spec) const; // recrossing factor
    double dyn_dev (int ward, int spec) const; // recrossing factor standard deviation
    double dyn_var (int ward, int spec) const; // recrossing factor variance

  };// FacetArray

  inline double FacetArray::min_ener () const throw(Error::General)
  {
    if(!samp_num()) {
      std::cout << "FacetArray::min_ener: minimal energy has not been initialized yet => STOP\n";
      throw Error::Logic();
    }
    return _min_ener; 
  }

  /*****************************************************************
   * Set of facet samplings arrays for a surface; Faset is defined *
   * as the surface part which separates two different species     *
   *****************************************************************/

  class SurArray : private std::map<DivSur::face_t, FacetArray>
  {
    int _fail;       // # of logical failures
    int _inner;      // # of species inner regions hits
    int _exclude;    // # of excluded region hits
    int _close;      // # of close atoms hits
    int _skip;       // # of skipped facet hits

  public:
    SurArray () : _fail(0), _inner(0), _exclude(0), _close(0), _skip(0) {} 

    int  fail_num () const { return _fail; }
    int inner_num () const { return _inner; }
    int  excl_num () const { return _exclude; }
    int close_num () const { return _close; }
    int  skip_num () const { return _skip; }

    void add_fail  () { ++_fail;    }
    void add_inner () { ++_inner;   }
    void add_excl  () { ++_exclude; }
    void add_close () { ++_close;   }
    void add_skip  () { ++_skip;    }

    int tot_smp_num () const;

    typedef std::map<DivSur::face_t, FacetArray>::iterator iterator; 
    typedef std::map<DivSur::face_t, FacetArray>::const_iterator const_iterator;
 
    FacetArray& operator[] (const DivSur::face_t face) 
    { return std::map<DivSur::face_t, FacetArray>::operator[](face); }

    iterator find (const DivSur::face_t& face) 
    { return std::map<DivSur::face_t, FacetArray>::find(face); }

    const_iterator find (const DivSur::face_t& face) const
    { return std::map<DivSur::face_t, FacetArray>::find(face); }

    const_iterator begin () const { return std::map<DivSur::face_t, FacetArray>::begin(); }
    const_iterator end   () const { return std::map<DivSur::face_t, FacetArray>::end(); }

    iterator begin () { return std::map<DivSur::face_t, FacetArray>::begin(); }
    iterator end   () { return std::map<DivSur::face_t, FacetArray>::end(); }

    int size () const { return std::map<DivSur::face_t, FacetArray>::size(); }

    // Facet volume fraction, its standard deviation, and relative error
    double vol_frac    (const_iterator sit) const;
    double vol_rmsd    (const_iterator sit) const;
    double vol_rel_var (const_iterator sit) const;

    double face_flux (const_iterator sit) const { return vol_frac(sit) * sit->second.flux_val(); }
  };// SurArray

  inline int SurArray::tot_smp_num() const
  {
    int res = _fail + _inner + _exclude + _close + _skip;

    for(const_iterator cit = begin(); cit != end(); ++cit)
      res += cit->second.face_num();
    
    return res;
}

  inline double SurArray::vol_frac (const_iterator sit) const
  {
    int itemp = sit->second.face_num();
    if(itemp)
      return (double)(itemp) / (double)(tot_smp_num());
    return 0.;
  }

  inline double SurArray::vol_rmsd (const_iterator sit) const
  {
    int itemp = sit->second.face_num();
    if(!itemp)
      return 0.;
    double r = (double)(itemp) / (double)(tot_smp_num());
    return std::sqrt(r * (1. - r));
  }

  inline double SurArray::vol_rel_var (const_iterator sit) const
  {
    int itemp = sit->second.face_num();
    if(itemp)
      return  1./(double)(itemp)- 1./(double)(tot_smp_num());    
    return 0.;
  }


  // samplings data for the multiple species dividing surface
  class MultiArray : public  std::vector<SurArray>
  {
      const DivSur::MultiSur _ms;
      Potential::Wrap _pot;

    int _sample (iterator, const std::set<DivSur::face_t>&);

    bool _work (Dynamic::CCP) throw(Error::General);

    // print samplings results
    void _print_sampling_results () const;

    // print trajectories results
    void _print_traj_results () const;
    void _analize () const;

    void _get_stat_flux(std::vector<double>&, 
			std::vector<double>&,
			std::vector<double>&) const;

    void _get_reac_flux (std::map<DivSur::face_t, double>&, 
			 std::map<DivSur::face_t, double>&,
			 std::map<DivSur::face_t, double>&) const;

  public:

    static int  min_sur_size; // minimum number of surface samplings to find facets
    static int  min_imp_size; // minimum importance samplings array size for the facet
    static int  max_imp_size; // maximum number of importance samplings for the facet; 
    static int  min_pot_size; // minimum number of facet samplings before the accuracy can be estimated
    static int  max_pot_size; // maximum number of facet samplings

    static double face_rel_tol;  // uniform relative tolerance for a facet statistical flux
    static double spec_rel_tol;  // relative tolerance for a species statistical flux
    static double reac_rel_tol;  // relative tolerance for a species reactive flux
    static double tran_rel_tol;  // relative tolerance for the reactive transition flux

    static std::vector<int> reactive_transition; // reactive transition for tolerance evaluation

    MultiArray(const DivSur::MultiSur&, Potential::Wrap, Dynamic::CCP) throw(Error::General);

    // run trajectories which have not been run
    void run_traj (Dynamic::CCP) throw(Error::General);
  };// MultiArray

}// CrossRate

#endif

