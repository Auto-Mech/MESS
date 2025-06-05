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

#ifndef MESS_HH
#define MESS_HH

#include <string>
#include <vector>
#include <set>
#include <map>

#include "error.hh"
#include "lapack.hh"
#include "model.hh"

namespace MasterEquation {
  //
  extern double temperature;

  extern double pressure;

  /*********************************************************************************************
   ******************************** USER DEFINED PARAMETERS ************************************
   *********************************************************************************************/
  
  // output
  //
  enum {TORR, BAR, ATM};                          // pressure units
  
  extern int pressure_unit;
  
  // multi-precision 
  //
  void set_precision (int);

  int get_precision ();

  extern int     use_mp;                          // use multi-precision for global matrix diagonalization

  // microcanonical rate maximum, 1/sec
  //
  extern double  micro_rate_max;

  // use kinetic matrix inversion for bimolecular-to-bimolecular rates calculation
  //
  extern    int  use_matrix_inversion;

  // well reduction
  //
  extern double  reduction_threshold;             // fast isomerization threshold

  extern double  reduction_tolerance;             // slow isomerization limit

  extern double  reduction_increment;             // minimal gap between the highest slow isomerization eigenvalue and the active one
  
  // chemical and relaxational subspaces
  //
  extern double  chemical_threshold;               // threshold separating chemical and relaxational eigenstates
  
  extern double  chemical_tolerance;               // smallest chemical eigenvalue

  // well partitioning
  //
  extern double  well_projection_threshold;        // single well projection threshold on the chemical subspace

  extern double  well_projection_tolerance;        // minimal single well projection on the chemical subspace

  // discretization and potential correction parameters
  //
  extern double  energy_step_over_temperature;

  extern double  excess_energy_over_temperature;

  extern double  dist_width_factor;                // thermal distribution width factor for excess energy  
  
  extern double  well_cutoff;                      // well cutoff parameter in temperature units

  extern double  well_extension;                   // well extension position
  
  extern double  well_ext_corr;                    // well extension correction parameter

  // time-propagation parameters in collision frequency units
  //
  extern double   time_step;

  extern double   time_limit;
  
  extern double   time_output;

  /****************************************************************************************************
   ************************************** CLASSES AND FUNCTIONS ***************************************
   ****************************************************************************************************/
  
  inline double energy_step () { return energy_step_over_temperature * temperature; }
  
  inline double thermal_factor (int e) { return std::exp((double)e * energy_step_over_temperature); }

  inline double thermal_factor_sqrt (int e) { return std::exp((double)e *  energy_step_over_temperature / 2.); }

  typedef std::set<int>            group_t;
  
  typedef std::vector<group_t> partition_t;

  inline group_t operator| (const group_t& g1, const group_t& g2)
  {
    group_t res = g2;

    for(group_t::const_iterator w1 = g1.begin(); w1 != g1.end(); ++w1)
      //
      res.insert(*w1);

    return res;
  }

  inline group_t operator+ (const group_t& g1, const group_t& g2) { return g1 | g2; }
  
  inline group_t operator& (const group_t& g1, const group_t& g2)
  {
    group_t res;

    for(group_t::const_iterator w1 = g1.begin(); w1 != g1.end(); ++w1)
      //
      if(g2.find(*w1) != g2.end())
	//
	res.insert(*w1);

    return res;
  }

  inline group_t operator- (const group_t& g1, const group_t& g2)
  {
    group_t res;

    for(group_t::const_iterator w1 = g1.begin(); w1 != g1.end(); ++w1)
      //
      if(g2.find(*w1) == g2.end())
	//
	res.insert(*w1);

    return res;
  }

  inline group_t& operator-= (group_t& g1, const group_t& g2)
  {
    int itemp = 0;

    while(itemp != g1.size()) {
      //
      itemp = g1.size();
      
      for(group_t::iterator w1 = g1.begin(); w1 != g1.end(); ++w1)
	//
	if(g2.find(*w1) != g2.end()) {
	  //
	  g1.erase(w1);

	  break;
	}
    }

    return g1;
  }

  /**************************************************************************************************
   ***************************************** REACTIVE COMPLEX ***************************************
   **************************************************************************************************/

  struct RateData;
  
  class ReactiveComplex: public Model::ChemGraph {
    //
    int _ener_index_begin;
    
    void _assert (const group_t&) const;
    
    std::map<int, int> _well_index_map;
    std::vector<int>   _index_well_map;

    std::map<int, int> _bim_index_map;
    std::vector<int>   _index_bim_map;

  public:
    //
    int kernel_bandwidth_max;

    double energy_bin (int e) const { return (_ener_index_begin - e) * energy_step(); }

    // constructor
    //
    ReactiveComplex (const ChemGraph&, int, int =0);
    
    /***************************************************************
     ********************** WELLS AND BARRIERS *********************
     ***************************************************************/
    
    // well and barrier base class
    //
    class PesObject {
      //
      int                   _type;
      
      int                  _index;
      
      Array<double>       _states;

      double              _weight;

      const ReactiveComplex*  _rc;
      
    public:
      //
      PesObject () : _weight(-1.) {}

      void init(int t, int i, const ReactiveComplex& rc) { _type = t; _index = i; _rc = &rc; }

      void set ();

      bool iset () const { if(_weight < 0.) return false; return true; }
      
      void assert (int) const;
      
      int    type () const { return _type;   }
      int   index () const { return _index;  }

      int   size () const { return _states.size(); }
      
      void  resize (int e)  { _states.resize(e); }

      double states (int e)  const    { assert(e); return _states[e]; }
      
      void   states (int e, double s) { assert(e);        _states[e] = s; }
      
      double states_sqrt   (int e) const { assert(e); return std::sqrt(_states[e]); }

      double      weight () const;

      double weight_sqrt () const { return std::sqrt(weight()); }

      double real_weight () const;

      double ground      () const { return _rc->energy_bin(_states.size()); }

      double real_ground () const;
      
      std::string name () const;

      const ReactiveComplex& reactive_complex () const { return *_rc; }
    };

    // well class
    //
    class Well: public PesObject {
      //
      double _collision_frequency;

      Lapack::Matrix _kernel;
      
      std::vector<double> _kernel_fraction;
  
    public:
      //
      Well () {}

      void init (int w, const ReactiveComplex& rc);

      void set ();
      
      double collision_frequency () const { return _collision_frequency; }

      double kernel (int i, int j) const { assert(i); assert(j); return _kernel(i, j); }

      int kernel_bandwidth;
    };

    // barrier class
    //
    class Barrier: public PesObject {
      //
    public:
      //
      Barrier () {}
      
      void init (int t, int b, const ReactiveComplex& rc) { PesObject::init(t, b, rc); }

      // number of states with microcanonical rate limit
      //
      double states (int e) const;

      void states (int e, double s) { PesObject::states(e, s); }
      
      std::pair<int, int> connect;
    };

    // well array
    //
    std::vector<Well>  well;

    int well_size () const { return well.size(); }
    
    int         well_to_index (int) const;
    group_t     well_to_index (const group_t&) const;
    partition_t well_to_index (const partition_t&) const;
    
    int         index_to_well (int) const;
    group_t     index_to_well (const group_t&) const;
    partition_t index_to_well (const partition_t&) const;

    std::string name (const     group_t& g) const;
    std::string name (const partition_t& p) const;

    const std::map<int, int>& well_index_map () { return _well_index_map; }
    const std::vector<int>&   index_well_map () { return _index_well_map; }
    
    // inner barrier array
    //
    std::vector<Barrier> inner_barrier;

    int inner_size () const { return inner_barrier.size(); }
    
    // outer barrier array
    //
    std::vector<Barrier> outer_barrier;

    int outer_size () const { return outer_barrier.size(); }
    
    // bimolecular-to-index maps
    //
    int bim_to_index (int) const;
    int index_to_bim (int) const;

    const std::map<int, int>& bim_index_map () { return _bim_index_map; }
    const std::vector<int>&   index_bim_map () { return _index_bim_map; }
    
    int bimolecular_size () const { return _index_bim_map.size(); }

    const Model::Bimolecular& bimolecular (int p) const { return Model::bimolecular(index_to_bim(p)); }

    /***********************************************************************************************
     **************************************** KINETIC BASIS ****************************************
     ***********************************************************************************************/

    // slow isomerization kinetic basis
    //
    struct KineticBasis {
      //
      int active_size;

      Lapack::Vector eigenvalue;
  
      Lapack::Matrix eigenvector;
      
      std::map<int, int> well_index_map;
  
      std::vector<int>   index_well_map;

      int size () const { return well_index_map.size(); }

      int     index_to_well (int)            const;
      int     well_to_index (int)            const;
      group_t index_to_well (const group_t&) const;
      group_t well_to_index (const group_t&) const;
    };

    int full_energy_range;
    
    // fast isomerization group
    //
    std::set<int>  xg;

    std::vector<KineticBasis> kinetic_basis;

    mutable std::vector<int>  well_shift;

    mutable int global_size;
    
    mutable Lapack::SymmetricMatrix kin_mat;

    // set kinetic basis and kinetic matrix
    //
    void  set_kinetic_matrix ();
    
    // get population in different wells in kinetic basis
    //
    Lapack::Vector population (double* state) const;

    // flux through bimolecular channels
    //
    Lapack::Vector escape_flux (Lapack::Vector) const;
    
    // populations as a result of state vector propagation in kinetic basis
    //
    Lapack::Vector propagate (Lapack::Vector state, Lapack::Vector escape = Lapack::Vector()) const;

    // rate partitioning
    //
    Lapack::Vector rate_partition ();

    /*****************************************************************************
     ************************** RATE CALCULATION METHODS **************************
     ******************************************************************************/

    void         well_reduction_method (RateData&, int flags =0) const;
  
    void direct_diagonalization_method (RateData&, int flags =0) const;

    void bound_species_rates (const Matrix<double>&, Lapack::Vector, Lapack::Matrix, RateData&, int) const;

    Lapack::SymmetricMatrix fast_bb_rate () const;
    
  };// reactive complex

  double default_energy_reference (const Model::ChemGraph&);
  
  /*********************************************************************************************
   ************************************* WELL PARTITIONING *************************************
   *********************************************************************************************/

  template<typename T>
  class WellSet {
    //
    Array<T> _weight;

    Array<T> _weight_sqrt;

    void _assert ()               const;
    void _assert (const group_t&) const;
    void _assert (const Matrix<T>&) const;
    
  public:
    //
    enum { INTERNAL = 1 };

    template<typename V>
    WellSet (const V&);
    
    int size () const { return _weight.size(); }
    
    T projection (const     group_t&, const Matrix<T>&) const;
    T projection (const partition_t&, const Matrix<T>&) const;

    T  weight (int          i) const { return _weight[i]; }
    T  weight (const group_t&) const;
    
    std::vector<double> weight (const partition_t&) const;

    Lapack::Vector basis (const     group_t&) const;
    Lapack::Matrix basis (const partition_t&) const;

    double threshold_well_partition (const Matrix<T>&, partition_t&, group_t&, int = 0) const;
  };
  
  template<typename T>
  template<typename V>
  inline WellSet<T>::WellSet (const V& v) : _weight(v)
  {
    _assert();

    _weight_sqrt.resize(_weight.size());

    for(int i = 0; i < _weight.size(); ++i)
      //
      _weight_sqrt[i] = sqrt(_weight[i]);
  }
    
  template<typename T>
  inline void WellSet<T>::_assert () const
  {
    const char funame [] = " MasterEquation::WellSet::assert: ";

    double dtemp;
    int    itemp;
  
    if(!size()) {
      //
      IO::log << IO::log_offset << funame << "empty set\n";

      throw Error::Init();
    }

    for(int i = 0; i < size(); ++i)
      //
      if(_weight[i] <= 0.) {
	//
	IO::log << IO::log_offset << funame << "non-positive weight: " << _weight[i] << "\n";

	throw Error::Range();
      }
  }

  template<typename T>
  inline void WellSet<T>::_assert (const group_t& g) const
  {
    if(!g.size()) {
      //
      IO::log << IO::log_offset << "group is empty\n";

      throw Error::Init();
    }
	
    if(*g.begin() < 0 || *g.rbegin() >= size()) {
      //
      IO::log << IO::log_offset << "well indices out of range:";

      for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
	//
	IO::log << "   " << *w;

      IO::log << "\n";
	
      throw Error::Range();
    }
  }

  template<typename T>
  inline void WellSet<T>::_assert (const Matrix<T>& pop_chem) const
  {
    if(pop_chem.size1() != size()) {
      //
      IO::log << IO::log_offset << "first dimension of the pop_chem matrix = " << pop_chem.size1()
	//
	      << " differs from the number of wells = " << size() << "\n";
    
      throw Error::Range();
    }
      
    if(pop_chem.size2() > size()) {
      //
      IO::log << IO::log_offset << "second dimension of the pop_chem matrix (chemically active species #) = " << pop_chem.size2()
	//
	      << " larger than the number of wells = " << size() << "\n";
    
      throw Error::Range();
    }      
  }

  template<typename T>
  inline T WellSet<T>::projection (const group_t& g, const Matrix<T>& pop_chem) const
  {
    const char funame [] = "MasterEquation::WellSet::projection: ";

    T dtemp;

    _assert(g);

    _assert(pop_chem);
  
    const int chem_size = pop_chem.size2();
  
    if(g.size() == 1)
      //
      return vdot(pop_chem.row(*g.begin()));

    T res = 0.;

    for(int l = 0; l < chem_size; ++l) {
      //
      dtemp = 0.;
    
      for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
	//
	dtemp += _weight_sqrt[*w] * pop_chem(*w, l);
    
      res += dtemp * dtemp;
    }

    // normalization
    //
    res /= weight(g);

    return res;
  }

  template<typename T>
  inline T WellSet<T>::projection (const partition_t& p, const Matrix<T>& pop_chem) const 
  {
    T res = 0.;
  
    for(int g = 0; g < p.size(); ++g)
    
      res += projection(p[g], pop_chem);

    return res;
  }

  template<typename T>
  inline T WellSet<T>::weight (const group_t& g) const
  {
    _assert(g);
  
    T res = 0.;
  
    for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
      //
      res += _weight[*w];

    return res;
  }

  template<typename T>
  inline std::vector<double> WellSet<T>::weight (const partition_t& p) const
  {
    std::vector<double> res(p.size());
  
    for(int g = 0; g < p.size(); ++g)
      //
      res[g] = (double)weight(p[g]);

    return res;
  }

  template<typename T>
  inline Lapack::Vector WellSet<T>::basis (const group_t& g) const
  {
    _assert(g);

    Lapack::Vector res(size(), 0.);
  
    if(g.size() == 1) {
      //
      res[*g.begin()] = 1.;
    
      return res;
    }
    
    double dtemp = std::sqrt((double)weight(g));

    for(group_t::const_iterator w = g.begin(); w != g.end(); ++w)
      //
      res[*w] = (double)_weight_sqrt[*w] / dtemp;
  
    return res;
  }

  template<typename T>
  inline Lapack::Matrix WellSet<T>::basis (const partition_t& p) const 
  {
    Lapack::Matrix res(size(), p.size());
  
    for(int g = 0; g < p.size(); ++g)
      //
      res.column(g) = basis(p[g]);
  
    return res;
  }

  // Generator of the partition of n items into m groups
  //
  class PartitionGenerator {
    //
    std::vector<int> _group_index;
    
    std::vector<int> _frame;
    
    bool             _end;
    
    const int        _partition_size;
    
  public:
    //
    PartitionGenerator (int m, int n);

    void operator++();
    
    void operator++(int) { operator++(); }
    
    bool end            ()      const { return _end; }
    
    int size            ()      const { return _group_index.size(); }
    
    int partition_size  ()      const { return _partition_size; }
    
    int operator[]      (int i) const { return _group_index[i]; }

    std::vector<std::set<int> > partition (const std::vector<int>& = std::vector<int> ()) const;
  };

  template<typename T>
  inline double WellSet<T>::threshold_well_partition (const Matrix<T>& pop_chem, partition_t& wp, group_t& bg, int flags) const
  {
    const char funame [] = "MasterEquation::WellSet::threshold_well_partition: ";

    int  itemp;
    T    dtemp;
    bool btemp;

    _assert(pop_chem);
  
    bg.clear();

    const int chem_size = pop_chem.size2();

    if(chem_size == size()) {
      //
      wp.resize(size());
    
      for(int i = 0; i < size(); ++i)
	//
	wp[i].insert(i);

      return (double)chem_size - (double)projection(wp, pop_chem);
    }

    if(well_projection_tolerance > 0. && flags & INTERNAL)
      //
      IO::log << IO::log_offset << funame << "WARNING: no exit channels: projection tolerance ignored\n";
  
    // large projection wells
    //
    std::vector<int> prim_well_map;

    std::set<int> exclude_group;
  
    for(int w = 0; w < size(); ++w) {
      //
      std::set<int> g;

      g.insert(w);

      dtemp = projection(g, pop_chem);

      if(dtemp > well_projection_threshold) {
	//
	prim_well_map.push_back(w);
      }
      else {
	//
	bg.insert(w);

	if(dtemp < well_projection_tolerance && !(flags & INTERNAL))
	  //
	  exclude_group.insert(w);
      }
    }
    if(prim_well_map.size() < chem_size) {
      //
      IO::log << IO::log_offset << funame << "primary wells # is less than chemical species #: decrease well partition threshold\n";

      throw Error::Range();
    }
  
    T proj = -1.;
  
    for(PartitionGenerator pg(chem_size, prim_well_map.size()); !pg.end(); ++pg) {
      //
      // new partition
      //
      partition_t ptemp = pg.partition(prim_well_map);

      // partition projection
      //
      dtemp = projection(ptemp, pop_chem);

      if(dtemp > proj) {
	//
	proj = dtemp;
      
	wp   = ptemp;
      }
    }

    while(bg.size()) {
      //
      int gi, wi;
    
      proj = -1.;
    
      for(int s = 0; s < chem_size; ++s) {
	//
	T pref = projection(wp[s], pop_chem);
	
	for(group_t::const_iterator w = bg.begin(); w !=bg.end(); ++w) {
	  //
	  if(exclude_group.find(*w) != exclude_group.end())
	    //
	    continue;
      
	  group_t g = wp[s];

	  g.insert(*w);
	
	  dtemp = projection(g, pop_chem) - pref;
	
	  if(dtemp > proj) {
	    //
	    proj = dtemp;
	  
	    gi = s;
	  
	    wi = *w;
	  }
	}
      }
    
      if(proj < 0. && !(flags & INTERNAL))
	//
	break;

      wp[gi].insert(wi);
    
      bg.erase(wi);
    }
  
    return (double)chem_size - (double)projection(wp, pop_chem);
  }
  
  /******************************************************************************
   ***************************** MANAGE RATE DATA *******************************
   ******************************************************************************/

  struct RateData {
    //
    std::map<std::set<int>, int> group_index_map;
    
    Array<int>          index_bim_map;

    partition_t        well_partition;

    Lapack::Matrix            ww_rate;

    Lapack::Matrix            wb_rate;

    Lapack::Matrix            bw_rate;
  
    Lapack::SymmetricMatrix   bb_rate;
  };

  // get reaction rates recursively
  //
  void get_rate_data (std::list<RateData>&, Model::ChemGraph = Model::ChemGraph(), int eref = 0, int flags = 0);

  void check_rate_data (const std::list<RateData>&);
  
  void print_rate_data (const std::list<RateData>&);
  
  void aggregate_rate_data (const std::list<RateData>&, RateData&);

  /**************************************************************
   *************************** HELPERS **************************
   **************************************************************/
  
  //  Generator of group of m elements from the pool of n elements
  //
  class GroupGenerator {
    //
    std::vector<int> _index;
    bool             _end;
    int              _size;
    
  public:
    //
    GroupGenerator (int m, int n);

    void operator++ ();
    void operator++ (int) { operator++(); }

    bool end       ()      const { return _end; }
    int size       ()      const { return _size; }
    int operator[] (int i) const { return _index[i]; }
  };

}

#endif
