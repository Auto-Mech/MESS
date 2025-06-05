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

#ifndef MESS_HH
#define MESS_HH

#include <string>
#include <vector>
#include <set>
#include <map>

#include "error.hh"
#include "lapack.hh"
#include "model.hh"

#if defined(WITH_MPACK) || defined(WITH_MPLAPACK)

#include "mpack_dd.hh"
#include "mplapack_dd.hh"

#define FLOAT  dd_real
#define LAPACK Mpack_dd
#define EXP    exp
#define SQRT   sqrt

#else

#define FLOAT double
#define LAPACK Lapack
#define EXP    std::exp
#define SQRT   std::sqrt

#endif


namespace MasterEquation {
  //
  class Group;
  class Partition;

  // eigenvector and eigenvalue output
  //
  extern std::ofstream eval_out;// eigenvalues output
  extern std::ofstream evec_out;// eigenvalues output
  extern int           evec_out_num;// number of relaxation eigenvalues to print

  extern std::string save_kinetic_matrix;
  
  // pressure units
  //
  enum {TORR, BAR, ATM};
  
  extern int pressure_unit;

  typedef FLOAT float_t;
  
  extern int use_mp;

  int get_precision ();

  void set_precision (int);

  // reduction of species
  //
  enum {DIAGONALIZATION, PROJECTION}; // possible reduction algorithms for low eigenvalue method
  extern int  red_out_num;// number of reduction schemes to print 
  extern int reduction_method;// reduction algorithm for low-eigenvalue method
  void set_default_partition (const std::vector<std::string>&) ;// default reduction scheme
  void set_default_chem_size (int);

  // capture probabilities
  //
  extern std::map<std::string, std::set<double> > hot_energy;

  // product energy distributions
  extern std::ofstream ped_out;
  void set_ped_pair(const std::vector<std::string>& ped_spec) ;

  // control parameters
  //
  extern double            reference_temperature;

  extern double            reference_pressure;
  
  extern bool              is_global_cutoff;

  extern double            lower_global_cutoff;  // lower global energy cutoff

  extern double            well_cutoff;// well cutoff parameter (lower part of the energy window)
  
  extern double            chemical_threshold;// threshold separating chemical and relaxational eigenstates
  
  extern double            chemical_tolerance;// smallest chemical eigenvalue
  
  extern double reduction_threshold;// relaxational eigenvalue to collision frequency threshold
  
  extern double            rate_max;// microcanonical rate maximum
  
  extern double      well_extension;// well extension parameter
  
  extern double       well_ext_corr;// well extension correction parameter
  
  extern double     energy_step_over_temperature;

  extern double   energy_window_over_temperature;

  extern double   excess_energy_over_temperature;// upper part of the energy window
  
  extern double   energy_cutoff_over_temperature;// lower part of the energy window

  extern double   time_propagation_step;// time-propagation step over collisional frequency time-scale
  
  extern double   time_propagation_limit;// time-propagation limit over collisional frequency time-scale
  
  extern double   time_propagation_interval; // time-interval between printouts
  
  double temperature            ();
  
  double pressure               ();
  
  double energy_step            ();
  
  double energy_reference       ();
  //double collision_frequency    ();

  void set_temperature      (double);
  void set_pressure         (double);
  void set_energy_reference (double);

  /**************************************************************************************************
   ***************************************** REACTIVE COMPLEX ***************************************
   **************************************************************************************************/

  // reactive complex
  //
  struct ReactiveComplex {
    //
    // global specifications
    //
    typedef std::map<std::set<int>, std::map<int, double> > bfdb_t;

    static void branching_fraction (double, double, const Model::ChemGraph&, const std::set<int>&, bfdb_t&, int = 0);

    typedef std::pair<std::list<std::set<int> >, std::set<int> > lump_t;

    static lump_t lumping_scheme (double temperature, double pressure, int f = 0);
  
    static lump_t lumping_scheme (double temperature, double pressure, const Model::ChemGraph& cg, std::pair<int, int>, int f = 0);
  
    // constructors
    //
    bool isinit;
    
    std::map<int, double> well_extension_cap;

    std::pair<int, int> control_barrier;
    
    ReactiveComplex () : isinit(false), control_barrier(std::make_pair((int)Model::UNKNOWN, -1)) {}
    
    void init       (double t, double p, double r, double c, const Model::ChemGraph& g, int f = 0);

    ReactiveComplex (double t, double p, double r, double c, const Model::ChemGraph& g, int f = 0);
    
    void init       (double t, double p, double r, double c, const Model::ChemGraph& g, std::pair<int, int> b, int f = 0);

    ReactiveComplex (double t, double p, double r, double c, const Model::ChemGraph& g, std::pair<int, int> b, int f = 0);

    void init       (double t, double p, double r, double c, const Model::ChemGraph& g, const std::map<int, double>& x, int f = 0);

    ReactiveComplex (double t, double p, double r, double c, const Model::ChemGraph& g, const std::map<int, double>& x, int f = 0);

    // well definition
    //
    struct Well {
      //
      bool isinit;
      
      Well () : isinit(false) {}
      
      void set (double t, double r, double c, const Model::Well& m, const std::map<int, double>& x, int f = 0);
      
      Well     (double t, double r, double c, const Model::Well& m, const std::map<int, double>& x, int f = 0);
      
      float_t weight;
      float_t weight_sqrt;
      
      std::vector<float_t> boltzman;
      std::vector<float_t> boltzman_sqrt;
      std::vector<float_t> state_density;

      int size () const { return state_density.size(); }

      // kernel stuff
      //
      typedef LAPACK::Matrix kernel_t;

      kernel_t kernel;

      int kernel_bandwidth;

      std::vector<float_t> kernel_escape;

      double collision_factor;

      double collision_frequency (double p) const { return collision_factor * p; }
    };

    // barrier definition
    //
    struct Barrier {
      //
      bool isinit;
      
      Barrier () : isinit(false) {}
      
      void set (double t, double r, int s, const Model::Species& m, int f = 0);

      Barrier  (double t, double r, int s, const Model::Species& m, int f = 0);
      
      float_t weight;
      
      std::vector<float_t> state_number;

      int size() const { return state_number.size(); }
    };

    Model::ChemGraph graph;
    
    double temperature;

    double pressure;

    double energy_reference;

    double energy_cutoff;

    double energy_step;

    // wells
    //
    std::vector<Well> _well;

    const Well& well (int i) const;
    
    int well_size () const { return _well.size(); }
    
    // well-to-index maps
    //
    std::map<int, int> _well_index_map;
    std::vector<int>   _index_well_map;

    int well_to_index (int) const;
    int index_to_well (int) const;

    std::set<int> well_to_index (const std::set<int>&) const;
    std::set<int> index_to_well (const std::set<int>&) const;
    
    std::vector<std::set<int> > well_to_index (const std::vector<std::set<int> >&) const;
    std::vector<std::set<int> > index_to_well (const std::vector<std::set<int> >&) const;

    const Model::Well& well_model (int i) const { return Model::well(index_to_well(i)); }

    // bimolecular-to-index maps
    //
    std::map<int, int> _bim_index_map;
    std::vector<int>   _index_bim_map;

    int bim_to_index (int) const;
    int index_to_bim (int) const;

    std::set<int> bim_to_index (const std::set<int>&) const;
    std::set<int> index_to_bim (const std::set<int>&) const;

    std::vector<std::set<int> > bim_to_index (const std::vector<std::set<int> >&) const;
    std::vector<std::set<int> > index_to_bim (const std::vector<std::set<int> >&) const;

    int bimolecular_size () const { return _index_bim_map.size(); }

    const Model::Bimolecular& bimolecular (int p) const { return Model::bimolecular(index_to_bim(p)); }

    // barriers
    //
    typedef std::map<std::set<int>,       Barrier> inner_t;
    typedef std::map<std::pair<int, int>, Barrier> outer_t;

    inner_t inner_barrier;
    outer_t outer_barrier;

    // dissociation channel
    //
    std::pair<int, Barrier> dissociation_channel;
    
    void set_dissociation_channel (std::pair<int, int>, int = 0);

    // thermal factor (Boltzman exponent)
    //
    static float_t thermal_factor (int);

    /*********************************************************************************************
     ************************************* WELL PARTITIONING *************************************
     *********************************************************************************************/

    typedef std::set<int> group_t;

    typedef std::vector<group_t> partition_t;

    void assert (const group_t&) const;

    float_t projection (const     group_t&, LAPACK::Matrix) const;
    float_t projection (const partition_t&, LAPACK::Matrix) const;

    float_t              group_weight (const     group_t&) const;
    std::vector<float_t> group_weight (const partition_t&) const;

    LAPACK::Vector basis (const     group_t&) const;
    LAPACK::Matrix basis (const partition_t&) const;

    float_t threshold_well_partition (LAPACK::Matrix pop_chem, partition_t&, group_t&, int flags = 0) const;

    enum {NO_WELL_EXTENSION = 1, NOPRINT = 2, USE_PROJECTION_TOLERANCE = 4};

    /***********************************************************************************************
     **************************************** KINETIC BASIS ****************************************
     ***********************************************************************************************/

    void set_kinetic_matrix (int = 0);
    
    struct KineticBasis {
      //
      int active_size;

      LAPACK::Vector eigenvalue;
  
      LAPACK::Matrix eigenvector;
      
      std::map<int, int> well_index_map;
  
      std::vector<int>   index_well_map;

      int size () const { return well_index_map.size(); }
    };

    std::vector<KineticBasis> kinetic_basis;

    // maximal energy index
    //
    int ener_index_max;

    // global index map
    //
    std::vector<int> well_shift;

    int global_size;

    // global kinetic relaxation matrix
    //
    LAPACK::Matrix kin_mat;

    /***********************************************************************************************
     ************************************* CALCULATION METHODS *************************************
     ***********************************************************************************************/

    // reaction type
    //
    typedef std::set<std::pair<int, int> > reac_t;

    struct RateData {
      //
      Model::ChemGraph        graph;
      
      std::map<int, int>      _bim_index_map;

      int bim_to_index (int) const;
      
      partition_t             well_partition;
 
      group_t                 bimolecular_group;
      
      Lapack::SymmetricMatrix bb_rate;

      Lapack::Matrix          ww_rate;

      Lapack::Matrix          wb_rate;

      Lapack::Matrix          bw_rate;

      bfdb_t                  bf;
      
      std::set<int>          missing_bf (const std::list<reac_t>&) const;

      std::set<int> add_rate_data (std::map<reac_t, double>& rate_data, const std::list<reac_t>& reactions, const bfdb_t& bfdb) const;
    };
      
    RateData well_reduction_method (int flags = 0) const;

    bool there_are_bound_groups () const;

    /*******************************************************************************************************
     ******************************** BRANCHING FRACTION BY TIME-PROPAGATION *******************************
     *******************************************************************************************************/
    
    std::map<std::pair<int, int>, double> branching_fraction (std::pair<int, int>, int = 0) const;

    std::map<int, double> branching_fraction (int = 0) const;
    
    std::map<std::pair<int, int>, double> propagate (const std::map<int, LAPACK::Vector>& pop_vector, int = 0) const;
    
    LAPACK::Vector escape_flux (LAPACK::Vector state) const;

    LAPACK::Vector  population (LAPACK::Vector state) const;

  };// reactive complex

  inline ReactiveComplex::Well::Well (double t, double r, double c, const Model::Well& m, const std::map<int, double>& x, int f)
	//
    : isinit(false) { set(t, r, c, m, x, f); }
      
  inline ReactiveComplex::Barrier::Barrier  (double t, double r, int s, const Model::Species& m, int f)
    //
    : isinit(false) { set(t, r, s, m, f); }
      
  inline void ReactiveComplex::init (double t, double p, double r, double c, const Model::ChemGraph& g, std::pair<int, int> b, int f)
  //
  { init(t, p, r, c, g, f); set_dissociation_channel(b, f); }
  
  inline void ReactiveComplex::init (double t, double p, double r, double c, const Model::ChemGraph& g, const std::map<int, double>& x, int f)
  //
  { well_extension_cap = x; init(t, p, r, c, g, f); }
  
  inline ReactiveComplex::ReactiveComplex (double t, double p, double r, double c, const Model::ChemGraph& g, int f)
  //
    : isinit(false) { init(t, p, r, c, g, f); set_kinetic_matrix(f); }
  
  inline ReactiveComplex::ReactiveComplex (double t, double p, double r, double c, const Model::ChemGraph& g, std::pair<int, int> b, int f)
  //
    : isinit(false) { init(t, p, r, c, g, b, f); set_kinetic_matrix(f); }
  
  inline ReactiveComplex::ReactiveComplex (double t, double p, double r, double c, const Model::ChemGraph& g, const std::map<int, double>& x, int f)
  //
    : isinit(false) { init(t, p, r, c, g, x, f); set_kinetic_matrix(f); }
  
  /************************************************************************************************/
  
  // set states densities, states numbers, etc.
  //
  enum {DEFAULT_EREF = 1};
  
  void set (std::map<std::pair<int, int>, double>& rate_data, std::map<int, double>& capture, int flags = 0);

  /********************************************************************************************
   ***************************************** WELL CLASS ***************************************
   ********************************************************************************************/

  // CRM = Collisional Relaxation Modes
  //
  // CER = Collisional Energy Relaxation
  //
  class Well {
    //
    double                  _weight;             // statistical weight
    double                  _weight_sqrt;        // square root of the statistical weight
    double                  _min_relax_eval;     // minimal collisional relaxation eigenvalue
    double                  _max_relax_eval;     // maximal collisional relaxation eigenvalue
    
    Lapack::Vector          _state_density;      // density of states on the grid
    Lapack::Vector          _boltzman;           // Boltzmann distribution
    Lapack::Vector          _boltzman_sqrt;      // Boltzmann distribution square root
    
    Lapack::Matrix          _crm_basis;          // CRM basis (ket)
    Lapack::Matrix          _crm_bra;            // CRM basis (bra)
    Lapack::Matrix          _kernel;             // energy relaxation kernel
    Lapack::SymmetricMatrix _crm_kernel;         // kernel in CRM basis
    //Lapack::Vector          _escape_rate;        // escape rate;

    double              _collision_factor;
    std::vector<double> _kernel_fraction;
    
    // radiational transitions
    Lapack::SymmetricMatrix     _radiation_rate;
    Lapack::SymmetricMatrix _crm_radiation_rate;

    void _set_state_density (const Model::Well&);
    void _set_kernel (const Model::Well&) ;
    void _set_crm_basis ();

  public:
    //
    explicit Well (const Model::Well&);

    int              size ()                const { return _state_density.size(); }
    double         weight ()                const { return               _weight; }
    double    weight_sqrt ()                const { return          _weight_sqrt; }
    double  state_density (int i)           const { return     _state_density[i]; }
    double       boltzman (int i)           const { return          _boltzman[i]; }
    double  boltzman_sqrt (int i)           const { return     _boltzman_sqrt[i]; }
    double         kernel (int i, int j)    const { return         _kernel(i, j); }

    double minimal_relaxation_eigenvalue () const { return       _min_relax_eval; }
    double maximal_relaxation_eigenvalue () const { return       _max_relax_eval; }

    int             crm_size ()             const { return      size() - 1;         }
    const double&  crm_basis (int i, int j) const { return      _crm_basis(i, j);   }
    double           crm_bra (int i, int j) const { return      _crm_bra(i, j);     }
    const double* crm_column (int i)        const { return     &_crm_basis(0, i);   }
    const double* crm_row    (int i)        const { return     &_crm_basis(i, 0);   }
    double crm_kernel (int i, int j) const { return      _crm_kernel(i, j);  }// kernel in CRM basis 

    //double escape_rate(int i)               const { return       _escape_rate[i]; }
    //const double* escape_rate()             const { return       _escape_rate; }

    // radiational transitions
    bool       radiation      ()             const { return     _radiation_rate.isinit(); }
    double     radiation_rate (int i, int j) const { return     _radiation_rate(i, j);    }
    double crm_radiation_rate (int i, int j) const { return _crm_radiation_rate(i, j);    }

    double collision_frequency     () const { return _collision_factor * pressure();    }
    double kernel_fraction    (int i) const { return _kernel_fraction[i];     }

    int kernel_bandwidth;
  };

  
  /********************************************************************************************
   **************************************** BARRIER CLASS *************************************
   ********************************************************************************************/

  class Barrier {
    //
    Lapack::Vector _state_number;
    double               _weight;

  public:
    //
    explicit Barrier  (const Model::Species&);
    
    int             size ()      const { return _state_number.size(); }
    double  state_number (int i) const { return _state_number[i]; }
    double& state_number (int i)       { return _state_number[i]; }
    double        weight ()      const { return _weight; }

    void      truncate (int) ;
  };

  /*********************************************************************************************/

  const Well&                            well (int w);
  const Barrier&                inner_barrier (int b);
  const Barrier&                outer_barrier (int b);

  /************************** RATE COEFFICIENTS CALCULATION METHODS ******************************/

  void divide_and_conquer_method (const Model::ChemGraph&);
  
  void low_eigenvalue_matrix (Lapack::SymmetricMatrix& k_11, Lapack::SymmetricMatrix& k_33,
			      //
			      Lapack::Matrix& k_13, Lapack::Matrix& l_21);

  typedef void             (*Method) (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  
  void         low_eigenvalue_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  
  void direct_diagonalization_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  
  void         well_reduction_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  
  void        high_pressure_analysis ();

  /************************** WELL PARTITION METHODS ******************************/

  void set_well_partition_method (const std::string&);

  extern double  well_projection_threshold;
  
  extern double  well_projection_tolerance;
  
  double   threshold_well_partition (Lapack::Matrix, Partition&, Group&);
  double        sort_well_partition (Lapack::Matrix, Partition&, Group&);
  double incremental_well_partition (Lapack::Matrix, Partition&, Group&);
  double  sequential_well_partition (Lapack::Matrix, Partition&, Group&);

  /******************************************* HELPERS *******************************************/

  // group of wells
  //
  class Group : public std::set<int> {
    //
    void _assert () const;
    
  public:
    
    Group () {}
    
    explicit Group (int w) { std::set<int>::insert(w); }

    bool insert  (const Group&);
    bool erase   (const Group&);
    bool contain (const Group&) const;
    
    bool insert (int i) { return std::set<int>::insert(i).second; }
    bool erase  (int i) { return std::set<int>::erase(i); }
    
    int    group_index ()       const;
    double      weight ()       const;
    double real_weight (double) const;
    double real_ground ()       const; 

    double projection (Lapack::Matrix) const;

    Lapack::Vector basis_vector () const;
  };

  class PartitionGenerator;
  
  class Partition : public std::vector<Group> {
    //
  public:
    //
    Partition () {}
    
    explicit Partition (int n) : std::vector<Group>(n) {}
    explicit Partition (const PartitionGenerator&, const std::vector<int>& = std::vector<int>());

    double projection (Lapack::Matrix) const;

    Lapack::Matrix      basis       ()       const;
    std::vector<int>    group_index ()       const;
    std::vector<double>      weight ()       const;
    std::vector<double> real_weight (double) const;
    std::vector<double> real_ground ()       const;
  };

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

  // Generator of the partition of n objects into m groups
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

  // description of the species as a group of equilibrated wells at high pressure
  //
  struct HPWell : public Group {
    double weight; // statistical weight
  };

  // description of the barrier which connect two species at high pressure
  //
  struct HPBarrier : public std::pair<int, int> {
    double weight; // statistical number of states (over 2 Pi)
  };

  // auxilliary eigenvalue output
  //
  extern std::ofstream eval_out; // eigenvalues output
}

#endif
