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
  class Group;
  class Partition;

  // eigenvector and eigenvalue output
  extern std::ofstream eval_out;// eigenvalues output
  extern std::ofstream evec_out;// eigenvalues output
  extern int           evec_out_num;// number of relaxation eigenvalues to print

  enum {TORR, BAR, ATM};
  extern int pressure_unit;

  //
  enum {DOUBLE, DD, QD};

  extern int float_type;
  
  // reduction of species
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

  extern double         well_cutoff;// well cutoff parameter
  extern double  chemical_threshold;// threshold separating chemical and relaxational eigenstates
  extern double       min_chem_eval;// smallest chemical eigenvalue
  extern double            rate_max;// microcanonical rate maximum
  extern double reduction_threshold;// maximal chemical eigenvalue to collision frequency ratio
  extern double      well_extension;
  extern double      well_ext_corr;

  extern bool well_reduction_correction;
  
  void set_global_cutoff (double);

  double temperature            ();
  double pressure               ();
  double energy_step            ();
  double energy_reference       ();
  //double collision_frequency    ();

  void set_temperature      (double);
  void set_pressure         (double);
  void set_energy_step      (double);
  void set_energy_reference (double);

  // set states densities, states numbers, etc.
  //
  enum {DEFAULT_EREF = 1};
  
  void set (std::map<std::pair<int, int>, double>& rate_data, std::map<int, double>& capture, int flags = 0) ;

  /********************************************************************************************
   ***************************************** WELL CLASS ***************************************
   ********************************************************************************************/

  // CRM = Collisional Relaxation Modes
  // CER = Collisional Energy Relaxation
  class Well {
    double                  _weight;             // statistical weight
    double                  _weight_sqrt;        // square root of the statistical weight
    double                  _real_weight;        // not truncated statistical weight
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
    explicit Well (const Model::Well&);

    int              size ()                const { return _state_density.size(); }
    double         weight ()                const { return               _weight; }
    double    real_weight ()                const { return          _real_weight; }
    double    weight_sqrt ()                const { return          _weight_sqrt; }
    double  state_density (int i)           const { return     _state_density[i]; }
    double       boltzman (int i)           const { return          _boltzman[i]; }
    double  boltzman_sqrt (int i)           const { return     _boltzman_sqrt[i]; }
    double         kernel (int i, int j)    const { return         _kernel(i, j); }

    double minimal_relaxation_eigenvalue () const { return       _min_relax_eval; }
    double maximal_relaxation_eigenvalue () const { return       _max_relax_eval; }

    int             crm_size ()             const { return      _crm_basis.size2(); }
    const double&  crm_basis (int i, int j) const { return      _crm_basis(i, j);   }
    double           crm_bra (int i, int j) const { return      _crm_bra(i, j);     }
    const double* crm_column (int i)        const { return     &_crm_basis(0, i);   }
    const double* crm_row    (int i)        const { return     &_crm_basis(i, 0);   }
    const double& crm_kernel (int i, int j) const { return      _crm_kernel(i, j);  }// kernel in CRM basis 

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
    Lapack::Vector _state_number;
    double               _weight;
    double          _real_weight;

  public:
    explicit Barrier  (const Model::Species&);
    int             size ()      const { return _state_number.size(); }
    double  state_number (int i) const { return _state_number[i]; }
    double& state_number (int i)       { return _state_number[i]; }
    double        weight ()      const { return _weight; }
    double   real_weight ()      const { return _real_weight; }

    void      truncate (int) ;
  };

  /********************************************************************************************
   ************************************ BIMOLECULAR CLASS *************************************
   ********************************************************************************************/

  class Bimolecular {
    double _weight;

  public:
    Bimolecular(const Model::Bimolecular& p);
    double   weight () const { return   _weight; }
  };

  /*********************************************************************************************/

  const Well&                            well (int w);
  const Bimolecular&              bimolecular (int p);
  const Barrier&                inner_barrier (int b);
  const Barrier&                outer_barrier (int b);

  /************************** RATE COEFFICIENTS CALCULATION METHODS ******************************/

  void low_eigenvalue_matrix (Lapack::SymmetricMatrix& k_11, Lapack::SymmetricMatrix& k_33, 
			      Lapack::Matrix& k_13, Lapack::Matrix& l_21) ;

  typedef void             (*Method) (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  void         low_eigenvalue_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  void direct_diagonalization_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  void         well_reduction_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  void         well_reduction_method_old (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);
  void             sequential_method (std::map<std::pair<int, int>, double>& rate_data, Partition& well_partition, int flags);

  void        high_pressure_analysis ();

  /************************** WELL PARTITION METHODS ******************************/

  void set_well_partition_method (const std::string&) ;

  extern double  well_projection_threshold;

  double   threshold_well_partition (const Lapack::Matrix&, Partition&, Group&, const std::vector<double>&)
    ;
  double        sort_well_partition (const Lapack::Matrix&, Partition&, Group&, const std::vector<double>&)
    ;
  double incremental_well_partition (const Lapack::Matrix&, Partition&, Group&, const std::vector<double>&)
    ;
  double  sequential_well_partition (const Lapack::Matrix&, Partition&, Group&, const std::vector<double>&)
    ;

  /******************************************* HELPERS *******************************************/

  // group of wells
  class Group : public std::set<int> {
  public:
    Group () {}
    explicit Group (int);

    Group insert (const Group&);
    Group erase  (const Group&);
    Group insert (int i) { std::set<int>::insert(i); return *this; }
    Group erase  (int i) { std::set<int>::erase(i); return *this; }
    bool contain (const Group&) const;

    double projection (const Lapack::Matrix&, const std::vector<double>& = std::vector<double>()) const
      ;
    double projection (const Lapack::Matrix&, int, const std::vector<Group>&, const std::vector<double>&) const 
      ;

    int operator[]    (int) const;
    Lapack::Vector basis () const;
    int    group_index () const;
    double      weight () const;
    double real_weight () const;
  };

  typedef Group::const_iterator Git;

  // partition of the wells
  //
  class PartitionGenerator;
  
  class Partition : public std::vector<Group> {
  public:
    explicit Partition (int n) : std::vector<Group>(n) {}
    explicit Partition (const PartitionGenerator&, const std::vector<int>& = std::vector<int>());
    Partition () {}

    double projection (const Lapack::Matrix&, const std::vector<double>& = std::vector<double>()) const;

    Lapack::Matrix basis () const;
    std::vector<int> group_index () const;
    std::vector<double> weight () const;
    std::vector<double> real_weight () const;
  };

  typedef Partition::const_iterator Pit;

  //  Generator of group of m elements from the pool of n elements
  class GroupGenerator {
    std::vector<int> _index;
    bool _end;
    int _size;
  public:
    GroupGenerator (int m, int n) ;

    void operator++ ();
    void operator++ (int) { operator++(); }

    bool end       ()      const { return _end; }
    int size       ()      const { return _size; }
    int operator[] (int i) const { return _index[i]; }
  };

  // Generator of the partition of n objects into m groups
  class PartitionGenerator {
    std::vector<int> _group_index;
    std::vector<int> _frame;
    bool             _end;
    const int        _partition_size;
  public:
    PartitionGenerator (int m, int n) ;

    void operator++();
    void operator++(int) { operator++(); }

    bool end            ()      const { return _end; }
    int size            ()      const { return _group_index.size(); }
    int partition_size  ()      const { return _partition_size; }
    int operator[]      (int i) const { return _group_index[i]; }
  };

  // description of the species as a group of equilibrated wells at high pressure
  struct HPWell : public Group {
    double weight; // statistical weight
  };

  // description of the barrier which connect two species at high pressure
  struct HPBarrier : public std::pair<int, int> {
    double weight; // statistical number of states (over 2 Pi)
  };

  // auxilliary output
  extern std::ofstream eval_out; // eigenvalues output
}

#endif
