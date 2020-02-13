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

#ifndef MODEL_HH
#define MODEL_HH

#include<vector>
#include<set>
#include<map>
//#include<multiset>

#include "graph_omp.hh"
#include "slatec.hh"
#include "shared.hh"
#include "lapack.hh"
#include "multindex.hh"
#include "atom.hh"
#include "math.hh"
#include "io.hh"
#include "system.hh"

namespace Model {

  enum {DENSITY, NUMBER, NOSTATES};// calculated properties
  enum {ANGSTROM, BOHR}; // distance units

  // minimal interatomic distance
  extern double atom_dist_min;

  // maximum energy to be used
  double  energy_limit () ;
  void set_energy_limit (double);
  bool is_energy_limit ();

  int escape_size       ();    // number of wells with escape channels
  int escape_well_index (int); // escape well index

  void shift_cm_to_zero(std::vector<Atom>&);
  Lapack::SymmetricMatrix inertia_moment_matrix(const std::vector<Atom>&);
  void read_geometry (IO::KeyBufferStream&, std::vector<Atom>&, int =ANGSTROM);

  /********************************************************************************************
   ******************************************** READERS ***************************************
   ********************************************************************************************/
  class Collision;
  class Kernel;
  class Species;
  class Bimolecular;
  class Rotor;
  class Core;
  class Tunnel;
  class Escape;

  SharedPointer<Collision>   new_collision      (IO::KeyBufferStream&);
  SharedPointer<Kernel>      new_kernel         (IO::KeyBufferStream&);
  SharedPointer<Tunnel>      new_tunnel         (IO::KeyBufferStream&);
  SharedPointer<Escape>      new_escape         (IO::KeyBufferStream&);
  SharedPointer<Bimolecular> new_bimolecular    (IO::KeyBufferStream&, const std::string&);
  SharedPointer<Species>     new_species        (IO::KeyBufferStream&, const std::string&, int);
  SharedPointer<Rotor>       new_rotor          (IO::KeyBufferStream&, const std::vector<Atom>&);
  SharedPointer<Core>        new_core           (IO::KeyBufferStream&, const std::vector<Atom>&, int);

  /********************************************************************************************
   ************************************* COLLISION MODEL **************************************
   ********************************************************************************************/

  class Collision {

  public:
    virtual ~Collision ();

    virtual double operator () (double temperature) const = 0;
  };

  class LennardJonesCollision : public Collision {
    double _frequency_factor;
    double _epsilon;

    double _omega_22_star (double) const;
  public:
    LennardJonesCollision (IO::KeyBufferStream&) ;
    ~LennardJonesCollision ();

    double operator () (double) const;
  };

  /********************************************************************************************
   ************************************ KERNEL CLASS FAMILY ***********************************
   ********************************************************************************************/

  // collisional energy transfer kernel
  class Kernel {
    static int _flags;

  public:
    Kernel () {}
    virtual ~Kernel ();
   
    enum {
      UP = 1,      // transition up probability form predefined
      DENSITY = 2, // transition probability is proportional to the final density of states
      NOTRUN  = 4  // no truncation even so the transition probability is negative
    };

    static int flags () { return _flags; }
    static void add_flag (int f) { _flags |= f; }

    virtual double    operator() (double ener, double temperature) const =0;
    virtual double cutoff_energy              (double temperature) const =0;
  };

  /********************************* EXPONENTIAL KERNEL MODEL **********************************/

  class ExponentialKernel : public Kernel {

    std::vector<double> _factor;
    std::vector<double> _power;
    std::vector<double> _fraction;

    double _cutoff;
    double _energy_down (int, double) const;
  public:

    ExponentialKernel (IO::KeyBufferStream&) ;
    ~ExponentialKernel ();

    double operator () (double, double) const;
    double cutoff_energy (double) const;
  };

  /**************************************************************************************
   ************************************* TUNNELING **************************************
   **************************************************************************************/

  class Tunnel {
    double    _wtol;// statistical weight tolerance
    static double _action_max; // maximum 

  protected:
    double _cutoff;// cutoff energy
    double   _freq;// imaginary frequency

    Tunnel(IO::KeyBufferStream&) ;
    
  public:
    virtual ~Tunnel ();

    double  cutoff () const { return _cutoff; }

    static double  action_max               () { return _action_max; }
    static void    set_action_max (double val) { _action_max = val; }

    double  factor (double) const; // tunneling factor
    double density (double) const; // energy derivative of tunneling factor
    double  weight (double) const; // statistical weight relative to cutoff energy
    void convolute (Array<double>&, double) const;// convolute number of states with tunneling density

    virtual double action (double, int =0) const =0; // semiclassical action
  };

  /**************************************************************************************
   ****************************** READ BARRIER TUNNELING ********************************
   **************************************************************************************/

  class ReadTunnel: public Tunnel {
    Slatec::Spline _action;

  public:
    ReadTunnel(IO::KeyBufferStream&) ;
    ~ReadTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /**************************************************************************************
   ************************** PARABOLIC BARRIER TUNNELING *******************************
   **************************************************************************************/

  class HarmonicTunnel: public Tunnel {
  public:
    HarmonicTunnel(IO::KeyBufferStream&) ;
    ~HarmonicTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /**************************************************************************************
   ************************** ECKART BARRIER TUNNELING ********************************
   **************************************************************************************/

  class EckartTunnel: public Tunnel {
    std::vector<double> _depth;
    double             _factor;

  public:
    EckartTunnel(IO::KeyBufferStream&) ;
    ~EckartTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /**************************************************************************************
   ************************** QUARTIC BARRIER TUNNELING ********************************
   **************************************************************************************/

  class QuarticTunnel: public Tunnel {
    double _vmin;// minimal well depth    
    double   _v3;// x^3 power expansion coefficient
    double   _v4;// x^4 power expansion coefficient

    double _potential (double x) { return x * x * (0.5  + _v3 * x + _v4 * x * x); }

    Slatec::Spline _action; // semiclassical action for energies below barrier 

    class XratioSearch : public Math::NewtonRaphsonSearch {
      double _vratio;

    public:
      XratioSearch(double v, double t) : _vratio(v) { tol = t; }
      double operator() (double, int) const;
    };

  public:
    QuarticTunnel(IO::KeyBufferStream&) ;
    ~QuarticTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /********************************************************************************************
   ******************************** INTERNAL ROTATION DEFINITION ******************************
   ********************************************************************************************/

  class InternalRotationBase {
    std::set<int>      _group; // moving group definition
    std::pair<int, int> _axis; // moving direction definition
    int             _symmetry; // internal rotation symmetry
    int                 _imax; // maximal atomic index in internal roation definition

    bool _isinit;

  protected:
    explicit InternalRotationBase (int s) : _isinit(false), _symmetry(s) {}
    InternalRotationBase (IO::KeyBufferStream&) ;
    virtual ~InternalRotationBase ();

  public:
    int symmetry () const { return _symmetry; }

    std::vector<Atom>            rotate (const std::vector<Atom>& atom, double angle) const;
    std::vector<D3::Vector> normal_mode (const std::vector<Atom>& atom, Lapack::Vector* =0) const;
  };

  /********************************************************************************************
   ************************* HINDERED ROTOR/UMBRELLA MODE CLASS FAMILY ************************
   ********************************************************************************************/

  class Rotor {

  protected:
    // controls
    int      _ham_size_max; // maximum dimension of the Hamiltonian
    int      _ham_size_min; // minimal dimension of the Hamiltonian
    int         _grid_size; // angular discretization size
    double  _therm_pow_max; // thermal exponent power maximum

    std::vector<Atom> _atom; // 3D structure

    Rotor ();
    Rotor (IO::KeyBufferStream&, const std::vector<Atom>&) ;

  public:
    virtual ~Rotor();

    virtual void   set (double) =0;// maximal energy relative to ground   

    virtual double ground        ()       const =0;// ground energy
    virtual double energy_level  (int)    const =0;// energy level relative to the ground
    virtual int    level_size    ()       const =0;// number of energy levels
    virtual double weight        (double) const =0;// statistical weight relative to the ground

    void convolute(Array<double>&, double) const;
  };

  // rotational constant for free and hindered rotors
  class RotorBase : public Rotor, public InternalRotationBase  {
    double  _rotational_constant;    

  protected:
    RotorBase (double r, int s) : InternalRotationBase(s), _rotational_constant(r) {} 
    RotorBase (IO::KeyBufferStream&, const std::vector<Atom>&) ;
    virtual ~RotorBase ();

  public:

    double rotational_constant () const { return _rotational_constant; }
  };

  /********************************************************************************************
   ****************************************** FREE ROTOR **************************************
   ********************************************************************************************/

  class FreeRotor : public RotorBase {
    int _level_size;

  public:
    FreeRotor (IO::KeyBufferStream&, const std::vector<Atom>&) ;
    ~FreeRotor ();

    void set (double);

    double ground        ()       const;
    double energy_level  (int)    const;
    int    level_size    ()       const;
    double weight        (double) const;

  };

  /********************************************************************************************
   **************************************** HINDERED ROTOR ************************************
   ********************************************************************************************/

  class HinderedRotor : public RotorBase {

    double                    _ground; // ground level energy
    std::vector<double> _energy_level; // energy levels relative to the ground
    std::map<int, double>   _pot_four; // potential fourier expansion

    std::vector<double>  _pot_grid; // potential-on-the-grid
    std::vector<double> _freq_grid; // frequency-on-the-grid
    double              _grid_step; // angular discretization step
    double                _pot_max; // global potential energy maximum
    double                _pot_min; // global potential energy minimum
    double               _freq_max; // maximum frequency
    double               _freq_min; // minimum frequency
    double            _harm_ground; // ground energy in harmonic approximation
    
    int _weight_output_temperature_step; // temperature step for statistical weight output
    int _weight_output_temperature_max;  // temperature maximum for statistical weight output
    int _weight_output_temperature_min;  // temperature maximum for statistical weight output
 
    void _set_energy_levels (int) ;
    void _read (IO::KeyBufferStream&);
    void _init ();

    bool _use_quantum_weight;

  public: 
    HinderedRotor (IO::KeyBufferStream&, const std::vector<Atom>&);
    HinderedRotor (const std::map<int, double>& p, double r, int s);
    ~HinderedRotor ();

    double                        potential           (double, int =0) const;
    int         semiclassical_states_number                   (double) const;
    Lapack::Vector real_space_energy_levels                         () const;
    double                   quantum_weight                   (double) const;
    int            get_semiclassical_weight (double, double&, double&) const;

    // integrate local energy distribution over angle
    void integrate (Array<double>&, double) const;
    double potential_minimum () { return _pot_min; }

    // virtual functions
    void set (double);

    double ground       ()       const;
    double energy_level (int)    const;
    int    level_size   ()       const;
    double weight       (double) const;
  };

  /********************************************************************************************
   ***************************************** UMBRELLA MODE ************************************
   ********************************************************************************************/

  class Umbrella : public Rotor {

    double _ground;
    std::vector<double> _energy_level;

    double _mass;
    Lapack::Vector _pot_coef;

    // discretization
    std::vector<double>  _pot_grid; // potential-on-the-grid
    std::vector<double> _freq_grid; // frequency-on-the-grid
    double                  _astep; // discretization step
    double                _pot_min; // global potential energy minimum
    
    static double _integral(int p, int n);// \int_0^1 dx x^p cos(n\pi x)
    void _set_energy_levels(int) ;

  public:
    Umbrella (IO::KeyBufferStream&, const std::vector<Atom>&) ;
    ~Umbrella ();

    double quantum_weight (double) const;
    int get_semiclassical_weight (double, double&, double&) const;
    double potential(double x, int der =0) const ;

    // virtual functions
    void set (double ener_max);

    double ground          () const; // ground energy
    double energy_level (int) const; // energy levels relative to the ground
    int    level_size      () const; // number of energy levels
    double weight    (double) const; // statistical weight relative to the ground
  };

  /**************************************************************************************
   ************************************* RRHO CORE **************************************
   **************************************************************************************/

  class Core {
    int _mode;
    Core ();

  protected:
    explicit Core(int m);

  public:
    virtual ~Core ();

    virtual double ground       () const =0;
    virtual double weight (double) const =0; // statistical weight relative to the ground
    virtual double states (double) const =0; // density or number of states relative to the ground

    int mode () const { return _mode; }
  };

  inline Core::Core (int m) : _mode(m) 
  {
    const char funame [] = "Model::Core::Core: ";

    if(mode() != NUMBER && mode() != DENSITY && mode() != NOSTATES) {
      std::cerr << funame << "wrong mode\n";
      throw Error::Logic();
    }
  }

  /**************************************************************************************
   *********************** PHASE SPACE THEORY NUMBER OF STATES **************************
   **************************************************************************************/

  class PhaseSpaceTheory : public Core {
    double _states_factor;
    double _weight_factor;
    double _power;

  public:
    PhaseSpaceTheory (IO::KeyBufferStream& from) ;
    ~PhaseSpaceTheory ();

    double ground       () const;
    double weight (double) const;
    double states (double) const;
  };

  
  /***************************************************************************************
   ************************************* UNHARMONICITY ************************************
   ***************************************************************************************/

  /**************************************************************************************
   ************************************* RIGID ROTOR ************************************
   **************************************************************************************/

  class RigidRotor : public Core {
    double            _factor;
    int               _rdim;
    double            _rofactor;
    
    double                                _ground;

    std::vector<double>                _frequency; // harmonic frequencies
    std::vector<int>                      _edegen; // electronic level degeneracies
    std::vector<int>                      _fdegen; // frequencies degeneracies
    Lapack::SymmetricMatrix               _anharm; // second order anharmonic energy level expansion
    std::vector<std::vector<double> >        _rvc; // rovibrational coupling

    // rovibrational expansion coefficients
    double _rovib (int ...) const;
    
    // interpolation
    double           _emax; // interpolation energy maximum
    double           _nmax; // extrapolation power value
    Slatec::Spline _states;

    double _core_states (double) const;
    double _core_weight (double) const;
    
  public:
    RigidRotor (IO::KeyBufferStream&, const std::vector<Atom>&, int) ;
    ~RigidRotor ();

    double ground ()       const;
    double weight (double) const;
    double states (double) const;
  };

  /*****************************************************************************************
   *************************** ROTATIONAL NUMBER OF STATES FROM ROTD ***********************
   *****************************************************************************************/
  
  class Rotd : public Core {

    Array<double>       _rotd_ener; // relative energy on the grid
    Array<double>       _rotd_nos; // density of states of the transitional modes on the grid

    double              _ground;
    Slatec::Spline      _rotd_spline;
    double              _rotd_emin, _rotd_emax;
    double              _rotd_nmin, _rotd_amin;
    double              _rotd_nmax, _rotd_amax;

  public:
    Rotd (IO::KeyBufferStream&, int) ;
    ~Rotd ();

    double ground       () const;
    double states (double) const; // density or the number of states relative to the ground
    double weight (double) const; // statistical weight relative to the ground
  };

  /********************************************************************************************
   ******************** INTERNAL ROTATION DEFINITION FOR MULTIROTOR CLASS *********************
   ********************************************************************************************/

  class InternalRotation : public InternalRotationBase {
    int                _msize; // mass fourier expansion size;
    int                _psize; // potential fourier expansion size
    int                _wsize; // angular grid size
    int                 _qmin; // minimum quantum state size
    int                 _qmax; // maximum quantum state size

  public:
    InternalRotation (IO::KeyBufferStream&) ;
    ~InternalRotation ();

    int      mass_fourier_size () const { return    _msize; }
    int potential_fourier_size () const { return    _psize; }
    int   weight_sampling_size () const { return    _wsize; }
    int       quantum_size_min () const { return     _qmin; }
    int       quantum_size_max () const { return     _qmax; }
  };

  /********************************************************************************************
   ******************************** INTERNAL ROTATION SAMPLING CLASS **************************
   ********************************************************************************************/
  
    class MultiRotorSampling: public Core {
      //
      std::vector<InternalRotation> _internal_rotation; // internal rotations description

      class Sampling {
	//
	double _weight_factor;

	double _pot_energy;

	std::vector<double> _frequency;

	int _dof;

      public:
	//
	Sampling (double wf, double pe, const std::vector<double>& f, int dof)
	  : _weight_factor(wf), _pot_energy(pe), _frequency(f), _dof(dof) {}

	double statistical_weight (double temperature) const;

	double states (double energy, int mode) const;
      };

      std::vector<Sampling> _sampling;

    };

  /********************************************************************************************
   ********************************* COUPLED INTERNAL ROTORS MODEL ****************************
   ********************************************************************************************/

  class MultiRotor: public Core {
    //
    std::vector<InternalRotation> _internal_rotation; // internal rotations description

    double  _external_symmetry; // external rotation symmetry number

    // internal & external mobilities fourier expansions
    //
    MultiIndexConvert                      _mass_index; // mass fourier expansion/sampling dimensions

    std::map<int, Lapack::SymmetricMatrix>     _imm_four; //           internal mobility matrix fourier expansion
    std::map<int, Lapack::SymmetricMatrix> _eff_imm_four; // effective internal mobility matrix fourier expansion

    std::map<int, double>     _ctf_four; // curvlinear transformation factor(dx/dx' determinant) fourier expansion
    std::map<int, double>     _irf_four; // internal rotation factor fourier expansion
    std::map<int, double>     _erf_four; // external rotation factor fourier expansion
    std::map<int, double> _eff_erf_four; // effective external rotation factor fourier expansion


    std::vector<Lapack::SymmetricMatrix> _internal_mobility_real;
    std::vector<Lapack::SymmetricMatrix> _external_mobility_real;
    std::vector<Lapack::Matrix>          _coriolis_coupling_real;
    Lapack::Vector                                     _ctf_real;

    std::map<int, Lapack::ComplexMatrix> _internal_mobility_fourier;// fourier transform
    std::map<int, Lapack::ComplexMatrix> _coriolis_coupling_fourier;// fourier transform
    std::map<int, Lapack::ComplexMatrix> _external_mobility_fourier;// fourier transform
    std::map<int, Lapack::complex>             _ctf_complex_fourier;

    // fourier expansion for potential and vibrational frequencies
    MultiIndexConvert                 _pot_four_index; // potential fourier expansion dimensions

    std::map<int, double>                   _pot_four; //           potential fourier expansion
    std::map<int, double>               _eff_pot_four; // effective potential fourier expansion

    std::vector<std::map<int, double> >     _vib_four; //           vibrational frequencies fourier expansion
    std::vector<std::map<int, double> > _eff_vib_four; // effective vibrational frequencies fourier expansion

    MultiIndexConvert                   _pot_index; // potential sampling dimensions

    Lapack::Vector                      _pot_real; // potential sampling data

    double                              _pot_shift; // constant part of the potential

    std::map<int, Lapack::complex> _pot_complex_fourier; // potential complex fourier transform

    // properties on the angular grid
    //
    MultiIndexConvert           _grid_index; // angular grid dimensions

    std::vector<double>          _pot_grid; // potential on the grid 
    std::vector<Lapack::Vector>  _vib_grid; // vibrational frequencies on the grid
    std::vector<Lapack::Vector> _freq_grid; // internal rotation frequencies on the grid
    std::vector<double>          _irf_grid; // square root of the internal mass determinant on the grid
    std::vector<double>          _erf_grid; // square root of the inertia moments product on the grid

    double                _angle_grid_cell; // angular grid cell volume
    std::vector<double>   _angle_grid_step; // angular grid step

    // quantum energy levels
    //
    std::vector<std::vector<double> > _energy_level; // quantum energy levels relative to the ground one

    std::vector<std::vector<double> >     _mean_erf; // state averaged external rotation factor [sqrt(I1*I2*I3)]

    double                                  _ground; // ground level energy

    // controls
    //
    bool   _with_ctf;               // include curvlinear transformation factor into hamiltonian
    bool   _with_ext_rot;           // include external rotation
    bool   _full_quantum_treatment; // all vibrational populational states calculated
    double _level_ener_max;         // maximum energy for quantum levels calculation relative to the potential minimum
    double _mtol;                   // mass matrix elements tolerance for pruning
    double _ptol;                   // potential matrix element tolerance for pruning
    double _vtol;                   // vibrational frequency matrix element tolerance for pruning
    double _extra_ener;             // maximal interpolation energy relative to the ground level
    double _extra_step;             // extrapolation logarithmic step
    double _ener_quant;             // energy descretization step
    int    _amom_max;               // angular momentum maximum

    // interpolation
    //
    double _pot_global_min;  // potential global minimum

    double _cstates_pow;     // extrapolation power value

    Slatec::Spline _cstates; // classical number/density of states relative to the potential minimum

    Slatec::Spline _qfactor; // quantum correction factor for number/density of states relative to the ground level

    void _set_qfactor ();

    void _set_states_base (Array<double>&, int =0) const;

    // estimates
    //
    std::vector<double> _mobility_parameter;

    Lapack::SymmetricMatrix _mobility_min;

    std::vector<double> _rotational_constant;
    
  public:
    
    MultiRotor(IO::KeyBufferStream&, const std::vector<Atom>&,  int = DENSITY) ;

    ~MultiRotor();

    int        internal_size ()      const { return _internal_rotation.size(); }

    int             symmetry (int i) const { return _internal_rotation[i].symmetry(); }

    double external_symmetry ()      const { return _external_symmetry; }

    double                  potential                (const std::vector<double>&, 
						      const std::map<int, int>& = std::map<int,int>()) const;
    Lapack::SymmetricMatrix mass                     (const std::vector<double>& angle)                const;
    Lapack::Vector          vibration                (const std::vector<double>& angle)                const;
    double                  external_rotation_factor (const std::vector<double>& angle)                const;
    double                  internal_rotation_factor (const std::vector<double>& angle)                const;
    double                         curvlinear_factor (const std::vector<double>& angle)                const;
    Lapack::Vector          frequencies              (const std::vector<double>& angle)                const;
    Lapack::SymmetricMatrix force_constant_matrix    (const std::vector<double>& angle)                const;
    Lapack::Vector          potential_gradient       (const std::vector<double>& angle)                const;
    void                    rotational_energy_levels ()                                                const;
 
    // statistical properties
    double quantum_weight           (double temperature)                         const;
    int    get_semiclassical_weight (double temperature, double& cw, double& pw) const;// classical & path integral
    void   quantum_states           (Array<double>&, double, int =0)             const;// relative to the ground

    // virtual functions
    double ground       () const;
    double states (double) const;// relative to the ground
    double weight (double) const;// relative to the ground
  };

  /********************************************************************************************
   *********** ABSTRACT CLASS REPRESENTING WELL, BARRIER, AND BIMOLECULAR FRAGMENT ************
   ********************************************************************************************/

  class Species {
    std::vector<Atom> _atom;
    std::string _name;
    int    _mode;

    Species ();

  protected:
    //
    double _ground;
    double _mass;
    double _print_min, _print_max, _print_step;

    void _print () const;

    Species (IO::KeyBufferStream&, const std::string&, int) ;
    Species (const std::string&, int);
    
  public:

    enum {CORE_WEIGHT,
	  HARMONIC_WEIGHT,
	  ANHARMONIC_WEIGHT,
	  HINDERED_ROTOR_WEIGHT,
	  ELECTRONIC_WEIGHT
    };
    
    virtual ~Species();

    virtual double states (double) const =0; // density or number of states of absolute energy
    virtual double weight (double) const =0; // weight relative to the ground

    double ground () const { return _ground; }
    virtual void shift_ground (double e) { _ground += e; }
    virtual double real_ground () const { return _ground; }
    virtual void init () {}
    
    double   mass () const ;

    const std::vector<Atom>& geometry () const { return _atom; }

    int mode () const { return _mode; }

    virtual double tunnel_weight (double) const;

    //void set_name (const std::string& n) { _name = n; }
    const std::string& name () const { return _name; }

    // radiational transitions
    virtual double   infrared_intensity (double, int) const;
    virtual double oscillator_frequency (int)         const;
    virtual int         oscillator_size ()            const;
  };

  inline Species::Species (const std::string& n, int m) : _name(n), _mode(m), _mass(-1.), _print_step(-1.), _ground(0.)
  {
    const char funame [] = "Model::Species::Species: ";
    
    if(m != DENSITY && m != NUMBER && m !=  NOSTATES) {
      std::cerr << funame << "wrong mode\n";
      throw Error::Logic();
    }
  }
  
  /***********************************************************************************************************
   ************************************* CRUDE MONTE-CARLO SAMPLING ******************************************
   ***********************************************************************************************************/

    // internal mode definition: bond distance, bond angle, or dihedral bond angle
    //
  class IntMod {
    //
    // atomic indices
    //
    std::vector<int> _atoms;

    std::set<int>     _pool;

    // no default constructor
    //
    IntMod ();
      
  public:
    //
    enum {
      DISTANCE = 2,
      ANGLE    = 3,
      DIHEDRAL = 4
    };

    IntMod (int molec_size, IO::KeyBufferStream&);
      
    double evaluate (Lapack::Vector          cart_pos,                 // cartesian coordinates of all atoms
		     const std::vector<int>& sign = std::vector<int>() // derivative signature
		     ) const;

    int type () const { return _atoms.size(); }
    
    static double increment;      
  };

  // explicitely sampled fluxional mode definition
  //
  class Fluxional : public IntMod {
    //
    double _span;

  public:

    Fluxional (int molec_size, IO::KeyBufferStream&);

    double span () const { return _span; }
  };

  class Constrain : public IntMod {
    //
    double _value;

  public:

    Constrain (int molec_size, IO::KeyBufferStream&);

    double value () const { return _value; }
  };
    
  extern "C" typedef void (*refp_t) (int& ifail, double& ener, const double* pos);

  extern "C" typedef void (*refw_t) (double& w, const double& t);

  class MonteCarlo : public Species {

    // atom mass square roots
    //
    std::vector<double> _mass_sqrt;

    double _total_mass_sqrt;

    double _atom_mass (int a) const { return _mass_sqrt[a] * _mass_sqrt[a]; }
    
    // explicitly sampled fluxional modes
    //
    std::vector<Fluxional> _fluxional;

    std::string _data_file;

    // reference energy
    //
    double _refen;

    // symmetry factor
    //
    double _symm_fac;

    // electronic energy levels
    //
    std::map<double, int> _elevel;
    
    // no quantum correction factor
    //
    bool _noqf;

    // no hessian data
    //
    bool _nohess;

    // no hessian curvlinear correction
    //
    bool _nocurv;

    // center-of-mass shift
    //
    bool _cmshift;

    // is transition state calculation
    //
    bool _ists;
    
    // non-fluxional modes frequencies
    //
    std::vector<double> _nm_freq;

    // minimal non-fluxional modes frequency 
    //
    static double nm_freq_min;
    
    // statistical weight prefactor including mass factors and quantum prefactor in local harmonic approximation
    //
    double _local_weight (double                  ener,       // energy
			  Lapack::Vector          cart_pos,   // cartesian coordinates    
			  Lapack::Vector          cart_grad,  // energy gradient in cartesian coordinates
			  Lapack::SymmetricMatrix cart_fc,    // cartesian force constant matrix
			  double                  temperature // temprature
			  ) const;

    // read data from the file
    //
    bool _read (std::istream&           from,       // data stream
	        double&                 ener,       // energy
		Lapack::Vector          cart_pos,   // cartesian coordinates    
		Lapack::Vector          cart_grad,  // energy gradient in cartesian coordinates
		Lapack::SymmetricMatrix cart_fc     // cartesian force constant matrix
		) const;

    // set reference energy to the minimal total energy including zero-point energy
    //
    void _set_reference_energy ();

    Lapack::SymmetricMatrix _inertia_matrix (Lapack::Vector pos) const;

    void _make_cm_shift (Lapack::Vector pos) const;
    
    // reference potential definition
    //
    class _RefPot {

      System::DynLib _lib;

      refp_t _pot;

      refw_t _weight;
      
    public:

      _RefPot () : _pot(0), _weight(0) {}

      void init (std::istream&);

      // energy as a function of fluxional coordinates
      //
      double operator() (const double*) const;

      operator bool () const { return _pot; }

      // reference statistical weight
      //
      double weight (double) const;

      class PotFail {};
    };

    // reference potential
    //
    _RefPot _ref_pot;

    // reference temperature
    //
    double  _ref_tem;
    
  public:
    //
    MonteCarlo (IO::KeyBufferStream&, const std::string&, int);
    
    ~MonteCarlo ();
    
    double states (double) const;
    
    double weight_with_error (double, double&) const;

    double weight (double temperature) const { double dtemp; return weight_with_error(temperature, dtemp); }

    int atom_size () const { return _mass_sqrt.size(); }
  };
  
  /***********************************************************************************************************
   ****************************** CRUDE MONTE-CARLO SAMPLING WITH DUMMY ATOMS ********************************
   ***********************************************************************************************************/

  class MonteCarloWithDummy : public Species {

    void _assert(ConstSlice<double> v) const;

    // dummy atoms indices
    //
    std::vector<int> _dummy_atom;

    // real atoms indices
    //
    std::vector<int> _real_atom;

    // atom sequence
    //
    std::vector<Atom> _atom_array;

    // vector operations with mass weighting
    //
    double _vdot(ConstSlice<double>, ConstSlice<double>) const;

    double _normalize(Slice<double>) const;
    
    void _orthogonalize (Slice<double>, ConstSlice<double>) const;
    
    void _orthogonalize (Lapack::Matrix m) const;

    // explicitly sampled fluxional modes
    //
    std::vector<Fluxional> _fluxional;

    std::vector<Constrain> _constrain;
    
    std::string _data_file;

    
    // statistical weight prefactor including mass factors and quantum prefactor in local harmonic approximation
    //
    double _local_weight (double                  ener,       // energy
			  Lapack::Vector          cart_pos,   // cartesian coordinates    
			  Lapack::Vector          cart_grad,  // energy gradient in cartesian coordinates
			  Lapack::SymmetricMatrix cart_fc,    // cartesian force constant matrix
			  double                  temperature // temprature
			  ) const;

  public:
    //
    MonteCarloWithDummy (IO::KeyBufferStream&, const std::string&, int);
    
    ~MonteCarloWithDummy ();
    
    double states (double) const;
    
    double weight (double) const;
  };
  
  /************************* RIGID ROTOR HARMONIC OSCILATOR MODEL ************************/

  class RRHO : public Species {
    //
    SharedPointer<Tunnel>                 _tunnel; // tunneling
    
    SharedPointer<Core>                     _core; // other degrees of freedom
    
    std::vector<SharedPointer<Rotor> >     _rotor; // hindered rotors
    
    std::vector<double>                _frequency; // real frequencies
    
    std::vector<double>                   _elevel; // electronic energy levels
    
    std::vector<int>                      _edegen; // electronic level degeneracies

    std::vector<int>                      _fdegen; // frequencies degeneracies
    
    Lapack::SymmetricMatrix               _anharm; // second order anharmonic energy level expansion

    double _sym_num;
    
    double _real_ground;

    // interpolation
    double           _emax; // interpolation energy maximum
    double           _nmax; // extrapolation power value
    Slatec::Spline _states;

    // radiative transitions
    std::vector<Slatec::Spline> _occ_num; // average occupation numbers for vibrational modes
    std::vector<double>     _occ_num_der; // occupation number derivatives (for extrapolation)
    std::vector<double>         _osc_int; // oscillator strength (infrared intensities)
    
    // graph perturbation theory
    Graph::Expansion  _graphex;
    void _init_graphex (std::istream&);

  public:
    RRHO (IO::KeyBufferStream&, const std::string&, int) ;
    ~RRHO ();

    double states (double) const; // density or number of states of absolute energy
    double weight (double) const; // weight relative to the ground

    double real_ground () const { return _real_ground; }
    void shift_ground (double e) { _ground += e; _real_ground += e; }

    double tunnel_weight (double) const;
    bool istunnel () const { return (bool)_tunnel; }

    // radiative transitions
    double   infrared_intensity (double, int) const;
    double oscillator_frequency (int)         const;
    int         oscillator_size ()            const;
  };

  /********************* READ DENSITY OF STATES FROM THE FILE AND INTERPOLATE ****************/
  
  class ReadSpecies : public Species {    
    Array<double>   _ener; // relative energy on the grid
    Array<double> _states; // density/number of states on the grid
    int _mode;

    Slatec::Spline _spline;
    double _emin, _emax;
    double _nmin, _amin;
    double _nmax, _amax;

    double _etol; // relative energy tolerance
    double _dtol; // absolute DoS    tolerance
    
  public:
    ReadSpecies (std::istream&, const std::string&, int) ;
    ~ReadSpecies ();

    double states (double) const;
    double weight (double) const;
  };
  
  
  /************************************* UNION OF SPECIES **************************************/

  class UnionSpecies : public Species {
    std::vector<SharedPointer<Species> > _species;
    typedef std::vector<SharedPointer<Species> >::const_iterator _Cit;

    double _real_ground;

    // radiational transitions
    std::vector<int> _osc_shift;
    std::vector<int> _osc_spec_index;

  public:

    UnionSpecies  (IO::KeyBufferStream&, const std::string&, int) ;
    ~UnionSpecies ();

    double states (double) const;
    double weight (double) const;

    void shift_ground (double);
    double real_ground () const { return _real_ground; }

    // radiative transitions
    double   infrared_intensity (double, int) const;
    double oscillator_frequency (int)         const;
    int         oscillator_size ()            const;
  };

  /****************************************** VARIATIONAL BARRIER MODEL ************************************/

  class VarBarrier : public Species {
    std::vector<SharedPointer<RRHO> > _rrho;
    SharedPointer<RRHO>               _outer;
    SharedPointer<Tunnel>             _tunnel;

    double _real_ground;

    double _ener_quant; // descretization step
    double       _emax; // extrapolation energy maximum
    double       _nmax; // extrapolation power value

    Array<double>  _stat_grid;
    Slatec::Spline _states;

    // two-transition-states model calculation method
    enum {STATISTICAL, DYNAMICAL};
    int _tts_method;

  public:
    VarBarrier (IO::KeyBufferStream& from, const std::string&) ;
    ~VarBarrier ();

    double states (double) const;
    double weight (double) const;

    double real_ground () const { return _real_ground; }
    void shift_ground (double e) { _ground += e; _real_ground += e; }

    double tunnel_weight (double) const;
  };

  

  /*********************************** ATOMIC FRAGMENT CLASS ***********************************/

  class AtomicSpecies : public Species {
    std::vector<double>                   _elevel; // electronic energy levels
    std::vector<int>                      _edegen; // electronic level degeneracies

  public:
    AtomicSpecies(IO::KeyBufferStream& from, const std::string&) ;
    ~AtomicSpecies ();

    double weight (double) const;
    double states (double) const;

    void shift_ground (double e);
  };

  /*********************************** ARRHENIUS  CLASS ***********************************/

  class Arrhenius : public Species {

    double _power;
    double _factor;
    double _ener;

    std::string _reactant;
    std::string _product;
    
    // interpolation
    double           _emax; // interpolation energy maximum
    double           _nmax; // extrapolation power value
    Slatec::Spline _states;

  public:
    Arrhenius(IO::KeyBufferStream& from, const std::string&) ;
    ~Arrhenius ();

    double weight (double) const;
    double states (double) const;

    void shift_ground (double e);

    void init ();
  };


  /********************************************************************************************
   ************************************* BIMOLECULAR CLASS ************************************
   ********************************************************************************************/

  // bimolecular products
  class Bimolecular {
    bool _dummy;
    std::vector<SharedPointer<Species> > _fragment;
    double _weight_fac;
    double _ground;
    std::string _name;

  public:
    Bimolecular (IO::KeyBufferStream&, const std::string&) ;
    ~Bimolecular ();

    bool     dummy       () const { return _dummy; }
    double  ground       () const;
    double  weight (double) const;
    void shift_ground (double);

    const std::string& fragment_name (int i) const { return _fragment[i]->name(); }
    double fragment_weight (int, double) const;

    //void set_name (const std::string& n) { _name = n; }
    const std::string& name () const { return _name; }
  };

  /********************************************************************************************
   *************************************** WELL ESCAPE ****************************************
   ********************************************************************************************/

  // abstract class for escape rate
  class Escape {
  public:
    virtual double rate (double) const =0;
    virtual void shift_ground (double) =0;
  };

  class ConstEscape : public Escape {
    double _rate;
  public:
    ConstEscape(IO::KeyBufferStream&) ;

    double rate (double) const {return _rate; }
    void shift_ground (double) {}
  };

  class FitEscape : public Escape {
    double _ground;
    Slatec::Spline _rate;

  public:
    FitEscape(IO::KeyBufferStream&) ;

    double rate (double) const;
    void shift_ground(double e) { _ground += e; }
  };

  /********************************************************************************************
   *********************** WELL = SPECIES + KERNEL + COLLISION + ESCAPE ***********************
   ********************************************************************************************/

  class ThermoChemistry {
  public:
    ThermoChemistry(std::istream&) ;

    void print (std::ostream&);
  };

  class Well {
    SharedPointer<Species>                  _species;
    std::vector<SharedPointer<Kernel> >     _kernel;
    std::vector<SharedPointer<Collision> >  _collision;
    SharedPointer<Escape>                   _escape;
    double                                  _extension;

  public:
    Well (IO::KeyBufferStream&, const std::string&) ;

    SharedPointer<Species>      species ()       { return _species; }
    ConstSharedPointer<Species> species () const { return _species; }

    ConstSharedPointer<Kernel>    kernel    (int i)  const { return _kernel[i]; }
    ConstSharedPointer<Collision> collision (int i)  const { return _collision[i]; }

    //void             set_name (const std::string&) ;
    const std::string&   name () const             ;
    double             ground () const             ;
    double             weight (double) const       ;
    double             states (double) const       ;
    double               mass () const             ;

    double        escape_rate (double ener) const { if(_escape) return _escape->rate(ener); else return 0.; }
    bool               escape () const { return (bool)_escape; } 

    void shift_ground (double) ;

    // radiation transitions
    double oscillator_frequency (int num) const;
    int         oscillator_size ()        const;

    // radiation down-transition probability
    double transition_probability (double ener, double temperature, int num) const;

    double dissociation_limit;
    double extension () const { return _extension; }
  }; // Well

  inline double Well::oscillator_frequency (int num) const 
  {
    const char funame [] = "Model::Well::oscillator_frequency: ";
    
    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->oscillator_frequency(num);
  }

  inline int Well::oscillator_size () const
  { 
    const char funame [] = "Model::Well::oscillator_size: ";
    
    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->oscillator_size();
  }

  inline void Well::shift_ground (double e) 
  {
    const char funame [] = "Model::Well::shift_ground: ";

    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    _species->shift_ground(e);

    if(_escape)
      _escape->shift_ground(e);
  }

  inline const std::string& Well::name () const 
  {
    const char funame [] = "Model::Well::name: ";

    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->name();
  }

  inline double Well::ground () const 
  {
    const char funame [] = "Model::Well::ground: ";

    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->ground();
  }

  inline double Well::mass () const 
  {
    const char funame [] = "Model::Well::mass: ";

    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->mass();
  }

  inline double Model::Well::weight (double t) const 
  {
    const char funame [] = "Model::Well::weight: ";

    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->weight(t);
  }

  inline double Well::states (double e) const 
  {
    const char funame [] = "Model::Well::states: ";

    if(!_species) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }

    return _species->states(e);
  }

  /********************************************************************************************
   **************************************** GLOBAL OBJECTS ************************************
   ********************************************************************************************/
  void init (IO::KeyBufferStream& from) ;
  bool isinit ();

  bool no_run ();
  
  int            well_size ();
  int     bimolecular_size ();
  int   inner_barrier_size ();
  int   outer_barrier_size ();

  int buffer_size ();
  double buffer_fraction (int);
  ConstSharedPointer<Collision>  collision (int);
  ConstSharedPointer<Kernel> default_kernel (int);

  const Well&                            well (int w);
  const Bimolecular&              bimolecular (int p);
  const Species&                inner_barrier (int b);
  const Species&                outer_barrier (int b);
  const std::pair<int, int>&    inner_connect (int b);
  const std::pair<int, int>&    outer_connect (int b);

  double  maximum_barrier_height ();

  // energy shift
  extern std::string reactant; // bimolecular species to use as an energy reference
  double energy_shift ();

  /********************************************************************************
   ****************************** TIME EVOLUTION **********************************
   ********************************************************************************/

  class TimeEvolution {
    double _excess;
    double _start;
    double _finish;
    double _step;
    int    _size;
    double _temperature;

    mutable int    _reactant;
    std::string    _reactant_name;

  public:
    TimeEvolution (IO::KeyBufferStream&) ;
    ~TimeEvolution () { out.close(); }

    void set_reactant () const;

    double  start () const { return _start; }
    double finish () const { return _finish; }
    double   step () const { return _step; }
    int      size () const { return _size; }
    int  reactant () const { if(_reactant < 0) set_reactant(); return _reactant; }
    double excess_reactant_concentration () const { return _excess; }
    double temperature () const { return _temperature; }

    std::ofstream out;
  };

  extern SharedPointer<TimeEvolution> time_evolution;

}// namespace Model

#endif
