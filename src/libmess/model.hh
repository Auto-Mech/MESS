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
#include<list>
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

  extern int log_precision;

  extern int out_precision;

  extern bool use_short_names;

  int name_size_max ();
  
  // minimal interatomic distance
  //
  extern double atom_dist_min;

  // maximum energy to be used
  //
  double  energy_limit () ;
  void set_energy_limit (double);
  bool is_energy_limit ();

  int escape_size       ();                 // number escape channels
  
  const std::pair<int, int>& escape_channel (int); // escape channel: <well, well escape index>

  std::string escape_name (int); // escape channel name

  void shift_cm_to_zero(std::vector<Atom>&);
  Lapack::SymmetricMatrix inertia_moment_matrix(const std::vector<Atom>&);
  void read_geometry (IO::KeyBufferStream&, std::vector<Atom>&, int =ANGSTROM);

  bool is_well (const std::string&);

  int well_by_name (const std::string&);
  
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
  SharedPointer<Species>     new_species        (IO::KeyBufferStream&, const std::string&, int, std::pair<bool, double> = std::make_pair(false, 0.));
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
    //
    double _frequency_factor;
    
    double _epsilon;

    double _omega_22_star (double) const;
    
  public:
    //
    LennardJonesCollision (IO::KeyBufferStream&);
    
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
      NOTRUN  = 4,  // no truncation even so the transition probability is negative
      DOWN = 8
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
    //
    double        _wtol;// statistical weight tolerance
    
    static double _action_max; // maximum 

    double        _freq; // imaginary frequency
    
  protected:
    //
    double        _cutoff; // cutoff energy
    
    double          _efac;// entanglement correction
    
    void read_freq(IO::KeyBufferStream&);

    void assert_freq() const;

    double frequency () const { return _freq; }
    
    Tunnel(IO::KeyBufferStream&);

    Tunnel (double c, double f) : _cutoff(c), _freq(f), _efac(0.) {}
    
  public:
    //
    virtual ~Tunnel ();

    double  cutoff () const { return _cutoff; }

    void set_cutoff (double c) { _cutoff = c; }

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
    ReadTunnel(IO::KeyBufferStream&);
    
    ~ReadTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /**************************************************************************************
   ************************** PARABOLIC BARRIER TUNNELING *******************************
   **************************************************************************************/

  class HarmonicTunnel: public Tunnel {
    //
  public:
    //
    HarmonicTunnel(IO::KeyBufferStream&);

    HarmonicTunnel (double c, double f) : Tunnel(c, f) {}
    
    ~HarmonicTunnel ();

    double action (double, int =0) const;
  };

  /**************************************************************************************
   ************************** POWER EXPANSION TUNNELING *******************************
   **************************************************************************************/

  class ExpTunnel: public Tunnel {

    std::map<int, double> _expansion;

  public:
    //
    ExpTunnel(IO::KeyBufferStream&);

    ~ExpTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /**************************************************************************************
   ************************** ECKART BARRIER TUNNELING ********************************
   **************************************************************************************/

  class EckartTunnel: public Tunnel {
    //
    std::vector<double> _depth;
    
    double             _factor;

  public:
    //
    EckartTunnel(IO::KeyBufferStream&);
    
    ~EckartTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /**************************************************************************************
   ************************** QUARTIC BARRIER TUNNELING ********************************
   **************************************************************************************/

  class QuarticTunnel: public Tunnel {
    //
    double _vmin;// minimal well depth
    
    double   _v3;// x^3 power expansion coefficient
    
    double   _v4;// x^4 power expansion coefficient

    double _potential (double x) { return x * x * (0.5  + _v3 * x + _v4 * x * x); }

    Slatec::Spline _action; // semiclassical action for energies below barrier 

    class XratioSearch : public Math::NewtonRaphsonSearch {
      //
      double _vratio;

    public:
      //
      XratioSearch(double v, double t) : _vratio(v) { tol = t; }
      
      double operator() (double, int) const;
    };

  public:
    //
    QuarticTunnel(IO::KeyBufferStream&);
    
    ~QuarticTunnel ();

    double action (double, int =0) const; // semiclassical action
  };

  /********************************************************************************************
   ******************************** INTERNAL ROTATION DEFINITION ******************************
   ********************************************************************************************/

  class InternalRotationDef {
    //
    std::set<int>      _group; // moving group definition
    
    std::pair<int, int> _axis; // moving direction definition
    
    int             _symmetry; // internal rotation symmetry
    
    int                 _imax; // maximal atomic index in internal roation definition

    bool _isinit;

  public:
    //
    void init(IO::KeyBufferStream&);
    
    InternalRotationDef()                           : _isinit(false), _symmetry(1), _axis(-1, -1) {}
    
    explicit InternalRotationDef(int s)             : _isinit(false), _symmetry(s), _axis(-1, -1) {}
    
    explicit InternalRotationDef (IO::KeyBufferStream& from) : _isinit(false), _symmetry(1), _axis(-1, -1) { init(from); }
    
    virtual ~InternalRotationDef ();

    int symmetry () const { return _symmetry; }

    std::vector<Atom>            rotate (const std::vector<Atom>& atom, double angle) const;
    
    std::vector<D3::Vector> normal_mode (const std::vector<Atom>& atom, Lapack::Vector* =0) const;
  };

  /********************************************************************************************
   ************************* HINDERED ROTOR/UMBRELLA MODE CLASS FAMILY ************************
   ********************************************************************************************/

  class Rotor {

  protected:
    //
    // controls

     // maximum dimension of the Hamiltonian
    //
    int      _ham_size_max;

     // minimal dimension of the Hamiltonian
    //
    int      _ham_size_min;

     // angular discretization size
    //
    int _grid_size;

     // thermal exponent power maximum
    //
    double  _therm_pow_max;

     // 3D structure
    //
    std::vector<Atom> _atom;

    Rotor ();
    
    Rotor (IO::KeyBufferStream&, const std::vector<Atom>&) ;

  public:
    
    virtual ~Rotor();

    int angular_grid_size () const { return _grid_size; }
    
    // set maximal energy relative to ground
    //
    virtual void   set (double) { std::cerr << "Model::Rotor::set: not defined\n"; throw Error::Init(); }

    // ground energy
    //
    virtual double ground        () const { std::cerr << "Model::Rotor::ground: not defined\n"; throw Error::Init(); }
    
    // energy level relative to the ground
    //
    virtual double energy_level  (int)  const { std::cerr << "Model::Rotor::energy_level: not defined\n"; throw Error::Init(); }

    // number of energy levels
    //
    virtual int    level_size    ()       const { std::cerr << "Model::Rotor::level_size: not defined\n"; throw Error::Init(); }
    
    // statistical weight relative to the ground
    //
    virtual double weight        (double) const { std::cerr << "Model::Rotor::weight: not defined\n"; throw Error::Init(); }

    virtual void convolute(Array<double>&, double) const;
  };

  // rotational constant for free and hindered rotors
  //
  class RotorBase : public Rotor, public InternalRotationDef  {
    //
    double  _rotational_constant;    

  public:

    RotorBase (double r, int s) : InternalRotationDef(s), _rotational_constant(r) {}
    
    RotorBase (IO::KeyBufferStream&, const std::vector<Atom>&);

    virtual ~RotorBase ();

    double rotational_constant () const { return _rotational_constant; }
  };

  /********************************************************************************************
   ****************************************** FREE ROTOR **************************************
   ********************************************************************************************/

  class FreeRotor : public RotorBase {
    //
    int _level_size;

  public:
    //
    FreeRotor (IO::KeyBufferStream&, const std::vector<Atom>&);
    
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
    //
    // hindered rotor potential
    //
    class _Potential {
      //
    public:
      //
      enum {AKIMA_SPLINE, C_SPLINE, FOURIER_EXPANSION, FOURIER_FIT };

    private:
      //
      // key-to-type map
      //
      class _KeyType : private std::map<std::string, int> {
	//
      public:
	//
	_KeyType () {
	  //
	  (*this)["akima_spline"     ] = AKIMA_SPLINE;
	  
	  (*this)["c_spline"         ] = C_SPLINE;
	  
	  (*this)["fourier_expansion"] = FOURIER_EXPANSION;
	  
	  (*this)["fourier_fit"      ] = FOURIER_FIT;
	}
	
	void print_keys (std::ostream&) const;

	int find (const std::string&) const;
      };

      static const _KeyType _key_type;
      
      int _type; // potential type

      int _symm; // symmetry
      
      Math::Spline  _spline;

      Array<double> _fourier;

      bool _isinit;

      int _print_pot_size;

      double _spline_value (double, int =0) const;

      double _fourier_value (double, int =0) const;
      
      //...

      // no copies
      //
      _Potential            (const _Potential&);
      
      _Potential& operator= (const _Potential&);

    public:
      //
      _Potential () :_isinit(false), _print_pot_size(0) {}

      void init (IO::KeyBufferStream&, int);
      
      _Potential (IO::KeyBufferStream& from, int s) : _print_pot_size(0) { init(from, s); }

      _Potential(const Array<double>& f, int s)
	//
	: _isinit(true), _type(FOURIER_EXPANSION), _symm(s), _fourier(f), _print_pot_size(0) {}

      bool isinit () const { return _isinit; }

      int type () const { return _type; }
      
      int symmetry () const { return _symm; }
      
      // value and derivatives
      //
      double operator() (double, int =0) const;

      // Fourier expansion
      //
      int fourier_size () const { return _fourier.size(); }
      
      double fourier (int i) const { return _fourier[i]; }
    };

    _Potential                _pot;

    static int _rotation_matrix_element (int, int, int, double&);
    
    double                    _ground; // ground level energy
    std::vector<double> _energy_level; // energy levels relative to the ground
    Array<double>           _pot_four; // potential fourier expansion

    std::vector<double>  _pot_grid; // potential-on-the-grid
    std::vector<double> _freq_grid; // frequency-on-the-grid
    double              _grid_step; // angular discretization step
    double                _pot_max; // global potential energy maximum
    double                _pot_min; // global potential energy minimum
    double               _freq_min; // frequency at minimum
    double               _freq_max; // maximal frequency
    double            _harm_ground; // ground energy in harmonic approximation
    
    double               _ener_max; // maximal level energy
    
    void _fourier_space_energy_levels (int);
    Lapack::Vector _real_space_energy_levels () const;
    
    void _read (IO::KeyBufferStream&);
    void _init ();

    int _flags;

    //Akima vs. C-spline
    //
    bool _use_akima;

    // Slatec vs GSL
    //
    bool _use_slatec;
    
  public:
    //
    enum {NOPRINT = 1};
    
    static int weight_output_temperature_step; // temperature step for statistical weight output
    static int weight_output_temperature_max;  // temperature maximum for statistical weight output
    static int weight_output_temperature_min;  // temperature maximum for statistical weight output

    // switch from quantum to semiclassical partition function calculation
    //
    static double weight_thres;

    // minimal sampling angle step, degrees
    //
    static double min_ang_step;

    HinderedRotor (IO::KeyBufferStream&, const std::vector<Atom>&, int f = 0);
    
    HinderedRotor (const RotorBase& r, const std::vector<double>& p, int f = 0);
    
    HinderedRotor (const std::vector<double>& p, double r, int s, int f = 0);
    
    ~HinderedRotor ();

    double                        potential           (double, int =0) const;
    int         semiclassical_states_number                   (double) const;
    double                   quantum_weight                   (double) const;
    int            get_semiclassical_weight (double, double&, double&) const;

    // integrate local energy distribution over angle
    //
    void integrate (Array<double>&, double) const;
    
    double potential_minimum () { return _pot_min; }

    // virtual functions
    //
    void set (double);

    double ground       ()       const;
    
    double energy_level (int)    const;
    
    int    level_size   ()       const;
    
    double weight       (double) const;
    
    void convolute(Array<double>&, double) const;
  };

  /*********************************************************************************************
   ******************* HINDERED ROTOR BUNDLE WITH ADIABATIC FREQUENCIES  ***********************
   ********************************************************************************************/

  class HinderedRotorBundle {

    double _ground; // ground energy

    double _ener_step; // discretization step
    
    Array<double> _qstates; //   quantum # of states per bin (low energies)

    Array<double> _cstates; // classical # of states per bin (high energies)

    int _flags;

  public:

    enum { NOPRINT = 1 };
    
    HinderedRotorBundle (IO::KeyBufferStream&, const std::vector<Atom>&);

    double ground ()       const { return _ground; }

    double energy_step () const { return _ener_step; }
    
    int size () const { return _cstates.size(); }

    double states (int e) const { if(e < _qstates.size()) return _qstates[e]; return _cstates[e]; }
    
    double weight (double) const;
  };
    
  /********************************************************************************************
   ***************************************** UMBRELLA MODE ************************************
   ********************************************************************************************/

  
  class Umbrella : public Rotor {

    double                    _ground;
    
    std::vector<double> _energy_level;

    double                      _mass;
    
    Lapack::Vector          _pot_coef;

    // discretization
    //
    std::vector<double>  _pot_grid; // potential-on-the-grid
    
    std::vector<double> _freq_grid; // frequency-on-the-grid
    
    double                  _astep; // discretization step
    
    double                _pot_min; // global potential energy minimum
    
    static double _integral(int p, int n);// \int_0^1 dx x^p cos(n\pi x)
    
    void _fourier_space_energy_levels(int) ;

  public:
    Umbrella (IO::KeyBufferStream&, const std::vector<Atom>&);
    
    ~Umbrella ();

    double quantum_weight (double) const;
    
    int get_semiclassical_weight (double, double&, double&) const;
    
    double potential(double x, int der =0) const ;

    // virtual functions
    //
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
    //
    int _mode;
    
    Core ();

  protected:
    //
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
      //
      std::cerr << funame << "wrong mode\n";
      
      throw Error::Logic();
    }
  }

  /**************************************************************************************
   ************************************* DUMMY CORE *************************************
   **************************************************************************************/

  class DummyCore : public Core {

    double _symmetry;

  public:

    DummyCore (IO::KeyBufferStream& from);

    ~DummyCore () {}

    double ground () const { return 0.; }

    double weight (double) const { return _symmetry; }

    double states (double) const;
  };
  
  /**************************************************************************************
   ******************************** PHASE SPACE THEORY **********************************
   **************************************************************************************/

  class PhaseSpaceTheory : public Core {
    //
    double _states_factor;
    
    double _weight_factor;
    
    double _power;

    // TST levels
    //
    enum {T_LEVEL, E_LEVEL, EJ_LEVEL};
    
  public:
    //
    PhaseSpaceTheory (IO::KeyBufferStream& from);
    
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
   ******************************** INTERNAL ROTATION SAMPLING CLASS **************************
   ********************************************************************************************/
  
    class MultiRotorWithBumps: public Core {
      //
      // internal rotations description
      //
      std::vector<InternalRotationDef> _internal_rotation;

      
      // hindered rotors potential
      //
      std::vector<std::map<int, double> > _hr_pot;

      // hydrogen bond potential
      //
      // ...

      // repulsion potential
      //
      // ...

      
    };

  /********************************************************************************************
   ********************************* COUPLED INTERNAL ROTORS MODEL ****************************
   ********************************************************************************************/

  class MultiRotor: public Core {
    //
    // Internal rotation for multirotor
    //
    class InternalRotation : public InternalRotationDef {
      int                _msize; // mass fourier expansion size;
      int                _psize; // potential fourier expansion size
      int                _wsize; // angular grid size
      int                 _qmin; // minimum quantum state size
      int                 _qmax; // maximum quantum state size

    public:
      //
      InternalRotation (IO::KeyBufferStream&);
    
      ~InternalRotation ();

      int      mass_fourier_size () const { return    _msize; }
      int potential_fourier_size () const { return    _psize; }
      int   weight_sampling_size () const { return    _wsize; }
      int       quantum_size_min () const { return     _qmin; }
      int       quantum_size_max () const { return     _qmax; }
    };

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
    //
    std::vector<Atom> _atom;
    
    std::string       _name;
    
    int               _mode;

    Species ();

  protected:
    //
    double _ground;
    
    double _mass;
    
    double _print_min, _print_max, _print_step;

    void _print () const;

    Species (IO::KeyBufferStream&, const std::string&, int);
    
    Species (const std::string&, int);
    
  public:

    enum {CORE_WEIGHT,
	  //
	  HARMONIC_WEIGHT,
	  //
	  ANHARMONIC_WEIGHT,
	  //
	  HINDERED_ROTOR_WEIGHT,
	  //
	  ELECTRONIC_WEIGHT
    };
    
    virtual ~Species();
    
    // density or number of states of absolute energy
    //
    virtual double states (double) const =0;

    // weight relative to the ground
    //
    virtual double weight (double) const =0;

    // free energy
    //
    double free_energy (double temperature) const { return -temperature * std::log(weight(temperature)); }
    
    // (relative) temperature differentiation step
    //
    static double temp_diff_step;
    
    // thermal parameters: energy, entropy, and thermal capacity
    //
    void esc_parameters (double temperature, double& e, double& s, double& c) const;

    // ground state energy
    //
    double ground () const { return _ground; }
    
    virtual void shift_ground (double e) { _ground += e; }
    
    virtual double real_ground () const { return _ground; }
    
    double thermal_energy (double temperature) const { double e, s, c; esc_parameters(temperature, e, s, c); return ground() + e; }
    
    virtual void init () {}
    
    double   mass () const ;

    const std::vector<Atom>& geometry () const { return _atom; }

    int mode () const { return _mode; }

    virtual double tunnel_weight (double) const;

    const std::string& name () const { return _name; }

    std::string short_name () const;

    // radiational transitions
    //
    virtual double   infrared_intensity (double, int) const;
    
    virtual double oscillator_frequency (int)         const;
    
    virtual int         oscillator_size ()            const;
  };

  inline Species::Species (const std::string& n, int m) : _name(n), _mode(m), _mass(-1.), _print_step(-1.), _ground(0.)
  {
    const char funame [] = "Model::Species::Species: ";
    
    if(m != DENSITY && m != NUMBER && m !=  NOSTATES) {
      //
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

    // correction factor
    //
    double _corr_fac;

    // electronic energy levels
    //
    std::map<double, int> _elevel;
    
    // no hessian data
    //
    bool _nohess;

    // no hessian curvlinear correction
    //
    bool _nocurv;

    // no Pitzer-Gwinn density of states correction
    //
    bool _nopg;

    // center-of-mass shift
    //
    bool _cmshift;

    // is transition state calculation
    //
    bool _ists;

    // sampling data
    //
    struct _Sampling {
      //
      double                  ener;      // potential energy
      
      Lapack::Vector          cart_pos;  // cartesian coordinates
      
      Lapack::Vector          cart_grad; // energy gradient in cartesian coordinates
      
      Lapack::SymmetricMatrix cart_fc;   // force constant matrix

      _Sampling (int s, std::istream& from, int flags =0);

      class End {};
    };

    std::vector<_Sampling> _sampling_data;
    
    // statistical weight prefactor including mass factors and quantum prefactor in local harmonic approximation
    //
    double _local_weight (const _Sampling&        s,
			  double                  temperature, // temprature
			  Array<double>&          local_dos,   // local density of states
			  int&                    flux_shift,  // fluxional modes energy shift
			  int                     flags = 0
			  ) const;

    // read data from the file
    //
    static bool _read (std::istream&           from,       // data stream
		       double&                 ener,       // energy
		       Lapack::Vector          cart_pos,   // cartesian coordinates    
		       Lapack::Vector          cart_grad,  // energy gradient in cartesian coordinates
		       Lapack::SymmetricMatrix cart_fc,    // cartesian force constant matrix
		       int                     flags = 0
		       );

    // internal motion basis
    //
    Lapack::Matrix _internal_basis (Lapack::Vector pos) const;
    
    // set reference energy to the minimal total energy including zero-point energy
    //
    int _set_reference_energy ();

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

    // internal rotation definition for calculations without constrain optimization
    //
    std::vector<InternalRotationDef> _internal_rotation;

    // reference energy
    //
    mutable double _refen;

    // DOS reference energy
    //
    mutable double _dosen;

    // minimal potential energy
    //
    mutable double _potmin;
    
    // reference point (true) frequencies
    //
    mutable std::vector<double> _ref_freq;

    // reference external rotation factor (square root of inertia moments product)
    //
    mutable double _rfac;

    // reference configuration symmetry
    //
    int _rsymm;

    Slatec::Spline _ref_states;
    
    // reference non-fluxional modes frequencies
    //
    mutable std::vector<double> _nm_freq;

    mutable std::vector<Atom> _atom;

    // descretization energy bin
    //
    double _ener_quant;

    // upper energy limit
    //
    double _ener_max;

    // shift density of states
    //
    static void _shift (int s, Array<double>& dos);

    Slatec::Spline _states;

    // inverse Laplace transform log number of states
    //
    Slatec::Spline _ilt_log;

    // temperature interpolation grid parameters
    //
    double _ilt_min, _ilt_max;

    int _ilt_num;

    // use inverse laplace transform states
    //
    bool _use_ilt;

    mutable bool _deep_tunnel;
    
  public:
    //
    enum { REF = 1, DOS = 2, NOHESS = 4 };
    
    MonteCarlo (IO::KeyBufferStream&, const std::string&, int, std::pair<bool, double> = std::make_pair(false, 0.));
    
    ~MonteCarlo ();
    
    double states (double) const;
    
    double weight_with_error (double, double&) const;

    double weight (double temperature) const { double dtemp; return weight_with_error(temperature, dtemp); }

    int atom_size () const { return _mass_sqrt.size(); }

    int fluxional_size() const { if (_fluxional.size()) return _fluxional.size(); return _internal_rotation.size(); }
    
    // minimal non-fluxional modes frequency 
    //
    double nm_freq_min;

    // high frequency threshold x = freq / 2 / T.
    //
    double high_freq_thres;

    // low frequency threshold x = freq / 2 / T.
    //
    double low_freq_thres;
  
    // maximal exponent argument
    //
    double exp_arg_max;

    // deep tunneling threshold
    //
    double deep_tunnel_thres;
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
    MonteCarloWithDummy (IO::KeyBufferStream&, const std::string&, int, std::pair<bool, double> = std::make_pair(false, 0.));
    
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

    std::vector<ConstSharedPointer<HinderedRotorBundle> > _hrb;
    
    std::vector<double>                _frequency; // real frequencies
    
    std::vector<double>                   _elevel; // electronic energy levels
    
    std::vector<int>                      _edegen; // electronic level degeneracies

    std::vector<int>                      _fdegen; // frequencies degeneracies
    
    Lapack::SymmetricMatrix               _anharm; // second order anharmonic energy level expansion

    double _sym_num;
    
    double _real_ground;

    double _ground_shift;

    // interpolation
    //
    double     _ener_quant; // interpolation energy step
    double           _emax; // interpolation energy maximum
    double           _nmax; // extrapolation power value

    Array<double>  _stat_grid;
    Slatec::Spline _states;

    // radiative transitions
    //
    std::vector<Slatec::Spline> _occ_num; // average occupation numbers for vibrational modes
    std::vector<double>     _occ_num_der; // occupation number derivatives (for extrapolation)
    std::vector<double>         _osc_int; // oscillator strength (infrared intensities)
    
    // graph perturbation theory
    //
    Graph::Expansion  _graphex;
    void _init_graphex (std::istream&);

  public:
    //
    RRHO (IO::KeyBufferStream&, const std::string&, int, std::pair<bool, double> = std::make_pair(false, 0.));
    
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
    //
    // relative energy on the grid
    //
    Array<double>   _ener;

    // density/number of states on the grid
    //
    Array<double> _states;
    
    //int           _mode;

    Slatec::Spline _spline;
    double _emin, _emax;
    double _nmin, _amin;
    double _nmax, _amax;

    double _etol; // relative energy tolerance
    double _dtol; // absolute DoS    tolerance
    
  public:
    //
    ReadSpecies (std::istream&, const std::string&, int, std::pair<bool, double> = std::make_pair(false, 0.));
    
    ~ReadSpecies ();

    double states (double) const;
    double weight (double) const;
  };
  
  
  /************************************* UNION OF SPECIES **************************************/

  class UnionSpecies : public Species {
    //
    std::vector<SharedPointer<Species> > _species;
    
    typedef std::vector<SharedPointer<Species> >::const_iterator _Cit;

    double _real_ground;

    // radiational transitions
    //
    std::vector<int> _osc_shift;
    std::vector<int> _osc_spec_index;

    void _read (IO::KeyBufferStream&, const std::string&, int, std::pair<bool, double>);

    void _set ();

  public:
    //
    UnionSpecies  (IO::KeyBufferStream&, const std::string&, int, std::pair<bool, double> = std::make_pair(false, 0.));
    
    UnionSpecies  (const std::vector<SharedPointer<Species> >&, const std::string&, int);
    
    ~UnionSpecies ();

    double states (double) const;
    double weight (double) const;

    void shift_ground (double);
    double real_ground () const { return _real_ground; }

    // radiative transitions
    //
    double   infrared_intensity (double, int) const;
    double oscillator_frequency (int)         const;
    int         oscillator_size ()            const;
  };

  /****************************************** VARIATIONAL BARRIER MODEL ************************************/

  class VarBarrier : public Species {
    //
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
    //
    enum {STATISTICAL, DYNAMICAL};
    
    int _tts_method;

  public:
    //
    VarBarrier (IO::KeyBufferStream& from, const std::string&, std::pair<bool, double>);
    
    ~VarBarrier ();

    double states (double) const;
    double weight (double) const;

    double real_ground () const { return _real_ground; }
    void shift_ground (double e) { _ground += e; _real_ground += e; }

    double tunnel_weight (double) const;
  };

  

  /*********************************** ATOMIC FRAGMENT CLASS ***********************************/

  class AtomicSpecies : public Species {
    //
    std::vector<double>                   _elevel; // electronic energy levels
    
    std::vector<int>                      _edegen; // electronic level degeneracies

  public:
    //
    AtomicSpecies(IO::KeyBufferStream& from, const std::string&);
    
    ~AtomicSpecies ();

    double weight (double) const;
    double states (double) const;

    void shift_ground (double e);
  };

  /*********************************** ARRHENIUS  CLASS ***********************************/

  class Arrhenius : public Species {
    //
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
    //
    Arrhenius(IO::KeyBufferStream& from, const std::string&);
    
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
    //
    bool _dummy;
    
    std::vector<SharedPointer<Species> > _fragment;
    
    double _weight_fac;
    
    double _ground;
    
    std::string _name;

  public:
    //
    Bimolecular (IO::KeyBufferStream&, const std::string&);
    
    ~Bimolecular ();

    bool     dummy       () const { return _dummy; }
    double  ground       () const;
    double  weight (double) const;
    void shift_ground (double);

    const std::string& fragment_name (int i) const { return _fragment[i]->name(); }
    
    double fragment_weight (int, double) const;

    const std::string& name () const { return _name; }

    std::string short_name () const;

    int mode (int f) const {return _fragment[f]->mode(); }
    
    double states (int f, double e) const { return _fragment[f]->states(e); }
  };

  /********************************************************************************************
   *************************************** WELL ESCAPE ****************************************
   ********************************************************************************************/

  // abstract class for escape rate
  //
  class Escape {
    //
    void _assert () const;
    
    ConstSharedPointer<Species> _spec;

  protected:

    std::string _name;
    
  public:
    //
    virtual double rate (double) const =0;

    double states (double ener) const { _assert(); return rate(ener) * _spec->states(ener); }
    
    virtual void shift_ground (double) =0;

    void set_spec (ConstSharedPointer<Species> s) { _spec = s; }

    const std::string& name () const { return _name; }
  };

  class ConstEscape : public Escape {
    //
    double _rate;
    
  public:
    //
    ConstEscape(IO::KeyBufferStream&);

    double rate (double) const {return _rate; }
    
    void shift_ground (double) {}
  };

  class FitEscape : public Escape {
    //
    double _ground;
    
    Slatec::Spline _rate;

  public:
    //
    FitEscape(IO::KeyBufferStream&) ;

    double rate (double) const;
    
    void shift_ground(double e) { _ground += e; }
  };

  /********************************************************************************************
   *********************** WELL = SPECIES + KERNEL + COLLISION + ESCAPE ***********************
   ********************************************************************************************/

  class ThermoChemistry {
    //
  public:
    //
    ThermoChemistry(std::istream&) ;

    void print (std::ostream&);
  };

  class Well {
    //
    SharedPointer<Species>                 _species;
    
    std::vector<SharedPointer<Kernel> >     _kernel;
    
    std::vector<SharedPointer<Collision> >  _collision;
    
    std::vector<SharedPointer<Escape> >     _escape;

    void _assert_spec() const;

  public:
    //
    Well (IO::KeyBufferStream&, const std::string&);

    Well (const std::vector<Well>&);

    const Kernel&       kernel (int i)   const { return *_kernel[i]; }
    
    const Collision& collision (int i)   const { return *_collision[i]; }

    std::string   short_name ()         const { _assert_spec(); return _species->short_name(); }
    
    const std::string&   name ()         const { _assert_spec(); return _species->name(); }
    
    double             ground ()         const { _assert_spec(); return _species->ground(); }
    
    double             weight (double t) const { _assert_spec(); return _species->weight(t); }
    
    double             states (double e) const { _assert_spec(); return _species->states(e); }
    
    double               mass ()         const { _assert_spec(); return _species->mass(); }
    
    void  esc_parameters (double t, double& e, double& s, double& c) const { _assert_spec(); return _species->esc_parameters(t, e, s, c); }
    
    double thermal_energy (double t) const { _assert_spec(); return _species->thermal_energy(t); }
    
    int  escape_size () const { return _escape.size(); }

    const std::string& escape_name (int i) const { return _escape[i]->name(); }
    
    double escape_rate (double) const;

    double escape_rate (double, int) const;

    void shift_ground (double);

    // radiation transitions
    //
    double oscillator_frequency (int num) const { _assert_spec(); return _species->oscillator_frequency(num); }
    
    int         oscillator_size ()        const { _assert_spec(); return _species->oscillator_size(); }

    // radiation down-transition probability
    //
    double transition_probability (double ener, double temperature, int num) const;

    double dissociation_limit;

    double well_ext_cap;
  };

  inline void Well::_assert_spec() const
  {
    const char funame [] = "Model::Well::_assert_spec: ";
	
    if(!_species) {
      //
      std::cerr << funame << "not initialized\n";
      
      throw Error::Init();
    }
  }

  inline void Well::shift_ground (double e) 
  {
    const char funame [] = "Model::Well::shift_ground: ";

    _assert_spec();
    
    _species->shift_ground(e);

    for(int i = 0; i < _escape.size(); ++i)
      //
      _escape[i]->shift_ground(e);
  }

  /********************************************************************************************
   **************************************** GLOBAL OBJECTS ************************************
   ********************************************************************************************/

  extern std::set<std::string> well_exclude_group;

  extern std::list<std::string>  lump_scheme;

  // main initialization
  //
  void init (IO::KeyBufferStream& from) ;
  bool isinit ();

  // reset the model on the base of the new lumping scheme
  //
  void reset (const std::vector<std::set<int> >& = std::vector<std::set<int> >());

  // print the wells, barriers, etc.
  //
  void pes_print ();
  
  void pf_print ();
  
  bool no_run ();
  
  int    buffer_size ();
  double buffer_fraction (int);
  
  ConstSharedPointer<Collision>  collision      (int);
  ConstSharedPointer<Kernel>     default_kernel (int);

  int            well_size ();
  int     bimolecular_size ();
  int   inner_barrier_size ();
  int   outer_barrier_size ();

  inline void assert_well (int w)
  {
    const char funame [] = "Model::assert_well: ";
    
    if(w >= 0 && w < well_size())
      //
      return;

    std::cerr << funame << "well index out of range: " << w << "\n";

    throw Error::Range();
  }
  
  inline void assert_inner (int b)
  {
    const char funame [] = "Model::assert_inner: ";
    
    if(b >= 0 && b < inner_barrier_size())
      //
      return;

    std::cerr << funame << "inner barrier index out of range: " << b << "\n";

    throw Error::Range();
  }
  
  inline void assert_outer (int b)
  {
    const char funame [] = "Model::assert_outer: ";
    
    if(b >= 0 && b < outer_barrier_size())
      //
      return;

    std::cerr << funame << "outer barrier index out of range: " << b << "\n";

    throw Error::Range();
  }
  
  inline void assert_bimolecular (int p)
  {
    const char funame [] = "Model::assert_bimolecular: ";
    
    if(p >= 0 && p < bimolecular_size())
      //
      return;

    std::cerr << funame << "bimolecular index out of range: " << p << "\n";

    throw Error::Range();
  }
  
  const Well&                            well (int w);
  const Bimolecular&              bimolecular (int p);
  const Species&                inner_barrier (int b);
  const Species&                outer_barrier (int b);
  const std::pair<int, int>&    inner_connect (int b);
  const std::pair<int, int>&    outer_connect (int b);

  enum {UNKNOWN, INNER, OUTER, WELL, BIMOLECULAR};

  typedef std::pair<int, int> spec_t;
  
  // well/bimolecular/inner/outer name output
  //
  inline std::string species_name (spec_t p)
  {
    const char funame [] = "Model::species_name: ";
    
    const int& type = p.first;

    const int& index = p.second;
    
    switch(type) {
      //
    case WELL:
      //
      return well(index).short_name();

    case BIMOLECULAR:
      //
      return bimolecular(index).short_name();
      
    case INNER:
      //
      return inner_barrier(index).short_name();
      
    case OUTER:
      //
      return outer_barrier(index).short_name();

    default:
      //
      std::cerr << funame << "unknown type: " << type << "\n";

      throw Error::Logic();
    }
  }
  
  inline  void assert_spec (spec_t p)
  {
    const char funame [] = "Model::assert_spec: ";
    
    const int& type = p.first;

    const int& index = p.second;
    
    switch(type) {
      //
    case WELL:
      //
      assert_well(index);

      break;
      //
    case BIMOLECULAR:
      //
      assert_bimolecular(index);

      break;
      //
    case INNER:
      //
      assert_inner(index);

      break;
      //
    case OUTER:
      //
      assert_outer(index);

      break;
      //
    default:
      //
      std::cerr << funame << "unknown type: " << type << "\n";

      throw Error::Logic();
    }
  }
  
  inline std::ostream& operator<< (std::ostream& to, spec_t p) { return to << species_name(p); }
  
  inline IO::LogOut&   operator<< (IO::LogOut&   to, spec_t p) { return to << species_name(p); }
  
  // well group state density
  //
  double state_density (const std::set<int>& g, double e);

  // well group weight
  //
  double weight        (const std::set<int>& g, double t);

  // well group ground energy
  //
  double ground        (const std::set<int>& g, int* = 0);

  // barrier thermal energy
  //
  inline double thermal_energy (spec_t b, double temperature)
  {
    const char funame [] = "Model::thermal_energy: ";
    
    switch(b.first) {
      //
    case Model::INNER:
      //
      return Model::inner_barrier(b.second).thermal_energy(temperature);
    
    case Model::OUTER:
      //
      return Model::outer_barrier(b.second).thermal_energy(temperature);
    
    default:
      //
      std::cerr << funame << "wrong dissociation barrier type: " << b.first << "\n";
    
      throw Error::Logic();
    }
  }
  
  // barrier (real) ground energy
  //
  inline double real_ground (spec_t b)
  {
    const char funame [] = "Model::real_ground: ";
    
    switch(b.first) {
      //
    case Model::INNER:
      //
      return Model::inner_barrier(b.second).real_ground();
    
    case Model::OUTER:
      //
      return Model::outer_barrier(b.second).real_ground();
    
    default:
      //
      std::cerr << funame << "wrong dissociation barrier type: " << b.first << "\n";
    
      throw Error::Logic();
    }
  }
  
  //                  energy - barrier map
  //
  //                    energy - barrier
  //                      |        |
  //                      |        |
  typedef std::multimap<double, spec_t > emap_t;

  emap_t  ground_energy_map ();

  emap_t thermal_energy_map (double);

  // reaction definition
  //
  typedef std::set<spec_t> reac_t;

  inline void assert_reac(const reac_t& r)
  {
    if(r.size() != 2) {
      //
      std::cerr << "Model::assert: wrong number of reactants: " << r.size() << "\n";

      throw Error::Logic();
    }
  }
    
  inline std::ostream& operator<< (std::ostream& to, const reac_t& r)
  {
    assert_reac(r);

    to << *r.begin() << " <---> " << *r.rbegin();

    return to;
  }

  inline std::ostream& operator<< (std::ostream& to, const std::list<reac_t>& reactions)
  {
    for(std::list<reac_t>::const_iterator r = reactions.begin(); r != reactions.end(); ++r)
      //
      to << IO::log_offset << std::setw(name_size_max()) << *r << "\n";

    return to;
  }
  
  inline IO::LogOut& operator<< (IO::LogOut& to, const reac_t&            r) { (std::ostream&)to << r; return to; }

  inline IO::LogOut& operator<< (IO::LogOut& to, const std::list<reac_t>& r) { (std::ostream&)to << r; return to; }

  //                 barrier - reaction map
  //
  // rate-limiting barrier ------ reactions list
  //                 |                 |
  //                 |                 |
  typedef std::map<spec_t, std::list<reac_t> > rmap_t;

  rmap_t  barrier_reaction_map (const emap_t& emap);

  inline rmap_t  barrier_reaction_map (double temperature) { return barrier_reaction_map(thermal_energy_map(temperature)); }

  inline rmap_t  barrier_reaction_map ()                   { return barrier_reaction_map(ground_energy_map()); }

  // chemical graph
  //
  struct ChemGraph {
    //
    std::set<int> well_set;

    std::set<int> inner_set;

    std::set<int> outer_set;

    std::list<ChemGraph> factorize () const;

    int split (double, std::list<ChemGraph>* = 0) const;

    ChemGraph operator+ (const ChemGraph&) const;

    ChemGraph subgraph (const std::set<int>&) const;

    bool does_include (const std::set<int>&) const;

    std::set<int> connected_to (int);
    
    void assert () const;

    void assert (spec_t s) const;

    bool is_connected () const { assert(); if(factorize().size() == 1) return true; return false; }
    
    ChemGraph& operator+= (const ChemGraph&);

    emap_t  ground_energy_map () const;
    
    emap_t thermal_energy_map (double t) const;

    ChemGraph () {}

    explicit ChemGraph (double);

    ChemGraph (double, double);

    explicit ChemGraph (const std::set<int>&);
  };

  std::ostream& operator<< (std::ostream&, const ChemGraph&);

  inline IO::LogOut& operator<< (IO::LogOut& to, const ChemGraph& g) { (std::ostream&)to << g; return to; }

  // kinetic landscape
  //
  typedef std::map<spec_t, std::list<ChemGraph> > landscape_t;

  landscape_t        kinetic_landscape (const emap_t& emap);
  
  inline landscape_t kinetic_landscape (double temperature) { return kinetic_landscape(thermal_energy_map(temperature)); }
  
  inline landscape_t kinetic_landscape ()                   { return kinetic_landscape(ground_energy_map()); }
  
  //                            bound groups of wells
  //
  //                        wells group --- dissociation barrier
  //                              |              |
  //                              |              |
  typedef std::list<std::pair<std::set<int>, spec_t > > bound_t;

  bound_t        bound_groups (const emap_t&);

  inline bound_t bound_groups (double temperature) { return bound_groups(thermal_energy_map(temperature)); }

  inline bound_t bound_groups ()                   { return bound_groups(ground_energy_map()); }

  std::vector<double> dissociation_energy_map (double temperature);
  
  double  maximum_barrier_height ();

  // energy shift
  //
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

  // well group output
  //
  std::ostream& operator<< (std::ostream&, const std::set<int>&);

  inline IO::LogOut& operator<< (IO::LogOut& to, const std::set<int>& g) { (std::ostream&)to << g; return to; }

  // well partition output
  //
  std::ostream& operator<< (std::ostream&, const std::vector<std::set<int> >&);
  
  inline IO::LogOut& operator<< (IO::LogOut& to, const std::vector<std::set<int> >& p) { (std::ostream&)to << p; return to; }

}// namespace Model

#endif
