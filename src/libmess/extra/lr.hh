

#ifndef LR_HH
#define LR_HH

#include "potential.hh"
#include "math.hh"

#include<iostream>

namespace LongRange {

  void init (std::istream&) throw(Error::General);
  void set_angular_grid (int g);

  enum sym_t { 
    SPHERICAL,     // sperically symmetric molecule: atom or spherical top
    LINEAR,        // linear molecule
    SYMMETRIC,     // symmetric top
    GENERIC        // asymmetric top
  };

  int   symmetric_fragment_size () throw(Error::General);
  sym_t symmetry_type   (int frag) throw(Error::General);  // fragment symmetry
  int   orientational_dimension () throw(Error::General);
  int   eff_dof                 () throw(Error::General); // effective number of degrees of freedom

  extern Potential::Wrap pot;

  // minimum energy search functions
  double minimum_energy (double dist);
  double gradient_search (Dynamic::Coordinates& dc, double ener_var);
  extern "C" void pos2drv (const double& t, const double* pos, double* drv, void*, void* par);
  

  // abstract class to represent the number of states orientational density
  class StatesNumberDensity {
    double _ener; // energy

  protected:

    double _density (const Dynamic::Coordinates& dc, int dim, double* mep =0) const;

  public:
    double energy () const { return _ener; }
    void set_ener (double e) { _ener = e; }

    // number of states (orientational integral)
    double integral (double dist, double* mep =0) const throw(Error::General);

    StatesNumberDensity (double e) : _ener(e) {}

    virtual double operator () (const Dynamic::Coordinates&, double* mep =0) const = 0;
    virtual double norm_factor () const = 0;
    virtual ~StatesNumberDensity () {}
  };

  
  
  // Number of states orientational density wrapper 
  class StatesNumberDensityWrap
  {
    ConstSharedPointer<StatesNumberDensity> _den;
	
  public:

    StatesNumberDensityWrap () {}
    StatesNumberDensityWrap (const StatesNumberDensity* p) : _den(p) {}
    StatesNumberDensityWrap (ConstSharedPointer<StatesNumberDensity> p) : _den(p) {}

    operator const StatesNumberDensity* () const { return _den; } 
    void   isinit     () const throw(Error::General);
    double operator() (const Dynamic::Coordinates&, double* mep =0) const throw(Error::General);
    double norm_factor () const throw(Error::General);
    ~StatesNumberDensityWrap () {}
  };

  inline void StatesNumberDensityWrap::isinit () const throw(Error::General)
  {
    const char funame [] = "LongRange::StatesNumberDensityWrap::isinit: ";

    if(!_den) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }
  }

  inline double StatesNumberDensityWrap::operator() (const Dynamic::Coordinates& dc, double* mep) const throw(Error::General)
  {
    isinit();
    return (*_den)(dc, mep);
  }

  inline double StatesNumberDensityWrap::norm_factor() const throw(Error::General)
  {
    isinit();
    return _den->norm_factor();
  }



  // Microcanonical number of states orientational density, spherically symmetric fragment degrees of freedom excluded
  class EDensity : public StatesNumberDensity  {

    static int    _dim;
    static double _nfac;

  public:
    static void init () throw(Error::General);

    EDensity (double e) : StatesNumberDensity(e) {}
    double operator() (const Dynamic::Coordinates& dc, double* mep =0) const;
    double norm_factor () const { return _nfac; }
    ~EDensity () {}
    
  };

  inline double EDensity::operator() (const Dynamic::Coordinates& dc, double* mep) const { return _density(dc, _dim, mep); }



  // EJ-resolved number of states orientational density, spherically symmetric fragment degrees of freedom excluded
  class JDensity : public StatesNumberDensity  {

    static int    _dim;
    static double _nfac;

  public:
    static void init () throw(Error::General);

    JDensity (double e) : StatesNumberDensity(e) {}
    double operator() (const Dynamic::Coordinates& dc, double* mep =0) const;
    double norm_factor () const { return _nfac; }
    ~JDensity () {}
    
  };

  inline double JDensity::operator() (const Dynamic::Coordinates& dc, double* mep) const { return _density(dc, _dim, mep); }



  // number of states orientational density with angular momentum projections on symmetry axes, if any, conserved
  // and spherical top degrees of freedom axcluded
  class MDensity : public StatesNumberDensity {

    static int    _dim;
    static double _nfac;

  public:
    static void init () throw(Error::General);

    MDensity (double e) : StatesNumberDensity(e) {}
    double operator() (const Dynamic::Coordinates& dc, double* mep =0) const;
    double norm_factor () const { return _nfac; }

    ~MDensity () {}

  };

  inline double MDensity::operator() (const Dynamic::Coordinates& dc, double* mep) const { return _density(dc, _dim, mep); }


  // number of states orientational density with angular momentum projections on symmetry axes, 
  // if any, conserved, and spherical top degrees of freedom axcluded
  // and angular momentum projection on the interfragment axis conserved
  // orbital (j-shift) and symmetry axes rotational energies (m-shift) being excluded from total energy 
  class KDensity : public StatesNumberDensity {

    static int    _dim;
    static double _nfac;

    double _m_proj [2];
    double _k_proj;
    
  public:
    static void init () throw(Error::General);
 
    void set_k_proj (double k) { _k_proj = k; }
    void set_m_proj (int i, double m);

    KDensity (double e, const std::vector<double>& m, double k) throw(Error::General);
    double operator() (const Dynamic::Coordinates&, double* mep =0) const;
    double norm_factor () const { return _nfac; }

    ~KDensity () {}
  };



  /*******************************************************************************************
   *************************************** Number of States Classes **************************
   *******************************************************************************************/

  // centrifugal energy shift
  double j_shift (double amom, double dist) throw(Error::General);

 // symmetric molecule inner rotation energy shift
  double m_shift (const std::vector<double>& m_proj) throw(Error::General);

  // centrifugal barrier
  double centrifugal_barrier (double amom);

  // base abstract class for the number of states
  class StatesNumber {
  public:
    virtual double operator() () = 0; // actual number of states calculation
    virtual double distance     () const =0;
    virtual void set_dist (double)       =0;
  };



  // abstract class representing E-resolved number of states
  class ENumber : public StatesNumber {
    double _ener;

  public:
    double ener () const { return _ener; }
    virtual void set_ener (double e) { _ener = e; }


    ENumber (double e =0.);
  };


  inline ENumber::ENumber (double e) : _ener(e) {}




  // abstract class representing J-resolved number of states
  class JNumber : public StatesNumber {
    double _ener;
    double _amom;

  public:
    double ener () const { return _ener; }
    double amom () const { return _amom; }

    virtual void set_ener (double e) { _ener = e; }
    virtual void set_amom (double j) { _amom = j; }

    JNumber (double e =0., double j =0.);
  };

  inline JNumber::JNumber (double e, double j) : _ener(e), _amom(j) {}



  // abstract class representing M-resolved number of states
  class MNumber : public StatesNumber {
    double _ener;
    double _amom;
    std::vector<double> _m_proj;

  protected:
    double _m_shift () const { return m_shift(_m_proj); }

  public:
    double ener () const { return _ener; }
    double amom () const { return _amom; }
    double m_proj (int i) const { return _m_proj[i]; }

    virtual void set_ener (double e) { _ener = e; }
    virtual void set_amom (double j) { _amom = j; }
    virtual void set_m_proj (int i, double m) { _m_proj[i] = m; }

    // minimal energy at m_proj=0 to set upper limits of m-integral
    virtual double ener_min () const = 0;

    MNumber (double e =0., double j =0., const std::vector<double>& m =
	     std::vector<double>(symmetric_fragment_size(), 0.)) throw(Error::General);
  };

  inline MNumber::MNumber (double e, double j, const std::vector<double>& m) throw(Error::General)
    : _ener(e), _amom(j), _m_proj(m) 
  {
    static const char funame [] = "LongRange::MNumber::MNumber: ";
    if(_m_proj.size() != symmetric_fragment_size()) {
      std::cerr << funame << "wrong number of symmetric fragments\n";
      throw Error::Init();
    }
  }



  // abstract class representing K-resolved number of states
  class KNumber : public StatesNumber {
    double _ener;
    double _amom;
    std::vector<double> _m_proj;
    double _k_proj;

  protected:
    double _m_shift () const { return m_shift(_m_proj); }

  public:
    double ener () const { return _ener; }
    double amom () const { return _amom; }
    double m_proj (int i) const { return _m_proj[i]; }
    double k_proj () const   { return _k_proj; }

    virtual void set_ener (double e) { _ener = e; }
    virtual void set_amom (double j) { _amom = j; }
    virtual void set_m_proj (int i, double m) { _m_proj[i] = m; }
    virtual void set_k_proj (double k) { _k_proj = k; }

    // minimal energy at k,m=0 to set limits for m-integral
    virtual double ener_min () const =0;

    KNumber (double e =0., double j =0., const std::vector<double>& m =
	     std::vector<double>(symmetric_fragment_size(), 0.), double k =0.) throw(Error::General);
  };

  inline KNumber::KNumber (double e, double j, const std::vector<double>& m, double k) throw(Error::General)
    : _ener(e), _amom(j), _m_proj(m), _k_proj(k) 
  {
    static const char funame [] = "LongRange::MNumber::MNumber: ";
    if(_m_proj.size() != symmetric_fragment_size()) {
      std::cerr << funame << "wrong number of symmetric fragments\n";
      throw Error::Init();
    }
  }



  /********************* N-INTEGRALS ****************************/


  // J-integral of J-resolved number of states
  class JIntegral : public ENumber {
    JNumber* _num;
    double _integral () const; // integral of J-resolved number of states over J

  public:
    static double step;

    mutable double arg_max;  // maximum value of J 
    mutable double arg_var;  // J variance 

    double distance ()         const { return _num->distance();  }
    void   set_dist (double d)       {        _num->set_dist(d); }

    JIntegral (JNumber* n, double e =0.);
    void set_ener (double e) { ENumber::set_ener(e); _num->set_ener(e); }

    double operator() () { return _integral(); }
  };

  inline JIntegral::JIntegral (JNumber* n, double e) : ENumber(e), _num(n), arg_max(-1.), arg_var(0.) { _num->set_ener(e); }



  // M-integral of M-resolved number of states
  class MIntegral : public JNumber {
    MNumber* _num;
    double _integral (int =0) const; // recursive integration over m projections (if any)

  public:
    static double step[2];      // integration step

    mutable double arg_max[2];  // maximum value of M
    mutable double arg_var[2];  // M variance 

    double distance ()         const { return _num->distance();  }
    void   set_dist (double d)       {        _num->set_dist(d); }

    MIntegral(MNumber* n, double e =0., double j=0.);
    
    void set_ener (double e) { JNumber::set_ener(e); _num->set_ener(e); }
    void set_amom (double j) { JNumber::set_amom(j); _num->set_amom(j); }
    
    double operator() () { return _integral(); } // integral of M-resolved number of states over M
  };

  inline MIntegral::MIntegral(MNumber* n, double e, double j) 
    : JNumber(e, j), _num(n)
  { _num->set_ener(e); _num->set_amom(j); arg_max[0] = arg_max[1] = -1.; arg_var[0] = arg_var[1] = 0.; }
    


  // K-integral of K-resolved number of states
  class KIntegral : public MNumber {
    KNumber* _num;
    double _integral () const; // integral of K-resolved number of states over K 
  public:
    static double step;

    mutable double arg_max;  // maximum value of K 
    mutable double arg_var;  // K variance 

    double distance ()         const { return _num->distance();  }
    void   set_dist (double d)       {        _num->set_dist(d); }

    KIntegral (KNumber* n, double e =0., double j =0., const std::vector<double>& m = 
	       std::vector<double>(symmetric_fragment_size(), 0.)) throw(Error::General);
    
    void set_ener          (double e) { MNumber::set_ener(e);      _num->set_ener(e); }
    void set_amom          (double j) { MNumber::set_amom(j);      _num->set_amom(j); }
    void set_m_proj (int i, double m) { MNumber::set_m_proj(i, m); _num->set_m_proj(i, m); }

    double ener_min () const { return _num->ener_min(); }

    double operator() () { return _integral(); }
  };

  inline KIntegral::KIntegral (KNumber* n, double e, double j, const std::vector<double>& m) 
    throw(Error::General) : MNumber(e, j, m), _num(n), arg_max(-1.), arg_var(0.) 
  {
    _num->set_ener(e); 
    _num->set_amom(j); 
    for(int i = 0; i < m.size(); ++i) 
      _num->set_m_proj(i, m[i]);
  }



  /********************* DENSITY INTEGRALS AT A FIXED DISTANCE ****************************/


  // E-resolved number of states orientational density integral
  class EDensityIntegral : public ENumber, private EDensity {
    double _dist;

  public:
    double ener     () const { return ENumber::ener(); }
    double distance () const { return _dist; }
    double integral () const { return StatesNumberDensity::integral(_dist); }
 
    EDensityIntegral(double d, double e =0.);

    void set_ener (double e) { ENumber::set_ener(e); StatesNumberDensity::set_ener(e); }
    virtual void set_dist (double d) { _dist = d; }

  };

  inline EDensityIntegral::EDensityIntegral(double d, double e) : ENumber(e), EDensity(e), _dist(d) {}


  
  // E-resolved number of states at a fixed distance (orientational density integral)
  class ENumberFix : public EDensityIntegral {
    double _vmin;

  public:
 
    ENumberFix(double d, double e =0.);

    void set_dist (double d) { EDensityIntegral::set_dist(d); _vmin = minimum_energy(d); }

    double operator() ();
  };

  inline ENumberFix::ENumberFix (double d, double e)
    : EDensityIntegral(d, e), _vmin(minimum_energy(d)) {}

  inline double ENumberFix::operator() ()
  { 
    if(ener() <= _vmin)
      return -1.; 
    return integral(); 
  }




  // J-resolved number of states orientational density integral
  class JDensityIntegral : public JNumber, private JDensity {
    double _dist;
    double _shift () const { return j_shift(amom(), _dist); }

  public:
    double shifted_energy () const { return StatesNumberDensity::energy(); }
    double ener           () const { return JNumber::ener(); }
    double distance       () const { return _dist; }
    double integral       () const { return  StatesNumberDensity::integral(_dist); }

    JDensityIntegral (double d, double e =0., double j =0.); 

    void         set_ener (double e) { JNumber::set_ener(e); 
      StatesNumberDensity::set_ener(ener() - _shift()); }
    virtual void set_amom (double j) { JNumber::set_amom(j); 
      StatesNumberDensity::set_ener(ener() - _shift()); }
    virtual void set_dist (double d) { _dist = d;
      StatesNumberDensity::set_ener(ener() - _shift()); } 

  };


  inline JDensityIntegral::JDensityIntegral (double d, double e, double j) 
    : JNumber(e, j), JDensity(e - j_shift(j, d)), _dist(d) {}



  // J-resolved number of states at a fixed distance (orientational density integral)
  class JNumberFix : public JDensityIntegral {
    double _vmin;

  public:

    JNumberFix(double d, double e =0., double j =0.); 

    void set_dist (double d) { JDensityIntegral::set_dist(d); _vmin = minimum_energy(d); } 

    double operator() ();
  };

  inline JNumberFix::JNumberFix (double d, double e, double j)
    : JDensityIntegral(d, e, j), _vmin(minimum_energy(d)) {}

  inline double JNumberFix::operator() ()
  { 
    if(shifted_energy() <= _vmin)
      return -1.; 
    return integral(); 
  }




  // M-resolved number of states orientational density integral 
  class MDensityIntegral : public MNumber, private MDensity {
    double _dist;
    double _shift () const { return _m_shift() + j_shift(amom(), _dist); }

  public:
    double shifted_energy () const { return StatesNumberDensity::energy(); }
    double ener           () const { return MNumber::ener(); }
    double distance       () const { return _dist; }
    double integral       () const { return  StatesNumberDensity::integral(_dist); }
 
    MDensityIntegral (double d, double e =0., double j =0., const std::vector<double>& m =
		      std::vector<double>(symmetric_fragment_size(), 0.)) throw(Error::General) ; 

    void         set_ener (double e) { MNumber::set_ener(e); StatesNumberDensity::set_ener(ener() - _shift()); }
    virtual void set_amom (double j) { MNumber::set_amom(j); StatesNumberDensity::set_ener(ener() - _shift()); }
    virtual void set_dist (double d) { _dist = d;            StatesNumberDensity::set_ener(ener() - _shift()); }

    void set_m_proj (int i, double m) { MNumber::set_m_proj(i, m); StatesNumberDensity::set_ener(ener() - _shift()); }
  };

  inline MDensityIntegral::MDensityIntegral (double d, double e, double j, const std::vector<double>& m) 
    throw(Error::General) : MNumber(e, j, m), MDensity(e - m_shift(m) - j_shift(j, d)), _dist(d) {}



  // M-resolved number of states at a fixed distance (orientation integral)
  class MNumberFix : public MDensityIntegral {
    double _vmin;

  public:
    double ener_min () const { return _vmin + j_shift(amom(), distance()); }
  
    MNumberFix(double d, double e =0., double j =0., const std::vector<double>& m =
	       std::vector<double>(symmetric_fragment_size(), 0.)) throw(Error::General); 

    void set_dist (double d) { MDensityIntegral::set_dist(d); _vmin = minimum_energy(d); }

    double operator() ();
  };

  inline MNumberFix::MNumberFix (double d, double e, double j, const std::vector<double>& m)
    throw(Error::General) : MDensityIntegral(d, e, j, m), _vmin(minimum_energy(d)) {}

  inline double MNumberFix::operator() () 
  { 
    if(shifted_energy() <= _vmin)
      return -1.;
    return integral(); 
  }



  // K-resolved number of states orientational density integral
  class KDensityIntegral : public KNumber, private KDensity {
    double _dist;
    double _shift () const { return _m_shift() + j_shift(amom(), _dist); }

  public:
    double shifted_energy () const { return StatesNumberDensity::energy(); }
    double ener           () const { return KNumber::ener(); }
    double distance       () const { return _dist; }
    double integral       () const { return  StatesNumberDensity::integral(_dist); }


    KDensityIntegral (double d, double e =0., double j =0., const std::vector<double>& m =
		      std::vector<double>(symmetric_fragment_size(), 0.), double k =0.) throw(Error::General);


    void         set_ener (double e) { KNumber::set_ener(e); StatesNumberDensity::set_ener(ener() - _shift()); }
    virtual void set_amom (double j) { KNumber::set_amom(j); StatesNumberDensity::set_ener(ener() - _shift()); }
    virtual void set_dist (double d) { _dist = d;            StatesNumberDensity::set_ener(ener() - _shift()); }

    void set_m_proj (int i, double m) { KNumber::set_m_proj(i, m); KDensity::set_m_proj(i, m); 
      StatesNumberDensity::set_ener(ener() - _shift()); }
    void set_k_proj        (double k) { KNumber::set_k_proj(k);    KDensity::set_k_proj(k); }
    
  };

  
  inline KDensityIntegral::KDensityIntegral (double d, double e, double j, const std::vector<double>& m, double k)
    throw(Error::General) : KNumber(e, j, m, k), KDensity(e - m_shift(m) - j_shift(j, d), m, k), _dist(d) {}



  // K-resolved number of states at a fixed distance (orientational density integral)
  class KNumberFix : public KDensityIntegral {
    double _vmin; // minimum potential energy

  public:
    double ener_min () const { return _vmin + j_shift(amom(), distance()); }

    KNumberFix (double d, double e =0., double j =0., const std::vector<double>& m = 
	       std::vector<double>(symmetric_fragment_size(), 0.), double k =0.) throw(Error::General);

    void set_dist (double d) { KDensityIntegral::set_dist(d); _vmin = minimum_energy(d); } 

    double operator() ();
  };

  inline KNumberFix::KNumberFix (double d, double e, double j, const std::vector<double>& m, double k) throw(Error::General)
    : KDensityIntegral(d, e, j, m, k), _vmin(minimum_energy(d)) {}


  inline double KNumberFix::operator() ()
  { 
    if(shifted_energy() <= _vmin)
      return -1.; 
    return integral(); 
  }



  /********************* DENSITY INTEGRALS VARIATIONALLY OPTIMIZED ****************************/



  // abstract class to optimize j-resolved number of states
  class Opt : private Math::MinimumSearch {
    int _opt_step_count; // number of times _fun has been called during the search (for debugging purposes)
    double _fun (double) throw (Math::MinimumSearch::Exception);

  public:

    static double xtol;
    static double ytol;
    static double step;

    virtual void   set_dist (double) =0;
    virtual double distance () const =0;
    virtual double integral () const =0;

    double opt_integral () throw(Error::General);

    Opt () : Math::MinimumSearch(xtol, ytol) {}
  };

  class Bar {
  protected:
    double _ebar; // centrifugal barrier energy
    //double _dbar; // centrifugal barrier distance

    void _set_bar (double j) { _ebar = centrifugal_barrier(j); } 
   
    Bar (double j) { _set_bar(j); }
  };

  // K-resolved number of states, optimized over distance (orientational density integral)
  class KNumberOpt : public KDensityIntegral, private Opt, private Bar {
  public:
    KNumberOpt (double d, double e =0., double j =0., const std::vector<double>& m = 
	       std::vector<double>(symmetric_fragment_size(), 0.), double k =0.) throw(Error::General);

    double distance   () const   { return KDensityIntegral::distance(); }
    double integral   () const   { return KDensityIntegral::integral(); }
    double ener_min   () const   { return _ebar; }
    
    void   set_amom   (double j) { KDensityIntegral::set_amom(j); _set_bar(j); }
    void   set_dist   (double d) { KDensityIntegral::set_dist(d); }

    double operator() ();
  }; 
    
  inline double KNumberOpt::operator() () { 
    if(ener() <= _ebar + _m_shift())
      return -1.;
    
    return opt_integral();
  }

  inline KNumberOpt::KNumberOpt (double d, double e, double j, const std::vector<double>& m, double k)
    throw(Error::General) : KDensityIntegral(d, e, j, m, k), Bar(j) {}



  // M-resolved number of states, optimized over distance (orientational density integral)
  class MNumberOpt : public MDensityIntegral, private Opt, private Bar {
  public:
    MNumberOpt (double d, double e =0., double j =0., const std::vector<double>& m = 
		std::vector<double>(symmetric_fragment_size(), 0.)) throw(Error::General); 

    double distance   () const   { return MDensityIntegral::distance(); }
    double integral   () const   { return MDensityIntegral::integral(); }
    double ener_min   () const   { return _ebar; }

    void   set_amom   (double j) { MDensityIntegral::set_amom(j); _set_bar(j); }
    void   set_dist   (double d) { MDensityIntegral::set_dist(d); }

    double operator() ();
  }; 
    
  inline double MNumberOpt::operator() () { 
    if(ener() <= _ebar)
      return -1.;
    
    return opt_integral();
  }

  inline MNumberOpt::MNumberOpt (double d, double e, double j, const std::vector<double>& m)
    throw(Error::General) : MDensityIntegral(d, e, j, m), Bar(j) {}



  // J-resolved number of states, optimized over distance (orientational density integral)
  class JNumberOpt : public JDensityIntegral, private Opt, private Bar {
  public:
    JNumberOpt (double d, double e =0., double j =0.);

    double distance     () const { return JDensityIntegral::distance(); }
    double integral     () const { return JDensityIntegral::integral(); }
    
    void   set_amom   (double j) { JDensityIntegral::set_amom(j); _set_bar(j); }
    void   set_dist   (double d) { JDensityIntegral::set_dist(d); }

    double operator() ();
  }; 
    
  inline double JNumberOpt::operator() () { 
    if(ener() <= _ebar)
      return -1.;
    
    return opt_integral();
  }

  inline JNumberOpt::JNumberOpt (double d, double e, double j)
    : JDensityIntegral(d, e, j), Bar(j) {}

  // E-resolved number of states, optimized over distance (orientational density integral)
  class ENumberOpt : public EDensityIntegral, private Opt {
  public:
    ENumberOpt (double d, double e =0.);

    double distance     () const { return EDensityIntegral::distance(); }
    double integral     () const { return EDensityIntegral::integral(); }
    
    void   set_dist   (double d) { EDensityIntegral::set_dist(d); }

    double operator() ();
  }; 
    
  inline double ENumberOpt::operator() () { 
    if(ener() <= 0.)
      return -1.;
    
    return opt_integral();
  }

  inline ENumberOpt::ENumberOpt (double d, double e)
    : EDensityIntegral(d, e) {}

}

#endif
