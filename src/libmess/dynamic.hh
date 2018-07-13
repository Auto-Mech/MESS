

#ifndef DYNAMIC_HH
#define DYNAMIC_HH

#include "structure.hh"
#include "array.hh"
#include "logical.hh"
#include "lapack.hh"

#include <iostream>
#include <string>

namespace Dynamic {

  /*********************************************************************
   *                    Dynamical variable derivatives                 *
   ********************************************************************/

  void set_dvd (const D3::Vector* torque, const double* dv, double* dvd);


  /*********************************************************************
   *                 Cartesian data for Coordinates class              *
   ********************************************************************/

  class CartData 
  {
    int _frag;
    std::vector<D3::Vector> _rel_pos; // relative positions of the fragment atoms in lab frame
    D3::Matrix* _mfo; // rotational matrix for a nonlinear fragment
    double* _ang;     // angular vector for a linear fragment    
    void _delete ();
    void _copy (const CartData&);

  public:
    CartData () : _mfo(0), _ang(0), _frag(-1) {}
    void init (int frag);
    CartData (int frag) : _mfo(0), _ang(0) { init(frag); }

    ~CartData () { _delete(); }
    CartData (const CartData& c) { _copy(c); }
    CartData& operator= (const CartData& c) { _delete(); _copy(c); return *this; }

    void mf2lf (const double* mf, double* lf) const throw(Error::General);
    void lf2mf (const double* lf, double* mf) const throw(Error::General);

    const std::vector<D3::Vector>& rel_pos () const { return _rel_pos; }

    void update_mfo(const double* ang) throw(Error::General);
    void update_rel()  throw(Error::General);

    // inertia moments matrix element
    double imm (int i, int j) const;

  };

  inline void CartData::_delete ()
  {
    if(_mfo) 
      delete _mfo; 
    if(_ang)
      delete[] _ang;
  }

  /*********************************************************************
   *                       Orientational coordinates                   *
   *********************************************************************/

  class Coordinates : private Array<double>
  {
    mutable CartData _fragment [2];
    double* _orb_pos;
    double* _ang_pos [2];

    double _length [2];

    void _init (); // set pointers

    void _normalize (int frag); // normalize fragment angular vector

    enum {MFO = 0, REL = 2, SIZE = 4};
    mutable bool _update [SIZE];
    
    void _init_update ()   const;

    void _update_mfo (int frag) const;
    void _update_rel (int frag) const;

  public:
    Coordinates ();
    Coordinates (const double* dv);
    Coordinates (const Coordinates& c);
    Coordinates& operator= (const Coordinates& c);

    static int size () { return Structure::pos_size(); }

    void get (const double*); // copy from
    void put (double*) const; // copy to
    double length (int frag) const { return _length[frag]; } // the length of original angular vector

    // cm-to-cm vector, 1->2
    double*       orb_pos  ()            { return _orb_pos; }
    const double* orb_pos  ()      const { return _orb_pos; }

    double&       orb_pos  (int i)       { return _orb_pos[i]; }
    double        orb_pos  (int i) const { return _orb_pos[i]; }

    double interfragment_distance () const { return vlength(_orb_pos, 3); }

    // orientational variables
    void    write_ang_pos (int frag, const double* pos) throw(Error::General);

    const double* ang_pos (int frag)              const throw(Error::General);
    double*       ang_pos (int frag)                    throw(Error::General);

    double        ang_pos (int frag, int i)       const throw(Error::General);
    double&       ang_pos (int frag, int i)             throw(Error::General);

    const std::vector<D3::Vector>& rel_pos (int frag)    const;
    void mf2lf (int frag, const double* mf, double* lf)  const;
    void lf2mf (int frag, const double* lf, double* mf)  const;
    Lapack::SymmetricMatrix imm () const;

    static double min_atom_dist;
    bool are_atoms_close () const;

    void print_geom (std::ostream&, const std::string&) const;

    void                        dipole (int frag,                  double* d) const throw(Error::General);
    void     quadrupole_vector_product (int frag, const double* v, double* q) const throw(Error::General);
    void polarizability_vector_product (int frag, const double* v, double* p) const throw(Error::General);

  };

  inline Coordinates::Coordinates (const double* dv) : Array<double>(size())
  {
    // set pointers and initialize cartesian data
    _init();

    // copy data and normalize angular vectors
    get(dv);

    // initialize update info
    _init_update();
  }

  inline const double* Coordinates::ang_pos  (int frag) const throw(Error::General)
  {
    const char funame [] = "Dynamic::Coordinates::ang_pos: ";
    
#ifdef DEBUG

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_pos[frag]; 
  }

  inline double* Coordinates::ang_pos  (int frag) throw(Error::General)
  {
    const char funame [] = "Dynamic::Coordinates::ang_pos: ";
    
#ifdef DEBUG

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_pos[frag]; 
  }

  inline double Coordinates::ang_pos  (int frag, int i) const throw(Error::General)
  {
    const char funame [] = "Dynamic::Coordinates::ang_pos: ";
    
#ifdef DEBUG

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_pos[frag][i]; 
  }

  inline double& Coordinates::ang_pos  (int frag, int i) throw(Error::General)
  {
    const char funame [] = "Dynamic::Coordinates::ang_pos: ";
    
#ifdef DEBUG

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_pos[frag][i]; 
  }

  inline void Coordinates::mf2lf(int frag, const double* mf, double* lf) const
  {
    if(_update[MFO + frag]) _update_mfo(frag);

    _fragment[frag].mf2lf(mf, lf); 
  }

  inline void Coordinates::lf2mf(int frag, const double* lf, double* mf) const
  {
    if(_update[MFO + frag]) _update_mfo(frag);

    _fragment[frag].lf2mf(lf, mf);
  }

  inline void Coordinates::_update_mfo (int frag) const 
  {
    _fragment[frag].update_mfo(ang_pos(frag));

    _update[MFO + frag] = false;
  }

  inline const std::vector<D3::Vector>& Coordinates::rel_pos (int frag) const
  {
    if(_update[REL + frag]) _update_rel(frag);

    return _fragment[frag].rel_pos();
  }

  inline void Coordinates::_update_rel (int frag) const
  {
    if(_update[MFO + frag]) _update_mfo(frag);
  
    _fragment[frag].update_rel();

    _update[REL + frag] = false;
  }

  /*************************************************************************
   ************************ Generalized Momenta ****************************
   *************************************************************************/

  class Momenta : private Array<double>
  {
    double* _orb_vel;
    double* _ang_vel [2];

    void _init(); // set pointers

  public:
    static int size () { return Structure::vel_size(); }

    Momenta () : Array<double>(size()) { _init(); }
    Momenta (const double* dv) : Array<double>(size()) { Array<double>::operator=(dv);  _init(); }
    Momenta (const Momenta& m) : Array<double>((const Array<double>&)m) { _init(); }
    Momenta& operator= (const Momenta& m) { Array<double>::operator=((const Array<double>&)m); return *this; }

    void get (const double* dv) { Array<double>::operator=(dv); }
    void put (double*) const;

    Momenta& operator*= (double d) { Array<double>::operator*=(d); return *this; }
    Momenta& operator/= (double d) { Array<double>::operator/=(d); return *this; }

    // rotational frequencies
    double*       ang_vel (int frag)            throw(Error::General);
    const double* ang_vel (int frag)      const throw(Error::General);
    double        ang_vel (int frag, int) const throw(Error::General);
    double&       ang_vel (int frag, int)       throw(Error::General);

    // cm-to-cm velocity
    double*       orb_vel ()            { return _orb_vel; }
    const double* orb_vel ()      const { return _orb_vel; }
    double        orb_vel (int i) const { return _orb_vel[i]; }
    double&       orb_vel (int i)       { return _orb_vel[i]; }

    double& operator[] (int i)       { return Array<double>::operator[](i); }
    double  operator[] (int i) const { return Array<double>::operator[](i); }

    double*       begin ()       { return Array<double>::begin(); }
    const double* begin () const { return Array<double>::begin(); }
    double*         end ()       { return Array<double>::end(); }
    const double*   end () const { return Array<double>::end(); }

    double orbital_kinetic_energy        () const;
    double total_kinetic_energy          () const;
    double fragment_rotational_energy (int) const;
  };

  inline const double* Momenta::ang_vel (int frag) const throw(Error::General)
  { 

#ifdef DEBUG

    const char funame [] = "Dynamic::Momenta::ang_vel: ";

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_vel[frag];
  }

  inline double* Momenta::ang_vel (int frag) throw(Error::General)
  { 
#ifdef DEBUG
    const char funame [] = "Dynamic::Momenta::ang_vel (int): ";

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }
#endif

    return _ang_vel[frag];
  }

  inline double Momenta::ang_vel (int frag, int i) const throw(Error::General)
  { 

#ifdef DEBUG

    const char funame [] = "Dynamic::Momenta::ang_vel: ";

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_vel[frag][i];
  }

  inline double& Momenta::ang_vel (int frag, int i)  throw(Error::General)
  { 

#ifdef DEBUG

    const char funame [] = "Dynamic::Momenta::ang_vel: ";

    if(Structure::fragment(frag).type() == Molecule::MONOATOMIC) {
      std::cerr << funame << "should not be called for monoatomic molecule\n";
      throw Error::Logic();
    }

#endif

    return _ang_vel[frag][i];
  }

  /*************************************************************************
   *************************** Dynamic Variables ***************************
   *************************************************************************/

  class Vars : public Coordinates, public Momenta
  {
  public:
    static int size () { return Structure::dv_size(); }

    Vars () {}
    Vars (const Coordinates& dc) : Coordinates(dc) {}
    Vars (const double* dv): Coordinates(dv), Momenta(dv + Coordinates::size()) {}

    void get (const double* dv) 
    { Coordinates::get(dv); Momenta::get(dv + Coordinates::size()); }

    void put (double* dv) const
    { Coordinates::put(dv); Momenta::put(dv + Coordinates::size()); }

    void    total_angular_momentum (double*)      const;
    void  orbital_angular_momentum (double*)      const;
    void fragment_angular_momentum (int, double*) const;

    D3::Vector    total_angular_momentum  ()    const;
    D3::Vector  orbital_angular_momentum  ()    const;
    D3::Vector fragment_angular_momentum  (int) const;

    double radial_kinetic_energy         ()    const;
    double total_angular_momentum_length ()    const;
    double angular_momentum_k_projection (int) const;
    double angular_momentum_m_projection (int) const throw(Error::General);
  };

  inline D3::Vector Vars::total_angular_momentum () const 
  { 
    D3::Vector res; 
    total_angular_momentum(res); 
    return res;
  }

  inline D3::Vector Vars::orbital_angular_momentum () const 
  { 
    D3::Vector res; 
    orbital_angular_momentum(res); 
    return res;
  }

  inline D3::Vector Vars::fragment_angular_momentum (int frag) const 
  { 
    D3::Vector res; 
    fragment_angular_momentum(frag, res); 
    return res;
  }

  /*************************************************************************
   ******************************** CLASSIFIER *****************************
   *************************************************************************/

  // abstract class to classify the configuration
  class Classifier {
  public:
    virtual int  classify (const Coordinates&)  const =0;
    virtual ~Classifier () {}
  };

  /*************************************************************************
   ******************************* CONDITION *******************************
   *************************************************************************/

  // abstract logical condition
  class Condition
  {
  public:
    virtual bool test (const Coordinates&) const =0;
    virtual void print (std::ostream&, const std::string&) const {}
    virtual ~Condition () {}
  };

  typedef ConstSharedPointer<Condition> CCP;

  // configuration region which is not sampled 
  extern CCP exclude_region;

  //distance condition
  class DistanceCondition : public Condition
  {
    double _dist;
  public:
    DistanceCondition (double d) : _dist(d) {}
    bool test (const Coordinates&) const;
    ~DistanceCondition () {}
  };

  inline bool DistanceCondition::test (const Coordinates& dc) const 
  {
    if(::vlength(dc.orb_pos(), 3) > _dist)
      return true;
    return false;
  }

  // compound conditions
  // binary condition
  class BinaryCondition : public Condition
  {
    Logical::BinExpr::Oper _op;
    CCP _c1;
    CCP _c2;

    BinaryCondition(CCP c1, CCP c2, Logical::BinExpr::Oper o) : _op(o), _c1(c1) , _c2(c2) {}

  public:
    bool test (const Coordinates&) const;

    ~BinaryCondition () {}

    friend CCP operator& (CCP c1, CCP c2); 
    friend CCP operator| (CCP c1, CCP c2); 
  };

  inline CCP operator& (CCP c1, CCP c2) { return CCP(new BinaryCondition(c1, c2, Logical::BinExpr::AND)); }
  inline CCP operator| (CCP c1, CCP c2) { return CCP(new BinaryCondition(c1, c2, Logical::BinExpr::OR)); }

  inline CCP& operator&= (CCP& c1, CCP c2) { c1 = c1 & c2; return c1; }
  inline CCP& operator|= (CCP& c1, CCP c2) { c1 = c1 | c2; return c1; }


  inline bool BinaryCondition::test (const Coordinates& dc) const
  {
    const char funame [] = "Dynamic::BinaryCondition::test: ";

    switch(_op) {
    case Logical::BinExpr::AND: 
      return  _c1->test(dc) && _c2->test(dc);
    case Logical::BinExpr::OR: 
      return  _c1->test(dc) || _c2->test(dc);
    default:
      std::cerr << funame << "unknown operation\n";
      throw Error::Logic();
    }
  }

  // negation condition
  class UnaryCondition : public Condition
  {
    CCP _c;

    UnaryCondition(CCP c) : _c(c) {}

  public:
    bool test (const Coordinates& dc) const { return !_c->test(dc); }

    ~UnaryCondition () {}

    friend CCP negate (CCP c);
  };

  inline CCP negate (CCP c) { return CCP(new UnaryCondition(c)); }

  /*************************************************************************
   ********************* CONFIGURATIONAL OBSERVABLE ************************
   *************************************************************************/


  class ObservBase {
  public:
    virtual double evaluate (const Coordinates&) const =0;
  };

  // atom-atom distance, atom-atom-atom angle, or atom-atom-atom-atom dihedral
  class ZmatObserv : public ObservBase {
    std::vector<int> _atom;

    static D3::Vector _vector(int, int, const Coordinates&);

  public:
    explicit ZmatObserv (const std::vector<int>&);
    double evaluate (const Coordinates&) const;
  };

  class ObservCondition : public Condition {
    ConstSharedPointer<ObservBase> _observ;
    std::pair<double, double>      _range;
    int                            _mode;

  public:
    enum {LOWER_BOUND, 
	  UPPER_BOUND, 
	  RANGE_BOUND
    };

    ObservCondition (ConstSharedPointer<ObservBase> p, double limit, int m);

    ObservCondition (ConstSharedPointer<ObservBase> p, double lower, double upper);

    bool test (const Coordinates&) const;
  };


  /*************************************************************************
   ********************************** WATCH ********************************
   *************************************************************************/

  // watch if some condition has changed along the trajectory
  class Watch
  {
    CCP _cond;
    bool _start;
    bool _changed;

  public:
    Watch (CCP c, const Coordinates& dv) : _cond(c), _start(c->test(dv)), _changed(false) {}
    Watch (CCP c, bool s) : _cond(c), _start(s), _changed(false) {}

    void test (const Coordinates&);
    bool is_changed () const { return _changed; }
  };

  inline void Watch::test (const Coordinates& c) 
  { 
    if(_changed) 
      return;
    if(_cond->test(c) != _start)
      _changed = true;
  }

  /*************************************************************************
   ********************************** ACTION *******************************
   *************************************************************************/

  class Action
  {
  public:
    virtual void execute (const Vars&)  =0;
  };

}// Dynamic

#endif
