

#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include "symmetry.hh"
#include "array.hh"
#include "shared.hh"
#include "lapack.hh"
#include "monom.hh"

#include <map>
#include <vector>
#include <set>

/*************************************************************************************************
 ***** CONFIGURATIONAL SPACE FOR INTERPOLATION OF THE POTENTIAL BETWEEN TWO RIGID FRAGMENTS ******
 *************************************************************************************************/

namespace Configuration {

/*******************************************************************************************
 *************************** CONFIGURATIONAL SPACE VECTOR LAYOUT ***************************
 *******************************************************************************************/

  // first index describes vector subspace dimension, 
  // other indices - diminsions of spheres of arbitrary dimensions
 
  class Layout {
    std::vector<int> _layout;
    std::vector<int> _shift;
    
    int   _linear_size;
    int  _simplex_size;
    int _fragment_type;

    void _isinit () const throw(Error::General);

  public:
    enum {ATOMIC, LINEAR, NONLINEAR};

    void set (const std::vector<int>&) throw(Error::General);

    Layout () {}
    Layout (const std::vector<int>& l) throw(Error::General) { set(l); }

    int simplex_size ()      const { _isinit(); return _simplex_size;   }    
    int  linear_size ()      const { _isinit(); return  _linear_size;   }
    int         size ()      const { _isinit(); return  _layout.size(); }
    int         size (int i) const { _isinit(); return  _layout[i];     }
    int        shift (int i) const { _isinit(); return  _shift[i];      }

    int orientation   () const;
    int radius_vector () const;

    bool operator== (const Layout& l) const { _isinit(); return _layout == l._layout; }
    bool operator!= (const Layout& l) const { _isinit(); return _layout != l._layout; }

    int  invariant_size () const { _isinit(); return size(0) == 1 ? 1 : 0; }
    int second_fragment () const { _isinit(); return _fragment_type; }
  };

  inline void Layout::_isinit () const throw(Error::General)
  {
    const char funame [] = "Configuration::Layout::_isinit: ";

    if(_layout.size())
      return;

    std::cerr << funame << "not initialized\n";
    throw Error::Init();
  }

  inline int Layout::radius_vector () const
  {
    _isinit();

    if(second_fragment() == ATOMIC)
      return *_shift.rbegin();
    else
      return *(_shift.rbegin() + 1);
  }

  inline int Layout::orientation () const
  {
    const char funame [] = "Configuration::Layout::orientation: ";

    _isinit();

    if(second_fragment() != ATOMIC)
      return *_shift.rbegin();

    std::cerr << funame << "wrong layout\n";
    throw Error::Logic();
  }

/*******************************************************************************************
 ******************************** CONFIGURATIONAL STATE VECTOR *****************************
 *******************************************************************************************/

  class State : public Array<double> {
  public:
    static Layout layout;

    State () : Array<double>(layout.linear_size()) {}

    double*       radius_vector ()       { return *this + layout.radius_vector(); }
    const double* radius_vector () const { return *this + layout.radius_vector(); }

    double*         orientation ()       { return *this + layout.orientation(); }
    const double*   orientation () const { return *this + layout.orientation(); }

    template <typename V> 
    State& operator= (const V& v) { Array<double>::operator=(v); return *this; }
  };

  inline std::ostream& operator << (std::ostream& to, const State& s) 
  {
    for(const double* it = s.begin(); it != s.end(); ++it)
      to << *it << " ";
    return to << "\n";
  }

  inline std::istream& operator >> (std::istream& from, State& s) 
  {
    for(double* it = s.begin(); it != s.end(); ++it)
      from >> *it;
    return from;
  }

/*******************************************************************************************
 ************ ABSTRACT CLASS FOR SYMMETRY OPERATIONS IN THE CONFIGURATIONAL SPACE **********
 *******************************************************************************************/

  class GroupBase : public Symmetry::GroupBase {
    
  public:
    GroupBase () {}
    GroupBase (const std::vector<Permutation>& p) : Symmetry::GroupBase(p) {}

    virtual void apply (int, const State&, State&) const = 0;
  };

/*******************************************************************************************
 **************** SYMMETRY GROUP FOR NONLINEAR + LINEAR OR ATOMIC FRAGMENTS ****************
 *******************************************************************************************/

  class SpaceGroup : public GroupBase { 
    bool _symmetric;
    std::vector<Symmetry::SpaceElement> _base;

  public:
    SpaceGroup (const Symmetry::SpaceGroup&, int =0) throw(Error::General);

    void apply (int, const State&, State&) const throw(Error::General);
  };

/*******************************************************************************************
 ****************** SYMMETRY GROUP FOR NONLINEAR + NONLINEAR FRAGMENTS *********************
 *******************************************************************************************/

  // for identical fragments symmetry operation assumed to have a form (q, p) * I, 
  // where q is the spatial symmetry operation on the first fragment (reference),
  // p is the spatial symmetry operation on the second (moving) fragment, and
  // I is the operation of exchange of idential fragments.
  // (q, p): (r, g) -> (q*r*q-1, q*g/p), where r is cm-to-cm radius-vector and 
  // g is the second fragment orientation vector (quaternion).
  // I: (r, g) -> (- 1/g * r * g, 1/g).

  class DoubleSpaceGroup : public GroupBase {

    // two fragments are identical
    bool _identical;

    // symmetry elements of individual fragments
    std::vector<std::vector<Symmetry::SpaceElement> > _symmetry_element;

    // map of the symmetry group index to the indices of the fragment symmetry groups, 
    // (optional) fragment exchange index, and quaternion sign change index
    std::vector<std::vector<int> >  _index_multi_map;

    enum {
      SIGN      = 2, // quaternion sign change map index
      EXCHANGE  = 3  // fragment exchange map index
   };

  public:
    DoubleSpaceGroup (const std::vector<Symmetry::SpaceGroup>&, int =0) throw(Error::General);

    void apply (int, const State&, State&) const throw(Error::General);

    Lapack::Matrix qmatrix (int) const throw(Error::General);
    Lapack::Matrix rmatrix (int) const throw(Error::General);
  };

/*******************************************************************************************
 ****************** PARTITIONING CONFIGURATIONAL SPACE INTO SIMPLEXES **********************
 *******************************************************************************************/

  class Simplexation {

    static double _det_tol;        // determinant tolerance
    static double _dis_tol;        // distance tolerance
    static int    _miss_count_max; // maximal number of misses at random initialization

    ConstSharedPointer<GroupBase> _symm_group;

    // data types
    enum {
      NOT_INIT,
      AB_INITIO,
      INTERPOLATE,
      NO_DATA
    };
    
    struct _Attribute {
      int attribute;
      _Attribute () : attribute(NOT_INIT) {}
    };

    struct _Double : public _Attribute {
      double value;
    };

    struct _Simplex {
      std::set<int> vertex;
      State         center;
      double        radius;
      _Double         data;
    };

    std::vector<State>       _vertex;
    std::vector<std::vector<int> >    _vertex_orbit;     // symmetry related vertices
    std::vector<std::pair<int, int> > _vertex_orbit_map; // vertex-to-vertex-orbit map

    std::vector<_Simplex>          _simplex;
    std::vector<std::vector<int> > _simplex_orbit; // symmetry related simplices

    std::map<std::set<int>, double> _vertex_vertex_distance;

    std::map<std::set<int>, std::set<int> > _facet_simplex_map;
    
    
    // ...

    int _iterate_simplex_center (const std::set<int>&, State&, double&) const;
    int    _find_simplex_center (const std::set<int>&, double, State&, double&, double&) const;

  public:
    explicit Simplexation (ConstSharedPointer<GroupBase> sg) : _symm_group(sg) {}

    void random_init (double, double = -1., double = -1.) throw(Error::General);

    // ...
  };

}

#endif
