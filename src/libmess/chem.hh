#ifndef CHEM_HH
#define CHEM_HH

#include <set>

#include "permutation.hh"
#include "atom.hh"

// chemical structure
//
namespace Chem {

  /*************************************************************************
   ************************* TOLERANCES AND LIMITS *************************
   *************************************************************************/

  extern double angle_tolerance, distance_tolerance;

  bool are_angles_equal (double, double);

  bool are_distances_equal (double, double);
  
  double max_bond_length(const AtomBase&, const AtomBase&);


  /**********************************************************************************
   *********************************** GRAPH ****************************************
   **********************************************************************************/

  class Graph: private std::set<std::set<int> > {

    int _vsize;// vertex size

    void _assert (int) const;

    void _assert ()    const;

    void _isinit ()    const;
    
    bool _is_ring ()    const; // is simple ring

    std::set<int> _find_neighbor (int)               const; // subset of nearest neigbors

    int _valence (int v) const { return _find_neighbor(v).size(); }  

    Graph _projection (const std::vector<int>&)      const; // graph projection on the subset of vertices

    std::set<Permutation> _default_symmetry_group () const; // brute force symmetry group calculation

  public:

    void init (int, const std::set<std::set<int> >&);

    Graph () : _vsize(0) {}

    explicit Graph (int s, const std::set<std::set<int> >& g) : _vsize(0) { init(s, g); }

    int vertex_size () const { return _vsize; }
    int   edge_size () const { return size(); }

    const std::set<std::set<int> >& base () const { return *this; }

    bool is_connected ()    const; // is graph connected

    bool is_bond (int, int) const;

    int distance (int, int) const;
    
    std::set<Permutation> symmetry_group () const;
  };
  
}

#endif
