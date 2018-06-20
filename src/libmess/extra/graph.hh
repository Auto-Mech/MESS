/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef GRAPH_HH
#define GRAPH_HH

#include<vector>
#include<set>
#include<map>
#include<cstdarg>
//#include<multiset>

#include "slatec.hh"
#include "shared.hh"
#include "lapack.hh"
#include "multindex.hh"
#include "atom.hh"
#include "math.hh"
#include "io.hh"

  /******************************************************************************************
   ********************** PARTITION FUNCTION GRAPH PERTURBATION THEORY **********************
   ******************************************************************************************/
  
class GraphExpansion {

  // graph connecting different vertices
  //
  class _graph_t: public std::multiset<std::multiset<int> > {
    //
    void _check () const;

  public:
    std::vector<std::multiset<int> > vertex_bond_map () const; // vertex to bond map
    bool                             is_connected    () const; // is graph connected
    int                              vertex_size     () const; // number of vertices
    int                              bond_size       () const; // number of bonds between different vertices
    int                              loop_size       () const; // number of bonds between same vertices
    int                              vertex_symmetry () const; // vertex permution symmetry factor
    int                              bond_symmetry   () const; // different points of the vertex symmetry (wick theorem)
    int                              symmetry_factor () const { return vertex_symmetry() * bond_symmetry(); }
    std::vector<int>                 bond_order      () const; // sequential bond orders (because bonds between the same vertices go sequentially)

    void  print (std::ostream&) const; // graph output
  };

  friend std::ostream& operator<< (std::ostream& to, const _graph_t& g);

  class _fg_t;

  typedef std::map<_fg_t, int> _mg_t;
    
  // frequency adapted graph
  //
  class _fg_t: public std::map<std::set<int>, std::multiset<int> > {
    //
    void _check_integrity (int = 0) const;
    //
    void _check_order     ()          const;

  public:
    //
    int                  vertex_size () const;

    int                  bond_size   () const;

    // factorize graph into the set of connected graphs
    //
    _mg_t                factorize   () const;
    
    // permutationally equivalent graphs
    //
    std::set<_fg_t>      perm_pool   (int* =0, int =0) const;

    // reduce graph with strongly coupled vertices
    //
    _fg_t  reduce (const std::vector<double>& freq,
		   //
		   double                     temperature,
		   //
		   const std::vector<double>& tanh_factor,
		   //
		   _mg_t&                     zpe_graph
		   //
		   )  const;
      
    // low  temperature / high frequency integral evaluation
    //
    double zpe_factor (const std::vector<double>& freq,
		       //
		       double                     temperature = -1.,
		       //
		       const std::vector<double>& tanh_factor = std::vector<double>()
		       //
		       ) const;
      
    // high temperature / low frequency integral evaluation
    //
    double fourier_sum (const std::vector<double>& freq, double temperature) const;
  };

  // frequency adapted graph converter (to save memory)
  //
  class _Convert {
    //
    int _vertex_size_max;
    
    int _freq_size;

  public:
    //
    typedef char         int_t;
    
    typedef Array<int_t> vec_t;

    static int mem_size (const vec_t& v) { return v.size() * sizeof(int_t); }
      
    _Convert () : _vertex_size_max(0), _freq_size(0) {}
      
    void init (int vertex_size_max, int freq_size);
      
    vec_t operator() (const _fg_t&) const;
  };
      
  _Convert _convert;
    
  // database format
  //
  struct _gmap_t : public std::map<_Convert::vec_t, double> {
    //
    long mem_size() const;
  };

  std::set<_graph_t>                     _graph_data;
  //
  std::vector<_graph_t>                  _sorted_graph;

  typedef std::map<std::multiset<int>, double> _potex_t;
  //
  _potex_t _potex;

  // reduced (and degenerate) frequencies
  //
  std::vector<double>                 _red_freq;

  // reduced frequency index to normal mode index map
  //
  std::vector<std::set<int> >     _red_freq_map;

  // normal mode index to reduced frequency index map
  //
  std::vector<int>              _red_freq_index;
  
  static std::set<_graph_t> _raw_graph_generator (std::vector<int>, int = -1);

  void           _set_frequencies (std::vector<double> freq);
  
  std::set<int>  _low_freq_set    (double temperature, std::vector<double>& tanh_factor) const;

public:
  //
  class DeepTunnel : public Exception::Base {
    
  public:
    //
    template <typename T>
    DeepTunnel& operator<< (const T& t) { (Exception::Base&)(*this) << t; return *this; } 
  };

  enum { KEEP_PERM = 1};

  void init (const std::vector<double>& freq, const _potex_t& potex);
    
  GraphExpansion () {}
    
  GraphExpansion (const std::vector<double>& freq, const _potex_t& potex) { init(freq, potex); }

  std::map<int, double>           correction (double temperature = -1.) const;
  //
  std::map<int, double>  centroid_correction (double temperature = -1.) const;
  //
  std::map<int, double>  centroid_correction (const std::map<std::multiset<int>, double>&, double temperature = -1.)  const;
  
  int size () const { return _graph_data.size(); }
  
  static int           mod_flag; // different modification flags
  //
  static int           bond_max; // maximum number of bonds in the graph
  //
  static double        freq_tol; // logarithmic tolerance for frequencies reduction
  //
  static double        four_cut; // fourier sum graph evaluation cutoff
  //
  static double      red_thresh; // low temperature / high frequency graph reduction threshold
  //
  static double low_freq_thresh; // low frequency threshold

  static void read_potex (const std::vector<double>& freq, std::istream& from, _potex_t& potex);
};

inline std::ostream& operator<< (std::ostream& to, const GraphExpansion::_graph_t& g) { g.print(to); return to; }

#endif  
