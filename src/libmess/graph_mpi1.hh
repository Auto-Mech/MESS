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

#ifndef GRAPH_MPI_1_HH 
#define GRAPH_MPI_1_HH

#include<vector>
#include<set>
#include<map>
#include<list>
#include<cstdarg>

#include<mpi.h>

#include "slatec.hh"
#include "shared.hh"
#include "lapack.hh"
#include "multindex.hh"
#include "atom.hh"
#include "math.hh"
#include "io.hh"
#include "array.hh"

/******************************************************************************************
 ********************** PARTITION FUNCTION GRAPH PERTURBATION THEORY **********************
 ******************************************************************************************/
  
class GraphExpansion {
  //
public:
  static int WORK_NODE;

  // node types
  enum { 
    MASTER,    // driver
    INT_SERV,  // integral database server 
    ZPE_SERV,  // zpe database server
    SUM_SERV  // sum database server
  };

  enum {GLOBAL, CENTROID};

private:
  // standard graph connecting different vertices
  class _graph_t: public std::multiset<std::multiset<int> > {
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
  class _fg_t: public std::map<std::set<int>, std::multiset<int> >
  {
    void _check_integrity (int = -1) const;
    void _check_order     ()        const;

  public:
    //_fg_t () {}
    //_fg_t (const std::map<std::set<int>, std::multiset<int> >& g) : std::map<std::set<int>, std::multiset<int> >(g) {}
      
    int              vertex_size () const;
    int              bond_size   () const;
    _mg_t            factorize   () const; // factorize graph into connected graphs in the standard form
    std::set<_fg_t>  perm_pool   (int* =0, int =0) const; // permutationally equivalent diagrams

    void print (std::ostream&) const;

    // reduce graph with strongly coupled vertices
    _fg_t  reduce (const std::vector<double>& freq,
		   double                     temperature,
		   const std::vector<double>& tanh_factor,
		   _mg_t&                     zpe_graph,
		   const std::set<int>&       low_freq = std::set<int>()
		   )  const;
      
    // low  temperature integral evaluation
    double zpe_factor (const std::vector<double>& freq,
		       double                     temperature = -1.,
		       const std::vector<double>& tanh_factor = std::vector<double>(),
		       const std::set<int>&       low_freq    = std::set<int>()
		       ) const;
      
    // high temperature integral evaluation
    double fourier_sum (const std::vector<double>& freq,
			double                     temperature,
			const std::set<int>&       low_freq    = std::set<int>()
			) const;
  };

  friend std::ostream& operator<< (std::ostream& to, const _fg_t& g);

  // frequency adapted graph converter (to save memory)
  class _Convert {
  private:
    int _vertex_size_max;
    int _freq_size;

    std::vector<std::set<int> > _index_map;

  public:
    typedef short              int_t;
    typedef Array<int_t>       vec_t;
      
    static int mem_size (const vec_t& v) { return v.size() * sizeof(int_t); }
      
    void init (int v, int f);
      
    _Convert () : _vertex_size_max(0), _freq_size(0) {}

    _Convert (int v, int f) { init(v, f); }
      
    vec_t operator() (const _fg_t&) const;
    _fg_t operator() (const vec_t&) const;
  };
      
  _Convert _convert;
    
  bool _is_work_node (int node) const {  return node >= WORK_NODE; }
  bool _is_driver    (int node) const {  return node > SUM_SERV && node < WORK_NODE; }
  bool _is_server    (int node) const {  return node && node <= SUM_SERV; }
  bool _is_master    (int node) const {  return !node; }
   
  // tags
  enum {
    END_TAG,        // end calcuation
    STD_TAG,        // standard graph request / result
    INT_TAG,        // int request / result
    SUM_TAG,        // sum result
    SUM_REQUEST,    // sum request to the sum server
    ZPE_TAG,        // zpe result
    ZPE_REQUEST,    // zpe request to the zpe server
    GINDEX_TAG,     // graph index request / result      
    GINDEX_REQUEST, // graph index request / result      
    NODE_TAG,       // working node request / result
    NODE_RELEASE,   // working node release request
    STAT_TAG        // statistical information
  };

  void   _global_driver                (double temperature)                                                   const;
  void _centroid_driver                (double temperature)                                                   const;
  void _centroid_driver_with_constrain (double temperature, const std::map<std::multiset<int>, double>& mmat) const;
    
  void _driver      (double temperature, int mode = GLOBAL, 
		     const std::map<std::multiset<int>, double>& mmat = 
		     std::map<std::multiset<int>, double>()) const; // driver

  void _master      (double temperature, int mode = GLOBAL) const; // master process: working nodes dispatcher and graph index distributor
  void _int_server  (double temperature)                    const; // integral    database server
  void _zpe_server  (double temperature)                    const; // zpe         database server
  void _sum_server  (double temperature)                    const; // fourier sum database server
  void _node_work   (double temperature)                    const; // node work

  struct _red_t {
    _mg_t           zpe_group;
    _fg_t           red_graph;
    std::set<_fg_t> raw_group;
    double          value;

    _red_t () : value(1.) {}
  };
    
  // base class for mpi communications
  class _gbase_t {
  protected:
    enum { 
      IFREE = 0,      // buffer is ifree 
      ISEND = 1,      // buffer is used by isend non-blocking communication
      IRECV = 2,      // buffer is used by irecv non-blocking communication
      BUFF_SIZE = 200 // buffer size in bytes
    };

    static _Convert _convert;
      
    mutable Array<char> _buff;

    mutable MPI::Request _request;

    mutable int _state; // sending, receiving, or free
      
    virtual int  _pack   (Array<char>&)  const = 0;
 
    virtual void _unpack (const Array<char>&)  = 0;

  public:
    _gbase_t () : _state(IFREE) {}

    void isend (int node, int tag, int flag = 0) const;
    void irecv (int node, int tag) const;
      
    void send (int node, int tag) const;
    void recv (int node, int tag);

    bool Test ();
    bool Test (MPI::Status&);

    void Wait ();
    void Wait (MPI::Status&);

    virtual ~_gbase_t ();
  };
    
  // frequency adapted graph itself
  class _gself_t : public _fg_t, public _gbase_t
  {
    int  _pack   (Array<char>&) const;
    void _unpack (const Array<char>&);
      
  public:
    _gself_t () {}
    _gself_t (const _fg_t& fg) : _fg_t(fg) {}
    _gself_t (int from, int tag) { recv(from, tag); }
  };

  // frequency adapted graph with the value (for mpi communication)
  class _gin_t : public std::pair<_fg_t, double>, public _gbase_t
  {
    int  _pack   (Array<char>&)  const;
    void _unpack (const Array<char>&);
      
  public:
    _gin_t () {}
    _gin_t (const _fg_t g, double d = 0.);
    _gin_t (int from, int tag) { recv(from, tag); }
  };

  // frequency adapted graph pair (for mpi communication)
  class _gpair_t : public std::pair<_fg_t, _fg_t>, public _gbase_t
  {
    int  _pack   (Array<char>&)  const;
    void _unpack (const Array<char>&);
      
  public:
    _gpair_t () {}
    _gpair_t (const _fg_t& rg, const _fg_t& sg) : std::pair<_fg_t, _fg_t>(rg, sg) {}
    _gpair_t (int from, int tag) { recv(from, tag); }
  };

  class _gval_t : public std::pair<int, double>, public _gbase_t {
    int  _pack   (Array<char>&)  const;
    void _unpack (const Array<char>&);
      
  public:
    _gval_t () {}
    _gval_t (int i, double v) : std::pair<int, double>(i, v) {}
    _gval_t (int from, int tag) { recv(from, tag); }
  };
    
  class _gmap_t : public std::pair<int, std::map<int, double> >, public _gbase_t {
    int  _pack   (Array<char>&)  const;
    void _unpack (const Array<char>&);
      
  public:
    _gmap_t () {}
    _gmap_t (int i, const std::map<int, double>& v) : std::pair<int, std::map<int, double> >(i, v) {}
    _gmap_t (int from, int tag) { recv(from, tag); }
  };
    
  class _gnode_t : public _gbase_t {
    int _value;

    int  _pack   (Array<char>&)  const;
    void _unpack (const Array<char>&);
      
  public:
    _gnode_t (int v) { _value = v; }
    _gnode_t (int from, int tag) { recv(from, tag); }

    operator int () const { return _value; }
  };

  static void _update_symbolic_map (std::map<int, std::map<_fg_t, double> >& sym_map,
				    const std::map<_fg_t, double>& received,
				    std::list<SharedPointer<_gbase_t> >* sending = 0);

  static void _update_symbolic_map (std::map<int, std::map<_mg_t, double> >& sym_map,
				    const std::map<_fg_t, double>& received,
				    std::list<SharedPointer<_gbase_t> >* sending = 0);

  static void _update_symbolic_map (std::map<int, std::map<int, std::map<_mg_t, double> > >& sym_map,
				    const std::map<_fg_t, double>& received,
				    std::list<SharedPointer<_gbase_t> >* sending = 0);

  // database format
  struct _db_t : public std::map<_Convert::vec_t, double> {
    long mem_size() const;
  };

  std::set<_graph_t>                     _graph_data;
  std::vector<_graph_t>                  _sorted_graph;

  typedef std::map<std::multiset<int>, double> _potex_t;
  //
  _potex_t _potex;

  std::vector<double>                 _red_freq; // reduced (and degenerate) frequencies
  std::vector<std::set<int> >     _red_freq_map; // reduced frequency index to normal mode index map
  std::vector<int>              _red_freq_index; // normal mode index to reduced frequency index map
  
  static std::set<_graph_t> _raw_graph_generator (std::vector<int>, int = -1);

  void           _set_frequencies (std::vector<double> freq);
  std::set<int>  _low_freq_set    (double temperature, std::vector<double>& tanh_factor) const;

  //std::vector<double> _adjust_frequencies (double temperature, std::vector<double>& tanh_factor) const;
  
public:
  void init (const std::vector<double>& freq, const PotentialExpansion& potex);

  GraphExpansion () {}
  GraphExpansion (const std::vector<double>& freq, const PotentialExpansion& potex) { init(freq, potex); }

  void correction (int mode = GLOBAL, double temperature = -1., const std::map<std::multiset<int>, double>& =
		   std::map<std::multiset<int>, double>()) const;
  
  int size () const { return _graph_data.size(); }
  
  static int           bond_max; // maximum number of bonds in the graph
  static double        freq_tol; // logarithmic tolerance for frequencies reduction
  static double        four_cut; // fourier sum graph evaluation cutoff
  static double      red_thresh; // low temperature / high frequency graph reduction threshold
  static double low_freq_thresh; // low frequency threshold
};

inline std::ostream& operator<< (std::ostream& to, const GraphExpansion::_graph_t& g) { g.print(to); return to; }
inline std::ostream& operator<< (std::ostream& to, const GraphExpansion::_fg_t&    g) { g.print(to); return to; }

#endif  
