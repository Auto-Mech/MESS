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

#ifndef GRAPH_COMMON_HH
#define GRAPH_COMMON_HH

#include<vector>
#include<set>
#include<map>
#include<iostream>

namespace Graph {

  // potential expansion type
  //
  typedef std::map<std::multiset<int>, double> potex_t;
  
  // maximal number of bonds in the graph
  //
  extern int  bond_max;

  // maximal potential expansion order
  //
  extern int potex_max;

  class  GenGraph;

  typedef std::vector<GenGraph>::const_iterator const_iterator;

  const_iterator begin ();

  int size ();

  void init ();

  bool isinit ();

  /*******************************************************************************
   ********************** GENERIC PERTURBATION THEORY GRAPH **********************
   *******************************************************************************/
  
  class GenGraph: public std::multiset<std::multiset<int> > {
    //
    void _check () const;

  public:
    std::vector<std::multiset<int> > vertex_bond_map () const; // vertex to bond map
    //
    bool                             is_connected    () const; // is graph connected
    //
    int                              vertex_size     () const; // number of vertices
    //
    int                              bond_size       () const; // number of bonds between different vertices
    //
    int                              loop_size       () const; // number of bonds between same vertices
    //
    int                              vertex_symmetry () const; // vertex permution symmetry factor
    //
    int                              bond_symmetry   () const; // different points of the vertex symmetry (wick theorem)
    //
    int                              symmetry_factor () const { return vertex_symmetry() * bond_symmetry(); }
    //
    std::vector<int>                 bond_order      () const; // sequential bond orders (because bonds between the same vertices go sequentially)

    void  print (std::ostream&) const; // graph output
  };

  inline std::ostream& operator<< (std::ostream& to, const GenGraph& g) { g.print(to); return to; }

  /*********************************************************************
   ********************** FREQUENCY ADAPTED GRAPH **********************
   *********************************************************************/
  
  class FreqGraph: public std::map<std::set<int>, std::multiset<int> > {
    //
    void _check_integrity (int = 0) const;
    //
    void _check_order     ()        const;

    typedef std::map<FreqGraph, int> _mg_t;

    // cyclic structure
    //
    class _Ring : private std::vector<int> {
      //
      void _assert(int) const;
      
    public:
      //
      explicit _Ring (const std::vector<int>&);

      int size () const { return std::vector<int>::size(); }
      
      int find (const std::set<int>&) const;

      std::set<int> edge (int) const;
    };

    class _NoRing {};
    
    _Ring _min_ring (const std::set<int>&) const;

    double _four_term (int                        mi,
		       double                     temperature,
		       const std::vector<double>& rfreq = std::vector<double>(),
		       const std::vector<double>& ifreq = std::vector<double>()) const;

  public:
    //
    int vertex_size () const;

    int   bond_size () const;

    bool is_connected (int, int) const;

    void print (std::ostream&) const;

    // factorize graph into the set of connected graphs
    //
    _mg_t factorize   () const;
    
    // permutationally equivalent graphs
    //
    std::set<FreqGraph> perm_pool   (int* =0, int =0) const;

    // reduce graph with strongly coupled vertices
    //
    FreqGraph reduce (const std::vector<double>& freq,
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

    // fourier sum cutoff
    //
    static int    four_cut;

    static double four_par;

    // graph reduction threshold
    //
    static double red_thresh;
  };

  inline std::ostream& operator<< (std::ostream& to, const FreqGraph& g) { g.print(to); return to; }

  void read_potex (const std::vector<double>& freq, std::istream& from, std::map<std::multiset<int>, double>& potex);
}

#endif  
