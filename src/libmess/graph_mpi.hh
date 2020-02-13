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

#ifndef GRAPH_MPI_HH
#define GRAPH_MPI_HH

#include "graph_common.hh"
#include "array.hh"

namespace Graph {

  /******************************************************************************************
   ********************** PARTITION FUNCTION GRAPH PERTURBATION THEORY **********************
   ******************************************************************************************/
  
  class Expansion {
    //
    typedef std::map<FreqGraph, int> _mg_t;

    // base class for mpi communications
    //
    class _gbase_t {
      //
      virtual int  _pack   (Array<char>&)  const = 0;

      virtual void _unpack (const Array<char>&)  = 0;

    public:
      //
      virtual ~_gbase_t () {}

      void send (int node, int tag) const;

      void recv (int node, int tag);
    };
    
    class _gmap_t : public std::map<int, double>, public _gbase_t {
      //
      int  _pack   (Array<char>&)  const;

      void _unpack (const Array<char>&);
      
    public:
      //
      _gmap_t () {}

      _gmap_t (const std::map<int, double>& m) : std::map<int, double>(m) {}

      _gmap_t (int from, int tag) { recv(from, tag); }
    };
    
    void   _global_work                (double temperature)                                                   const;
    void _centroid_work                (double temperature)                                                   const;
    void _centroid_work_with_constrain (double temperature, const std::map<std::multiset<int>, double>& mmat) const;
    
    void _work   (int mode, double temperature, const std::map<std::multiset<int>, double>& mmat) const;
    void _master (int mode, double temperature, const std::map<std::multiset<int>, double>& mmat) const;

    // potential expansion
    //
    potex_t  _potex;

    // reduced frequencies
    //
    std::vector<double> _red_freq;

    // reduced frequency index to normal mode index map
    //
    std::vector<std::set<int> > _red_freq_map;

    // normal mode index to reduced frequency index map
    //
    std::vector<int>  _red_freq_index;
  
    void           _set_frequencies (std::vector<double> freq);

    std::set<int>  _low_freq_set    (double temperature, std::vector<double>& tanh_factor) const;

  public:
    //
    enum {MASTER, GLOBAL, CENTROID};

    // tags
    //
    enum {
      END_TAG, 
      WORK_TAG,
      STAT_TAG
    };
    
    void init (const std::vector<double>& freq, const potex_t& potex);

    Expansion () {}

    Expansion (const std::vector<double>& freq, const potex_t& potex) { init(freq, potex); }

    void correction (int                                         mode = GLOBAL, 
		     double                                      temperature = -1., 
		     const std::map<std::multiset<int>, double>& mmat = std::map<std::multiset<int>, double>()
		     ) const;
  
    // logarithmic tolerance for frequencies reduction
    //
    static double freq_tol;

    // low frequency threshold
    //
    static double low_freq_thresh;
  };
}

#endif  
