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

#ifndef GRAPH_OMP_HH
#define GRAPH_OMP_HH

#include "graph_common.hh"
#include "array.hh"

  /******************************************************************************************
   ********************** PARTITION FUNCTION GRAPH PERTURBATION THEORY **********************
   ******************************************************************************************/

namespace Graph { 
 
  class Expansion {

    typedef std::map<FreqGraph, int> _mg_t;
    
    // frequency adapted graph converter (to save memory)
    //
    class _Convert {
      //
      int _vertex_size_max;
    
      int _freq_size;

      std::vector<std::set<int> > _index_map;

    public:
      //
      typedef char         int_t;
    
      typedef Array<int_t> vec_t;

      static int mem_size (const vec_t& v) { return v.size() * sizeof(int_t); }
      
      _Convert () : _vertex_size_max(0), _freq_size(0) {}
      
      void init (int vertex_size_max, int freq_size);
      
      vec_t operator() (const FreqGraph&) const;

      FreqGraph operator() (const vec_t&) const;
    };
      
    _Convert _convert;
    
    // database format
    //
    struct _gmap_t : public std::map<_Convert::vec_t, double> {
      //
      long mem_size() const;
    };

    // potential expansion
    //
    potex_t _potex;

    // reduced frequencies
    //
    std::vector<double> _red_freq;

    // reduced frequency index to normal mode index map
    //
    std::vector<std::set<int> > _red_freq_map;

    // normal mode index to reduced frequency index map
    //
    std::vector<int> _red_freq_index;
  
    void _set_frequencies (std::vector<double> freq);
  
    std::set<int> _low_freq_set (double temperature, std::vector<double>& tanh_factor) const;

  public:
    //
    enum { KEEP_PERM = 1};

    void init (const std::vector<double>& freq, const potex_t& potex);
    
    Expansion () {}
    
    Expansion (const std::vector<double>& freq, const potex_t& potex) { init(freq, potex); }

    std::map<int, double>           correction (double temperature = -1.) const;
    //
    std::map<int, double>  centroid_correction (double temperature = -1.) const;
    //
    std::map<int, double>  centroid_correction (const std::map<std::multiset<int>, double>&, double temperature = -1.)  const;
  
    // different modification flags
    //
    static int mod_flag;

    // logarithmic tolerance for frequencies reduction
    //
    static double freq_tol;
  
    // low frequency threshold
    //
    static double low_freq_thresh;
  };
}

#endif  
