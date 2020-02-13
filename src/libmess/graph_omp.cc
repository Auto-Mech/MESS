#include "graph_omp.hh"
#include "key.hh"
#include "io.hh"
#include "permutation.hh"
#include "multindex.hh"
#include "units.hh"

#include <cmath>

// different modification flags
//
int Graph::Expansion::mod_flag = Graph::Expansion::KEEP_PERM;

// reduced frequency logarithmic tolerance
//
double Graph::Expansion::freq_tol = 0.1;

// low frequency threshold
//
double Graph::Expansion::low_freq_thresh = 0.2;

/********************************************************************************************
 ************************ PERTURBATION THEORY GRAPH EXPANSION *******************************
 ********************************************************************************************/

std::map<int, double> Graph::Expansion::correction (double temperature) const
{
  const char funame [] = "Graph::Expansion::correction: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;

  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "temperature(K) = " << temperature / Phys_const::kelv << "\n\n";
  }
  else {
    //
    IO::log << IO::log_offset << "zero-point energy correction calculation:" << "\n\n";
    
    for(int i = 0.; i < _red_freq.size(); ++i)
      //
      if(_red_freq[i] <= 0.) {
	ErrOut err_out;
	err_out << funame << "for zero-point energy calculations all frequencies should be real";
      }
  }

  IO::log << IO::log_offset << std::setw(5) << "#";

  if(temperature > 0.) {
    //
    IO::log << std::setw(15) << "Value";
  }
  else {
    //
    IO::log << std::setw(15) << "Value, 1/cm";
  }

  IO::log << std::setw(5) << "V#"
	  << std::setw(5) << "B#"
	  << std::setw(5) << "L#" 
	  << std::setw(7) << "Time"
	  << "   "        << "Graph"   
	  << std::endl;
  
    
  std::vector<double> tanh_factor;
  //
  std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);
  
  int sum_calc_tot = 0;
  int zpe_calc_tot = 0;
  int int_calc_tot = 0;
    
  long sum_read_tot = 0;
  long zpe_read_tot = 0;
  long int_read_tot = 0;
    
  int sum_miss_tot = 0;
  int zpe_miss_tot = 0;
  int int_miss_tot = 0;
    
  _gmap_t int_data;
  _gmap_t sum_data;
  _gmap_t zpe_data;

  std::map<int, double> corr;

#ifndef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_calc_tot, zpe_calc_tot, int_calc_tot, sum_read_tot, zpe_read_tot, int_read_tot, sum_miss_tot, zpe_miss_tot, int_miss_tot) private(itemp, dtemp, btemp) schedule(dynamic)
#endif
        
  for(int gindex = 0; gindex < Graph::size(); ++gindex) {
    //
    Graph::const_iterator graphit = Graph::begin() + gindex;

    std::time_t  start_time = std::time(0);

    const std::vector<std::multiset<int> > vertex_map = graphit->vertex_bond_map();
    const int vertex_size = vertex_map.size();

    double gvalue = 0.;

    int sum_calc = 0;
    int zpe_calc = 0;
    int int_calc = 0;
    
    long sum_read = 0;
    long zpe_read = 0;
    long int_read = 0;
    
    int sum_miss = 0;
    int zpe_miss = 0;
    int int_miss = 0;
    
    MultiIndexConvert corr_multi_index(graphit->size(), _red_freq_index.size());

#ifdef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_calc, zpe_calc, int_calc, sum_read, zpe_read, int_read, sum_miss, zpe_miss, int_miss, gvalue) private(itemp, dtemp, btemp) schedule(dynamic)
#endif

    for(long corr_lin = 0; corr_lin < corr_multi_index.size(); ++corr_lin) {
      //
      std::vector<int> corrin = corr_multi_index(corr_lin);

      double gfactor = 1.;

      btemp = false;
      //
      for(std::vector<std::multiset<int> >::const_iterator mit = vertex_map.begin(); mit != vertex_map.end(); ++mit) {
	//
	std::multiset<int> potex_sign;
	//
	for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it) {
	  //
	  potex_sign.insert(corrin[*it]);
	}
	
	potex_t::const_iterator pexit = _potex.find(potex_sign);
	
	if(pexit != _potex.end()) {
	  //
	  gfactor *= pexit->second;
	}
	else {
	  //
	  btemp = true;
	  //
	  break;
	}
      }
      
      if(btemp) {
	//
	continue;
      }

      // frequency adapted graph
      //
      FreqGraph mod_graph;
      
      itemp = 0;
      //
      for(GenGraph::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {
	//
	int ci = corrin[itemp]; // correlator (normal mode) index
	
	int fi = _red_freq_index[ci]; // reduced frequency index

	// low frequency correlator is a constant
	//
	if(low_freq.find(fi) != low_freq.end()) {
	  //
	  if(_red_freq[fi] > 0.) {
	    //
	    gfactor *=  temperature / _red_freq[fi] / _red_freq[fi];
	  }
	  else {
	    //
	    gfactor *= -temperature / _red_freq[fi] / _red_freq[fi];
	  }

	  continue;
	}

	std::set<int> bond;
	//
	for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
	  //
	  bond.insert(*it);

	// bond loop
	//
	if(bond.size() == 1) {
	  //
	  if(temperature > 0.) {
	    //
	    gfactor /= 2. * _red_freq[fi] * tanh_factor[fi];
	  }
	  else {
	    //
	    gfactor /= 2. * _red_freq[fi];
	  }
	}
	// add frequency index to the graph
	//
	else {
	  //
	  mod_graph[bond].insert(fi);
	}
      }

      itemp = vertex_size - mod_graph.vertex_size();
      //
      if(itemp && temperature > 0.) {
	//
	gfactor /= std::pow(temperature, (double)itemp);
      }

      // frequency adapted graph avaluation
      //
      if(mod_graph.size()) {
	//
	// zero temperature integral (zpe factor) evaluation
	//
	if(temperature <= 0.) {
	  //
	  _Convert::vec_t mod_graph_conv;

	  if(mod_flag & KEEP_PERM) {
	    //
	    mod_graph_conv = _convert(mod_graph);
	  }
	  else {
	    //
	    mod_graph_conv = _convert(*mod_graph.perm_pool().begin());
	  }

#pragma omp critical(zpe_critical)
	  {
	    itemp = zpe_data.find(mod_graph_conv) != zpe_data.end();
	  }
	
	  // read graph value from the database
	  //
	  if(itemp) {
	    //
	    ++zpe_read;
	  
#pragma omp critical(zpe_critical)
	    {
	      gfactor *= zpe_data[mod_graph_conv];
	    }
	  }
	  // zero temperature integral (zpe factor) calculation
	  //
	  else {
	    //
	    ++zpe_calc;

	    dtemp = mod_graph.zpe_factor(_red_freq);
	    //
	    gfactor *= dtemp;

	    // save zero temperature integral value in the database
	    //
#pragma omp critical(zpe_critical)
	    {
	      if(zpe_data.find(mod_graph_conv) != zpe_data.end()) {
		//
		++zpe_miss;
	      }
	      else if(mod_flag & KEEP_PERM) {
		//
		std::set<FreqGraph> pool = mod_graph.perm_pool();

		for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		  //
		  zpe_data[_convert(*pit)] = dtemp;
	      }
	      else {
		//
		zpe_data[mod_graph_conv] = dtemp;
	      }
	    }
	    //
	    //
	  } // zero temperature integral (zpe factor) calculation
	  //
	  //
	} // zero temperature integral (zpe factor) evaluation
	//
	// positive temperature
	//
	else {
	  //
	  // graph factorization into connected graphs
	  //
	  _mg_t fac_graph = mod_graph.factorize();

	  // factorized graph cycle
	  //
	  for(_mg_t::const_iterator fgit = fac_graph.begin(); fgit != fac_graph.end(); ++fgit) {
	    //
	    // whole integral evaluation
	    //
	    _Convert::vec_t fac_graph_conv;

	    if(mod_flag & KEEP_PERM) {
	      //
	      fac_graph_conv = _convert(fgit->first);
	    }
	    else {
	      //
	      fac_graph_conv = _convert(*fgit->first.perm_pool().begin());
	    }

#pragma omp critical(int_critical)
	    {
	      itemp = int_data.find(fac_graph_conv) != int_data.end();
	    }

	    // read whole integral value from the database
	    //
	    if(itemp) {
	      //
	      ++int_read;

#pragma omp critical(int_critical)
	      {
		dtemp = int_data[fac_graph_conv];
	      }
	      
	      for(int i = 0; i < fgit->second; ++i)
		gfactor *= dtemp;
	    }
	    //
	    // whole integral calculation
	    //
	    else {
	      //
	      ++int_calc;
	      
	      double int_val = 1.;

	      // graph reduction
	      //
	      _mg_t zpe_graph;
	      FreqGraph red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);

	      // zpe graph cycle
	      //
	      for(_mg_t::const_iterator zgit = zpe_graph.begin(); zgit != zpe_graph.end(); ++zgit) {
		//
		// low temperature integral (zpe factor) evaluation
		//
		_Convert::vec_t zpe_graph_conv;
	      
		if(mod_flag & KEEP_PERM) {
		  //
		  zpe_graph_conv = _convert(zgit->first);
		}
		else {
		  //
		  zpe_graph_conv = _convert(*zgit->first.perm_pool().begin());
		}
		  
#pragma omp critical(zpe_critical)
		{
		  itemp = zpe_data.find(zpe_graph_conv) != zpe_data.end();
		}

		// read low temperature integral value from the database
		//
		if(itemp) {
		  ++zpe_read;
		  
#pragma omp critical(zpe_critical)
		  {
		    dtemp = zpe_data[zpe_graph_conv];
		  }
		
		  for(int i = 0; i < zgit->second; ++i)
		    int_val *= dtemp;
		}
		// low temperature integral (zpe factor) calculation
		//
		else {
		  ++zpe_calc;
		  
		  dtemp = zgit->first.zpe_factor(_red_freq, temperature, tanh_factor);

		  for(int i = 0; i < zgit->second; ++i)
		    int_val *= dtemp;

		  // save low temperature integral value in the database
		  //
#pragma omp critical(zpe_critical)
		  {
		    if(zpe_data.find(zpe_graph_conv) != zpe_data.end()) {
		      ++zpe_miss;
		    }
		    else if(mod_flag & KEEP_PERM) {
		      //
		      // permutationally equivalent configurations
		      //
		      std::set<FreqGraph> pool = zgit->first.perm_pool();

		      for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
			zpe_data[_convert(*pit)] = dtemp;
		    }
		    else {
		      zpe_data[zpe_graph_conv] = dtemp;
		    }
		  }
		  //
		  //
		} // low temperature integral (zpe factor) calculation
		//
		//
	      } // zpe graph cycle

	      // reduced graph fourier sum evaluation
	      //
	      if(red_graph.size()) {
		//
		_Convert::vec_t red_graph_conv;

		if(mod_flag & KEEP_PERM) {
		  //
		  red_graph_conv = _convert(red_graph);
		}
		else {
		  //
		  red_graph_conv = _convert(*red_graph.perm_pool().begin());
		}
		
#pragma omp critical(sum_critical)
		{
		  itemp = sum_data.find(red_graph_conv) != sum_data.end();
		}

		// read reduced graph fourier sum value from the database
		//
		if(itemp) {
		  ++sum_read;
		
#pragma omp critical(sum_critical)
		  {
		    int_val *= sum_data[red_graph_conv];
		  }
		}
		// reduced graph fourier sum calculation
		//
		else {
		  //
		  ++sum_calc;
		
		  dtemp = red_graph.fourier_sum(_red_freq, temperature);
		  int_val *= dtemp;

		  // save fourier sum calculation result in the database
		  //
#pragma omp critical(sum_critical)
		  {
		    if(sum_data.find(red_graph_conv) != sum_data.end()) {
		      ++sum_miss;
		    }
		    else if(mod_flag & KEEP_PERM) {
		      std::set<FreqGraph> pool = red_graph.perm_pool();
		  
		      for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
			sum_data[_convert(*pit)] = dtemp;
		    }
		    else {
		      sum_data[red_graph_conv] = dtemp;
		    }
		  }
		  //
		  //
		} // reduced graph fourier sum calculation
		//
		//
	      } // reduced graph fourier sum evaluation
	      //
	      else {
		int_val /= temperature;
	      }

	      for(int i = 0; i < fgit->second; ++i)
		gfactor *= int_val;

	      // save whole integral calculation result in the database
	      //
#pragma omp critical(int_critical)
	      {
		if(int_data.find(fac_graph_conv) != int_data.end())
		  ++int_miss;
		else if(mod_flag & KEEP_PERM) {
		  std::set<FreqGraph> pool = fgit->first.perm_pool();

		  for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		    int_data[_convert(*pit)] = int_val;
		}
		else {
		  int_data[fac_graph_conv] = int_val;
		}
	      }
	      //
	      //
	    } // whole integral calculation
	    //
	    //
	  } // factorized graph cycle
	  //
	  //
	} // positive temperature
	//
	//
      } // frequency adapted graph evaluation

      gvalue += gfactor;
      //
      //
    } // normal mode indices cycle
    
    zpe_calc_tot += zpe_calc;
    int_calc_tot += int_calc;
    sum_calc_tot += sum_calc;
    
    zpe_read_tot += zpe_read;
    int_read_tot += int_read;
    sum_read_tot += sum_read;
    
    zpe_miss_tot += zpe_miss;
    int_miss_tot += int_miss;
    sum_miss_tot += sum_miss;
    
    gvalue /= (double)graphit->symmetry_factor();
   
    // odd number of vertices has minus sign
    //
    if(vertex_size % 2)
      gvalue = -gvalue;
     
    // zero-point energy has an opposite sign
    //
    if(temperature < 0.)
      gvalue = -gvalue;

#ifndef INNER_CYCLE_PARALLEL
#pragma omp critical
#endif
    {
      corr[graphit->size()] += gvalue;

      IO::log << IO::log_offset << std::setw(5) << gindex;

      if(temperature > 0.) {
	//
	IO::log << std::setw(15) << gvalue;
      }
      else {
	//
	IO::log << std::setw(15) << gvalue / Phys_const::incm;
      }

      IO::log << std::setw(5) << graphit->vertex_size()
	      << std::setw(5) << graphit->bond_size()
	      << std::setw(5) << graphit->loop_size()
	      << std::setw(7) << std::time(0) - start_time 
	      << "   "        << *graphit 
	      << std::endl;

    }
    //
    //
  } // graph cycle
  
  IO::log << "\n";

  IO::log << IO::log_offset << "Statistics:\n\n";

  IO::log << IO::log_offset;

  if(zpe_calc_tot)
    //
    IO::log << std::setw(10) << "ZPE Calc";

  if(zpe_miss_tot)
    //
    IO::log << std::setw(10) << "ZPE Miss";

  if(zpe_read_tot)
    //
    IO::log << std::setw(15) << "ZPE Read";

  if(zpe_data.size())
    //
    IO::log << std::setw(15) << "ZPE Data";

  if(temperature > 0.) {
    //
    if(int_calc_tot)
      //
      IO::log << std::setw(10) << "Int Calc";

    if(int_miss_tot)
      //
      IO::log<< std::setw(10) << "Int Miss";

    if(int_read_tot)
      //
      IO::log<< std::setw(15) << "Int Read";

    if(int_data.size())
      //
      IO::log << std::setw(15) << "Int Data";

    if(sum_calc_tot)
      //
      IO::log << std::setw(10) << "Sum Calc";
    
    if(sum_miss_tot)
      //
      IO::log << std::setw(10) << "Sum Miss";

    if(sum_read_tot)
      //
      IO::log << std::setw(15) << "Sum Read";

    if(sum_data.size())
      //
      IO::log << std::setw(15) << "Sum Data";
  }
  IO::log << "\n";
  
  IO::log << IO::log_offset;

  if(zpe_calc_tot)
    //
    IO::log << std::setw(10) << zpe_calc_tot;

  if(zpe_miss_tot)
    //
    IO::log << std::setw(10) << zpe_miss_tot;

  if(zpe_read_tot)
    //
    IO::log << std::setw(15) << zpe_read_tot;

  if(zpe_data.size())
    //
    IO::log << std::setw(15) << zpe_data.size();

  if(temperature > 0.) {
    //
    if(int_calc_tot)
      //
      IO::log << std::setw(10) << int_calc_tot;

    if(int_miss_tot)
      //
      IO::log << std::setw(10) << int_miss_tot;

    if(int_read_tot)
      //
      IO::log << std::setw(15) << int_read_tot;

    if(int_data.size())
      //
      IO::log << std::setw(15) << int_data.size();

    if(sum_calc_tot)
      //
      IO::log << std::setw(10) << sum_calc_tot;

    if(sum_miss_tot)
      //
      IO::log << std::setw(10) << sum_miss_tot;

    if(sum_read_tot)
      //
      IO::log << std::setw(15) << sum_read_tot;

    if(sum_data.size())
      //
      IO::log << std::setw(15) << sum_data.size();
  }
  IO::log << "\n\n";
  
  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "anharmonic correction:\n";
  }
  else {
    //
    IO::log << IO::log_offset << "zero-point energy anharmonic correction, 1/cm:\n";
  }

  IO::log << IO::log_offset 
	  << std::setw(2) << "BO" 
	  << std::setw(15) << "Value"
	  << "\n";

  std::map<int, double> res;

  double curr_val = 0.;
  //
  for(std::map<int, double>::const_iterator mit = corr.begin(); mit != corr.end(); ++mit) {
    //
    curr_val += mit->second;

    if(temperature > 0.) {
      //
      res[mit->first] = std::exp(curr_val);
    }
    else {
      //
      res[mit->first] = curr_val;
    }

    dtemp = curr_val;
    //
    if(temperature < 0.) {
      //
      dtemp /= Phys_const::incm;
    }

    IO::log << IO::log_offset 
	    << std::setw(2) << mit->first 
	    << std::setw(15) << dtemp 
	    << "\n";
  }

  IO::log << "\n";

  return res;
}

std::map<int, double> Graph::Expansion::centroid_correction (double temperature) const
{
  const char funame [] = "Graph::Expansion::centroid_correction: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;

  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "temperature, K = " << temperature / Phys_const::kelv << "\n\n";
  }
  else {
    //
    IO::log << IO::log_offset << "low temperature expansion calculation:" << "\n\n";
    
    for(int i = 0.; i < _red_freq.size(); ++i) {
      //
      if(_red_freq[i] <= 0.) {
	ErrOut err_out;
	err_out << funame << "for low temperature expansion calculation all frequencies should be real";
      }
    }
  }
    
  IO::log << IO::log_offset << std::setw(5) << "#";

  if(temperature > 0.)
    //
    IO::log << std::setw(15) << "Value";
	
  IO::log << std::setw(5) << "V#"
	  << std::setw(5) << "B#"
	  << std::setw(5) << "L#" 
	  << std::setw(7) << "Time"
	  << "   "        << "Graph"
	  << std::endl;
  
  std::vector<double> tanh_factor;
  //
  std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

  int sum_calc_tot = 0;
  int zpe_calc_tot = 0;
  int int_calc_tot = 0;
    
  int sum_miss_tot = 0;
  int zpe_miss_tot = 0;
  int int_miss_tot = 0;
    
  long sum_read_tot = 0;
  long zpe_read_tot = 0;
  long int_read_tot = 0;
    
  _gmap_t sum_data;
  _gmap_t zpe_data;
  _gmap_t int_data;

  std::map<int, double> corr;
  std::map<int, std::map<int, double> > zpe;

#ifndef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_calc_tot, zpe_calc_tot, int_calc_tot, sum_read_tot, zpe_read_tot, int_read_tot, sum_miss_tot, zpe_miss_tot, int_miss_tot) private(itemp, dtemp, btemp) schedule(dynamic)
#endif
        
  for(int gindex = 0; gindex < Graph::size(); ++gindex) {
    //
    Graph::const_iterator graphit = Graph::begin() + gindex;

    std::time_t  start_time = std::time(0);

    const std::vector<std::multiset<int> > vertex_map = graphit->vertex_bond_map();
    const int vertex_size = vertex_map.size();

    double gvalue = 0.;
    std::map<int, double> gze;

    int sum_calc = 0;
    int zpe_calc = 0;
    int int_calc = 0;
    
    int sum_miss = 0;
    int zpe_miss = 0;
    int int_miss = 0;
    
    long sum_read = 0;
    long zpe_read = 0;
    long int_read = 0;
    
    MultiIndexConvert corr_multi_index(graphit->size(), _red_freq_index.size());

#ifdef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_calc, zpe_calc, int_calc, sum_read, zpe_read, int_read, sum_miss, zpe_miss, int_miss, gvalue) private(itemp, dtemp, btemp) schedule(dynamic)
#endif

    for(long corr_li = 0; corr_li < corr_multi_index.size(); ++corr_li) {

      std::vector<int> corrin = corr_multi_index(corr_li);

      double potfac = 1.;

      btemp = false;
      //
      for(std::vector<std::multiset<int> >::const_iterator mit = vertex_map.begin(); mit != vertex_map.end(); ++mit) {
	//
	std::multiset<int> potex_sign;
	  //
	for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
	    //
	  potex_sign.insert(corrin[*it]);

	potex_t::const_iterator pexit = _potex.find(potex_sign);
	
	if(pexit != _potex.end()) {
	  //
	  potfac *= pexit->second;
	}
	else {
	  //
	  btemp = true;
	  //
	  break;
	}
      }
      
      if(btemp) {
	//
	continue;
      }

      potfac /= (double)graphit->symmetry_factor();
   
      // odd number of vertices term has negative sign
      //
      if(vertex_size % 2)
	//
	potfac = -potfac;

      double fvalue = 0.;
      std::map<int, double> fze;
      //
      // centroid correction mask cycle
      //
      for(MultiIndex cmask(graphit->size(), 1); !cmask.end(); ++cmask) {
	//
	double gfactor = potfac;
	//
	int    t_count = 0;
	
	// frequency adapted graph
	//
	FreqGraph mod_graph;
	
	itemp = 0;
	//
	btemp = false;
	//
	for(GenGraph::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {
	  //
	  // normal mode index
	  //
	  int ci = corrin[itemp];

	  // reduced frequency index
	  //
	  int fi = _red_freq_index[ci];

	  // quantum correlator 
	  //
	  if(cmask[itemp]) {
	    //
	    std::set<int> bond;
	    //
	    for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
	      //
	      bond.insert(*it);
	
	    // bond loop
	    //
	    if(bond.size() == 1) {
	      //
	      // centroid-constrained low frequency correlator value
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		gfactor /= 12. * temperature;
	      }
	      // positive temperature correlator value
	      //
	      else if(temperature > 0.) {
		//
		gfactor /= 2. * _red_freq[fi] * tanh_factor[fi];
	      }
	      // zero temperature correlator value
	      //
	      else {
		//
		gfactor /= 2. * _red_freq[fi];
	      }
	    }
	    // add frequency index to the graph
	    //
	    else {
	      //
	      // low frequency
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		mod_graph[bond].insert(-1);
	      }
	      else {
		//
		mod_graph[bond].insert(fi);
	      }
	    }
	  }
	  // low frequency centroid correction incorporated into the correlator
	  //
	  else if(low_freq.find(fi) != low_freq.end()) {
	    //
	    btemp = true;
	    //
	    break;
	  }
	  // centroid correction
	  //
	  else {
	    //
	    if(temperature > 0.) {
	      //
	      gfactor *= temperature;
	    }
	    else
	      //
	      ++t_count;

	    if(_red_freq[fi] > 0.) {
	      //
	      gfactor /= -_red_freq[fi] * _red_freq[fi];
	    }
	    else {
	      //
	      gfactor /=  _red_freq[fi] * _red_freq[fi];
	    }
	  }
	}

	if(btemp)
	  //
	  continue;
	
	itemp = vertex_size - mod_graph.vertex_size();
	//
	if(itemp) {
	  //
	  if(temperature > 0.) {
	    //
	    gfactor /= std::pow(temperature, (double)itemp);
	  }
	  else {
	    //
	    t_count -= itemp;
	  }
	}

	// frequency adapted graph evaluation
	//
	if(mod_graph.size()) {
	  //
	  // graph factorization into connected graphs
	  //
	  _mg_t fac_graph = mod_graph.factorize();

	  // factorized graph cycle
	  //
	  for(_mg_t::const_iterator fgit = fac_graph.begin(); fgit != fac_graph.end(); ++fgit) {
	    //
	    _Convert::vec_t fac_graph_conv;

	    if(mod_flag & KEEP_PERM) {
	      //
	      fac_graph_conv = _convert(fgit->first);
	    }
	    else {
	      //
	      fac_graph_conv = _convert(*fgit->first.perm_pool().begin());
	    }

	    // zero temperature integral evaluation
	    //
	    if(temperature <= 0.) {
	      //
	      t_count -= fgit->second;

#pragma omp critical(zpe_critical)
	      {
		itemp = zpe_data.find(fac_graph_conv) != zpe_data.end();
	      }

	      // read zero temperature integral value from the database
	      //
	      if(itemp) {
		//
		++zpe_read;
		  
#pragma omp critical(zpe_critical)
		{
		  dtemp = zpe_data[fac_graph_conv];
		}
		
		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= dtemp;
	      }
	      // zero temperature integral calculation
	      //
	      else {
		//
		++zpe_calc;
		  
		dtemp = fgit->first.zpe_factor(_red_freq);

		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= dtemp;

		// save zero temperature integral calculation result in the database
		//
#pragma omp critical(zpe_critical)
		{
		  if(zpe_data.find(fac_graph_conv) != zpe_data.end()) {
		    //
		    ++zpe_miss;
		  }
		  else if(mod_flag & KEEP_PERM) {
		    //
		    // permutationally equivalent configurations
		    //
		    std::set<FreqGraph> pool = fgit->first.perm_pool();
		    
		    for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		      zpe_data[_convert(*pit)] = dtemp;
		  }
		  else {
		    //
		    zpe_data[fac_graph_conv] = dtemp;
		  }
		}
		//
		//
	      } // zero temperature integral calculation
	      //
	      //
	    } // zero temperature integral evaluation
	    //
	    //
	    // thermal whole integral evaluation
	    //
	    else {
	      //
#pragma omp critical(int_critical)
	      {
		itemp = int_data.find(fac_graph_conv) != int_data.end();
	      }

	      // read whole integral value from the database
	      //
	      if(itemp) {
		//
		++int_read;

#pragma omp critical(int_critical)
		{
		  dtemp = int_data[fac_graph_conv];
		}
	      
		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= dtemp;
	      }
	      //
	      // whole integral calculation
	      //
	      else {
		//
		++int_calc;
	      
		double int_val = 1.;

		// graph reduction
		//
		_mg_t zpe_graph;
		//
		FreqGraph red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);

		// low temperature graphs cycle
		//
		for(_mg_t::const_iterator zgit = zpe_graph.begin(); zgit != zpe_graph.end(); ++zgit) {
		  //
		  // low temperature integral (zpe factor) evaluation
	      
		  _Convert::vec_t zpe_graph_conv;
	      
		  if(mod_flag & KEEP_PERM) {
		    //
		    zpe_graph_conv = _convert(zgit->first);
		  }
		  else {
		    //
		    zpe_graph_conv = _convert(*zgit->first.perm_pool().begin());
		  }
		  
#pragma omp critical(zpe_critical)
		  {
		    itemp = zpe_data.find(zpe_graph_conv) != zpe_data.end();
		  }

		  // read low temperature integral value from the database
		  if(itemp) {
		    //
		    ++zpe_read;
		  
#pragma omp critical(zpe_critical)
		    {
		      dtemp = zpe_data[zpe_graph_conv];
		    }
		
		    for(int i = 0; i < zgit->second; ++i)
		      //
		      int_val *= dtemp;
		  }
		  // low temperature integral calculation
		  //
		  else {
		    //
		    ++zpe_calc;
		  
		    dtemp = zgit->first.zpe_factor(_red_freq, temperature, tanh_factor);

		    for(int i = 0; i < zgit->second; ++i)
		      //
		      int_val *= dtemp;

		    // save calculation result in the database
		    //
#pragma omp critical(zpe_critical)
		    {
		      if(zpe_data.find(zpe_graph_conv) != zpe_data.end()) {
			//
			++zpe_miss;
		      }
		      else if(mod_flag & KEEP_PERM) {
			//
			std::set<FreqGraph> pool = zgit->first.perm_pool();

			for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
			  //
			  zpe_data[_convert(*pit)] = dtemp;
		      }
		      else {
			//
			zpe_data[zpe_graph_conv] = dtemp;
		      }
		    }
		    //
		    //
		  } // low temperature integral (zpe factor) calculation
		  //
		  //
		} // zpe factor calculation cycle

		// reduced graph fourier sum evaluation
		//
		if(red_graph.size()) {
		  //
		  _Convert::vec_t red_graph_conv;

		  if(mod_flag & KEEP_PERM) {
		    //
		    red_graph_conv = _convert(red_graph);
		  }
		  else {
		    //
		    red_graph_conv = _convert(*red_graph.perm_pool().begin());
		  }

#pragma omp critical(sum_critical)
		  {
		    itemp = sum_data.find(red_graph_conv) != sum_data.end();
		  }

		  // read reduced graph fourier sum from the database
		  //
		  if(itemp) {
		    //
		    ++sum_read;
		
#pragma omp critical(sum_critical)
		    {
		      int_val *= sum_data[red_graph_conv];
		    }
		  }
		  // reduced graph fourier sum calculation
		  //
		  else {
		    //
		    ++sum_calc;
		
		    dtemp = red_graph.fourier_sum(_red_freq, temperature);
		    //
		    int_val *= dtemp;

		    // save calculation result in the database
		    //
#pragma omp critical(sum_critical)
		    {
		      if(sum_data.find(red_graph_conv) != sum_data.end()) {
			//
			++sum_miss;
		      }
		      else if(mod_flag & KEEP_PERM) {
			//
			std::set<FreqGraph> pool = red_graph.perm_pool();

			for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
			  //
			  sum_data[_convert(*pit)] = dtemp;
		      }
		      else {
			//
			sum_data[red_graph_conv] = dtemp;
		      }
		    }
		    //
		    //
		  } // reduced graph fourier sum calculation
		  //
		  //
		} // reduced graph fourier sum evaluation
		//
		else {
		  //
		  int_val /= temperature;
		}

		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= int_val;

		// save whole integral value in the database
		//
#pragma omp critical(int_critical)
		{
		  if(int_data.find(fac_graph_conv) != int_data.end()) {
		    //
		    ++int_miss;
		  }
		  else if(mod_flag & KEEP_PERM) {
		    //
		    std::set<FreqGraph> pool = fgit->first.perm_pool();
	      
		    for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		      //
		      int_data[_convert(*pit)] = int_val;
		  }
		  else {
		    //
		    int_data[fac_graph_conv] = int_val;
		  }
		}
		//
		//
	      } // whole integral calculation 
	      //
	      //
	    } // thermal whole integral evaluation
	    //
	    //
	  } // factorized graph cycle
	  //
	  //
	} // modified graph evaluation
	
	if(temperature > 0.) {
	  //
	  fvalue += gfactor;
	}
	else {
	  //
	  fze[t_count] += gfactor;
	}
	//
	//
      } // centroid correction mask cycle

      if(temperature > 0.) {
	//
	gvalue += fvalue;
      }
      else {
#ifdef INNER_CYCLE_PARALLEL
#pragma omp critical
#endif
	{
	  for(std::map<int, double>::const_iterator fzit = fze.begin(); fzit != fze.end(); ++fzit)
	    //
	    gze[fzit->first] += fzit->second;
	}
      }
    }// normal mode indices cycle

    zpe_calc_tot += zpe_calc;
    int_calc_tot += int_calc;
    sum_calc_tot += sum_calc;
    
    zpe_miss_tot += zpe_miss;
    int_miss_tot += int_miss;
    sum_miss_tot += sum_miss;
    
    zpe_read_tot += zpe_read;
    int_read_tot += int_read;
    sum_read_tot += sum_read;
    
#ifndef INNER_CYCLE_PARALLEL
#pragma omp critical
#endif
    {
      if(temperature > 0.) {
	//
	corr[graphit->size()] += gvalue;
      }
      else
	for(std::map<int, double>::const_iterator gzit = gze.begin(); gzit != gze.end(); ++gzit) {
	  //
	  zpe[graphit->size()][gzit->first] += gzit->second; // opposite to zero-point energy correction sign
	}

      IO::log << IO::log_offset 
	      << std::setw(5) << gindex;

      if(temperature > 0.)
	IO::log << std::setw(15) << gvalue;  
	
      IO::log << std::setw(5) << graphit->vertex_size()
	      << std::setw(5) << graphit->bond_size()
	      << std::setw(5) << graphit->loop_size()
	      << std::setw(7) << std::time(0) - start_time 
	      << "   "        << *graphit 
	      << std::endl;
    }
  } // graph cycle

  IO::log << "\n";

  if(temperature <= 0.) {
    //
    for(std::map<int, std::map<int, double> >::const_iterator zit = zpe.begin(); zit != zpe.end(); ++zit) {
      //
      if(zit->second.size() && zit->second.begin()->first < -1)	{
	ErrOut err_out;
	err_out << funame << "temperature expansion has term of 1/T^" << -zit->second.begin()->first << " order";
      }
    }
  }

  IO::log << IO::log_offset
	  << "Statistics:\n\n";

  IO::log << IO::log_offset;

  if(zpe_calc_tot)
    //
    IO::log << std::setw(10) << "ZPE Calc";

  if(zpe_miss_tot)
    //
    IO::log << std::setw(10) << "ZPE Miss";

  if(zpe_read_tot)
    //
    IO::log << std::setw(15) << "ZPE Read";

  if(zpe_data.size())
    //
    IO::log << std::setw(15) << "ZPE Data";

  if(temperature > 0.) {
    //
    if(int_calc_tot)
      //
      IO::log << std::setw(10) << "Int Calc";

    if(int_miss_tot)
      //
      IO::log<< std::setw(10) << "Int Miss";

    if(int_read_tot)
      //
      IO::log<< std::setw(15) << "Int Read";

    if(int_data.size())
      //
      IO::log << std::setw(15) << "Int Data";

    if(sum_calc_tot)
      //
      IO::log << std::setw(10) << "Sum Calc";

    if(sum_miss_tot)
      //
      IO::log << std::setw(10) << "Sum Miss";

    if(sum_read_tot)
      //
      IO::log << std::setw(15) << "Sum Read";

    if(sum_data.size())
      //
      IO::log << std::setw(15) << "Sum Data";
  }
  IO::log << "\n";
  
  IO::log << IO::log_offset;

  if(zpe_calc_tot)
    //
    IO::log << std::setw(10) << zpe_calc_tot;

  if(zpe_miss_tot)
    //
    IO::log << std::setw(10) << zpe_miss_tot;
  
  if(zpe_read_tot)
    //
    IO::log << std::setw(15) << zpe_read_tot;

  if(zpe_data.size())
    //
    IO::log << std::setw(15) << zpe_data.size();

  if(temperature > 0.) {
    //
    if(int_calc_tot)
      //
      IO::log << std::setw(10) << int_calc_tot;

    if(int_miss_tot)
      //
      IO::log << std::setw(10) << int_miss_tot;

    if(int_read_tot)
      //
      IO::log << std::setw(15) << int_read_tot;

    if(int_data.size())
      //
      IO::log << std::setw(15) << int_data.size();

    if(sum_calc_tot)
      //
      IO::log << std::setw(10) << sum_calc_tot;

    if(sum_miss_tot)
      //
      IO::log << std::setw(10) << sum_miss_tot;

    if(sum_read_tot)
      //
      IO::log << std::setw(15) << sum_read_tot;

    if(sum_data.size())
      //
      IO::log << std::setw(15) << sum_data.size();
  }
  IO::log << "\n\n";
  
  std::map<int, double> res;
 
  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "centroid-constrained anharmonic correction:\n";

    IO::log << IO::log_offset 
	    << std::setw(2) << "BO" 
	    << std::setw(15) << "Value" 
	    << "\n";

    double curr_val = 0.;
    //
    for(std::map<int, double>::const_iterator mit = corr.begin(); mit != corr.end(); ++mit) {
      //
      curr_val += mit->second;

      res[mit->first] = std::exp(curr_val);

      IO::log << IO::log_offset 
	      << std::setw(2) << mit->first 
	      << std::setw(15) << curr_val 
	      << "\n";
    }
    IO::log << "\n";
  }
  // low temperature expansion
  //
  else {
    int print_term_max = 3;

    IO::log << IO::log_offset << "centroid-constrained low temperature expansion:\n";

    IO::log << IO::log_offset 
	    << std::setw(2)  << "BO"
	    << std::setw(15) << "1/T(K) term"
	    << std::setw(15) << "T^0 term";

    for(int i = 1; i <= print_term_max; ++i) {
      //
      IO::log << std::setw(9) << "T(K)^" << i << " term";
    }
    IO::log << "\n";

    for(std::map<int, std::map<int, double> >::const_iterator zit = zpe.begin(); zit != zpe.end(); ++zit) {
      //
      IO::log << IO::log_offset 
	      << std::setw(2) << zit->first;
      
      for(std::map<int, double>::const_iterator it = zit->second.begin(); it != zit->second.end(); ++it)
	//
	res[it->first] += it->second;

      for(int i = -1 ; i <= print_term_max; ++i) {
	//
	IO::log << std::setw(15);

	if(res.find(i) == res.end()) {
	  //
	  IO::log << "0";
	}
	else {
	  //
	  dtemp = res[i];

	  switch(i){
	  case 0:
	    break;

	  default:
	    //
	    dtemp *= std::pow(Phys_const::kelv, (double)i);
	    break;
	  }

	  IO::log << dtemp;
	}
      }

      IO::log << "\n";
    }
    IO::log << "\n";
  }

  return res;
}

std::map<int, double> Graph::Expansion::centroid_correction (const std::map<std::multiset<int>, double>& mmat, double temperature) const
{
  const char funame [] = "Graph::Expansion::centroid_correction: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;

  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "temperature(K) = " << temperature / Phys_const::kelv << "\n";
  }
  else {
    //
    IO::log << IO::log_offset << "low temperature expansion calculation:\n";
    
    for(int i = 0.; i < _red_freq.size(); ++i)
      //
      if(_red_freq[i] <= 0.) {
	ErrOut err_out;
	err_out << funame << "for zero-point energy calculations all frequencies should be real";
      }
  }
    
  IO::log << IO::log_offset 
	  << std::setw(5) << "#";

  if(temperature > 0.)
    IO::log << std::setw(15) << "Value";
	
  IO::log << std::setw(5) << "V#"
	  << std::setw(5) << "B#"
	  << std::setw(5) << "L#" 
	  << std::setw(7) << "Time"
	  << "   Graph"   << std::endl;
  
  std::vector<double> tanh_factor;
  //
  std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);
  
  int sum_calc_tot = 0;
  int zpe_calc_tot = 0;
  int int_calc_tot = 0;
    
  int sum_miss_tot = 0;
  int zpe_miss_tot = 0;
  int int_miss_tot = 0;
    
  long sum_read_tot = 0;
  long zpe_read_tot = 0;
  long int_read_tot = 0;
    
  _gmap_t sum_data;
  _gmap_t zpe_data;
  _gmap_t int_data;

  std::map<int, double> corr;
  //
  std::map<int, std::map<int, double> > zpe;

#ifndef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_miss_tot, zpe_miss_tot, int_miss_tot, sum_calc_tot, zpe_calc_tot, int_calc_tot, sum_read_tot, zpe_read_tot, int_read_tot) private(itemp, dtemp, btemp) schedule(dynamic)
#endif
        
  for(int gindex = 0; gindex < Graph::size(); ++gindex) {
    //
    Graph::const_iterator graphit = Graph::begin() + gindex;

    std::time_t  start_time = std::time(0);

    const int vertex_size = graphit->vertex_size();

    std::vector<std::set<int> > vertex_map(vertex_size);
    //
    itemp = 0;
    //
    for(GenGraph::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {
      //
      if(git->size() != 2) {
	std::cerr << funame << "bond should connect exactly two vertices: " << git->size() << "\n";
	throw Error::Logic();
      }
      
      int v = 0;
      //
      for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit, ++v)
	//
	if(!vertex_map[*bit].insert(2 * itemp + v).second) {
	  std::cerr << funame << "duplicated frequency index: " << 2 * itemp + v << "\n";
	  throw Error::Logic();
	}
    }

    double gvalue = 0.;         // thermal graph value
    //
    std::map<int, double> gze;  // graph temperature expansion at low temperatures

    int sum_calc = 0;
    int zpe_calc = 0;
    int int_calc = 0;
    
    int sum_miss = 0;
    int zpe_miss = 0;
    int int_miss = 0;
    
    long sum_read = 0;
    long zpe_read = 0;
    long int_read = 0;
    
    MultiIndexConvert corr_multi_index(graphit->size() * 2, _red_freq_index.size());

#ifdef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_miss, zpe_miss, int_miss, sum_calc, zpe_calc, int_calc, sum_read, zpe_read, int_read, gvalue) private(itemp, dtemp, btemp) schedule(dynamic)
#endif

    for(long corr_li = 0; corr_li < corr_multi_index.size(); ++corr_li) {

      std::vector<int> corrin = corr_multi_index(corr_li);

      double potfac = 1.;

      btemp = false;
      //
      for(std::vector<std::set<int> >::const_iterator mit = vertex_map.begin(); mit != vertex_map.end(); ++mit) {
	//
	std::multiset<int> potex_sign;
	//
	for(std::set<int>::const_iterator it = mit->begin(); it != mit->end(); ++it) {
	    //
	  potex_sign.insert(corrin[*it]);
	}
	      
	potex_t::const_iterator pexit = _potex.find(potex_sign);
	
	if(pexit != _potex.end()) {
	  //
	  potfac *= pexit->second;
	}
	else {
	  //
	  btemp = true;
	  //
	  break;
	}
      }
      
      if(btemp) {
	//
	continue;
      }

      potfac /= (double)graphit->symmetry_factor();
   
      // odd number of vertices term has negative sign
      //
      if(vertex_size % 2)
	//
	potfac = -potfac;

      double fvalue = 0.;
      //
      std::map<int, double> fze;

      // centroid correction mask cycle
      //
      for(MultiIndex cmask(graphit->size(), 1); !cmask.end(); ++cmask) {
	//
	double gfactor = potfac;
	//
	int t_count = 0;

	// frequency adapted graph
	//
	FreqGraph mod_graph;

	itemp = 0;
	//
	btemp = false;
	//
	for(GenGraph::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {
	  //
	  // individual bond normal mode indices
	  //
	  int ci[2];
	  //
	  for(int i = 0; i < 2; ++i)
	    //
	    ci[i] = corrin[2 * itemp + i];
	
	  // quantum correlator
	  //
	  if(cmask[itemp]) {
	    //
	    // cross correlator terms provided by centroid correction only
	    //
	    if(ci[0] != ci[1]) {
	      //
	      btemp = true;
	      //
	      break;
	    }

	    // reduced frequency index
	    //
	    const int fi = _red_freq_index[ci[0]];
	    
	    std::set<int> bond;
	    //
	    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit) {
	      //
	      bond.insert(*bit);
	    }
	
	    // bond loop
	    //
	    if(bond.size() == 1) {
	      //
	      // centroid-constrained low frequency correlator value
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		gfactor /= 12. * temperature;
	      }
	      // positive temperature correlator value
	      //
	      else if(temperature > 0.) {
		//
		gfactor /= 2. * _red_freq[fi] * tanh_factor[fi];
	      }
	      // zero temperature correlator value
	      //
	      else {
		//
		gfactor /= 2. * _red_freq[fi];
	      }
	    }
	    // add frequency index to the graph
	    //
	    else {
	      //
	      // low frequency
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		mod_graph[bond].insert(-1);
	      }
	      else {
		//
		mod_graph[bond].insert(fi);
	      }
	    }
	  }
	  // centroid correction
	  //
	  else {
	    //
	    std::multiset<int> mi;
	    //
	    for(int i = 0; i < 2; ++i)
	      //
	      mi.insert(ci[i]);

	    std::map<std::multiset<int>, double>::const_iterator mit = mmat.find(mi);

	    // correlator cross terms
	    //
	    if(ci[0] != ci[1]) {
	      //
	      if(mmat.end() != mit) {
		//
		gfactor *= -mit->second;
	      }
	      else {
		//
		btemp = true;

		break;
	      }
	    
	      for(int i = 0; i < 2; ++i) {
		//
		const int fi = _red_freq_index[ci[i]]; // reduced frequency index
	      
		if(_red_freq[fi] > 0.) {
		  //
		  gfactor /= _red_freq[fi] * _red_freq[fi];
		}
		else {
		  //
		  gfactor /= -_red_freq[fi] * _red_freq[fi];
		}
	      }
	    }
	    // diagonal terms
	    //
	    else {
	      //
	      const int fi = _red_freq_index[ci[0]]; // reduced frequency index

	      if(_red_freq[fi] > 0.) {
		//
		dtemp = _red_freq[fi] * _red_freq[fi];
	      }
	      else {
		//
		dtemp = -_red_freq[fi] * _red_freq[fi];
	      }

	      // low frequency centroid correction residue
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		if(mmat.end() != mit) {
		  //
		  gfactor *= (1. - mit->second / dtemp) / dtemp;
		}
		else {
		  //
		  ErrOut err_out;

		  err_out << funame << "low frequency " << fi << " does not have counterpart in the M matrix";
		}
	      }
	      // diagonal centroid correction value
	      //
	      else if(mmat.end() != mit) {
		//
		gfactor *=  -mit->second / dtemp / dtemp;
	      }
	      else {
		//
		btemp = true;
		//
		break;
	      }
	    }

	    if(temperature > 0.) {
	      //
	      gfactor *= temperature;
	    }
	    else {
	      //
	      ++t_count;
	    }
	    //
	    //
	  } // centroid correction
	  //
	  //
	} // graph cycle

	if(btemp)
	  //
	  continue;
	
	itemp = vertex_size - mod_graph.vertex_size();
	//
	if(itemp) {
	  //
	  if(temperature > 0.)
	    //
	    gfactor /= std::pow(temperature, (double)itemp);
	  else
	    //
	    t_count -= itemp;
	}

	// frequency adapted graph evaluation
	//
	if(mod_graph.size()) {
	  //
	  // graph factorization into connected graphs
	  //
	  _mg_t fac_graph = mod_graph.factorize();

	  // factorized graph cycle
	  //
	  for(_mg_t::const_iterator fgit = fac_graph.begin(); fgit != fac_graph.end(); ++fgit) {
	    //
	    _Convert::vec_t fac_graph_conv;

	    if(mod_flag & KEEP_PERM) {
	      //
	      fac_graph_conv = _convert(fgit->first);
	    }
	    else {
	      //
	      fac_graph_conv = _convert(*fgit->first.perm_pool().begin());
	    }
	    
	    // zero temperature integral (zpe factor) evaluation
	    //
	    if(temperature <= 0.) {
	      //
	      t_count -= fgit->second;

#pragma omp critical(zpe_critical)
	      {	      
		itemp = zpe_data.find(fac_graph_conv) != zpe_data.end();
	      }

	      // read zero temperature integral (zpe factor) value from the database
	      //
	      if(itemp) {
		//
		++zpe_read;
		  
#pragma omp critical(zpe_critical)
		{
		  dtemp = zpe_data[fac_graph_conv];
		}
		
		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= dtemp;

	      }
	      // zero temperature integral (zpe factor) calculation
	      //
	      else {
		//
		++zpe_calc;
		  
		dtemp = fgit->first.zpe_factor(_red_freq);
		
		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= dtemp;

		// save zero temperature integral (zpe factor) value in the database
		//
#pragma omp critical(zpe_critical)
		{	      
		  if(zpe_data.find(fac_graph_conv) != zpe_data.end()) {
		    //
		    ++zpe_miss;
		  }
		  else if(mod_flag & KEEP_PERM) {
		    //
		    std::set<FreqGraph> pool = fgit->first.perm_pool();
		  
		    for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		      //
		      zpe_data[_convert(*pit)] = dtemp;
		  }
		  else {
		    //
		    zpe_data[fac_graph_conv] = dtemp;
		  }
		}
		//
		//
	      } // zero temperature integral (zpe factor) calculation
	      //
	      //
	    } // zero temperature integral (zpe factor) evaluation
	    //
	    // thermal whole integral evaluation
	    //
	    else {
	      //
#pragma omp critical(int_critical)
	      {	      
		itemp = int_data.find(fac_graph_conv) != int_data.end();
	      }

	      // read whole integral value from the database
	      //
	      if(itemp) {
		//
		++int_read;
	      
#pragma omp critical(int_critical)
		{	      
		  dtemp = int_data[fac_graph_conv];
		}

		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= dtemp;

	      }
	      // whole integral calculation
	      //
	      else {
		//
		++int_calc;
	      
		double int_val = 1.;

		// graph reduction
		//
		_mg_t zpe_graph;
		//
		FreqGraph red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);
	      
		// zpe graph cycle
		//
		for(_mg_t::const_iterator zgit = zpe_graph.begin(); zgit != zpe_graph.end(); ++zgit) {
		  //
		  // low temperature integral (zpe factor) evaluation
		  //
		  _Convert::vec_t zpe_graph_conv;

		  if(mod_flag & KEEP_PERM) {
		    //
		    zpe_graph_conv = _convert(zgit->first);
		  }
		  else {
		    //
		    zpe_graph_conv = _convert(*zgit->first.perm_pool().begin());
		  }
	      
#pragma omp critical(zpe_critical)
		  {	      
		    itemp = zpe_data.find(zpe_graph_conv) != zpe_data.end();
		  }

		  // read low temperature integral (zpe factor) value from the database
		  //
		  if(itemp) {
		    //
		    ++zpe_read;
		  
#pragma omp critical(zpe_critical)
		    {	      
		      dtemp = zpe_data[zpe_graph_conv];
		    }
		
		    for(int i = 0; i < zgit->second; ++i)
		      //
		      int_val *= dtemp;
		  }
		  // low temperature integral (zpe factor) calculation
		  //
		  else {
		    //
		    ++zpe_calc;
		  
		    dtemp = zgit->first.zpe_factor(_red_freq, temperature, tanh_factor);
		
		    for(int i = 0; i < zgit->second; ++i)
		      //
		      int_val *= dtemp;

		    // save low temperature integral (zpe factor) value in the database
		    //
#pragma omp critical(zpe_critical)
		    {	      
		      if(zpe_data.find(zpe_graph_conv) != zpe_data.end()) {
			//
			++zpe_miss;
		      }
		      else if(mod_flag & KEEP_PERM) {
			//
			std::set<FreqGraph> pool = zgit->first.perm_pool();
	  
			for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
			  //
			  zpe_data[_convert(*pit)] = dtemp;
		      }
		      else {
			//
			zpe_data[zpe_graph_conv] = dtemp;
		      }
		    }
		    //
		    //
		  } // low temperature integral (zpe factor) calculation
		  //
		  //
		} // zpe graph cycle

		// reduced graph fourier sum evaluation
		//
		if(red_graph.size()) {
		  //
		  _Convert::vec_t red_graph_conv;

		  if(mod_flag & KEEP_PERM) {
		    //
		    red_graph_conv = _convert(red_graph);
		  }
		  else {
		    //
		    red_graph_conv = _convert(*red_graph.perm_pool().begin());
		  }

#pragma omp critical(sum_critical)
		  {	      
		    itemp = sum_data.find(red_graph_conv) != sum_data.end();
		  }

		  // read reduced graph fourier sum from the database
		  //
		  if(itemp) {
		    //
		    ++sum_read;
		
#pragma omp critical(sum_critical)
		    {	      
		      int_val *= sum_data[red_graph_conv];
		    }
		  }
		  // reduced graph fourier sum calculation
		  //
		  else {
		    //
		    ++sum_calc;
		
		    dtemp = red_graph.fourier_sum(_red_freq, temperature);
		    //
		    int_val *= dtemp;

		    // save reduced graph fourier sum in the database
		    //
#pragma omp critical(sum_critical)
		    {	      
		      if(sum_data.find(red_graph_conv) != sum_data.end()) {
			//
			++sum_miss;
		      }
		      else if(mod_flag & KEEP_PERM) {
			//
			std::set<FreqGraph> pool = red_graph.perm_pool();
		
			for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
			  //
			  sum_data[_convert(*pit)] = dtemp;
		      }
		      else {
			//
			sum_data[red_graph_conv] = dtemp;
		      }
		    }
		    //
		    //
		  } // reduced graph fourier sum calculation
		  //
		  //
		} // reduced graph fourier sum evaluation
		//
		else {
		  //
		  int_val /= temperature;
		}

		for(int i = 0; i < fgit->second; ++i)
		  //
		  gfactor *= int_val;

		// save whole integral value in the database
		//
#pragma omp critical(int_critical)
		{	      
		  if(int_data.find(fac_graph_conv) != int_data.end()) {
		    //
		    ++int_miss;
		  }
		  else if(mod_flag & KEEP_PERM) {
		    //
		    std::set<FreqGraph> pool = fgit->first.perm_pool();

		    for(std::set<FreqGraph>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		      //
		      int_data[_convert(*pit)] = int_val;
		  }
		  else {
		    //
		    int_data[fac_graph_conv] = int_val;
		  }
		}
		//
		//
	      } // whole integral calculation
	      //
	      //
	    } // whole integral evaluation
	    //
	    //
	  } // factorized graph cycle
	  //
	  //
	} // frequency adapted graph evaluation

	if(temperature > 0.) {
	  //
	  fvalue += gfactor;
	}
	else {
	  //
	  fze[t_count] += gfactor;
	}

      } // centroid correction mask cycle

      if(temperature > 0.) {
	//
	gvalue += fvalue;
      }
      else {
	//
#ifdef INNER_CYCLE_PARALLEL
#pragma omp critical
#endif
	{
	  for(std::map<int, double>::const_iterator fzit = fze.begin(); fzit != fze.end(); ++fzit)
	    //
	    gze[fzit->first] += fzit->second;
	}
      }
      //
      //
    } // normal mode indices cycle

    zpe_calc_tot += zpe_calc;
    int_calc_tot += int_calc;
    sum_calc_tot += sum_calc;
    
    zpe_miss_tot += zpe_miss;
    int_miss_tot += int_miss;
    sum_miss_tot += sum_miss;
    
    zpe_read_tot += zpe_read;
    int_read_tot += int_read;
    sum_read_tot += sum_read;
    
#ifndef INNER_CYCLE_PARALLEL
#pragma omp critical
#endif
    {
      if(temperature > 0.) {
	//
	corr[graphit->size()] += gvalue;
      }
      // zero temperature
      //
      else {
	//
	for(std::map<int, double>::const_iterator gzit = gze.begin(); gzit != gze.end(); ++gzit) {
	  //
	  zpe[graphit->size()][gzit->first] += gzit->second; // opposite to zero-point energy sign
	}
      }
      
      IO::log << IO::log_offset << std::setw(5) << gindex;

      if(temperature > 0.)
	IO::log << std::setw(15) << gvalue;  
	
      IO::log << std::setw(5) << graphit->vertex_size()
	      << std::setw(5) << graphit->bond_size()
	      << std::setw(5) << graphit->loop_size()
	      << std::setw(7) << std::time(0) - start_time 
	      << "   "        << *graphit
	      << std::endl;
    }
    //
    //
  } // graph cycle

  IO::log << "\n";

  if(temperature <= 0) {
    //
    for(std::map<int, std::map<int, double> >::const_iterator zit = zpe.begin(); zit != zpe.end(); ++zit) {
      //
      if(zit->second.size() && zit->second.begin()->first < -1)	{
	//
	ErrOut err_out;

	err_out << funame << "temperature expansion has term of 1/T^" << -zit->second.begin()->first << " order";
      }
    }
  }
  
  IO::log << IO::log_offset << "Statistics:\n\n";

  IO::log << IO::log_offset;

  if(zpe_calc_tot)
    //
    IO::log << std::setw(10) << "ZPE Calc";

  if(zpe_miss_tot)
    //
    IO::log << std::setw(10) << "ZPE Miss";

  if(zpe_read_tot)
    //
    IO::log << std::setw(15) << "ZPE Read";

  if(zpe_data.size())
    //
    IO::log << std::setw(15) << "ZPE Data";

  if(temperature > 0.) {
    //
    if(int_calc_tot)
      //
      IO::log << std::setw(10) << "Int Calc";

    if(int_miss_tot)
      //
      IO::log<< std::setw(10) << "Int Miss";

    if(int_read_tot)
      //
      IO::log<< std::setw(15) << "Int Read";

    if(int_data.size())
      //
      IO::log << std::setw(15) << "Int Data";

    if(sum_calc_tot)
      //
      IO::log << std::setw(10) << "Sum Calc";

    if(sum_miss_tot)
      //
      IO::log << std::setw(10) << "Sum Miss";

    if(sum_read_tot)
      //
      IO::log << std::setw(15) << "Sum Read";

    if(sum_data.size())
      //
      IO::log << std::setw(15) << "Sum Data";
  }
  IO::log << "\n";
  
  IO::log << IO::log_offset;

  if(zpe_calc_tot)
    //
    IO::log << std::setw(10) << zpe_calc_tot;

  if(zpe_miss_tot)
    //
    IO::log << std::setw(10) << zpe_miss_tot;

  if(zpe_read_tot)
    //
    IO::log << std::setw(15) << zpe_read_tot;

  if(zpe_data.size())
    //
    IO::log << std::setw(15) << zpe_data.size();

  if(temperature > 0.) {
    //
    if(int_calc_tot)
      //
      IO::log << std::setw(10) << int_calc_tot;

    if(int_miss_tot)
      //
      IO::log << std::setw(10) << int_miss_tot;

    if(int_read_tot)
      //
      IO::log << std::setw(15) << int_read_tot;

    if(int_data.size())
      //
      IO::log << std::setw(15) << int_data.size();

    if(sum_calc_tot)
      //
      IO::log << std::setw(10) << sum_calc_tot;

    if(sum_miss_tot)
      //
      IO::log << std::setw(10) << sum_miss_tot;

    if(sum_read_tot)
      //
      IO::log << std::setw(15) << sum_read_tot;

    if(sum_data.size())
      //
      IO::log << std::setw(15) << sum_data.size();
  }
  IO::log << "\n\n";

  std::map<int, double> res;

  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "centroid-constrained anharmonic correction:\n";

    IO::log << IO::log_offset 
	    << std::setw(2) << "BO" 
	    << std::setw(15) << "Value" 
	    << "\n";

    double curr_val = 0.;
    //
    for(std::map<int, double>::const_iterator mit = corr.begin(); mit != corr.end(); ++mit) {
      //
      curr_val += mit->second;

      res[mit->first] = std::exp(curr_val);

      IO::log << IO::log_offset 
	      << std::setw(2) << mit->first 
	      << std::setw(15) << curr_val 
	      << "\n";
    }
    IO::log << "\n";
  }
  // low temperature expansion
  //
  else {
    //
    int print_term_max = 3;

    IO::log << IO::log_offset << "centroid-constrained low temperature expansion:\n";

    IO::log << IO::log_offset 
	    << std::setw(2)  << "BO"
	    << std::setw(15) << "1/T(K) term"
	    << std::setw(15) << "T^0 term";
    
    for(int i = 1; i <= print_term_max; ++i) {
      //
      IO::log << std::setw(9) << "T(K)^" << i << " term";
    }
    IO::log << "\n";

    for(std::map<int, std::map<int, double> >::const_iterator zit = zpe.begin(); zit != zpe.end(); ++zit) {
      //
      IO::log << IO::log_offset 
	      << std::setw(2) << zit->first;
      
      for(std::map<int, double>::const_iterator it = zit->second.begin(); it != zit->second.end(); ++it) {
	//
	res[it->first] += it->second;
      }

      for(int i = -1 ; i <= print_term_max; ++i) {
	//
	IO::log << std::setw(15);

	if(res.find(i) == res.end()) {
	  //
	  IO::log << "0";
	}
	else {
	  //
	  dtemp = res[i];

	  switch(i){
	    //
	  case 0:

	    break;

	  default:
	    //
	    dtemp *= std::pow(Phys_const::kelv, (double)i);
	    break;
	  }

	  IO::log << dtemp;
	}
      }

      IO::log << "\n";
    }
    
    IO::log << "\n";
  }

  return res;
}

void Graph::Expansion::init (const std::vector<double>& freq, const potex_t& pex)
{
  const char funame [] = "Graph::Expansion::init: ";

  int itemp;

  if(!freq.size()) {
    //
    ErrOut err_out;

    err_out << funame << "no frequencies";
  }
  
  // initialize frequencies
  //
  _set_frequencies(freq);

  // initialize potential expansion
  //
  _potex = pex;

  // maximal number of vertices
  //
  itemp = 2 * Graph::bond_max / 3;

  // initialize frequency adapted graph converter
  //
  _convert.init(itemp, _red_freq.size() + 1);
}

#include "graph_include.cc"

/*************************************************************************************************
 ****************************** FREQUENCY ADAPTED GRAPH CONVERTER ********************************
 *************************************************************************************************/

void Graph::Expansion::_Convert::init (int v, int f)
{
  const char funame [] = "Graph::Expansion::_Convert::init: ";

  int itemp;
  
  _vertex_size_max = v;

  _freq_size       = f;

  if(v <= 0 || f <= 0) {
    //
    ErrOut err_out;

    err_out << funame << "maximum number of vertices and/or number of frequencies out of range: " << v << ", " << f;
  }

  // check if the used integer type can accommodate the graph data
  //
  if(sizeof(int_t) < sizeof(int)) {
    //
    itemp = 1;

    itemp <<= sizeof(int_t) * 8 - 1;
    
    if(itemp < f * v * (v - 1) / 2) {
      //
      ErrOut err_out;

      err_out << funame << "integer type is too small to accommodate graph data";
    }
  }

  _index_map.resize(v * (v - 1) / 2);

  itemp = 0;
  //
  for(int i = 1; i < v; ++i) {
    //
    for(int j = 0; j < i; ++j, ++itemp) {
      //
      _index_map[itemp].insert(i);
   
      _index_map[itemp].insert(j);
    }
  }
}

Graph::FreqGraph Graph::Expansion::_Convert::operator() (const vec_t& gconv) const
{
  const char funame [] = "Graph::Expansion::_Convert::operator(): ";
  
  int itemp;

  if(!_vertex_size_max) {
    //
    ErrOut err_out;

    err_out << funame << "frequency graph data converter not initialized";
  }

  FreqGraph res;

  if(!gconv.size())
    return res;

  for(int i = 0; i < gconv.size(); ++i) {
    //
    if(gconv[i] < 0) {
      //
      ErrOut err_out;

      err_out << funame << "negative index";
    }

    itemp = gconv[i] / _freq_size;
    
    if(itemp >= _vertex_size_max * (_vertex_size_max - 1) / 2) {
      //
      ErrOut err_out;
      //
      err_out << funame << "vertex index convert out of range";
    }
    
    res[_index_map[itemp]].insert(gconv[i] % _freq_size - 1);
  }

  return res;
}
  
Graph::Expansion::_Convert::vec_t Graph::Expansion::_Convert::operator() (const FreqGraph& freq_graph) const
{
  const char funame [] = "Graph::Expansion::_Convert::operator(): ";
  
  int itemp;

  if(!_vertex_size_max) {
    //
    ErrOut err_out;

    err_out << funame << "frequency graph converter not initialized";
  }
  
  if(!freq_graph.size())
    //
    return vec_t();
  
  if(freq_graph.vertex_size() > _vertex_size_max) {
    //
    ErrOut err_out;

    err_out << funame << "number of vertices exceeds the maximum: " << freq_graph.vertex_size();
  }
  
  vec_t res(freq_graph.bond_size());
  
  int bond_index = 0;
  //
  for(FreqGraph::const_iterator git = freq_graph.begin(); git != freq_graph.end(); ++git) {
    //
    itemp  = *git->first.rbegin();
    //
    itemp  = *git->first.begin() + itemp * (itemp - 1) / 2;
    //
    itemp *= _freq_size;
    
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit, ++bond_index) {
      //
      if(*fit < -1 || *fit >= _freq_size - 1) {
	//
	ErrOut err_out;

	err_out << funame << "frequency index out of range: " << *fit;
      }

      res[bond_index] = itemp + *fit + 1;
    }
  }
  return res;  
}

long Graph::Expansion::_gmap_t::mem_size () const
{
  long res = 0;
  //
  for(const_iterator dit = begin(); dit != end(); ++dit)
    //
    res += (long)_Convert::mem_size(dit->first);

  res += (long)size() * long(sizeof(_Convert::vec_t) + 32); // 32 stands for three pointers and the double
  
  return res;
}

