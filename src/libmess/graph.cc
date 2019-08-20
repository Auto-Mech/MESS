#include "graph.hh"
#include "key.hh"

#include <cstdarg>
//#include <omp.h>

//const int max_thread_num = omp_get_max_threads();

int    GraphExpansion::mod_flag        = GraphExpansion::KEEP_PERM; // different modification flags
int    GraphExpansion::bond_max        = 6;    // maximal number of graph bonds (edges)
double GraphExpansion::freq_tol        = 0.1;  // reduced frequency tolerance
double GraphExpansion::four_cut        = 10.;  // fourier sum graph evaluation cutoff
double GraphExpansion::red_thresh      = 10.;  // low temperature / high frequency graph reduction threshold
double GraphExpansion::low_freq_thresh = 0.2;  // low frequency threshold (relative to the temperature)

/********************************************************************************************
 ************************ PERTURBATION THEORY GRAPH EXPANSION *******************************
 ********************************************************************************************/

std::map<int, double> GraphExpansion::correction (double temperature) const
{
  const char funame [] = "GraphExpansion::correction: ";

  IO::Marker funame_marker(funame);

  int    itemp;
  double dtemp;
  bool   btemp;

  if(temperature > 0.) {
    //
    IO::log << IO::log_offset << "temperature(K) = " << temperature / Phys_const::kelv << std::endl;
  }
  else {
    //
    IO::log << IO::log_offset << "zero-point energy correction calculation:" << std::endl;
    
    for(int i = 0.; i < _red_freq.size(); ++i)
      //
      if(_red_freq[i] <= 0.) {
	ErrOut err_out;
	err_out << funame << "for zero-point energy calculations all frequencies should be real";
      }
  }
    
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
        
  for(int gindex = 0; gindex < _sorted_graph.size(); ++gindex) {
    std::vector<_graph_t>::const_iterator graphit = _sorted_graph.begin() + gindex;

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
	
	_potex_t::const_iterator pexit = _potex.find(potex_sign);
	
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
      _fg_t mod_graph;
      
      itemp = 0;
      //
      for(_graph_t::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {
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
	    ++zpe_read;
	  
#pragma omp critical(zpe_critical)
	    {
	      gfactor *= zpe_data[mod_graph_conv];
	    }
	  }
	  // zero temperature integral (zpe factor) calculation
	  //
	  else {
	    ++zpe_calc;

	    dtemp = mod_graph.zpe_factor(_red_freq);
	    gfactor *= dtemp;

	    // save zero temperature integral value in the database
	    //
#pragma omp critical(zpe_critical)
	    {
	      if(zpe_data.find(mod_graph_conv) != zpe_data.end()) {
		++zpe_miss;
	      }
	      else if(mod_flag & KEEP_PERM) {
		std::set<_fg_t> pool = mod_graph.perm_pool();

		for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
		  zpe_data[_convert(*pit)] = dtemp;
	      }
	      else {
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
	      _fg_t red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);

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
		      std::set<_fg_t> pool = zgit->first.perm_pool();

		      for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
		      std::set<_fg_t> pool = red_graph.perm_pool();
		  
		      for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
		  std::set<_fg_t> pool = fgit->first.perm_pool();

		  for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
   
    if(vertex_size % 2)
      gvalue = -gvalue;
      
    if(temperature < 0.)
      gvalue = -gvalue;

#ifndef INNER_CYCLE_PARALLEL
#pragma omp critical
#endif
    {
      corr[graphit->size()] += gvalue;

      IO::log << IO::log_offset 
	      << "graph # = "            << std::setw(5)  << gindex; 

      if(temperature > 0.) 
	IO::log << " graph value = "       << std::setw(15) << gvalue;
      else
	IO::log << " graph value[1/cm] = " << std::setw(15) << gvalue / Phys_const::incm;
      
      IO::log << " vertex # = "          << std::setw(3)  << graphit->vertex_size()
	      << " bond # = "            << std::setw(3)  << graphit->bond_size()
	      << " loop # = "            << std::setw(3)  << graphit->loop_size()
	      << " elapsed time[sec] = " << std::setw(7)  << std::time(0) - start_time; 

      IO::log << "   " << *graphit << std::endl;
    }
  } // graph cycle

  IO::log << IO::log_offset
	  << "Statistics:\n";

  IO::log << IO::log_offset;
  if(zpe_calc_tot)
    IO::log << std::setw(10) << "ZPE Calc";
  if(zpe_miss_tot)
    IO::log << std::setw(10) << "ZPE Miss";
  if(zpe_read_tot)
    IO::log << std::setw(15) << "ZPE Read";

  if(mod_flag & KEEP_PERM && zpe_data.size())
    IO::log << std::setw(15) << "ZPE Data";

  if(temperature > 0.) {
    if(int_calc_tot)
      IO::log << std::setw(10) << "Int Calc";
    if(int_miss_tot)
      IO::log<< std::setw(10) << "Int Miss";
    if(int_read_tot)
      IO::log<< std::setw(15) << "Int Read";

    if(mod_flag & KEEP_PERM && int_data.size())
      IO::log << std::setw(15) << "Int Data";

    if(sum_calc_tot)
      IO::log << std::setw(10) << "Sum Calc";
    if(sum_miss_tot)
      IO::log << std::setw(10) << "Sum Miss";
    if(sum_read_tot)
      IO::log << std::setw(15) << "Sum Read";

    if(mod_flag & KEEP_PERM && sum_data.size())
      IO::log << std::setw(15) << "Sum Data";
  }
  IO::log << "\n";
  
  IO::log << IO::log_offset;
  if(zpe_calc_tot)
    IO::log << std::setw(10) << zpe_calc_tot;
  if(zpe_miss_tot)
    IO::log << std::setw(10) << zpe_miss_tot;
  if(zpe_read_tot)
    IO::log << std::setw(15) << zpe_read_tot;

  if(mod_flag & KEEP_PERM && zpe_data.size())
    IO::log << std::setw(15) << zpe_data.size();

  if(temperature > 0.) {
    if(int_calc_tot)
      IO::log << std::setw(10) << int_calc_tot;
    if(int_miss_tot)
      IO::log << std::setw(10) << int_miss_tot;
    if(int_read_tot)
      IO::log << std::setw(15) << int_read_tot;

    if(mod_flag & KEEP_PERM && int_data.size())
      IO::log << std::setw(15) << int_data.size();

    if(sum_calc_tot)
      IO::log << std::setw(10) << sum_calc_tot;
    if(sum_miss_tot)
      IO::log << std::setw(10) << sum_miss_tot;
    if(sum_read_tot)
      IO::log << std::setw(15) << sum_read_tot;

    if(mod_flag & KEEP_PERM && sum_data.size())
      IO::log << std::setw(15) << sum_data.size();
  }
  IO::log << "\n\n";
  
  if(temperature > 0.)
    IO::log << IO::log_offset << "anharmonic correction:\n";
  else
    IO::log << IO::log_offset << "zero-point energy anharmonic correction, 1/cm:\n";

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

std::map<int, double> GraphExpansion::centroid_correction (double temperature) const
{
  const char funame [] = "GraphExpansion::centroid_correction: ";

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
    
  IO::log << IO::log_offset 
	  << std::setw(5) << "#";

  if(temperature > 0.)
    IO::log << std::setw(15) << "value";
	
  IO::log << std::setw(5) << "V#"
	  << std::setw(5) << "B#"
	  << std::setw(5) << "L#" 
	  << std::setw(7) << "time"
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
  std::map<int, std::map<int, double> > zpe;

#ifndef INNER_CYCLE_PARALLEL
#pragma omp parallel for default(shared) reduction(+: sum_calc_tot, zpe_calc_tot, int_calc_tot, sum_read_tot, zpe_read_tot, int_read_tot, sum_miss_tot, zpe_miss_tot, int_miss_tot) private(itemp, dtemp, btemp) schedule(dynamic)
#endif
        
  for(int gindex = 0; gindex < _sorted_graph.size(); ++gindex) {
    std::vector<_graph_t>::const_iterator graphit = _sorted_graph.begin() + gindex;

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

	_potex_t::const_iterator pexit = _potex.find(potex_sign);
	
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
	_fg_t mod_graph;
	
	itemp = 0;
	//
	btemp = false;
	//
	for(_graph_t::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {
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
		    std::set<_fg_t> pool = fgit->first.perm_pool();
		    
		    for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
		_fg_t red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);

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
			std::set<_fg_t> pool = zgit->first.perm_pool();

			for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
			std::set<_fg_t> pool = red_graph.perm_pool();

			for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
		    std::set<_fg_t> pool = fgit->first.perm_pool();
	      
		    for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
	      << "   " << *graphit << std::endl;
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
	  << "Statistics:\n";

  IO::log << IO::log_offset;
  if(zpe_calc_tot)
    IO::log << std::setw(10) << "ZPE Calc";
  if(zpe_miss_tot)
    IO::log << std::setw(10) << "ZPE Miss";
  if(zpe_read_tot)
    IO::log << std::setw(15) << "ZPE Read";

  if(mod_flag & KEEP_PERM && zpe_data.size())
    IO::log << std::setw(15) << "ZPE Data";

  if(temperature > 0.) {
    if(int_calc_tot)
      IO::log << std::setw(10) << "Int Calc";
    if(int_miss_tot)
      IO::log<< std::setw(10) << "Int Miss";
    if(int_read_tot)
      IO::log<< std::setw(15) << "Int Read";

    if(mod_flag & KEEP_PERM && int_data.size())
      IO::log << std::setw(15) << "Int Data";

    if(sum_calc_tot)
      IO::log << std::setw(10) << "Sum Calc";
    if(sum_miss_tot)
      IO::log << std::setw(10) << "Sum Miss";
    if(sum_read_tot)
      IO::log << std::setw(15) << "Sum Read";

    if(mod_flag & KEEP_PERM && sum_data.size())
      IO::log << std::setw(15) << "Sum Data";
  }
  IO::log << "\n";
  
  IO::log << IO::log_offset;
  if(zpe_calc_tot)
    IO::log << std::setw(10) << zpe_calc_tot;
  if(zpe_miss_tot)
    IO::log << std::setw(10) << zpe_miss_tot;
  if(zpe_read_tot)
    IO::log << std::setw(15) << zpe_read_tot;

  if(mod_flag & KEEP_PERM && zpe_data.size())
    IO::log << std::setw(15) << zpe_data.size();

  if(temperature > 0.) {
    if(int_calc_tot)
      IO::log << std::setw(10) << int_calc_tot;
    if(int_miss_tot)
      IO::log << std::setw(10) << int_miss_tot;
    if(int_read_tot)
      IO::log << std::setw(15) << int_read_tot;

    if(mod_flag & KEEP_PERM && int_data.size())
      IO::log << std::setw(15) << int_data.size();

    if(sum_calc_tot)
      IO::log << std::setw(10) << sum_calc_tot;
    if(sum_miss_tot)
      IO::log << std::setw(10) << sum_miss_tot;
    if(sum_read_tot)
      IO::log << std::setw(15) << sum_read_tot;

    if(mod_flag & KEEP_PERM && sum_data.size())
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

std::map<int, double> GraphExpansion::centroid_correction (const std::map<std::multiset<int>, double>& mmat, double temperature) const
{
  const char funame [] = "GraphExpansion::centroid_correction: ";

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
    IO::log << std::setw(15) << "value";
	
  IO::log << std::setw(5) << "V#"
	  << std::setw(5) << "B#"
	  << std::setw(5) << "L#" 
	  << std::setw(7) << "time"
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
        
  for(int gindex = 0; gindex < _sorted_graph.size(); ++gindex) {
    //
    std::vector<_graph_t>::const_iterator graphit = _sorted_graph.begin() + gindex;

    std::time_t  start_time = std::time(0);

    const int vertex_size = graphit->vertex_size();

    std::vector<std::set<int> > vertex_map(vertex_size);
    //
    itemp = 0;
    //
    for(_graph_t::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {
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
	      
	_potex_t::const_iterator pexit = _potex.find(potex_sign);
	
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
	_fg_t mod_graph;

	itemp = 0;
	//
	btemp = false;
	//
	for(_graph_t::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {
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
		    std::set<_fg_t> pool = fgit->first.perm_pool();
		  
		    for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
		_fg_t red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);
	      
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
			std::set<_fg_t> pool = zgit->first.perm_pool();
	  
			for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
			std::set<_fg_t> pool = red_graph.perm_pool();
		
			for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
		    std::set<_fg_t> pool = fgit->first.perm_pool();

		    for(std::set<_fg_t>::const_iterator pit = pool.begin(); pit != pool.end(); ++pit)
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
      
      IO::log << IO::log_offset 
	      << std::setw(5) << gindex;

      if(temperature > 0.)
	IO::log << std::setw(15) << gvalue;  
	
      IO::log << std::setw(5) << graphit->vertex_size()
	      << std::setw(5) << graphit->bond_size()
	      << std::setw(5) << graphit->loop_size()
	      << std::setw(7) << std::time(0) - start_time 
	      << "   " << *graphit << std::endl;
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
  
  IO::log << IO::log_offset
	  << "Statistics:\n";

  IO::log << IO::log_offset;
  if(zpe_calc_tot)
    IO::log << std::setw(10) << "ZPE Calc";
  if(zpe_miss_tot)
    IO::log << std::setw(10) << "ZPE Miss";
  if(zpe_read_tot)
    IO::log << std::setw(15) << "ZPE Read";

  if(mod_flag & KEEP_PERM && zpe_data.size())
    IO::log << std::setw(15) << "ZPE Data";

  if(temperature > 0.) {
    if(int_calc_tot)
      IO::log << std::setw(10) << "Int Calc";
    if(int_miss_tot)
      IO::log<< std::setw(10) << "Int Miss";
    if(int_read_tot)
      IO::log<< std::setw(15) << "Int Read";

    if(mod_flag & KEEP_PERM && int_data.size())
      IO::log << std::setw(15) << "Int Data";

    if(sum_calc_tot)
      IO::log << std::setw(10) << "Sum Calc";
    if(sum_miss_tot)
      IO::log << std::setw(10) << "Sum Miss";
    if(sum_read_tot)
      IO::log << std::setw(15) << "Sum Read";

    if(mod_flag & KEEP_PERM && sum_data.size())
      IO::log << std::setw(15) << "Sum Data";
  }
  IO::log << "\n";
  
  IO::log << IO::log_offset;
  if(zpe_calc_tot)
    IO::log << std::setw(10) << zpe_calc_tot;
  if(zpe_miss_tot)
    IO::log << std::setw(10) << zpe_miss_tot;
  if(zpe_read_tot)
    IO::log << std::setw(15) << zpe_read_tot;

  if(mod_flag & KEEP_PERM && zpe_data.size())
    IO::log << std::setw(15) << zpe_data.size();

  if(temperature > 0.) {
    if(int_calc_tot)
      IO::log << std::setw(10) << int_calc_tot;
    if(int_miss_tot)
      IO::log << std::setw(10) << int_miss_tot;
    if(int_read_tot)
      IO::log << std::setw(15) << int_read_tot;

    if(mod_flag & KEEP_PERM && int_data.size())
      IO::log << std::setw(15) << int_data.size();

    if(sum_calc_tot)
      IO::log << std::setw(10) << sum_calc_tot;
    if(sum_miss_tot)
      IO::log << std::setw(10) << sum_miss_tot;
    if(sum_read_tot)
      IO::log << std::setw(15) << sum_read_tot;

    if(mod_flag & KEEP_PERM && sum_data.size())
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

// perturbation graph theory initializer
//
void GraphExpansion::init (const std::vector<double>& freq, const _potex_t& pex)
{
  const char funame [] = "GraphExpansion::init: ";

  IO::Marker funame_marker(funame);

  int itemp;

  if(!freq.size()) {
    ErrOut err_out;
    err_out << funame << "no frequencies";
  }
  
  // initialize frequencies
  //
  _set_frequencies(freq);

  // initialize potential expansion
  //
  _potex = pex;

  std::set<int> rank_pool;
  //
  for(_potex_t::const_iterator pit = _potex.begin(); pit != _potex.end(); ++pit)
    //
    rank_pool.insert(pit->first.size());

  std::vector<int> potex_rank;
  for(std::set<int>::const_iterator pit = rank_pool.begin(); pit != rank_pool.end(); ++pit)
    if(*pit > 2)
      potex_rank.push_back(*pit);
  
  if(!potex_rank.size()) {
    ErrOut err_out;
    err_out << funame << "anharmonic terms in the potential expansion do not exist";
  }

  std::vector<int> glimit(potex_rank.size());
  //
  for(int i = 0; i < potex_rank.size(); ++i)
    //
    glimit[i] = 2 * bond_max / potex_rank[i];

  // initialize frequency adapted graph converter (glimit[0] - maximal number of vertices)
  //
  _convert.init(glimit[0], _red_freq.size() + 1);
  
  int unconnected = 0;
  //
  int graph_count = 0;
  //
  int   raw_count = 0;
  
  for(MultiIndex mi(glimit); !mi.end(); ++mi) {
    //
    itemp = 0;
    //
    for(int i = 0; i < potex_rank.size(); ++i)
      //
      itemp += potex_rank[i] * mi.base()[i];

    if(!itemp || itemp % 2 || itemp / 2 > bond_max)
      //
      continue;

    std::map<int, int> ginit;
    //
    for(int i = 0; i < potex_rank.size(); ++i) {
      //
      if(mi.base()[i])
	//
	ginit[potex_rank[i]] = mi.base()[i];
    }

    IO::log << IO::log_offset << "perturbation term:";
    //
    for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
      //
      IO::log << "  " << it->first << "(" << it->second << ")";
    //
    IO::log << std::endl;
	     
    std::vector<int>               vertex_order;
    //
    std::vector<int>               multi_perm_init;
    //
    std::vector<std::vector<int> > vertex_perm_base;

    for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it) {
      //
      // MultiPerm dimensions
      //
      multi_perm_init.push_back(it->second);

      // vertex permutation base & vertex order
      //
      vertex_perm_base.push_back(std::vector<int>(it->second));
      //
      for(int i = 0; i < it->second; ++i) {
	//
	vertex_perm_base.back()[i] = vertex_order.size();
	//
	vertex_order.push_back(it->first);
      }
    }
  
    std::set<_graph_t> raw_graph = _raw_graph_generator(vertex_order);

    raw_count += raw_graph.size();
  
    itemp = 1;
    //
    for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
      //
      itemp *= Math::factorial(it->second);

    const int vertex_perm_size = itemp;

    while(1) {
      //
      // erase not-connected graphs
      //
      while(raw_graph.size() && !raw_graph.begin()->is_connected()) {
	//
	raw_graph.erase(raw_graph.begin());
	//
	++unconnected;
      }

      if(!raw_graph.size())
	//
	break;
    
      ++graph_count;

      _graph_t ancor = *raw_graph.begin();
    
      if(!_graph_data.insert(ancor).second) {
	//
	std::cerr << funame << "graph already in the pool\n";
	//
	throw Error::Logic();
      }

      // remove all permutationally related graphs from the raw graph pool
      //
      int perm_equal = 0;
      //
      int perm_diff  = 0;
      //
      for(MultiPerm multi_perm(multi_perm_init); !multi_perm.end(); ++multi_perm) {
	//
	std::vector<int> vertex_perm(vertex_order.size());
	//
	for(int i = 0; i < vertex_perm_base.size(); ++i) {
	  //
	  for(int j = 0; j < vertex_perm_base[i].size(); ++j) {
	    //
	    vertex_perm[vertex_perm_base[i][j]] = vertex_perm_base[i][multi_perm[i][j]];
	  }
	}
      
	_graph_t perm_graph;
	//
	for(_graph_t::const_iterator pit = ancor.begin(); pit != ancor.end(); ++pit) {
	  //
	  std::multiset<int> bond;
	  //
	  for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it)
	    //
	    bond.insert(vertex_perm[*it]);
	  
	  perm_graph.insert(bond);
	}

	if(perm_graph == ancor)
	  //
	  ++perm_equal;

	std::set<_graph_t>::iterator git = raw_graph.find(perm_graph);
	//
	if(git != raw_graph.end()) {
	  //
	  ++perm_diff;
	  //
	  raw_graph.erase(git);
	}
      }
    
      if(perm_diff * perm_equal != vertex_perm_size) {
	//
	ErrOut err_out;
	//
	err_out << funame << "numbers of permutationally equal and different graphs inconsistent: "
		<< perm_equal << ", "<< perm_diff;
      }

      if(perm_equal != ancor.vertex_symmetry()) {
	//
	ErrOut err_out;
	//
	err_out << funame << "permutational symmetry factors differ";
      }
    }
  }
  IO::log << "\n";
  
  // sorting graph according to certain criterion (vertex size, for example)
  //
  std::multimap<int, _graph_t> sort_map;
  //
  for(std::set<_graph_t>::const_iterator graphit = _graph_data.begin(); graphit != _graph_data.end(); ++graphit)
    //
    sort_map.insert(make_pair(graphit->vertex_size(), *graphit));

  _sorted_graph.resize(sort_map.size());
    
  itemp = 0;
  //
  for(std::multimap<int, _graph_t>::const_reverse_iterator sit = sort_map.rbegin(); sit != sort_map.rend(); ++sit, ++itemp)
    //
    _sorted_graph[itemp] = sit->second;

  IO::log << IO::log_offset << "number of raw graphs = " << raw_count << "\n";
  //
  IO::log << IO::log_offset << "number of unconnected graphs = " << unconnected << "\n";
  //
  IO::log << IO::log_offset << "number of permutationally distinct graphs = " << graph_count << "\n\n";

  IO::log << IO::log_offset << "permutationally distinct graphs(" << _graph_data.size() << "):\n";
  IO::log << IO::log_offset
	  << std::setw(5) << "#"
	  << std::setw(5) << "VS"
	  << std::setw(5) << "BS"
	  << std::setw(5) << "V#"
	  << std::setw(5) << "B#"
	  << std::setw(3) << " "
	  << "graph" << "\n";

  itemp = 0;
  for(std::vector<_graph_t>::const_iterator it = _sorted_graph.begin(); it != _sorted_graph.end(); ++it, ++itemp)
    IO::log << IO::log_offset
	    << std::setw(5) << itemp
	    << std::setw(5) << it->vertex_symmetry()
	    << std::setw(5) << it->bond_symmetry()
	    << std::setw(5) << it->vertex_size()
	    << std::setw(5) << it->bond_size()
	    << std::setw(3) << " "
	    << *it << "\n";
  
  IO::log << std::endl;
}

std::set<GraphExpansion::_graph_t> GraphExpansion::_raw_graph_generator (std::vector<int> vertex_order, int root)
{
  const char funame [] = "GraphExpansion::_raw_graph_generator: ";
  
  int itemp;

  if(::sum(vertex_order) % 2) {
    //
    ErrOut err_out;
    //
    err_out << funame << "odd number of connction points";
  }
  
  std::set<_graph_t> res;

  if(root < 0) {
    //
    std::set<_graph_t> add;
    
    for(int i = 0; i < vertex_order.size(); ++i) {
      //
      if(vertex_order[i] > 0) {
	//
	add = _raw_graph_generator(vertex_order, i);
	//
	root = i;
	//
	break;
      }
    }

    if(!add.size())
      //
      return res;

    for(std::set<_graph_t>::const_iterator at = add.begin(); at != add.end(); ++at) {
      //
      std::vector<int> new_order = vertex_order;
      //
      for(_graph_t::const_iterator pit = at->begin(); pit != at->end(); ++pit) {
	//
	for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it)
	  //
	  --new_order[*it];
      }

      for(int i = 0; i < new_order.size(); ++i) {
	//
	if(new_order[i] < 0 || i <= root && new_order[i]) {
	  //
	  ErrOut err_out;
	  //
	  err_out << funame << "wrong number of connection points";
	}
      }
      
      std::set<_graph_t> g  = _raw_graph_generator(new_order);

      if(g.size()) {
	//
	for(std::set<_graph_t>::const_iterator git = g.begin(); git != g.end(); ++git) {
	  //
	  _graph_t gval = *git;
	  //
	  for(_graph_t::const_iterator pit = at->begin(); pit != at->end(); ++pit)
	    //
	    gval.insert(*pit);
	  
	  res.insert(gval);
	}
      }
      else
	//
	res.insert(*at);
    }
  }
  // front addition
  //
  else {
    //
    if(!vertex_order[root])
      return res;

    --vertex_order[root];
  
    for(int i = 0; i < vertex_order.size(); ++i) {
      //
      if(vertex_order[i] < 0 || i < root && vertex_order[i]) {
	//
	ErrOut err_out;
	//
	err_out << funame << "wrong number of connection points";
      }
      
      if(vertex_order[i] > 0) {
	//
	--vertex_order[i];
	
	std::set<_graph_t> g  = _raw_graph_generator(vertex_order, root);
	
	++vertex_order[i];

	std::multiset<int> bond;
	//
	bond.insert(root);
	//
	bond.insert(i);

	if(g.size()) {
	  //
	  for(std::set<_graph_t>::const_iterator git = g.begin(); git != g.end(); ++git) {
	    //
	    _graph_t gval = *git;
	    //
	    gval.insert(bond);
	    
	    res.insert(gval);
	  }
	}
	else {
	  //
	  _graph_t gval;
	  //
	  gval.insert(bond);
	  
	  res.insert(gval);
	} 
      }
    }
  }
  
  return res;
}
  
std::set<int> GraphExpansion::_low_freq_set (double temperature, std::vector<double>& tanh_factor) const
{
  const char funame [] = "GraphExpansion::_low_freq_set: ";

  std::set<int> res;

  if(temperature > 0.) {
    //
    tanh_factor.resize(_red_freq.size());

    for(int i = 0; i < _red_freq.size(); ++i)
      //
      if(_red_freq[i] > 0.) {
	//
	tanh_factor[i] = std::tanh(_red_freq[i] / temperature / 2.);
      
	if(_red_freq[i] / temperature < low_freq_thresh)
	  //
	  res.insert(i);
      }
      else {
	//
	tanh_factor[i] = std::tan(-_red_freq[i] / temperature / 2.);
	
	if(_red_freq[i] / temperature > -low_freq_thresh)
	  //
	  res.insert(i);
      }
  }
  // zero temperature
  //
  else {
    //
    tanh_factor.clear();
  }

  return res;
}

void GraphExpansion::_set_frequencies (std::vector<double> freq)
{
  const char funame [] = "GraphExpansion::_set_frequencies: ";
  
  static const double min_freq = Phys_const::incm;
  
  double dtemp;
  int    itemp;

  _red_freq_map.clear();

  // removing very low frequencies
  for(int f = 0; f < freq.size(); ++f) {
    //
    if(freq[f] >= 0. && freq[f] < min_freq) {
      //
      freq[f] =  min_freq;
    }
    else if(freq[f] <= 0. && freq[f] > -min_freq) {
      //
      freq[f] = -min_freq;
    }
  }
  
  for(int f = 0; f < freq.size(); ++f) {
    //
    itemp = 1;
    //
    for(int i = 0; i < _red_freq_map.size(); ++i) {
      //
      dtemp = freq[f] / freq[*_red_freq_map[i].begin()];

      if(1. - freq_tol < dtemp && dtemp < 1. + freq_tol) {
	//
	_red_freq_map[i].insert(f);

	itemp = 0;
	
	break;
      }
    }
    
    if(itemp) {
      //
      _red_freq_map.push_back(std::set<int>());
      //
      _red_freq_map.back().insert(f);
    }
  }

  _red_freq.resize(_red_freq_map.size());
  //
  for(int i = 0; i < _red_freq_map.size(); ++i) {
    //
    dtemp = 0.;
    
    itemp =  1;
    //
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it) {
      //
      if(it == _red_freq_map[i].begin() && freq[*it] < 0.)
	//
	itemp = 0;

      if(itemp && freq[*it] < 0. || !itemp && freq[*it] > 0.) {
	//
	ErrOut err_out;
	//
	err_out << funame << "frequencies of different signs in one group";
      }
      
      if(itemp) {
	//
	dtemp += std::log(freq[*it]);
      }
      else {
	//
	dtemp += std::log(-freq[*it]);
      }
    }
    
    dtemp /= (double)_red_freq_map[i].size();

    if(itemp) {
      //
      _red_freq[i] =  std::exp(dtemp);
    }
    else {
      //
      _red_freq[i] = -std::exp(dtemp);
    }
  }

  _red_freq_index.resize(freq.size());
  //
  for(int i = 0; i < _red_freq_map.size(); ++i) {
    //
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it)
      //
      _red_freq_index[*it] = i;
  }
  
  IO::log << IO::log_offset << "logarithmic frequency tolerance  = " << freq_tol << "\n"
	  << IO::log_offset << "reduced frequencies, 1/cm:";
  
  for(int i = 0; i < _red_freq.size(); ++i) {
    //
    IO::log << "   " << _red_freq[i] / Phys_const::incm;
  }
  
  IO::log << std::endl;

  IO::log << IO::log_offset << "reduced frequencies groups: ";
  //
  for(int i = 0; i < _red_freq_map.size(); ++i) {
    //
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it) {
      //
      if(it != _red_freq_map[i].begin())
	IO::log << ", ";
      else
	IO::log << "(";
      
      IO::log << *it;
    }
    
    IO::log << ")";
  }
  
  IO::log << "\n\n";
}

/*************************************************************************************************
 ********************************************* GRAPH *********************************************
 *************************************************************************************************/

std::vector<std::multiset<int> >  GraphExpansion::_graph_t::vertex_bond_map () const
{
  const char funame [] = "GraphExpansion::_graph_t::vertex_bond_map: ";
  
  _check();

  std::vector<std::multiset<int> > res(vertex_size());
  
  int bi = 0;
  for(const_iterator git = begin(); git != end(); ++git, ++bi)
    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
      res[*bit].insert(bi);

  return res;
}

bool GraphExpansion::_graph_t::is_connected () const
{
  const char funame [] = "GraphExpansion::_graph_t::is_connected: ";

  if(!size())
    return true;
  
  _check();

  std::set<int> pool, new_pool;
  new_pool.insert(*begin()->begin());

  while(pool.size() != new_pool.size()) {
    pool = new_pool;
    for(const_iterator git = begin(); git != end(); ++git)
      for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
	if(new_pool.find(*bit) != new_pool.end()) {
	  if(bit == git->begin())
	    new_pool.insert(*git->rbegin());
	  else
	    new_pool.insert(*git->begin());

	  break;
	}
  }

  if(pool.size() == vertex_size())
    return true;

  return false;
}

int GraphExpansion::_graph_t::vertex_symmetry () const
{
  const char funame [] = "GraphExpansion::_graph_t::vertex_symmetry: ";

  int itemp;

  if(!size())
    return 1;

  _check();
  
  std::vector<int> vertex_order(vertex_size());

  for(const_iterator git = begin(); git != end(); ++git)
    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
      ++vertex_order[*bit];

  std::map<int, std::vector<int> > order_map;
  for(int i = 0; i < vertex_order.size(); ++i)
    order_map[vertex_order[i]].push_back(i);

  if(order_map.begin()->first < 3) {
    ErrOut err_out;
    err_out << funame << "vertex order less then 3: " << order_map.begin()->first;
  }

  itemp = 0;
  std::vector<int> order_size(order_map.size());
  for(std::map<int, std::vector<int> >::const_iterator it = order_map.begin(); it != order_map.end(); ++it, ++itemp)
    order_size[itemp] = it->second.size();

  int res = 0;

  // vertex permutation symmetry factor
  for(MultiPerm multi_perm(order_size); !multi_perm.end(); ++multi_perm) {
    itemp = 0;
    std::vector<int> vertex_perm(vertex_order.size());
    for(std::map<int, std::vector<int> >::const_iterator it = order_map.begin(); it != order_map.end(); ++it, ++itemp)
      for(int i = 0; i < it->second.size(); ++i)
	vertex_perm[it->second[i]] = it->second[multi_perm[itemp][i]];

    _graph_t new_graph;

    for(const_iterator git = begin(); git != end(); ++git) {
      std::multiset<int> bond;
      for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
	bond.insert(vertex_perm[*bit]);
      new_graph.insert(bond);
    }
    
    if(new_graph == *this)
      ++res;
  }

  return res;
}

int GraphExpansion::_graph_t::bond_symmetry () const
{
  const char funame [] = "GraphExpansion::_graph_t::bond_symmetry: ";
  
  int itemp;

  if(!size())
    return 1;

  int res = 1;
  
  // loop symmetry
  res <<= loop_size();

  // high order bond symmetry
  std::vector<int> bov = bond_order();
  for(int i = 0; i < bov.size(); ++i)
    res *= Math::factorial(bov[i]);
  
  return res;
}

std::vector<int> GraphExpansion::_graph_t::bond_order () const
{
  int itemp;

  std::map<std::multiset<int>, int> bond_order_map;
  for(const_iterator git = begin(); git != end(); ++git)
    ++bond_order_map[*git];

  std::vector<int> res(bond_order_map.size());

  itemp = 0;
  for(std::map<std::multiset<int>, int>::const_iterator it = bond_order_map.begin(); it != bond_order_map.end(); ++it, ++itemp)
    res[itemp] = it->second;

  return res;
}

void GraphExpansion::_graph_t::_check () const
{
  const char funame [] = "GraphExpansion::_graph_t::_check: ";

  if(!size())
    return;

  std::set<int> vertex_pool;
   for(const_iterator git = begin(); git != end(); ++git) {
    if(git->size() != 2) {
      ErrOut err_out;
      err_out << funame << "there should be exactly two points for a bond";
    }

    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
      vertex_pool.insert(*bit);
  }

  if(*vertex_pool.begin() || *vertex_pool.rbegin() >= vertex_pool.size()) {
    ErrOut err_out;
    err_out << funame << "vertex indices not ordered";
  }
}

int GraphExpansion::_graph_t::vertex_size () const
{
  const char funame [] = "GraphExpansion::_graph_t::vertex_size: ";

  if(!size())
    return 0;

  std::set<int> res;
  for(const_iterator mit = begin(); mit != end(); ++mit) {
    if(mit->size() != 2) {
      ErrOut err_out;
      err_out << funame << "there should be exactly two points for a bond";
    }

    for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
      res.insert(*it);
  }

  if(*res.begin() || *res.rbegin() >= res.size()) {
    ErrOut err_out;
    err_out << funame << "vertex indices not ordered";
  }

  return res.size();
}

int GraphExpansion::_graph_t::bond_size () const
{
  const char funame [] = "GraphExpansion::_graph_t::bond_size: ";

  int res = 0;
  for(const_iterator git = begin(); git != end(); ++git) {
    if(git->size() != 2) {
      ErrOut err_out;
      err_out << funame << "there should be exactly two points for a bond";
    }

    if(*git->begin() != *git->rbegin())
      ++res;
  }

  return res;
}

int GraphExpansion::_graph_t::loop_size () const
{
  const char funame [] = "GraphExpansion::_graph_t::bond_size: ";

  int res = 0;
  for(const_iterator git = begin(); git != end(); ++git) {
    if(git->size() != 2) {
      ErrOut err_out;
      err_out << funame << "there should be exactly two points for a bond";
    }

    if(*git->begin() == *git->rbegin())
      ++res;
  }

  return res;
}

void GraphExpansion::_graph_t::print (std::ostream& to) const
{
  for(const_iterator pit = begin(); pit != end(); ++pit) {
    for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it) {
      if(it == pit->begin())
	to << "(";
      else
	to << ", ";

      to << *it;
    }

    to << ")";
  }
}

/*************************************************************************************************
 ****************************** FREQUENCY ADAPTED GRAPH CONVERTER ********************************
 *************************************************************************************************/

void GraphExpansion::_Convert::init (int v, int f)
{
  const char funame [] = "GraphExpansion::_Convert::init: ";

  int itemp;
  
  _vertex_size_max = v;
  _freq_size       = f;

  if(v <= 0 || f <= 0) {
    ErrOut err_out;
    err_out << funame << "maximum number of vertices and/or number of frequencies out of range: " << v << ", " << f;
  }

  // check if the used integer type can accommodate the graph data
  if(sizeof(int_t) < sizeof(int)) {
  
    itemp = 1;
    itemp <<= sizeof(int_t) * 8 - 1;
    
    if(itemp < f * v * (v - 1) / 2) {
      ErrOut err_out;
      err_out << funame << "integer type is too small to accommodate graph data";
    }
  }
}

GraphExpansion::_Convert::vec_t GraphExpansion::_Convert::operator() (const _fg_t& freq_graph) const
{
  const char funame [] = "GraphExpansion::_Convert::operator(): ";
  
  int itemp;

  if(!_vertex_size_max) {
    ErrOut err_out;
    err_out << funame << "frequency graph data converter not initialized";
  }
  
  if(!freq_graph.size())
    return vec_t();
  
  if(freq_graph.vertex_size() > _vertex_size_max) {
    ErrOut err_out;
    err_out << funame << "number of vertices exceeds the maximum: " << freq_graph.vertex_size();
  }
  
  vec_t res(freq_graph.bond_size());
  
  int bond_index = 0;
  for(_fg_t::const_iterator git = freq_graph.begin(); git != freq_graph.end(); ++git) {
    
    itemp  = *git->first.rbegin();
    //
    itemp  = *git->first.begin() + itemp * (itemp - 1) / 2;
    //
    itemp *= _freq_size;
    
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit, ++bond_index) {
      //
      if(*fit < -1 || *fit >= _freq_size - 1) {
	ErrOut err_out;
	err_out << funame << "frequency index out of range: " << *fit;
      }

      res[bond_index] = itemp + *fit + 1;
    }
  }
  return res;  
}

long GraphExpansion::_gmap_t::mem_size () const
{
  long res = 0;
  for(const_iterator dit = begin(); dit != end(); ++dit)
    res += (long)_Convert::mem_size(dit->first);

  res += (long)size() * long(sizeof(_Convert::vec_t) + 32); // 32 stands for three pointers and the double
  
  return res;
}

/*************************************************************************************************
 ****************************** FREQUENCY ADAPTED GRAPH ******************************************
 *************************************************************************************************/

int GraphExpansion::_fg_t::bond_size () const
{
  _check_integrity();
  
  int res = 0;
  for(const_iterator git = begin(); git != end(); ++git)
    res += git->second.size();

  return res;
}

// permutationally equivalent graphs
//
std::set<GraphExpansion::_fg_t> GraphExpansion::_fg_t::perm_pool (int* symm, int flag) const
{
  const char funame [] = "GraphExpansion::_fg_t::perm_pool: ";

  _check_integrity();

  int itemp;

  std::set<_fg_t> res;
  if(symm)
    *symm = 0;
  
  if(!size())
    return res;
  
  std::set<int> vertex_pool;
  for(_fg_t::const_iterator git = begin(); git != end(); ++git)
    for(std::set<int>::const_iterator it = git->first.begin(); it != git->first.end(); ++it)
      vertex_pool.insert(*it);

  // arbitrary indices
  //
  if(flag) {
    std::vector<int> vertex;
    for(std::set<int>::const_iterator it = vertex_pool.begin(); it != vertex_pool.end(); ++it)
      vertex.push_back(*it);
    
    for(Permutation perm(vertex_pool.size()); !perm.end(); ++perm) {

      _fg_t perm_graph;

      std::vector<int> perm_vertex = perm(vertex);
      
      std::map<int, int> vertex_map;
      for(int i = 0; i < vertex.size(); ++i)
	vertex_map[vertex[i]] = perm_vertex[i];
      
      for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
	
	std::set<int> bond;
	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	  bond.insert(vertex_map[*bit]);
	      
	perm_graph[bond] = git->second;
      }
      
      if(symm && perm_graph == *this)
	++(*symm);
      
      res.insert(perm_graph);
    }
  }
  // ordered vertices
  //
  else {
    if(*vertex_pool.begin() || *vertex_pool.rbegin() >= vertex_pool.size()) {
      ErrOut err_out;
      err_out << funame << "graph not in the standard form";
    }
    
    for(Permutation perm(vertex_pool.size()); !perm.end(); ++perm) {

      _fg_t perm_graph;
      for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
	std::set<int> bond;
	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	  bond.insert(perm[*bit]);
	      
	perm_graph[bond] = git->second;
      }

      if(symm && perm_graph == *this)
	++(*symm);
      
      res.insert(perm_graph);
    }
  }
  
  return res;
}

void GraphExpansion::_fg_t::_check_order () const
{
  const char funame [] = "GraphExpansion::_fg_t::_check_order: ";

  if(!size())
    return;
  
  std::set<int> vertex_pool;
  for(_fg_t::const_iterator git = begin(); git != end(); ++git)
    for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
      vertex_pool.insert(*bit);

  if(!vertex_pool.size() || *vertex_pool.begin() || *vertex_pool.rbegin() >= vertex_pool.size()) {
    ErrOut err_out;
    err_out << funame << "graph not in a standard form";
  }
}

void GraphExpansion::_fg_t::_check_integrity (int freq_size) const
{
  const char funame [] = "GraphExpansion::_fg_t::_check_integrity: ";

  for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
    //
    if(git->first.size() != 2) {
      //
      ErrOut err_out;
      //
      err_out << funame << "bond should have exactly two different ends:";
      //
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	//
	err_out << " " << *bit;
    }
    
    if(!git->second.size()) {
      //
      ErrOut err_out;
      //
      err_out << funame << "zero-order bond";
    }

    if(freq_size > 0) {
      //
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
	//
	if(*fit >= freq_size || *fit < -1) {
	  //
	  ErrOut err_out;
	  //
	  err_out << funame << "frequency index out of order: " << *fit;
	}
      }
    }
  }
}
 
int GraphExpansion::_fg_t::vertex_size () const
{
  const char funame [] = "GraphExpansion::_fg_t::vertex_size: ";

  _check_integrity();

  std::set<int> res;
  for(_fg_t::const_iterator git = begin(); git != end(); ++git) 
    for(std::set<int>::const_iterator vit = git->first.begin(); vit != git->first.end(); ++vit)
      res.insert(*vit);

  return res.size();
}

// factorize the general graph into a set of connected ones (there may be identical ones) 
//
GraphExpansion::_mg_t  GraphExpansion::_fg_t::factorize () const
{
  const char funame [] = "GraphExpansion::_fg_t::factorize: ";

  _check_integrity();

  int itemp;

  std::vector<_fg_t> fac;
  for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
    std::set<int> connect;
    for(std::set<int>::const_iterator vit = git->first.begin(); vit != git->first.end(); ++vit)
      for(int i = 0; i < fac.size(); ++i) {
	itemp = 0;
	for(_fg_t::const_iterator mit = fac[i].begin(); mit != fac[i].end(); ++mit)
	  if(mit->first.find(*vit) != mit->first.end()) {
	    connect.insert(i);
	    itemp = 1;
	    break;
	  }
	if(itemp)
	  break;
      }
    
    switch(connect.size()) {
    case 0:
      fac.push_back(_fg_t());
      fac.back()[git->first] = git->second;
      break;

    case 1:
      fac[*connect.begin()][git->first] = git->second;
      break;

    case 2:
      fac[*connect.begin()][git->first] = git->second;

      for(_fg_t::const_iterator mit = fac[*connect.rbegin()].begin(); mit != fac[*connect.rbegin()].end(); ++mit)
	fac[*connect.begin()][mit->first] = mit->second;
      
      fac.erase(fac.begin() + *connect.rbegin());
      break;

    default:
      std::cerr << funame << "wrong connect size: " << connect.size() << "\n";
      throw Error::Logic();
    }
  }
  
  _mg_t res;
  for(int g = 0; g < fac.size(); ++g) {
    std::set<int> vset;
    for(_fg_t::const_iterator git = fac[g].begin(); git != fac[g].end(); ++git)
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	vset.insert(*bit);

    std::map<int, int> vmap;
    itemp = 0;
    for(std::set<int>::const_iterator it = vset.begin(); it != vset.end(); ++it, ++itemp)
      vmap[*it] = itemp;
 
    _fg_t perm_graph;
    for(_fg_t::const_iterator git = fac[g].begin(); git != fac[g].end(); ++git) {
      std::set<int> bond;
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	bond.insert(vmap[*bit]);
      
      perm_graph[bond] = git->second;
    }
    
    ++res[perm_graph];
  }

  return res;
}

GraphExpansion::_fg_t  GraphExpansion::_fg_t::reduce (const std::vector<double>& freq,
						      double                     temperature,  
						      const std::vector<double>& tanh_factor,
						      _mg_t&                     zpe_graph
						      ) const
{
  const char funame [] = "GraphExpansion::_fg_t::reduce: ";

  _check_integrity();
  //
  _check_order();

  int    itemp;
  double dtemp;

  const int vsize = vertex_size();

  // reduced graph
  //
  _fg_t red_graph;

  std::vector<std::set<int> > red_group(vsize);
  //
  for(int v = 0; v < vsize; ++v)
    //
    red_group[v].insert(v);

  while(1) {// reduction loop
    //
    // new effective graph
    //
    red_graph.clear();

    for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
      //
      std::set<int> bond;
      //
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit) {
	//
	for(int v = 0; v < red_group.size(); ++v) {
	  //
	  if(red_group[v].find(*bit) != red_group[v].end()) {
	    //
	    bond.insert(v);
	    //
	    break;
	  }
	}
      }
      
      if(bond.size() == 2) {
	//
	for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit)
	  //
	  red_graph[bond].insert(*fit);
      }
    }

    // which groups to merge
    std::vector<std::set<int> > merge_group;
    //
    for(_fg_t::const_iterator git = red_graph.begin(); git != red_graph.end(); ++git) {
      //
      dtemp = 0.;
      //
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
	//
	if(*fit >= 0 && freq[*fit] > 0.)
	  //
	  dtemp += freq[*fit] * tanh_factor[*fit];
      }
      
      if(dtemp > temperature * red_thresh) {
	//
	std::multiset<int> merge;
	//
	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit) {
	  //
	  for(int v = 0; v < merge_group.size(); ++v) {
	    //
	    if(merge_group[v].find(*bit) != merge_group[v].end()) {
	      //
	      merge.insert(v);
	      
	      break;
	    }
	  }
	}

	switch(merge.size()) {
	case 0:
	  // create new group
	  //
	  merge_group.push_back(git->first);
	  
	  break;

	case 1:
	  // add free vertex to the group
	  //
	  for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	    //
	    merge_group[*merge.begin()].insert(*bit);
	  
	  break;

	case 2:
	  // vertices already belong to the same group: do nothing
	  //
	  if(*merge.begin() == *merge.rbegin())
	    //
	    break;

	  // merge two groups
	  //
	  for(std::set<int>::const_iterator mit = merge_group[*merge.begin()].begin(); mit != merge_group[*merge.begin()].end(); ++mit)
	    //
	    merge_group[*merge.rbegin()].insert(*mit);
	  
	  merge_group.erase(merge_group.begin() + *merge.begin());
	  
	  break;
	}
      }
    }
	    
    if(!merge_group.size())
      //
      break;
	    
    std::set<int> pool;
    
    std::vector<std::set<int> > new_red_group(merge_group.size());
    //
    for(int v = 0; v < merge_group.size(); ++v) {
      //
      for(std::set<int>::const_iterator mit = merge_group[v].begin(); mit != merge_group[v].end(); ++mit) {
	//
	pool.insert(*mit);

	for(std::set<int>::const_iterator vit = red_group[*mit].begin(); vit != red_group[*mit].end(); ++vit)
	  //
	  new_red_group[v].insert(*vit);
      }
    }
    
    for(int v = 0; v < red_group.size(); ++v) {
      //
      if(pool.find(v) == pool.end())
	//
	new_red_group.push_back(red_group[v]);
    }

    red_group = new_red_group;
  } // reduction loop

  // zpe graph
  //
  for(int v = 0; v < red_group.size(); ++v) {// zpe factor calculation cycle
    //
    if(red_group[v].size() == 1)
      //
      continue;

    std::map<int, int> vertex_map;
    //
    itemp = 0;
    //
    for(std::set<int>::const_iterator it = red_group[v].begin(); it != red_group[v].end(); ++it, ++itemp)
      //
      vertex_map[*it] = itemp;
	    
    _fg_t new_graph;
    //
    for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
      //
      std::set<int> bond;
      //
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit) {
	//
	if(red_group[v].find(*bit) != red_group[v].end())
	  //
	  bond.insert(vertex_map[*bit]);

	if(bond.size() == 2)
	  //
	  new_graph[bond] = git->second;
      }
    }
    
    ++zpe_graph[new_graph];
  }

  return red_graph;
}

double GraphExpansion::_fg_t::zpe_factor  (const std::vector<double>& freq,
					   double                     temperature,
					   const std::vector<double>& tanh_factor
					   ) const
{
 const char funame [] = "GraphExpansion::_fg_t::zpe_factor: ";

 _check_integrity(freq.size());
 //
 _check_order();

  int    itemp;
  double dtemp;

  if(!size())
    return 1.;

  const int vsize = vertex_size();

  int symm_fac = 0;
  //
  std::set<_fg_t> perm_set = perm_pool(&symm_fac);

  double res = 0.;

  for(std::set<_fg_t>::const_iterator pit = perm_set.begin(); pit != perm_set.end(); ++pit) {
    //
    std::vector<double> freq_band(vsize - 1);
    //
    for(_fg_t::const_iterator git = pit->begin(); git != pit->end(); ++git) {
      //
      dtemp = 0.;
      //
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
	//
	// low frequency
	//
	if(*fit < 0) {
	  //
	  dtemp += 6. * temperature;
	}
	else if(temperature > 0.) {
	  //
	  if(freq[*fit] > 0.)
	    //
	    dtemp += freq[*fit] * tanh_factor[*fit];
	}
	// zero temperature
	//
	else {
	  //
	  dtemp += freq[*fit];
	}
      }
      
      for(int i = *git->first.begin(); i < *git->first.rbegin(); ++i)
	//
	freq_band[i] += dtemp;
    }

    if(temperature > 0.) {
      //
      itemp = 0;
      //
      for(int i = 0; i < freq_band.size(); ++i) {
	//
	dtemp = freq_band[i] / temperature;
	//
	if(dtemp < red_thresh) {
	  //
	  itemp = 1;
	  //
	  break;
	}
      }
      
      if(itemp)
	//
	std::cerr << funame << "WARNING: high frequency cutoff condition not met: " << dtemp << "\n";
    }
    
    dtemp = 1.;
    //
    for(int i = 0; i < freq_band.size(); ++i)
      //
      dtemp /= freq_band[i];
    
    res += dtemp;
  }
  
  res *= (double)symm_fac;

  for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
    //
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      //
      // low frequency
      //
      if(*fit < 0) {
	//
	res /= 12. * temperature;
      }
      else if(temperature > 0.) {
	//
	res /= 2. * freq[*fit] * tanh_factor[*fit];
      }
      else {
	//
	res /= 2. * freq[*fit];
      }
    }
  }
  
  return res;
}

double GraphExpansion::_fg_t::fourier_sum (const std::vector<double>& freq, double temperature) const
{
  const char funame [] = "GraphExpansion::_fg_t::fourier_sum: ";
  
  int    itemp;
  double dtemp;

  if(!size())
    //
    return 1. / temperature;

  _check_integrity(freq.size());
  //
  _check_order();

  const int vsize = vertex_size();

  // frequency map
  //
  std::vector<int> freq_map;
  //
  for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
    //
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      //
      freq_map.push_back(*fit);
    }
  }
  
  std::vector<int> index_shift(freq_map.size());
  //
  std::vector<int> index_range(freq_map.size());
  //
  std::vector<double>     xval(freq_map.size());

  for(int i = 0; i < freq_map.size(); ++i) {
    //
    itemp = freq_map[i];

    if(itemp >= 0) {
      //
      dtemp = freq[itemp] / temperature / 2. / M_PI;
    }
    // low frequency
    //
    else {
      //
      dtemp = 0.;
    }
      
    if(dtemp <= -0.9) {
      throw DeepTunnel() <<funame << "deep tunneling regime: "
			 << "frequency[1/cm] = "	<< freq[itemp] / Phys_const::incm << ",  "
			 << "Temperature[K]  = " << temperature / Phys_const::kelv;
    }

    if(dtemp < 0.) {
      //
      xval[i] = -dtemp * dtemp;
      //
      dtemp = -dtemp;
    }
    else {
      //
      xval[i] = dtemp * dtemp;
    }

    if(dtemp > 1.) {
      //
      itemp = (int)std::ceil(four_cut * dtemp);
    }
    else {
      //
      itemp = (int)std::ceil(four_cut);
    }
    
    index_shift[i] = itemp;
    
    index_range[i] = 2 * itemp + 1;
  }
    
  double res = 0.;

  MultiIndexConvert harmonic_index(index_range);

  for(long ml = 0; ml < harmonic_index.size(); ++ml) {
    //
    std::vector<int> mi = harmonic_index(ml);

    for(int i = 0; i < mi.size(); ++i)
      //
      mi[i] -= index_shift[i];
    
    std::vector<int> constrain(vsize);

    itemp = 0;
    //
    for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
      //
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit, ++itemp) {
	//
	int sign = -1;
	//
	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit, sign += 2) {
	  //
	  constrain[*bit] += sign * mi[itemp];
	}
      }
    }
    
    int test = 0;
    //
    for(std::vector<int>::const_iterator it = constrain.begin(); it != constrain.end(); ++it) {
      //
      if(*it) {
	//
	test = 1;
	//
	break;
      }
    }
    
    if(test)
      //
      continue;

    dtemp = 1.;
    //
    for(int i = 0; i < freq_map.size(); ++i) {
      //
      itemp = mi[i];

      // low frequency
      //
      if(freq_map[i] < 0) {
	//
	if(!itemp) {
	  //
	  test = 1;
	  //
	  break;
	}
	else {
	  //
	  dtemp /= double(itemp * itemp);
	}
      }
      //
      else {
	//
	dtemp /= xval[i] + double(itemp * itemp);
      }
    }

    if(test)
      //
      continue;
    
    res += dtemp;
  }
  
  res /= std::pow(temperature, vsize) * std::pow(4. * M_PI * M_PI * temperature, freq_map.size());

  return res;
}

void GraphExpansion::read_potex (const std::vector<double>& freq, std::istream& from, _potex_t& potex)
{
  const std::string funame [] = "GraphExpansion::read_potex: ";

  int    itemp;
  double dtemp;
  
  if(!freq.size()) {
    ErrOut err_out;
    err_out << funame << "no frequencies";
  }

  std::vector<double> freq_sqrt(freq.size());
  
  for(int i = 0; i < freq.size(); ++i) {
    //
    if(freq[i] <= 0.) {
      ErrOut err_out;
      err_out << funame << "negative frequency[1/cm]: " << freq[i] / Phys_const::incm;
    }
    else {
      //
      freq_sqrt[i] = std::sqrt(freq[i]);
    }
  }
  
  while(from) {
    //
    IO::LineInput lin(from);
    
    IO::String stemp;
    
    std::vector<IO::String> string_value;
    //
    while(lin >> stemp) {
      //
      string_value.push_back(stemp);
    }

    if(!string_value.size())
      //
      continue;
    
    if(string_value[0] == IO::end_key())
      //
      return;

    if(string_value.size() < 4) {
      ErrOut err_out;
      err_out << funame << "not enough values (linear and quadratic terms not allowed)";
    }

    double value = double(string_value.back()) * Phys_const::incm;

    const int imax = string_value.size() - 1;
    
    std::multiset<int> index;
    //
    for(int i = 0; i < imax; ++i) {
      itemp = int(string_value[i]) - 1;

      if(itemp < 0 || itemp >= freq.size()) {
	ErrOut err_out;
	err_out << funame << "index out of range: " << itemp;
      }

      index.insert(itemp);
      value *= freq_sqrt[itemp];
    }

    if(potex.find(index) != potex.end()) {
      ErrOut err_out;
      err_out << funame << "term already in the expansion";
    }
    
    potex[index] = value;
  }

  if(!from) {
    ErrOut err_out;
    err_out << funame << "corrupted";
  }
}

