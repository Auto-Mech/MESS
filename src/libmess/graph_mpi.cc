#include "graph_mpi.hh"

#include "key.hh"
#include "multindex.hh"
#include "permutation.hh"
#include "io.hh"
#include "units.hh"

#include <mpi.h>
#include <cmath>

#ifndef __MPI
#define __MPI
#endif

// reduced frequency tolerance
//
double Graph::Expansion::freq_tol = 0.1;

// low frequency threshold
//
double Graph::Expansion::low_freq_thresh = 0.2;

/********************************************************************************************
 ************************ PERTURBATION THEORY GRAPH EXPANSION *******************************
 ********************************************************************************************/

void Graph::Expansion::correction (const int mode, const double temperature, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "Graph::Expansion::correction: ";
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  //
  // master process
  //
  if(!mpi_rank) {
    //
    IO::Marker funame_marker(funame);
    
    _master(mode, temperature, mmat);
  }
  //
  // work processes
  //
  else {
    //
    _work(mode, temperature, mmat);
  }
}
  
void Graph::Expansion::_work (const int mode, const double temperature, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "Graph::Expansion::_driver: ";

  if(GLOBAL == mode) {
    //
    _global_work(temperature);
  }
  else if(CENTROID == mode) {
    //
    if(mmat.size()) {
      //
      _centroid_work_with_constrain(temperature, mmat);
    }
    else {
      //
      _centroid_work(temperature);
    }
  }
  else {
    //
    ErrOut err_out;

    err_out << funame << "unknown mode: " << mode;
  }
}

/*******************************************************************************************
 ************************************** MASTER PROCESS *************************************
 *******************************************************************************************/

void Graph::Expansion::_master (const int mode, const double temperature, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "Graph::Expansion::_master: ";

  try {
    int    itemp;
    double dtemp;

    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  
    if(mpi_rank) {
      //
      ErrOut err_out;

      err_out << funame << "not a master";
    }

    std::vector<double> tanh_factor;
    //
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    if(temperature > 0.) {
      //
      IO::log << IO::log_offset << "temperature(K) = " << temperature / Phys_const::kelv << "\n\n";

      if(low_freq.size()) {
	//
	IO::log << IO::log_offset << "low frequencies[1/cm]:\n";
	
	for(std::set<int>::const_iterator it = low_freq.begin(); it != low_freq.end(); ++it) {
	  //
	  IO::log << IO::log_offset << std::setw(15) << _red_freq[*it] / Phys_const::incm << "\n";
	}
	
	IO::log << "\n";
      }
    }
    // zero temperature
    else {
      //
      for(int i = 0.; i < _red_freq.size(); ++i)
	//
	if(_red_freq[i] <= 0.) {
	  //
	  ErrOut err_out;

	  err_out << funame << "for zero temperature calculation all frequencies should be real";
	}

      if(GLOBAL == mode) {
	//
	IO::log << IO::log_offset << "zero-point energy correction calculation:\n";
      }
      else {
	//
	IO::log << IO::log_offset << "centroid-constrained low temperature expansion calculation:\n";
      }
    }

    std::map<int, double> graph_value;

    std::map<int, std::map<int, double> > graphex_value; // low temperature expansion

    int work_node = 1;

    MPI::Status stat;

    std::time_t start_time = std::time(0);

    //
    // main loop
    //
    for(int gi = 0; gi < Graph::size(); ++gi) {

      IO::log << IO::log_offset
	      << "graph index = " << std::setw(6) << gi
	      << "   elapsed time[s] = " << std::setw(8) << std::time(0) - start_time
	      << std::endl;

      MultiIndexConvert multi;
      
      if(CENTROID == mode && mmat.size()) {
	//
	multi.resize((Graph::begin() + gi)->size() * 2, _red_freq_index.size());
      }
      else {
	//
	multi.resize((Graph::begin() + gi)->size(),     _red_freq_index.size());
      }

      //
      // normal mode indices cycle
      //
      for(long li = 0; li < multi.size(); ++li) {
	//
	if(work_node < mpi_size) {
	  //
	  MPI::COMM_WORLD.Send(&gi, 1, MPI::INT,  work_node, WORK_TAG);

	  MPI::COMM_WORLD.Send(&li, 1, MPI::LONG, work_node, WORK_TAG);
	  ++work_node;
	}
	else {
	  //
	  int gindex;

	  MPI::COMM_WORLD.Recv(&gindex, 1, MPI::INT,  MPI::ANY_SOURCE, MPI::ANY_TAG, stat);
	  
	  const int tag  = stat.Get_tag();
	  const int from = stat.Get_source();

	  if(tag != WORK_TAG) {
	    //
	    ErrOut err_out;

	    err_out << funame << "wrong tag: " << tag;
	  }
	  
	  if(temperature > 0. || GLOBAL == mode) {
	    //
	    // graph value
	    //
	    MPI::COMM_WORLD.Recv(&dtemp, 1, MPI::DOUBLE, from, tag);

	    graph_value[gindex] += dtemp;
	  }
	  else {
	    //
	    // graph expansion value
	    //
	    _gmap_t gmap(from, tag);

	    for(_gmap_t::const_iterator gmit = gmap.begin(); gmit != gmap.end(); ++gmit) {
	      //
	      graphex_value[gindex][gmit->first] += gmit->second;
	    }
	  }
	  
	  MPI::COMM_WORLD.Send(&gi, 1, MPI::INT,  from, WORK_TAG);
	  MPI::COMM_WORLD.Send(&li, 1, MPI::LONG, from, WORK_TAG);
	}
	//
	//
      } // normal mode indices cycle
      
      //
      //
    } // main loop

    // finish with the rest of receives
    //
    for(int node = 1; node < work_node; ++node) {
      //
      int gindex;

      MPI::COMM_WORLD.Recv(&gindex, 1, MPI::INT,  MPI::ANY_SOURCE, MPI::ANY_TAG, stat);
	  
      const int tag  = stat.Get_tag();
      const int from = stat.Get_source();

      if(tag != WORK_TAG) {
	//
	ErrOut err_out;

	err_out << funame << "wrong tag: " << tag;
      }
      
      // graph value
      //
      if(temperature > 0. || GLOBAL == mode) {
	//
	MPI::COMM_WORLD.Recv(&dtemp, 1, MPI::DOUBLE, from, tag);

	graph_value[gindex] += dtemp;
      }
      // low temperature expansion graph value
      //
      else {
	//
	_gmap_t gmap(from, tag);

	for(_gmap_t::const_iterator gmit = gmap.begin(); gmit != gmap.end(); ++gmit) {
	  //
	  graphex_value[gindex][gmit->first] += gmit->second;
	}
      }
    }

    // break execution loop on working nodes
    //
    for(int node = 1; node < mpi_size; ++node) {
      //
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, node, END_TAG);
    }

    double calc_time = 0.;

    long zpe_count = 0; // zpe calculation #
    long fac_count = 0; // factorization #
    long red_count = 0; // reduction #
    long sum_count = 0; // fourier sum calculation #

    for(int node = 1; node < mpi_size; ++node) {
      //
      long ltemp;

      MPI::COMM_WORLD.Recv(&dtemp, 1, MPI::DOUBLE, node, STAT_TAG);
      //
      calc_time += dtemp;

      MPI::COMM_WORLD.Recv(&ltemp, 1, MPI::LONG,   node, STAT_TAG);
      //
      zpe_count += ltemp;

      if(CENTROID == mode) {
	//
	MPI::COMM_WORLD.Recv(&ltemp, 1, MPI::LONG, node, STAT_TAG);
	//
	fac_count += ltemp;
      }

      if(temperature > 0.) {
	//
	MPI::COMM_WORLD.Recv(&ltemp, 1, MPI::LONG, node, STAT_TAG);
	//
	red_count += ltemp;

	MPI::COMM_WORLD.Recv(&ltemp, 1, MPI::LONG, node, STAT_TAG);
	//
	sum_count += ltemp;
      }
    }

    IO::log << "\n";

    IO::log << IO::log_offset
	    << "Statistics:\n";

    IO::log << IO::log_offset
	    << std::setw(15) << "Calc time[sec]";
    if(zpe_count)
      //
      IO::log << std::setw(15) << "ZPE #";

    if(fac_count)
      //
      IO::log << std::setw(15) << "Factor #";

    if(red_count)
      //
      IO::log << std::setw(15) << "Reduct #";

    if(sum_count)
      //
      IO::log << std::setw(15) << "FourSum #";

    IO::log << "\n";

    IO::log << IO::log_offset
	    << std::setw(15) << calc_time;

    if(zpe_count)
      //
      IO::log << std::setw(15) << zpe_count;

    if(fac_count)
      //
      IO::log << std::setw(15) << fac_count;

    if(red_count)
      //
      IO::log << std::setw(15) << red_count;

    if(sum_count)
      //
      IO::log << std::setw(15) << sum_count;

    IO::log << "\n\n";

    if(temperature > 0.) {
      //
      std::map<int, double> corr;
      //
      for(int i = 0; i < Graph::size(); ++i) {
	//
	corr[(Graph::begin() + i)->size()] += graph_value[i];
      }

      if(GLOBAL == mode) {
	//
	IO::log << IO::log_offset << "anharmonic correction:\n";
      }
      else {
	//
	IO::log << IO::log_offset << "centroid-constrained anharmonic correction:\n";
      }

      IO::log << IO::log_offset 
	      << std::setw(2) << "BO" 
	      << std::setw(15) << "Value" 
	      << "\n";

      double res = 0.;
      for(std::map<int, double>::const_iterator cit = corr.begin(); cit != corr.end(); ++cit) {
	//
	res += cit->second;

	IO::log << IO::log_offset 
		<< std::setw(2) << cit->first 
		<< std::setw(15) << res 
		<< "\n";
      }
      IO::log << "\n";
    }
    // zero temperature
    //
    else {
      //
      // global mode (no constrains)
      //
      if(GLOBAL == mode) {
	//
	std::map<int, double> zpe;
	//
	for(int i = 0; i < Graph::size(); ++i) {
	  //
	  zpe[(Graph::begin() + i)->size()] += graph_value[i];
	}

	IO::log << IO::log_offset << "zero-point energy correction, 1/cm:\n";
	//
	IO::log << IO::log_offset 
		<< std::setw(2) << "BO" 
		<< std::setw(15) << "Value" 
		<< "\n";

	double res = 0.;
	for(std::map<int, double>::const_iterator gzit = zpe.begin(); gzit != zpe.end(); ++gzit) {
	  //
	  res += gzit->second;

	  IO::log << IO::log_offset 
		  << std::setw(2) << gzit->first 
		  << std::setw(15) << -res / Phys_const::incm
		  << "\n";
	}
	IO::log << "\n";
      }
      //
      // centroid mode
      //
      else {
	//
	int print_term_max = 3;

	std::map<int, std::map<int, double> > graphex;

	for(int i = 0; i < Graph::size(); ++i) {
	  //
	  for(std::map<int, double>::const_iterator gxit = graphex_value[i].begin(); gxit != graphex_value[i].end(); ++gxit) {
	    //
	    graphex[(Graph::begin() + i)->size()][gxit->first] += gxit->second;
	  }
	}
	
	IO::log << IO::log_offset << "centroid-constrained low temperature expansion:\n";
	//
	IO::log << IO::log_offset 
		<< std::setw(2)  << "BO"
		<< std::setw(15) << "1/T(K) term"
		<< std::setw(15) << "T^0 term";

	for(int i = 1; i <= print_term_max; ++i) {
	  //
	  IO::log << std::setw(9) << "T(K)^" << i << " term";
	}
	IO::log << "\n";

	std::map<int, double> res;
	for(std::map<int, std::map<int, double> >::const_iterator gxit = graphex.begin(); gxit != graphex.end(); ++gxit) {
	  //
	  for(std::map<int, double>::const_iterator it = gxit->second.begin(); it != gxit->second.end(); ++it) {
	    //
	    res[it->first] += it->second;
	  }

	  itemp = -res.begin()->first;
	  if(itemp > 1) {
	    ErrOut err_out;
	    err_out << funame << "low temperature expansion begins with 1/T^" << itemp << "terms";
	  }

	  IO::log << IO::log_offset << std::setw(2) << gxit->first;
	
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
    }
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

//
// non-constrain graph perturbation theory
//
void Graph::Expansion::_global_work (const double temperature) const
{
  const char funame [] = "Graph::Expansion::_global_work: ";

  try {
    //
    int    itemp;
    double dtemp;
    bool   btemp;

    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    MPI::Status stat;

    std::clock_t calc_time = std::clock();

    long fac_count = 0;
    long zpe_count = 0;
    long red_count = 0;
    long sum_count = 0;

    //
    // main loop
    //
    while(1) {
      // graph index
      int gindex;
      MPI::COMM_WORLD.Recv(&gindex, 1, MPI::INT, MASTER, MPI::ANY_TAG, stat);

      const int tag = stat.Get_tag();
 
      if(END_TAG == tag) {
	//
	break;
      }
      else if(WORK_TAG != tag) {
	ErrOut err_out;
	err_out << funame << "wrong tag: " << tag;
      }

      Graph::const_iterator graphit = Graph::begin() + gindex;

      // normal mode indices converted to a single linear index
      //
      long li;
      MPI::COMM_WORLD.Recv(&li, 1, MPI::LONG, MASTER, WORK_TAG);
      
      std::vector<int> corrin = MultiIndexConvert(graphit->size(), _red_freq_index.size())(li);
      //
      const std::vector<std::multiset<int> > vertex_map = graphit->vertex_bond_map();
      //
      const int vertex_size = vertex_map.size();
      //
      double gfactor = 1.;

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
	dtemp = 0.;
	MPI::COMM_WORLD.Send(&gindex, 1, MPI::INT,    MASTER, WORK_TAG);
	MPI::COMM_WORLD.Send(&dtemp,  1, MPI::DOUBLE, MASTER, WORK_TAG);

	continue;
      }

      gfactor /= (double)graphit->symmetry_factor();

      if(vertex_size % 2)
	//
	gfactor = -gfactor;

      // frequency adapted graph
      //
      FreqGraph mod_graph;
      
      itemp = 0;
      for(GenGraph::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {
	//
	int ci = corrin[itemp]; // correlator index
	//
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
	  // positive temperature correlator value
	  //
	  if(temperature > 0.) {
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
	// insert frequency index into the frequency adapted graph
	//
	else {
	  //
	  mod_graph[bond].insert(fi);
	}
      }

      // zero temperature integral (zpe factor) calculation
      //
      if(temperature <= 0.) {
	//
	if(mod_graph.size()) {
	  //
	  ++zpe_count;
	
	  gfactor *= mod_graph.zpe_factor(_red_freq);
	}
      }
      // thermal whole integral evaluation
      //
      else {
	//
	itemp = vertex_size - mod_graph.vertex_size();
	//
	if(itemp)
	  //
	  gfactor /= std::pow(temperature, (double)itemp);
	
	if(mod_graph.size()) {
	  //
	  ++fac_count;

	  // factorize frequency adapted graph into several connected graphs
	  //
	  _mg_t fac_graph = mod_graph.factorize();

	  for(_mg_t::const_iterator fgit = fac_graph.begin(); fgit != fac_graph.end(); ++fgit) {
	    //
	    ++red_count;

	    // graph reduction
	    //
	    _mg_t zpe_graph;
	    //
	    FreqGraph red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);
	  
	    double gf = 1.;

	    if(!red_graph.size()) {
	      //
	      gf /= temperature;
	    }
	    // reduced graph fourier sum calculation
	    //
	    else {
	      //
	      ++sum_count;

	      gf *= red_graph.fourier_sum(_red_freq, temperature);
	    }

	    // zpe graph cycle
	    //
	    for(_mg_t::const_iterator zgit = zpe_graph.begin(); zgit != zpe_graph.end(); ++zgit) {
	      //
	      // low temperature integral (zpe factor) calculation
	      //
	      ++zpe_count;

	      dtemp = zgit->first.zpe_factor(_red_freq, temperature, tanh_factor);
	      
	      for(int i = 0; i < zgit->second; ++i)
		//
		gf *= dtemp;
	    }
	    
	    for(int i = 0; i < fgit->second; ++i)
	      //
	      gfactor *= gf;
	  }
	}
      }
      
      MPI::COMM_WORLD.Send(&gindex,  1, MPI::INT,    MASTER, WORK_TAG);
      MPI::COMM_WORLD.Send(&gfactor, 1, MPI::DOUBLE, MASTER, WORK_TAG);
      //
      //
    } // main loop

    calc_time = std::clock() - calc_time;

    dtemp = calc_time / CLOCKS_PER_SEC;
    
    MPI::COMM_WORLD.Send(&dtemp, 1, MPI::DOUBLE, MASTER, STAT_TAG);

    MPI::COMM_WORLD.Send(&zpe_count, 1, MPI::LONG,   MASTER, STAT_TAG);

    if(temperature > 0.) {
      MPI::COMM_WORLD.Send(&red_count, 1, MPI::LONG, MASTER, STAT_TAG);
      MPI::COMM_WORLD.Send(&sum_count, 1, MPI::LONG, MASTER, STAT_TAG);
    }
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void Graph::Expansion::_centroid_work (const double temperature) const
{
  const char funame [] = "Graph::Expansion::_centroid_work: ";

  try {
    int    itemp;
    double dtemp;
    bool   btemp;

    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    MPI::Status stat;

    std::clock_t calc_time = std::clock();

    long zpe_count = 0;
    long fac_count = 0;
    long red_count = 0;
    long sum_count = 0;

    //
    // main loop
    //
    while(1) {
      //
      int gindex;
      MPI::COMM_WORLD.Recv(&gindex, 1, MPI::INT, MASTER, MPI::ANY_TAG, stat);

      const int tag = stat.Get_tag();

      if(END_TAG == tag) {
	//
	break;
      }
      else if(WORK_TAG != tag) {
	ErrOut err_out;
	err_out << funame << "wrong tag: " << tag;
      }

      long li;
      MPI::COMM_WORLD.Recv(&li, 1, MPI::LONG, MASTER, WORK_TAG);
      //
      Graph::const_iterator graphit = Graph::begin() + gindex;
      //
      std::vector<int> corrin = MultiIndexConvert(graphit->size(), _red_freq_index.size())(li);
      //
      const std::vector<std::multiset<int> > vertex_map = graphit->vertex_bond_map();
      //
      const int vertex_size = vertex_map.size();
      //
      double gvalue = 0.;
      //
      _gmap_t gmap;
      //
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

	MPI::COMM_WORLD.Send(&gindex, 1, MPI::INT, MASTER, WORK_TAG);
	
	if(temperature > 0.) {
	  //
	  MPI::COMM_WORLD.Send(&gvalue, 1, MPI::DOUBLE, MASTER, WORK_TAG);
	}
	else {
	  //
	  gmap.send(MASTER, WORK_TAG);
	}

	continue;
      }

      potfac /= (double)graphit->symmetry_factor();
      //
      if(vertex_size % 2)
	potfac = -potfac;

      //
      // centroid correction mask cycle
      //
      for(MultiIndex cmask(graphit->size(), 1); !cmask.end(); ++cmask) {
	//
	double gfactor = potfac;

	int    t_count = 0;
	
	FreqGraph mod_graph;

	itemp = 0;

	btemp = false;
	//
	for(GenGraph::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {
	  //
	  int ci = corrin[itemp]; // normal mode index

	  int fi = _red_freq_index[ci];  // reduced frequency index

	  if(cmask[itemp]) {
	    //
	    std::set<int> bond;

	    for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
	      //
	      bond.insert(*it);
	
	    // bond loop
	    //
	    if(bond.size() == 1) {
	      //
	      // low frequency correlator incorporates centroid correction
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		gfactor /= 12. * temperature;
	      }
	      // thermal correlator value
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
	    // include frequency index into the graph
	    //
	    else {
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
	  // centroid correction for low frequency correlator incorporated into it
	  //
	  else if(low_freq.find(fi) != low_freq.end()) {
	    //
	    btemp = true;

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
	  else
	    t_count -= itemp;
	}

	// frequency adapted graph calculation
	//
	if(mod_graph.size()) {
	  //
	  ++fac_count;

	  // graph factorization into the set of connected graphs
	  //
	  _mg_t fac_graph = mod_graph.factorize();

	  // factorized graph cycle
	  //
	  for(_mg_t::const_iterator fgit = fac_graph.begin(); fgit != fac_graph.end(); ++fgit) {
	    //
	    // zero temperature integral (zpe factor) calculation
	    //
	    if(temperature <= 0.) {
	      //
	      ++zpe_count;

	      dtemp = fgit->first.zpe_factor(_red_freq);

	      for(int i = 0; i < fgit->second; ++i)
		//
		gfactor *= dtemp;
	    
	      t_count -= fgit->second;
	    }
	    // thermal whole integral calculation
	    //
	    else {
	      //
	      ++red_count;
	    
	      _mg_t zpe_graph;
	      //
	      FreqGraph red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);
	  
	      double gf = 1.;

	      if(!red_graph.size()) {
		//
		gf /= temperature;
	      }
	      else {
		//
		++sum_count;

		gf *= red_graph.fourier_sum(_red_freq, temperature);
	      }

	      // zpe graph cycle
	      //
	      for(_mg_t::const_iterator zgit = zpe_graph.begin(); zgit != zpe_graph.end(); ++zgit) {
		//
		++zpe_count;

		dtemp = zgit->first.zpe_factor(_red_freq, temperature, tanh_factor);

		for(int i = 0; i < zgit->second; ++i)
		  //
		  gf *= dtemp;
		//
		//
	      } // zpe graph cycle

	      for(int i = 0; i < fgit->second; ++i)
		//
		gfactor *= gf;
	      //
	      //
	    } // whole integral calculation
	    //
	    //
	  } // factorized graph cycle
	  //
	  //
	} // frequency adapted graph evaluation
	  
	if(temperature > 0.) {
	  //
	  gvalue += gfactor;
	}
	else {
	  //
	  gmap[t_count] += gfactor;
	}
	//
	//
      } // centroid correction mask cycle
	
      MPI::COMM_WORLD.Send(&gindex, 1, MPI::INT, MASTER, WORK_TAG);

      if(temperature > 0.) {
	//
	MPI::COMM_WORLD.Send(&gvalue, 1, MPI::DOUBLE, MASTER, WORK_TAG);
      }
      else {
	//
	gmap.send(MASTER, WORK_TAG);
      }
      //
      //
    } // main loop

    calc_time = std::clock() - calc_time;

    dtemp = calc_time / CLOCKS_PER_SEC;
    
    MPI::COMM_WORLD.Send(&dtemp, 1, MPI::DOUBLE, MASTER, STAT_TAG);

    MPI::COMM_WORLD.Send(&zpe_count, 1, MPI::LONG,   MASTER, STAT_TAG);

    MPI::COMM_WORLD.Send(&fac_count, 1, MPI::LONG,   MASTER, STAT_TAG);

    if(temperature > 0.) {
      MPI::COMM_WORLD.Send(&red_count, 1, MPI::LONG, MASTER, STAT_TAG);

      MPI::COMM_WORLD.Send(&sum_count, 1, MPI::LONG, MASTER, STAT_TAG);
    }
  }
  catch(Error::General) {
    //
    std::cerr << funame << "Oops\n";

    throw;
  }
}

void Graph::Expansion::_centroid_work_with_constrain (const double temperature, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "Graph::Expansion::_centroid_work_with_constrain: ";

  try {
    int    itemp;
    double dtemp;
    bool   btemp;

    // mpi stuff
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    MPI::Status stat;

    std::clock_t calc_time = std::clock();

    long zpe_count = 0;
    long fac_count = 0;
    long red_count = 0;
    long sum_count = 0;

    //
    // main loop
    //
    while(1) {
      int gindex;
      MPI::COMM_WORLD.Recv(&gindex, 1, MPI::INT, MASTER, MPI::ANY_TAG, stat);

      const int tag = stat.Get_tag();

      if(END_TAG == tag) {
	//
	break;
      }
      else if(WORK_TAG != tag) {
	ErrOut err_out;
	err_out << funame << "wrong tag: " << tag;
      }

      Graph::const_iterator graphit = Graph::begin() + gindex;

      long li;
      MPI::COMM_WORLD.Recv(&li, 1, MPI::LONG, MASTER, WORK_TAG);
      //
      std::vector<int> corrin = MultiIndexConvert(graphit->size() * 2, _red_freq_index.size())(li);
      //
      const int vertex_size = graphit->vertex_size();
      //
      std::vector<std::set<int> > vertex_map(vertex_size);
      //
      itemp = 0;
      for(GenGraph::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {
	//
	if(git->size() != 2) {
	  ErrOut err_out;
	  err_out << funame << "bond should connect exactly two vertices: " << git->size();
	}
      
	int v = 0;
	for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit, ++v) {
	  //
	  if(!vertex_map[*bit].insert(2 * itemp + v).second) {
	    ErrOut err_out;
	    err_out << funame << "duplicated frequency index: " << 2 * itemp + v;
	  }
	}
      }
      //
      double gvalue = 0.;
      //
      _gmap_t gmap;
      //
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
	MPI::COMM_WORLD.Send(&gindex, 1, MPI::INT, MASTER, WORK_TAG);
	
	if(temperature > 0.) {
	  //
	  MPI::COMM_WORLD.Send(&gvalue, 1, MPI::DOUBLE, MASTER, WORK_TAG);
	}
	else {
	  //
	  gmap.send(MASTER, WORK_TAG);
	}
	
	continue;
      }

      potfac /= (double)graphit->symmetry_factor();
      //
      if(vertex_size % 2)
	//
	potfac = -potfac;

      // centroid correction mask cycle
      //
      for(MultiIndex cmask(graphit->size(), 1); !cmask.end(); ++cmask) {
	//
	double gfactor = potfac;
	
	int    t_count = 0;

	// frequency adapted graph
	//
	FreqGraph mod_graph;

	itemp = 0;

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
	  if(cmask[itemp]) {
	    //
	    // cross term contribution comes from centroid correction only
	    //
	    if(ci[0] != ci[1]) {
	      //
	      btemp = true;

	      break;
	    }

	    const int fi = _red_freq_index[ci[0]];
	    
	    std::set<int> bond;
	    //
	    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
	      //
	      bond.insert(*bit);
	
	    // bond loop
	    //
	    if(bond.size() == 1) {
	      //
	      // second order expansion term in low frequency correlator; zero order constant term included into centroid correction
	      //
	      if(low_freq.find(fi) != low_freq.end()) {
		//
		gfactor /= 12. * temperature;
	      }
	      // thermal correlator value
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
	    // add frequency index to frequency adapted graph
	    //
	    else {
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

	    // centroid cross-term contribution
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
	    // diagonal centroid correction
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

	      // low frequency zero order expansion term included into centroid correction
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

		  err_out << funame << "low frequency " << fi << " does not have corresponding value in M matrix";
		}
	      }
	      // centroid correction value
	      //
	      else if(mmat.end() != mit) {
		//
		gfactor *= - mit->second / dtemp / dtemp;
	      }
	      else {
		//
		btemp = true;

		break;
	      }

	      if(temperature > 0.) {
		//
		gfactor *= temperature;
	      }
	      else
		//
		++t_count;
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
	  else
	    //
	    t_count -= itemp;
	}

	// frequency adapted graph evaluation
	//
	if(mod_graph.size()) {
	  //
	  ++fac_count;

	  // graph factorization into connected ones
	  //
	  _mg_t fac_graph = mod_graph.factorize();

	  // factorized graph cycle
	  //
	  for(_mg_t::const_iterator fgit = fac_graph.begin(); fgit != fac_graph.end(); ++fgit) {
	    //
	    // zero temperature integral (zpe factor) calculation
	    //
	    if(temperature <= 0.) {
	      //
	      ++zpe_count;

	      dtemp = fgit->first.zpe_factor(_red_freq);

	      for(int i = 0; i < fgit->second; ++i)
		//
		gfactor *= dtemp;
	    
	      t_count -= fgit->second;
	    }
	    // thermal whole integral calculation
	    //
	    else {
	      //
	      ++red_count;

	      // graph reduction
	      //
	      _mg_t zpe_graph;
	      //
	      FreqGraph red_graph = fgit->first.reduce(_red_freq, temperature, tanh_factor, zpe_graph);
	  
	      double gf = 1.;

	      if(!red_graph.size()) {
		//
		gf /= temperature;
	      }
	      // reduced graph fourier sum calculation
	      //
	      else {
		//
		++sum_count;

		gf *= red_graph.fourier_sum(_red_freq, temperature);
	      }

	      // zpe graph cycle
	      //
	      for(_mg_t::const_iterator zgit = zpe_graph.begin(); zgit != zpe_graph.end(); ++zgit) {
		//
		// low temperature / high frequency integral (zpe factor) calculation
		//
		++zpe_count;

		dtemp = zgit->first.zpe_factor(_red_freq, temperature, tanh_factor);

		for(int i = 0; i < zgit->second; ++i)
		  //
		  gf *= dtemp;

		//
		//
	      } // zpe graph cycle

	      for(int i = 0; i < fgit->second; ++i)
		//
		gfactor *= gf;
	      //
	      //
	    } // thermal whole integral calculation
	    //
	    //
	  } // factorized graph cycle
	  //
	  //
	} // frequency adapted graph evaluation
	  
	if(temperature > 0.) {
	  //
	  gvalue += gfactor;
	}
	else {
	  //
	  gmap[t_count] += gfactor;
	}
	//
	//
      } // centroid correction mask cycle
	
      MPI::COMM_WORLD.Send(&gindex, 1, MPI::INT, MASTER, WORK_TAG);
	
      if(temperature > 0.) {
	//
	MPI::COMM_WORLD.Send(&gvalue, 1, MPI::DOUBLE, MASTER, WORK_TAG);
      }
      else {
	//
	gmap.send(MASTER, WORK_TAG);
      }
      //
      //
    } // main loop

    calc_time = std::clock() - calc_time;

    dtemp = calc_time / CLOCKS_PER_SEC;
    
    MPI::COMM_WORLD.Send(&dtemp, 1, MPI::DOUBLE, MASTER, STAT_TAG);

    MPI::COMM_WORLD.Send(&zpe_count, 1, MPI::LONG,   MASTER, STAT_TAG);

    MPI::COMM_WORLD.Send(&fac_count, 1, MPI::LONG,   MASTER, STAT_TAG);

    if(temperature > 0.) {
      MPI::COMM_WORLD.Send(&red_count, 1, MPI::LONG, MASTER, STAT_TAG);

      MPI::COMM_WORLD.Send(&sum_count, 1, MPI::LONG, MASTER, STAT_TAG);
    }
  }
  catch(Error::General) {
    //
    std::cerr << funame << "Oops\n";

    throw;
  }
}

// perturbation graph theory initializer
//
void Graph::Expansion::init (const std::vector<double>& freq, const potex_t& pex)
{
  const char funame [] = "Graph::Expansion::init: ";

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

}

#include "graph_include.cc"

/*******************************************************************************************
 ******************************** COMMUNICAIION CLASS BASE *********************************
 *******************************************************************************************/

void Graph::Expansion::_gbase_t::send (int node, int tag) const
{
  Array<char> buff;
  
  int pack_size = _pack(buff);

  MPI::COMM_WORLD.Send(&pack_size, 1,         MPI::INT,    node, tag);

  MPI::COMM_WORLD.Send(buff,       pack_size, MPI::PACKED, node, tag);
}

void Graph::Expansion::_gbase_t::recv (int node, int tag)
{
  int pack_size;

  MPI::COMM_WORLD.Recv(&pack_size, 1,   MPI::INT,    node, tag);
  
  Array<char> buff(pack_size);
  //
  MPI::COMM_WORLD.Recv(buff, pack_size, MPI::PACKED, node, tag);

  _unpack(buff);
}


/*******************************************************************************************
 ************************* COMMUNICAIION CLASS FOR MAP<INT, DOUBLE> ************************
 *******************************************************************************************/

int Graph::Expansion::_gmap_t::_pack (Array<char>& buff) const
{
  int itemp;

  itemp = (size() + 1) * MPI::INT.Pack_size(1, MPI::COMM_WORLD) + size() * MPI::DOUBLE.Pack_size(1, MPI::COMM_WORLD);

  buff.resize(itemp);
  
  int pos = 0;

  itemp = size();

  MPI::INT.Pack(&itemp, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  for(std::map<int, double>::const_iterator mit = this->begin(); mit != this->end(); ++mit) {
    //
    MPI::INT.Pack(&mit->first, 1, buff, buff.size(), pos, MPI::COMM_WORLD);
    
    MPI::DOUBLE.Pack(&mit->second, 1, buff, buff.size(), pos, MPI::COMM_WORLD);
  }
  
  return pos;
}

void Graph::Expansion::_gmap_t::_unpack (const Array<char>& buff)
{
  const char funame [] = "Graph::Expansion::_gmap_t::_unpack: ";
  
  int    itemp;
  double dtemp;

  int pos = 0;

  int map_size;
  //
  MPI::INT.Unpack(buff, buff.size(), &map_size, 1, pos, MPI::COMM_WORLD);
  
  for(int i = 0; i < map_size; ++i) {
    //
    MPI::INT.Unpack(buff, buff.size(), &itemp, 1, pos, MPI::COMM_WORLD);
    
    MPI::DOUBLE.Unpack(buff, buff.size(), &dtemp, 1, pos, MPI::COMM_WORLD);

    if(!this->insert(std::make_pair(itemp, dtemp)).second) {
      //
      ErrOut err_out;
      
      err_out << funame << "cannot insert int-double pair: " << itemp << ", " << dtemp << ": map values:\n";
      
      for(std::map<int, double>::const_iterator mit = this->begin(); mit != this->end(); ++mit)
	//
	err_out << mit->first << "  " << mit->second << "\n";
      
      err_out << "\n";
    }
  }
}

