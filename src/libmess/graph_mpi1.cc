#include "graph_mpi1.hh"
#include "key.hh"

#include <cstdarg>

int    GraphExpansion::WORK_NODE       = GraphExpansion::SUM_SERV + 2;// first warking node (default to one driver)

int    GraphExpansion::bond_max        = 6;      // maximal number of graph bonds (edges)
double GraphExpansion::freq_tol        = 0.1;    // reduced frequency tolerance
double GraphExpansion::four_cut        = 5.;     // fourier sum graph evaluation cutoff
double GraphExpansion::red_thresh      = 6.;     // low temperature / high frequency graph reduction threshold
double GraphExpansion::low_freq_thresh = 1.e-2;  // low frequency threshold (relative to the temperature)

GraphExpansion::_Convert  GraphExpansion::_gbase_t::_convert(20, 100); // 20 indices and 100 frequencies are allowed

/********************************************************************************************
 ************************ PERTURBATION THEORY GRAPH EXPANSION *******************************
 ********************************************************************************************/

void GraphExpansion::correction (int mode, double temperature, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "GraphExpansion::correction: ";
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  // master process
  if(!mpi_rank) {
    IO::Marker funame_marker(funame);
    
    _master(temperature, mode);
  }
  //
  // integral server
  else if(mpi_rank == INT_SERV) {
    _int_server(temperature);
  }
  // zpe database server
  else if(mpi_rank == ZPE_SERV) {
    _zpe_server(temperature);
  }
  //
  // fourier sum database server
  else if(mpi_rank == SUM_SERV) {
    _sum_server(temperature);
  }
  //
  // working nodes
  else if(_is_work_node(mpi_rank)) {
    _node_work(temperature);
  }
  else if(_is_driver(mpi_rank)) {
    _driver(temperature, mode, mmat);
  }
}
  
void GraphExpansion::_driver (double temperature, int mode, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "GraphExpansion::_driver: ";

  if(GLOBAL == mode) {
    _global_driver(temperature);
  }
  else if(CENTROID == mode) {
    if(mmat.size())
      _centroid_driver_with_constrain(temperature, mmat);
    else
      _centroid_driver(temperature);
  }
  else {
    ErrOut err_out;
    err_out << funame << "unknown mode: " << mode;
  }
}

void GraphExpansion::_global_driver (double temperature) const
{
  const char funame [] = "GraphExpansion::_global_driver: ";

  try {
    int    itemp;
    double dtemp;
    bool   btemp;

    // mpi stuff
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    // posted sends
    typedef std::list<SharedPointer<_gbase_t> >  sent_t;
    sent_t sending;

    //posted receives
    std::list<_gin_t> receiving;

    // symbolic map
    std::map<int, std::map<_mg_t, double> > sym_map;

    if(!_is_driver(mpi_rank)) {
      ErrOut err_out;
      err_out << funame << "not a driver";
    }
    
    // calculation stuff
    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    while(1) {
      // mpi stuff
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, MASTER, GINDEX_REQUEST);

      int gindex;
      gindex = _gnode_t(MASTER, GINDEX_TAG);

      //std::cout << "gindex = " << gindex << std::endl;

      if(_graph_data.size() == gindex)
	break;
    
      std::vector<_graph_t>::const_iterator graphit = _sorted_graph.begin() + gindex;

      // calculation stuff
      //
      const std::vector<std::multiset<int> > vertex_map = graphit->vertex_bond_map();
      //
      const int vertex_size = vertex_map.size();
     
      MultiIndexConvert corr_multi_index(graphit->size(), _red_freq_index.size());

      for(long corr_lin = 0; corr_lin < corr_multi_index.size(); ++corr_lin) {

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

	gfactor /= (double)graphit->symmetry_factor();
	if(vertex_size % 2)
	  gfactor = -gfactor;

	// frequency adapted graph setting
	_fg_t mod_graph;
      
	itemp = 0;
	for(_graph_t::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {

	  int ci = corrin[itemp]; // correlator index
	
	  int fi = _red_freq_index[ci]; // reduced frequency index

	  if(low_freq.find(fi) != low_freq.end()) {
	    if(_red_freq[fi] > 0.)
	      gfactor *=  temperature / _red_freq[fi] / _red_freq[fi];
	    else
	      gfactor *= -temperature / _red_freq[fi] / _red_freq[fi];

	    continue;
	  }

	  std::set<int> bond;
	  for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
	    bond.insert(*it);

	  if(bond.size() == 1) {
	    // bond loop
	    if(temperature > 0.)
	      gfactor /= 2. * _red_freq[fi] * tanh_factor[fi];
	    else
	      gfactor /= 2. * _red_freq[fi];

	    continue;
	  }

	  mod_graph[bond].insert(fi);
	}

	itemp = vertex_size - mod_graph.vertex_size();
	if(itemp && temperature > 0.) {
	  gfactor /= std::pow(temperature, (double)itemp);
	}
	
	// graph factorization into connected graphs
	_mg_t sym = mod_graph.factorize();
	
	// send graphs (symbol components) to integral server
	for(_mg_t::const_iterator sit = sym.begin(); sit != sym.end(); ++sit) {// modified graph factor cycle
	  sending.push_back(SharedPointer<_gbase_t>(new _gself_t(sit->first)));
	  sending.back()->isend(INT_SERV, INT_TAG);
	  
	  receiving.push_back(_gin_t());
	  receiving.back().irecv(INT_SERV, INT_TAG);
	}
	  
	sym_map[gindex][sym] += gfactor;
      }// correlator indices cycle

      // clean completed sends
      for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	if((*sit)->Test())
	  sit = sending.erase(sit);
	else
	  ++sit;
      
      // extract and clean completed receives
      std::map<_fg_t, double> received;
      for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
	if(pit->Test()) {
	  received[pit->first] = pit->second;
          
	  pit = receiving.erase(pit);
	}
	else
	  ++pit;
      }

      // update symbolic map
      _update_symbolic_map(sym_map, received, &sending);

    } // graph cycle

    // finish with receives
    std::map<_fg_t, double> received;
    for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
      pit->Wait();

      received[pit->first] = pit->second;
      
      pit = receiving.erase(pit);
    }
    
    // update symbolic map
    _update_symbolic_map(sym_map, received, &sending);

    // finish with sends
    for(sent_t::iterator sit = sending.begin(); sit != sending.end();) {
      (*sit)->Wait();
      sit = sending.erase(sit);
    }

    if(sym_map.size()) {
      ErrOut err_out;
      err_out << funame << "some unresolved graphs are still in symbolic map";
    }
    MPI::COMM_WORLD.Send(0, 0, MPI::INT, MASTER, END_TAG);

    MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MASTER, END_TAG);
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void GraphExpansion::_centroid_driver (double temperature) const
{
  const char funame [] = "GraphExpansion::centroid_correction: ";

  try {
    int    itemp;
    double dtemp;
    bool   btemp;

    // mpi stuff
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    if(!_is_driver(mpi_rank)) {
      ErrOut err_out;
      err_out << funame << "not a driver";
    }

    typedef std::list<SharedPointer<_gbase_t> >  sent_t;

    sent_t            sending;
    std::list<_gin_t> receiving;
    
    std::map<int, std::map<_mg_t, double> >                 sym_val_map;
    std::map<int, std::map<int, std::map<_mg_t, double> > > sym_zpe_map;

    // calculation stuff
    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    // graph cycle
    while(1) {
      // graph index request
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, MASTER, GINDEX_REQUEST);

      // graph index
      const int gindex = _gnode_t(MASTER, GINDEX_TAG);

      if(_graph_data.size() == gindex)
	break;
    
      std::vector<_graph_t>::const_iterator graphit = _sorted_graph.begin() + gindex;

      // calculation stuff
      const std::vector<std::multiset<int> > vertex_map = graphit->vertex_bond_map();
      const int vertex_size = vertex_map.size();

      MultiIndexConvert corr_multi_index(graphit->size(), _red_freq_index.size());

      //
      // normal mode indices cycle
      //
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
	  //
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

	if(btemp)
	  continue;


	potfac /= (double)graphit->symmetry_factor();
	//
	if(vertex_size % 2)
	  //
	  potfac = -potfac;

	//
	// classical correlator mask cycle
	//
	for(MultiIndex cmask(graphit->size(), 1); !cmask.end(); ++cmask) {
	
	  double gfactor = potfac;
	  int    t_count = 0;
	  
	  _fg_t mod_graph;

	  itemp = 0;
	  btemp = false;
	  for(_graph_t::const_iterator mit = graphit->begin(); mit != graphit->end(); ++mit, ++itemp) {

	    int ci = corrin[itemp]; // normal mode index
	    int fi = _red_freq_index[ci];  // reduced frequency index

	    if(cmask[itemp]) {
	      std::set<int> bond;
	      for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
		bond.insert(*it);
	
	      // bond loop
	      if(bond.size() == 1) {
		// low frequency correlator value
		if(low_freq.find(fi) != low_freq.end()) {
		  gfactor /= 12. * temperature;
		}
		// positive temperature correlator value
		else if(temperature > 0.) {
		  gfactor /= 2. * _red_freq[fi] * tanh_factor[fi];
		}
		// zero temperature correlator value
		else {
		  gfactor /= 2. * _red_freq[fi];
		}
	      }
	      else {
		// insert frequency index into the graph
		mod_graph[bond].insert(fi);
	      }
	    }
	    // low frequency centroid correction already included into correlator
	    else if(low_freq.find(fi) != low_freq.end()) {
	      btemp = true;
	      break;
	    }
	    // centroid correction
	    else {
	      if(temperature > 0.)
		gfactor *= temperature;
	      else
		++t_count;

	      if(_red_freq[fi] > 0.)
		gfactor /= -_red_freq[fi] * _red_freq[fi];
	      else
		gfactor /=  _red_freq[fi] * _red_freq[fi];
	    }
	  }

	  if(btemp)
	    continue;
	
	  itemp = vertex_size - mod_graph.vertex_size();
	  if(itemp) {
	    if(temperature > 0.)
	      gfactor /= std::pow(temperature, (double)itemp);
	    else
	      t_count -= itemp;
	  }

	  // graph factorization into connected graphs
	  _mg_t sym = mod_graph.factorize();
	  
	  if(temperature <=  0.)
	    for(_mg_t::const_iterator mit = sym.begin(); mit != sym.end(); ++mit)
	      t_count -= mit->second;

	  // send graphs (symbol components) to integral server
	  for(_mg_t::const_iterator sit = sym.begin(); sit != sym.end(); ++sit) {// modified graph factor cycle
	    sending.push_back(SharedPointer<_gbase_t>(new _gself_t(sit->first)));
	    sending.back()->isend(INT_SERV, INT_TAG);
	  
	    receiving.push_back(_gin_t());
	    receiving.back().irecv(INT_SERV, INT_TAG);
	  }
	  
	  if(temperature > 0.)
	    sym_val_map[gindex][sym] += gfactor;
	  else
	    sym_zpe_map[gindex][t_count][sym] += gfactor;
	  //
	  //
	} // correlator mask cycle
	
	// clean completed sends
	for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	  if((*sit)->Test())
	    sit = sending.erase(sit);
	  else
	    ++sit;
	  
      // extract and clean completed receives
	std::map<_fg_t, double> received;
	for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
	  if(pit->Test()) {
	    received[pit->first] = pit->second;
	    
	    pit = receiving.erase(pit);
	  }
	  else
	    ++pit;
	}

	//update symbolic maps (but not send the results)
	if(temperature > 0.)
	  _update_symbolic_map(sym_val_map, received);
	else
	  _update_symbolic_map(sym_zpe_map, received);
	//
	//
      } // normal mode indices cycle

      // extract and clean completed receives
      std::map<_fg_t, double> received;
      for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
	if(pit->Test()) {
	  received[pit->first] = pit->second;
	    
	  pit = receiving.erase(pit);
	}
	else
	  ++pit;
      }

      //update symbolic maps and send the results to the master node
      if(temperature > 0.)
	_update_symbolic_map(sym_val_map, received, &sending);
      else
	_update_symbolic_map(sym_zpe_map, received, &sending);
      //
      //
    } // graph cycle

    // complete rest of receives
    std::map<_fg_t, double> received;
    for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
      pit->Wait();
	
      received[pit->first] = pit->second;
      
      pit = receiving.erase(pit);
    }

    //update symbolic maps
    if(temperature > 0.)
      _update_symbolic_map(sym_val_map, received, &sending);
    else
      _update_symbolic_map(sym_zpe_map, received, &sending);
    
    // complete rest of sends
    for(sent_t::iterator sit = sending.begin(); sit != sending.end(); ) {
      (*sit)->Wait();

      sit = sending.erase(sit);
    }
    
    if(temperature > 0. && sym_val_map.size() || temperature <= 0. && sym_zpe_map.size()) {
      ErrOut err_out;
      err_out << "some graphs in the symbolic map not resolved";
    } 

    MPI::COMM_WORLD.Send(0, 0, MPI::INT, MASTER, END_TAG);

    MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MASTER, END_TAG);
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void GraphExpansion::_centroid_driver_with_constrain (double temperature, const std::map<std::multiset<int>, double>& mmat) const
{
  const char funame [] = "GraphExpansion::centroid_driver_with_constrain: ";

  try {
    int    itemp;
    double dtemp;
    bool   btemp;

    // mpi stuff
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    if(!_is_driver(mpi_rank)) {
      ErrOut err_out;
      err_out << funame << "not a driver";
    }

    typedef std::list<SharedPointer<_gbase_t> >  sent_t;

    sent_t            sending;
    std::list<_gin_t> receiving;
    
    std::map<int, std::map<_mg_t, double> >                  sym_val_map;
    std::map<int, std::map<int, std::map<_mg_t, double> > >  sym_zpe_map;

    // calculation stuff
    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);

    while(1) {
      // mpi stuff
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, MASTER, GINDEX_REQUEST);

      const int gindex = _gnode_t(MASTER, GINDEX_TAG);

      if(_graph_data.size() == gindex)
	break;
    
      std::vector<_graph_t>::const_iterator graphit = _sorted_graph.begin() + gindex;

      const int vertex_size = graphit->vertex_size();

      std::vector<std::set<int> > vertex_map(vertex_size);
      //
      itemp = 0;
      //
      for(_graph_t::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {
	//
	if(git->size() != 2) {
	  ErrOut err_out;
	  err_out << funame << "bond should connect exactly two vertices: " << git->size();
	}
      
	int v = 0;
	//
	for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit, ++v)
	  //
	  if(!vertex_map[*bit].insert(2 * itemp + v).second) {
	    ErrOut err_out;
	    err_out << funame << "duplicated frequency index: " << 2 * itemp + v;
	  }
      }

      MultiIndexConvert corr_multi_index(graphit->size() * 2, _red_freq_index.size());

      for(long corr_li = 0; corr_li < corr_multi_index.size(); ++corr_li) {

	std::vector<int> corrin = corr_multi_index(corr_li);

	double potfac = 1.;

	btemp = false;
	//
	for(std::vector<std::set<int> >::const_iterator mit = vertex_map.begin(); mit != vertex_map.end(); ++mit) {
	  //
	  std::multiset<int> potex_sign;
	  //
	  for(std::set<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
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
	
	if(btemp)
	  continue;

	potfac /= (double)graphit->symmetry_factor();
	//
	if(vertex_size % 2)
	  //
	  potfac = -potfac;

	for(MultiIndex cmask(graphit->size(), 1); !cmask.end(); ++cmask) {// classical correlator mask cycle
	  //
	  double gfactor = potfac;
	  int    t_count = 0;
	
	  _fg_t mod_graph;

	  itemp = 0;
	  btemp = false;
	  for(_graph_t::const_iterator git = graphit->begin(); git != graphit->end(); ++git, ++itemp) {

	    int ci[2]; // individual bond normal mode indices
	    for(int i = 0; i < 2; ++i)
	      ci[i] = corrin[2 * itemp + i];
	
	    // quantum correlator
	    if(cmask[itemp]) {

	      if(ci[0] != ci[1]) {
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
	      if(bond.size() == 1) {
		//
		if(low_freq.find(fi) != low_freq.end()) {
		  //
		  gfactor /= 12. * temperature;
		}
		else if(temperature > 0.) {
		  //
		  gfactor /= 2. * _red_freq[fi] * tanh_factor[fi];
		}
		else {
		  //
		  gfactor /= 2. * _red_freq[fi];
		}
	      }
	      // add frequency index to the graph
	      else {
		//
		mod_graph[bond].insert(fi);
	      }	      
	    }
	    // centroid correction
	    else {
	      std::multiset<int> mi;
	      for(int i = 0; i < 2; ++i)
		mi.insert(ci[i]);

	      std::map<std::multiset<int>, double>::const_iterator mit = mmat.find(mi);
	    
	      if(ci[0] != ci[1]) {
		if(mmat.end() != mit) {
		  gfactor *= -mit->second;
		}
		else {
		  btemp = true;
		  break;
		}
	    
		for(int i = 0; i < 2; ++i) {
		  const int fi = _red_freq_index[ci[i]]; // reduced frequency index
	      
		  if(_red_freq[fi] > 0.)
		    gfactor /= _red_freq[fi] * _red_freq[fi];
		  else
		    gfactor /= -_red_freq[fi] * _red_freq[fi];
		}
	      }
	      else {
		const int fi = _red_freq_index[ci[0]]; // reduced frequency index

		if(_red_freq[fi] > 0.)
		  dtemp = _red_freq[fi] * _red_freq[fi];
		else
		  dtemp = -_red_freq[fi] * _red_freq[fi];

		if(low_freq.find(fi) != low_freq.end()) {
		  if(mmat.end() != mit) {
		    gfactor *= (1. - mit->second / dtemp) / dtemp;
		  }
		  else {
		    ErrOut err_out;
		    err_out << funame << "low frequency " fi << " does not have corresponding value in M matrix";
		  }
		}
		else if(mmat.end() != mit) {
		  gfactor *=  -mit->second / dtemp / dtemp;
		}
		else {
		  btemp = true;
		  break;
		}
	      }

	      if(temperature > 0.)
		gfactor *= temperature;
	      else
		++t_count;
	    }
	  }

	  if(btemp)
	    continue;
	
	  itemp = vertex_size - mod_graph.vertex_size();
	  if(itemp) { 
	    if(temperature > 0.)
	      gfactor /= std::pow(temperature, (double)itemp);
	    else
	      t_count -= itemp;
	  }

	  // graph factorization into connected graphs
	  _mg_t sym = mod_graph.factorize();
	  
	  if(temperature <=  0.)
	    for(_mg_t::const_iterator mit = sym.begin(); mit != sym.end(); ++mit)
	      t_count -= mit->second;

	  for(_mg_t::const_iterator sit = sym.begin(); sit != sym.end(); ++sit) {// modified graph factor cycle
	    sending.push_back(SharedPointer<_gbase_t>(new _gself_t(sit->first)));
	    sending.back()->isend(INT_SERV, INT_TAG);
	  
	    receiving.push_back(_gin_t());
	    receiving.back().irecv(INT_SERV, INT_TAG);
	  }
	  
	  if(temperature > 0.)
	    sym_val_map[gindex][sym] += gfactor;
	  else
	    sym_zpe_map[gindex][t_count][sym] += gfactor;
	} // correlator mask cycle
	
	for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	  if((*sit)->Test())
	    sit = sending.erase(sit);
	  else
	    ++sit;
      
	// extract and clean completed receives
	std::map<_fg_t, double> received;
	for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
	  if(pit->Test()) {
	    received[pit->first] = pit->second;
	      
	    pit = receiving.erase(pit);
	  }
	  else
	    ++pit;
	}

	//update symbolic maps
	if(temperature > 0.)
	  _update_symbolic_map(sym_val_map, received);
	else
	  _update_symbolic_map(sym_zpe_map, received);
	  
      }// normal mode indices cycle

      // extract and clean completed receives
      std::map<_fg_t, double> received;
      for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
	if(pit->Test()) {
	  received[pit->first] = pit->second;
	      
	  pit = receiving.erase(pit);
	}
	else
	  ++pit;
      }

      //update symbolic maps
      if(temperature > 0.)
	_update_symbolic_map(sym_val_map, received, &sending);
      else
	_update_symbolic_map(sym_zpe_map, received, &sending);
	  
    }  // graph cycle

    // complete rest of receives
    std::map<_fg_t, double> received;
    for(std::list<_gin_t>::iterator pit = receiving.begin(); pit != receiving.end(); ) {
      pit->Wait();
    
      received[pit->first] = pit->second;
      
      pit = receiving.erase(pit);
    }

    //update symbolic maps
    if(temperature > 0.)
      _update_symbolic_map(sym_val_map, received, &sending);
    else
      _update_symbolic_map(sym_zpe_map, received, &sending);
    
    // complete rest of sends
    for(sent_t::iterator sit = sending.begin(); sit != sending.end(); ) {
      (*sit)->Wait();
    
      sit = sending.erase(sit);
    }
  
    if(temperature > 0. && sym_val_map.size() || temperature <= 0. && sym_zpe_map.size()) {
      ErrOut err_out;
      err_out << "some graphs in the symbolic map not resolved";
    }

    MPI::COMM_WORLD.Send(0, 0, MPI::INT, MASTER, END_TAG);

    MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MASTER, END_TAG);
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

// perturbation graph theory initializer
void GraphExpansion::init (const std::vector<double>& freq, const _potex_t& pex)
{
  const char funame [] = "GraphExpansion::init: ";

  int itemp;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  if(!freq.size()) {
    ErrOut err_out;
    err_out << funame << "no frequencies";
  }
  
  if(pex.index_range() != freq.size()) {
    ErrOut err_out;
    err_out << funame << "potential expansion index range and frequencies number mismatch: " << pex.index_range() << ", " << freq.size();
  }
  
  // initialize frequencies
  _set_frequencies(freq);

  // initialize potential expansion
  _potex = pex;

  std::set<int> pex_rank_pool;
  //
  for(_potex_t::const_iterator pit = _potex.begin(); pit != _potex.end(); ++pit)
    pex_rank_pool.insert(pit->first.size());

  std::vector<int> pex_rank;
  for(std::set<int>::const_iterator pit = pex_rank_pool.begin(); pit != pex_rank_pool.end(); ++pit)
    if(*pit > 2)
      pex_rank.push_back(*pit);
  
  if(!pex_rank.size()) {
    ErrOut err_out;
    err_out << funame << "anharmonic terms in the potential expansion do not exist";
  }

  std::vector<int> glimit(pex_rank.size());
  for(int i = 0; i < pex_rank.size(); ++i)
    glimit[i] = 2 * bond_max / pex_rank[i];

  // initialize frequency adapted graph converter (glimit[0] - maximal number of vertices)
  _convert.init(glimit[0], _red_freq.size());

  IO::Marker funame_marker(funame);

  int unconnected = 0;
  int graph_count = 0;
  int   raw_count = 0;
  
  for(MultiIndex mi(glimit); !mi.end(); ++mi) {

    itemp = 0;
    for(int i = 0; i < pex_rank.size(); ++i)
      itemp += pex_rank[i] * mi.base()[i];

    if(!itemp || itemp % 2 || itemp / 2 > bond_max)
      continue;

    std::map<int, int> ginit;
    for(int i = 0; i < pex_rank.size(); ++i)
      if(mi.base()[i])
	ginit[pex_rank[i]] = mi.base()[i];

    if(!mpi_rank) {
      IO::log << IO::log_offset << "perturbation term:";
      for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
	IO::log << "  " << it->first << "(" << it->second << ")";
      IO::log << std::endl;
    }
    
    std::vector<int>               vertex_order;
    std::vector<int>               multi_perm_init;
    std::vector<std::vector<int> > vertex_perm_base;

    for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it) {
      // MultiPerm dimensions
      multi_perm_init.push_back(it->second);

      // vertex permutation base & vertex order
      vertex_perm_base.push_back(std::vector<int>(it->second));
      for(int i = 0; i < it->second; ++i) {
	vertex_perm_base.back()[i] = vertex_order.size();
	vertex_order.push_back(it->first);
      }
    }
  
    std::set<_graph_t> raw_pool = _raw_graph_generator(vertex_order);

    raw_count += raw_pool.size();
  
    itemp = 1;
    for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
      itemp *= Math::factorial(it->second);
    const int vertex_perm_size = itemp;

    while(1) {
      // erase not-connected graphs
      while(raw_pool.size() && !raw_pool.begin()->is_connected()) {
	raw_pool.erase(raw_pool.begin());
	++unconnected;
      }

      if(!raw_pool.size())
	break;
    
      ++graph_count;

      _graph_t ancor = *raw_pool.begin();
    
      if(!_graph_data.insert(ancor).second) {
	ErrOut err_out;
	err_out << funame << "graph already in the pool";
      }

      // remove all permutationally related graphs from the raw graph pool
      int perm_equal = 0;
      int perm_diff  = 0;
      for(MultiPerm multi_perm(multi_perm_init); !multi_perm.end(); ++multi_perm) {

	std::vector<int> vertex_perm(vertex_order.size());

	for(int i = 0; i < vertex_perm_base.size(); ++i)
	  for(int j = 0; j < vertex_perm_base[i].size(); ++j)
	    vertex_perm[vertex_perm_base[i][j]] = vertex_perm_base[i][multi_perm[i][j]];
      
	_graph_t perm_graph;
	for(_graph_t::const_iterator pit = ancor.begin(); pit != ancor.end(); ++pit) {
	  std::multiset<int> bond;
	  for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it)
	    bond.insert(vertex_perm[*it]);
	  perm_graph.insert(bond);
	}

	if(perm_graph == ancor)
	  ++perm_equal;

	std::set<_graph_t>::iterator git = raw_pool.find(perm_graph);
	if(git != raw_pool.end()) {
	  ++perm_diff;
	  raw_pool.erase(git);
	}
      }
    
      if(perm_diff * perm_equal != vertex_perm_size) {
	ErrOut err_out;
	err_out << funame << "numbers of permutationally equal and different graphs inconsistent: "
		<< perm_equal << ", "<< perm_diff;
      }

      if(perm_equal != ancor.vertex_symmetry()) {
	ErrOut err_out;
	err_out << funame << "permutational symmetry factors differ";
      }
    }
  }

  // sorting graph according to certain criterion (vertex size, for example)
  std::multimap<int, _graph_t> sort_map;
  for(std::set<_graph_t>::const_iterator graphit = _graph_data.begin(); graphit != _graph_data.end(); ++graphit)
    sort_map.insert(make_pair(graphit->vertex_size(), *graphit));

  _sorted_graph.resize(sort_map.size());
    
  itemp = 0;
  for(std::multimap<int, _graph_t>::const_reverse_iterator sit = sort_map.rbegin(); sit != sort_map.rend(); ++sit, ++itemp)
    _sorted_graph[itemp] = sit->second;


  if(!mpi_rank) {
    IO::log << "\n";
  
    IO::log << IO::log_offset << "number of raw graphs = " << raw_count << "\n";
    IO::log << IO::log_offset << "number of unconnected graphs = " << unconnected << "\n";
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
}

std::set<GraphExpansion::_graph_t> GraphExpansion::_raw_graph_generator (std::vector<int> vertex_order, int root)
{
  const char funame [] = "GraphExpansion::_raw_graph_generator: ";
  
  int itemp;

  if(::sum(vertex_order) % 2) {
    ErrOut err_out;
    err_out << funame << "odd number of connction points";
  }
  
  std::set<_graph_t> res;

  if(root < 0) {
    std::set<_graph_t> add;
    
    for(int i = 0; i < vertex_order.size(); ++i)
      if(vertex_order[i] > 0) {
	add = _raw_graph_generator(vertex_order, i);
	root = i;
	break;
      }

    if(!add.size())
      return res;

    for(std::set<_graph_t>::const_iterator at = add.begin(); at != add.end(); ++at) {
      
      std::vector<int> new_order = vertex_order;
      for(_graph_t::const_iterator pit = at->begin(); pit != at->end(); ++pit)
	for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it)
	  --new_order[*it];

      for(int i = 0; i < new_order.size(); ++i)
	if(new_order[i] < 0 || i <= root && new_order[i]) {
	  ErrOut err_out;
	  err_out << funame << "wrong number of connection points";
	}
      
      std::set<_graph_t> g  = _raw_graph_generator(new_order);

      if(g.size())
	for(std::set<_graph_t>::const_iterator git = g.begin(); git != g.end(); ++git) {
	  _graph_t gval = *git;
	  for(_graph_t::const_iterator pit = at->begin(); pit != at->end(); ++pit)
	    gval.insert(*pit);
	  res.insert(gval);
	}
      else
	res.insert(*at);
    }
  }
  // front addition
  else {
    if(!vertex_order[root])
      return res;

    --vertex_order[root];
  
    for(int i = 0; i < vertex_order.size(); ++i) {
      if(vertex_order[i] < 0 || i < root && vertex_order[i]) {
	ErrOut err_out;
	err_out << funame << "wrong number of connection points";
      }
      
      if(vertex_order[i] > 0) {
	--vertex_order[i];
	std::set<_graph_t> g  = _raw_graph_generator(vertex_order, root);
	++vertex_order[i];

	std::multiset<int> bond;
	bond.insert(root);
	bond.insert(i);

	if(g.size())
	  for(std::set<_graph_t>::const_iterator git = g.begin(); git != g.end(); ++git) {
	    _graph_t gval = *git;
	    gval.insert(bond);
	    res.insert(gval);
	  }
	else {
	  _graph_t gval;
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
    tanh_factor.resize(_red_freq.size());

    for(int i = 0; i < _red_freq.size(); ++i)
      if(_red_freq[i] > 0.) {
	tanh_factor[i] = std::tanh(_red_freq[i] / temperature / 2.);
      
	if(_red_freq[i] / temperature < low_freq_thresh)
	  res.insert(i);
      }
      else {
	tanh_factor[i] = std::tan(-_red_freq[i] / temperature / 2.);
	
	if(_red_freq[i] / temperature > -low_freq_thresh)
	  res.insert(i);
      }
  }
  else 
    tanh_factor.clear();

  return res;
}

/*
std::vector<double> GraphExpansion::_adjust_frequencies (double temperature, std::vector<double>& tanh_factor) const
{
  const char funame [] = "GraphExpansion::_adjust_frequencies: ";
  
  static const double eps = 1.e-3;

  std::vector<double> adj_freq(_red_freq.size());
  
  if(temperature <= 0.) {
    adj_freq = _red_freq;
  }
  else {
    for(int i = 0; i < _red_freq.size(); ++i)
      if(_red_freq[i] >= 0. && _red_freq[i] < eps * temperature)
	adj_freq[i] = eps * temperature;
      else if(_red_freq[i] <= 0. && _red_freq[i] > -eps * temperature)
	adj_freq[i] = -eps * temperature;
      else
	adj_freq[i] = _red_freq[i];

    tanh_factor.resize(_red_freq.size());
  
    for(int i = 0; i < adj_freq.size(); ++i)
      if(adj_freq[i] > 0.)
	tanh_factor[i] = std::tanh(adj_freq[i] / temperature / 2.);
      else
	tanh_factor[i] = std::tan(-adj_freq[i] / temperature / 2.);
  }

  return adj_freq;
}
*/

void GraphExpansion::_set_frequencies (std::vector<double> freq)
{
  const char funame [] = "GraphExpansion::_set_frequencies: ";
  
  static const double min_freq = Phys_const::incm;
  
  double dtemp;
  int    itemp;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  _red_freq_map.clear();

  // removing very low frequencies
  for(int f = 0; f < freq.size(); ++f)
    if(freq[f] >= 0. && freq[f] < min_freq)
      freq[f] =  min_freq;
    else if(freq[f] <= 0. && freq[f] > -min_freq)
      freq[f] = -min_freq;
    
  for(int f = 0; f < freq.size(); ++f) {
    itemp = 1;
    for(int i = 0; i < _red_freq_map.size(); ++i) {
      dtemp = freq[f] / freq[*_red_freq_map[i].begin()];

      if(1. - freq_tol < dtemp && dtemp < 1. + freq_tol) { 
	_red_freq_map[i].insert(f);
	itemp = 0;
	break;
      }
    }
    if(itemp) {
      _red_freq_map.push_back(std::set<int>());
      _red_freq_map.back().insert(f);
    }
  }

  _red_freq.resize(_red_freq_map.size());
  for(int i = 0; i < _red_freq_map.size(); ++i) {
    dtemp = 0.;
    itemp =  1;
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it) {
      if(it == _red_freq_map[i].begin() && freq[*it] < 0.)
	itemp = 0;

      if(itemp && freq[*it] < 0. || !itemp && freq[*it] > 0.) {
	ErrOut err_out;
	err_out << funame << "frequencies of different signs in one group";
      }
      
      if(itemp)
	dtemp += std::log(freq[*it]);
      else
	dtemp += std::log(-freq[*it]);
    }  
    dtemp /= (double)_red_freq_map[i].size();

    if(itemp)
      _red_freq[i] =  std::exp(dtemp);
    else
      _red_freq[i] = -std::exp(dtemp);
  }

  _red_freq_index.resize(freq.size());
  for(int i = 0; i < _red_freq_map.size(); ++i)
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it)
      _red_freq_index[*it] = i;

  if(!mpi_rank) {
  
    IO::log << IO::log_offset << "logarithmic frequency tolerance  = " << freq_tol << "\n"
	    << IO::log_offset << "reduced frequencies, 1/cm:";
    for(int i = 0; i < _red_freq.size(); ++i)
      IO::log << "   " << _red_freq[i] / Phys_const::incm;
    IO::log << std::endl;

    IO::log << IO::log_offset << "reduced frequencies groups: ";
    for(int i = 0; i < _red_freq_map.size(); ++i) {
      for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it) {
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
}

void GraphExpansion::_update_symbolic_map (std::map<int, std::map<_fg_t, double> >& sym_map,
					   //
					   const std::map<_fg_t, double>&           received,
					   //
					   std::list<SharedPointer<_gbase_t> >*     sending
					   )
{
  const char funame [] = "GraphExpansion::_update_symbolic_map: ";

  int itemp;

  typedef std::map<_fg_t, double> sv_t;
  typedef std::map<int, sv_t> sm_t;

  itemp = 0;
  for(sm_t::const_iterator smit = sym_map.begin(); smit != sym_map.end(); ++smit)
    itemp += smit->second.size();

  const int sym_map_size = itemp;

  sm_t new_sym_map;

  for(sm_t::const_iterator smit = sym_map.begin(); smit != sym_map.end(); ++smit) {
    sv_t sym_val;
    
    for(sv_t::const_iterator svit = smit->second.begin(); svit != smit->second.end(); ++svit) {
      _fg_t sym;
      double val = svit->second;

      sv_t::const_iterator fit = received.find(svit->first);
      if(fit != received.end())
	val  *= fit->second;
      else
	sym = svit->first;

      sym_val[sym] += val;
    }
      
    if(sending && sym_val.size() == 1) {

      sending->push_back(SharedPointer<_gbase_t>(new _gval_t(smit->first, sym_val.begin()->second)));
	
      sending->back()->isend(MASTER, GINDEX_TAG);

      sending->push_back(SharedPointer<_gbase_t>(new _gnode_t(sym_map_size)));
      
      sending->back()->isend(MASTER, GINDEX_TAG, 1);
    }
    else {
      new_sym_map[smit->first] = sym_val;
    }
  }
  sym_map = new_sym_map;
}

void GraphExpansion::_update_symbolic_map (std::map<int, std::map<_mg_t, double> >& sym_map,
					   //
					   const std::map<_fg_t, double>&           received,
					   //
					   std::list<SharedPointer<_gbase_t> >*     sending
					   )
{
  const char funame [] = "GraphExpansion::_update_symbolic_map: ";

  int itemp;

  typedef std::map<_mg_t, double> sv_t;
  typedef std::map<int, sv_t>    sm_t;

  itemp = 0;
  for(sm_t::const_iterator smit = sym_map.begin(); smit != sym_map.end(); ++smit)
    itemp += smit->second.size();

  const int sym_map_size = itemp;

  sm_t new_sym_map;
	    
  for(sm_t::const_iterator smit = sym_map.begin(); smit != sym_map.end(); ++smit) {
    sv_t sym_val;
    for(sv_t::const_iterator svit = smit->second.begin(); svit != smit->second.end(); ++svit) {
      _mg_t sym;
      double val = svit->second;
	      
      for(_mg_t::const_iterator mgit = svit->first.begin(); mgit != svit->first.end(); ++mgit) {
	std::map<_fg_t, double>::const_iterator fit = received.find(mgit->first);
		
	if(fit != received.end())
	  for(int i = 0; i < mgit->second; ++i)
	    val *=  fit->second;
	else
	  sym.insert(*mgit);
      }
	      
      sym_val[sym] += val;
    }

    if(sending && sym_val.size() == 1) {
      sending->push_back(SharedPointer<_gbase_t>(new _gval_t(smit->first, sym_val.begin()->second)));
      
      sending->back()->isend(MASTER, GINDEX_TAG);

      sending->push_back(SharedPointer<_gbase_t>(new _gnode_t(sym_map_size)));
      
      sending->back()->isend(MASTER, GINDEX_TAG, 1);
    }
    else {
      new_sym_map[smit->first] = sym_val;
    }
  }
  sym_map = new_sym_map;
}

void GraphExpansion::_update_symbolic_map (std::map<int, std::map<int, std::map<_mg_t, double> > >& sym_map,
					   //
					   const std::map<_fg_t, double>&                           received,
					   //
					   std::list<SharedPointer<_gbase_t> >*                     sending
					   )
{
  const char funame [] = "GraphExpansion::_update_symbolic_map: ";

  int  itemp;
  bool btemp;

  typedef std::map<_mg_t, double> sv_t;
  typedef std::map<int, sv_t>    is_t;
  typedef std::map<int, is_t>    sm_t;

  itemp = 0;
  for(sm_t::const_iterator smit = sym_map.begin(); smit != sym_map.end(); ++smit)
    for(is_t::const_iterator isit = smit->second.begin(); isit != smit->second.end(); ++isit)
      itemp += isit->second.size();

  const int sym_map_size = itemp;

  sm_t new_sym_map;
	    
  for(sm_t::const_iterator smit = sym_map.begin(); smit != sym_map.end(); ++smit) {
    is_t sym_val;

    for(is_t::const_iterator isit = smit->second.begin(); isit != smit->second.end(); ++isit) {
      for(sv_t::const_iterator svit = isit->second.begin(); svit != isit->second.end(); ++svit) {
	_mg_t sym;
	double val = svit->second;
	      
	for(_mg_t::const_iterator mgit = svit->first.begin(); mgit != svit->first.end(); ++mgit) {
	  std::map<_fg_t, double>::const_iterator fit = received.find(mgit->first);
		
	  if(fit != received.end())
	    for(int i = 0; i < mgit->second; ++i)
	      val *=  fit->second;
	  else
	    sym.insert(*mgit);
	}
	      
	sym_val[isit->first][sym] += val;
      }
    }
    
    btemp = true;
    std::map<int, double> res;

    if(sending) {
      for(is_t::const_iterator isit = sym_val.begin(); isit != sym_val.end(); ++isit) {
	if(isit->second.size() != 1) {
	  btemp = false;
	  break;
	}

	res[isit->first] = isit->second.begin()->second;
      }
    }
	
    if(sending && btemp) {

      sending->push_back(SharedPointer<_gbase_t>(new _gmap_t(smit->first, res)));
      
      sending->back()->isend(MASTER, GINDEX_TAG);

      sending->push_back(SharedPointer<_gbase_t>(new _gnode_t(sym_map_size)));
      
      sending->back()->isend(MASTER, GINDEX_TAG, 1);
    }
    else {
      new_sym_map[smit->first] = sym_val;
    }
  }
  sym_map = new_sym_map;
}

/*******************************************************************************************
 ********************************** COMMUNICAIION CLASS BASE *******************************
 *******************************************************************************************/

GraphExpansion::_gbase_t::~_gbase_t () 
{
  const char funame [] = "GraphExpansion::_gbase_t::~_gbase_t: ";
      
  if(IFREE != _state) {
    ErrOut err_out;
    err_out << funame << "mpi communication is not complete: ";

    if(_state == IRECV) {
      err_out << "waiting to receive";
    }
    else if(_state == ISEND) {
      err_out << "waiting to send";
    }
    else {
      err_out << "unknown state: " << _state;
    }
  }
}
    
void GraphExpansion::_gbase_t::isend (int node, int tag, int flag) const
{
  const char funame [] = "GraphExpansion::_gbase_t::isend: ";

  if(IFREE != _state) {
    ErrOut err_out;
    err_out << funame << "buffer is used for communication: " << _state;
  }

  if(!_buff.size())
    _buff.resize(BUFF_SIZE);

  int pos = _pack(_buff);

  if(!flag)
    MPI::COMM_WORLD.Isend(0, 0, MPI::INT, node, tag).Free();
  
  _request = MPI::COMM_WORLD.Isend(_buff, pos, MPI::PACKED, node, tag);

  _state = ISEND;
}

void GraphExpansion::_gbase_t::send (int node, int tag) const
{
  const char funame [] = "GraphExpansion::_gbase_t::send: ";

  if(IFREE != _state) {
    ErrOut err_out;
    err_out << funame << "buffer is used for communication: " << _state;
  }

  Array<char> alt_buff((int)BUFF_SIZE);
  
  int pos = _pack(alt_buff);
  
  MPI::COMM_WORLD.Send(alt_buff, pos, MPI::PACKED, node, tag);
}

void GraphExpansion::_gbase_t::irecv (int node, int tag) const
{
  const char funame [] = "GraphExpansion::_gbase_t::irecv: ";

  if(IFREE != _state) {
    ErrOut err_out;
    err_out << funame << "buffer is used for communication: " << _state;
  }

  if(!_buff.size())
    _buff.resize(BUFF_SIZE);

  _request = MPI::COMM_WORLD.Irecv(_buff, _buff.size(), MPI::PACKED, node, tag);

  _state = IRECV;
}
      
void GraphExpansion::_gbase_t::recv (int node, int tag)
{
  const char funame [] = "GraphExpansion::_gbase_t::recv: ";

  if(IFREE != _state) {
    ErrOut err_out;
    err_out << funame << "buffer is used for communication: " << _state;
  }

  Array<char> alt_buff((int)BUFF_SIZE);
  
  MPI::COMM_WORLD.Recv(alt_buff, alt_buff.size(), MPI::PACKED, node, tag);

  _unpack(alt_buff);
}

bool GraphExpansion::_gbase_t::Test ()
{
  const char funame [] = "GraphExpansion::_gbase_t::Test: ";
      
  if(IFREE == _state) {
    ErrOut err_out;
    err_out << funame << "test on ifree buffer";
  }
      
  if(_request.Test()) {
    if(IRECV == _state)
      _unpack(_buff);

    _state = IFREE;
    return true;
  }

  return false;
}
    
bool GraphExpansion::_gbase_t::Test (MPI::Status& stat)
{
  const char funame [] = "GraphExpansion::_gbase_t::Test: ";
      
  if(IFREE == _state) {
    ErrOut err_out;
    err_out << funame << "test on ifree buffer";
  }
      
  if(_request.Test(stat)) {
    if(IRECV == _state)
      _unpack(_buff);

    _state = IFREE;
    return true;
  }

  return false;
}
    
void GraphExpansion::_gbase_t::Wait ()
{
  const char funame [] = "GraphExpansion::_gbase_t::Wait: ";
      
  if(IFREE == _state) {
    ErrOut err_out;
    err_out << funame << "wait on ifree buffer";
  }
      
  _request.Wait();
      
  if(IRECV == _state)
    _unpack(_buff);
      
  _state = IFREE;
}
    
void GraphExpansion::_gbase_t::Wait (MPI::Status& stat)
{
  const char funame [] = "GraphExpansion::_gbase_t::Wait: ";
      
  if(IFREE == _state) {
    ErrOut err_out;
    err_out << funame << "wait on ifree buffer";
  }
      
  _request.Wait(stat);
  
  if(IRECV == _state)
    _unpack(_buff);
      
  _state = IFREE;
}

/*******************************************************************************************
 ***************************** COMMUNICAIION CLASS FOR GRAPH ITSELF ************************
 *******************************************************************************************/

int GraphExpansion::_gself_t::_pack (Array<char>& buff) const
{
  int itemp;
  
  _Convert::vec_t gconv = _convert((const _fg_t&)*this);

  itemp = gconv.size();

  int pos = 0;

  MPI::INT.Pack(&itemp, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  MPI::SHORT.Pack(gconv, gconv.size(), buff, buff.size(), pos, MPI::COMM_WORLD);

  return pos;
}

void GraphExpansion::_gself_t::_unpack (const Array<char>& buff)
{
  int itemp;
  
  _Convert::vec_t gconv;

  int pos = 0;

  MPI::INT.Unpack(buff, buff.size(), &itemp, 1, pos, MPI::COMM_WORLD);
    
  gconv.resize(itemp);
  MPI::SHORT.Unpack(buff, buff.size(), gconv, gconv.size(), pos, MPI::COMM_WORLD);

  (_fg_t&)*this = _convert(gconv);
}

/*******************************************************************************************
 *********************** COMMUNICAIION CLASS FOR GRAPH AND ITS VALUE ***********************
 *******************************************************************************************/

GraphExpansion::_gin_t::_gin_t (const _fg_t g, double d) : std::pair<_fg_t, double>(g, d) {}

int GraphExpansion::_gin_t::_pack (Array<char>& buff) const
{
  int itemp;
  
  _Convert::vec_t gconv = _convert(this->first);
  
  itemp = gconv.size();

  int pos = 0;

  MPI::INT.Pack(&itemp, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  MPI::SHORT.Pack(gconv, gconv.size(), buff, buff.size(), pos, MPI::COMM_WORLD);

  MPI::DOUBLE.Pack(&this->second, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  return pos;
}

void GraphExpansion::_gin_t::_unpack (const Array<char>& buff)
{
  int    itemp;
  double dtemp;

  _Convert::vec_t gconv;

  int pos = 0;

  MPI::INT.Unpack(buff, buff.size(), &itemp, 1, pos, MPI::COMM_WORLD);
    
  gconv.resize(itemp);
  MPI::SHORT.Unpack(buff, buff.size(), gconv, gconv.size(), pos, MPI::COMM_WORLD);

  this->first = _convert(gconv);

  MPI::DOUBLE.Unpack(buff, buff.size(), &this->second, 1, pos, MPI::COMM_WORLD);
}

/*******************************************************************************************
 ************************* COMMUNICAIION CLASS FOR TWO GRAPHS ******************************
 *******************************************************************************************/

int GraphExpansion::_gpair_t::_pack (Array<char>& buff) const
{
  int itemp;
  
  _Convert::vec_t gconv;
  
  int pos = 0;

  gconv = _convert(this->first);
  itemp = gconv.size();

  MPI::INT.Pack(&itemp, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  MPI::SHORT.Pack(gconv, gconv.size(), buff, buff.size(), pos, MPI::COMM_WORLD);

  gconv = _convert(this->second);
  itemp = gconv.size();

  MPI::INT.Pack(&itemp, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  MPI::SHORT.Pack(gconv, gconv.size(), buff, buff.size(), pos, MPI::COMM_WORLD);

  return pos;
}

void GraphExpansion::_gpair_t::_unpack (const Array<char>& buff)
{
  int    itemp;

  _Convert::vec_t gconv;

  int pos = 0;

  MPI::INT.Unpack(buff, buff.size(), &itemp, 1, pos, MPI::COMM_WORLD);
    
  gconv.resize(itemp);
  MPI::SHORT.Unpack(buff, buff.size(), gconv, gconv.size(), pos, MPI::COMM_WORLD);

  this->first = _convert(gconv);

  MPI::INT.Unpack(buff, buff.size(), &itemp, 1, pos, MPI::COMM_WORLD);
    
  gconv.resize(itemp);
  MPI::SHORT.Unpack(buff, buff.size(), gconv, gconv.size(), pos, MPI::COMM_WORLD);

  this->second = _convert(gconv);
}

/*******************************************************************************************
 ************************* COMMUNICAIION CLASS FOR INDEX-VALUE *****************************
 *******************************************************************************************/

int GraphExpansion::_gval_t::_pack (Array<char>& buff) const
{
  int pos = 0;

  MPI::INT.Pack   (&this->first,  1, buff, buff.size(), pos, MPI::COMM_WORLD);
  MPI::DOUBLE.Pack(&this->second, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  return pos;
}

void GraphExpansion::_gval_t::_unpack (const Array<char>& buff)
{
  int pos = 0;

  MPI::INT.Unpack   (buff, buff.size(), &this->first,  1, pos, MPI::COMM_WORLD);
  MPI::DOUBLE.Unpack(buff, buff.size(), &this->second, 1, pos, MPI::COMM_WORLD);
}

/*******************************************************************************************
 ************************* COMMUNICAIION CLASS FOR INTEGER *********************************
 *******************************************************************************************/

int GraphExpansion::_gnode_t::_pack (Array<char>& buff) const
{
  int pos = 0;

  MPI::INT.Pack(&_value,  1, buff, buff.size(), pos, MPI::COMM_WORLD);

  return pos;
}

void GraphExpansion::_gnode_t::_unpack (const Array<char>& buff)
{
  int pos = 0;

  MPI::INT.Unpack(buff, buff.size(), &_value,  1, pos, MPI::COMM_WORLD);
}

/*******************************************************************************************
 ******************* COMMUNICAIION CLASS ZERO TEMPERATURE GRAPH EXPANSION ******************
 *******************************************************************************************/

int GraphExpansion::_gmap_t::_pack (Array<char>& buff) const
{
  int itemp;
  
  int pos = 0;

  MPI::INT.Pack(&this->first, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  itemp = this->second.size();

  MPI::INT.Pack(&itemp, 1, buff, buff.size(), pos, MPI::COMM_WORLD);

  for(std::map<int, double>::const_iterator mit = this->second.begin(); mit != this->second.end(); ++mit) {
    MPI::INT.Pack(&mit->first, 1, buff, buff.size(), pos, MPI::COMM_WORLD);
    
    MPI::DOUBLE.Pack(&mit->second, 1, buff, buff.size(), pos, MPI::COMM_WORLD);
  }
  
  return pos;
}

void GraphExpansion::_gmap_t::_unpack (const Array<char>& buff)
{
  const char funame [] = "GraphExpansion::_gmap_t::_unpack: ";
  
  int    itemp;
  double dtemp;

  this->second.clear();
  
  int pos = 0;

  MPI::INT.Unpack(buff, buff.size(), &this->first,  1, pos, MPI::COMM_WORLD);

  MPI::INT.Unpack(buff, buff.size(), &itemp,        1, pos, MPI::COMM_WORLD);

  const int map_size = itemp;
  
  for(int i = 0; i < map_size; ++i) {
    MPI::INT.Unpack(buff, buff.size(), &itemp,    1, pos, MPI::COMM_WORLD);
    
    MPI::DOUBLE.Unpack(buff, buff.size(), &dtemp, 1, pos, MPI::COMM_WORLD);

    if(!this->second.insert(std::make_pair(itemp, dtemp)).second) {
      ErrOut err_out;
      err_out << funame << "cannot insert int-double pair: " << itemp << ", " << dtemp << ": available values:\n";
      for(std::map<int, double>::const_iterator mit = this->second.begin(); mit != this->second.end(); ++mit)
	err_out << mit->first << "  " << mit->second << "\n";
      err_out << "\n";
    }
  }
}

/*******************************************************************************************
 ************************************ WORKING NODE PROCESS *********************************
 *******************************************************************************************/

void GraphExpansion::_node_work (double temperature) const
{
  const char funame [] = "GraphExpansion::_node_work: ";

  try {
    int    itemp;
    double dtemp;

    double std_time = 0;
    double zpe_time = 0.;
    double sum_time = 0.;

    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);
  
    typedef std::list<SharedPointer<_gbase_t> > sent_t;
    sent_t sending;
  
    // main loop
    while(1) {
      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);

      const int tag  = stat.Get_tag();
      const int from = stat.Get_source();

      // clean completed sends
      for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	if((*sit)->Test())
	  sit = sending.erase(sit);
	else
	  ++sit;

      // end of work
      if(tag ==  END_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "end request should come from master: " << from;
	}

	break;
      }

      const _fg_t graph = _gself_t(from, tag);

      // zpe calculation
      if(tag == ZPE_TAG) {
	if(temperature > 0. && from != ZPE_SERV || temperature <= 0. && from != INT_SERV) {
	  ErrOut err_out;
	  err_out << funame << "zpe factor calculation request comes not from zpe server: " << from;
	}

	double start_time = std::clock();

	dtemp = graph.zpe_factor(_red_freq, temperature, tanh_factor, low_freq);

	zpe_time += std::clock() - start_time;

	sending.push_back(SharedPointer<_gbase_t>(new _gin_t(graph, dtemp)));
      
	sending.back()->isend(from, tag);

	MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_RELEASE).Free();
      }
      // fourier sum calculation
      else if(tag == SUM_TAG) {
	if(from != SUM_SERV) {
	  ErrOut err_out;
	  err_out << funame << "fourier sum calculation request does not come from  the server: " << from;
	}

	if(temperature <= 0.) {
	  ErrOut err_out;
	  err_out << funame << "fourier sum calculation request for zero temperature";
	}
      
	double start_time = std::clock();

	dtemp = graph.fourier_sum(_red_freq, temperature, low_freq);

	sum_time += std::clock() - start_time;

	sending.push_back(SharedPointer<_gbase_t>(new _gin_t(graph, dtemp)));
      
	sending.back()->isend(from, tag);

	MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_RELEASE).Free();
      }
      // graph reduction to standard form
      else if(tag == STD_TAG) {
	if(temperature <= 0. && from != INT_SERV || temperature > 0. && !_is_server(from)) {
	  ErrOut err_out;
	  err_out << funame << "standard graph calculation request comes not from a server: " << from;
	}


	double start_time = std::clock();

	_fg_t std_graph = *graph.perm_pool().begin();

	std_time += std::clock() - start_time;

	sending.push_back(SharedPointer<_gbase_t>(new _gpair_t(graph, std_graph)));
      
	sending.back()->isend(from, tag);

	MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_RELEASE).Free();
      }
      // wrong tag
      else {
	ErrOut err_out;
	err_out << funame << "unknown tag: " << tag;
      }
      //
      //
    } // main loop

    // some checking
    if(sending.size()) {
      ErrOut err_out;
      err_out << funame << "some sends not completed";
    }

    MPI::COMM_WORLD.Send(&std_time, 1, MPI::DOUBLE, MASTER, STAT_TAG);
    MPI::COMM_WORLD.Send(&zpe_time, 1, MPI::DOUBLE, MASTER, STAT_TAG);

    if(temperature > 0.)
      MPI::COMM_WORLD.Send(&sum_time, 1, MPI::DOUBLE, MASTER, STAT_TAG);

  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void GraphExpansion::_int_server (double temperature) const
{
  const char funame [] = "GraphExpansion::_int_server: ";

  try {
    int    itemp;
    double dtemp;
  
    std::vector<double> tanh_factor;
    std::set<int> low_freq = _low_freq_set(temperature, tanh_factor);
  
    _db_t int_data;

    typedef std::map<_fg_t, std::multiset<int> > raw_t;
    typedef std::map<_fg_t, _red_t>              std_t;

    raw_t raw_pool; // <raw graph, send #>
    std_t std_pool; // <std graph, set<raw graph>>

    std::set<_fg_t> zpe_pool;
    std::set<_fg_t> red_pool;
  
    std::list<std::pair<_fg_t, int> >  waiting; // waiting queue; arguments: <graph, job type>

    typedef std::list<SharedPointer<_gbase_t> > sent_t;
    sent_t sending; // messages posted to send
    
    long calc_count = 0;
    long read_count = 0;
    long conv_count = 0;

    while(1) {
      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);

      const int tag  = stat.Get_tag();
      const int from = stat.Get_source();

      // remove completed sends
      for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	if((*sit)->Test())
	  sit = sending.erase(sit);
	else
	  ++sit;
    
      // end of work
      if(tag ==  END_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "end request should come from master: " << from;
	}
      
	break;
      }

      //
      // getting working node
      //
      if(tag == NODE_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "received node service from unexpected source: " << from; 
	}
      
	if(!waiting.size()) {
	  ErrOut err_out;
	  err_out << funame << "no job for working node";
	}
	
	const int node = _gnode_t(from, tag);
      
	sending.push_back(SharedPointer<_gbase_t>(new _gself_t(waiting.begin()->first)));
      
	sending.back()->isend(node, waiting.begin()->second);
      
	waiting.erase(waiting.begin());
      }
      //
      // integral request
      //
      else if(tag == INT_TAG) {
      
	const _fg_t graph = _gself_t(from, tag);

	if(raw_pool.find(graph) == raw_pool.end()) {
	  ++conv_count;

	  waiting.push_back(std::make_pair(graph, STD_TAG));
	
	  // node request
	  MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_TAG).Free();
	}

	raw_pool[graph].insert(from);
      }
      //
      // standard graph result
      //
      else if(tag == STD_TAG) {
	if(!_is_work_node(from)) {
	  ErrOut err_out;
	  err_out << funame << "standard graph should come from working node: " << from;
	}
	  
	_gpair_t graph_pair(from, tag);

	const _fg_t& raw_graph = graph_pair.first;
	const _fg_t& std_graph = graph_pair.second;
	
	if(raw_pool.find(raw_graph) == raw_pool.end()) {
	  ErrOut err_out;
	  err_out << funame << "raw graph is not in the pool";
	}

	if(!std_pool[std_graph].raw_group.insert(raw_graph).second) {
	  ErrOut err_out;
	  err_out << funame << "raw graph already has corresponding standard one";
	}

	std_t::iterator sgit = std_pool.find(std_graph);

	if(sgit->second.raw_group.size() == 1) {
	
	  _db_t::const_iterator idit = int_data.find(_convert(std_graph));

	  // graph in the database
	  if(idit != int_data.end()) {
	    ++read_count;

	    const std::set<_fg_t>& rg = sgit->second.raw_group;
	  
	    for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {

	      raw_t::iterator rit = raw_pool.find(*git);

	      if(rit == raw_pool.end()) {
		ErrOut err_out;
		err_out << funame << "raw graph not in the pool";
	      }
	  
	      for(std::multiset<int>::const_iterator nit = rit->second.begin(); nit != rit->second.end(); ++nit) {
		sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, idit->second)));
	  
		sending.back()->isend(*nit, INT_TAG, 1);
	      }
	  
	      raw_pool.erase(rit);
	    }
	    std_pool.erase(sgit);
	  }
	  // graph not in the database at zero temperature
	  else if(temperature <= 0.) {
	    ++calc_count;

	    waiting.push_back(std::make_pair(std_graph, ZPE_TAG));
	
	    // node request
	    MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_TAG).Free();
	  }
	  // non-zero temperature
	  else {
	    ++calc_count;

	    _mg_t zpe_group;
	    _fg_t red_graph = std_graph.reduce(_red_freq, temperature, tanh_factor, zpe_group, low_freq);

	    sgit->second.red_graph = red_graph;
	    sgit->second.zpe_group = zpe_group;

	    if(!red_graph.size()) {
	      sgit->second.value /= temperature;
	    }
	    else if(red_pool.insert(red_graph).second) {
	      sending.push_back(SharedPointer<_gbase_t>(new _gself_t(red_graph)));
	      
	      sending.back()->isend(SUM_SERV, SUM_REQUEST);
	    }
	    
	    for(_mg_t::const_iterator zgit = zpe_group.begin(); zgit != zpe_group.end(); ++zgit) {
	      if(zpe_pool.insert(zgit->first).second) {
		sending.push_back(SharedPointer<_gbase_t>(new _gself_t(zgit->first)));

		sending.back()->isend(ZPE_SERV, ZPE_REQUEST);
	      }
	    }
	  }
	}
      }
      //
      // zpe result
      //
      else if(tag == ZPE_TAG) {
	if(temperature <= 0.) {
	  if(!_is_work_node(from)) {
	    ErrOut err_out;
	    err_out << funame << "at zero temperature zpe result should come from working node: " << from;
	  }
	
	  _gin_t gin(from, tag);

	  const _fg_t& std_graph = gin.first;

	  std_t::iterator sgit = std_pool.find(std_graph);
      
	  if(sgit == std_pool.end()) {
	    ErrOut err_out;
	    err_out << funame << "standard graph is not in the pool: " << std_graph;
	  }
	
	  _Convert::vec_t gconv = _convert(std_graph);
	
	  if(int_data.find(gconv) != int_data.end()) {
	    ErrOut err_out;
	    err_out << funame << "graph value is already in the database";
	  }

	  int_data[gconv] = gin.second;

	  const std::set<_fg_t>& rg = sgit->second.raw_group;
	
	  for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {

	    raw_t::iterator rit = raw_pool.find(*git);
	
	    if(rit == raw_pool.end()) {
	      ErrOut err_out;
	      err_out << funame << "raw graph not in the pool";
	    }
	  
	    for(std::multiset<int>::const_iterator nit = rit->second.begin(); nit != rit->second.end(); ++nit) {
	      sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, gin.second)));
	    
	      sending.back()->isend(*nit, INT_TAG, 1);
	    }
	    raw_pool.erase(rit);
	  }  
	  std_pool.erase(sgit);
	}
	// positive temperature      
	else {
	  if(from != ZPE_SERV) {
	    ErrOut err_out;
	    err_out << funame << "at positive temperature zpe result should come from the server: " << from;
	  }

	  _gin_t gin(from, tag);

	  const _fg_t& zpe_graph = gin.first;

	  //std::cout << funame << "zpe value = " << gin.second << std::endl;

	  if(zpe_pool.find(zpe_graph) == zpe_pool.end()) {
	    ErrOut err_out;
	    err_out << funame << "zpe graph is not found in zpe pool";
	  }
	  zpe_pool.erase(zpe_graph);

	  std_t new_std_pool;
	
	  for(std_t::iterator sgit = std_pool.begin(); sgit != std_pool.end(); ++sgit) {
	    std::map<_fg_t, int>::iterator zgit = sgit->second.zpe_group.find(zpe_graph);

	    // resolving zpe graph symbol and removing it from the group
	    if(zgit != sgit->second.zpe_group.end()) {
	      for(int i = 0; i < zgit->second; ++i)
		sgit->second.value *= gin.second;
	    
	      sgit->second.zpe_group.erase(zgit);
	    }

	    // sending the result to the driver and removing the standard graph and associated with it group of raw graphs
	    if(!sgit->second.zpe_group.size() && !sgit->second.red_graph.size()) {
	    
	      _Convert::vec_t gconv = _convert(sgit->first);
	
	      if(int_data.find(gconv) != int_data.end()) {
		ErrOut err_out;
		err_out << funame << "graph value is already in the database";
	      }

	      int_data[gconv] = sgit->second.value;

	      const std::set<_fg_t>& rg = sgit->second.raw_group;
	
	      for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {
	      
		raw_t::iterator rit = raw_pool.find(*git);
	
		if(rit == raw_pool.end()) {
		  ErrOut err_out;
		  err_out << funame << "raw graph not in the pool";
		}
	  
		for(std::multiset<int>::const_iterator nit = rit->second.begin(); nit != rit->second.end(); ++nit) {
		  sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, sgit->second.value)));
	    
		  sending.back()->isend(*nit, INT_TAG, 1);
		}
		raw_pool.erase(rit);
	      }  
	    }
	    else
	      new_std_pool.insert(*sgit);
	  }
	  std_pool = new_std_pool;
	}
      }
      //
      // fourier sum result
      //
      else if(tag == SUM_TAG) {
	if(temperature <= 0.) {
	  ErrOut err_out;
	  err_out << funame << "fourier sum calculation only at positive temperature";
	}

	if(from != SUM_SERV) {
	  ErrOut err_out;
	  err_out << funame << "fourier sum result should come from the server: " << from;
	}
      
	_gin_t gin(from, tag);
      
	const _fg_t& red_graph = gin.first;

	if(red_pool.find(red_graph) == red_pool.end()) {
	  ErrOut err_out;
	  err_out << funame << "red graph is not found in the pool";
	}
	red_pool.erase(red_graph);

	std_t new_std_pool;
	
	for(std_t::iterator sgit = std_pool.begin(); sgit != std_pool.end(); ++sgit) {

	  // resolving reduced graph symbol and clearing it
	  if(red_graph == sgit->second.red_graph) {
	    sgit->second.value *= gin.second;
	    
	    sgit->second.red_graph.clear();
	  }

	  // sending the result to the driver and removing the standard graph and associated with it group of raw graphs
	  if(!sgit->second.zpe_group.size() && !sgit->second.red_graph.size()) {
	    
	    _Convert::vec_t gconv = _convert(sgit->first);
	
	    if(int_data.find(gconv) != int_data.end()) {
	      ErrOut err_out;
	      err_out << funame << "graph value is already in the database";
	    }

	    int_data[gconv] = sgit->second.value;

	    const std::set<_fg_t>& rg = sgit->second.raw_group;
	
	    for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {
	      
	      raw_t::iterator rit = raw_pool.find(*git);
	
	      if(rit == raw_pool.end()) {
		ErrOut err_out;
		err_out << funame << "raw graph not in the pool";
	      }
	  
	      for(std::multiset<int>::const_iterator nit = rit->second.begin(); nit != rit->second.end(); ++nit) {
		sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, sgit->second.value)));
	      
		sending.back()->isend(*nit, INT_TAG, 1);
	      }

	      raw_pool.erase(rit);
	    }  
	  }
	  else
	    new_std_pool.insert(*sgit);
	}
	std_pool = new_std_pool;
      }
      //
      // wrong tag
      //
      else {
	ErrOut err_out;
	err_out << funame << "unknown tag: " << tag;
      }
    }

    // some checking
    if(sending.size()) {
      ErrOut err_out;
      err_out << funame << "some sends not completed";
    }

    if(raw_pool.size()) {
      ErrOut err_out;
      err_out << funame << "raw graph pool is not empty";
    }

    if(std_pool.size()) {
      ErrOut err_out;
      err_out << funame << "standard graph pool is not empty";
    }

    if(zpe_pool.size()) {
      ErrOut err_out;
      err_out << funame << "there are some unresolved zpe graphs";
    }

    if(red_pool.size()) {
      ErrOut err_out;
      err_out << funame << "there are some unresolved fourier sum graphs";
    }

    if(waiting.size()) {
      ErrOut err_out;
      err_out << funame << "there are pending jobs";
    }

    MPI::COMM_WORLD.Send(&calc_count, 1, MPI::LONG, MASTER, STAT_TAG);
    MPI::COMM_WORLD.Send(&read_count, 1, MPI::LONG, MASTER, STAT_TAG);
    MPI::COMM_WORLD.Send(&conv_count, 1, MPI::LONG, MASTER, STAT_TAG);
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void GraphExpansion::_zpe_server (double temperature) const
{
  const char funame [] = "GraphExpansion::_zpe_server: ";

  try {
    int    itemp;
    double dtemp;
  
    _db_t zpe_data;

    typedef std::map<_fg_t, std::multiset<int> > raw_t;
    typedef std::map<_fg_t, std::set<_fg_t> >    std_t;

    raw_t raw_pool; // <raw graph, send #>
    std_t std_pool; // <std graph, <raw graph>>

    std::list<std::pair<_fg_t, int> >  waiting; // queue; arguments: <graph, job type>

    typedef std::list<SharedPointer<_gbase_t> > sent_t;
    sent_t sending;
    
    long calc_count = 0;
    long read_count = 0;
    long conv_count = 0;

    while(1) {
      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);

      const int tag  = stat.Get_tag();
      const int from = stat.Get_source();

      // remove completed sends
      for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	if((*sit)->Test())
	  sit = sending.erase(sit);
	else
	  ++sit;
    
      // end of work
      if(tag ==  END_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "end request should come from master: " << from;
	}
      
	break;
      }

      if(temperature <= 0.) {
	ErrOut err_out;
	err_out << funame << "at zero temperature integral server gets zpe value directly from the working node";
      }
    
      //
      // getting working node
      //
      if(tag == NODE_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "received node service from unexpected source: " << from; 
	}
      
	if(!waiting.size()) {
	  ErrOut err_out;
	  err_out << funame << "no job for working node";
	}
	
	const int node = _gnode_t(from, tag);
      
	sending.push_back(SharedPointer<_gbase_t>(new _gself_t(waiting.begin()->first)));

	sending.back()->isend(node, waiting.begin()->second);

	waiting.erase(waiting.begin());
      }
      //
      // zpe request
      //
      else if(tag == ZPE_REQUEST) {
	if(from != INT_SERV) {
	  ErrOut err_out;
	  err_out << funame << "request should come from the integral server: " << from;
	}

	_gself_t gself(from, tag);
      
	const _fg_t& graph = gself;
      
	if(raw_pool.find(graph) == raw_pool.end()) {
	  ++conv_count;

	  waiting.push_back(std::make_pair(graph, STD_TAG));
      
	  // node request
	  MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_TAG).Free();
	}

	raw_pool[graph].insert(from);
      }
      //
      // standard graph result
      //
      else if(tag == STD_TAG) {
	if(!_is_work_node(from)) {
	  ErrOut err_out;
	  err_out << funame << "standard graph should come from working node: " << from;
	}
	  
	_gpair_t graph_pair(from, tag);

	const _fg_t& raw_graph = graph_pair.first;
	const _fg_t& std_graph = graph_pair.second;
	
	if(raw_pool.find(raw_graph) == raw_pool.end()) {
	  ErrOut err_out;
	  err_out << funame << "raw graph is not in the pool";
	}

	if(!std_pool[std_graph].insert(raw_graph).second) {
	  ErrOut err_out;
	  err_out << funame << "raw graph already has corresponding standard one";
	}
      
	if(std_pool[std_graph].size() == 1) {

	  _db_t::const_iterator zdit = zpe_data.find(_convert(std_graph));

	  // graph in the database
	  if(zdit != zpe_data.end()) {
	    ++read_count;

	    std::set<_fg_t>& rg = std_pool[std_graph];
	  
	    for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {
	    
	      raw_t::iterator rit = raw_pool.find(*git);

	      if(rit == raw_pool.end()) {
		ErrOut err_out;
		err_out << funame << "raw graph is not in the pool";
	      }
	    
	      for(std::multiset<int>::const_iterator mit = rit->second.begin();
		  mit != rit->second.end(); ++mit) {
		sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, zdit->second)));
	    
		sending.back()->isend(*mit, ZPE_TAG);
	      }
	      raw_pool.erase(rit);
	    }
	    std_pool.erase(std_graph);
	  }
	  // graph not in the database
	  else {
	    ++calc_count;

	    waiting.push_back(std::make_pair(std_graph, ZPE_TAG));
	
	    // node request
	    MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_TAG).Free();
	  }
	}
      }
      //
      // zpe result
      //
      else if(tag == ZPE_TAG) {
	if(!_is_work_node(from)) {
	  ErrOut err_out;
	  err_out << funame << "zpe result should come from working node: " 
		  << from;
	}

	_gin_t gin(from, tag);

	const _fg_t& graph = gin.first;

	std_t::iterator sgit = std_pool.find(graph);
      
	if(sgit == std_pool.end()) {
	  ErrOut err_out;
	  err_out << funame << "standard graph is not in the pool";
	}
	
	_Convert::vec_t gconv = _convert(graph);
	
	if(zpe_data.find(gconv) != zpe_data.end()) {
	  ErrOut err_out;
	  err_out << funame << "graph value is already in the database";
	}

	zpe_data[gconv] = gin.second;

	const std::set<_fg_t>& rg = sgit->second;
	
	for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {
	
	  raw_t::iterator rit = raw_pool.find(*git);
	
	  if(rit == raw_pool.end()) {
	    ErrOut err_out;
	    err_out << funame << "raw graph is not in the pool";
	  }
	    
	  for(std::multiset<int>::const_iterator mit = rit->second.begin(); 
	      mit != rit->second.end(); ++mit) {
	    sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, gin.second)));
	    
	    sending.back()->isend(*mit, ZPE_TAG);
	  }
	  raw_pool.erase(rit);
	}
	std_pool.erase(sgit);
      }
      //
      // wrong tag
      //
      else {
	ErrOut err_out;
	err_out << funame << "unknown tag: " << tag;
      }
    }

    // some checking
    if(sending.size()) {
      ErrOut err_out;
      err_out << funame << "some sends not completed";
    }

    if(raw_pool.size()) {
      ErrOut err_out;
      err_out << funame << "raw graph pool is not empty";
    }

    if(std_pool.size()) {
      ErrOut err_out;
      err_out << funame << "standard graph pool is not empty";
    }

    if(waiting.size()) {
      ErrOut err_out;
      err_out << funame << "there are pending jobs";
    }

    if(temperature > 0.) {
      MPI::COMM_WORLD.Send(&calc_count, 1, MPI::LONG, MASTER, STAT_TAG);
      MPI::COMM_WORLD.Send(&read_count, 1, MPI::LONG, MASTER, STAT_TAG);
      MPI::COMM_WORLD.Send(&conv_count, 1, MPI::LONG, MASTER, STAT_TAG);
    }
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void GraphExpansion::_sum_server (double temperature) const
{
  const char funame [] = "GraphExpansion::_sum_server: ";

  try {
    int    itemp;
    double dtemp;
  
    _db_t sum_data;

    typedef std::map<_fg_t, std::multiset<int> > raw_t;
    typedef std::map<_fg_t, std::set<_fg_t> >    std_t;

    raw_t raw_pool; // <raw graph, send #>
    std_t std_pool; // <std graph, <raw graph>>

    std::list<std::pair<_fg_t, int> >  waiting; // queue; arguments: <graph, job type>

    typedef std::list<SharedPointer<_gbase_t> > sent_t;
    sent_t sending;
    
    long calc_count = 0;
    long read_count = 0;
    long conv_count = 0;

    while(1) {
      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);

      const int tag  = stat.Get_tag();
      const int from = stat.Get_source();

      // remove completed sends
      for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	if((*sit)->Test())
	  sit = sending.erase(sit);
	else
	  ++sit;
    
      // end of work
      if(tag ==  END_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "should come from master: " << from;
	}
      
	break;
      }

      if(temperature <= 0.) {
	ErrOut err_out;
	err_out << funame << "zero temperature fourier sum request";
      }
      
      //
      // getting working node
      //
      if(tag == NODE_TAG) {
	if(from != MASTER) {
	  ErrOut err_out;
	  err_out << funame << "received node service is not from node server: " << from; 
	}
      
	if(!waiting.size()) {
	  ErrOut err_out;
	  err_out << funame << "no job for working node";
	}
	
	const int node = _gnode_t(from, tag);
      
	sending.push_back(SharedPointer<_gbase_t>(new _gself_t(waiting.begin()->first)));

	sending.back()->isend(node, waiting.begin()->second);
      
	waiting.erase(waiting.begin());
      }
      //
      // sum request
      //
      else if(tag == SUM_REQUEST) {
	if(from != INT_SERV) {
	  ErrOut err_out;
	  err_out << funame << "request should come from the integral server: " << from;
	}

	_gself_t gself(from, tag);
      
	const _fg_t& graph = gself;
      
	if(raw_pool.find(graph) == raw_pool.end()) {
	  ++conv_count;

	  waiting.push_back(std::make_pair(graph, STD_TAG));
      
	  // node request
	  MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_TAG).Free();
	}

	raw_pool[graph].insert(from);
      }
      //
      // standard graph result
      //
      else if(tag == STD_TAG) {
	if(!_is_work_node(from)) {
	  ErrOut err_out;
	  err_out << funame << "standard graph result should come from working node: " << from;
	}
	  
	_gpair_t graph_pair(from, tag);

	const _fg_t& raw_graph = graph_pair.first;
	const _fg_t& std_graph = graph_pair.second;
	
	if(raw_pool.find(raw_graph) == raw_pool.end()) {
	  ErrOut err_out;
	  err_out << funame << "raw graph is not in the pool";
	}

	if(!std_pool[std_graph].insert(raw_graph).second) {
	  ErrOut err_out;
	  err_out << funame << "raw graph already has corresponding standard one";
	}
      
	if(std_pool[std_graph].size() == 1){

	  _db_t::const_iterator sdit = sum_data.find(_convert(std_graph));

	  // graph in the database
	  if(sdit != sum_data.end()) {
	    ++read_count;

	    std::set<_fg_t>& rg = std_pool[std_graph];
	  
	    for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {
	    
	      raw_t::iterator rit = raw_pool.find(*git);

	      if(rit == raw_pool.end()) {
		ErrOut err_out;
		err_out << funame << "raw graph is not in the pool";
	      }
	    
	      for(std::multiset<int>::const_iterator mit = rit->second.begin(); 
		  mit != rit->second.end(); ++mit) {
		sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, sdit->second)));
	    
		sending.back()->isend(*mit, SUM_TAG);
	      }
	      raw_pool.erase(rit);
	    }
	    std_pool.erase(std_graph);
	  }
	  // graph not in the database
	  else {
	    ++calc_count;

	    waiting.push_back(std::make_pair(std_graph, SUM_TAG));
	
	    // node request
	    MPI::COMM_WORLD.Isend(0, 0, MPI::INT, MASTER, NODE_TAG).Free();
	  }
	}
      }
      //
      // sum result
      //
      else if(tag == SUM_TAG) {
	if(!_is_work_node(from)) {
	  ErrOut err_out;
	  err_out << funame << "fourier sum result should come from working node: " << from;
	}

	_gin_t gin(from, tag);

	const _fg_t& graph = gin.first;

	std_t::iterator sgit = std_pool.find(graph);
      
	if(sgit == std_pool.end()) {
	  ErrOut err_out;
	  err_out << funame << "standard graph is not in the pool";
	}
	
	_Convert::vec_t gconv = _convert(graph);
	
	if(sum_data.find(gconv) != sum_data.end()) {
	  ErrOut err_out;
	  err_out << funame << "graph value is already in the database";
	}

	sum_data[gconv] = gin.second;

	const std::set<_fg_t>& rg = sgit->second;
	
	for(std::set<_fg_t>::const_iterator git = rg.begin(); git != rg.end(); ++git) {
	
	  raw_t::iterator rit = raw_pool.find(*git);
	
	  if(rit == raw_pool.end()) {
	    ErrOut err_out;
	    err_out << funame << "raw graph is not in the pool";
	  }
	    
	  for(std::multiset<int>::const_iterator mit = rit->second.begin(); 
	      mit != rit->second.end(); ++mit) {
	    sending.push_back(SharedPointer<_gbase_t>(new _gin_t(*git, gin.second)));
	    
	    sending.back()->isend(*mit, SUM_TAG);
	  }
	  raw_pool.erase(rit);
	}
	std_pool.erase(sgit);
      }
      //
      // wrong tag
      //
      else {
	ErrOut err_out;
	err_out << funame << "unknown tag: " << tag;
      }
    }

    // some checking
    if(sending.size()) {
      ErrOut err_out;
      err_out << funame << "some sends not completed";
    }

    if(raw_pool.size()) {
      ErrOut err_out;
      err_out << funame << "raw graph pool is not empty";
    }

    if(std_pool.size()) {
      ErrOut err_out;
      err_out << funame << "standard graph pool is not empty";
    }

    if(waiting.size()) {
      ErrOut err_out;
      err_out << funame << "there are pending jobs";
    }

    if(temperature > 0.) {
      MPI::COMM_WORLD.Send(&calc_count, 1, MPI::LONG, MASTER, STAT_TAG);
      MPI::COMM_WORLD.Send(&read_count, 1, MPI::LONG, MASTER, STAT_TAG);
      MPI::COMM_WORLD.Send(&conv_count, 1, MPI::LONG, MASTER, STAT_TAG);
    }
  }
  catch(Error::General) {
    std::cerr << funame << "Oops\n";
    throw;
  }
}

void GraphExpansion::_master (double temperature, int mode) const
{
  const char funame [] = "GraphExpansion::_master: ";

  try {
    int    itemp;
    double dtemp;

    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  
    if(mpi_rank) {
      ErrOut err_out;
      err_out << funame << "not a master";
    }

    if(temperature > 0.) {
      IO::log << IO::log_offset << "temperature(K) = " << temperature / Phys_const::kelv << "\n";
    }
    // zero temperature
    else {
      for(int i = 0.; i < _red_freq.size(); ++i)
	if(_red_freq[i] <= 0.) {
	  ErrOut err_out;
	  err_out << funame << "for zero temperature calculation all frequencies should be real";
	}

      if(GLOBAL == mode)
	IO::log << IO::log_offset << "zero-point energy calculation:\n";
      else
	IO::log << IO::log_offset << "low temperature expansion calculation:\n";
    }

    std::set<int>    free_nodes;
    for(int i = WORK_NODE; i < mpi_size; ++i)
      free_nodes.insert(i);

    std::list<int>   waiting;

    std::map<int, double> graph_value;

    std::map<int, std::map<int, double> > graphex_value; // low temperature expansion

    int gindex = 0;  // current graph index

    int end_num = 0; // number of stopped drivers

    typedef std::list<SharedPointer<_gbase_t> > sent_t;
    sent_t sending;
  
    //
    // main loop
    //
    while(graph_value.size() < _graph_data.size() && graphex_value.size() < _graph_data.size() || end_num < WORK_NODE - SUM_SERV - 1) {
      MPI::Status stat;
      MPI::COMM_WORLD.Recv(0, 0, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);
    
      const int tag  = stat.Get_tag();
      const int from = stat.Get_source();

      //
      // graph index result
      //
      if(tag == GINDEX_TAG) {
	if(!_is_driver(from)) {
	  ErrOut err_out;
	  err_out << funame << "graph index request should come from a driver: " << from;
	}

	_gval_t gval;
	_gmap_t gmap;

	int gind;
	if(temperature <= 0. && mode == CENTROID) {
	  gmap.recv(from, tag);

	  if(!graphex_value.insert(gmap).second) {
	    ErrOut err_out;
	    err_out << funame << "graph value is already in the map: " << gmap.first;
	  }

	  itemp = graphex_value.size();

	  gind = gmap.first;

	  if(gind < 0 || gind >= _graph_data.size()) {
	    ErrOut err_out;
	    err_out << funame << "graph index out of range: " << gmap.first;
	  }
	}
	else {
	  gval.recv(from, tag);

	  if(!graph_value.insert(gval).second) {
	    ErrOut err_out;
	    err_out << funame << "graph value is already in the map: " << gval.first;
	  }

	  itemp = graph_value.size();

	  gind = gval.first;

	  if(gind < 0 || gind >= _graph_data.size()) {
	    ErrOut err_out;
	    err_out << funame << "graph index out of range: " << gval.first;
	  }
	}

	const int sym_map_size = _gnode_t(from, tag);
      
	IO::log << IO::log_offset 
		<< "gindex = "               << std::setw(6) << gind 
		<< "   from  = "             << std::setw(6) << from
		<< "   symbolic map size = " << std::setw(6) << sym_map_size 
		<< "   received graphs # = " << std::setw(6) << itemp
		<< std::endl;
      }
      //
      // graph index request
      //
      else if(tag == GINDEX_REQUEST) {
	if(!_is_driver(from)) {
	  ErrOut err_out;
	  err_out << funame << "graph index request should come from a driver: " << from;
	}

	sending.push_back(SharedPointer<_gbase_t>(new _gnode_t(gindex)));

	sending.back()->isend(from, GINDEX_TAG, 1);
      
	if(gindex < _graph_data.size())
	  ++gindex;
      }
      //
      // driver stopped working
      //
      else if(tag == END_TAG) {
	if(!_is_driver(from)) {
	  ErrOut err_out;
	  err_out << funame << "end of work sygnal should come from the driver: " << from;
	}

	++end_num;
      }
      //
      // woriking node release request
      //
      else if(tag == NODE_RELEASE) {
	if(!_is_work_node(from)) {
	  ErrOut err_out;
	  err_out << funame << "working node release request should come from a working node: " << from;
	}
	
	if(free_nodes.find(from) != free_nodes.end()) {
	  ErrOut err_out;
	  err_out << funame << "free node " << from << " should not request for release";
	}

	free_nodes.insert(from);
      }
      //
      // work node request
      //
      else if(tag == NODE_TAG) {
	if(!_is_server(from)) {
	  ErrOut err_out;
	  err_out << funame << "working node request should come from a server: " << from;
	}
	
	waiting.push_back(from);
      }
      //
      // wrong tag
      //
      else {
	ErrOut err_out;
	err_out << funame << "unknown tag: " << tag;
      }

      // respond to working node requests
      while(free_nodes.size() && waiting.size()) {
      
	sending.push_back(SharedPointer<_gbase_t>(new _gnode_t(*free_nodes.begin())));

	sending.back()->isend(*waiting.begin(), NODE_TAG);
      
	free_nodes.erase(free_nodes.begin());
	waiting.erase(waiting.begin());
      }

      // clean the completed sends
      for(sent_t::iterator sit = sending.begin(); sit != sending.end(); )
	if((*sit)->Test())
	  sit = sending.erase(sit);
	else
	  ++sit;
      //
      //
    } // main loop

    for(int node = 1; node < mpi_size; ++node)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, node, END_TAG);

    // some checking
    if(waiting.size()) {
      ErrOut err_out;
      err_out << funame << "there are pending requests";
    }

    if(sending.size()) {
      ErrOut err_out;
      err_out << funame << "some sends not completed";
    }

    long     conv_count = 0;

    long int_calc_count = 0;
    long int_read_count = 0;

    long zpe_calc_count = 0;
    long zpe_read_count = 0;

    long sum_calc_count = 0;
    long sum_read_count = 0;

    MPI::COMM_WORLD.Recv(&int_calc_count, 1, MPI::LONG, INT_SERV, STAT_TAG);
    MPI::COMM_WORLD.Recv(&int_read_count, 1, MPI::LONG, INT_SERV, STAT_TAG);
    MPI::COMM_WORLD.Recv(    &conv_count, 1, MPI::LONG, INT_SERV, STAT_TAG);

    if(temperature > 0.) {
      long count;

      MPI::COMM_WORLD.Recv(&zpe_calc_count, 1, MPI::LONG, ZPE_SERV, STAT_TAG);
      MPI::COMM_WORLD.Recv(&zpe_read_count, 1, MPI::LONG, ZPE_SERV, STAT_TAG);
      MPI::COMM_WORLD.Recv(         &count, 1, MPI::LONG, ZPE_SERV, STAT_TAG);
      conv_count += count;

      MPI::COMM_WORLD.Recv(&sum_calc_count, 1, MPI::LONG, SUM_SERV, STAT_TAG);
      MPI::COMM_WORLD.Recv(&sum_read_count, 1, MPI::LONG, SUM_SERV, STAT_TAG);
      MPI::COMM_WORLD.Recv(         &count, 1, MPI::LONG, SUM_SERV, STAT_TAG);
      conv_count += count;
    }
    
    IO::log << "\n";
    IO::log << IO::log_offset
	    << "Statistics:\n";

    IO::log << IO::log_offset;
    if(conv_count)
      IO::log << std::setw(12) << "Conv #";
    if(int_calc_count)
      IO::log << std::setw(12) << "Int Calc #";
    if(int_read_count)
      IO::log << std::setw(12) << "Int Read #";
    if(zpe_calc_count)
      IO::log << std::setw(12) << "ZPE Calc #";
    if(zpe_read_count)
      IO::log << std::setw(12) << "ZPE Read #";
    if(sum_calc_count)
      IO::log << std::setw(12) << "Sum Calc #";
    if(sum_read_count)
      IO::log << std::setw(12) << "Sum Read #";
    IO::log << "\n";

    IO::log << IO::log_offset;
    if(conv_count)
      IO::log << std::setw(12) << conv_count;
    if(int_calc_count)
      IO::log << std::setw(12) << int_calc_count;
    if(int_read_count)
      IO::log << std::setw(12) << int_read_count;
    if(zpe_calc_count)
      IO::log << std::setw(12) << zpe_calc_count;
    if(zpe_read_count)
      IO::log << std::setw(12) << zpe_read_count;
    if(sum_calc_count)
      IO::log << std::setw(12) << sum_calc_count;
    if(sum_read_count)
      IO::log << std::setw(12) << sum_read_count;
    IO::log << "\n\n";

    double std_time = 0.;
    double zpe_time = 0.;
    double sum_time = 0.;

    for(int node = WORK_NODE; node < mpi_size; ++node) {

      MPI::COMM_WORLD.Recv(&dtemp, 1, MPI::DOUBLE, node, STAT_TAG);
      std_time += dtemp;

      MPI::COMM_WORLD.Recv(&dtemp, 1, MPI::DOUBLE, node, STAT_TAG);
      zpe_time += dtemp;

      if(temperature > 0.) {
	MPI::COMM_WORLD.Recv(&dtemp, 1, MPI::DOUBLE, node, STAT_TAG);
	sum_time += dtemp;
      }
    }

    IO::log << IO::log_offset << "Working time[sec]:\n";
    
    IO::log << IO::log_offset
	    << std::setw(15) << "Conv"
	    << std::setw(15) << "ZPE";
    if(temperature > 0.)
      IO::log << std::setw(15) << "Sum";
    IO::log << "\n";

    IO::log << IO::log_offset
	    << std::setw(15) << std_time / CLOCKS_PER_SEC
	    << std::setw(15) << zpe_time / CLOCKS_PER_SEC;
    if(temperature > 0.)
      IO::log << std::setw(15) << sum_time / CLOCKS_PER_SEC;
    IO::log << "\n\n";

    if(temperature > 0.) {
      std::map<int, double> corr;
      for(int i = 0; i < _sorted_graph.size(); ++i)
	corr[_sorted_graph[i].size()] += graph_value[i];

      IO::log << IO::log_offset << "total correction:\n";

      IO::log << IO::log_offset 
	      << std::setw(2) << "BO" 
	      << std::setw(15) << "Value" 
	      << "\n";

      double res = 0.;
      for(std::map<int, double>::const_iterator cit = corr.begin(); cit != corr.end(); ++cit) {

	res += cit->second;

	IO::log << IO::log_offset 
		<< std::setw(2) << cit->first 
		<< std::setw(15) << res 
		<< "\n";
      }
      IO::log << "\n";
    }
    // zero temperature
    else {
      // global mode
      if(GLOBAL == mode) {
	std::map<int, double> zpe;

	for(int i = 0; i < _sorted_graph.size(); ++i)
	  zpe[_sorted_graph[i].size()] += graph_value[i];

	IO::log << IO::log_offset << "zero-point energy correction, 1/cm:\n";
      
	IO::log << IO::log_offset 
		<< std::setw(2) << "BO" 
		<< std::setw(15) << "Value" 
		<< "\n";

	double res = 0.;
	for(std::map<int, double>::const_iterator gzit = zpe.begin(); gzit != zpe.end(); ++gzit) {

	  res += gzit->second;

	  IO::log << IO::log_offset 
		  << std::setw(2) << gzit->first 
		  << std::setw(15) << -res / Phys_const::incm
		  << "\n";
	}
	IO::log << "\n";
      }
      // centroid mode
      else {
	int print_term_max = 3;

	std::map<int, std::map<int, double> > graphex;

	for(int i = 0; i < _sorted_graph.size(); ++i)
	  for(std::map<int, double>::const_iterator gxit = graphex_value[i].begin(); gxit != graphex_value[i].end(); ++gxit)
	    graphex[_sorted_graph[i].size()][gxit->first] += gxit->second;

	IO::log << IO::log_offset << "low temperature expansion:\n";

	IO::log << IO::log_offset 
		<< std::setw(2)  << "BO"
		<< std::setw(15) << "1/T(K) term"
		<< std::setw(15) << "T^0 term";

	for(int i = 1; i <= print_term_max; ++i)
	  IO::log << std::setw(9) << "T(K)^" << i << " term";
	IO::log << "\n";

	std::map<int, double> res;
	for(std::map<int, std::map<int, double> >::const_iterator gxit = graphex.begin(); gxit != graphex.end(); ++gxit) {
	
	  for(std::map<int, double>::const_iterator it = gxit->second.begin(); it != gxit->second.end(); ++it)
	    res[it->first] += it->second;

	  itemp = -res.begin()->first;
	  if(itemp > 1) {
	    ErrOut err_out;
	    err_out << funame << "low temperature expansion begins with 1/T^" << itemp << "terms";
	  }

	  IO::log << IO::log_offset << std::setw(2) << gxit->first;
	
	  for(int i = -1 ; i <= print_term_max; ++i) {
	    IO::log << std::setw(15);

	    if(res.find(i) == res.end()) {
	      IO::log << "0";
	    }
	    else {
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

  _vertex_size_max = v;
  _freq_size       = f;

  _index_map.resize(v * (v - 1) / 2);

  itemp = 0;
  for(int i = 1; i < v; ++i)
    for(int j = 0; j < i; ++j, ++itemp) {
      _index_map[itemp].insert(i);
      _index_map[itemp].insert(j);
    }
}

GraphExpansion::_fg_t GraphExpansion::_Convert::operator() (const vec_t& gconv) const
{
  const char funame [] = "GraphExpansion::_Convert::operator(): ";
  
  int itemp;

  if(!_vertex_size_max) {
    ErrOut err_out;
    err_out << funame << "frequency graph data converter not initialized";
  }

  _fg_t res;

  if(!gconv.size())
    return res;

  for(int i = 0; i < gconv.size(); ++i) {
    if(gconv[i] < 0) {
      ErrOut err_out;
      err_out << funame << "negative index";
    }

    itemp = gconv[i] / _freq_size;
    
    if(itemp >= _vertex_size_max * (_vertex_size_max - 1) / 2) {
      ErrOut err_out;
      err_out << funame << "vertex index convert out of range";
    }
    
    res[_index_map[itemp]].insert(gconv[i] % _freq_size);
  }

  return res;
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
    itemp  = *git->first.begin() + itemp * (itemp - 1) / 2;
    itemp *= _freq_size;
    
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit, ++bond_index) {
      if(*fit < 0 || *fit >= _freq_size) {
	ErrOut err_out;
	err_out << funame << "frequency index out of range: " << *fit;
      }

      res[bond_index] = itemp + *fit;
    }
  }
  return res;  
}

long GraphExpansion::_db_t::mem_size () const
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

void GraphExpansion::_fg_t::print (std::ostream& to) const
{
  for(const_iterator git = begin(); git != end(); ++git) {
    if(git != begin())
      to << "...";
    to << "(";
    for(std::set<int>::const_iterator vit = git->first.begin(); vit != git->first.end(); ++vit) {
      if(vit != git->first.begin())
	to << ", ";
      to << *vit;
    }
    to << ")->(";

    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      if(fit != git->second.begin())
	to << ", ";
      to << *fit;
    }
    to << ")";
  }
}

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
    if(git->first.size() != 2) {
      ErrOut err_out;
      err_out << funame << "bond should have exactly two different ends:";
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	err_out << " " << *bit;
    }
    
    if(!git->second.size()) {
      ErrOut err_out;
      err_out << funame << "zero-order bond";
    }

    if(freq_size >= 0)
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit)
	if(*fit >= freq_size || *fit < 0) {
	  ErrOut err_out;
	  err_out << funame << "frequency index out of order: " << *fit;
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

  int itemp;

  _mg_t res;
  if(!size())
    return res;

  _check_integrity();

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
						      _mg_t&                     zpe_graph,
						      const std::set<int>&       low_freq
						      ) const
{
  const char funame [] = "GraphExpansion::_fg_t::reduce: ";

  _check_integrity();
  _check_order();

  int  itemp;
  double dtemp;

  const int vsize = vertex_size();

  // reduced graph
  _fg_t red_graph;

  std::vector<std::set<int> > red_group(vsize);
  for(int v = 0; v < vsize; ++v)
    red_group[v].insert(v);

  while(1) {// reduction loop

    // new effective graph
    red_graph.clear();
    for(_fg_t::const_iterator git = begin(); git != end(); ++git) {

      std::set<int> bond;
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	for(int v = 0; v < red_group.size(); ++v)
	  if(red_group[v].find(*bit) != red_group[v].end()) {
	    bond.insert(v);
	    break;
	  }

      if(bond.size() == 2)
	for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit)
	  red_graph[bond].insert(*fit);
    }

    // which groups to merge
    std::vector<std::set<int> > merge_group;
    for(_fg_t::const_iterator git = red_graph.begin(); git != red_graph.end(); ++git) {
      dtemp = 0.;
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit)
	if(low_freq.find(*fit) == low_freq.end())
	  dtemp += freq[*fit] * tanh_factor[*fit];
      
      if(dtemp > temperature * red_thresh) {
	std::multiset<int> merge;
	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit) {
	  for(int v = 0; v < merge_group.size(); ++v)
	    if(merge_group[v].find(*bit) != merge_group[v].end()) {
	      merge.insert(v);
	      break;
	    }
	}

	switch(merge.size()) {
	case 0:
	  // create new group
	  merge_group.push_back(git->first);
	  break;

	case 1:
	  // add free vertex to the group
	  for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	    merge_group[*merge.begin()].insert(*bit);
	  break;

	case 2:
	  // vertices already belong to the same group: do nothing
	  if(*merge.begin() == *merge.rbegin())
	    break;

	  // merge two groups
	  for(std::set<int>::const_iterator mit = merge_group[*merge.begin()].begin(); mit != merge_group[*merge.begin()].end(); ++mit)
	    merge_group[*merge.rbegin()].insert(*mit);
	  merge_group.erase(merge_group.begin() + *merge.begin());
	  break;
	}
      }
    }
	    
    if(!merge_group.size())
      break;
	    
    std::set<int> pool;
    std::vector<std::set<int> > new_red_group(merge_group.size());
    for(int v = 0; v < merge_group.size(); ++v) {

      for(std::set<int>::const_iterator mit = merge_group[v].begin(); mit != merge_group[v].end(); ++mit) {
	pool.insert(*mit);

	for(std::set<int>::const_iterator vit = red_group[*mit].begin(); vit != red_group[*mit].end(); ++vit)
	  new_red_group[v].insert(*vit);
      }
    }
    for(int v = 0; v < red_group.size(); ++v)
      if(pool.find(v) == pool.end())
	new_red_group.push_back(red_group[v]);

    red_group = new_red_group;
  } // reduction loop

  // zpe graph
  for(int v = 0; v < red_group.size(); ++v) {// zpe factor calculation cycle
    if(red_group[v].size() == 1)
      continue;

    std::map<int, int> vertex_map;
    itemp = 0;
    for(std::set<int>::const_iterator it = red_group[v].begin(); it != red_group[v].end(); ++it, ++itemp)
      vertex_map[*it] = itemp;
	    
    _fg_t new_graph;
    for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
      std::set<int> bond;
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit) {
	if(red_group[v].find(*bit) != red_group[v].end())
	  bond.insert(vertex_map[*bit]);

	if(bond.size() == 2)
	  new_graph[bond] = git->second;
      }
    }
    ++zpe_graph[new_graph];
  }

  return red_graph;
}

double GraphExpansion::_fg_t::zpe_factor (const std::vector<double>& freq, 
					  double                     temperature, 
					  const std::vector<double>& tanh_factor,
					  const std::set<int>&       low_freq
					  ) const
{
 const char funame [] = "GraphExpansion::_fg_t::zpe_factor: ";

 _check_integrity(freq.size());
 _check_order();

  int    itemp;
  double dtemp;

  if(!size())
    return 1.;

  const int vsize = vertex_size();

  int symm_fac = 0;
  std::set<_fg_t> perm_set = perm_pool(&symm_fac);

  double res = 0.;

  //  if(perm_graph.size() < max_thread_num)
  // omp_set_num_threads(perm_graph.size());

  //  #pragma omp parallel for default(shared) private(itemp, dtemp) reduction(+: res) schedule(dynamic)

  for(std::set<_fg_t>::const_iterator pit = perm_set.begin(); pit != perm_set.end(); ++pit) {

    std::vector<double> freq_band(vsize - 1);
    for(_fg_t::const_iterator git = pit->begin(); git != pit->end(); ++git) {

      dtemp = 0.;
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
	if(low_freq.find(*fit) != low_freq.end()) {
	  dtemp += 6. * temperature;
	}
	else if(temperature > 0.) {
	  dtemp += freq[*fit] * tanh_factor[*fit];
	}
	else {
	  dtemp += freq[*fit];
	}
      }

      for(int i = *git->first.begin(); i < *git->first.rbegin(); ++i)
	freq_band[i] += dtemp;
    }

    if(temperature > 0.) {
      itemp = 0;
      for(int i = 0; i < freq_band.size(); ++i) {
	dtemp = freq_band[i] / temperature;
	if(dtemp < red_thresh) {
	  itemp = 1;
	  break;
	}
      }
      if(itemp)
	std::cout << funame << "WARNING: high frequency cutoff condition not met: " << dtemp << "\n";
    }
    
    dtemp = 1.;
    for(int i = 0; i < freq_band.size(); ++i)
      dtemp /= freq_band[i];
    
    res += dtemp;
  }

  //  if(perm_graph.size() < max_thread_num)
  //  omp_set_num_threads(max_thread_num);

  res *= (double)symm_fac;

  for(_fg_t::const_iterator git = begin(); git != end(); ++git) {
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      if(low_freq.find(*fit) != low_freq.end()) {
	res /= 12. * temperature;
      }
      else if(temperature > 0.) {
	res /= tanh_factor[*fit] * 2. * freq[*fit];
      }
      else {
	res /= 2. * freq[*fit];
      }
    }
  }

  return res;
}

double GraphExpansion::_fg_t::fourier_sum (const std::vector<double>& freq, 
					   double                     temperature,
					   const std::set<int>&       low_freq
					   ) const
{
  const char funame [] = "GraphExpansion::_fg_t::fourier_sum: ";
  
  int    itemp;
  double dtemp;

  if(!size())
    return 1. / temperature;

  _check_integrity(freq.size());
  _check_order();

  const int vsize = vertex_size();

  std::vector<int> freq_map;
  for(_fg_t::const_iterator git = begin(); git != end(); ++git)
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit)
      freq_map.push_back(*fit);
  
  std::vector<int> index_shift(freq_map.size());
  std::vector<int> index_range(freq_map.size());
  std::vector<double>     xval(freq_map.size());

  for(int i = 0; i < freq_map.size(); ++i) {
    itemp = freq_map[i];
    dtemp = freq[itemp] / temperature / 2. / M_PI;

    if(dtemp <= -1.) {
      ErrOut err_out;
      err_out <<funame << "deep tunneling regime: "
	      << "frequency[1/cm] = "	<< freq[freq_map[i]] / Phys_const::incm << ",  "
	      << "Temperature[K]  = " << temperature / Phys_const::kelv;
    }
    
    if(dtemp < 0.) {
      xval[i] = - dtemp * dtemp;
      dtemp = -dtemp;
    }
    else {
      xval[i] = dtemp * dtemp;
    }
    
    if(dtemp > 1.)
      itemp = (int)std::ceil(four_cut * dtemp);
    else
      itemp = (int)std::ceil(four_cut);

    index_shift[i] = itemp;
    index_range[i] = 2 * itemp + 1;
  }   
    
  double res = 0.;

  MultiIndexConvert harmonic_index(index_range);

  for(long ml = 0; ml < harmonic_index.size(); ++ml) {

    std::vector<int> mi = harmonic_index(ml);

    for(int i = 0; i < mi.size(); ++i)
      mi[i] -= index_shift[i];
    
    std::vector<int> constrain(vsize);

    itemp = 0;
    for(_fg_t::const_iterator git = begin(); git != end(); ++git)
      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit, ++itemp) {
	int sign = -1;
	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit, sign += 2)
	  constrain[*bit] += sign * mi[itemp];
      }
    
    int test = 0;
    for(std::vector<int>::const_iterator it = constrain.begin(); it != constrain.end(); ++it)
      if(*it) {
	test = 1;
	break;
      }
    
    if(test)
      continue;

    dtemp = 1.;
    for(int i = 0; i < freq_map.size(); ++i) {
      itemp = mi[i];

      if(itemp || low_freq.find(freq_map[i]) == low_freq.end()) {
	dtemp /= xval[i] + double(itemp * itemp);
      }
      else {
	test = 1;
	break;
      }
    }

    if(test)
      continue;

    res += dtemp;
  }

  //if(harmonic_index.size() < max_thread_num)
  //  omp_set_num_threads(max_thread_num);

  res /= std::pow(temperature, (double)vsize) * std::pow(4. * M_PI * M_PI * temperature, freq_map.size());

  return res;
}
