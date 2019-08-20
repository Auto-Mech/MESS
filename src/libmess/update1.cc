// update graph database for certain perturbation terms
void GraphExpansion::update (const PotentialExpansion& pex)
{
  const char funame [] = "GraphExpansion::update: ";

  IO::Marker funame_marker(funame);

  int itemp;

  _used.insert(std::map<int, int>());
  
  if(_potex.find(pex.rank()) != _potex.end()) {
    std::cerr << funame << "graph expansion of " << pex.rank() << " order already inserted\n";
    throw Error::Logic();
  }

  _potex[pex.rank()] = pex;

  int unconnected = 0;
  int graph_count = 0;
  int   raw_count = 0;
  
  std::set<std::map<int, int> > new_used = _used;
  int pex_count = 0;
  
  while(1) {
    ++pex_count;
    if(pex_count * pex.rank() > 2 * bond_max)
      break;
    for(std::set<std::map<int, int> >::const_iterator pexit = _used.begin(); pexit != _used.end(); ++pexit) {
      std::map<int, int> ginit = *pexit;
      ginit[pex.rank()] = pex_count;
      
      itemp = pex_count * pex.rank();
      for(std::map<int, int>::const_iterator it = pexit->begin(); it != pexit->end(); ++it)
	itemp += it->first * it->second;
    
      if(itemp % 2 || itemp > 2 * bond_max)
	continue;
      
      new_used.insert(ginit);

      for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
	if(it->first < 3 || it->second < 1) {
	  std::cerr << funame << "potential expansion initializer out of range: " << it->first << ", " << it->second << "\n";
	  throw Error::Range();
	}

      IO::log << IO::log_offset << "perturbation term:";
      for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
	IO::log << "  " << it->first << "(" << it->second << ")";
      IO::log << std::endl;
	     
      // odd number of connection points
      itemp = 0;
      for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
	itemp += it->first * it->second;
      if(!itemp || itemp % 2)
	return;

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
  
      std::set<_graph_t> raw_graph = _raw_graph_generator(vertex_order);

      raw_count += raw_graph.size();
  
      itemp = 1;
      for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
	itemp *= Math::factorial(it->second);
      const int vertex_perm_size = itemp;

      while(1) {
	// erase not-connected graphs
	while(raw_graph.size() && !raw_graph.begin()->is_connected()) {
	  raw_graph.erase(raw_graph.begin());
	  ++unconnected;
	}

	if(!raw_graph.size())
	  break;
    
	++graph_count;

	_graph_t ancor = *raw_graph.begin();
    
	if(!_graph.insert(ancor).second) {
	  std::cerr << funame << "graph already in the pool\n";
	  throw Error::Logic();
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

	  std::set<_graph_t>::iterator git = raw_graph.find(perm_graph);
	  if(git != raw_graph.end()) {
	    ++perm_diff;
	    raw_graph.erase(git);
	  }
	}
    
	if(perm_diff * perm_equal != vertex_perm_size) {
	  std::cerr << funame << "numbers of permutationally equal and different graphs inconsistent: "
		    << perm_equal << ", "<< perm_diff << "\n";
	  throw Error::Logic();
	}

	if(perm_equal != ancor.vertex_symmetry()) {
	  std::cerr << funame << "permutational symmetry factors differ\n";
	  throw Error::Logic();
	}
      }
    }
  }

  _used = new_used;

  IO::log << IO::log_offset << "number of raw graphs = " << raw_count << "\n";
  IO::log << IO::log_offset << "number of unconnected graphs = " << unconnected << "\n";
  IO::log << IO::log_offset << "number of permutationally distinct graphs = " << graph_count << "\n";
}
