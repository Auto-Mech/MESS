#include "graph_common.hh"
#include "key.hh"
#include "permutation.hh"
#include "multindex.hh"
#include "io.hh"
#include "units.hh"

#include <list>
#include <cmath>
#include <complex>

namespace Graph {

  // fourier sum graph evaluation cutoff
  //
  double FreqGraph::four_cut  = 10.;

  int    FreqGraph::four_max  = 10;
  
  // low temperature / high frequency graph reduction threshold
  //
  double FreqGraph::red_thresh = 10.;

  // maximal potential expansion order  
  //
  int potex_max  =  4;

  // maximal number of bonds in the graph
  //
  int bond_max   =  6;

  // graph sorted according to certain criteria
  //
  std::vector<GenGraph> _sorted_graph;

  int size () { return _sorted_graph.size(); }

  bool isinit () { return size(); }

  void _check_init ();

  const_iterator begin () { _check_init(); return _sorted_graph.begin(); }  

  std::set<GenGraph> _raw_graph_generator (std::vector<int> vertex_order, int root = -1);
}

void Graph::_check_init () 
{
  const char funame [] = "Graph::_check_init(): ";
      
  if(isinit())
    //
    return;
      
  ErrOut err_out;

  err_out << funame << "not initialized";
}

// perturbation graph theory initializer
//
void Graph::init ()
{
  const char funame [] = "Graph::init: ";

  IO::Marker funame_marker(funame);

  int itemp;

  if(potex_max < 3) {
    //
    ErrOut err_out;

    err_out << "no anharmonic terms in the potential";
  }

  std::vector<int> potex_rank(potex_max - 2);
  //
  for(int i = 0; i < potex_rank.size(); ++i)
    //
    potex_rank[i] = i + 3;

  std::vector<int> glimit(potex_rank.size());
  //
  for(int i = 0; i < potex_rank.size(); ++i)
    //
    glimit[i] = 2 * bond_max / potex_rank[i];


  int unconnected = 0;
  //
  int graph_count = 0;
  //
  int   raw_count = 0;

  std::set<GenGraph> graph_set;

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

    if(!IO::mpi_rank) {
      //
      IO::log << IO::log_offset << "perturbation term:";

      for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
	//
	IO::log << "  " << it->first << "(" << it->second << ")";

      IO::log << std::endl;
    }
   
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

      for(int i = 0; i < it->second; ++i) {
	//
	vertex_perm_base.back()[i] = vertex_order.size();

	vertex_order.push_back(it->first);
      }
    }
  
    std::set<GenGraph> raw_graph = _raw_graph_generator(vertex_order);

    raw_count += raw_graph.size();
  
    itemp = 1;
    //
    for(std::map<int, int>::const_iterator it = ginit.begin(); it != ginit.end(); ++it)
      //
      itemp *= Permutation::factorial(it->second);

    const int vertex_perm_size = itemp;

    while(1) {
      //
      // erase not-connected graphs
      //
      while(raw_graph.size() && !raw_graph.begin()->is_connected()) {
	//
	raw_graph.erase(raw_graph.begin());

	++unconnected;
      }

      if(!raw_graph.size())
	//
	break;
    
      ++graph_count;

      GenGraph ancor = *raw_graph.begin();
    
      if(!graph_set.insert(ancor).second) {
	//
	std::cerr << funame << "graph already in the pool\n";

	throw Error::Logic();
      }

      int perm_equal = 0;

      int perm_diff  = 0;

      // remove all permutationally related graphs from the raw graph pool
      //
      for(MultiPerm multi_perm(multi_perm_init); !multi_perm.end(); ++multi_perm) {
	//
	std::vector<int> vertex_perm(vertex_order.size());

	for(int i = 0; i < vertex_perm_base.size(); ++i) {
	  //
	  for(int j = 0; j < vertex_perm_base[i].size(); ++j) {
	    //
	    vertex_perm[vertex_perm_base[i][j]] = vertex_perm_base[i][multi_perm[i][j]];
	  }
	}
      
	GenGraph perm_graph;
	//
	for(GenGraph::const_iterator pit = ancor.begin(); pit != ancor.end(); ++pit) {
	  //
	  std::multiset<int> bond;

	  for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it)
	    //
	    bond.insert(vertex_perm[*it]);
	  
	  perm_graph.insert(bond);
	}

	if(perm_graph == ancor)
	  //
	  ++perm_equal;

	std::set<GenGraph>::iterator git = raw_graph.find(perm_graph);
	//
	if(git != raw_graph.end()) {
	  //
	  ++perm_diff;

	  raw_graph.erase(git);
	}
      }
    
      if(perm_diff * perm_equal != vertex_perm_size) {
	//
	ErrOut err_out;

	err_out << funame << "numbers of permutationally equal and different graphs inconsistent: "
		<< perm_equal << ", "<< perm_diff;
      }

      if(perm_equal != ancor.vertex_symmetry()) {
	//
	ErrOut err_out;

	err_out << funame << "permutational symmetry factors differ";
      }
    }
  }

  // sorting graph according to certain criterion (vertex size, for example)
  //
  std::multimap<int, GenGraph> sort_map;
  //
  for(std::set<GenGraph>::const_iterator gsit = graph_set.begin(); gsit != graph_set.end(); ++gsit)
    //
    sort_map.insert(make_pair(gsit->vertex_size(), *gsit));

  _sorted_graph.resize(sort_map.size());
    
  itemp = 0;
  //
  for(std::multimap<int, GenGraph>::const_reverse_iterator sit = sort_map.rbegin(); sit != sort_map.rend(); ++sit, ++itemp)
    //
    _sorted_graph[itemp] = sit->second;

  if(!IO::mpi_rank) {
    //
    IO::log << "\n";
  
    IO::log << IO::log_offset << "number of raw graphs = " << raw_count << "\n";

    IO::log << IO::log_offset << "number of unconnected graphs = " << unconnected << "\n";

    IO::log << IO::log_offset << "number of permutationally distinct graphs = " << graph_count << "\n\n";

    IO::log << IO::log_offset << "permutationally distinct graphs(" << _sorted_graph.size() << "):\n";

    IO::log << IO::log_offset
	    << std::setw(5) << "#"
	    << std::setw(5) << "VS"
	    << std::setw(5) << "BS"
	    << std::setw(5) << "V#"
	    << std::setw(5) << "B#"
	    << std::setw(3) << " "
	    << "graph" << "\n";

    itemp = 0;
    for(std::vector<GenGraph>::const_iterator it = _sorted_graph.begin(); it != _sorted_graph.end(); ++it, ++itemp)
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

/*************************************************************************************
 ************************ GENERIC PERTURBATION THEORY GRAPH **************************
 *************************************************************************************/

std::vector<std::multiset<int> >  Graph::GenGraph::vertex_bond_map () const
{
  const char funame [] = "Graph::GenGraph::vertex_bond_map: ";
  
  _check();

  std::vector<std::multiset<int> > res(vertex_size());
  
  int bi = 0;
  //
  for(const_iterator git = begin(); git != end(); ++git, ++bi)
    //
    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
      //
      res[*bit].insert(bi);

  return res;
}

bool Graph::GenGraph::is_connected () const
{
  const char funame [] = "Graph::GenGraph::is_connected: ";

  if(!size())
    //
    return true;
  
  _check();

  std::set<int> pool, new_pool;
  //
  new_pool.insert(*begin()->begin());

  while(pool.size() != new_pool.size()) {
    //
    pool = new_pool;

    for(const_iterator git = begin(); git != end(); ++git) {
      //
      for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit) {
	//
	if(new_pool.find(*bit) != new_pool.end()) {
	  //
	  if(bit == git->begin()) {
	    //
	    new_pool.insert(*git->rbegin());
	  }
	  else {
	    //
	    new_pool.insert(*git->begin());
	  }

	  break;
	}
      }
    }
  }

  if(pool.size() == vertex_size())
    //
    return true;

  return false;
}

int Graph::GenGraph::vertex_symmetry () const
{
  const char funame [] = "Graph::GenGraph::vertex_symmetry: ";

  int itemp;

  if(!size())
    //
    return 1;

  _check();
  
  std::vector<int> vertex_order(vertex_size());

  for(const_iterator git = begin(); git != end(); ++git)
    //
    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
      //
      ++vertex_order[*bit];

  std::map<int, std::vector<int> > order_map;
  //
  for(int i = 0; i < vertex_order.size(); ++i)
    //
    order_map[vertex_order[i]].push_back(i);

  if(order_map.begin()->first < 3) {
    //
    ErrOut err_out;

    err_out << funame << "vertex order less then 3: " << order_map.begin()->first;
  }

  itemp = 0;
  //
  std::vector<int> order_size(order_map.size());
  //
  for(std::map<int, std::vector<int> >::const_iterator it = order_map.begin(); it != order_map.end(); ++it, ++itemp)
    //
    order_size[itemp] = it->second.size();

  int res = 0;

  // vertex permutation symmetry factor
  //
  for(MultiPerm multi_perm(order_size); !multi_perm.end(); ++multi_perm) {
    //
    itemp = 0;

    std::vector<int> vertex_perm(vertex_order.size());
    //
    for(std::map<int, std::vector<int> >::const_iterator it = order_map.begin(); it != order_map.end(); ++it, ++itemp)
      //
      for(int i = 0; i < it->second.size(); ++i)
	//
	vertex_perm[it->second[i]] = it->second[multi_perm[itemp][i]];

    GenGraph new_graph;

    for(const_iterator git = begin(); git != end(); ++git) {
      //
      std::multiset<int> bond;

      for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
	//
	bond.insert(vertex_perm[*bit]);

      new_graph.insert(bond);
    }
    
    if(new_graph == *this)
      //
      ++res;
  }

  return res;
}

int Graph::GenGraph::bond_symmetry () const
{
  const char funame [] = "Graph::GenGraph::bond_symmetry: ";
  
  int itemp;

  if(!size())
    //
    return 1;

  int res = 1;
  
  // loop symmetry
  //
  res <<= loop_size();

  // high order bond symmetry
  //
  std::vector<int> bov = bond_order();
  //
  for(int i = 0; i < bov.size(); ++i)
    //
    res *= Permutation::factorial(bov[i]);
  
  return res;
}

std::vector<int> Graph::GenGraph::bond_order () const
{
  int itemp;

  std::map<std::multiset<int>, int> bond_order_map;
  //
  for(const_iterator git = begin(); git != end(); ++git)
    //
    ++bond_order_map[*git];

  std::vector<int> res(bond_order_map.size());

  itemp = 0;
  //
  for(std::map<std::multiset<int>, int>::const_iterator it = bond_order_map.begin(); it != bond_order_map.end(); ++it, ++itemp)
    //
    res[itemp] = it->second;

  return res;
}

void Graph::GenGraph::_check () const
{
  const char funame [] = "Graph::GenGraph::_check: ";

  if(!size())
    //
    return;

  std::set<int> vertex_pool;
  //
   for(const_iterator git = begin(); git != end(); ++git) {
     //
    if(git->size() != 2) {
      //
      ErrOut err_out;

      err_out << funame << "there should be exactly two points for a bond";
    }

    for(std::multiset<int>::const_iterator bit = git->begin(); bit != git->end(); ++bit)
      //
      vertex_pool.insert(*bit);
  }

  if(*vertex_pool.begin() || *vertex_pool.rbegin() >= vertex_pool.size()) {
    //
    ErrOut err_out;

    err_out << funame << "vertex indices not ordered";
  }
}

int Graph::GenGraph::vertex_size () const
{
  const char funame [] = "Graph::GenGraph::vertex_size: ";

  if(!size())
    //
    return 0;

  std::set<int> res;
  //
  for(const_iterator mit = begin(); mit != end(); ++mit) {
    //
    if(mit->size() != 2) {
      //
      ErrOut err_out;

      err_out << funame << "there should be exactly two points for a bond";
    }

    for(std::multiset<int>::const_iterator it = mit->begin(); it != mit->end(); ++it)
      //
      res.insert(*it);
  }

  if(*res.begin() || *res.rbegin() >= res.size()) {
    //
    ErrOut err_out;

    err_out << funame << "vertex indices not ordered";
  }

  return res.size();
}

int Graph::GenGraph::bond_size () const
{
  const char funame [] = "Graph::GenGraph::bond_size: ";

  int res = 0;
  //
  for(const_iterator git = begin(); git != end(); ++git) {
    //
    if(git->size() != 2) {
      //
      ErrOut err_out;

      err_out << funame << "there should be exactly two points for a bond";
    }

    if(*git->begin() != *git->rbegin())
      //
      ++res;
  }

  return res;
}

int Graph::GenGraph::loop_size () const
{
  const char funame [] = "Graph::GenGraph::bond_size: ";

  int res = 0;
  //
  for(const_iterator git = begin(); git != end(); ++git) {
    //
    if(git->size() != 2) {
      //
      ErrOut err_out;

      err_out << funame << "there should be exactly two points for a bond";
    }

    if(*git->begin() == *git->rbegin())
      //
      ++res;
  }

  return res;
}

void Graph::GenGraph::print (std::ostream& to) const
{
  for(const_iterator pit = begin(); pit != end(); ++pit) {
    //
    for(std::multiset<int>::const_iterator it = pit->begin(); it != pit->end(); ++it) {
      //
      if(it == pit->begin())
	//
	to << "(";
      else
	//
	to << ", ";

      to << *it;
    }

    to << ")";
  }
}

/*************************************************************************************************
 ****************************** FREQUENCY ADAPTED GRAPH ******************************************
 *************************************************************************************************/

void Graph::FreqGraph::print (std::ostream& to) const
{
  for(const_iterator git = begin(); git != end(); ++git) {
    //
    if(git != begin())
      //
      to << "...";

    to << "(";
    //
    for(std::set<int>::const_iterator vit = git->first.begin(); vit != git->first.end(); ++vit) {
      //
      if(vit != git->first.begin())
	//
	to << ", ";
      
      to << *vit;
    }

    to << ")->(";

    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      //
      if(fit != git->second.begin())
	//
	to << ", ";

      to << *fit;
    }

    to << ")";
  }
}

int Graph::FreqGraph::bond_size () const
{
  _check_integrity();
  
  int res = 0;
  //
  for(const_iterator git = begin(); git != end(); ++git)
    //
    res += git->second.size();

  return res;
}

// permutationally equivalent graphs
//
std::set<Graph::FreqGraph> Graph::FreqGraph::perm_pool (int* symm, int flag) const
{
  const char funame [] = "Graph::FreqGraph::perm_pool: ";

  _check_integrity();

  int itemp;

  std::set<FreqGraph> res;

  if(symm)
    //
    *symm = 0;
  
  if(!size())
    //
    return res;
  
  std::set<int> vertex_pool;
  //
  for(const_iterator git = begin(); git != end(); ++git)
    //
    for(std::set<int>::const_iterator it = git->first.begin(); it != git->first.end(); ++it)
      //
      vertex_pool.insert(*it);

  // arbitrary indices
  //
  if(flag) {
    //
    std::vector<int> vertex;

    for(std::set<int>::const_iterator it = vertex_pool.begin(); it != vertex_pool.end(); ++it)
      //
      vertex.push_back(*it);
    
    for(Permutation perm(vertex_pool.size()); !perm.end(); ++perm) {
      //
      FreqGraph perm_graph;

      std::vector<int> perm_vertex = perm(vertex);
      
      std::map<int, int> vertex_map;
      //
      for(int i = 0; i < vertex.size(); ++i)
	//
	vertex_map[vertex[i]] = perm_vertex[i];
      
      for(const_iterator git = begin(); git != end(); ++git) {
	//
	std::set<int> bond;

	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	  //
	  bond.insert(vertex_map[*bit]);
	      
	perm_graph[bond] = git->second;
      }
      
      if(symm && perm_graph == *this)
	//
	++(*symm);
      
      res.insert(perm_graph);
    }
  }
  // ordered vertices
  //
  else {
    //
    if(*vertex_pool.begin() || *vertex_pool.rbegin() >= vertex_pool.size()) {
      //
      ErrOut err_out;

      err_out << funame << "graph not in the standard form";
    }
    
    for(Permutation perm(vertex_pool.size()); !perm.end(); ++perm) {
      //
      FreqGraph perm_graph;

      for(const_iterator git = begin(); git != end(); ++git) {
	//
	std::set<int> bond;

	for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	  //
	  bond.insert(perm[*bit]);
	      
	perm_graph[bond] = git->second;
      }

      if(symm && perm_graph == *this)
	//
	++(*symm);
      
      res.insert(perm_graph);
    }
  }
  
  return res;
}

void Graph::FreqGraph::_check_order () const
{
  const char funame [] = "Graph::FreqGraph::_check_order: ";

  if(!size())
    //
    return;
  
  std::set<int> vertex_pool;
  //
  for(const_iterator git = begin(); git != end(); ++git)
    //
    for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
      //
      vertex_pool.insert(*bit);

  if(!vertex_pool.size() || *vertex_pool.begin() || *vertex_pool.rbegin() >= vertex_pool.size()) {
    //
    ErrOut err_out;

    err_out << funame << "graph not in a standard form";
  }
}

void Graph::FreqGraph::_check_integrity (int freq_size) const
{
  const char funame [] = "Graph::FreqGraph::_check_integrity: ";

  for(const_iterator git = begin(); git != end(); ++git) {
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

    if(freq_size) {
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
 
int Graph::FreqGraph::vertex_size () const
{
  const char funame [] = "Graph::FreqGraph::vertex_size: ";

  _check_integrity();

  std::set<int> res;
  //
  for(const_iterator git = begin(); git != end(); ++git) 
    //
    for(std::set<int>::const_iterator vit = git->first.begin(); vit != git->first.end(); ++vit)
      //
      res.insert(*vit);

  return res.size();
}

// factorize the general graph into a set of connected ones (there may be identical ones) 
//
Graph::FreqGraph::_mg_t  Graph::FreqGraph::factorize () const
{
  const char funame [] = "Graph::FreqGraph::factorize: ";

  _check_integrity();

  int itemp;

  std::vector<FreqGraph> fac;
  //
  for(const_iterator git = begin(); git != end(); ++git) {
    //
    std::set<int> connect;
    //
    for(std::set<int>::const_iterator vit = git->first.begin(); vit != git->first.end(); ++vit) {
      //
      for(int i = 0; i < fac.size(); ++i) {
	//
	itemp = 0;
	//
	for(const_iterator mit = fac[i].begin(); mit != fac[i].end(); ++mit) {
	  //
	  if(mit->first.find(*vit) != mit->first.end()) {
	    //
	    connect.insert(i);

	    itemp = 1;

	    break;
	  }
	}

	if(itemp)
	  //
	  break;
      }
    }

    switch(connect.size()) {
    case 0:
      //
      fac.push_back(FreqGraph());

      fac.back()[git->first] = git->second;

      break;

    case 1:
      //
      fac[*connect.begin()][git->first] = git->second;

      break;

    case 2:
      //
      fac[*connect.begin()][git->first] = git->second;

      for(const_iterator mit = fac[*connect.rbegin()].begin(); mit != fac[*connect.rbegin()].end(); ++mit)
	//
	fac[*connect.begin()][mit->first] = mit->second;
      
      fac.erase(fac.begin() + *connect.rbegin());

      break;

    default:
      //
      ErrOut err_out;

      err_out << funame << "wrong connect size: " << connect.size();
    }
  }
  
  _mg_t res;
  //
  for(int g = 0; g < fac.size(); ++g) {
    //
    std::set<int> vset;
    //
    for(const_iterator git = fac[g].begin(); git != fac[g].end(); ++git)
      //
      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	//
	vset.insert(*bit);

    std::map<int, int> vmap;

    itemp = 0;
    //
    for(std::set<int>::const_iterator it = vset.begin(); it != vset.end(); ++it, ++itemp)
      //
      vmap[*it] = itemp;
 
    FreqGraph perm_graph;
    //
    for(const_iterator git = fac[g].begin(); git != fac[g].end(); ++git) {
      //
      std::set<int> bond;

      for(std::set<int>::const_iterator bit = git->first.begin(); bit != git->first.end(); ++bit)
	//
	bond.insert(vmap[*bit]);
      
      perm_graph[bond] = git->second;
    }
    
    ++res[perm_graph];
  }

  return res;
}

Graph::FreqGraph  Graph::FreqGraph::reduce (const std::vector<double>& freq,
					    double                     temperature,  
					    const std::vector<double>& tanh_factor,
					    _mg_t&                     zpe_graph
					    ) const
{
  const char funame [] = "Graph::FreqGraph::reduce: ";

  _check_integrity();
  //
  _check_order();

  if(temperature <= 0.) {
    //
    ErrOut err_out;
    
    err_out << funame << "negative temperature, K: " << temperature / Phys_const::kelv;
  }

  zpe_graph.clear();
  
  int    itemp;
  double dtemp;

  const int vsize = vertex_size();

  // reduced graph
  //
  FreqGraph red_graph;

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

    for(const_iterator git = begin(); git != end(); ++git) {
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
    //
    std::vector<std::set<int> > merge_group;
    //
    for(const_iterator git = red_graph.begin(); git != red_graph.end(); ++git) {
      //
      dtemp = 0.;

      for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
	//
	if(*fit >= 0 && freq[*fit] > 0.)
	  //
	  dtemp += freq[*fit] * tanh_factor[*fit];
      }
      
      if(dtemp > temperature * red_thresh) {
	//
	std::multiset<int> merge;

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
	  //
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

    itemp = 0;
    //
    for(std::set<int>::const_iterator it = red_group[v].begin(); it != red_group[v].end(); ++it, ++itemp)
      //
      vertex_map[*it] = itemp;
	    
    FreqGraph new_graph;
    //
    for(const_iterator git = begin(); git != end(); ++git) {
      //
      std::set<int> bond;

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

double Graph::FreqGraph::zpe_factor  (const std::vector<double>& freq,
				      double                     temperature,
				      const std::vector<double>& tanh_factor
				      ) const
{
 const char funame [] = "Graph::FreqGraph::zpe_factor: ";

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
  std::set<FreqGraph> perm_set = perm_pool(&symm_fac);

  double res = 0.;

  for(std::set<FreqGraph>::const_iterator pit = perm_set.begin(); pit != perm_set.end(); ++pit) {
    //
    std::vector<double> freq_band(vsize - 1);
    //
    for(const_iterator git = pit->begin(); git != pit->end(); ++git) {
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

  for(const_iterator git = begin(); git != end(); ++git) {
    //
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      //
      // low frequency correlator value
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

// high temperature / low frequency integral evaluation
//
Graph::FreqGraph::_Ring::_Ring (const std::vector<int>& base) : std::vector<int>(base)
{
  const char funame [] = "Graph::FreqGraph::_Ring::_Ring: ";

  if(size() < 3) {
    //
    ErrOut err_out;
    
    err_out << funame << "ring size out of range: " << size();
  }

  std::set<int> test;

  for(const_iterator it = begin(); it != end(); ++it)
    //
    test.insert(*it);

  if(test.size() != size()) {
    //
    ErrOut err_out;
    
    err_out << funame << "identical vertices:";

    for(const_iterator it = begin(); it != end(); ++it)
      //
      err_out << "   " << *it;
  }  
}

void Graph::FreqGraph::_Ring::_assert (int i) const
{
  const char funame [] = "Graph::FreqGraph::_Ring::_assert: ";

  if(i < 0 || i >= size()) {
    //
    ErrOut err_out;

    err_out << funame << "index out of range: " << i;
  }
}

std::set<int> Graph::FreqGraph::_Ring::edge (int i) const
{
  _assert(i);
  
  std::set<int> res;

  res.insert((*this)[i]);

  res.insert((*this)[(i + 1) % size()]);
  
  return res;
}

int Graph::FreqGraph::_Ring::find (const std::set<int>& e) const
{
  const char funame [] = "Graph::FreqGraph::_Ring::find: ";

  if(e.size() != 2) {
    //
    ErrOut err_out;

    err_out << funame << "wrong edge size: " << e.size();
  }

  for(int i = 0; i < size(); ++i) {
    //
    int v0 = (*this)[i];
    
    int v1 = (*this)[(i + 1) % size()];

    std::set<int> test;

    test.insert(v0);

    test.insert(v1);

    if(e == test)
      //
      return v0 < v1 ? 1 : -1;
  }

  return 0;
}

bool Graph::FreqGraph::is_connected (int i0, int i1) const
{
  std::set<int> test;

  test.insert(i0);
  
  test.insert(i1);

  if(find(test) == end())
    //
    return false;

  return true;
}

Graph::FreqGraph::_Ring Graph::FreqGraph::_min_ring (const std::set<int>& edge) const
{
  const char funame [] = "Graph::FreqGraph::_min_ring: ";
  
  int    itemp;
  double dtemp;
  bool   btemp;

  if(find(edge) == end()){
    //
    ErrOut err_out;
    
    err_out << funame << "edge does not exist: (";
    
    for(std::set<int>::const_iterator sit = edge.begin(); sit != edge.end(); ++sit) {
      //
      if(sit != edge.begin())
	//
	err_out << ", ";

      err_out << *sit;
    }
    
    err_out << ")";
  }

  std::list<int> vertex_pool;
  
  for(const_iterator git = begin(); git != end(); ++git)
    //
    for(std::set<int>::const_iterator sit = git->first.begin(); sit != git->first.end(); ++sit)
      //
      if(edge.find(*sit) == edge.end())
	//
	vertex_pool.push_back(*sit);

  std::list<std::vector<int> > ring_pool;

  ring_pool.push_back(std::vector<int>(1, *edge.begin()));

  while(1) {
    //
    std::list<std::vector<int> > new_ring_pool;
    
    for(std::list<int>::iterator pit = vertex_pool.begin(); pit != vertex_pool.end();) {
      //
      btemp = true;
      
      for(std::list<std::vector<int> >::iterator lit = ring_pool.begin(); lit != ring_pool.end(); ++lit)
	//
	if(is_connected(lit->back(), *pit)) {
	  //
	  btemp = false;

	  new_ring_pool.push_back(*lit);
	  
	  new_ring_pool.back().push_back(*pit);

	  break;
	}

      if(btemp) {
	//
	++pit;
      }
      else {
	//
	pit = vertex_pool.erase(pit);
      }
    }

    if(!new_ring_pool.size())
      //
      throw _NoRing();

    ring_pool = new_ring_pool;

    for(std::list<std::vector<int> >::const_iterator lit = ring_pool.begin(); lit != ring_pool.end(); ++lit)
      //
      if(is_connected(lit->back(), *edge.rbegin())) {
	//
	std::vector<int> base = *lit;
	  
	base.push_back(*edge.rbegin());

	return _Ring(base);
      }
  }
}

double Graph::FreqGraph::_four_term (int                        mi,
				     double                     temperature,
				     const std::vector<double>& rfreq,
				     const std::vector<double>& ifreq) const
{
  const char funame [] = "Graph::FreqGraph::_four_term: ";

  static const double eps = 1.e-3;
  
  int    itemp;
  double dtemp;
  bool   btemp;

  if(mi < 0)
    //
    mi = -mi;

  const double mipi2 = M_PI * M_PI * double(mi * mi);

  // zero frequency
  //
  if(!rfreq.size() && !ifreq.size()) {
    //
    if(!mi)
      //
      return 0.;

    return 0.25 / temperature / mipi2;
  }

  double res = 0.;
  
  if(!rfreq.size() || !ifreq.size()) {
    //
    const std::vector<double>* fp = &rfreq;

    if(ifreq.size())
      //
      fp = &ifreq;
    
    const int nf = 1 << (fp->size() - 1);

    for(int fi = 0; fi < nf; ++fi) {
      //
      double ff = (*fp)[0];

      itemp = fi;
      
      for(int i = 1; i < fp->size(); ++i, itemp >>= 1) {
	//
	int bit = itemp % 2;

	if(bit) {
	  //
	  ff += (*fp)[i];
	}
	else
	  //
	  ff -= (*fp)[i];
      }

      ff /= 2. * temperature;

      if(rfreq.size()) {
	//
	dtemp = ff * ff + mipi2;

	if(!mi && dtemp < eps * eps) {
	  //
	  res += 1.;
	}
	else
	  //
	  res += std::sinh(ff) * ff / dtemp;
      }
      else {
	//
	dtemp = ff * ff - mipi2;

	if(!mi && dtemp < eps * eps) {
	  //
	  res += 1.;
	}
	else if(mi && dtemp < eps && dtemp > -eps) {
	  //
	  if(mi % 2) {
	    //
	    res -= 0.5;
	  }
	  else
	    //
	    res += 0.5;
	}
	else
	  //
	  res += std::sin(ff) * ff / dtemp;
      }
    }

    res /= (double)nf;

    return res;
  }
  // complex frequencies
  //
  else {
    //
    const int nf = 1 << (rfreq.size() + ifreq.size() - 2);

    for(int fi = 0; fi < nf; ++fi) {
      //
      std::complex<double> ff(rfreq[0], ifreq[0]);

      itemp = fi;
      
      for(int i = 1; i < rfreq.size(); ++i, itemp >>= 1) {
	//
	int bit = itemp % 2;

	if(bit) {
	  //
	  ff += rfreq[i];
	}
	else
	  //
	  ff -= rfreq[i];
      }

      for(int i = 1; i < ifreq.size(); ++i, itemp >>= 1) {
	//
	int bit = itemp % 2;

	if(bit) {
	  //
	  ff += std::complex<double>(0., ifreq[i]);
	}
	else
	  //
	  ff -= std::complex<double>(0., ifreq[i]);
      }

      ff /= 2. * temperature;

      std::complex<double> ctemp = ff * ff + mipi2;

      dtemp = std::abs(ctemp);
      
      if(!mi && dtemp < eps * eps) {
	//
	res += 1.;
      }
      else if(mi && dtemp < eps) {
	//
	if(mi % 2) {
	  //
	  res -= 0.5;
	}
	else
	  //
	  res += 0.5;
      }
      else {
	//
	ctemp = std::sinh(ff) * ff / ctemp;

	res += ctemp.real();
      }
    }

    res /= (double)nf;

    return res;
  }
}

double Graph::FreqGraph::fourier_sum (const std::vector<double>& freq, double temperature) const
{
  const char funame [] = "Graph::FreqGraph::fourier_sum: ";
  
  int    itemp;
  double dtemp;
  bool   btemp;

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

  // real frequencies
  //
  std::vector<std::vector<double> > rfreq(size());
  
  //imaginary frequencies
  //
  std::vector<std::vector<double> > ifreq(size());

  // low frequencies
  //
  std::vector<int>                  nfreq(size());

  itemp = 0;
  //
  for(const_iterator git = begin(); git != end(); ++git, ++itemp) {
    //
    if(!git->second.size()) {
      //
      ErrOut err_out;

      err_out << funame << "no frequencies at " << itemp << " edge: ";

      print(err_out);
    }
    
    for(std::multiset<int>::const_iterator fit = git->second.begin(); fit != git->second.end(); ++fit) {
      //
      freq_map.push_back(*fit);

      if(*fit < 0) {
	//
	++nfreq[itemp];
      }
      else {
	//
	dtemp = freq[*fit];
      
	if(dtemp > 0.) {
	  //
	  rfreq[itemp].push_back(dtemp);
	}
	else
	  //
	  ifreq[itemp].push_back(-dtemp);
      }
    }
  }
  
  double res = 0.;
  
  std::set<std::set<int> > edge_pool;

  for(const_iterator git = begin(); git != end(); ++git)
    //
    edge_pool.insert(git->first);

  // set of rings
  //
  std::list<_Ring>  ring_set;
  //
  while(edge_pool.size()) {
    //
    try {
      //
      _Ring ring = _min_ring(*edge_pool.begin());
    
      for(int i = 0; i < ring.size(); ++i)
	//
	edge_pool.erase(ring.edge(i));
    
      ring_set.push_back(ring);
    }
    catch(_NoRing) {
      //
      edge_pool.erase(edge_pool.begin());
    }
  }

  // purging set of rings
  //
  for(std::list<_Ring>::iterator rit0 = ring_set.begin(); rit0 != ring_set.end();) {
    //
    for(int i = 0; i < rit0->size(); ++i) {
      //
      btemp = true;
      
      for(std::list<_Ring>::iterator rit1 = ring_set.begin(); rit1 != ring_set.end(); ++rit1) {
	//
	if(rit1 != rit0 && rit1->find(rit0->edge(i))) {
	  //
	  btemp = false;
	  
	  break;
	}
      }
      
      if(btemp)
	//
	break;
    }

    if(btemp) {
      //
      ++rit0;
    }
    else
      //
      rit0 = ring_set.erase(rit0);
  }

  std::vector<std::map<int, int> > index_map;

  for(const_iterator git = begin(); git != end(); ++git) {
    //
    std::map<int, int> im;
    
    itemp = 0;
    //
    for(std::list<_Ring>::const_iterator rit = ring_set.begin(); rit != ring_set.end(); ++rit, ++itemp) {
      //
      int s = rit->find(git->first);

      if(s)
	//
	im[itemp] = s;
    }

    index_map.push_back(im);
  }

  itemp = ring_set.size();

  for(int i = 0; i < size(); ++i) {
    //
    if(rfreq[i].size() || ifreq[i].size())
      //
      ++itemp;

    itemp += nfreq[i] - 1;
  }

  const int idim = itemp;
    
  if(!idim) {
    //
    res = 1.;
    
    for(int i = 0; i < size(); ++i)
      //
      res *= _four_term(0, temperature, rfreq[i], ifreq[i]);
  }
  else {
    //
    for(MultiIndex multi(idim, 2 * four_max); !multi.end(); ++multi) {
      //
      double gvalue = 1.;

      int ci = size();
	
      for(int i = 0; i < size(); ++i) {
	//
	int mindex = 0;
	
	for(std::map<int, int>::const_iterator mit = index_map[i].begin(); mit != index_map[i].end(); ++mit)
	  //
	  mindex += mit->second * (multi[mit->first] - four_max);

	itemp = -1;
	  
	if(rfreq[i].size() || ifreq[i].size())
	  //
	  ++itemp;

	itemp += nfreq[i];

	const int nf = itemp;
	  
	for(int j = 0; j < nf; ++j, ++ci) {
	  //
	  itemp = multi[ci] - four_max;
	    
	  mindex -= itemp ;

	  gvalue *= _four_term(itemp, temperature);
	}

	gvalue *= _four_term(mindex, temperature, rfreq[i], ifreq[i]);
      }

      res += gvalue;
    }
  }
  // normalization
  //
  for(int i = 0; i < size(); ++i) {
    //
    for(int f = 0; f < rfreq[i].size(); ++f)
      //
      res /= 2. * rfreq[i][f] * std::sinh(rfreq[i][f] / 2. / temperature);

    for(int f = 0; f < ifreq[i].size(); ++f)
      //
      res /= -2. * ifreq[i][f] * std::sin(ifreq[i][f] / 2. / temperature);
  }

  res /= std::pow(temperature, (double)vsize);
    
  // old method
  //
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
      
    if(dtemp <= -1.) {
      //
      ErrOut err_out;

      err_out << funame << "deep tunneling regime: "
	      << "frequency[1/cm] = " << freq[itemp] / Phys_const::incm << ",  "
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
    
  double old_res = 0.;

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
    for(const_iterator git = begin(); git != end(); ++git) {
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
    
    old_res += dtemp;
  }
  
  old_res /= std::pow(temperature, vsize) * std::pow(4. * M_PI * M_PI * temperature, freq_map.size());

  if(old_res != 0.) {
    //
    dtemp = res / old_res - 1.;

    if(dtemp > 0.1 || dtemp < -0.1) {
      //
      std::cerr << funame << "WARNING: Graph = ";

      print(std::cerr);

      std::cerr << "\n"
		<< funame << "low frequency fourier sum = " << std::left << std::setw(13) << old_res
		<< "    standard fourier sum = " << std::setw(13) << res << "\n";
    }
  }
    
  return res;
}

void Graph::read_potex (const std::vector<double>& freq, std::istream& from, std::map<std::multiset<int>, double>& potex)
{
  const std::string funame = "Graph::read_potex: ";

  int    itemp;
  double dtemp;
  
  if(!freq.size()) {
    //
    ErrOut err_out;
    
    err_out << funame << "no frequencies";
  }

  std::vector<double> freq_sqrt(freq.size());
  
  for(int i = 0; i < freq.size(); ++i) {
    //
    if(freq[i] <= 0.) {
      //
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
      //
      ErrOut err_out;

      err_out << funame << "not enough values (linear and quadratic terms not allowed)";
    }

    double value = double(string_value.back()) * Phys_const::incm;

    const int imax = string_value.size() - 1;
    
    std::multiset<int> index;
    //
    for(int i = 0; i < imax; ++i) {
      //
      itemp = int(string_value[i]) - 1;

      if(itemp < 0 || itemp >= freq.size()) {
	//
	ErrOut err_out;

	err_out << funame << "index out of range: " << itemp;
      }

      index.insert(itemp);

      value *= freq_sqrt[itemp];
    }

    if(potex.find(index) != potex.end()) {
      //
      ErrOut err_out;

      err_out << funame << "term already in the expansion";
    }
    
    potex[index] = value;
  }

  if(!from) {
    //
    ErrOut err_out;

    err_out << funame << "corrupted";
  }
}

std::set<Graph::GenGraph> Graph::_raw_graph_generator (std::vector<int> vertex_order, int root)
{
  const char funame [] = "Graph::_raw_graph_generator: ";
  
  int itemp;

  itemp = 0;
  //
  for(int i = 0; i < vertex_order.size(); ++i)
    //
    itemp += vertex_order[i];
  
  if(itemp % 2) {
    //
    ErrOut err_out;
    //
    err_out << funame << "odd number of connction points";
  }
  
  std::set<GenGraph> res;

  if(root < 0) {
    //
    std::set<GenGraph> add;
    
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

    for(std::set<GenGraph>::const_iterator at = add.begin(); at != add.end(); ++at) {
      //
      std::vector<int> new_order = vertex_order;
      //
      for(GenGraph::const_iterator pit = at->begin(); pit != at->end(); ++pit) {
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
      
      std::set<GenGraph> g  = _raw_graph_generator(new_order);

      if(g.size()) {
	//
	for(std::set<GenGraph>::const_iterator git = g.begin(); git != g.end(); ++git) {
	  //
	  GenGraph gval = *git;
	  //
	  for(GenGraph::const_iterator pit = at->begin(); pit != at->end(); ++pit)
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
	
	std::set<GenGraph> g  = _raw_graph_generator(vertex_order, root);
	
	++vertex_order[i];

	std::multiset<int> bond;
	//
	bond.insert(root);
	//
	bond.insert(i);

	if(g.size()) {
	  //
	  for(std::set<GenGraph>::const_iterator git = g.begin(); git != g.end(); ++git) {
	    //
	    GenGraph gval = *git;
	    //
	    gval.insert(bond);
	    
	    res.insert(gval);
	  }
	}
	else {
	  //
	  GenGraph gval;
	  //
	  gval.insert(bond);
	  
	  res.insert(gval);
	} 
      }
    }
  }
  
  return res;
}
  
