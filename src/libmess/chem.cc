#include "chem.hh"
#include "io.hh"

#include <map>
#include <set>

/*************************************************************************
 ************************* TOLERANCES AND LIMITS *************************
 *************************************************************************/

double Chem::angle_tolerance = 5.0;

double Chem::distance_tolerance = 0.05;

bool Chem::are_angles_equal (double a1, double a2)
{
  double da = a2 - a1;
  
  if(da < angle_tolerance && da > -angle_tolerance) {
    //
    return true;
  }
  
  return false;
}

bool Chem::are_distances_equal (double a1, double a2)
{
  double da = a2 - a1;
  
  if(da < distance_tolerance && da > -distance_tolerance) {
    //
    return true;
  }
  
  return false;
}

double Chem::max_bond_length(const AtomBase& a1, const AtomBase& a2)
{
  if(a1.number() == AtomBase::HYDROGEN || a2.number() == AtomBase::HYDROGEN)
    //
    return 2.5;
  
  return 3.5;
}


/**********************************************************************************
 *********************************** GRAPH ****************************************
 **********************************************************************************/


// brute force symmetry group generator
//
std::set<Permutation> Chem::Graph::_default_symmetry_group () const
{
  const char funame [] = "Chem::Graph::_default_symmetry_group: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  std::vector<int> ivec;

  _isinit();
  
  if(!is_connected()) {
    //
    ErrOut err_out;
    
    err_out << funame << "graph is not connected";
  }

  std::set<Permutation> res;

  // permutations which leave vertex valences intact
  //
  std::map<int, std::set<int> > valence_map;

  std::map<int, std::set<int> >::const_iterator vmit;

  for(int v = 0; v < _vsize; ++v)
    //
    valence_map[_valence(v)].insert(v);

  ivec.resize(valence_map.size());

  itemp = 0;
  
  for(vmit = valence_map.begin(); vmit != valence_map.end(); ++vmit, ++itemp)
    //
    ivec[itemp] = vmit->second.size();

  // multi-permutation of valence groups
  //
  for(MultiPerm multi_perm(ivec); !multi_perm.end(); ++multi_perm) {
    //
    std::vector<int> perm(_vsize);

    int vgi = 0; // valence group index
    
    for(vmit = valence_map.begin(); vmit != valence_map.end(); ++vmit, ++vgi) {// valence map cycle
      //
      std::vector<int> forw(vmit->second.size());

      itemp = 0;
      
      for(std::set<int>::const_iterator sit = vmit->second.begin(); sit != vmit->second.end(); ++sit, ++itemp)
	//
	forw[itemp] = *sit;
      
      for(int i = 0; i < vmit->second.size(); ++i)
	//
	perm[forw[i]] = forw[multi_perm[vgi][i]];
      //
    }// valence map cycle

    // check if the permutation belongs to the symmetry group
    //
    std::set<std::set<int> > test;
    
    std::set<std::set<int> >::const_iterator git;

    for(git = begin(); git != end(); ++git) {
      //
      std::set<int> edge;

      for(std::set<int>::const_iterator sit = git->begin(); sit != git->end(); ++sit)
	//
	edge.insert(perm[*sit]);

      test.insert(edge);
    }

    if(test == *this)
      //
      res.insert(Permutation(perm));
    //
  }// multi-permutation cycle

  return res;
}

// graph projection on the subset of vertices
//
Chem::Graph Chem::Graph::_projection (const std::vector<int>& vset) const
{
  const char funame [] = "Chem::Graph::_projection: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  _isinit();
  
  std::map<int, int> back;

  std::set<int> vpool;
  
  for(int i = 0; i < vset.size(); ++i) {
    //
    itemp = vset[i];

    _assert(itemp);
    
    if(!vpool.insert(itemp).second) {
      //
      ErrOut err_out;

      err_out << funame << "identical vertices in the subset: " << itemp;
    }
    
    back[itemp] = i;
  }

  std::set<std::set<int> > res;
  
  std::set<std::set<int> >::const_iterator git;

  for(git = begin(); git != end(); ++git) {
    //
    std::set<int> edge;
    
    for(std::set<int>::const_iterator sit = git->begin(); sit != git->end(); ++sit) {
      //
      std::map<int, int>::const_iterator mit = back.find(*sit);

      if(mit != back.end())
	//
	edge.insert(mit->second);
      //
    }//

    if(edge.size() == 2)
      //
      res.insert(edge);
  }

  return Graph(vset.size(), res);
}

// find subset of nearest neighbors
//
std::set<int> Chem::Graph::_find_neighbor (int v) const
{
  const char funame [] = "Chem::Graph::_find_neighbor: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  _isinit();
  
  std::set<int> res;

  _assert(v);

  std::set<std::set<int> >::const_iterator git;

  for(git = begin(); git != end(); ++git) {
    //
    std::set<int>::const_iterator sit = git->find(v);
    
    if(sit != git->end()) {
      //
      if(sit != git->begin()) {
	//
	res.insert(*git->begin());
      }
      else
	//
	res.insert(*git->rbegin());
      //
    }//
    //
  }//
    
  return res;
}

void Chem::Graph::_isinit () const
{
  const char funame [] = "Chem::Graph::isinit: ";
  
  if(_vsize)
    //
    return;

  ErrOut err_out;

  err_out << funame << "not initialized";
}

// is graph connected
//
bool Chem::Graph::is_connected () const
{
  const char funame [] = "Chem::Graph::is_connected: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  _isinit();
  
  std::set<int> test;

  test.insert(0);

  while(1) {
    //
    std::set<int> next;// next connected layer
    
    for(std::set<int>::const_iterator sit = test.begin(); sit != test.end(); ++sit) {
      //
      std::set<int> near = _find_neighbor(*sit);

      for(std::set<int>::const_iterator nit = near.begin(); nit != near.end(); ++nit)
	//
	if(test.find(*nit) == test.end())
	  //
	  next.insert(*nit);
    }

    if(!next.size()) {
      //
      if(test.size() < _vsize) {
	//
	return false;
      }
      else
	//
	return true;
    }

    // merge with the rest
    //
    for(std::set<int>::const_iterator sit = next.begin(); sit != next.end(); ++sit)
      //
      test.insert(*sit);
  }
}

// bond distance between two vertices
//
int Chem::Graph::distance (int v0, int v1) const
{
  const char funame [] = "Chem::Graph::distance: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  _isinit();

  _assert(v0);

  _assert(v1);
  
  std::set<int> test;

  test.insert(v0);

  int res = 0;
  
  while(1) {
    //
    if(test.find(v1) != test.end())
      //
      return res;
    
    std::set<int> next;// next connected layer
    
    for(std::set<int>::const_iterator sit = test.begin(); sit != test.end(); ++sit) {
      //
      std::set<int> near = _find_neighbor(*sit);

      for(std::set<int>::const_iterator nit = near.begin(); nit != near.end(); ++nit)
	//
	if(test.find(*nit) == test.end())
	  //
	  next.insert(*nit);
    }

    if(!next.size()) {
      //
      if(test.size() < _vsize) {
	//
	IO::log << IO::log_offset << funame << "WARNING: graph is not connected";

	return -1;
      }
      else {
	//
	ErrOut err_out;

	err_out << funame << "logical error";
      }
    }
    
    // merge with the rest
    //
    for(std::set<int>::const_iterator sit = next.begin(); sit != next.end(); ++sit)
      //
      test.insert(*sit);

    ++res;
  }
}

// is graph a ring
//
bool Chem::Graph::_is_ring () const
{
  const char funame [] = "Chem::Graph::_is_ring: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  _isinit();

  for(int v = 0; v < _vsize; ++v)
    //
    if(_valence(v) != 2)
      //
      return false;
  
  return true;
}

void Chem::Graph::_assert (int v) const
{
  const char funame [] = "Chem::Graph::_assert: ";

  if(v < 0 || v >= _vsize) {
    //
    ErrOut err_out;

    err_out << funame << "vertex index out of range: " << v;
  }
}
 
// graph initialization
// 
void Chem::Graph::init (int v, const std::set<std::set<int> >& base)
{
  const char funame [] = "Chem::Graph::init: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  if(_vsize) {
    //
    ErrOut err_out;

    err_out << funame << "already initialized";
  }

  if(v <= 0) {
    //
    ErrOut err_out;

    err_out << funame << "vertex size out of range: " << v;
  }

  _vsize = v;
  
  (std::set<std::set<int> >&)*this = base;
  
  _assert();
}

 void Chem::Graph::_assert () const
 {
   const char funame [] = "Chem::Graph::_assert: ";

  std::set<std::set<int> >::const_iterator git;

  for(git = begin(); git != end(); ++git) {
    //
    if(git->size() != 2) {
      //
      ErrOut err_out;

      err_out << funame << "edge size out of range: " << git->size();
    }
  
    for(std::set<int>::const_iterator sit = git->begin(); sit != git->end(); ++sit)
      //
      _assert(*sit);
  }
}
  
// recursive linear branches removal
// 
std::set<Permutation> Chem::Graph::symmetry_group () const
{
  const char funame [] = "Chem::Graph::symmetry_group: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  _isinit();

  std::set<Permutation> res;

  if(size() == 1) {
    //
    res.insert(Permutation(1));

    return res;
  }
  
  if(!is_connected()) {
    //
    ErrOut err_out;
    
    err_out << funame << "graph is not connected";
  }

  std::set<int> vpool;
  
  for(int v = 0; v < _vsize; ++v)
    //
    vpool.insert(v);

  std::map<int, std::map<int, std::vector<std::vector<int> > > > tree;

  std::map<int, std::map<int, std::vector<std::vector<int> > > >::const_iterator ctit;// tree iterator

  std::map<int, std::vector<std::vector<int> > >::const_iterator cbit;// bush iterator

  // tree cycle
  //
  while(1) {
    //
    int root;

    std::vector<int> branch;

    // find the head of the branch
    //
    for(std::set<int>::iterator sit = vpool.begin(); sit != vpool.end(); ++sit) {
      //
      std::set<int> next = _find_neighbor(*sit);

      if(!next.size()) {
	//
	ErrOut err_out;

	err_out << funame << "logic error: there is a dangling vertex: " << *sit;
      }
      
      if(next.size() == 1) {
	//
	branch.push_back(*sit);

	vpool.erase(sit);
	
	root = *next.begin();

	break;
      }
    }

    if(!branch.size())
      //
      break;

    // branch cycle
    //
    while(1) {
      //
      std::set<int> next = _find_neighbor(root);

      // root of the branch found
      //
      if(next.size() > 2)
	//
	break;

      // linear graph
      //
      if(next.size() == 1) {
	//
	branch.push_back(root);

	vpool.erase(root);

	if(branch.size() != _vsize) {
	  //
	  ErrOut err_out;

	  err_out << funame << "logical error: branch size for linear graph differs from the total number of vertices: "
		  << branch.size() << " vs " << _vsize;
	}

	res.insert(Permutation(branch.size()));
	
	std::vector<int> perm(branch.size());

	for(int i = 0; i < perm.size(); ++i)
	  //
	  perm[branch[i]] = branch[perm.size() - 1 - i];

	res.insert(Permutation(perm));

	return res;
	//
      }// linear graph

      next.erase(branch.back());
      
      branch.push_back(root);

      vpool.erase(root);

      root = *next.begin();
      //
    }// branch cycle

    tree[root][branch.size()].push_back(branch);
    //
  }// tree cycle

  // graph is the rings system
  //
  if(!tree.size()) {
    //
    // graph is a simple ring
    //
    if(_is_ring()) {
      //
      std::vector<int> ring(_vsize);

      ring[0] = 0;

      std::set<int> next = _find_neighbor(0);

      for(int i = 1; i < _vsize; ++i) {
	//
	ring[i] = *next.begin();

	next = _find_neighbor(ring[i]);

	next.erase(ring[i - 1]);
	//
      }//

      std::vector<int> perm(_vsize);

      // rotation permutation
      //
      for(int i = 0; i < _vsize; ++i)
	//
	if(!i) {
	  //
	  perm[ring[i]] = ring.back();
	}
	else
	  //
	  perm[ring[i]] = ring[i - 1];

      res.insert(Permutation(perm));

      // reflection permutation
      //
      for(int i = 0; i < _vsize; ++i)
	//
	perm[ring[i]] = ring[_vsize - 1 - i];

      res.insert(Permutation(perm));

      return permutation_group(res);
      //
    }// graph is a simple ring

    return _default_symmetry_group();
    //
  }// graph is a system of rings

  // bush signature to root map
  //
  std::map<std::map<int, int>, std::set<int> > sign_map;
  
  std::map<std::map<int, int>, std::set<int> >::const_iterator smit;
  
  for(ctit = tree.begin(); ctit != tree.end(); ++ctit) {
    //
    std::map<int, int> sign;
    
    for(cbit = ctit->second.begin(); cbit != ctit->second.end(); ++cbit)
      //
      sign[cbit->first] = cbit->second.size();

    sign_map[sign].insert(ctit->first);
  }

  itemp = 0;

  std::vector<int>   forw(vpool.size());

  std::map<int, int> back;

  for(std::set<int>::const_iterator sit = vpool.begin(); sit != vpool.end(); ++sit, ++itemp) {
    //
    forw[itemp] = *sit;

    back[*sit] = itemp;
  }

  // symmetry group of the leftover subset of vertices
  //
  std::set<Permutation> sg = _projection(forw).symmetry_group();

  std::set<Permutation>::const_iterator sgit;

  // symmetry group cycle
  //
  for(sgit = sg.begin(); sgit != sg.end(); ++sgit) {
    //
    btemp = false;

    std::vector<int> perm(_vsize);

    for(std::set<int>::const_iterator sit = vpool.begin(); sit != vpool.end(); ++sit)
      //
      perm[*sit] = forw[(*sgit)[back[*sit]]];

    for(smit = sign_map.begin(); smit != sign_map.end(); ++smit) {
      //
      for(std::set<int>::const_iterator sit = smit->second.begin(); sit != smit->second.end(); ++sit) { 
	//
	int new_root = forw[(*sgit)[back[*sit]]];

	// the root does not have the right signuture
	//
	if(smit->second.find(new_root) == smit->second.end()) {
	  //
	  btemp = true;

	  break;
	}

	for(std::map<int, int>::const_iterator mit = smit->first.begin(); mit != smit->first.end(); ++mit) {
	  //
	  for(int n = 0; n < mit->second; ++n) {
	    //
	    for(int m = 0; m < mit->first; ++m) {
	      //
	      perm[tree[*sit][mit->first][n][m]] = tree[new_root][mit->first][n][m];
	    }//
	    //
	  }//
	  //
	}//
	//
      }//

      if(btemp)
	//
	break;

    }
    
    if(btemp)
      //
      continue;
    
    res.insert(Permutation(perm));
  }

  // add self-permutations for each bush
  //
  for(ctit = tree.begin(); ctit != tree.end(); ++ctit) {// bush root cycle
    //
    for(cbit = ctit->second.begin(); cbit != ctit->second.end(); ++cbit) {// branch length cycle
      //
      // permutations of branches of the same length with the same root
      //
      for(Permutation branch_perm(cbit->second.size()); !branch_perm.end(); ++branch_perm) {
	//
	std::vector<int> perm(_vsize);

	for(int i = 0; i < _vsize; ++i)
	  //
	  perm[i] = i;

	for(int n = 0; n < cbit->second.size(); ++n)
	  //
	  for(int m = 0; m < cbit->first; ++m) 
	    //
	    perm[tree[ctit->first][cbit->first][n][m]] = tree[ctit->first][cbit->first][branch_perm[n]][m];

	res.insert(Permutation(perm));
	//
      }// permutation of branches of the same length
      //
    }// branch length cycle
    //
  }// bush root cycle

  return permutation_group(res);
}
