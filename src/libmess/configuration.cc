/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
*/

#include "configuration.hh"
#include "random.hh"
#include "lapack.hh"
#include "io.hh"

namespace Configuration {
  
  Layout State::layout;

  double Simplexation::_det_tol        = 1.e-5;  // simplex determinant tolerance
  double Simplexation::_dis_tol        = 1.e-10; // distance tolerance
  int    Simplexation::_miss_count_max = 10000;  // maximal number of random misses

}

/****************************************************************************************************************************
 ********************************************* CONFIGURATIONAL SPACE VECTOR LAYOUT ******************************************
 ****************************************************************************************************************************/

// first size describes the linear subspace dimension;
// other sizes correspond to 3- and 4-dimensional spheres

void Configuration::Layout::set (const std::vector<int>& lt) 
{
  const char funame [] = "Configuration::Layout::Layout: ";

  _layout = lt;

  int itemp = 0;
  for(int i = 0; i < _layout.size(); itemp += _layout[i++])
    _shift.push_back(itemp); 

  _linear_size  = itemp;
  _simplex_size = itemp - size() + 2;

  // second fragment atomic
  if(size() == 1 && size(0) == 3 ||
     size() == 2 && size(0) == 1 &&  size(1) == 3 ||
     size() == 2 && size(0) == 0 &&  size(1) == 3)
    _fragment_type = ATOMIC;
  // second fragment linear
  else if(size() == 2 && size(0) == 3 &&  size(1) == 3 ||
	  size() == 3 && size(0) == 1 &&  size(1) == 3 && size(2) == 3 ||
	  size() == 3 && size(0) == 0 &&  size(1) == 3 && size(2) == 3)
    _fragment_type = LINEAR;
  // second fragment nonlinear
  else if(size() == 2 && size(0) == 3 &&  size(1) == 4 ||
	  size() == 3 && size(0) == 1 &&  size(1) == 3 && size(2) == 4 ||
	  size() == 3 && size(0) == 0 &&  size(1) == 3 && size(2) == 4)
    _fragment_type = NONLINEAR;
  // unknown layout
  else {
    std::cerr << funame << "unknown layout\n";
    throw Error::Range();
  }
}

/***************************************************************************************************************************
 ***************************** SYMMETRY GROUP FOR NONLINEAR + LINEAR AND/OR ATOMIC FRAGMENTS *******************************
 ***************************************************************************************************************************/

Configuration::SpaceGroup::SpaceGroup (const Symmetry::SpaceGroup& sg, int symmetric)
 
  : _base(sg.base()), _symmetric(symmetric)
{
  const char funame [] = "Configuration::SpaceGroup::SpaceGroup: ";

  int itemp;

  // linear fragment with inversion symmetry element
  if(_symmetric) {
    std::vector<int> perm(sg.size() * 2);
    std::vector<Permutation> table;

    for(int e = 0; e < sg.size(); ++e) {
      for(int g = 0; g < sg.size(); ++g) {
	itemp = sg.product(e, g);
	perm[g]             = itemp;
	perm[g + sg.size()] = itemp + sg.size();
      }
      table.push_back(Permutation(perm));
    }

    for(int e = 0; e < sg.size(); ++e) {
      for(int g = 0; g < sg.size(); ++g) {
	itemp = sg.product(e, g);
	perm[g]             = itemp + sg.size();
	perm[g + sg.size()] = itemp;
      }
      table.push_back(Permutation(perm));
    }

    init(table);
  }
  // linear or atomic fragments symmetry group
  else
    init(sg.multiplication_table());
}
    
void Configuration::SpaceGroup::apply (int g, const State& v, State& res) const 
{
  const char funame [] = "Configuration::SpaceGroup::apply: ";

  int itemp;

  // linear fragment with inversion symmetry element
  if(_symmetric) {

    // check layout
    if(State::layout.second_fragment() != Layout::LINEAR) {
      std::cerr << funame << "wrong layout\n";
      throw Error::Logic();
    }

    itemp = g % _base.size();
    _base[itemp].apply(v.orientation(),   res.orientation());
    _base[itemp].apply(v.radius_vector(), res.radius_vector()); 

    // second fragment identical atoms exchange operation
    if(g / _base.size()) {
      double* end  = res.orientation() + 3;
      for(double* it = res.orientation(); it != end; ++it)
	*it = - *it;
    }    
  }
  // symmetry group for atomic and linear fragmens
  else
    switch(State::layout.second_fragment()) {
    case Layout::ATOMIC:

    _base[g].apply(v.radius_vector(), res.radius_vector());

    break;
    case Layout::LINEAR:

      _base[g].apply(v.orientation(),   res.orientation()); 
      _base[g].apply(v.radius_vector(), res.radius_vector()); 

      break;
    default:
      std::cerr << funame << "wrong layout\n";
      throw Error::Logic();
    }

  // copy invariant part of the state vector
  for(int i = 0; i < State::layout.invariant_size(); ++i)
    res[i] = v[i];
}

/***************************************************************************************************************************
 ************************************** SYMMETRY GROUP FOR NONLINEAR + NONLINEAR FRAGMENTS *********************************
 ***************************************************************************************************************************/

Configuration::DoubleSpaceGroup::DoubleSpaceGroup (const std::vector<Symmetry::SpaceGroup>& sg, int identical)
  : _identical(identical)
{
  const char funame [] = "Configuration::DoubleSpaceGroup::DoubleSpaceGroup: ";

  int    itemp;
  double dtemp;

  if(sg.size() != 2) {
    std::cerr << funame << "fragment symmetry groups number\n";
    throw Error::Logic();
  }

  for(int frag = 0; frag < 2; ++frag)
    _symmetry_element.push_back(sg[frag].base());

  //check if the fragments symmetry groups are identical
  if(_identical && sg[0] != sg[1]) {
    std::cerr << funame << "identical fragments should have identical symmetry groups\n";
    throw Error::Logic();
  }

  // multi-index to symmetry group index map
  std::map<std::vector<int>, int> multi_index_map;

  // identical fragments
  if(_identical) {
    std::vector<int> multi(4);

    // symmetry group index to multi-index map
    for(int g0 = 0; g0 < sg[0].size(); ++g0)
      for(int g1 = 0; g1 < sg[0].size(); ++g1)
	if(sg[0][g0].inversion() &&  sg[0][g1].inversion() ||
	  !sg[0][g0].inversion() && !sg[0][g1].inversion()) {

	  // fragments symmetry groups indices
	  multi[0] = g0; 
	  multi[1] = g1;

	  for(int i = 0; i < 2; ++i) {
	    // orientational quaternion of the second fragment sign change 
	    multi[SIGN] = 1 - 2 * i; 

	    for(int x = 0; x < 2; ++x) {
	      // identical fragments exchange
	      multi[EXCHANGE] = x;

	      // add multi-index to the map
	      multi_index_map[multi] = _index_multi_map.size();
	      _index_multi_map.push_back(multi);
	    }
	  }
	}

    // multiplication table
    std::vector<Permutation> table;
    for(int g = 0; g < _index_multi_map.size(); ++g) {
      std::vector<int> perm(_index_multi_map.size());
      for(int f = 0; f < _index_multi_map.size(); ++f) {

	//  multi-index of the product of the symmetry group elements  
	switch(_index_multi_map[g][EXCHANGE]) {
	case 0: // (q1, p1) * (q2, p2) * X = (q1*q2, p1*p2) * X; (q1, p1) * (q2, p2) = (q1*q2, p1*p2)

	  // fragment symmetry group elements products indices
	  for(int i = 0; i < 2; ++i)
	    multi[i] = sg[0].product(_index_multi_map[g][i], _index_multi_map[f][i]);

	  // sign change
	  itemp = _index_multi_map[g][SIGN] * _index_multi_map[f][SIGN];
	  for(int i = 0; i < 2; ++i)
	    itemp *= sg[0].sign(_index_multi_map[g][i], _index_multi_map[f][i]);
	  multi[SIGN] = itemp;
	
	  // fragment exchange
	  multi[EXCHANGE] = _index_multi_map[f][EXCHANGE];
	  break;

	case 1: // (q1, p1) * X * (q2, p2) * X = (q1*p2, p1*q2);  (q1, p1) * X * (q2, p2) = (q1*p2, p1*q2) * X

	  // fragment symmetry group elements products indices
	  for(int i = 0; i < 2; ++i)
	    multi[i] = sg[0].product(_index_multi_map[g][i], _index_multi_map[f][1 - i]);

	  // sign change
	  itemp = _index_multi_map[g][SIGN] * _index_multi_map[f][SIGN];
	  for(int i = 0; i < 2; ++i)
	    itemp *= sg[0].sign(_index_multi_map[g][i], _index_multi_map[f][1 - i]);
	  multi[SIGN] = itemp;
	
	  // fragment exchange
	  multi[EXCHANGE] = 1 - _index_multi_map[f][EXCHANGE];
	  break;
	default:
	  std::cerr << funame << "identical fragments exchange index out of range\n";
	  throw Error::Range();
	}

	// add index to multiplication table
	std::map<std::vector<int>, int>::const_iterator mit = multi_index_map.find(multi);
	if(mit == multi_index_map.end()) {
	  std::cerr << funame << "did not find index in the map\n";
	  throw Error::Logic();
	}
	perm[f] = mit->second;
      }
      table.push_back(Permutation(perm));
    }

    init(table);
  } 
  // different fragments
  else {
    std::vector<int> multi(3);

    // symmetry group index to multi-index map 
    for(int g0 = 0; g0 < sg[0].size(); ++g0)
      for(int g1 = 0; g1 < sg[1].size(); ++g1)
	if(sg[0][g0].inversion() &&  sg[1][g1].inversion() ||
	  !sg[0][g0].inversion() && !sg[1][g1].inversion()) {

	  // fragment symmetry groups indices
	  multi[0]    = g0;
	  multi[1]    = g1;

	  for(int i = 0; i < 2; ++i) {
	    // second fragment quaternion sign change
	    multi[SIGN] = 1 - 2 * i;

	    // add multi-index to the map
	    multi_index_map[multi] = _index_multi_map.size();
	    _index_multi_map.push_back(multi);
	  }
	}

    // multiplication table
    std::vector<Permutation> table;
    for(int g = 0; g < _index_multi_map.size(); ++g) {

      // g * f = P(f) permutation
      std::vector<int> perm(_index_multi_map.size());
      for(int f = 0; f < _index_multi_map.size(); ++f) {

	// fragment symmetry group elements products indices
	for(int i = 0; i < 2; ++i)
	  multi[i] = sg[i].product(_index_multi_map[g][i], _index_multi_map[f][i]);

	// sign change
	itemp = _index_multi_map[g][SIGN] * _index_multi_map[f][SIGN];
	for(int i = 0; i < 2; ++i)
	  itemp *= sg[i].sign(_index_multi_map[g][i], _index_multi_map[f][i]);
	multi[SIGN] = itemp;
	
	// add index to multiplication table
	std::map<std::vector<int>, int>::const_iterator mit = multi_index_map.find(multi);
	if(mit == multi_index_map.end()) {
	  std::cerr << funame << "did not find index in the map\n";
	  throw Error::Logic();
	}
	perm[f] = mit->second;
      }
      table.push_back(Permutation(perm));
    }

    init(table);
  }

  // check that it is actually the group
  /*
  Lapack::Matrix test;
  for(int g0 = 0; g0 < size(); ++g0)
    for(int g1 = 0; g1 < size(); ++g1) {

      // qmatrix test
      test = qmatrix(g0) * qmatrix(g1);
      test -= qmatrix(product(g0, g1));

      for(int i = 0; i < test.size(); ++i)
	for(int j = 0; j < test.size(); ++j) {
	  dtemp = test(i, j);
	  dtemp = dtemp >= 0. ? dtemp : -dtemp;
	  if(dtemp > 1.e-10) {
	    std::cerr << funame << "qmatrix group inconsistency\n";
	    throw Error::Logic();
	  }
	}

      // rmatrix test
      test = rmatrix(g0) * rmatrix(g1);
      test -= rmatrix(product(g0, g1));

      for(int i = 0; i < test.size(); ++i)
	for(int j = 0; j < test.size(); ++j) {
	  dtemp = test(i, j);
	  dtemp = dtemp >= 0. ? dtemp : -dtemp;
	  if(dtemp > 1.e-10) {
	    std::cerr << funame << "rmatrix group inconsistency\n";
	    throw Error::Logic();
	  }
	}
    }
  */
}

void Configuration::DoubleSpaceGroup::apply (int g, const State& v, State& res) const 
{
  const char funame [] = "Configuration::DoubleSpaceGroup::apply: ";

  //check layout
  if(State::layout.second_fragment() != Layout::NONLINEAR) {
    std::cerr << funame << "wrong layout\n";
    throw Error::Logic();
  }

  // identical fragments and exchange operation is on
  if(_identical && _index_multi_map[g][EXCHANGE]) {

    // fragment exchange symmetry operation on the second fragment orientation (quaternion)
    Quaternion q(v.orientation());
    for(Quaternion::iterator it = q.begin() + 1; it != q.end(); ++it)
      *it = - *it;
   
    // spatial symmetry group operation on the second fragment orientation (quaternion)
    Quaternion qtemp = (const Quaternion&)(_symmetry_element[0][_index_multi_map[g][0]]) * q
      / (const Quaternion&)(_symmetry_element[0][_index_multi_map[g][1]]);
    
    // (optional) sign change of the second fragment orientation (quaternion)
    double* rit =  res.orientation();
    switch(_index_multi_map[g][SIGN]) {
    case  1:
      for(const double* it = qtemp.begin(); it != qtemp.end(); ++it, ++rit)
	*rit = *it;
      break;
    case -1:
      for(const double* it = qtemp.begin(); it != qtemp.end(); ++it, ++rit)
	*rit = -*it;
      break;
    default:
      std::cerr << funame << "wrong sign in the index map\n";
      throw Error::Logic();
    }
    
    // fragment exchange symmetry operation on the radius-vector
    Quaternion r;
    const double* vit = v.radius_vector();
    for(double* it = r.begin() + 1; it != r.end(); ++it, ++vit)
      *it = - *vit;
    
    r = q * r / q;

    // spatial symmetry group operation on the radius-vector
    _symmetry_element[0][_index_multi_map[g][0]].apply(r.begin() + 1, res.radius_vector()); 
  }
  // no fragments exchange operation
  else {

    // spatial symmetry operation on the second fragment orientation vector (quaternion)
    Quaternion qtemp = (const Quaternion&)(_symmetry_element[0][_index_multi_map[g][0]]) * v.orientation()
      / (const Quaternion&)(_symmetry_element[1][_index_multi_map[g][1]]);
    
    // (optional) sign change of the second fragment orientation vector (quaternion)
    double* rit =  res.orientation();
    switch(_index_multi_map[g][SIGN]) {
    case  1:
      for(const double* it = qtemp.begin(); it != qtemp.end(); ++it, ++rit)
	*rit = *it;
      break;
    case -1:
      for(const double* it = qtemp.begin(); it != qtemp.end(); ++it, ++rit)
	*rit = -*it;
      break;
    default:
      std::cerr << funame << "wrong sign in the index map\n";
      throw Error::Logic();
    }
    
    // spatial symmetry group operation on the center-of-mass vector
    _symmetry_element[0][_index_multi_map[g][0]].apply(v.radius_vector(), res.radius_vector()); 
  }

  // copy the rest of variables
  for(int i = 0; i < State::layout.invariant_size(); ++i)
    res[i] = v[i];
}

// radius-vector transformation matrix
Lapack::Matrix Configuration::DoubleSpaceGroup::rmatrix (int g) const 
{
  const char funame [] = "Configuration::DoubleSpaceGroup::rmatrix: ";

  if(_identical) {
    std::cerr << funame << "not applicable to identical fragments\n";
    throw Error::Logic();
  }

  if(g < 0 || g >= size()) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  Lapack::Matrix res(3);

  for(int i = 0; i < 3; ++i)
    res.row(i) = _symmetry_element[0][_index_multi_map[g][0]].rmatrix().column(i);

  return res;
}

// quaternion orientation transformation matrix
Lapack::Matrix Configuration::DoubleSpaceGroup::qmatrix (int g) const 
{
  const char funame [] = "Configuration::DoubleSpaceGroup::qmatrix: ";

  if(_identical) {
    std::cerr << funame << "not applicable to identical fragments\n";
    throw Error::Logic();
  }

  if(g < 0 || g >= size()) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  Lapack::Matrix res(4);
  res = 0.;
  
  int k, m, n;
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j) {
      double qp = _symmetry_element[0][_index_multi_map[g][0]][i]
	*         _symmetry_element[1][_index_multi_map[g][1]][j];
      // i = 0 and j = 0
      if(!i && !j) {
	for(int k = 0; k < 4; ++k)
	  res(k, k) += qp;
      }
      // i = j != 0
      else if(i == j) {
	for(int k = 0; k < 4; ++k)
	  if(k != 0 && k != i)
	    res(k, k) -= qp;
	  else
	    res(k, k) += qp;
      }
      // i = 0 or j = 0
      else if(!i || !j) {
	res(i, j) += qp;
	res(j, i) -= qp;

	k = i + j;
	m = k % 3 + 1;
	n = m % 3 + 1;
	res(m, n) -= qp;
	res(n, m) += qp;
      }
      // i != 0 and j != 0
      else {
	res(i, j) += qp;
	res(j, i) += qp;

	k = i % 3 + 1;
	if(j != k) {
	  res(0, k) += qp;
	  res(k, 0) += qp;
	}
	else {
	  k = j % 3 + 1;
	  res(0, k) -= qp;
	  res(k, 0) -= qp;
	}
      }
    }

  if(_index_multi_map[g][SIGN] == -1)
    res *= -1.;

  return res;
}

/***************************************************************************************************************************
 ************************************ PARTITIONING CONFIGURATIONAL SPACE INTO SIMPLEXES ************************************
 ***************************************************************************************************************************/

int Configuration::Simplexation::_find_simplex_center (const std::set<int>& simplex, double rmax, State& center, 
						       double& radius, double& det) const
{
  const char funame [] = "Configuration::Simplexation::_find_simplex_center: ";

  int    itemp;
  double dtemp;

  State vtemp;

  if(simplex.size() != State::layout.simplex_size()) {
    std::cerr << funame << "wrong simplex size\n";
    throw Error::Logic();
  }

  vtemp = center;
  
  if(_iterate_simplex_center(simplex, center, det))
    return 1;

  radius = vdistance(center, _vertex[*simplex.begin()]);

  if(radius > rmax)
    return 1;

  if(State::layout.size() == 1 || State::layout.size() == 2 && State::layout.size(0) == 0)
    return 0;

  while (vdistance(vtemp, center) > _dis_tol) {
    vtemp = center;
    if(_iterate_simplex_center(simplex, center, det)) 
      return 1;
  }

  radius = vdistance(center, _vertex[*simplex.begin()]);

  return 0;
}

int Configuration::Simplexation::_iterate_simplex_center (const std::set<int>& simplex, State& center, double& det) const
{
  const char funame [] = "Configuration::Simplexation::_iterate_simplex_center: ";

  int    itemp;
  double dtemp;

  if(simplex.size() != State::layout.simplex_size()) {
    std::cerr << funame << "wrong simplex size\n";
    throw Error::Logic();
  }

  // linear equations for the simplex center 
  Lapack::Matrix lhs_mat(State::layout.linear_size());
  Lapack::Vector rhs_vec(State::layout.linear_size());
  
  double rhs_shift;
  int row_index = -1;
  // the conditions of equal distances between the simplex center and the vertices
  for(std::set<int>::const_iterator vit = simplex.begin(); vit != simplex.end(); ++vit, ++row_index) {

    const double* cdp = _vertex[*vit].begin();

    dtemp = 0.;
    for(const double* end = cdp + State::layout.size(0); cdp != end; ++cdp)
      dtemp += *cdp * *cdp;
    dtemp /= 2.;

    if(vit == simplex.begin())
      rhs_shift = dtemp;

    else {
      rhs_vec[row_index] = dtemp - rhs_shift;

      Slice<double> row(lhs_mat.row(row_index));
      row  = _vertex[*vit];
      row -= _vertex[*simplex.begin()];
      
      // normalize
      rhs_vec[row_index] /= normalize(row);
    }
  }

  // initial guess orthogonality conditions
  for(int s = 1; s < State::layout.size(); ++s, ++row_index) {

    Slice<double> row(lhs_mat.row(row_index));
    row = 0.;

    itemp = State::layout.shift(s) + State::layout.size(s);
    for(int i = State::layout.shift(s); i < itemp; ++i)
      row[i] = center[i];

    rhs_vec[row_index] = 1.;
  }

  try {
    Lapack::LU lu(lhs_mat);
    
    dtemp = lu.det();
    det = dtemp >= 0. ? dtemp : -dtemp;

    if(det < _det_tol)
      return 1;

    Lapack::Vector res = lu.invert(rhs_vec);

    for(int s = 1; s < State::layout.size(); ++s, ++row_index)
      normalize(res.begin() + State::layout.shift(s), State::layout.size(s));

    center = res;
  }
  catch(Error::General) {
    return 1;
  }

  return 0;
}

void Configuration::Simplexation::random_init (double dist_min, double rmax, double rmin) 
{
  const char funame [] = "Configuration::Simplexation::random_init: ";

  int    itemp;
  double dtemp, dmin, dmax;

  switch(State::layout.size(0)) {
  case 0:
    // do nothing
    break;
  case 1:
    if(rmax <= 0.) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }
    break;
  case 3:
    if(rmin <= 0. || rmax <= 0. || rmin >= rmax) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }
    break;
  default:
    std::cerr << funame << "unknown linear subspace dimension: " << State::layout.size(0) << "\n";
    throw Error::Logic();
  }

  int miss_count = 0;

  // vertex best separated from its symmetric images
  double vertex_dist_max = 0.;
  int    vertex_indx_max = -1;

  // random guess
  State guess;
  std::vector<State> guess_orbit(_symm_group->size(), guess);
    
  while(miss_count < _miss_count_max) {

    
    // linear subspace
    switch(State::layout.size(0)) {
    case 0:
      // do nothing
      break;
    case 1:
      guess[0] = Random::flat() * rmax;
      break;
    case 3:
      Random::spherical_layer(guess.begin(), 3, rmin, rmax);
      break;
    default:
      std::cerr << funame << "unknown linear subspace dimension: " << State::layout.size(0) << "\n";
      throw Error::Logic();
    }

    // spherical subspaces
    for(int s = 1;  s < State::layout.size(); ++s)
      Random::orient(guess.begin() + State::layout.shift(s), State::layout.size(s));

    bool miss = false;

    // check the distances to vertices
#pragma omp parallel for default(shared) schedule(static)	

    for(int v = 0;  v < _vertex.size(); ++v)
      if(!miss && vdistance(guess, _vertex[v]) < dist_min)
	miss = true;
 
    if(miss) {
      miss_count++;
      continue;
    }

    // generate guess orbit and check the distances to symmetric configurations
    for(int g = 1; g < _symm_group->size(); ++g) {

      _symm_group->apply(g, guess, guess_orbit[g]);

      dtemp = vdistance(guess, guess_orbit[g]);
      if(dtemp < dist_min) {
	  miss = true;
	  break;
      }
      if(g == 1 || dtemp < dmin)
	dmin = dtemp;
    }

    if(miss) {
      miss_count++;
      continue;
    }

    if(dmin > vertex_dist_max) {
      vertex_dist_max  = dmin;
      vertex_indx_max = _vertex.size();
    }
    
    // update vertex and vertex orbit arrays and vertex-orbit map;
    std::vector<int> orbit(_symm_group->size());
    for(int v = 0; v < _symm_group->size(); ++v) {
      orbit[v] = _vertex.size();
      _vertex_orbit_map.push_back(std::make_pair(_vertex_orbit.size(), v));
      _vertex.push_back(guess_orbit[v]);
    }
    // add orbit 
    _vertex_orbit.push_back(orbit);

    miss_count = 0;
  }

}

