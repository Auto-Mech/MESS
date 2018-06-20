/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include "symmetry.hh"
#include <set>

/***************************************************************************************************
 ************************************** ABSTRACT GROUP CLASS ***************************************
 ***************************************************************************************************/

void Symmetry::GroupBase::_isinit () const throw(Error::General)
{
  const char funame [] = "Symmetry::GroupBase::_isinit: ";

  if(_init)
    return;

  std::cerr << funame << "not initialized\n";
  throw Error::Init();
}

void Symmetry::GroupBase::init (const std::vector<Permutation>& table) throw(Error::General)
{
  const char funame [] = "Symmetry::GroupBase::init: ";

  _init = true;
  std::vector<Permutation>::operator=(table);

  if(!table.size()) {
    std::cerr << funame << "no table\n";
    throw Error::Init();
  }

  // check dimensions
  for(int e = 0; e < table.size(); ++e)
    if(table[e].size() != table.size()) {
      std::cerr << funame << "dimensions mismatch\n";
      throw Error::Logic();
    }

  //check that the first element is unity
  if(table[0] != Permutation(table.size())) {
    std::cerr << funame << "first element is not unity\n";
    throw Error::Logic();
  }

  // check associativity
  for(int e = 0; e < table.size(); ++e)
    for(int g = 0; g < table.size(); ++g)
      for(int f = 0; f < table.size(); ++f)
	if(table[table[e][g]][f] != table[e][table[g][f]]) {
	  std::cerr << funame << "associativity check failed\n";
	  throw Error::Logic();
	}
  // check solution uniqueness
  for(int e = 0; e < table.size(); ++e) {
    std::set<int> pool;
    for(int g = 0; g < table.size(); ++g)
      pool.insert(table[g][e]);
    if(pool.size() != table.size()) {
      std::cerr << funame << "solution is not unique\n";
      throw Error::Logic();
    }
  }
  
  // find inverse
  _inverse.resize(table.size(), -1);
  for(int e = 0; e < table.size(); ++e) {
    for(int g = 0; g < table.size(); ++g)
      if(table[e][g] == 0) {
	if(_inverse[e] >= 0) {
	  std::cerr << funame << "inverse already defined\n";
	  throw Error::Logic();
	}
	_inverse[e] = g;
      }

    if(_inverse[e] < 0) {
      std::cerr << funame << "inverse not found\n";
      throw Error::Logic();
    }
  } 
}

/**************************************************************************************************
 ********************************** SPATIAL SYMMETRY OPERATIONS ***********************************
 **************************************************************************************************/

Symmetry::SpaceElement::SpaceElement (const Quaternion& q, bool i, int flags) throw(Error::General)
  : Quaternion(q), _inversion(i), _rmatrix(q)
{
  const char funame [] = "Symmetry::SpaceElement::SpaceElement: ";
  
  if(_inversion)
    _rmatrix *= -1.;

  if(flags & Quaternion::NOCHECK)
    return;

  double dtemp = vdot(q) - 1.;
  dtemp = dtemp >= 0. ? dtemp : -dtemp;

  if(dtemp > Quaternion::tolerance) {
    std::cerr << funame << "quaternion not normalized\n";
    throw Error::Init();
  }
}

Symmetry::SpaceElement::SpaceElement (const D3::Matrix& m, int flags) throw(Error::General)
  : _inversion(m.inversion()), _rmatrix(m) 
{
  if(_inversion)
    mat2quat(-m, *this, flags);
}

Symmetry::SpaceElement Symmetry::SpaceElement::operator* (const SpaceElement& e) const
{
  return SpaceElement(Quaternion::operator*(e),
		       inversion() && !e.inversion() || 
		      !inversion() &&  e.inversion(), Quaternion::NOCHECK);  
}

int Symmetry::SpaceElement::compare (const SpaceElement& e) const
{
  if(inversion() != e.inversion())
    return  0;

  if(std::sqrt(vdot(*this - e)) < Quaternion::tolerance)
    return  1;

  if(std::sqrt(vdot(*this + e)) < Quaternion::tolerance)
    return -1;

  return 0;
}

void Symmetry::SpaceElement::apply (const double* v, double* res) const
{
  const char funame [] = "Symmetry::SpaceElement::apply: ";

  D3::vprod(v, _rmatrix, res);
}

/************************************************************************************************
 *********************************** SPATIAL SYMMETRY GROUP *************************************
 ************************************************************************************************/

Symmetry::SpaceGroup::SpaceGroup (const std::string& nm, const std::vector<SpaceElement>& group_base, int size_max) throw(Error::General)
  : _name(nm), std::vector<SpaceElement>(group_base), _sign(group_base.size(), group_base.size())
{
  const char funame [] = "Symmetry::SpaceGroup::SpaceGroup: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  // trivial group
  if(!std::vector<SpaceElement>::size()) {
    std::vector<SpaceElement>::push_back(unity());
    init(std::vector<Permutation>(1, Permutation(1)));
    return;
  }

  // check that the first element is unity
  if(*begin() != unity()) {
    std::cerr << funame << "firgst element of the group base should be unity\n";
    throw Error::Logic();
  }

  // check that all elements are different
  for(const_iterator i = begin(); i != end(); ++i)
    for(const_iterator j = begin(); j != i; ++j)
      if(compare(*i, *j)) {
	std::cerr << funame << "identical elements\n";
	throw Error::Logic();
      }

  // generate missing elements
  if(size_max) {
    do {
      // check the group size
      if(std::vector<SpaceElement>::size() > size_max) {
	std::cerr << funame << "number of elements in the group exceeded the limit\n";
	throw Error::Range();
      }

      // check if the product is in the group and add it if needed
      for(const_iterator i = begin(); i != end(); ++i) {
	for(const_iterator j = begin(); j != end(); ++j) {
	  // elements product
	  SpaceElement ij = *i * *j;
	  btemp = true;

	  // compare product with available elements
	  for(const_iterator k = begin(); k != end(); ++k)
	    if(compare(ij, *k)) {
	      btemp = false;
	      break;
	    }

	  // add not-found product to the group
	  if(btemp) {
	    std::vector<SpaceElement>::push_back(ij);
	    break;
	  }
	}
	// go to the next cycle
	if(btemp)
	  break;
      }
    } while(btemp);

    _sign.resize(std::vector<SpaceElement>::size(), std::vector<SpaceElement>::size());
  }

  // make a multiplication table
  std::vector<Permutation> table; 
  for(const_iterator e = begin(); e != end(); ++e) {
    int ee = e - begin();
    std::vector<int> perm(std::vector<SpaceElement>::size(), -1);

    for(const_iterator f = begin(); f != end(); ++f) {
      int ff = f - begin();
      // elements product
      SpaceElement ef = *e * *f;
      
      for(const_iterator g = begin(); g != end(); ++g) {
	int gg = g - begin();
	itemp = compare(ef, *g);
	if(itemp) {

	  if(perm[ff] >= 0) {
	    std::cerr << funame << "identical elements\n";
	    throw Error::Logic();
	  }
	  
	  // update multiplication table
	  perm[ff]      = gg;
	  _sign(ee, ff) = itemp;
	}
      }

      if(perm[ff] < 0) {
	std::cerr << funame << "product not found\n";
	throw Error::Logic();
      }
    }
    
    table.push_back(Permutation(perm));
  }

  // initialize the group
  init(table);
}
