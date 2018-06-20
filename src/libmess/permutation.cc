/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include <map>
#include <cstdlib>

#include "permutation.hh"
#include "io.hh"

/************************************************************************************************
 ************************************ INDEX PERMUTATION *****************************************
 ************************************************************************************************/

int Permutation::factorial (int n)
{
  const char funame [] = "Permutation::factorial: ";

  if(n > 10) {
    std::cerr << funame << "argument is too large: " << n << "\n";
    throw Error::Range();
  }
  
  if(n < 0) {
    std::cerr << funame << "negative argument: " << n << "\n";
    throw Error::Range();
  }
  
  switch(n) {
  case  0: return 1;
  case  1: return 1;
  case  2: return 2;
  case  3: return 6;
  case  4: return 24;
  default: return n * factorial(n - 1);
  }
}
  
Permutation::Permutation (int n)
  : _end(false)
{
  Exception::Base funame = "Permutation::Permutation: ";

  if(n < 0)
    throw funame << "out of range";

  if(!n) {
    _end = true;
    return;
  }

  resize(n);
  for(int i = 0; i < n; ++i)
    std::vector<int>::operator[](i) = i;
}

Permutation::Permutation (const std::vector<int>& p, int flags)
  : std::vector<int>(p), _end(false)
{
  Exception::Base funame = "Permutation::Permutation: ";

  if(!p.size()) {
    _end = true;
    return;
  }

  if(flags & NOCHECK) 
    return;

  std::set<int> pool;
  for(int i = 0; i < p.size(); ++i)
    if(!pool.insert(p[i]).second) {
      funame << "in permutation (";
      for(int j = 0; j < p.size(); ++j) {
	if(j)
	  funame << ", ";
	funame << p[j];
      }	    
      funame << ") " << i <<"-th element (" << p[i] << ") is repeated";
      throw funame;
    }

  if(*pool.begin() < 0 || *pool.rbegin() >= p.size())
    throw funame << "out of range";
}

// two elements permutation
Permutation::Permutation (int i1, int i2, int s)
  : std::vector<int>(s), _end(false)
{
  if(s < 2 || i1 < 0 || i2 < 0 || i1 >= s || i2 >= s || i1 == i2)
    throw Exception::Base() << "Permutation::Permutation: indices out of range: " << i1 << ", " <<  i2 << ", " << s;

  for(int i = 0; i < size(); ++i)
    if(i == i1)
      std::vector<int>::operator[](i) = i2;
    else if(i == i2)
      std::vector<int>::operator[](i) = i1;
    else
      std::vector<int>::operator[](i) = i;
}

// permutation from orbit
Permutation::Permutation (const std::set<std::vector<int> >& orbit, int s)
  : _end(false)
{
  Exception::Base funame = "Permutation::Permutation: ";
  
  if(s <= 0)
    throw funame << "non-positive size";

  resize(s);

  // elements pool
  std::set<int> pool;
  for(std::set<std::vector<int> >::const_iterator ot = orbit.begin(); ot != orbit.end(); ++ot) {
    if(!ot->size())
      throw funame << "zero orbit";

    for(const_iterator it = ot->begin(); it != ot->end(); ++it)
      if(!pool.insert(*it).second)
	throw funame << "duplicate elements";
  }

  if(pool.size()) {

    if(*pool.begin() < 0 || *pool.rbegin() >= size())
      throw funame << "orbit element(s) out of range";

    for(std::set<std::vector<int> >::const_iterator ot = orbit.begin(); ot != orbit.end(); ++ot)
      for(const_iterator it = ot->begin(); it != ot->end(); ++it) {
	const_iterator jt = it + 1;
	if(jt == ot->end())
	  jt = ot->begin();
	std::vector<int>::operator[](*it) = *jt;
      }
  }

  // update with identical permutations
  for(int i = 0; i < size(); ++i)
    if(pool.find(i) == pool.end())
      std::vector<int>::operator[](i) = i;
}

Permutation Permutation::operator* (const Permutation& p) const
{
  Exception::Base funame = "Permutation::operator*: ";

  if(p.size() != size())
    throw funame << "dimensions mismatch";
    
  std::vector<int> res(size());
  for(int i = 0; i < size(); ++i)
    res[i] = (*this)[p[i]];

  return Permutation(res);
}

Permutation Permutation::invert () const 
{
  std::vector<int> res(size());
  for(int i = 0; i < size(); ++i)
    res[(*this)[i]] = i;

  return Permutation(res);   
}

int Permutation::_compare (const Permutation& p) const
{ 
  Exception::Base funame = "Permutation::_compare: ";

  if(size() != p.size())
    throw funame << "sizes mismatch";

  std::vector<int>::const_iterator pit = p.begin();
  for(std::vector<int>::const_iterator it = begin(); it != std::vector<int>::end(); ++it, ++pit)
    if(*it < *pit)
      return -1;
    else if(*it > *pit)
      return  1;

  return 0;
}

void Permutation::operator++ ()
{
  const char funame [] = "Permutation::operator++: ";
  if(_end)
    return;

  std::set<int> pool;
  std::set<int>::const_iterator p;

  for(reverse_iterator i = rbegin(); i != rend(); ++i) {
    p = pool.insert(*i).first;

    if(++p != pool.end()) {
      *i = *p;
      pool.erase(p);
      std::set<int>::const_reverse_iterator  q = pool.rbegin();
      for(reverse_iterator j = rbegin(); j != i; ++j, ++q)
	*j = *q;
      return;
    }
  }

  std::set<int>::const_reverse_iterator  q = pool.rbegin();
  for(reverse_iterator j = rbegin(); j != rend(); ++j, ++q)
    *j = *q;

  _end = true;
}

std::set<std::vector<int> > Permutation::orbit () const
{
  std::set<std::vector<int> > res;
  
  int          itemp;
  std::set<int> done;

  for(const_iterator pit = begin();  pit != std::vector<int>::end(); ++pit)
    if(done.find(*pit) == done.end()) {

      std::vector<int> orb;

      itemp = *pit;
      do {
	orb.push_back(itemp);
	done.insert(itemp);
	itemp = (*this)[itemp];
      } while(itemp != *pit);

      if(orb.size() > 1)
	res.insert(orb);
    }

  return res;
}

std::vector<int> Permutation::orbit_max () const
{
  std::set<std::vector<int> > orb = orbit();

  std::vector<int> res;
  for(std::set<std::vector<int> >::const_iterator ot = orb.begin(); ot != orb.end(); ++ot)
    if(ot == orb.begin() || ot->size() > res.size())
      res = *ot;

  return res;    
}

std::set<std::vector<int> > Permutation::read_orbit(std::istream& from)
{
  Exception::Base funame = "Permutation::read_orbit: ";
  
  std::set<std::vector<int> > res;

  char next = IO::skip_space(from);
  switch(next) {
  case EOF:
    throw Exception::Eof(funame << "end of file");
    
  case '(':
    from.get(next);
    break;

  default:
    throw funame << "unknown opening symbol: " << next;
  }

  while(1)
    switch(next = IO::skip_space(from)) {
    case EOF:
      throw funame << "unexpected end of input";

    case '(':
      from.get(next);
      res.insert(read_cycle(from));
      break;

    case ')':
      from.get(next);
      if(!res.size())
	throw funame << "empty orbit";
      return res;

    default:
      throw funame << "unknown symbol: " << next;
    }
}

std::vector<int> Permutation::read_cycle(std::istream& from)
{
  Exception::Base funame = "Permutation::read_cycle: ";

  int itemp;
  std::vector<int> res;
  char next;

  while(1)
    switch(next = IO::skip_space(from)) {
    case EOF:
      throw funame << "unexpected end of input";

    case ')':
      from.get(next);
      if(!res.size())
	throw funame << "empty cycle";
      return res;
    
    default:
      if(!(from >> itemp))
	throw funame << ": corrupted; current char = " << next;
      res.push_back(itemp);
    }
}

std::set<Permutation> permutation_group (const std::set<Permutation>& base) 
{
  Exception::Base funame = "permutation_group: ";

  std::set<Permutation> res = base;

  bool btemp;

  if(!res.size() || !res.begin()->size())
    throw funame << "no permutation";

  // checking sizes
  for(std::set<Permutation>::const_iterator i = res.begin(); i != res.end(); ++i)
    if(i->size() != res.begin()->size())
      throw funame << "sizes mismatch";

  do {
    btemp = false;
    for(std::set<Permutation>::const_iterator i = res.begin(); i != res.end(); ++i) {
      for(std::set<Permutation>::const_iterator j = res.begin(); j != res.end(); ++j) {
	Permutation ij = *i * *j;
	if(res.find(ij) == res.end()) {
	  res.insert(ij);
	  btemp = true;
	  break;
	}
      }
      if(btemp)
	break;
    }
  } while(btemp);

  return res;
}

// full permutation group
std::set<Permutation> permutation_group (int n) 
{
  Exception::Base funame = "permutation_group: ";

  if (n <= 0)
    throw funame << "out of range";

  if (n == 1) {
    std::set<Permutation> res;
    res.insert(Permutation(1));
    return res;
  }

  std::set<Permutation> res;
  const std::set<Permutation> base = permutation_group(n - 1);
  for(std::set<Permutation>::const_iterator p = base.begin(); p != base.end(); ++p)
    for(int i = 0; i < n; ++i) {
      std::vector<int> perm = p->base();
      perm.insert(perm.begin() + i, n - 1);
      res.insert(Permutation(perm));
    }
  
  return res;
}

/************************************************************************************************
 ***************************** ITERATION OVER MULTIPLE PERMUTATIONS *****************************
 ************************************************************************************************/

MultiPerm::MultiPerm (const std::vector<int>& dim) throw(Error::General) : _end(false)
{
  const char funame [] = "MultiPerm::MultiPerm: ";

  if(!dim.size()) {
    std::cerr << funame << "no permutations";
    throw Error::Range();
  }
  
  for(std::vector<int>::const_iterator it = dim.begin(); it != dim.end(); ++it)
    if(*it > 0)
      std::vector<Permutation>::push_back(Permutation(*it));
    else {
      std::cerr << funame << "dimension out of range\n";
      throw Error::Range();
    }
}

void MultiPerm::operator++ ()
{
  _end = true;
  for(std::vector<Permutation>::reverse_iterator i = rbegin(); i != rend(); ++i) {
    ++(*i);
    if(!i->end()) {
      _end = false;
      break;
    }
    else
      *i = Permutation(i->size());
  }
}

