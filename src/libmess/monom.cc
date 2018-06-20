/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include "monom.hh"
#include "multindex.hh"

Monom::iterator::iterator (const std::vector<int>& val) throw(Error::General)
  : std::vector<int>(val), _end(false) 
{
  const char funame [] = "Monom::iterator::iterator: ";

  int itemp = 0;
  for(std::vector<int>::const_iterator it = begin(); it != std::vector<int>::end(); ++it)
    if(*it < 0) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }
    else
      itemp += *it;

  if(!itemp) {
    std::cerr << funame << "zero index\n";
    throw Error::Range();
  }
}

void Monom::iterator::operator++ () 
{
  int itemp;

  std::vector<int>::reverse_iterator fin = rend() - 1;
  for(std::vector<int>::reverse_iterator i = rbegin(); i != fin; ++i)
    if(*i) {
      itemp = *i;
      *i++ = 0;
      ++(*i);
      *rbegin() = itemp - 1;
      return;
    }
  
  _end = true;
}

Monom::Monom (int rn, int sz) throw(Error::General) : _rank(rn), _size(sz)
{
  const char funame [] = "Monom::Monom: ";

  if(rn < 1 || sz < 1) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  double dtemp;
  int    itemp;

  // iterator begin 
  std::vector<int> vtemp(_size, 0);
  *vtemp.rbegin() = _rank;
  _begin = iterator(vtemp);

  for(iterator it = _begin; !it.end(); ++it) {
    _multi_index_map[*it] = _index_multi_map.size();
    _index_multi_map.push_back(*it);
  }  

  //std::cout << funame << "expansion size = " << linear_size() << "\n";
}


Lapack::Matrix Monom::operator() (const Lapack::Matrix& a) const throw(Error::General)
{
  const char funame [] = "Monom::operator(): ";

  if(a.size() != size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

  int    itemp;
  double dtemp;

  // test multiindex
  /*
  itemp = 0;
  for(MultiIndex tj(rank(), size() - 1); !tj.end(); ++tj)
    ++itemp;
  std::cout << funame << "multiindex size = " << itemp << "\n";
  */

  Lapack::Matrix res(linear_size());
  res = 0.;

#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)      

  for(int li = 0; li < linear_size(); ++li) {
    // monomial power
    std::vector<int> mi = _index_multi_map[li];
    // tensor index
    std::vector<int> ti;
    for(int i = 0; i < size(); ++i)
      for(int j = 0; j < mi[i]; ++j)	
	ti.push_back(i);

    for(MultiIndex tj(rank(), size() - 1); !tj.end(); ++tj) {
      std::vector<int> mj(size(), 0);
      for(int i = 0; i < rank(); ++i)
	++mj[tj[i]];
      
      std::map<std::vector<int>, int>::const_iterator mit = _multi_index_map.find(mj);
      if(mit == _multi_index_map.end()) {
	std::cerr << funame << "monomial index is not in the map\n";
	throw Error::Range();
      }
      int lj = mit->second;

      dtemp = 1.;
      for(int k = 0; k < rank(); ++k)
	dtemp *= a(ti[k], tj[k]);

      res(li, lj) += dtemp;
    }
  }

  return res;
}
