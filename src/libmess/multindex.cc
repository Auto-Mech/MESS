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

#include "multindex.hh"
#include <iostream>

/************************************************************************************************
 ********************************** MULTIDIMENSIONAL INDEX **************************************
 ************************************************************************************************/

void MultiIndex::set (const std::vector<int>& gr, const std::vector<int>& gs) 
{
  const char funame [] = "MultiIndex::set: ";

  int itemp;

  _clear();

  if(!gr.size() && !gs.size())
    return;
  
  if(gr.size() != gs.size()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

  for(int g = 0; g < gr.size(); g++)
    if(gr[g] < 1 || gs[g] < 0) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }

  _end = false;

  _group_rank = gr;
  _group_size = gs;
  _group_end.resize(gr.size());

  itemp = 0;
  for(int g = 0; g < gr.size(); g++) {
    itemp += gr[g];
    _group_end[g] = itemp;
  }

  resize(itemp, 0);
    
  int group = 0;
  _group_map.resize(size());

  for(int i = 0; i < size(); ++i) {
    if(i == _group_end[group])
      ++group;

    itemp = size() - i - 1;
    _group_map[itemp] = group;
  }
}

void MultiIndex::operator++ () 
{
  const char funame [] = "MultiIndex::operator++: ";

  int itemp;

  for(reverse_iterator i = rbegin(); i != rend(); ++i) {
    int group = _group_map[i - rbegin()];
    if(*i < _group_size[group]) {
      ++(*i);

      reverse_iterator g = rbegin() + size() - _group_end[group];
      for(reverse_iterator j = rbegin(); j != g; ++j)
	*j = 0;
      for(reverse_iterator j = g; j != i; ++j)
	*j = *i;
      
      return;
    }
  }

  _end = true;
}

std::vector<std::multiset<int> >  MultiIndex::convert_multiset () const
{
  int itemp;

  std::vector<std::multiset<int> > res;

  res.resize(_group_rank.size());
  itemp = 0;
  for(int g = 0; g < _group_rank.size(); itemp += _group_rank[g++])
    for(int i = 0; i < _group_rank[g]; ++i)
      res[g].insert((*this)[itemp + i]);
    
  return res;
}

int MultiIndex::cumulative_index () const
{
  const char funame [] = "MultiIndex::cumulative_index: ";

  int res = 0;
  for(const_iterator it = begin(); it != std::vector<int>::end(); ++it)
    res += *it;
  return res;

}

/************************************************************************************************
 **************************** MULTI-DIMENSIONAL INDEX CONVERTER *********************************
 ************************************************************************************************/

void MultiIndexConvert::resize (const std::vector<int>& s) 
{
  static const char funame [] = "MultiIndexConvert::resize: ";

  std::vector<int>::operator=(s);

  if(!rank()) {
    _linear_size = 0;
    return;
  }
  
  _linear_size = 1;
  for(const_iterator i = begin(); i != end(); _linear_size *= (long)*i++)
    if(*i < 1) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }  
}

std::vector<int> MultiIndexConvert::operator() (long lindex) const 
{
  static const char funame [] = "MultiIndexConvert::operator(): ";

  if(lindex < 0 || lindex >= _linear_size) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  std::vector<int> res(rank());
  iterator j = res.begin();
  for(const_iterator i = begin(); i != end(); lindex /= (long)*i++, ++j)
    *j = lindex % (long)*i;

  return res;
}

long MultiIndexConvert::operator() (const std::vector<int>& v) const 
{
  static const char funame [] = "MultiIndexConvert::operator(): ";

  if(v.size() != rank()) {
    std::cerr << funame << "sizes mismatch\n";
    throw Error::Range();
  }

  long res = 0;
  long dim = 1;

  const_iterator j = v.begin();
  for(const_iterator i = begin(); i != end(); dim *= (long)*i++, ++j) {
    if(*j < 0 || *j >= *i) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }

    res += (long)*j * dim;
  }

  return res;
}

long MultiIndexConvert::conjugate (long lindex) const
{
  static const char funame [] = "MultiIndexConvert::conjugate: ";
  
  std::vector<int> v = (*this)(lindex);
  iterator j = v.begin();
  for(const_iterator i = begin(); i != end(); ++i, ++j)
    if(*j)
      *j = *i - *j;
  
  return (*this)(v);
}


