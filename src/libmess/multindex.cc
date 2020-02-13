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
#include <cmath>

/********************************************************************************************
 ******************************* GROWING SYMMETRIC TENSOR INDEX ***************************** 
 ********************************************************************************************/

void SymIndex::operator++ () 
{
  int itemp;

  for(int i = 0; i < size(); ++i) {
    //
    itemp = ++(*this)[i];
    //
    if(itemp < _range) {
      //
      for(int j = 0; j < i; ++j) {
	//
	(*this)[j] = itemp;
      }
      //
      return;
    }
  }

  itemp = size() + 1;
  //
  clear();
  //
  resize(itemp, 0);
}

SymIndex::operator std::multiset<int> () const 
{
  std::multiset<int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    res.insert(*it);
  }

  return res;
}

SymIndex::operator std::map<int, int> () const 
{
  std::map<int, int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    ++res[*it];
  }

  return res;
}

/*********************************************************************************************
 ******************************* FIXED SIZE GENERIC TENSOR INDEX ***************************** 
 *********************************************************************************************/

void GenIndex::operator++ () 
{
  const char funame [] = "GenIndex::operator++: ";

  if(_fin) {
    ErrOut err_out;
    err_out << funame << "already finished";
  }
  
  int itemp;

  for(int i = 0; i < size(); ++i) {
    //
    itemp = ++(*this)[i];
    //
    if(itemp < _range) {
      //
      for(int j = 0; j < i; ++j) {
	//
	(*this)[j] = 0;
      }
      //
      return;
    }
  }
  
  _fin = true;
}

GenIndex::operator std::multiset<int> () const 
{
  std::multiset<int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    res.insert(*it);
  }

  return res;
}

GenIndex::operator std::map<int, int> () const 
{
  std::map<int, int> res;
  //
  for(const_iterator it = begin(); it != end(); ++it) {
    //
    ++res[*it];
  }

  return res;
}

/***********************************************************************************************
***************************** NUMERICAL DIFFERENTIATION MAPPER *********************************
************************************************************************************************/

int NumDer::order (const der_t& der)
{
  int res = 0;
  
  for(der_t::const_iterator dit = der.begin(); dit != der.end(); ++dit)
    //
    res += dit->second;

  return res;
}

void NumDer::_assert (const der_t& der) const
{
  const char funame [] = "NumDer::_assert: ";
   
  for(der_t::const_iterator dit = der.begin(); dit != der.end(); ++dit)
    //
    if(dit->first < 0 || dit->first >= size() || dit->second <= 0) {
      //
      ErrOut err_out;

      err_out << funame << "derivative index and/or order out of range:" << dit->first << ", " << dit->second;
    }
}

// numerical derivative from the function-on-the-grid database
//
double NumDer::operator() (const der_t& der, const fdata_t& fdata, const std::vector<double>& step) const
{
  const char funame [] = "NumDer::operator(): ";

  double dtemp;
  
  int    itemp;

  _check();

  _assert(der);

  if(step.size() != size()) {
    //
    ErrOut err_out;

    err_out << funame << "wrong differentiation steps number: " << step.size();
  }

  for(int i = 0; i < step.size(); ++i)
    //
    if(step[i] <= 0.) {
      //
      ErrOut err_out;

      err_out << funame << i + 1 << "-th step our of range: " << step[i];
    }
  
  cmap_t cmap = convert(der);

  double res = 0.;
  //
  for(cmap_t::const_iterator cit = cmap.begin(); cit != cmap.end(); ++cit) {
    //
    fdata_t::const_iterator fit = fdata.find(cit->first);

    if(fit == fdata.end()) {
      //
      ErrOut err_out;
	 
      err_out << funame << "cannot find the grid point (";

      for(int i = 0; i < cit->first.size(); ++i) {
	//
	if(i)
	  //
	  err_out << ", ";
	  
	err_out << cit->first[i];
      }
	
      err_out <<  ") in the function-on-the-grid database";
    }
    else
      //
      res += (double)cit->second * fit->second;
  }

  // normalization
  //
  for(der_t::const_iterator dit = der.begin(); dit != der.end(); ++dit) {
    //
    res /= std::pow(step[dit->first], (double)dit->second);
    
    // correction for odd derivatives
    //
    if(dit->second % 2)
      //
      res /= 2.;
  }

  return res;
}

// converts derivative signature into the displaced configurations mapping
//
NumDer::cmap_t NumDer::convert (const der_t& der) const
{
  const char funame [] = "NumDer::operator(): ";
  
  _check();

  cmap_t cmap;

  if(!der.size()) {
    //
    cmap[std::vector<int>(_size)] = 1;
    //
    return cmap;
  }

  if(der.begin()->first < 0 || der.begin()->first >= _size) {
    //
    ErrOut err_out;
    //
    err_out << funame << "wrong derivative index: " << der.begin()->first;
  }
    
  if(der.begin()->second <= 0) {
    //
    ErrOut err_out;
    //
    err_out << funame << "wrong derivative order: " << der.begin()->second;
  }

  der_t new_der = der;

  if(der.begin()->second <= 2) {
    //
    new_der.erase(new_der.begin());
  }
  else {
    //
    new_der.begin()->second -= 2;
  }

  cmap_t new_cmap = convert(new_der);

  for(cmap_t::const_iterator cmit = new_cmap.begin(); cmit != new_cmap.end(); ++cmit) {
    //
    std::vector<int> conf = cmit->first;
	
    // first derivative
    //
    if(der.begin()->second == 1) {
      //
      conf[der.begin()->first] += 1;
      //
      cmap[conf] += cmit->second;

      conf[der.begin()->first] -= 2;
      //
      cmap[conf] -= cmit->second;
    }
    // second derivative
    //
    else {
      //
      cmap[conf] -= 2 * cmit->second;

      conf[der.begin()->first] += 1;
      //
      cmap[conf] += cmit->second;
	
      conf[der.begin()->first] -= 2;
      //
      cmap[conf] += cmit->second;
    }
  }

  // clean configuration map
  //
  cmap_t res;
  //
  for(cmap_t::const_iterator cmit = cmap.begin(); cmit != cmap.end(); ++cmit) {
    //
    if(cmit->second) {
      //
      res.insert(*cmit);
    }
  }
  return res;
}

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


