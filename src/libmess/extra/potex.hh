/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef POTEX_HH
#define POTEX_HH

#include<vector>
#include<set>
#include<map>
#include<iostream>

/******************************************************************************************
 ***************************** UNHARMONIC POTENTIAL EXPANSION *****************************
 ******************************************************************************************/
  
typedef short _potex_int_t;

class PotentialExpansion : private std::map<std::multiset<_potex_int_t>, double> {
public:
  typedef _potex_int_t               int_t;
  typedef std::multiset<int_t>     index_t;
  typedef std::map<index_t, double> base_t;
    
private:
  template<typename T>
  static index_t _convert (const std::multiset<T>&);
    
  template<typename T> 
  static index_t _convert (const std::vector<T>&);

  int_t _index_range;

  void _check_index(const index_t&) const;

public:
  class NotExist {};
  class    Exist {};

  int size          () const { return base_t::size(); }

  int_t index_range () const { return _index_range; }
    
  std::set<int> rank_pool () const;
    
  void read_from_file (const std::vector<double>& freq, std::istream& from);

  PotentialExpansion () : _index_range(0) {}

  template<typename T>
  PotentialExpansion (int_t, const std::map<std::multiset<T>, double>&);
    
  PotentialExpansion (const std::vector<double>& freq, std::istream& from) : _index_range(0) { read_from_file (freq, from); }

  template<typename T>
  bool find (const std::vector<T>&   i)   const { return find(_convert(i)); }

  template<typename T>
  bool find (const std::multiset<T>& i)   const { return find(_convert(i)); }

  template<typename T>
  void assign (const std::vector<T>&   in, double val) throw(Exist) { assign(_convert(in), val); }
    
  template<typename T>
  void assign (const std::multiset<T>& in, double val) throw(Exist) { assign(_convert(in), val); }

  template<typename T>
  void erase (const std::vector<T>&   in) { erase(_convert(in)); }
    
  template<typename T>
  void erase (const std::multiset<T>& in) { erase(_convert(in)); }

  template<typename T>
  double operator() (const std::vector<T>&   in) const throw(NotExist) { return (*this)(_convert(in)); }

  template<typename T>
  double operator() (const std::multiset<T>& in) const throw(NotExist) { return (*this)(_convert(in)); }

  template<typename T>
  double operator() (int ...)  const throw(NotExist);
};

template<typename T> 
PotentialExpansion::index_t PotentialExpansion::_convert (const std::vector<T>& in)
{
  index_t res;
  for(typename std::vector<T>::const_iterator it = in.begin(); it != in.end(); ++it)
    res.insert((int_t)*it);

  return res;
}
    
template<typename T> 
PotentialExpansion::index_t PotentialExpansion::_convert (const std::multiset<T>& in)
{
  index_t res;
  for(typename std::multiset<T>::const_iterator it = in.begin(); it != in.end(); ++it)
    res.insert((int_t)*it);

  return res;
}
    
template<> 
inline PotentialExpansion::index_t PotentialExpansion::_convert (const index_t& in) { return in; }

template<typename T>
PotentialExpansion::PotentialExpansion (int_t r, const std::map<std::multiset<T>, double>& d)
  : _index_range(r)
{
  const char funame [] = "PotentialExpansion::PotentialExpansion: ";
    
  if(index_range() <= 0) {
    ErrOut err_out;
    err_out << funame << "index range out of range: " << index_range();
  }
    
  for(typename std::map<std::multiset<T>, double>::const_iterator mit = d.begin(); mit != d.end(); ++mit) {
    index_t in = _convert(mit->first);

    _check_index(in);
      
    (*this)[in] = mit->second;
  }
}

template<>
inline PotentialExpansion::PotentialExpansion (int_t r, const base_t& d)
  : _index_range(r), base_t(d)
{
  const char funame [] = "PotentialExpansion::PotentialExpansion: ";

  if(index_range() <= 0) {
    ErrOut err_out;
    err_out << funame << "index range out of range: " << index_range();
  }
    
  for(const_iterator it = begin(); it != end(); ++it)
    _check_index(it->first);
}

template<typename T>
double PotentialExpansion::operator() (int r ...) const throw(NotExist)
{
  const std::string funame [] = "PotentialExpansion::operator(): ";

  va_list ap;
  va_start(ap, r);

  index_t in;
  for(int i = 0; i < r; ++i)
    in.insert((int_t)va_arg(ap, T));

  va_end(ap);

  _check_index(in);
    
  return (*this)(in);
}

template<>
inline bool PotentialExpansion::find (const index_t& in) const
{
  _check_index(in);
    
  return base_t::find(in) != end();
}
    
template<>
inline void PotentialExpansion::assign (const index_t& in, double val) throw(Exist)
{
  _check_index(in);
    
  if(find(in))
    throw Exist();

  (*this)[in] = val;
}

template<>
inline void PotentialExpansion::erase (const index_t& in)
{
  _check_index(in);
    
  base_t::erase(in);
}

template<>
inline double PotentialExpansion::operator() (const index_t& in) const throw(NotExist)
{
  _check_index(in);
    
  const_iterator it = base_t::find(in);
  
  if(it != end())
    return it->second;
  else
    throw NotExist();
}

#endif  
