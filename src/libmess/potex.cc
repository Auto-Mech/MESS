#include "potex.hh"

/********************************************************************************************
 ************************************ POTENTIAL EXPANSION ***********************************
 ********************************************************************************************/

void PotentialExpansion::read_from_file (const std::vector<double>& freq, std::istream& from)
{
  const std::string funame [] = "PotentialExpansion::read_from_file: ";

  int    itemp;
  double dtemp;
  
  if(!freq.size()) {
    ErrOut err_out;
    err_out << funame << "no frequencies";
  }

  if(index_range() && index_range() != freq.size()) {
    ErrOut err_out;
    err_out << funame << "index range changed";
  }

  _index_range = freq.size();
  
  std::vector<double> freq_sqrt(freq.size());
  
  for(int i = 0; i < freq.size(); ++i)
    if(freq[i] == 0.){
      ErrOut err_out;
      err_out << funame << "zero frequency";
    }
    else if(freq[i] < 0.)
      freq_sqrt[i] = std::sqrt(-freq[i]);
    else
      freq_sqrt[i] = std::sqrt(freq[i]);
  
  while(from) {
    IO::LineInput lin(from);
    
    IO::String stemp;
    
    std::vector<IO::String> string_value;
    while(lin >> stemp)
      string_value.push_back(stemp);

    if(!string_value.size())
      continue;
    
    if(string_value[0] == IO::end_key())
      return;

    if(string_value.size() < 4) {
      std::cerr << funame << "not enough values (linear and quadratic terms not allowed)\n";
      throw Error::Input();
    }

    double value = double(string_value.back()) * Phys_const::incm;

    const int imax = string_value.size() - 1;
    
    index_t index;
    for(int i = 0; i < imax; ++i) {
      itemp = int(string_value[i]) - 1;

      if(itemp < 0 || itemp >= index_range()) {
	ErrOut err_out;
	err_out << funame << "index out of range: " << itemp;
      }

      index.insert((int_t)itemp);
      value *= freq_sqrt[itemp];
    }

    assign(index, value);
  }

  if(!from) {
    ErrOut err_out;
    err_out << funame << "corrupted";
  }
}

void PotentialExpansion::_check_index (const index_t& i) const
{
  const char funame [] = "PotentialExpansion::_check_index: ";

  if(index_range() <= 0) {
    ErrOut err_out;
    err_out << funame << "index range out of range: " << index_range();
  }
  
  if(!i.size())
    return;
    
  if(*i.begin() < 0 || *i.rbegin() >= index_range()) {
    ErrOut err_out;
    err_out << funame << "index out of range:";
    for(index_t::const_iterator it = i.begin(); it != i.end(); ++it)
      err_out << " " << (int)*it;
  }
}
  
std::set<int> PotentialExpansion::rank_pool () const
{
  std::set<int> res;
  for(const_iterator pit = begin(); pit != end(); ++pit)
    res.insert(pit->first.size());

  return res;
}
