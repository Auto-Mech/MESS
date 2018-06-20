std::set<int> Graph::Expansion::_low_freq_set (double temperature, std::vector<double>& tanh_factor) const
{
  const char funame [] = "Graph::Expansion::_low_freq_set: ";

  std::set<int> res;

  if(temperature > 0.) {
    //
    tanh_factor.resize(_red_freq.size());

    for(int i = 0; i < _red_freq.size(); ++i)
      //
      if(_red_freq[i] > 0.) {
	//
	tanh_factor[i] = std::tanh(_red_freq[i] / temperature / 2.);
      
	if(_red_freq[i] / temperature < low_freq_thresh)
	  //
	  res.insert(i);
      }
      else {
	//
	tanh_factor[i] = std::tan(-_red_freq[i] / temperature / 2.);
	
	if(_red_freq[i] / temperature > -low_freq_thresh)
	  //
	  res.insert(i);
      }
  }
  // zero temperature
  //
  else {
    //
    tanh_factor.clear();
  }

  return res;
}

void Graph::Expansion::_set_frequencies (std::vector<double> freq)
{
  const char funame [] = "Graph::Expansion::_set_frequencies: ";
  
  static const double min_freq = Phys_const::incm;
  
  double dtemp;
  int    itemp;

  _red_freq_map.clear();

  // removing very low frequencies
  //
  for(int f = 0; f < freq.size(); ++f) {
    //
    if(freq[f] >= 0. && freq[f] < min_freq) {
      //
      freq[f] =  min_freq;
    }
    else if(freq[f] <= 0. && freq[f] > -min_freq) {
      //
      freq[f] = -min_freq;
    }
  }
  
  for(int f = 0; f < freq.size(); ++f) {
    //
    itemp = 1;
    //
    for(int i = 0; i < _red_freq_map.size(); ++i) {
      //
      dtemp = freq[f] / freq[*_red_freq_map[i].begin()];

      if(1. - freq_tol < dtemp && dtemp < 1. + freq_tol) {
	//
	_red_freq_map[i].insert(f);

	itemp = 0;
	
	break;
      }
    }
    
    if(itemp) {
      //
      _red_freq_map.push_back(std::set<int>());
      //
      _red_freq_map.back().insert(f);
    }
  }

  _red_freq.resize(_red_freq_map.size());
  //
  for(int i = 0; i < _red_freq_map.size(); ++i) {
    //
    dtemp = 0.;
    
    itemp =  1;
    //
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it) {
      //
      if(it == _red_freq_map[i].begin() && freq[*it] < 0.)
	//
	itemp = 0;

      if(itemp && freq[*it] < 0. || !itemp && freq[*it] > 0.) {
	//
	ErrOut err_out;
	//
	err_out << funame << "frequencies of different signs in one group";
      }
      
      if(itemp) {
	//
	dtemp += std::log(freq[*it]);
      }
      else {
	//
	dtemp += std::log(-freq[*it]);
      }
    }
    
    dtemp /= (double)_red_freq_map[i].size();

    if(itemp) {
      //
      _red_freq[i] =  std::exp(dtemp);
    }
    else {
      //
      _red_freq[i] = -std::exp(dtemp);
    }
  }

  _red_freq_index.resize(freq.size());
  //
  for(int i = 0; i < _red_freq_map.size(); ++i) {
    //
    for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it)
      //
      _red_freq_index[*it] = i;
  }
  
  if(!IO::mpi_rank) {
    //
    IO::log << IO::log_offset << "logarithmic frequency tolerance  = " << freq_tol << "\n";

    IO::log << IO::log_offset << "reduced frequencies, 1/cm:";
  
    for(int i = 0; i < _red_freq.size(); ++i) {
      //
      IO::log << "   " << _red_freq[i] / Phys_const::incm;
    }
  
    IO::log << "\n";

    IO::log << IO::log_offset << "reduced frequencies groups: ";
    //
    for(int i = 0; i < _red_freq_map.size(); ++i) {
      //
      for(std::set<int>::const_iterator it = _red_freq_map[i].begin(); it != _red_freq_map[i].end(); ++it) {
	//
	if(it != _red_freq_map[i].begin())
	  IO::log << ", ";
	else
	  IO::log << "(";
      
	IO::log << *it;
      }
    
      IO::log << ")";
    }
  
    IO::log << "\n\n";
  }
}
