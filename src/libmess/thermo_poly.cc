




void MuratPot::read_potex (const std::vector<double>& freq, std::istream& from, std::map<std::multiset<int>, double>& potex)
{
  const std::string funame [] = "MuratPot::read_potex: ";

  int    itemp;
  double dtemp;
  
  if(!freq.size()) {
    //
    ErrOut err_out;
    
    err_out << funame << "no frequencies";
  }

  std::vector<double> freq_sqrt(freq.size());
  
  for(int i = 0; i < freq.size(); ++i) {
    //
    if(freq[i] <= 0.) {
      //
      ErrOut err_out;

      err_out << funame << "negative frequency[1/cm]: " << freq[i] / Phys_const::incm;
    }
    else {
      //
      freq_sqrt[i] = std::sqrt(freq[i]);
    }
  }
  
  while(from) {
    //
    IO::LineInput lin(from);
    
    IO::String stemp;
    
    std::vector<IO::String> string_value;
    //
    while(lin >> stemp) {
      //
      string_value.push_back(stemp);
    }

    if(!string_value.size())
      //
      continue;
    
    if(string_value[0] == IO::end_key())
      //
      return;

    if(string_value.size() < 4) {
      //
      ErrOut err_out;

      err_out << funame << "not enough values (linear and quadratic terms not allowed)";
    }

    double value = double(string_value.back()) * Phys_const::incm;

    const int imax = string_value.size() - 1;
    
    std::multiset<int> index;
    //
    for(int i = 0; i < imax; ++i) {
      //
      itemp = int(string_value[i]) - 1;

      if(itemp < 0 || itemp >= freq.size()) {
	//
	ErrOut err_out;

	err_out << funame << "index out of range: " << itemp;
      }

      index.insert(itemp);

      value *= freq_sqrt[itemp];
    }

    if(potex.find(index) != potex.end()) {
      //
      ErrOut err_out;

      err_out << funame << "term already in the expansion";
    }
    
    potex[index] = value;
  }

  if(!from) {
    //
    ErrOut err_out;

    err_out << funame << "corrupted";
  }
}

