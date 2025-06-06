void MasterEquation::branching_fraction (double temperature, double pressure, const std::set<int>& wg,
					 //
					 std::map<std::set<int>, std::map<int, double> >& bfdb)
{
  const char funame [] = "MasterEquation::branching_fraction: ";

  int    itemp;
  double dtemp;
  
  if(!wg.size()) {
    //
  }

  std::map<std::set<int>, std::map<int, double> >::const_iterator bfi = bfdb.find(wg);
  
  if(bfi != bfdb.end())
    //
    return;
  
  std::map<int, double> res;

  if(wg.size() == 1) {
    //
    res[*wg.begin()] = 1.;

    bfdb[wg] = res;
    
    return;
  }

  IO::Marker funame_marker(funame);
  
  // chemical graph
  //
  IO::log << IO::log_offset << "the chemical group:  ";
			       
  for(std::set<int>::const_iterator w = wg.begin(); w != wg.end(); ++w) {
    //
    if(w != wg.begin())
      //
      IO::log << "+";

    IO::log << Model::well(*w).name();
  }

  IO::log << "\n";

    
  Model::ChemGraph cg(wg);

  std::list<Model::ChemGraph> gl;

  // splitting barrier
  //
  const int sb = cg.split(temperature, &gl);

  const double be =  Model::inner_barrier(sb).dist_ener_max(temperature);

  IO::log << IO::log_offset << "well-splitting barrier " << Model::inner_barrier(sb).name()
    //
	  << ": distribution maximum energy = " << std::ceil(be / Phys_const::kcal * 10.) / 10. << " kcal/mol\n";
  
  double er = be + excess_energy_over_temperature * temperature;

  double ec = be - energy_cutoff_over_temperature * temperature;

  ReactiveComplex rc(temperature, pressure, er, ec, cg, ReactiveComplex::BRANCHING);

  itemp = 1;

  ReactiveComplex::RateData rate_data;
  
  try {
    //
    rate_data = rc.well_reduction_method();
    
    if(!rate_data.well_partition.size()) {
      //
      std::cerr << funame << "no bound groups\n";

      throw Error::Run();
    }

    itemp = rate_data.well_partition.size() + rate_data.bimolecular_group.size();
  }
  catch(Error::Lapack) {
    //
    IO::log << "\n" << IO::log_offset << "WARNING: lapack diagonalization failed, will use instead state densities for branching\n\n";
  }
  catch(BranchingException) {
    //
    IO::log << IO::log_offset << "WARNING: branching calculation failed, will use state densities instead\n\n";
  }
  
  // setting branching fraction manually
  //
  if(itemp == 1) {
    //
    IO::log << IO::log_offset << "getting the branching fraction for the two groups from their state densities by\n";

    IO::log << IO::log_offset << "splitting the group into two: ";

    for(std::set<int>::const_iterator w = gl.front().well_set.begin(); w != gl.front().well_set.end(); ++w) {
    //
    if(w != gl.front().well_set.begin())
      //
      IO::log << "+";

    IO::log << Model::well(*w).name();
  }

  IO::log << " and ";
    
  for(std::set<int>::const_iterator w = gl.back().well_set.begin(); w != gl.back().well_set.end(); ++w) {
    //
    if(w != gl.back().well_set.begin())
      //
      IO::log << "+";
    
    IO::log << Model::well(*w).name();
  }

  IO::log << ", at the " << Model::inner_barrier(sb).name() << " barrier distribution maximum energy, "
      //
	  << std::ceil(be / Phys_const::kcal * 10.) / 10. <<" kcal/mol\n";
	       
    int count = 0;

    double nfac = 0.;

    std::vector<double> vtemp(gl.size());
    
    for(std::list<Model::ChemGraph>::const_iterator li = gl.begin(); li != gl.end(); ++li, ++count) {
      //
      dtemp = Model::state_density(li->well_set, be);

      vtemp[count] = dtemp;

      nfac += dtemp;
    }

    if(nfac == 0.) {
      //
      IO::log << IO::log_offset << "WARNING: the group has zero state density\n";
      
      for(int i = 0; i < gl.size(); ++i)
	//
	vtemp[i] = 1. / vtemp.size();
    }
    else {
      //
      for(int i = 0; i < gl.size(); ++i)
	//
	vtemp[i] /= nfac;
    }

    IO::log << IO::log_offset << "groups branching fractions  based on their state densities:  "
      //
	    << vtemp[0] << " and " << vtemp[1] << "\n\n";

    count = 0;
      
    for(std::list<Model::ChemGraph>::const_iterator li = gl.begin(); li != gl.end(); ++li, ++count) {
      //
      branching_fraction(temperature, pressure, li->well_set, bfdb);

      std::map<int, double>& bf = bfdb[li->well_set];
      
      for(std::map<int, double>::const_iterator w = bf.begin(); w != bf.end(); ++w)
	//
	res[w->first] = w->second * vtemp[count];
    }

    IO::log << IO::log_offset << "branching fractions for the ";

    for(std::set<int>::const_iterator w = wg.begin(); w != wg.end(); ++w) {
      //
      if(w != wg.begin())
	//
	IO::log << "+";

      IO::log << Model::well(*w).name();
    }

    IO::log << " group:";
      
    for(std::map<int, double>::const_iterator wi = res.begin(); wi!= res.end(); ++wi)
      //
      IO::log << "   " << wi->second << "/" << Model::well(wi->first).name();

    IO::log << "\n";

    bfdb[wg] = res;
    
    return;
  }

  for(int g = 0; g < rate_data.well_partition.size(); ++g) {
    //
    branching_fraction(temperature, pressure, rate_data.well_partition[g], bfdb);

    std::map<int, double>& bf = bfdb[rate_data.well_partition[g]];

    for(std::map<int, double>::const_iterator w = bf.begin(); w != bf.end(); ++w)
      //
      res[w->first] = w->second * rate_data.branching_fraction[g];
  }
  
  for(std::set<int>::const_iterator w = rate_data.bimolecular_group.begin(); w != rate_data.bimolecular_group.end(); ++w)
    //
    res[*w] = 0.;

  IO::log << IO::log_offset << "branching fractions for the ";

  for(std::set<int>::const_iterator w = wg.begin(); w != wg.end(); ++w) {
    //
    if(w != wg.begin())
      //
      IO::log << "+";

    IO::log << Model::well(*w).name();
  }

  IO::log << " group:";

  for(std::map<int, double>::const_iterator wi = res.begin(); wi != res.end(); ++wi)
    //
    IO::log << "   " << wi->second << "/" << Model::well(wi->first).name();

  IO::log << "\n";
  
  bfdb[wg] = res;
  
  return;
}
