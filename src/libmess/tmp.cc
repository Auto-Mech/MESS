MasterEquation::ReactiveComplex::RateData MasterEquation::ReactiveComplex::well_reduction_method (int flags) const
{
  const char funame [] = "MasterEquation::ReactiveComplex::well_reduction_method: ";

  using Model::operator<<;

  if(!isinit) {
    //
    std::cerr << funame << "reactive complex has not been initialized\n";
    
    throw Error::Init();
  }

  if(!kinetic_basis.size()) {
    //
    std::cerr << funame << "kinetic matrices have not been initialized\n";
    
    throw Error::Init();
  }
    
  IO::Marker funame_marker(funame);
  
  int                  itemp;
  double               dtemp;
  bool                 btemp;
  std::string          stemp;
  float_t            dd_temp;
  
  RateData rate_data;

  rate_data._bim_index_map = _bim_index_map;
  
  IO::log << IO::log_offset << "collision frequencies[1/sec] / well: ";

  for(int w = 0; w < well_size(); ++w)
    //
    IO::log << "   " << well(0).collision_frequency(pressure) / Phys_const::herz << "/" << well_model(w).name();

  IO::log << "\n";
  
  if(global_size < well_size() + 1) {
    //
    if(mode() == BRANCHING) {
      //
      IO::log << IO::log_offset << "WARNING: there are not enough kinetically active states for branching analysis\n";

      throw BranchingException();
    }

    group_t owg;
    
    for(int i = 0; i < well_size(); ++i)
      //
      owg.insert(index_to_well(i));

    rate_data.bimolecular_group = owg;

    if(mode() == LUMPING) {
      //
      return rate_data;
    }
    
    if(!bimolecular_size()) {
      //
      std::cerr << "no kinetically active states & no bimolecular channels... Hm...\n";

      throw Error::Logic();
    }

    // bimolecular-to-bimolecular rates
    //
    LAPACK::SymmetricMatrix mtemp(bimolecular_size());
    
    mtemp = 0.;
    
    for(int e = 0; e < ener_index_max; ++e) {
      //
      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	//
	LAPACK::Vector vtemp(bimolecular_size(), 0.);
    
	for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	  //
	  const int w = bit->first.first;
      
	  const int p = bit->first.second;

	  if(e >= bit->second.size())
	    //
	    continue;

	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "well size assiciated with the barrier smaller than the barrier size\n";

	    throw Error::Logic();
	  }
	  
	  vtemp[p] += kinetic_basis[e].eigenvector(i->second, k)
	    //
	    * bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e] / 2. / M_PI;
	}

	for(int p = 0; p < bimolecular_size(); ++p)
	  //
	  for(int q = 0; q < p; ++q)
	    //
	    mtemp(p, q) += vtemp[p] * vtemp[q] / kinetic_basis[e].eigenvalue[k];
      }
    }

    // normalization
    //
    mtemp *= energy_step / std::exp(energy_reference / temperature);

    rate_data.bb_rate.resize(bimolecular_size());

    for(int p = 0; p < bimolecular_size(); ++p)
      //
      for(int q = 0; q < p; ++q)
	//
	rate_data.bb_rate(p, q) = CONVERT_DD(mtemp(p, q));
    
    return rate_data;
  }
  
  /***********************************************************************
   ************************* GLOBAL MATRICES *****************************
   ***********************************************************************/

  LAPACK::Vector eigenval;
  
  LAPACK::Matrix global_eigen(global_size);

  {
    IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);

    eigenval = kin_mat.eigenvalues(&global_eigen);
  }

  const float_t& relax_eval_min = eigenval[well_size()];
  
  const float_t& relax_eval_max = eigenval.back();

  IO::log << IO::log_offset << "minimal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_min / well(0).collision_frequency(pressure) << "\n";
  
  IO::log << IO::log_offset << "maximal relaxation eigenvalue over collision frequency = "
    //
	  << relax_eval_max / well(0).collision_frequency(pressure) << "\n";

  // eigenvectors projection on thermal subspace;
  //
  Lapack::Matrix eigen_pop(global_size, well_size());

  eigen_pop = 0.;

#pragma omp parallel for default(shared) private(dd_temp, itemp) schedule(static)
  
  for(int l = 0; l < global_size; ++l) {
    //
    for(int w = 0; w < well_size(); ++w) {
      //
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	if(i == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "eigen_pop: well is not in the kinetic basis space\n";

	  throw Error::Logic();
	}

	dd_temp = 0.;
	
	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dd_temp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	dd_temp *= well(w).boltzman_sqrt[e] / well(w).weight_sqrt;

	eigen_pop(l, w) += CONVERT_DD(dd_temp);
      }
      //
    }//
    //
  }//

  /************************************* EIGENVECTOR OUTPUT *****************************************/

#ifdef DEBUG
  
  IO::log << IO::log_offset << "orthogonality test:";

  for(int l = 0; l < well_size(); ++l) {
    //
    float_t res = 0.;
    
    for(int w = 0; w < well_size(); ++w) {
      //
      for(int e = 0; e < well(w).size(); ++e) {
	//
	std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	if(i == kinetic_basis[e].well_index_map.end()) {
	  //
	  std::cerr << funame << "eigen_pop: well is not in the kinetic basis space\n";

	  throw Error::Logic();
	}


	dd_temp = 0.;
	
	for(int k = 0; k < kinetic_basis[e].active_size; ++k)
	  //
	  dd_temp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);

	res += dd_temp * dd_temp;
      }
    }

    IO::log << "   " << res;
  }

  IO::log << "\n";

#endif

  std::vector<double> relaxation_projection(well_size());
  
  for(int l = 0; l < well_size(); ++l)
    //
    relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));

 
  IO::log << IO::log_offset << "eigenstate populations:\n"
    //
	  << IO::log_offset
    //
	  << std::setw(5)  << "L"
    //
	  << std::setw(10) << "*R"
    //
	  << std::setw(10) << "*P";
  
  for(int w = 0; w < well_size(); ++w)
    //
    IO::log << std::setw(10) << well_model(w).name();
  
  IO::log << "\n";

  for(int l = 0; l < well_size(); ++l) {
    //
    double nfac;
    
    for(int w = 0; w < well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * CONVERT_DD(well(w).weight_sqrt);

      dtemp = dtemp < 0. ? -dtemp : dtemp;

      if(!w || dtemp > nfac)
	//
	nfac = dtemp;
    }
    
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(10) << eigenval[l] / relax_eval_min
      //
	    << std::setw(10) << relaxation_projection[l];
    
    for(int w = 0; w < well_size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * CONVERT_DD(well(w).weight_sqrt) / nfac;
      
      IO::log << std::setw(10);
      
      if(dtemp > .01 || dtemp < -.01) {
	//
	IO::log << dtemp;
      }
      else
	//
	IO::log << 0;
    }
    
    IO::log << "\n";
  }

  IO::log << IO::log_offset
    //
	  << std::setw(5)  << "*Z"
    //
	  << std::setw(10) << "---"
    //
	  << std::setw(10) << "---";

  for(int w = 0; w < well_size(); ++w)
    //
    if(!w || dtemp < well(w).weight_sqrt)
      //
      dtemp = CONVERT_DD(well(w).weight_sqrt);
  
  
  for(int w = 0; w < well_size(); ++w)
    //
    IO::log << std::setw(10) << well(w).weight_sqrt / dtemp;
  
  IO::log << "\n";

  IO::log << IO::log_offset
    //
	  << "*R - eigenvalue over the relaxation limit\n"
    //
	  << IO::log_offset
    //
	  << "*P - eigenvector projection squared on the relaxation subspace (= 1 -F_ne)\n"
    //
	  << IO::log_offset
    //
	  << "*Z - well partition function square root (normalized)\n\n";
  
  /***************************************************************
   **************** CHEMICAL SUBSPACE DIMENSION ******************
   ***************************************************************/
  
  // number of low eigenvalues
  //
  if(chemical_tolerance > 0.) {
    //
  for(itemp = 0; itemp < well_size(); ++itemp)
    //
    if(eigenval[itemp] > chemical_tolerance * relax_eval_min)
      //
      break;
  //
    else
      //
      eigenval[itemp] = 0.;
  }
  else
    //
    itemp = 0;
  
  int low_size = 0;
  
  if((!bimolecular_size() && itemp > 1 || bimolecular_size() && itemp) && mode() != LUMPING) {
    //
    low_size = itemp;

    if(itemp == 1) {
      //
      IO::log << IO::log_offset << "WARNING: there is one low eigenvalue\n\n";
    }
    else
      //
      IO::log << IO::log_offset << "WARNING: there are " << itemp << " low eigenvalues\n\n";
  }

  // number of chemical eigenvalues
  //
  if(chemical_threshold > 1.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(relaxation_projection[itemp] * chemical_threshold > 1.)
	//
	break;
  }
  // absolute eigenvalue threshold
  //
  else if(chemical_threshold < 1. && chemical_threshold > 0.) {
    //
    for(itemp = 0; itemp < well_size(); ++itemp)
      //
      if(eigenval[itemp] > chemical_threshold * relax_eval_min)
	//
	break;
  }
  // relative eigenvalue threshold
  //
  else if(chemical_threshold < 0. && chemical_threshold > -1.) {
    //
    for(itemp = well_size(); itemp > 0; --itemp)
      //
      if(eigenval[itemp] <= 0. || eigenval[itemp - 1] / eigenval[itemp] < -chemical_threshold)
	//
	break;
  }
  //
  else {
    //
    std::cerr << funame << "chemical threshold has not been initialized properly: " << chemical_threshold << "\n";
    //
    throw Error::Logic();
  }
  
  const int chem_size = itemp;

  if(chem_size < low_size) {
    //
    std::cerr << funame << "the number of chemical eigenvalues cannot be less than the number of low eigenvalues\n";

    throw Error::Logic();
  }

  IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n\n";

  // well partitioning infrastructure initialization
  //
  Lapack::Matrix pop_chem;
    
  partition_t well_partition;
  
  group_t bimolecular_group;

  if(chem_size == well_size()) {
    //
    // no partitioning
    //
    pop_chem.resize(well_size());

    pop_chem = 0.;

    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

    well_partition.resize(well_size());
      
    for(int w = 0; w < well_size(); ++w)
      //
      well_partition[w].insert(w);

    IO::log << IO::log_offset << "partition projection error = "
      //
	    << (double)chem_size - projection(well_partition, pop_chem) << "\n\n";
  }

  /*******************************************************************************************
   ********************************* LOW EIGENVALUE REGIME ***********************************
   *******************************************************************************************/

  //low eigenvalue regime
  //
  if(low_size) {
    //
    if(mode() == BRANCHING) {
      //
      IO::log << IO::log_offset << "in low eigenvalue regime branching is impossible\n";
	
      throw BranchingException();
    }

    IO::Marker low_eigenvalue_marker("low eigenvalue regime");

    // low eigenvalue wells partition
    //
    if(low_size != well_size()) {
      //
      pop_chem.resize(well_size(), low_size);
    
      for(int l = 0; l < low_size; ++l)
	//
	pop_chem.column(l) = eigen_pop.row(l);

      // partitioning wells into equilibrated groups
      //
      dtemp = threshold_well_partition(pop_chem, well_partition, bimolecular_group);

      IO::log << IO::log_offset << "well partition: " << index_to_well(well_partition) << "\n";

      if(bimolecular_group.size())
	//
	IO::log << IO::log_offset << "bimolecular group: " << index_to_well(bimolecular_group) << "\n";

      IO::log << IO::log_offset << "partition projection error = " << dtemp << "\n\n";
    }

    std::vector<std::set<int> > owp = index_to_well(well_partition);

    std::set<int>               owg = index_to_well(bimolecular_group);

    // rate-limiting barrier
    //
    const int& type = control_barrier.first;

    const int& b    = control_barrier.second;

    switch(type) {
      //
    case Model::INNER:
      //
      if(graph.inner_set.find(b) == graph.inner_set.end()) {
	//
	std::cerr << funame << "(inner) rate-limiting barrier is not in the graph\n";

	throw Error::Logic();
      }
      
      break;
      //
    case Model::OUTER:
      //
      if(graph.outer_set.find(b) == graph.outer_set.end()) {
	//
	std::cerr << funame << "(outer) rate-limiting barrier is not in the graph\n";

	throw Error::Logic();
      }
      
      break;
      //
    default:
      //
      std::cerr << "rate-limiting barrier is not initialized\n";

      throw Error::Init();
    }
    
    // kinetic landscape
    //
    Model::emap_t ener_map;
    
    if(flags & THERMAL_ORDER) {
      //
      ener_map = graph.thermal_energy_map(temperature);
    }
    else
      //
      ener_map = graph.ground_energy_map();

    Model::landscape_t landscape = Model::kinetic_landscape(ener_map);
    
    Model::landscape_t::const_iterator li = landscape.find(control_barrier);

    if(li == landscape.end()) {
      //
      std::cerr << funame << "control barrier ";
	
      switch(control_barrier.first) {
	//
      case Model::INNER:
	//
	std::cerr << Model::inner_barrier(control_barrier.second).name();

	break;
	//
      case Model::OUTER:
	//
	std::cerr << Model::outer_barrier(control_barrier.second).name();
      }
	  
      std::cerr << " not in the landscape:";

      for(li = landscape.begin(); li != landscape.end(); ++li)
	//
	switch(li->first.first) {
	  //
	case Model::INNER:
	  //
	  std::cerr << "  " << Model::inner_barrier(li->first.second).name();

	  break;
	  //
	case Model::OUTER:
	  //
	  std::cerr << "  " << Model::outer_barrier(li->first.second).name();
	}

      std::cerr << "\n";
      
      throw Error::Logic();
    }
    
    // inner barrier
    //
    if(type == Model::INNER) {
      //
      if(li->second.size() != 2) {
	//
	std::cerr << funame << "wrong reactive groups number: " << li->second.size() << "\n";

	throw Error::Logic();
      }

      Model::ChemGraph rg1 = li->second.front();

      Model::ground(rg1.well_set, &itemp);
      
      const int w1 = itemp;
	
      Model::ChemGraph rg2 = li->second.back();

      Model::ground(rg2.well_set, &itemp);
      
      const int w2 = itemp;

      IO::log << IO::log_offset << "rate limiting barrier " << Model::inner_barrier(b).name()
	//
	      << ": " <<  rg1 << "(" << Model::well(w1).name() << ")" << "<-->" << Model::inner_barrier(b).name()
	//
	      << "<-->" << rg2 << "(" << Model::well(w2).name() << ")" << "\n";
    
      // bound groups associated with the rate-limiting barrier
      //
      int g1 = -1, g2 = -1;
	
      for(int g = 0; g < owp.size(); ++g) {
	//
	if(owp[g].find(w1) != owp[g].end())
	  //
	  g1 = g;

	if(owp[g].find(w2) != owp[g].end())
	  //
	  g2 = g;
      }

      // controlling barrier belongs to bimolecular group
      //
      if(g1 < 0 && g2 < 0) {
	//
	IO::log << IO::log_offset << "rate-limiting barrier is an inner barrier of the bimolecular group " << owg
	  //
		<< ": no low eigenvalue treatment is necessary\n";
      }
      // controlling barrier is internal barrier of the group
      //
      else if(g1 == g2) {
	//
	IO::log << IO::log_offset << "rate-limiting barrier is an inner barrier of the bound group " << owp[g1]
	  //
		<< ": no low eigenvalue treatment is necessary\n";
      }
      // control barrier is actually rate-limiting barrier, which connects two different groups
      //
      else {
	//
	if(g1 >= 0 && g2 >= 0) {
	//
	  IO::log << IO::log_offset << "rate-limiting barrier connects two bound groups: " << owp[g1]
	    //
		  << " and " << owp[g2] << "\n";

	  double ener = Model::outer_barrier(b).dist_ener_max(temperature);

	  dtemp = std::min(Model::state_density(owp[g1], ener), Model::state_density(owp[g2], ener));

	  IO::log << IO::log_offset << "the microscopic rate constant at the distribution energy maximum over relaxation time-scale = "
	    //
		  << Model::outer_barrier(b).states(ener) / dtemp / 2. / M_PI / relax_eval_min
	    //
		  << "\n";
	}
	else if(g1 >= 0) {
	  //
	  IO::log << IO::log_offset << "rate-limiting barrier connects the bound group " << owp[g1]
	  
		  << " and the bimolecular group " << owg << "\n";

	  double ener = Model::outer_barrier(b).dist_ener_max(temperature);

	  IO::log << IO::log_offset << "the microscopic rate constant at the distribution energy maximum over relaxation time-scale = "
	    //
		  << Model::outer_barrier(b).states(ener) / Model::state_density(owp[g1], ener) / 2. / M_PI / relax_eval_min
	    //
		  << "\n";
	}
	else {
	  //
	  IO::log << IO::log_offset << "rate-limiting barrier connects the bound group " << owp[g2]
	  
		  << " and the bimolecular group " << owg << "\n";

	  double ener = Model::outer_barrier(b).dist_ener_max(temperature);

	  IO::log << IO::log_offset << "the microscopic rate at the distribution energy maximum over relaxation time-scale = "
	    //
		  << Model::outer_barrier(b).states(ener) / Model::state_density(owp[g2], ener) / 2. / M_PI / relax_eval_min
	    //
		  << "\n";
	}
	
	// first subsystem
	//
	std::map<std::pair<int, int>, double> bf1;

	if(rg1.well_set.size() == 1 && !rg1.outer_set.size()) {
	  //
	  bf1[std::make_pair((int)Model::WELL, *rg1.well_set.begin())] = 1.;
	}
	else {
	  //
	  ReactiveComplex rc1(temperature, pressure, energy_reference, energy_cutoff, rg1);

	  bf1 = rc1.branching_fraction(control_barrier);
	}
	
	// second subsystem
	//
	std::map<std::pair<int, int>, double> bf2;
	
	if(rg2.well_set.size() == 1 && !rg2.outer_set.size()) {
	  //
	  bf2[std::make_pair((int)Model::WELL, *rg2.well_set.begin())] = 1.;
	}
	else {
	  //
	  ReactiveComplex rc2(temperature, pressure, energy_reference, energy_cutoff, rg2);

	  bf2 = rc2.branching_fraction(control_barrier);
	}
	
	// bimolecular and wells to internal indices mappings
	//
	std::map<int, int> b2i;

	std::map<int, int> w2i;

	int bi = 0, wi = 0;

	rate_data.well_partition.clear();
	
	rate_data.well_partition.resize(rg1.well_set.size() + rg2.well_set.size());
	  
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf1.begin(); i != bf1.end(); ++i)
	  //
	  switch(i->first.first) {
	    //
	  case Model::WELL:
	    //
	    rate_data.well_partition[wi].insert(i->first.second);

	    w2i[i->first.second] = wi++;

	    break;
	    //
	  case Model::BIMOLECULAR:
	    //
	    if(b2i.find(i->first.second) == b2i.end())
	      //
	      b2i[i->first.second] = bi++;

	    break;
	    //
	  default:
	    //
	    std::cerr << funame << "wrong reactant type: " << i->first.first << "\n";

	    throw Error::Logic();
	  }
	      
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf2.begin(); i != bf2.end(); ++i)
	  //
	  switch(i->first.first) {
	    //
	  case Model::WELL:
	    //
	    rate_data.well_partition[wi].insert(i->first.second);

	    w2i[i->first.second] = wi++;

	    break;
	    //
	  case Model::BIMOLECULAR:
	    //
	    if(b2i.find(i->first.second) == b2i.end())
	      //
	      b2i[i->first.second] = bi++;

	    break;
	    //
	  default:
	    //
	    std::cerr << funame << "wrong reactant type: " << i->first.first << "\n";

	    throw Error::Logic();
	  }

	rate_data.bimolecular_group.clear();
	
	for(std::set<int>::const_iterator w = graph.well_set.begin(); w != graph.well_set.end(); ++w)
	  //
	  if(rg1.well_set.find(*w) == rg1.well_set.end() && rg2.well_set.find(*w) == rg1.well_set.end())
	    //
	    rate_data.bimolecular_group.insert(*w);

	rate_data._bim_index_map = b2i;

	rate_data.ww_rate.resize(wi);

	rate_data.ww_rate = 0.;

	rate_data.bb_rate.resize(bi);

	rate_data.bb_rate = 0.;

	rate_data.wb_rate.resize(wi, bi);

	rate_data.wb_rate = 0.;

	dtemp = Model::inner_barrier(b).weight(temperature) / std::exp(Model::inner_barrier(b).ground() / temperature)

	  * temperature / 2. / M_PI;
	  
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf1.begin(); i != bf1.end(); ++i) {
	  //
	  for(std::map<std::pair<int, int>, double>::const_iterator j = bf2.begin(); j != bf2.end(); ++j) {
	    //
	    if(i->first.first == Model::BIMOLECULAR && j->first.first == Model::BIMOLECULAR) {
	      //
	      rate_data.bb_rate(b2i[i->first.second], b2i[j->first.second]) += dtemp * i->second * j->second;
	    }
	    else if(i->first.first == Model::WELL && j->first.first == Model::BIMOLECULAR) {
	      //
	      rate_data.wb_rate(w2i[i->first.second], b2i[j->first.second]) = dtemp * i->second * j->second;
	    }
	    else if(i->first.first == Model::WELL && j->first.first == Model::WELL) {
	      //
	      rate_data.ww_rate(w2i[i->first.second], w2i[j->first.second]) = dtemp * i->second * j->second;
	    }
	    else if(i->first.first == Model::BIMOLECULAR && j->first.first == Model::WELL) {
	      //
	      rate_data.wb_rate(w2i[j->first.second], b2i[i->first.second]) = dtemp * i->second * j->second;
	    }
	  }
	}

	rate_data.bw_rate = rate_data.wb_rate.transpose();

	return rate_data;
	//
      }// flux through the rate-limiting barrier
      //
    }// inner barrier
    //
    // outer barrier
    //
    else if(type == Model::OUTER) {
      //
      if(li->second.size() != 1) {
	//
	std::cerr << funame << "wrong reactive groups number: " << li->second.size() << "\n";

	throw Error::Logic();
      }

      Model::ChemGraph rg = li->second.front();

      Model::ground(rg.well_set, &itemp);

      const int w = itemp;
      
      IO::log << IO::log_offset << "rate-limiting barrier " <<  Model::outer_barrier(b).name() << ": "
	//
	      << rg << "(" << Model::well(w).name() << ")" << "<-->" << Model::outer_barrier(b).name() << "<-->"
	  //
		<< Model::bimolecular(Model::outer_connect(b).second).name() << "\n";
      
      int g;
	
      for(g = 0; g < owp.size(); ++g)
	//
	if(owp[g].find(w) != owp[g].end())
	  //
	  break;

      // barrier connects to bimolecular group
      //
      if(g == owp.size()) {
	//
	IO::log << IO::log_offset << "rate-limiting outer barrier conects to bimolecular group " << owg
	  //
		<< ": no low igenvalue treatment is necessary\n";
      }
      // barrier connects to the bound group
      //
      else {
	//
	IO::log << IO::log_offset << "rate-limiting barrier connects to the bound group: " << owp[g] << "\n";

	dtemp = Model::outer_barrier(b).dist_ener_max(temperature);

	IO::log << IO::log_offset << "the microscopic rate constant at the distribution energy maximum over relaxation time-scale = "
	  //
		<< Model::outer_barrier(b).states(dtemp) / Model::state_density(owp[g], dtemp) / 2. / M_PI / relax_eval_min
	  //
		<< "\n";
	
	std::map<std::pair<int, int>, double> bf;
	
	if(rg.well_set.size() == 1 && !rg.outer_set.size()) {
	  //
	  bf[std::make_pair((int)Model::WELL, *rg.well_set.begin())] = 1.;
	}
	else {
	  //
	  ReactiveComplex rc(temperature, pressure, energy_reference, energy_cutoff, rg);

	  bf = rc.branching_fraction(control_barrier);
	}

	std::map<int, int> b2i;

	std::map<int, int> w2i;

	int bi = 0, wi = 0;

	rate_data.well_partition.clear();
	
	rate_data.well_partition.resize(well_size());
	  
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i)
	  //
	  switch(i->first.first) {
	    //
	  case Model::WELL:
	    //
	    rate_data.well_partition[wi].insert(i->first.second);

	    w2i[i->first.second] = wi++;

	    break;
	    //
	  case Model::BIMOLECULAR:
	    //
	    b2i[i->first.second] = bi++;
	      
	    break;
	    //
	  default:
	    //
	    std::cerr << funame << "wrong reactant type: " << i->first.first << "\n";

	    throw Error::Logic();
	  }
	  
	if(b2i.find(Model::outer_connect(b).second) == b2i.end())
	  //
	  b2i[Model::outer_connect(b).second] = bi++;

	rate_data.bimolecular_group.clear();
	
	for(std::set<int>::const_iterator w = graph.well_set.begin(); w != graph.well_set.end(); ++w)
	  //
	  if(rg.well_set.find(*w) == rg.well_set.end())
	    //
	    rate_data.bimolecular_group.insert(*w);

	rate_data._bim_index_map = b2i;
	  
	rate_data.ww_rate.resize(wi);

	rate_data.ww_rate = 0.;

	rate_data.bb_rate.resize(bi);

	rate_data.bb_rate = 0.;

	rate_data.wb_rate.resize(wi, bi);

	rate_data.wb_rate = 0.;
	  
	dtemp = Model::outer_barrier(b).weight(temperature) / std::exp(Model::outer_barrier(b).ground() / temperature)

	  * temperature / 2. / M_PI;
	  
	for(std::map<std::pair<int, int>, double>::const_iterator i = bf.begin(); i != bf.end(); ++i) {
	  //
	  if(i->first.first == Model::BIMOLECULAR) {
	    //
	    rate_data.bb_rate(b2i[i->first.second], b2i[Model::outer_connect(b).second]) += dtemp * i->second;
	  }
	  else if(i->first.first == Model::WELL) {
	    //
	    rate_data.wb_rate(w2i[i->first.second], b2i[Model::outer_connect(b).second]) = dtemp * i->second;
	  }
	}

	rate_data.bw_rate = rate_data.wb_rate.transpose();

	return rate_data;
	//
      }// bound group
      //
    }// outer barrier
    //
  }// low eigenvalue regime
  
  /*********************************************************************
   ***************** EIGENSTATES-BIMOLECULAR COUPLING  *****************
   *********************************************************************/

  // eigenvectors projection on the bimolecular subspace;
  //
  LAPACK::Matrix eigen_bim;

  if(mode() != LUMPING) {
    //
    if(mode() == BRANCHING) {
      //
      eigen_bim.resize(global_size, 1);

      eigen_bim = 0.;

#pragma omp parallel for default(shared) private(dd_temp, itemp) schedule(static)
  
      for(int l = 0; l < global_size; ++l)
	//
	for(int w = 0; w < well_size(); ++w)
	  //
	  for(int e = 0; e < well(w).kernel_escape.size(); ++e) {
	    //
	    std::map<int, int>::const_iterator wit = kinetic_basis[e].well_index_map.find(w);

	    if(wit == kinetic_basis[e].well_index_map.end()) {
	      //
	      std::cerr << funame << "eigen_bim (escape): well is not in the kinetic basis space\n";

	      throw Error::Logic();
	    }

	    dd_temp = 0.;
	
	    for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	      //
	      dd_temp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(wit->second, k);
	    }

	    dd_temp *= well(w).kernel_escape[e] * well(w).boltzman_sqrt[e];

	    eigen_bim(l, 0) += dd_temp * pressure;
	  }
    }
    else if(bimolecular_size()) {
      //
      eigen_bim.resize(global_size, bimolecular_size());

      eigen_bim = 0.;

#pragma omp parallel for default(shared) private(dd_temp, itemp) schedule(static)
  
      for(int l = 0; l < global_size; ++l) {
	//
	for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	  //
	  const int w = bit->first.first;
      
	  const int p = bit->first.second;

	  for(int e = 0; e < bit->second.size(); ++e) {
	    //
	    std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	    if(i == kinetic_basis[e].well_index_map.end()) {
	      //
	      std::cerr << funame << "eigen_bim: well is not in the kinetic basis space\n";

	      throw Error::Logic();
	    }

	    dd_temp = 0.;
	
	    for(int k = 0; k < kinetic_basis[e].active_size; ++k) {
	      //
	      dd_temp += global_eigen(well_shift[e] + k, l) * kinetic_basis[e].eigenvector(i->second, k);
	    }

	    dd_temp *= bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e];

	    eigen_bim(l, p) += dd_temp / 2. / M_PI;
	  }
	}
      }
    }// bimolecular

    // high eigenstate correction
    //
    for(int l = 0; l < chem_size; ++l) {
      //
      // kernel escape
      //
      if(mode() == BRANCHING) {
	//
	for(int w1 = 0; w1 < well_size(); ++w1) {
	  //
	  float_t dd_val = 0.;

#pragma omp parallel for default(shared) private(dd_temp, itemp) reduction(+: dd_val) schedule(dynamic)

	  for(int e1 = 0; e1 < well(w1).size(); ++e1) {
	    //
	    std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);

	    if(i11 == kinetic_basis[e1].well_index_map.end()) {
	      //
	      std::cerr << funame << "high eigenstate correction: well w1 is not in the e1 map\n";

	      throw Error::Logic();
	    }
	
	    dd_temp = 0.;
	
	    for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
	      //
	      dd_temp += global_eigen(well_shift[e1] + k, l) * kinetic_basis[e1].eigenvector(i11->second, k);
	 
	    const float_t proj = dd_temp;

	    itemp = e1 + well(w1).kernel_bandwidth;

	    const int e2_max = itemp < well(w1).size() ? itemp : well(w1).size();
	  
	    itemp = e1 - well(w1).kernel_bandwidth + 1;

	    const int e2_min = itemp > 0 ? itemp : 0;
	    
	    for(int e2 = e2_min; e2 < e2_max; ++e2) {
	      //
	      std::map<int, int>::const_iterator i21 = kinetic_basis[e2].well_index_map.find(w1);

	      if(i21 == kinetic_basis[e2].well_index_map.end()) {
		//
		std::cerr << funame << "high eigenstate correction: well w1 is not in the e2 map\n";

		throw Error::Logic();
	      }

	      for(int w2 = 0; w2 < well_size(); ++w2) {
		//
		if(e2 >= well(w2).kernel_escape.size())
		  //
		  continue;
	      
		std::map<int, int>::const_iterator i22 = kinetic_basis[e2].well_index_map.find(w2);
	      
		if(i22 == kinetic_basis[e2].well_index_map.end()) {
		  //
		  std::cerr << funame << "high eigenstate eigen_bim correction (kernel escape): well w2 is not in the e2 map\n";

		  throw Error::Logic();
		}
	
		dd_temp = 0.;
	      
		for(int k = kinetic_basis[e2].active_size; k < kinetic_basis[e2].size(); ++k)
		  //
		  dd_temp += kinetic_basis[e2].eigenvector(i21->second, k)
		    //
		    * kinetic_basis[e2].eigenvector(i22->second, k)
		    //
		    / kinetic_basis[e2].eigenvalue[k];

		dd_temp *= proj * well(w1).kernel(e1, e2)
		  //
		  * well(w1).boltzman_sqrt[e1] / well(w1).boltzman_sqrt[e2]
		  //
		  * well(w2).kernel_escape[e2] * well(w2).boltzman_sqrt[e2];
	    
		dd_val -= dd_temp;
		//
	      }// kernel escape
	      //
	    }// e2 cycle
	    //
	  }// e1 cycle

	  dd_val *= pressure * pressure;
	
	  eigen_bim(l, 0) += dd_val;
	  //
	}// w1 cycle
	//
      }// kernel escape
      //
      // bimolecular escape
      //
      else {
	//
	for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	  //
	  const int w = bit->first.first;
	
	  const int p  = bit->first.second;

	  float_t dd_val = 0.;
	
#pragma omp parallel for default(shared) private(dd_temp, itemp) reduction(+: dd_val) schedule(dynamic)
	
	  for(int e = 0; e < bit->second.size(); ++e) {
	    //
	    std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);
		  
	    if(i == kinetic_basis[e].well_index_map.end()) {
	      //
	      std::cerr << funame << "high eigenstate eigen_bim correction: well w is not in the e map\n";

	      throw Error::Logic();
	    }

	    for(int w1 = 0; w1 < well_size(); ++w1) {
	      //
	      std::map<int, int>::const_iterator i1 = kinetic_basis[e].well_index_map.find(w1);
		  
	      if(i1 == kinetic_basis[e].well_index_map.end())
		//
		continue;
	    
	      float_t proj = 0.;
	      
	      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k)
		//
		proj += kinetic_basis[e].eigenvector(i1->second, k)
		  //
		  * kinetic_basis[e].eigenvector(i->second, k)
		  //
		  / kinetic_basis[e].eigenvalue[k];

	      proj *= bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e];

	      itemp = e + well(w1).kernel_bandwidth;

	      const int e1_max = itemp < well(w1).size() ? itemp : well(w1).size();
	  
	      itemp = e - well(w1).kernel_bandwidth + 1;

	      const int e1_min = itemp > 0 ? itemp : 0;
	    
	      for(int e1 = e1_min; e1 < e1_max; ++e1) {
	      
		std::map<int, int>::const_iterator i11 = kinetic_basis[e1].well_index_map.find(w1);
		  
		if(i11 == kinetic_basis[e1].well_index_map.end()) {
		  //
		  std::cerr << funame << "high eigenstate eigen_bim correction: well w1 is not in the e1 map\n";

		  throw Error::Logic();
		}

		dd_temp = 0.;
	
		for(int k = 0; k < kinetic_basis[e1].active_size; ++k)
		  //
		  dd_temp += global_eigen(well_shift[e1] + k, l) * kinetic_basis[e1].eigenvector(i11->second, k);
	 
		dd_temp *= proj * well(w1).kernel(e1, e) * well(w1).boltzman_sqrt[e1] / well(w1).boltzman_sqrt[e];
	    
		dd_val -= dd_temp;
		//
	      }// e1 cycle 
	      //
	    } // w1 cycle
	    //
	  }// e cycle

	  dd_val *=  pressure / 2. / M_PI;
	
	  eigen_bim(l, p) += dd_val;
	  //
	}// outer barrier cycle
	//
      }// bimolecular escape
      //
    }// eigenstate cycle
    //
  }// no low eigenvalue regime or lumping
  
  /***********************************************************************************************
   ************************************* CHEMICAL SUBSPACE ***************************************
   ***********************************************************************************************/

  if(!chem_size) {
    //
    // no kinetically distinct species
    //
    group_t  owg;

    for(int i = 0; i < well_size(); ++i)
      //
      owg.insert(index_to_well(i));

    rate_data.bimolecular_group = owg;

    if(mode() == LUMPING)
      //
      return rate_data;

    if(mode() == BRANCHING) {
      //
      IO::log << IO::log_offset << "WARNING: no bound groups\n";
      
      throw BranchingException();
    }
  }
  if(chem_size) {
    //
    // projection of the chemical eigenvectors onto the thermal subspace
    //
    if(chem_size != well_size())  {
      //
      pop_chem.resize(well_size(), chem_size);
    
      for(int l = 0; l < chem_size; ++l)
	//
	pop_chem.column(l) = eigen_pop.row(l);

#ifdef DEBUG

      IO::log << IO::log_offset << "orthogonality check:\n";
      
      IO::log_offset.increase();

      double proj = 0.;
      
      for(int l = 0; l < chem_size; ++l) {
	//
	// normalize
	//
	for(int m = 0; m < l; ++m) {
	  //
	  dtemp = vdot(pop_chem.column(l), pop_chem.column(m));
	  
	  dtemp = dtemp >= 0. ? dtemp : -dtemp;
	  
	  proj = dtemp > proj ? dtemp : proj;
	}
      }

      IO::log << IO::log_offset << "maximal scalar product of different chemical eigenvectors = "
	//
	      << proj << "\n";

      for(int l = 0; l < chem_size; ++l) {
	//
	dtemp = vdot(pop_chem.column(l));
	
	if(!l || dtemp < proj)
	  //
	  proj = dtemp;
      }

      IO::log << IO::log_offset << "minimal chemical eigenvector square = "
	//
	      << proj << "\n";

      for(int l = 0; l < chem_size; ++l) {
	//
	dtemp = vdot(pop_chem.column(l));
	
	if(!l || dtemp > proj)
	  //
	  proj = dtemp;
      }

      IO::log << IO::log_offset << "maximal chemical eigenvector square = "
	//
	      << proj << "\n";

      IO::log_offset.decrease();
      
      IO::log << IO::log_offset << "orthogonality check done\n\n";

#endif

      // partitioning wells into equilibrated groups
      //
      dtemp = threshold_well_partition(pop_chem, well_partition, bimolecular_group);

      // convert chemical eigenvectors in the new basis
      //
      pop_chem = basis(well_partition).transpose() * pop_chem;
    }

    // well partition in original well representation
    //
    partition_t owp(well_partition.size());

    for(int g = 0; g < well_partition.size(); ++g)
      //
      for(group_t::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w)
	//
	owp[g].insert(index_to_well(*w));

    group_t  owg;

    for(group_t::const_iterator w = bimolecular_group.begin(); w != bimolecular_group.end(); ++w)
      //
      owg.insert(index_to_well(*w));

    if(chem_size != well_size()) {
      //
      IO::log << IO::log_offset << "well partition:";

      for(int g = 0; g < owp.size(); ++g) {
	//
	IO::log << "   ";
      
	for(group_t::const_iterator w = owp[g].begin(); w != owp[g].end(); ++w) {
	  //
	  if(w != owp[g].begin())
	    //
	    IO::log << "+";
	
	  IO::log << Model::well(*w).name();
	}

	IO::log << "/" << g;
      }

      IO::log << "\n";
    }

    if(bimolecular_group.size()) {
      //
      IO::log << IO::log_offset << "bimolecular group:";

      for(group_t::const_iterator w = owg.begin(); w != owg.end(); ++w)
	//
	IO::log << "   " << Model::well(*w).name();

      IO::log << "\n";
    }
    
    rate_data.well_partition = owp;

    rate_data.bimolecular_group = owg;

    // for well lumping only well partitioning data are needed
    //
    if(mode() == LUMPING)
      //
      return rate_data;
    
    Lapack::Matrix  m_direct = pop_chem;
  
    Lapack::Matrix m_inverse = m_direct.invert();

    std::vector<double>  weight = group_weight(well_partition);

    // branching fraction data
    //
    if(mode() == BRANCHING) {
      //
      rate_data.branching_fraction.resize(chem_size);

      double pval = 0.;

      double nval = 0.;
      
      for(int g = 0; g < chem_size; ++g) {
	//
	dtemp = 0.;
	  
	for(int l = 0; l < chem_size; ++l)
	  //
	  dtemp += m_inverse(l, g) * CONVERT_DD(eigen_bim(l, 0));
	  
	dtemp *= std::sqrt(weight[g]); // * vdot(m_inverse.column(g), &eigen_bim(0, 0))

	if(dtemp < 0.) {
	  //
	  nval -= dtemp;
	  
	  dtemp = 0.;
	}
	
	rate_data.branching_fraction[g] = dtemp;

	pval += dtemp;
      }

      if(pval == 0.) {
	//
	IO::log << IO::log_offset << "WARNING: all positive branching fractions are zero\n";

	throw BranchingException();
      }

      dtemp = nval / (pval + nval);

      if(dtemp > 0.01) {
	//
	IO::log << IO::log_offset << "WARNING: negative branching fraction = " << dtemp << "\n";

	throw BranchingException();
      }
      
      for(int g = 0; g < chem_size; ++g)
	//
	rate_data.branching_fraction[g] /= pval;

      if(chem_size > 1) {
	//
	IO::log << IO::log_offset << "branching fraction/group index:";
    
	for(int g = 0; g < chem_size; ++g)
	  //
	  IO::log << "  " << rate_data.branching_fraction[g] << "/" << g;
    
	IO::log << "\n";
      }

      return rate_data;
    }
  
#ifdef DEBUG

    /*
      IO::log << IO::log_offset << "partition groups weights:";

      for(int g = 0; g < chem_size; ++g)
      //
      IO::log << "   " << weight[g] << "/" << g;

      IO::log << "\n";
    
      Lapack::Matrix one(m_direct * m_inverse);
    
      one.diagonal() -= 1.;
    
      double val_max = -1.;
    
      for(int i = 0; i < one.size1(); ++i)
      //
      for(int j = 0; j < one.size2(); ++j) {
      //
      dtemp = one(i, j);
	
      dtemp = dtemp < 0. ? -dtemp : dtemp;
	
      if(dtemp > epsilon && dtemp > val_max)
      //
      val_max = dtemp;
      }
    
      if(val_max > 0.)
      //
      IO::log << IO::log_offset << funame
      //
      << "WARNING: matrix inversion error = " << val_max
      //
      << " exceeds numerical accuracy = " << epsilon
      //
      << "\n";
    */
    
#endif
    
    /***************************************************************************
     ********************** WELL-TO-WELL RATE COEFFICIENTS *********************
     ***************************************************************************/

    Lapack::Matrix ww_rate(chem_size);

    if(mode() != BRANCHING) {
      //
      ww_rate = 0.;
    
      for(int i = 0; i < chem_size; ++i)
	//
	for(int j = 0; j < chem_size; ++j) {
	  //
	  for(int l = 0; l < chem_size; ++l)
	    //
	    ww_rate(i, j) += m_direct(j, l) * m_inverse(l, i) * CONVERT_DD(eigenval[l]);

	  ww_rate(i, j) *= -std::sqrt(weight[i] * weight[j])
	    //
	    * energy_step / std::exp(energy_reference / temperature);
	}

      rate_data.ww_rate = ww_rate;

#ifdef DEBUG
      
      IO::log << IO::log_offset << "isomerization rate constants ratios (G - well group index):\n"
	//
	      << IO::log_offset << std::setw(3) << "G\\G";
    
      for(int i = 0; i < chem_size; ++i)
	//
	IO::log << std::setw(10) << i;
    
      IO::log << "\n";
    
      for(int j = 0; j < chem_size; ++j) {
	//
	IO::log << IO::log_offset << std::setw(3) << j;
      
	for(int i = 0; i < chem_size; ++i)
	  //
	  if(i != j) {
	    //
	    if(ww_rate(j, i) != 0.) {
	      //
	      IO::log << std::setw(10) << ww_rate(i, j) / ww_rate(j, i);
	    }
	    else
	      //
	      IO::log << std::setw(10) << "***";
	  }
	  else
	    //
	    IO::log << std::setw(10) << "1";
      
	IO::log << "\n";
      }

#endif

    }
  
    
    /******************************************************************************************************************
     ********************** WELL-TO-BIMOLECULAR RATE COEFFICIENTS / WELL-TO-WELL BRANCHING RATIOS *********************
     ******************************************************************************************************************/

    if(bimolecular_size()) {
      //
      Lapack::Matrix wb_rate(chem_size, bimolecular_size()), bw_rate(bimolecular_size(), chem_size);
    
      for(int g = 0; g < chem_size; ++g)
	//
	for(int p = 0; p < bimolecular_size(); ++p) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_inverse(l, g) * CONVERT_DD(eigen_bim(l, p));
	  
	  wb_rate(g, p) = dtemp //vdot(m_inverse.column(g), &eigen_bim(0, p))
	    //
	    * std::sqrt(weight[g]) * energy_step / std::exp(energy_reference / temperature);
	}
      
      for(int p = 0; p < bimolecular_size(); ++p)
	//
	for(int g = 0; g < chem_size; ++g) {
	  //
	  dtemp = 0.;
	  
	  for(int l = 0; l < chem_size; ++l)
	    //
	    dtemp += m_direct(g, l) * CONVERT_DD(eigen_bim(l, p));
	  
	  bw_rate(p, g) = dtemp //vdot(m_direct.row(g), &eigen_bim(0, p))
	    //
	    * std::sqrt(weight[g]) * energy_step / std::exp(energy_reference / temperature);
	}
      
      rate_data.bw_rate = bw_rate;

      rate_data.wb_rate = wb_rate;

#ifdef DEBUG

      IO::log << IO::log_offset << "w->b/b->w ratios:\n";

      for(int p = 0; p < bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(6) << bimolecular(p).name();
	
	for(int g = 0; g < chem_size; ++g)
	  //
	  if(bw_rate(p, g) != 0.) {
	    //
	    IO::log << std::setw(10) << wb_rate(g, p) / bw_rate(p, g);
	  }
	  else
	    //
	    IO::log << std::setw(10) << "***";

	IO::log << "\n";
      }//

      /*
	IO::log << IO::log_offset << "well-to-bimolecular rate constants:\n";

	IO::log << IO::log_offset << std::setw(3) << "W\\P";

	for(int p = 0; p < bimolecular_size(); ++p)
	//
	IO::log << std::setw(10) << bimolecular(p).name();

	IO::log << "\n";

	for(int g = 0; g < chem_size; ++g) {
	//
	IO::log << IO::log_offset << std::setw(3) << g;

	for(int p = 0; p < bimolecular_size(); ++p)
	//
	IO::log << std::setw(10) << wb_rate(g, p) / Model::weight(owp[g], temperature) / Phys_const::herz;

	IO::log << "\n";
	}
      */
      
#endif
	
    }// well-bimolecular rates
    //
  }// bound species

  /*****************************************************************************************
   ************************** BIMOLECULAR-TO-BIMOLECULAR RATES *****************************
   *****************************************************************************************/

  if(bimolecular_size()) {
    //
    const int relax_size = global_size - chem_size;
        
    Lapack::SymmetricMatrix bb_rate;
  
    bb_rate.resize(bimolecular_size());

    for(int p = 0; p < bimolecular_size(); ++p)
      //
      for(int q = 0; q < p; ++q) {
	//
	dd_temp = 0.;
	
	for(int l = chem_size; l < global_size; ++l)
	  //
	  dd_temp += eigen_bim(l, p) * eigen_bim(l, q) / eigenval[l];
	  
	bb_rate(p, q) = CONVERT_DD(dd_temp);
      }

    // high eigenvalue contribution
    //
    for(int e = 0; e < ener_index_max; ++e)
      //
      for(int k = kinetic_basis[e].active_size; k < kinetic_basis[e].size(); ++k) {
	//
	LAPACK::Vector vtemp(bimolecular_size(), 0.);
	
	for(outer_t::const_iterator bit = outer_barrier.begin(); bit != outer_barrier.end(); ++bit) {
	  //
	  const int w = bit->first.first;
	  
	  const int p = bit->first.second;

	  if(e >= bit->second.size())
	    //
	    continue;
	  
	  std::map<int, int>::const_iterator i = kinetic_basis[e].well_index_map.find(w);

	  if(i == kinetic_basis[e].well_index_map.end()) {
	    //
	    std::cerr << funame << "size of the well associated with the barrier is smaller than the barrier size\n";

	    throw Error::Logic();
	  }
	  
	  vtemp[p] += kinetic_basis[e].eigenvector(i->second, k)
	    //
	    * bit->second.state_number[e] / well(w).state_density[e] * well(w).boltzman_sqrt[e] / 2. / M_PI;
	}
	
	for(int p = 0; p < bimolecular_size(); ++p)
	  //
	  for(int q = 0; q < p; ++q) {
	    //
	    dd_temp = vtemp[p] * vtemp[q] / kinetic_basis[e].eigenvalue[k];
	    
	    bb_rate(p, q) += CONVERT_DD(dd_temp);
	  }
      }
    
    bb_rate *= energy_step / std::exp(energy_reference / temperature);

    rate_data.bb_rate = bb_rate;
    //
  }// bimolecular-to-bimolecular rates

  return rate_data;
}
