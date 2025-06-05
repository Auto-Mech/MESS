#indclude "mess_2d.hh"

void Mess::direct_diagonalization_method (std::map<std::pair<int, int>, double>& rate_data,
						    //
						    Partition& well_partition, int flags)
  
{
  const char funame [] = "Mess2::direct_diagonalization_method: ";

  // bimolecular rate units
  //
  const double bru = Phys_const::cm * Phys_const::cm * Phys_const::cm * Phys_const::herz;

  if(!isset()) {
    //
    std::cerr << funame << "reactive complex is not set\n";
    
    throw Error::Init();
  }

  IO::Marker funame_marker(funame);

  // pressure output
  //
  IO::log << IO::log_offset << "Pressure = ";
  
  switch(pressure_unit) {
  case BAR:
    IO::log << pressure() / Phys_const::bar << " bar";
    break;
  case TORR:
    IO::log << pressure() / Phys_const::tor << " torr";
    break;
  case ATM:
    IO::log << pressure() / Phys_const::atm << " atm";
    break;
  }

  // temperature output
  //
  IO::log << "\t Temperature = "
    //
	  << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";

  IO::aux << "Pressure = " << pressure() / Phys_const::bar << " bar"
    //
	  << "\t Temperature = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << " K\n";
  
  // collisional frequency output
  //
  IO::log << IO::log_offset
	  << std::setw(Model::log_precision + 7) << "Well"
	  << std::setw(20) << "Collision, 1/sec"
	  << "\n";
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    IO::log << IO::log_offset
	    << std::setw(Model::log_precision + 7) << Model::well(w).short_name()
	    << std::setw(20) << well(w).collision_frequency() / Phys_const::herz
	    << "\n";

  rate_data.clear();

  int                 itemp;
  double              dtemp;
  bool                btemp;
  std::string         stemp;

  Lapack::Vector      vtemp;
  Lapack::Matrix      mtemp;

  double proj;

  // total kinetic relaxation matrix dimension and index shifts for individual wells
  //
  std::vector<int> well_shift(well.size());
  
  itemp = 0;
  
  for(int w = 0; w < well.size(); itemp += well(w++).size())
    //
    well_shift[w] = itemp;

  const int global_size = itemp;

  IO::log << IO::log_offset << "global relaxation matrix dimension = " << global_size << "\n";

  /********************************* SETTING GLOBAL MATRICES *********************************/

  // kinetic relaxation matrix
  //
  Lapack::SymmetricMatrix kin_mat(global_size); // kinetic relaxation matrix
  
  kin_mat = 0.;

  // bimolecular product vectors
  //
  Lapack::Matrix global_bim;
  
  if(Model::bimolecular_size()) {
    //
    global_bim.resize(global_size,  Model::bimolecular_size());
    
    global_bim = 0.;
  }

  // Boltzmann distributions
  //
  Lapack::Matrix global_pop(global_size, well.size());
  
  global_pop = 0.;

  {
    IO::Marker set_marker("setting global matrices", IO::Marker::ONE_LINE);

    // kin_mat initialization
    
    //  inner barrier isomerization contribution
    //
    for(int b = 0; b < inner_barrier.size(); ++b) {
      //
      const int w1 = Model::inner_connect(b).first;
      
      const int w2 = Model::inner_connect(b).second;
      
      for(int i = 0; i < inner_barrier[b].size(); ++i) {
	//
	const int e = inner_barrier[b].ener_index(i);

	const int a = inner_barrier[b].amom_index(i);

        const int i1 = well[w1].linear_index(e, a);

        const int i2 = well[w2].linear_index(e, a);

	if(i1 < 0 || i2 < 0)
	  //
	  continue;

	dtemp = inner_barrier[b].states(i) / 2. / M_PI;
	
	kin_mat(i1 + well_shift[w1], i2 + well_shift[w2]) += -dtemp / std::sqrt(well[w1].states(i1) * well[w2].states(i2));

	kin_mat(i1 + well_shift[w1], i1 + well_shift[w1]) += dtemp / well[w1].states(i1);

	kin_mat(i2 + well_shift[w2], i2 + well_shift[w2]) += dtemp / well[w2].states(i2);
    }

    //  outer barrier isomerization contribution and bimolecular product vectors
    //
    for(int b = 0; b < outer_barrier.size(); ++b) {
      //
      const int w = Model::outer_connect(b).first;  // well index
      
      const int p = Model::outer_connect(b).second; // product index

      for(int i = 0; i < outer_barrier[b].size(); ++i) {
	//
	const int e = outer_barrier[b].ener_index(i);

	const int a = outer_barrier[b].amom_index(i);

        const int wi = well[w].linear_index(e, a);

	if(wi < 0)
	  //
	  continue;

	dtemp = outer_barrier[b].states(i) / 2. / M_PI;

	// outer barrier kinetic matrix contribution
	//
	kin_mat(wi + well_shift[w], wi + well_shift[w]) += dtemp / well[w].states(wi);

	// bimolecular product vectors
	//
	global_bim(wi + well_shift[w], p) = dtemp * std::sqrt(thermal_factor(e) / well[w].states(wi));
      }

    // collision relaxation contribution
    //
    for(int w = 0; w < well.size(); ++w) {
      //
      for(int i = 0; i < well[w].size(); ++i) {
	//
	for(int j = i; j < well[w].size(); ++j)
	  //
	  // symmetrized collisional relaxation kernel
	  //
	  kin_mat(i + well_shift[w], j + well_shift[w]) +=  well[w].collision_frequency() * well[w].kernel(i, j);
      }
    }

    // thermal distributions
    //
    for(int w = 0; w < well.size(); ++w)
      //
      for(int i = 0; i < well[w].size(); ++i) {
	//
	const int e = well[w].ener_index(i);
	
	global_pop(i + well_shift[w], w) = std::sqrt(well[w].states(i) * boltzmann_factor(e) / well[w].weight());
      }
    //
  }// global matrices

  /******************** DIAGONALIZING THE GLOBAL KINETIC RELAXATION MATRIX ********************/

  Lapack::Vector eigenval;
  
  Lapack::Matrix eigen_global(global_size);

  { IO::Marker solve_marker("diagonalizing global relaxation matrix", IO::Marker::ONE_LINE);

    if(Mpack::mp_type == Mpack::DOUBLE) {
      //
      eigenval = Offload::eigenvalues(kin_mat, &eigen_global);
    }
    else if(use_mp) {
      //
      eigenval = Mpack::eigenvalues(kin_mat, &eigen_global);
    }
    else {
      //
      eigenval = kin_mat.eigenvalues(&eigen_global);
    }

    eigen_global = eigen_global.transpose();
  }

  const double min_relax_eval = eigenval[well.size()];
  
  const double max_relax_eval = eigenval.back();

  Lapack::Matrix eigen_well(global_size, well.size());
  
  for(int l = 0; l < global_size; ++l)
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      eigen_well(l, w) = vlength(&eigen_global(l, well_shift[w]), well[w].size(), global_size);
  
  // projection of the  eigenvectors onto the thermal subspace
  //
  Lapack::Matrix eigen_pop = eigen_global * global_pop;

  // eigenvector to bimolecular vector projection
  //
  Lapack::Matrix eigen_bim;
  
  if(Model::bimolecular_size())
    //
    eigen_bim = eigen_global * global_bim;

  std::vector<double> relaxation_projection(well.size());
  
  for(int l = 0; l < well.size(); ++l)
    //
    relaxation_projection[l] = 1. - vdot(eigen_pop.row(l));

 
  /************************************* EIGENVECTOR OUTPUT *****************************************/

  IO::log  << IO::log_offset << "eigenvector populations normalized:\n"
    //
	   << IO::log_offset
    //
	   << std::setw(5)  << "L"
    //
	   << std::setw(Model::log_precision + 7) << "*R"
    //
	   << std::setw(Model::log_precision + 7) << "*P";
  
  for(int w = 0; w < well.size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
  IO::log << "\n";

  // maximal population
  //
  for(int l = 0; l < well.size(); ++l) {
    //
    double pos_pop = 0.;
    
    double neg_pop = 0.;
    
    for(int w = 0; w < well.size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * std::sqrt(well[w].weight());
      
      if(dtemp > 0.) {
	//
	pos_pop += dtemp;
      }
      if(dtemp < 0.)
	//
	neg_pop += dtemp;
    }
    
    double max_pop = pos_pop > -neg_pop ? pos_pop : neg_pop;
    
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(Model::log_precision + 7) << eigenval[l] / min_relax_eval
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < well.size(); ++w) {
      //
      dtemp = eigen_pop(l, w) * std::sqrt(well[w].weight()) / max_pop;
      
      IO::log << std::setw(Model::log_precision + 7);
      
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
	  << "*R - eigenvalue over the relaxation limit\n"
    //
	  << IO::log_offset
    //
	  << "*P - eigenvector projection squared on the relaxation subspace (= 1 -F_ne)\n";
  
  IO::log << IO::log_offset << "eigenvector projections:\n"
    //
	  << IO::log_offset << std::setw(5)  << "L" << std::setw(Model::log_precision + 7) << "*Q" << std::setw(Model::log_precision + 7) << "*P";

  for(int w = 0; w < well.size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::well(w).short_name();
  
  IO::log << "\n";

  for(int l = 0; l < well.size(); ++l) {
    //
    IO::log << IO::log_offset
      //
	    << std::setw(5)  << l
      //
	    << std::setw(Model::log_precision + 7) << eigenval[l] / well.front().collision_frequency()
      //
	    << std::setw(Model::log_precision + 7) << relaxation_projection[l];
    
    for(int w = 0; w < well.size(); ++w)
      //
      IO::log << std::setw(Model::log_precision + 7) << eigen_pop(l,w);
    
    IO::log << "\n";
  }
  
  IO::log << IO::log_offset
    //
	  << std::setw(5) << "*Z"
    //
	  << std::setw(Model::log_precision + 7) << "---"
    //
	  << std::setw(Model::log_precision + 7) << "---";
  
  for(int w = 0; w < well.size(); ++w)
    //
    IO::log << std::setw(Model::log_precision + 7) << std::sqrt(well[w].weight());
  
  IO::log << "\n";

  IO::log << IO::log_offset
    //
	  << "*Q - eigenvalue over the collision frequency in first well\n"
    //
	  << IO::log_offset
    //
	  << "*P - eigenvector projection squared on the relaxation subspace (= 1 - F_ne)\n"
    //
	  << IO::log_offset
    //
	  << "*Z - well partition function square root\n";
  
  // eigenvalues output
  //
  if(eval_out.is_open()) {
    //
    eval_out << std::setw(13) << temperature() / Phys_const::kelv
      //
	     << std::setw(13);
    
    switch(pressure_unit) {
      //
    case BAR:
      //
      eval_out << pressure() / Phys_const::bar;
      break;
      
    case TORR:
      //
      eval_out << pressure() / Phys_const::tor;
      break;
      
    case ATM:
      //
      eval_out << pressure() / Phys_const::atm;
      break;
    }
    
    eval_out << std::setw(13) << well(0).collision_frequency() / Phys_const::herz
      //
	     << std::setw(13) << min_relax_eval / well(0).collision_frequency();
    
    int eval_max = Model::well_size() + evec_out_num;
    
    for(int l = 0; l < eval_max; ++l)
      //
      eval_out << std::setw(13) << eigenval[l] / well(0).collision_frequency()
	//
	       << std::setw(13) <<  1. - vdot(eigen_pop.row(l));
    
    eval_out << "\n";
  }

  // eigenvector output
  //
  if(evec_out.is_open()) {
    //
    evec_out << "EIGENVECTORS:\n";
    
    int evec_max = Model::well_size() + evec_out_num;
    
    for(int l = 0; l < evec_max; ++l) {
      //
      evec_out << "l = " << l << "\n"
	//
	       << "eigenvalue / collision frequency = "
	//
	       << eigenval[l] / well(0).collision_frequency()
	//
	       << "\n";

      evec_out << std::setw(13) << "well length";
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << eigen_well(l, w);
      
      evec_out << "\n";

      evec_out << std::setw(13) << "E, kcal/mol";
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << Model::well(w).short_name();
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	evec_out << std::setw(13) << Model::well(w).short_name();
      
      evec_out << "\n";

      for(int i = 0; i < well_size_max; ++i) {
	//
	evec_out << std::setw(13) << energy_bin(i) / Phys_const::kcal;
	
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  if(i < well(w).size()) {
	    //
	    // micropopulational distribution
	    //
	    dtemp = eigen_global(l, well_shift[w] + i) * well(w).boltzman_sqrt(i);
	    
	    evec_out << std::setw(13) << dtemp;
	  }
	  else
	    //
	    evec_out << std::setw(13) << 0;
	
	for(int w = 0; w < Model::well_size(); ++w)
	  //
	  if(i < well(w).size()) {
	    //
	    // ratio to the thermal distribution
	    //
	    dtemp = eigen_global(l, well_shift[w] + i) / well(w).boltzman_sqrt(i)
	      //
	      * well(w).weight_sqrt();
	    
	    evec_out << std::setw(13) << dtemp;
	  }
	  else
	    //
	    evec_out << std::setw(13) << 0;
	
	evec_out << "\n";
      }
    }
    
    evec_out << "\n";
  }

  /********************************** CHEMICAL SUBSPACE DIMENSION ****************************************/

  //
  // predefined chemical subspace dimension
  //
  if(_default_chem_size >= 0) {
    //
    itemp = _default_chem_size;
  }
  //
  // default partitioning scheme
  //
  else if(default_partition.size()) {
    //
    itemp = default_partition.size();
  }
  //
  // absolute eigenvalue threshold
  //
  else if(chemical_threshold > 1.) {
    //
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(min_relax_eval / eigenval[itemp]  < chemical_threshold)
	//
	break;
  }
  //
  // relaxation projection threshold
  //
  else if(chemical_threshold < 1. && chemical_threshold > 0.) {
    //
    for(itemp = 0; itemp < Model::well_size(); ++itemp)
      //
      if(relaxation_projection[itemp] > chemical_threshold)
	//
	break;
  }
  //
  // relative eigenvalue threshold
  //
  else if(chemical_threshold < -1. ) {
    //
    for(itemp = Model::well_size(); itemp > 0; --itemp)
      //
      if(eigenval[itemp] / eigenval[itemp - 1] > - chemical_threshold)
	//
	break;
  }
  //
  else
    //
    itemp = Model::well_size();
	
  const int chem_size = itemp;

  if(chem_size != Model::well_size())
    //
    IO::log << IO::log_offset << "dimension of the chemical subspace = " << chem_size << "\n";

  /***** PARTITIONING THE GLOBAL PHASE SPACE INTO THE CHEMICAL AND COLLISIONAL SUBSPACES *****/

  // collisional relaxation eigenvalues and eigenvectors
  //
  const int relax_size = global_size - chem_size;
  
  Lapack::Vector relax_lave(relax_size);
  
  for(int r = 0; r < relax_size; ++r) {
    //
    itemp = r + chem_size;

    relax_lave[r] = 1. / eigenval[itemp];
  }

  // kinetic matrix modified

#pragma omp parallel for default(shared) private(dtemp) schedule(dynamic)
	
  for(int i = 0; i < global_size; ++i) {
    //
    for(int j = i; j < global_size; ++j) {
      //
      dtemp = 0.;
      
      for(int l = 0; l < chem_size; ++l)
	//
	dtemp += eigen_global(l, i) * eigen_global(l, j);
      
      kin_mat(i, j) += dtemp * well(0).collision_frequency();
    }
  }

  Lapack::Matrix proj_bim = global_bim.copy();
  
  for(int p = 0; p < Model::bimolecular_size(); ++p)
    //
    for(int l = 0; l < chem_size; ++l)
      //
      parallel_orthogonalize(&proj_bim(0, p), &eigen_global(l, 0), global_size, 1, global_size);

  Lapack::Matrix inv_proj_bim;
  
  if(Model::bimolecular_size())
    //
    inv_proj_bim = Lapack::Cholesky(kin_mat).invert(proj_bim);

  Lapack::Matrix proj_pop = global_pop.copy();
  
  for(int w = 0; w < Model::well_size(); ++w)
    //
    for(int l = 0; l < chem_size; ++l)
      //
      parallel_orthogonalize(&proj_pop(0, w), &eigen_global(l, 0), global_size, 1, global_size);
  
  // kappa matrix
  //
  if(Model::bimolecular_size()) {
    //
    Lapack::Matrix kappa = proj_pop.transpose() * inv_proj_bim;

    IO::log << IO::log_offset << "isomers-to-bimolecular equilibrium coefficients (kappa matrix):\n"
      //
	    << IO::log_offset << std::setw(5) << "W\\P";
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::bimolecular(p).short_name();
    
    IO::log << "\n";
   
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(w).short_name();
      
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	//dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w), relax_lave, relax_size) 
	// / well(w).weight_sqrt();
	//
	dtemp = kappa(w, p) / well(w).weight_sqrt();
	
	IO::log << std::setw(Model::log_precision + 7);
	
	if(dtemp < 0.05 && dtemp > -0.05) {
	  //
	  IO::log << "0";
	}
	else
	  //
	  IO::log << dtemp;
      }
      
      IO::log << "\n";
    }

    //    IO::log << std::setprecision(6);
  }

  // bimolecular-to-bimolecular rate coefficients
  //
  if(Model::bimolecular_size()) {
    //
    Lapack::SymmetricMatrix bb_rate = Lapack::SymmetricMatrix(proj_bim.transpose() * inv_proj_bim);
    
    //  Lapack::SymmetricMatrix bb_rate(Model::bimolecular_size());
    //  for(int i = 0; i < Model::bimolecular_size(); ++i)
    //    for(int j = i; j < Model::bimolecular_size(); ++j)
    //	bb_rate(i, j) = triple_product(&eigen_bim(chem_size, i), &eigen_bim(chem_size, j), 
    //				       relax_lave, relax_size);

    for(int i = 0; i < Model::bimolecular_size(); ++i)
      //
      if(bimolecular(i).weight() > 0.)
	//
	for(int j = 0; j < Model::bimolecular_size(); ++j) {
	  //
	  // bimolecular reactant loss
	  //
	  if(i == j) {
	    //
	    dtemp = 0.;
	    
	    for(int b = 0; b < Model::outer_barrier_size(); ++b)
	      //
	      if(Model::outer_connect(b).second == i)
		//
		//dtemp += temperature() / 2. / M_PI * outer_barrier(b).weight();
		for(int e = 0; e < outer_barrier(b).size(); ++e)
		  //
		  dtemp += outer_barrier(b).state_number(e) * thermal_factor(e);
	    
	    dtemp /= 2. * M_PI;
	    
	    dtemp -= bb_rate(i, i);
	  }
	  // crossrate
	  //
	  else
	    //
	    dtemp = bb_rate(i, j);

	  dtemp *= energy_step() / bimolecular(i).weight() / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + i, Model::well_size() + j)] = dtemp; 
	  //arr_out;
	}
    
    // bimolecular-to-escape rate coefficients
    //
    if(Model::escape_size()) {
      /*
	Lapack::Matrix proj_escape = global_escape.copy();
	for(int count = 0; count < Model::escape_size(); ++count)
	for(int l = 0; l < chem_size; ++l)
	parallel_orthogonalize(&proj_escape(0, count), &eigen_global(l, 0), global_size, 1, global_size);
    
	Lapack::Matrix escape_bim = proj_escape.transpose() * inv_proj_bim;
      */
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	if(!Model::bimolecular(p).dummy())
	  //
	  for(int e = 0; e < Model::escape_size(); ++e) {
	    //
	    dtemp = triple_product(&eigen_bim(chem_size, p), &eigen_escape(chem_size, e),
				   //
				   relax_lave, relax_size) * energy_step() / bimolecular(p).weight() / bru;
	    
	    //dtemp = escape_bim(e, p) * energy_step() / bimolecular(p).weight() / bru;
	    
	    rate_data[std::make_pair(Model::well_size() + p, Model::well_size() +
				     //
				     Model::bimolecular_size() + e)] = dtemp; 
	  }
    }
  }

  // product energy distributions
  //
  if(ped_out.is_open()) {
    //
    switch(pressure_unit) {
      //
    case BAR:
      //
      ped_out << "pressure[bar]        = " << pressure() / Phys_const::bar << "\n";
      break;
      
    case TORR:
      //
      ped_out << "pressure[torr]       = " << pressure() / Phys_const::tor << "\n";
      break;
      
    case ATM:
      //
      ped_out << "pressure[atm]        = " << pressure() / Phys_const::atm << "\n";
      break;
    }
    
    ped_out << "temperature[K]           = " << (int)std::floor(temperature() / Phys_const::kelv + 0.5) << "\n"
	    << "energy step[1/cm]        = " << energy_step() / Phys_const::incm << "\n"
	    << "maximum energy[kcal/mol] = " << energy_reference() / Phys_const::kcal << "\n\n";

    if(!Model::bimolecular_size()) {
      //
      std::cerr << funame << "no bimolecular products\n";

      throw Error::Logic();
    }
    
    int ener_index_max;

    // bimolecular-to-bimolecular PEDs
    //
    if(ped_pair.size()) {
      //
      ped_out << "bimolecular PEDs:\n";

      // dimensions
      //
      itemp = 0;
      
      for(int pi = 0; pi < ped_pair.size(); ++ pi)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	  //
	  if(Model::outer_connect(b).second == ped_pair[pi].second && outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

      Lapack::Matrix ped(itemp, ped_pair.size());

      ped = 0.;
      
      // PED
      //
      for(std::vector<std::pair<int, int> >::const_iterator pi = ped_pair.begin(); pi != ped_pair.end(); ++pi)
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {// outer barrier cycle
	  //
	  const int& w = Model::outer_connect(b).first;
	  
	  const int& p = Model::outer_connect(b).second;
	    
	  if(p != pi->second)
	    //
	    continue;
	    
	  for(int e = 0; e < outer_barrier(b).size(); ++e)
	    //
	    ped(e, pi - ped_pair.begin()) += global_bim(e + well_shift[w], p) *
	      //
	      triple_product(&eigen_bim(chem_size, pi->first),
			     //
			     &eigen_global(chem_size, e + well_shift[w]),
			     //
			     relax_lave, relax_size);
	
	}// outer barrier cycle
      
      // output
      //
      ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";
      
      std::vector<int> ped_name_size(ped_pair.size());
      
      for(int pi = 0; pi < ped_pair.size(); ++ pi) {
	//
	stemp = " " + Model::bimolecular(ped_pair[pi].first).short_name() + "->"
	  //
	  + Model::bimolecular(ped_pair[pi].second).short_name();

	itemp = Model::ped_precision + 7;
	
	ped_name_size[pi] = stemp.size() > itemp ? stemp.size() : itemp;
	
	ped_out << std::setw(ped_name_size[pi]) << stemp;
      }
      
      ped_out << "\n";

      for(int e = 0; e < ped.size1(); ++e) {
	//
	ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	for(int pi = 0; pi < ped_pair.size(); ++pi)
	  //
	  ped_out << std::setw(ped_name_size[pi]) << ped(e, pi);
	
	ped_out << "\n";
      }
      
      ped_out << "\n";
      //
    }// bimolecular PEDs
    
    // escape PEDs
    //
    if(Model::escape_size()) {
      //
      // hot energies-to-escape PEDs
      //
      if(hot_index.size()) {
	//
	ped_out << "hot-to-escape PEDs:\n";

	std::vector<int> escape_name_size(Model::escape_size());

	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  itemp = Model::ped_precision + 7;

	  escape_name_size[s] = itemp > Model::escape_name(s).size() + 1 ? itemp : Model::escape_name(s).size() + 1;
	}

	// dimensions
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  const int& w = Model::escape_channel(s).first;
	  
	  if(!s || itemp < well(w).size())
	    //
	    itemp = well(w).size();
	}
	
	Lapack::Matrix ped(itemp, Model::escape_size());

	// hot energies cycle
	//
	for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit) {
	  //
	  const int hw = hit->first;
	      
	  const int he = hit->second;

	  ped = 0.;

	  // PED
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    const int& w = Model::escape_channel(s).first;
	    
	    const int& c = Model::escape_channel(s).second;

	    for(int e = 0; e < well(w).size(); ++e)
	      //
	      ped(e, s) = triple_product(&eigen_global(chem_size, he + well_shift[hw]),
					 //
					 &eigen_global(chem_size, e + well_shift[w]),
					 //
					 relax_lave, relax_size)
		//
		* Model::well(w).escape_rate(energy_bin(e), c) * well(w).boltzman_sqrt(e);
	  }
	  
	  // output
	  //
	  ped_out << "well: " << Model::well(hw).short_name() << "  hot energy = "
	      //
		  << energy_bin(he) / Phys_const::kcal << " kcal/mol\n";

	  ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	  for(int s = 0; s < Model::escape_size(); ++s)
	    //
	    ped_out << std::setw(escape_name_size[s]) << Model::escape_name(s);

	  ped_out << "\n";
	
	  for(int e = 0; e < ped.size1(); ++e) {
	      //
	    ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	    for(int s = 0; s < Model::escape_size(); ++s)
	      //
	      ped_out << std::setw(escape_name_size[s]) << ped(e, s);

	    ped_out << "\n";
	  }

	  ped_out << "\n";
	  //
	}// hot energies cycle
	//
      }// hot energies-to-escape PEDs
	
      // bimolecular-to-escape PEDs
      //
      ped_out << "bimolecular-to-escape PEDs:\n";

      // dimensions
      //
      for(int s = 0; s < Model::escape_size(); ++s) {
	//
	const int& w  = Model::escape_channel(s).first;
	  
	if(!s || well(w).size() > itemp)
	  //
	  itemp = well(w).size();
      }
	
      Lapack::Matrix ped(itemp, Model::bimolecular_size() * Model::escape_size());

      ped = 0.;

      // PED
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  const int& w = Model::escape_channel(s).first;

	  const int& c = Model::escape_channel(s).second;

	  for(int e = 0; e < well(w).size(); ++e) {
	    //
	    itemp = s + Model::escape_size() * p;
	    
	    ped(e, itemp) = triple_product(&eigen_bim(chem_size, p),
					   //
					   &eigen_global(chem_size, e + well_shift[w]),
					   //
					   relax_lave, relax_size)
	      //
	      * Model::well(w).escape_rate(energy_bin(e), c) * well(w).boltzman_sqrt(e);
	  }
	}

      // output
      //
      ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

      std::vector<int> ped_name_size(Model::bimolecular_size() * Model::escape_size());
	
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    stemp = " " + Model::bimolecular(p).short_name() + "->" + Model::escape_name(s);

	    itemp = Model::ped_precision + 7;
	      
	    itemp = stemp.size() > itemp ? stemp.size() : itemp;
	      
	    ped_name_size[s + Model::escape_size() * p] = itemp;
	    
	    ped_out << std::setw(itemp) << stemp;
	  }
	    
      ped_out << "\n";

      for(int e = 0; e < ped.size1(); ++e) {
	//
	ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	for(int i = 0; i < ped.size2(); ++i)
	  //
	  ped_out << std::setw(ped_name_size[i]) << ped(e, i);
	  
	ped_out << "\n";
      }
      
      ped_out << "\n";
      //
    }// escape PEDs
  
    // hot energies-to-bimolecular PEDs
    //
    if(hot_index.size()) {
      //
      ped_out << "hot energies PEDs:\n";

      std::vector<int> bim_name_size(Model::bimolecular_size());
	  
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	itemp = Model::bimolecular(p).short_name().size() + 1;

	bim_name_size[p] = itemp > Model::ped_precision + 7 ? itemp : Model::ped_precision + 7;
      }
	    
      // dimensions
      //
      for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	if(!b || outer_barrier(b).size() > itemp)
	  //
	  itemp = outer_barrier(b).size();
      
      Lapack::Matrix ped(itemp, Model::bimolecular_size());

      // PED
      //
      int count = 0;
	
      for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
	//
	ped = 0.;
	
	for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	  //
	  const int& w = Model::outer_connect(b).first;

	  const int& p = Model::outer_connect(b).second;
	  
	  for(int e = 0; e < outer_barrier(b).size(); ++e)
	    //
	    ped(e, p) += triple_product(&eigen_global(chem_size, e + well_shift[w]),
					//
					&eigen_hot(chem_size, count),
					//
					relax_lave, relax_size)

	      * global_bim(e + well_shift[w], p);
	}
	  
	//output
	//
	ped_out << "well: "<< Model::well(hit->first).short_name()
	  //
		<< "   hot energy = " << energy_bin(hit->second) / Phys_const::kcal
	  //
		<< " kcal/mol\n";
	
	ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  ped_out << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
	    
	ped_out << "\n";
	
	for(int e = 0; e < ped.size1(); ++e) {
	  //
	  ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(bim_name_size[p]) << ped(e, p);
	      
	  ped_out << "\n";
	}

	ped_out << "\n";
	//
      }// hot energy cycle
      //
    }// hot energies PEDs
    //
  }// PEDs output

  /*************************************** BOUND SPECIES *******************************************/

  std::vector<int>  group_index;
  Lapack::Matrix m_direct;
  std::vector<double>      weight;
  std::vector<double> real_weight;

  if(chem_size) {
    //
    // projection of the chemical eigenvectors onto the thermal subspace
    //
    Lapack::Matrix pop_chem(Model::well_size(), chem_size);
    
    for(int l = 0; l < chem_size; ++l)
      //
      pop_chem.column(l) = eigen_pop.row(l);

#ifdef DEBUG

    IO::log << IO::log_offset << "orthogonality check starts\n";
      
    IO::log_offset.increase();

    proj = 0.;
      
    for(int l = 0; l < chem_size; ++l) {
      //
      // normalize
      //
      //normalize(&pop_chem(0, l), Model::well_size());
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
      
    IO::log << IO::log_offset << "orthogonality check done\n";
  
    
#endif

    // partitioning wells into equilibrated groups
    //
    if(default_partition.size()) {
      //
      // default reduction scheme
      //
      well_partition = default_partition;

      IO::log << IO::log_offset << "using default reduction scheme, projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    
      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }
    else if(chem_size == Model::well_size()) {
      //
      // no partitioning
      //
      well_partition.resize(Model::well_size());
      
      for(int w = 0; w < Model::well_size(); ++w)
	//
	well_partition[w].insert(w);

      IO::log << IO::log_offset << "projection error = "
	//
	      << (double)chem_size - well_partition.projection(pop_chem) << "\n";
    }
    else {
      //
      // well partitioning
      //
      Group bimolecular_group;
      
      dtemp = well_partition_method(pop_chem, well_partition, bimolecular_group);

      // convert chemical eigenvectors in the new basis
      //
      pop_chem = well_partition.basis().transpose() * pop_chem;
    }

    group_index = well_partition.group_index();
    weight      = well_partition.weight();
    real_weight = well_partition.real_weight();

    // auxiliary output
    //
    IO::aux << "number of species = " << well_partition.size() << "\n";
    
    for(int g = 0; g < well_partition.size(); ++g) {
      //
      std::multimap<double, int> ww_map;
      
      for(Group::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w)
	//
	ww_map.insert(std::make_pair(well(*w).real_weight(),*w));

      for(std::multimap<double, int>::const_reverse_iterator mit = ww_map.rbegin(); mit != ww_map.rend(); ++mit)
	//
	IO::aux << std::setw(15) << Model::well(mit->second).short_name();

      IO::aux << "\n";
      
      for(std::multimap<double, int>::const_reverse_iterator mit = ww_map.rbegin(); mit != ww_map.rend(); ++mit)
	//
	IO::aux << std::setw(15) << mit->first;
      
      IO::aux << "\n";
    }
    
    // output
    //
    if(chem_size != Model::well_size()) {
      //
      IO::log << IO::log_offset << "combined species:\n"
	//
	      << IO::log_offset << std::setw(2) << "#"  << std::setw(Model::log_precision + 7)
	//
	      << "new name" << IO::first_offset
	//
	      << "group\n";
      
      for(int g = 0; g < well_partition.size(); ++g) {
	//
	IO::log << IO::log_offset << std::setw(2) << g << std::setw(Model::log_precision + 7)
	  //
		<< Model::well(group_index[g]).short_name() << IO::first_offset;
	
	for(Group::const_iterator w = well_partition[g].begin(); w != well_partition[g].end(); ++w) {
	  //
	  if(w != well_partition[g].begin())
	    //
	    IO::log << "+";
	  
	  IO::log << Model::well(*w).short_name();
	}
	
	IO::log << "\n";
      }
    }

    m_direct = pop_chem;
    
    Lapack::Matrix m_inverse = m_direct.invert();

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
  
    // well-to-well rate coefficients
    //
    Lapack::Matrix ww_rate(chem_size);
    
    ww_rate = 0.;
    
    //std::cout << funame << "well-to-well rate contributions:\n";

    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j) {
	//
	//std::cout << i << " <---> " << j << ":";
	
	for(int l = 0; l < chem_size; ++l) {
	  //
	  dtemp = m_direct(j, l) * m_inverse(l, i) * eigenval[l];

	  //std::cout << " " << dtemp;
	  
	  ww_rate(i, j) += dtemp;
	}
	//std::cout << "\n";
      }
    
    // well-to-bimolecular rate coefficients
    //
    Lapack::Matrix wb_rate, bw_rate;
    
    if(Model::bimolecular_size()) {
      //
      wb_rate.resize(chem_size, Model::bimolecular_size());
      
      for(int w = 0; w < chem_size; ++w)
	//
	for(int p = 0; p < Model::bimolecular_size(); ++p)
	  //
	  wb_rate(w, p) = vdot(m_inverse.column(w), &eigen_bim(0, p));

      // bimolecular-to-well rate coefficients
      //
      bw_rate.resize(Model::bimolecular_size(), chem_size);
      
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	for(int w = 0; w < chem_size; ++w)
	  //
	  bw_rate(p, w) = vdot(m_direct.row(w), &eigen_bim(0, p));
    }

    // output
    //
    IO::log << IO::log_offset << "Wa->Wb/Wb->Wa rate constants ratios:\n"
      //
	    << IO::log_offset << std::setw(5) << "Wb\\Wa";
    
    for(int i = 0; i < chem_size; ++i)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::well(group_index[i]).short_name();
    
    IO::log << "\n";
    
    for(int j = 0; j < chem_size; ++j) {
      //
      IO::log << IO::log_offset << std::setw(5) << Model::well(group_index[j]).short_name();
      
      for(int i = 0; i < chem_size; ++i)
	//
	if(i != j) {
	  //
	  if(ww_rate(j, i) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << ww_rate(i, j) / ww_rate(j, i);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	}
	else
	  //
	  IO::log << std::setw(Model::log_precision + 7) << "1";
      
      IO::log << "\n";
    }
    
    if(Model::bimolecular_size()) {
      //
      IO::log << IO::log_offset << "W->P/P->W rate constants ratios:\n"
	//
	      << IO::log_offset << std::setw(5) << "P\\W";
      
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(Model::log_precision + 7) << Model::well(group_index[w]).short_name();
      
      IO::log << "\n";
    
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	IO::log << IO::log_offset << std::setw(5) << Model::bimolecular(p).short_name();
	
	for(int w = 0; w < chem_size; ++w)
	  //
	  if(bw_rate(p, w) != 0.) {
	    //
	    IO::log << std::setw(Model::log_precision + 7) << wb_rate(w, p) / bw_rate(p, w);
	  }
	  else
	    //
	    IO::log << std::setw(Model::log_precision + 7) << "***";
	
	IO::log << "\n";
      }
    }

    //IO::log << std::setprecision(6);

    // well-to-well rate coefficients
    //
    for(int i = 0; i < chem_size; ++i)
      //
      for(int j = 0; j < chem_size; ++j) {
	//
	dtemp = ww_rate(i, j) * std::sqrt(weight[i] * weight[j]) * energy_step() / real_weight[i] / Phys_const::herz;
	
	if(i != j)
	  //
	  dtemp = -dtemp;
	
	rate_data[std::make_pair(group_index[i], group_index[j])] = dtemp; 	
      }

    // well-to-escape rate coefficients
    //
    if(Model::escape_size())
      //
      for(int w = 0; w < chem_size; ++w)
	//
	for(int e = 0; e < Model::escape_size(); ++e) {
	  //
	  dtemp = vdot(m_inverse.column(w), &eigen_escape(0, e))
	    //
	    * std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	  
	  rate_data[std::make_pair(group_index[w], Model::well_size() +
				   //
				   Model::bimolecular_size() + e)] = dtemp;
	}

    // output of steady state distributions associated with individual wells
    //
    IO::aux << "steady state distributions:\n";
    
    for(int g = 0; g < chem_size; ++g) {
      //
      IO::aux << "Initial well: " << Model::well(group_index[g]).short_name() << "\n";

      IO::aux  << std::setw(13) << "E, kcal/mol";

      for(int w = 0; w < Model::well_size(); ++w)
	//
	IO::aux << std::setw(13) << Model::well(w).short_name();

      IO::aux << "\n";

      for(int i = 0; i < well_size_max; ++i) {
	//
	IO::aux << std::setw(13) << energy_bin(i) / Phys_const::kcal;

	for(int w = 0; w < Model::well_size(); ++w) {
	  //
	  if(i < well(w).size()) {
	    //
	    IO::aux << std::setw(13) << vdot(m_inverse.column(g), &eigen_global(0, i + well_shift[w]))
	      //
	      * well(w).boltzman_sqrt(i) / std::sqrt(weight[g]);
	  }
	  else
	    //
	    IO::aux<< std::setw(13) << "0";
	}

	IO::aux << "\n";
      }

      IO::aux << "\n";
    }
    
    // well-to-bimolecular rate coefficients
    //
    for(int w = 0; w < chem_size; ++w)
      //
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	dtemp = wb_rate(w, p) * std::sqrt(weight[w]) * energy_step() / real_weight[w] / Phys_const::herz;
	
	rate_data[std::make_pair(group_index[w], Model::well_size() + p)] = dtemp;
      }
    
    // bimolecular-to-well rate coefficients
    //
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      if(bimolecular(p).weight() > 0.)
	//
	for(int w = 0; w < chem_size; ++w) {
	  //
	  dtemp = bw_rate(p, w) * std::sqrt(weight[w]) * energy_step() / bimolecular(p).weight() / bru;
	  
	  rate_data[std::make_pair(Model::well_size() + p, group_index[w])] = dtemp;
	}

    // product energy distributions
    //
    if(ped_out.is_open()) {
      //
      std::vector<int> bim_name_size(Model::bimolecular_size());
	  
      for(int p = 0; p < Model::bimolecular_size(); ++p) {
	//
	itemp = Model::ped_precision + 7;

	bim_name_size[p] = Model::bimolecular(p).short_name().size() + 1 > itemp ?  Model::bimolecular(p).short_name().size() + 1 : itemp;
      }
      
      // well-to-bimolecular distributions
      //
      if(Model::bimolecular_size()) {
	//
	ped_out << "well PEDs:\n";

	// dimensions
	//
	for(int b = 0; b < Model::outer_barrier_size(); ++b)
	//
	  if(!b || outer_barrier(b).size() > itemp)
	    //
	    itemp = outer_barrier(b).size();

	Lapack::Matrix ped(itemp, Model::bimolecular_size());

	// PED
	//
	for(int c = 0; c < chem_size; ++c) {
	  //
	  ped = 0.;
	  
	  for(int b = 0; b < Model::outer_barrier_size(); ++b) {
	    //
	    const int& w = Model::outer_connect(b).first;
	      
	    const int& p = Model::outer_connect(b).second;
	      
	    for(int e = 0; e < outer_barrier(b).size(); ++e) 
	    //
	      ped(e, p) += vdot(m_inverse.column(c), &eigen_global(0, e + well_shift[w]))
		//
		* global_bim(e + well_shift[w], p);
	  }

	  // output
	  //
	  ped_out << "well: " << Model::well(group_index[c]).short_name() << "\n";
	  
	  ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	  for(int p = 0; p < Model::bimolecular_size(); ++p)
	    //
	    ped_out << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
	  
	  ped_out << "\n";

	  for(int e = 0; e < ped.size1(); ++e) {
	    //
	    ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	    for(int p = 0; p < Model::bimolecular_size(); ++p)
	      //
	      ped_out << std::setw(bim_name_size[p]) << ped(e, p);
	    
	    ped_out << "\n";
	  }
	  
	  ped_out << "\n";
	  //
	}//well cycle
	//
      }// well PEDs   

      // well-to-escape PEDs
      //
      if(Model::escape_size()) {
	//
	ped_out << "well-to-escape PEDs:\n";

	// dimensions
	//
	for(int s = 0; s < Model::escape_size(); ++s) {
	  //
	  const int& w = Model::escape_channel(s).first;
	  
	  if(!s || well(w).size() > itemp)
	    //
	    itemp = well(w).size();
	}
	
	Lapack::Matrix ped(itemp, chem_size * Model::escape_size());

	ped = 0.;
	
	// PED
	//
	for(int c = 0; c < chem_size; ++c)
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    const int& w = Model::escape_channel(s).first;
	      
	    const int& d = Model::escape_channel(s).second;
	      
	    for(int e = 0; e < well(w).size(); ++e) {
	      //
	      itemp = s + Model::escape_size() * c;
	  
	      ped(e, itemp) = vdot(m_inverse.column(c), &eigen_global(0, e + well_shift[w]))
		//
		* Model::well(w).escape_rate(energy_bin(e), d) * well(w).boltzman_sqrt(e);
	    }
	  }

	// output
	//
	ped_out << std::setw(Model::ped_precision + 7) << "E, kcal";

	std::vector<int> ped_name_size(chem_size * Model::escape_size());
	
	for(int c = 0; c < chem_size; ++c)
	  //
	  for(int s = 0; s < Model::escape_size(); ++s) {
	    //
	    stemp = " " + Model::well(group_index[c]).short_name() + "->" + Model::escape_name(s);

	    itemp = Model::ped_precision + 7;

	    itemp = itemp > stemp.size() ? itemp : stemp.size();

	    ped_name_size[s + Model::escape_size() * c] = itemp;
	    
	    ped_out << std::setw(itemp) << stemp;
	  }
	
	ped_out << "\n";
      
	for(int e = 0; e < ped.size1(); ++e) {
	  //
	  ped_out << std::setw(Model::ped_precision + 7) << energy_bin(e) / Phys_const::kcal;

	  for(int i = 0; i < ped.size2(); ++i)
	    //
	    ped_out << std::setw(ped_name_size[i]) << ped(e, i);
	  
	  ped_out << "\n";
	}
	
	ped_out << "\n";
	//
      }// well-to-escape PEDs
      //
    }// product energy distributions
    //
  }// bound species

  // hot distribution branching ratios
  //
  if(hot_index.size()) {
    //
    int well_name_size_max;

    std::vector<int> well_name_size(Model::well_size());
    
    for(int w = 0; w < Model::well_size(); ++w) {
      //
      itemp = Model::well(w).short_name().size() + 1;

      if(!w || itemp > well_name_size_max)
	//
	well_name_size_max = itemp;
    
      if(itemp < Model::log_precision + 7)
	//
	itemp = Model::log_precision + 7;

      well_name_size[w] = itemp;
    }

    if(well_name_size_max < 5)
      //
      well_name_size_max = 5;

    IO::log << IO::log_offset << "hot energies branching fractions:\n"
      //
	    << IO::log_offset //<< std::setprecision(6)
      //
	    << std::setw(well_name_size_max)  << "Well"
      //
	    << std::setw(Model::log_precision + 7) << "E, kcal";
    
    for(int w = 0; w < chem_size; ++w)
      //
      IO::log << std::setw(well_name_size[group_index[w]]) << Model::well(group_index[w]).short_name();

    std::vector<int> bim_name_size(Model::bimolecular_size());
    
    for(int p = 0; p < Model::bimolecular_size(); ++p) {
      //
      itemp = Model::bimolecular(p).short_name().size() + 1;

      if(itemp < Model::log_precision + 7)
	//
	itemp = Model::log_precision + 7;

      bim_name_size[p] = itemp;
      
      IO::log << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
    }

    for(int e = 0; e < Model::escape_size(); ++e)
      //
      IO::log << std::setw(Model::log_precision + 7) << Model::escape_name(e);
    
    IO::log << "\n";
    
    int count = 0;
    
    for(std::set<std::pair<int, int> >::const_iterator hit = hot_index.begin(); hit != hot_index.end(); ++hit, ++count) {
      //
      IO::log << IO::log_offset
	//
	      << std::setw(well_name_size_max)  << Model::well(hit->first).short_name()
	//
	      << std::setw(Model::log_precision + 7) << energy_bin(hit->second) / Phys_const::kcal;
	
      for(int w = 0; w < chem_size; ++w)
	//
	IO::log << std::setw(well_name_size[group_index[w]]) << vdot(m_direct.row(w), &eigen_hot(0, count)) * std::sqrt(weight[w]);
	
      for(int p = 0; p < Model::bimolecular_size(); ++p)
	//
	IO::log << std::setw(bim_name_size[p]) << triple_product(&eigen_bim(chem_size, p), &eigen_hot(chem_size, count),
								 //
								 relax_lave, relax_size);
	
      for(int e = 0; e < Model::escape_size(); ++e)
	//
	IO::log << std::setw(Model::log_precision + 7) << triple_product(&eigen_escape(chem_size, e), &eigen_hot(chem_size, count),
									 //
									 relax_lave, relax_size);
	
      IO::log << "\n";
    }
  }

  // prompt isomerization
  //
  IO::log << IO::log_offset << "prompt isomerization/dissociation:\n";

  std::vector<int> well_name_size(Model::well_size());

  int well_name_size_max;

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    itemp = Model::well(w).short_name().size() + 1;

    if(!w || itemp > well_name_size_max)
      //
      well_name_size_max = itemp;
    
    if(itemp < Model::log_precision + 7)
      //
      itemp = Model::log_precision + 7;

    well_name_size[w] = itemp;
  }

  if(well_name_size_max < 7)
    //
    well_name_size_max = 7;

  IO::log << IO::log_offset << std::setw(well_name_size_max) << "W\\W,P,E";

  if(chem_size)
    //
    for(int w = 0; w < Model::well_size(); ++w)
      //
      IO::log << std::setw(well_name_size[w]) << Model::well(w).short_name();

  IO::log << std::setw(Model::log_precision + 7) << "Total";

  std::vector<int> bim_name_size(Model::bimolecular_size());
  
  for(int p = 0; p < Model::bimolecular_size(); ++p) {
    //
    itemp = Model::bimolecular(p).short_name().size() + 1;

    if(itemp < Model::log_precision + 7)
      //
      itemp = Model::log_precision + 7;

    bim_name_size[p] = itemp;
    
    IO::log << std::setw(bim_name_size[p]) << Model::bimolecular(p).short_name();
  }
  
  for(int e = 0; e < Model::escape_size(); ++e)
    //
    IO::log << std::setw(Model::log_precision + 7) << Model::escape_name(e);
    
  IO::log << "\n";

  for(int w = 0; w < Model::well_size(); ++w) {
    //
    IO::log << IO::log_offset << std::setw(well_name_size_max) << Model::well(w).short_name();

    double diss = 1.;

    if(chem_size)
      //
      for(int v = 0; v < Model::well_size(); ++v) {
	//
	dtemp = 0.;

	for(int l = 0; l < chem_size; ++l)
	  //
	  dtemp += eigen_pop(l, w) * eigen_pop(l, v);

	// renormalization
	//
	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	dtemp *= well(v).weight_sqrt() * well(w).weight_sqrt() * energy_step() / well(w).real_weight();

	if(v == w)
	  //
	  dtemp = 1. - dtemp;
      
	diss -= dtemp;

	IO::log << std::setw(well_name_size[v]) << dtemp;
      }

    IO::log << std::setw(Model::log_precision + 7) << diss;
    
    for(int p = 0; p < Model::bimolecular_size(); ++p)
      //
      IO::log << std::setw(bim_name_size[p]) << triple_product(&eigen_bim(chem_size, p), &eigen_pop(chem_size, w),
							       //
							       relax_lave, relax_size) * energy_step()
	//
	* well(w).weight_sqrt() / well(w).real_weight();

    for(int e = 0; e < Model::escape_size(); ++e)
      //
      IO::log << std::setw(Model::log_precision + 7) << triple_product(&eigen_escape(chem_size, e), &eigen_pop(chem_size, w),
								       //
								       relax_lave, relax_size) * energy_step()
	//
	* well(w).weight_sqrt() / well(w).real_weight();
      
    IO::log << "\n";
  }

  IO::log << "\n";
}

