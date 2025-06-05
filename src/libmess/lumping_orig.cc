lump_t MasterEquation::lumping_scheme (double temperature, double pressure, const Model::ChemGraph& graph, std::pair<int, int> diss)
{
  const char funame [] = "MasterEquation::lumping_scheme: ";

  using Model::operator<<;
  
  int    itemp;
  double dtemp;
  bool   btemp;

  lump_t res;

  if(graph.well_set.size() == 1 && diss.first == Model::INIT) {
    //
    res.first.push_back(graph.well_set);

    return res;
  }

  if(!graph.is_connected()) {
    //
    std::cerr << funame << "graph is not connected: " << graph << "\n";

    throw Error::Logic();
  }
  
  IO::Marker funame_marker(funame);
  
  IO::log << IO::log_offset << "well group: " << graph.well_set << "\n";

  double be;

  int sb = -1;
  
  switch(diss.first) {
    //
  case Model::INNER:
    //
    be = Model::inner_barrier(diss.second).dist_ener_max(temperature);
    
    IO::log << IO::log_offset << "(inner) dissociation barrier " << Model::inner_barrier(diss.second).name();
    
    break;
    //
  case Model::OUTER:
    //
    be = Model::outer_barrier(diss.second).dist_ener_max(temperature);
    
    IO::log << IO::log_offset << "(outer) dissociation barrier " << Model::outer_barrier(diss.second).name();
    
    break;
    //
  case Model::INIT:
    //
    sb = graph.split(temperature);
    
    be = Model::inner_barrier(sb).dist_ener_max(temperature);
    
    IO::log << IO::log_offset << "well-splitting barrier " << Model::inner_barrier(itemp).name();
    
    break;
    //
  default:
    //
    std::cerr << funame << "wrong dissociation barrier type: " << diss.first << "\n";
    
    throw Error::Logic();
  }
  
  IO::log << ": distribution maximum energy = " << std::ceil(be / Phys_const::kcal * 10.) / 10. << " kcal/mol\n";
													    
  double er = be + excess_energy_over_temperature * temperature;
    
  double ec = be - energy_cutoff_over_temperature * temperature;
    
  ReactiveComplex rc;

  rc.init(temperature, pressure, er, ec, graph, ReactiveComplex::LUMPING);

  if(diss.first != Model::INIT)
    //
    rc.set_dissociation_channel(diss);

  rc.set_kinetic_matrix();
  
  ReactiveComplex::RateData rate_data = rc.well_reduction_method();

  const std::vector<std::set<int> >& wp = rate_data.well_partition;

  const std::set<int>& bg = rate_data.bimolecular_group;

  if(wp.size() == 1 && !bg.size() && diss.first == Model::INIT) {
    //
    res.first.push_back(wp[0]);

    return res;
  }
  
  for(int g = 0; g < wp.size(); ++g) {
    //
    if(wp[g].size() == 1) {
      //
      res.first.push_back(wp[g]);

      continue;
    }

    if(sb >= 0 &&
       //
       wp[g].find(Model::inner_connect(sb).first) != wp[g].end() &&
       //
       wp[g].find(Model::inner_connect(sb).second) != wp[g].end()) {

      res.first.push_back(wp[g]);

      continue;
    }

    Model::ChemGraph cg(wp[g]);

    if(cg.is_connected()) {
      //
      lump_t ls = lumping_scheme(temperature, pressure, cg);

      for(std::list<std::set<int> >::const_iterator li = ls.first.begin(); li != ls.first.end(); ++li)
	//
	res.first.push_back(*li);
      
      for(std::set<int>::const_iterator w = ls.second.begin(); w != ls.second.end(); ++w)
	//
	res.second.insert(*w);

      continue;
    }
      
    cg = graph;

    btemp = true;

    while(btemp) {
      //
      std::list<Model::ChemGraph> gl;

      cg.split(temperature, &gl);
	
      btemp = false;

      for(std::list<Model::ChemGraph>::const_iterator li = gl.begin(); li != gl.end(); ++li)
	//
	if(li->does_include(wp[g])) {
	  //
	  cg = *li;

	  btemp = true;

	  break;
	}
    }

    if(cg.well_set.size() == graph.well_set.size() && diss.first == Model::INIT) {
      //
      IO::log << IO::log_offset << "WARNING: it seems that we are in the inifinite loop: assuming that the "
	//
	      << wp[g] << " patition group of the " << graph << " reactive complex is unsplitable\n";
								  
      res.first.push_back(wp[g]);

      continue;
    }

    /*
      btemp = false;
      
      for(std::set<int>::const_iterator w = cg.well_set.begin(); w != cg.well_set.end(); ++w) {
      //
      if(wp[g].find(*w) == wp[g].end()) {
      //
      btemp = true;
	  
      IO::log << IO::log_offset << "WARNING: " << std::setw(4) << Model::well(*w).name()
      //
      << " well is not found in the original group";
											
      if(bg.find(*w) == bg.end())
      //
      IO::log << " and in the bimolecular group";

      IO::log << "\n";
      }
      }

      if(btemp)
    */
    IO::log << IO::log_offset << "reactive complex:  " << graph << "\n"
      //
	    << IO::log_offset << "bimolecular group: " << bg << "\n"
      //
	    << IO::log_offset << "original partition group: " << wp[g] << "\n"
      //
	    << IO::log_offset << "extended partition group: " << cg.well_set << "\n";

    lump_t ls = lumping_scheme(temperature, pressure, cg);

    for(std::list<std::set<int> >::const_iterator li = ls.first.begin(); li != ls.first.end(); ++li) {
      //
      std::set<int> wg;

      for(std::set<int>::const_iterator w = li->begin(); w != li->end(); ++w)
	//
	if(wp[g].find(*w) != wp[g].end())
	  //
	  wg.insert(*w);

      if(wg.size())
	//
	res.first.push_back(wg);
    }
      
    for(std::set<int>::const_iterator w = ls.second.begin(); w != ls.second.end(); ++w)
      //
      if(wp[g].find(*w) != wp[g].end())
	//
	res.second.insert(*w);
  }
  
  for(std::set<int>::const_iterator w = bg.begin(); w != bg.end(); ++w)
    //
    res.second.insert(*w);
  
  return res;
}
