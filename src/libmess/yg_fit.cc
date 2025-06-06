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

#include "key.hh"
#include "io.hh"
#include "system.hh"
#include "structure.hh"
#include "potential.hh"

#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <list>

#include <fstream>
#include <sstream>

double cap_ener = -1.;

void ener_cap (double& ener)
{
  double dtemp;

  if(-ener > cap_ener) {
    //
    dtemp = -ener / cap_ener - 1.;

    if(dtemp > 0.1)
      //
      ener = -cap_ener * (1. + dtemp * (1. - std::exp(-1. / dtemp)));
  }
  else if(ener > cap_ener) {
    //
    dtemp = ener / cap_ener - 1.;

    if(dtemp > 0.1)
      //
      ener = cap_ener * (1. + dtemp * (1. - std::exp(-1. / dtemp)));
  }
}

class InputError {};

class Index : private std::vector<int> {
  //
  int _order;

  bool _end;

  void _assert () const;
  
public:
  //
  Index (int, int);

  bool end () const { return _end; }

  void operator++ ();

  void operator++ (int) { operator++(); }

  vector<int> value () const; 
};

Index::Index (int n, int ord) : std::vector<int>(n - 1, 0), _order(ord), _end(false)
  {
  const char funame [] = "Index::Index: ";

  if(n < 2 || ord < 1) {
    //
    std::cerr << funame << "out of range: " << n << ", " << ord <<"\n";

    throw Error::Range();
  }
}

void Index::_assert () const
{
  const char funame [] = "Index::_assert: ";
  
  if(_end) {
    //
    std::cerr << funame << "out of range\n";

    throw Error::Range();
  }
}
  
void Index::operator++ ()
{
  int  itemp;

  _assert();
  
  for(int i = 0; i < size(); ++i) {
    //
    if(++(*this)[i] == _order + 1)
      //
      continue;

    for(int j = 0; j < i; ++i)
      //
      (*this)[j] = (*this)[i];

    return;
  }

  _end = true;
}

std::vector<int> Index::value () const
{
  _assert();

  std::vector<int> res(size() + 1);

  for(int i = 0; i < res.size(); ++i)
    //
    if(!i) {
      //
      res[i] = _order - (*this)[i];
    }
    else if(i == size()) {
      //
      res[i] = (*this)[i - 1];
    }
    else
      //
      res[i] =(*this)[i - 1] - (*this)[i];

  return res;
}      
  
int main (int argc, char* argv[])
{
  const char funame [] = "yg_fit: ";

  double         dtemp;
  int            itemp;
  int            btemp;
  std::string    stemp;
  Lapack::Vector vtemp;
  Lapack::Matrix mtemp;

  Structure::mute = 0;

  std::cout.precision(3);
  
  // usage
  //
  if(argc != 2) {
    //
    std::cerr << "usage: yg_fit input_file\n";

    return 0;
  }

  std::ifstream from(argv[1]);

  if(!from) {
    //
    std::cerr << funame << "cannot open " << argv[1] << " file\n";

    return 0;
  }

  /*********************************************************************************
   ******************************* INPUT PARAMETERS ********************************
   *********************************************************************************/
  //
  std::string smp_file, out_file, scat_file;

  int order_min = -1, order_max = -1;
  
  double svd_prec = -1., bf_tol = 1.e-7, ener_decr = -1.;

  double anchor_dist = -1., ener_min = -1., ener_tol = -1.;
  
  KeyGroup YGfitingGroup;

  Key   struc_key("Structure"        );
  Key     pot_key("Potential"        );
  Key     ord_key("OrderRange"       );
  Key     smp_key("SmpFile"          );
  Key     out_key("OutFile"          );
  Key    scat_key("ScatterFile"      );
  Key   danch_key("AnchorDist[bohr]" );
  Key     svd_key("SVDPrecision"     );
  Key     bft_key("BasisFunTol"      );
  Key    edec_key("EnerDec[kcal/mol]");
  Key    emin_key("EnerMin[kcal/mol]");
  Key    ecap_key("EnerCap[kcal/mol]");
  Key    etol_key("EnerTol"          );
  
  std::string token, comment;

  IO::LineInput lin;
      
  while(from >> token) {
    //
    // molecular structure initialization
    //
    if(struc_key == token) {
      //
      token += ": ";

      std::getline(from, comment);

      if(Structure::isinit()) {
	//
	std::cerr << funame << token << "already initialized\n";
	
	throw InputError();
      }
      
      try {
	//
	Structure::init(from);
      }
      catch(Error::General) {
	//
	throw InputError();
      }
    }
    // potential initialization
    //
    else if(pot_key == token) {
      //
      token += ": ";

      if(Potential::YG::pot) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }
      
      try {
	//
	Potential::YG::init(from);
      }
      catch(Error::General) {
	//
	throw InputError();
      }
    }
    // output file
    //
    else if(out_key == token) {
      //
      token += ": ";

      if(out_file.size()) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> out_file)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }
    }
    // scatter plot file
    //
    else if(scat_key == token) {
      //
      token += ": ";

      if(scat_file.size()) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> scat_file)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }
    }
    // sampling file
    //
    else if(smp_key == token) {
      //
      token += ": ";

      if(smp_file.size()) {
	//
	std::cerr << funame << token << "already initialized\n";
	
	throw InputError();
      }

      lin.read_line(from);
      
      if(!(lin >> smp_file)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }
    }
    // cap energy
    //
    else if(ecap_key == token) {
      //
      token += ": ";

      if(cap_ener > 0.) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      IO::LineInput lin(from);
      
      if(!(lin >> cap_ener)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(cap_ener <= 0.) {
	//
	std::cerr << funame << token << "should be positive\n";
	
	throw InputError();
      }

      //cap_ener *= Phys_const::kcal;
    }
    // test energy minimum
    //
    else if(emin_key == token) {
      //
      token += ": ";

      if(ener_min > 0.) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> ener_min)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(ener_min <= 0.) {
	//
	std::cerr << funame << token << "out of range\n";

	throw InputError();
      }
    }
    // energy fitting tolerance
    //
    else if(etol_key == token) {
      //
      token += ": ";

      if(ener_tol > 0.) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> ener_tol)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(ener_tol <= 0.) {
	//
	std::cerr << funame << token << "out of range\n";

	throw InputError();
      }
    }
    // scatter plot distance
    //
    else if(danch_key == token) {
      //
      token += ": ";

      if(anchor_dist > 0.) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> anchor_dist)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(anchor_dist <= 0.) {
	//
	std::cerr << funame << token << "out of range\n";

	throw InputError();
      }
    }
    // SVD precision
    //
    else if(svd_key == token) {
      //
      token += ": ";

      if(svd_prec > 0.) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> svd_prec)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(svd_prec <= 0.) {
	//
	std::cerr << funame << token << "out of range\n";

	throw InputError();
      }
    }
    // basis function colinearity tolerance
    //
    else if(bft_key == token) {
      //
      token += ": ";

      lin.read_line(from);

      if(!(lin >> bf_tol)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(bf_tol <= 0.) {
	//
	std::cerr << funame << token << "out of range\n";

	throw InputError();
      }
    }
    // basis function energy decrement minimum
    //
    else if(edec_key == token) {
      //
      token += ": ";

      if(ener_decr > 0.) {
	//
	std::cerr << funame << token << "already initialized\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> ener_decr)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }

      if(ener_decr <= 0.) {
	//
	std::cerr << funame << token << "should be positive\n";

	throw InputError();
      }

      //ener_decr *= Phys_const::kcal;
    }
    // order range of the basis functions
    //
    else if(ord_key == token) {
      //
      token += ": ";

      if(!Potential::YG::pot) {
	//
	std::cerr << funame << token << "potential should be initialized first\n";

	throw InputError();
      }

      lin.read_line(from);

      if(!(lin >> order_min >> order_max)) {
	//
	std::cerr << funame << token << "corrupted\n";
	
	throw InputError();
      }
      
      if(order_min <= 0 || order_min > order_max) {
	//
	std::cerr << funame << token << "out of range: " << order_min << ", " << order_max << "\n";
	
	throw InputError();
      }
    }
    // unknown keyword
    //
    else if(IO::skip_comment(token, from)) {
      //
      std::cerr << funame << "unknown keyword " << token << "\n";
      
      Key::show_all(std::cerr);
      
      std::cerr << "\n";

      throw InputError();
    }//
    //
  }// input cycle

  from.close();
  from.clear();
  
  if(!out_file.size()) {
    //
    std::cerr << funame << "coefficients output file not initialized\n";
    
    throw InputError();
  }

  if(!smp_file.size()) {
    //
    std::cerr << funame << "training data file not initialized\n";

    throw InputError();
  }

  // order range
  //
  if(order_min < 0) {
    //
    std::cerr << funame << "order range not initialized\n";

    throw InputError();
  }
      
  // cap energy
  //
  if(cap_ener < 0.) {
    //
    std::cerr << funame << "energy cap not initialized\n";

    throw InputError();
  }
      
  // angular vector dimension
  //
  switch(Structure::type(1)) {
    //
  case Molecule::MONOATOMIC:
    //
    itemp = 0;

    break;
    //
  case Molecule::LINEAR:
    //
    itemp = 3;

    break;
    //
  case Molecule::NONLINEAR:
    //
    itemp = 4;
  }

  // angular vector
  //
  Array<double> ang_pos(3 + itemp, 0.);

  /******************************************************************************************
   ******************************* TRAINING SAMPLING FILE ***********************************
   ******************************************************************************************/
  //
  from.open(smp_file.c_str());

  if(!from) {
    //
    std::cerr << funame << "cannot open sampling file: " << smp_file << "\n";

    throw InputError();
  }

  // header line
  //
  int ray_size, dist_size, ang_size;
  
  lin.read_line(from);

  if(!(lin >> ray_size >> dist_size >> ang_size)) {
    //
    std::cerr << funame << "cannot read sampling dimensions\n";

    throw InputError();
  }

  if(ray_size <= 0 || dist_size <= 0 || ang_size != ang_pos.size()) {
    //
    std::cerr << funame << "sampling dimensions our of range: " << ray_size << ", "
      //
	      << dist_size << ", " << ang_size << " vs " << ang_pos.size() << "\n";

    throw InputError();
  }
    
  // distances data section
  //
  std::vector<double> dist_data(dist_size);

  itemp = -1;
  
  for(int d = 0; d < dist_data.size(); ++d) {
    //
    if(!(from >> dtemp)) {
      //
      std::cerr << funame << d << "-th distance: reading failed\n";

      throw InputError();
    }
      
    if(dtemp <= 0.) {
      //
      std::cerr << funame << d << "-th distance: out of range: " << dtemp << "\n";

      throw InputError();
    }

    if(itemp < 0 && dtemp < anchor_dist)
      //
      itemp = d;
    
    dist_data[d] = dtemp;
  }

  if(itemp < 0)
    //
    itemp = dist_data.size() - 1;

  const int anchor_index = itemp;

  // angular vectors and energies
  //
  std::vector<Array<double> > ang_data(ray_size,   ang_pos);
  Lapack::Matrix             ener_data(ray_size, dist_size);
  std::vector<int>           ener_size(ray_size, -1);

  int count = 0;

  int ray_index;
  
  while(from >> ray_index) {
    //
    if(ray_index < 0 || ray_index >= ray_size) {
      //
      std::cerr <<  "ray sampling index out of range: " << ray_index << " vs " << ray_size << "\n";

      throw InputError();
    }
    
    count++;
    
    lin.read_line(from);

    for(int i = 0; i < ang_size; ++i)
      //
      lin >> ang_data[ray_index][i];

    if(!lin) {
      //
      std::cerr << ray_index << "-th ray sampling: reading angular position failed\n";

      throw InputError();
    }

    lin.read_line(from);

    if(!(lin >> itemp)) {
      //
      std::cerr << ray_index << "-th ray sampling: cannot read energies #\n";

      throw InputError();
    }

    if(itemp < 0 || itemp > dist_size) {
      //
      std::cerr << ray_index << "-th ray sampling: energies # out of range: " << itemp << "\n";

      throw InputError();
    }

    ener_size[ray_index] = itemp;
    
    // energies
    //
    for(int e = 0; e < itemp; ++e) {
      //
      if(!(from >> dtemp)) {
	//
	std::cerr << ray_index << "-th ray sampling: " << e << "-the energy reading failed\n";

	throw InputError();
      }

      ener_cap(dtemp);
      
      ener_data(ray_index, e) = dtemp;
    }
  }
  
  from.close();
  from.clear();

  if(count != ray_size) {
    //
    std::cerr << "missing ray samplings: " << count << " vs " << ray_size << "\n";

    throw InputError();
  }
  
  /*************************************************************************************
   ********************************* BASIS FUNCTIONS ***********************************
   *************************************************************************************/
  //
  const std::vector<int> index_range = Potential::YG::pot->index_range();

  itemp = 0;
  
  for(int i = 0; i < index_range.size(); ++i)
    //
    if(index_range[i] <= 0)
      //
      ++itemp;

  const int unbound_size = itemp;

  if(unbound_size < 2) {
    //
    std::cerr << funame << "# of unbound indices out of range: " << unbound_size << "\n";

    throw Error::Range();
  }
  
  std::vector<int> curr_bfi(index_range.size(), 0);

  // fixed order bisis function index set
  //
  std::set<std::vector<int> > fix_bfi_set; 

  std::set<std::vector<int> > zero_bfi_set; 

  // total basis function index set
  //
  std::set<std::vector<int> > tot_bfi_set;

  Lapack::Matrix ang_bfi_mat;

  Lapack::Vector ener(ray_size);

  itemp = 0;
  
  for(int r = 0; r < ray_size; ++r)
    //
    if(ener_size[r] > anchor_index)
      //
      ener[itemp++] = ener_data(r, anchor_index);

  ener.resize(itemp);
  
  vtemp.resize(itemp);

  double smp_size_sqrt = std::sqrt((double)itemp);

  std::cout << "test dist [bohr] = " << dist_data[anchor_index] <<  ",  samplings size = " << itemp << "\n";

  // basis functions
  //
  for(int order = 0; order < order_max; ++order) {
    //
    // zero-th order basis function indices
    //
    if(!order) {
      //
      itemp = 0;

      while(itemp < index_range.size()) {
	//
	fix_bfi_set.insert(curr_bfi);

	zero_bfi_set.insert(curr_bfi);
	
	for(itemp = 0; itemp < index_range.size(); ++itemp)
	  //
	  if(index_range[itemp] > 0) {
	    //
	    if(++curr_bfi[itemp] == index_range[itemp])
	      //
	      continue;

	    for(int i = 0; i < itemp; ++i) 
	      //
	      curr_bfi[i] = 0;

	    break;
	  }
      }
    }
    // non-zero fixed order basis function indices
    //
    else {
      //
      fix_bfi_set.clear();
      
      for(Index ind(unbound_size, order); !ind.end(); ++ind) {
	//
	std::vector<int> ival = ind.value();	
	
	for(std::set<std::vector<int> >::const_iterator bit = zero_bfi_set.begin(); bit != zero_bfi_set.end(); ++bit) {
	  //
	  curr_bfi = *bit;

	  itemp = 0;
	  
	  for(int i = 0; i < index_range.size(); ++i)
	    //
	    if(index_range[i] <= 0)
	      //
	      curr_bfi[i] = ival[itemp++];

	  fix_bfi_set.insert(curr_bfi);
	}
      }
    }
    
    // basis function test
    //
    std::list<double> er_list;

    for(std::set<std::vector<int> >::const_iterator bit = fix_bfi_set.begin(); bit != fix_bfi_set.end(); ++bit) {
      //
      itemp = 0;
      
      for(int r = 0; r < ray_size; ++r)
	//
	if(ener_size[r] > anchor_index)
	  //
	  vtemp[itemp++] = Potential::YG::pot->basis_function(*bit, ang_data[r]);

      for(int i = 0; i < ang_bfi_mat.size2(); ++i)
	//
	::orthogonalize(vtemp, &ang_bfi_mat(0, i), vtemp.size());

      dtemp = ::normalize(vtemp);

      if(tot_bfi_set.size() && dtemp < bf_tol) {
	//
	std::cout << "linear-dependent basis function:";

	for(int i = 0; i < bit->size(); ++i)
	  //
	  std::cout << "  " << (*bit)[i];

	std::cout << "\n";
	
	continue;
      }
      
      if(order >= order_min && std::fabs(vdot(vtemp, ener)) / smp_size_sqrt < ener_decr)
	//
	continue;
      
      orthogonalize(ener, vtemp);

      er_list.push_back(vlength(ener) / smp_size_sqrt);

      // add basis function to the set
      //
      mtemp.resize(ener.size(), tot_bfi_set.size() + 1);

      for(int i = 0; i < tot_bfi_set.size(); ++i)
	//
	mtemp.column(i) = ang_bfi_mat.column(i);

      mtemp.column(tot_bfi_set.size()) = vtemp;

      ang_bfi_mat = mtemp.copy();

      tot_bfi_set.insert(*bit);
    }
    
    std::cout << "order = " << std::setw(3) << order << ", rmsd [kcal/mol]:";

    if(er_list.size()) {
      //
      for(std::list<double>::const_iterator rit = er_list.begin(); rit != er_list.end(); ++rit)
	//
	std::cout << "  " << *rit;
    }
    else
      //
      std::cout << "none";

    std::cout << std::endl;
    //
  }// order cycle

  /************************************************************************************************
   ********************************** BASIS FUNCTION COEFFICIENTS *********************************
   ************************************************************************************************/
  //
  // scattering plot data stream
  //
  std::ofstream scat_out;
    
  if(scat_file.size())
    //
    scat_out.open(scat_file.c_str());

  // large energy deviation ray samplings
  //
  std::map<double, std::pair<int, int> > check_ray;

  std::vector<double> ener_rmsd(dist_size);

  // distance cycle
  //
  for(int d = 0; d < dist_size; ++d) {
    //
    ener.resize(ray_size);

    std::vector<int> n2o(ray_size);
    
    itemp = 0;
    
    for(int r = 0; r < ray_size; ++r)
      //
      if(ener_size[r] > d) {
	//
	n2o[itemp] = r;
	
	ener[itemp++] = ener_data(r, d);
      }

    n2o.resize(itemp);
    
    ener.resize(itemp);

    ang_bfi_mat.resize(itemp, tot_bfi_set.size());

    count = 0;
    
    for(std::set<std::vector<int> >::const_iterator bit = tot_bfi_set.begin(); bit != tot_bfi_set.end(); ++bit, ++count) {
      //
      itemp = 0;
      
      for(int r = 0; r < ray_size; ++r)
	//
	if(ener_size[r] > d)
	  //
	  ang_bfi_mat(itemp++, count) = Potential::YG::pot->basis_function(*bit, ang_data[r]);
    }
    
    Lapack::Vector coef = Lapack::svd_solve(ang_bfi_mat, ener, svd_prec);

    Lapack::Vector fit = ang_bfi_mat * coef;

    // large energy deviation ray samplings
    //
    if(ener_min > 0. && ener_tol > 0.)
      //
      for(int i = 0; i < ener.size(); ++i)
	//
	if(std::fabs(ener[i]) > ener_min) {
	  //
	  dtemp = std::fabs(fit[i] / ener[i] - 1.);

	  if(dtemp > ener_tol) {
	    //
	    if(check_ray.find(dtemp) != check_ray.end())
	      //
	      std::cerr << "WARNING: " << check_ray[dtemp].first << "-th ray sampling has the same relative deviation as " << n2o[i] << "-th one\n";
	    
	    check_ray[dtemp] = std::make_pair(n2o[i], d);
	  }
	}
    
    // scatter plot
    //
    if(scat_out) {
      //
      std::cout << "distance [bohr] = " << dist_data[d] << " energies # = " << ener.size() << "\n";
      
      for(int i = 0; i < ener.size(); ++i)
	//
	scat_out << ener[i] << "  " << fit[i] << "\n";
    }

    // rmsd
    //
    double res = 0.;
    
    for(int i = 0; i < ener.size(); ++i) {
      //
      dtemp = fit[i] - ener[i];
	
      res += dtemp * dtemp;
    }
    
    res /= (double)ener.size();
    
    ener_rmsd[d] = std::sqrt(res);
  }
  
  // energy RMSD output
  //
  std::cout << "\nEnergy RMSD [kcal/mol]:\n" << std::setw(7) << "R, bohr";

  for(int d = 0; d < dist_size; ++d)
    //
    std::cout << std::setw(7) << dist_data[d];

  std::cout << "\n" << std::setw(7) << "ERMSD";
    
  for(int d = 0; d < dist_size; ++d) {
    //
    std::cout << std::setw(7);

    dtemp = ener_rmsd[d];
    
    if(dtemp > 0.1) {
      //
      std::cout << dtemp;
    }
    else
      //
      std::cout << "<0.1";
  }
  std::cout << "\n";
  
  if(check_ray.size()) {
    //
    std::cout << "\nlarge deviation ray samplings:";

    for(std::map<double, std::pair<int, int> >::const_reverse_iterator rit = check_ray.rbegin(); rit != check_ray.rend(); ++rit)
      //
      std::cout << "  (" << rit->second.first << ", " << rit->second.second << ")";

    std::cout << "\n";
  }
  
  return 0;
}
