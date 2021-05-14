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

#include "atom.hh"
#include "units.hh"
#include "io.hh"

#include <sstream>
#include <cstdlib>
#include <vector>

/************************** Atom description ****************************/

AtomBase::DataBase AtomBase::_data;

AtomBase::DataBase::DataBase ()
{
  const char funame [] = "AtomBase::DataBase::DataBase: ";

  (*this)[DUMMY].name    = "X";
  (*this)[DUMMY].valence = 0;
  (*this)[DUMMY].isotope = 0;

  (*this)[HYDROGEN].name    = "H";
  (*this)[HYDROGEN].valence = 1;
  (*this)[HYDROGEN].isotope = 1;
  (*this)[HYDROGEN].mass[1] = 1.007825;
  (*this)[HYDROGEN].mass[2] = 2.014;
  (*this)[HYDROGEN].mass[3] = 3.01605;

  (*this)[HELIUM].name    = "He";
  (*this)[HELIUM].valence = 0;
  (*this)[HELIUM].isotope = 4;
  (*this)[HELIUM].mass[4] = 4.;

  (*this)[CARBON].name     = "C";
  (*this)[CARBON].valence  = 4;
  (*this)[CARBON].isotope  = 12;
  (*this)[CARBON].mass[12] = 12.;
  (*this)[CARBON].mass[13] = 13.00335;

  (*this)[NITROGEN].name     = "N";
  (*this)[NITROGEN].valence  = 3;
  (*this)[NITROGEN].isotope  = 14;
  (*this)[NITROGEN].mass[14] = 14.00307;
  (*this)[NITROGEN].mass[15] = 15.00011;

  (*this)[OXYGEN].name     = "O";
  (*this)[OXYGEN].valence  = 2;
  (*this)[OXYGEN].isotope  = 16;
  (*this)[OXYGEN].mass[16] = 15.99491;
  (*this)[OXYGEN].mass[17] = 17.;
  (*this)[OXYGEN].mass[18] = 18.;

  (*this)[FLUORINE].name     = "F";
  (*this)[FLUORINE].valence  = 1;
  (*this)[FLUORINE].isotope  = 19;
  (*this)[FLUORINE].mass[19] = 18.9984;

  (*this)[SODIUM].name     = "Na";
  (*this)[SODIUM].valence  = 1;
  (*this)[SODIUM].isotope  = 23;
  (*this)[SODIUM].mass[23] = 22.9898;

  (*this)[SILICON].name     = "Si";
  (*this)[SILICON].valence  = 4;
  (*this)[SILICON].isotope  = 28;
  (*this)[SILICON].mass[28] = 27.97693;
  (*this)[SILICON].mass[29] = 28.97649;
  (*this)[SILICON].mass[30] = 29.97376;

  (*this)[PHOSPHORUS].name     = "P";
  (*this)[PHOSPHORUS].valence  = 3;
  (*this)[PHOSPHORUS].isotope  = 31;
  (*this)[PHOSPHORUS].mass[31] = 30.97376;

  (*this)[SULFUR].name     = "S";
  (*this)[SULFUR].valence  = 2;
  (*this)[SULFUR].isotope  = 32;
  (*this)[SULFUR].mass[32] = 31.9720;
  (*this)[SULFUR].mass[33] = 32.97146;
  (*this)[SULFUR].mass[34] = 33.96786;
  (*this)[SULFUR].mass[36] = 35.96709;

  (*this)[CHLORINE].name     = "Cl";
  (*this)[CHLORINE].valence  = 1;
  (*this)[CHLORINE].isotope  = 35;
  (*this)[CHLORINE].mass[35] = 34.96885;
  (*this)[CHLORINE].mass[37] = 37.;

  (*this)[TITANIUM].name     = "Ti";
  (*this)[TITANIUM].valence  = 1;
  (*this)[TITANIUM].isotope  = 49;
  (*this)[TITANIUM].mass[46] = 45.95263;
  (*this)[TITANIUM].mass[47] = 46.95176;
  (*this)[TITANIUM].mass[48] = 47.947947;
  (*this)[TITANIUM].mass[49] = 48.947871;
  (*this)[TITANIUM].mass[50] = 49.944792;

  (*this)[BROMINE].name     = "Br";
  (*this)[BROMINE].valence  = 1;
  (*this)[BROMINE].isotope  = 79;
  (*this)[BROMINE].mass[79] = 78.9183;
  (*this)[BROMINE].mass[81] = 80.9163;

  (*this)[URANIUM].name     = "U";
  (*this)[URANIUM].valence  = 6;
  (*this)[URANIUM].isotope  = 238;
  (*this)[URANIUM].mass[238] = 238.029;

  // name-to-number map
  for(std::map<int, Data>::const_iterator it = this->begin(); it != this->end(); ++it)
    _name_num_map[it->second.name] = it->first;
}

int AtomBase::DataBase::number (const std::string& n) const 
{
  const char funame [] = "AtomBase::DataBase::number: ";
  
  std::map<std::string, int>::const_iterator nit = _name_num_map.find(n);
  if(nit != _name_num_map.end())
    return nit->second;

  std::cerr << funame << n << ": not implemented\n";
  throw Error::Init();
}

const std::string&  AtomBase::DataBase::name (int n)  const 
{
  const char funame [] = "AtomBase::DataBase::name: ";

  std::map<int, Data>::const_iterator nit = find(n);
  if(nit != end())
    return nit->second.name;

  std::cerr << funame << "not implemented\n";
  throw Error::Init();
}

int AtomBase::DataBase::default_isotope (int n)  const 
{
  const char funame [] = "AtomBase::DataBase::default_isotope: ";

  std::map<int, Data>::const_iterator nit = find(n);
  if(nit != end())
    return nit->second.isotope;

  std::cerr << funame << "not implemented\n";
  throw Error::Init();
}

int AtomBase::DataBase::valence (int n)  const 
{
  const char funame [] = "AtomBase::DataBase::valence: ";

  std::map<int, Data>::const_iterator nit = find(n);
  if(nit != end())
    return nit->second.valence;

  std::cerr << funame << "not implemented\n";
  throw Error::Init();
}

double AtomBase::DataBase::mass (int n, int i) const 
{
  const char funame [] = "AtomBase::DataBase::mass: ";

  std::map<int, Data>::const_iterator nit = find(n);
  if(nit != end()) {
    // dummy atom mass can be anything
    if(nit->first == DUMMY)
      return (double)i * Phys_const::amu;

    std::map<int, double>::const_iterator it = nit->second.mass.find(i);
    if(it != nit->second.mass.end())
      return it->second * Phys_const::amu;
  }

  std::cerr << funame << "not implemented\n";
  throw Error::Init();
}

void Atom::_read (std::istream& from) 
{
    const char funame [] = "Atom::_read: ";

    int         itemp;
    double      dtemp;
    std::string stemp;
    bool        btemp;

    IO::skip_space(from);

    IO::LineInput lin(from);

    // find the number of items on line
    stemp = lin.str();
    //std::cout << stemp << "\n";
    std::istringstream iss(stemp);
    int item_num = 0;
    while(iss >> stemp) {
      // check for comments
      bool is_comment = false;
      for(std::string::const_iterator it = IO::comment_symbol().begin(); it != IO::comment_symbol().end(); ++it)
	if(*it == stemp[0]) {
	  is_comment = true;
	  break;
	}
      if(is_comment)
	break;

      ++item_num;
    }

    // gaussian input
    if(item_num == 6) {
      // atomic index
      lin >> itemp;
 
      // atomic number
      lin >> itemp;
      set(itemp);

      // atomic type
      lin >> itemp;

      for(int i = 0; i < 3; ++i)
	lin >> (*this)[i];

      if(!lin) {
	std::cerr << funame << "assumed gaussian input is corrupted\n";
	throw Error::Input();
      }
      
      return;
    }

    if(!(lin >> stemp)) {
      std::cerr << funame << "input stream is corrupted\n";
      throw Error::Input();
    }
    set(stemp);
    
    switch(item_num) {
    case 1: // atom name only
      return;

    case 2: // atom name and isotope number
      if(!(lin >> itemp)) {
	std::cerr << funame << "corrupted\n";
	throw Error::Input();
      }

      if(itemp <= 0) {
	std::cerr << funame << "out of range\n";
	throw Error::Range();
      }

      set_isotope(itemp);
      return;

    case 4: // atom name and coordinates
      for(int i = 0; i < 3; ++i) {
	if(!(lin >> dtemp)) {
	  std::cerr << funame << "corrupted\n";
	  throw Error::Input();
	}
	(*this)[i] = dtemp;
      }
      return;

    case 5: // atom name, isotope number, and coordinates
      if(!(lin >> itemp)) {
	std::cerr << funame << "corrupted\n";
	throw Error::Input();
      }

      if(itemp <= 0) {
	std::cerr << funame << "out of range\n";
	throw Error::Range();
      }

      set_isotope(itemp);

      for(int i = 0; i < 3; ++i) {
	if(!(lin >> dtemp)) {
	  std::cerr << funame << "corrupted\n";
	  throw Error::Input();
	}
	(*this)[i] = dtemp;
      }
      return;

    default:
      std::cerr << funame << "wrong number of items, " << item_num << "\n";
      throw Error::Input();
    }
}

std::ostream& operator<< (std::ostream& out , const Atom& a)
{
  out << std::setw(2) << a.name() 
      << std::fixed; // fixed format
  for(int i = 0; i < 3; ++i)
    out << std::setw(11) << a[i];

  out.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  return out;
}

bool are_equal (double r1, double r2, double tol) 
{
  double dtemp = r1 > r2 ? r1 - r2 : r2 - r1;

  if(dtemp < tol)
    return true;

  return false;
}

bool are_equal (const Atom& a, const Atom& b, double tol, int flags)
{
  if(flags & IGNORE_ISOTOPE) {
    if(a.number() != b.number())
      return false;
  }
  else if((const AtomBase&)a != (const AtomBase&)b)
    return false;

  if(vdistance(a, b) > tol)
    return false;

  return true;
}

bool are_equal (const std::vector<Atom>& m, const std::vector<Atom>& n, double tol, int flags)
{
  if(m.size() != n.size())
    return false;

  std::vector<Atom>::const_iterator  nit = n.begin();
  for(std::vector<Atom>::const_iterator mit = m.begin(); mit != m.end(); ++mit, ++nit)
    if(!are_equal(*mit, *nit, tol, flags))
      return false;

  return true;
}

Permutation is_symmetric (const std::vector<Atom>& mol, const Symmetry::SpaceElement& symmel, double tolerance, int flags)
{
  const char funame [] = "is_symmetric: ";
  
  if(mol.size() < 2) {
    std::cerr << funame << "not a molecule\n";
    throw Error::Range();
  }

  if(tolerance <= 0. || tolerance >= 1.) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  bool btemp;

  std::vector<Atom> test = mol;
  symmel.apply(mol, test);
  
  std::vector<int> perm(mol.size());
  for(std::vector<Atom>::const_iterator mat = mol.begin(); mat != mol.end(); ++mat) {
    btemp = false;
    for(std::vector<Atom>::const_iterator tat = test.begin(); tat != test.end(); ++tat) 
      if(are_equal(*mat, *tat, tolerance, flags)) {
	if(btemp) {
	  std::cerr << funame << "identical atoms\n";
	  throw Error::Logic();
	}
	btemp = true;
	perm[mat - mol.begin()]  = tat - test.begin();	
      }
    if(!btemp) {
      return Permutation();
    }
  }

  return Permutation(perm);
}

std::set<Permutation> identical_atoms_permutation_symmetry_group (const Lapack::SymmetricMatrix& dist, double tol) 
{
  const char funame [] = "identical_atoms_permutation_symmetry_group: ";

  bool btemp;

  std::set<Permutation> res;
  if(!dist.size())
    return res;

  if(dist.size() == 1) {
    res.insert(Permutation(1));
    return res;
  }

  std::vector<std::vector<int> > perm(1);

  for(int n = 0; n < dist.size(); ++n) {
    std::vector<std::vector<int> > new_perm;
    for(int p = 0; p < perm.size(); ++p) {
      std::set<int> used;
      for(int i = 0; i < n; ++i)
	used.insert(perm[p][i]);

      for(int ii = 0; ii < dist.size(); ++ii) {
	if(used.find(ii) != used.end())
	  continue;

	btemp = false;
	for(int i = 0; i < n; ++i)
	  if(!are_equal(dist(i, n), dist(perm[p][i], ii), tol)) {
	    btemp = true;
	    break;
	  }
	if(btemp)
	  continue;

	new_perm.push_back(perm[p]);
	new_perm.back().push_back(ii);
      }
    }
    perm = new_perm;
  }

  for(int i = 0; i < perm.size(); ++i)
    res.insert(Permutation(perm[i], Permutation::NOCHECK));

  // check that the result is the group
  for(std::set<Permutation>::const_iterator i = res.begin(); i != res.end(); ++i)
    for(std::set<Permutation>::const_iterator j = res.begin(); j != res.end(); ++j)
      if(res.find(*i * *j) == res.end()) {
	std::cerr << funame << "not a group\n";
	throw Error::Logic();
      }

  return res;
}

std::set<Permutation> permutation_symmetry_group (const std::vector<Atom>& molecule, double tolerance, int flags) 
{
  const char funame [] = "permutation_symmetry_group: ";
  
  int    itemp;
  double dtemp;
  bool   btemp;

  if(tolerance <= 0. || tolerance >= 1.) {
    std::cerr << funame << "tolerance out of range\n";
    throw Error::Range();
  }

  // no atoms
  if(molecule.size() < 2) {
    std::cerr << funame << "no molecule\n";
    throw Error::Range();
  }

  // no dummies
  for(int i = 0; i < molecule.size(); ++i)
    if(molecule[i].number() == AtomBase::DUMMY) {
      std::cerr << funame << "not implemented for dummies\n";
      throw Error::Run();
    }

  // distance matrix
  Lapack::SymmetricMatrix dist(molecule.size());
  for(int i = 0; i < molecule.size(); ++i) 
    for(int j = 0; j < i; ++j)
      dist(i, j) = vdistance(molecule[i], molecule[j]);

  // identical atoms groups
  std::vector<std::vector<int> > atom_group;
  if(flags & IGNORE_ISOTOPE) {
    std::map<int, std::vector<int> > ag;
    for(int i = 0; i < molecule.size(); ++i)
      ag[molecule[i].number()].push_back(i);

    for(std::map<int, std::vector<int> >::const_iterator it = ag.begin(); it != ag.end(); ++it)
      atom_group.push_back(it->second);
  }
  else {
    std::map<AtomBase, std::vector<int> > ag;
    for(int i = 0; i < molecule.size(); ++i)
      ag[molecule[i]].push_back(i);
    
    for(std::map<AtomBase, std::vector<int> >::const_iterator it = ag.begin(); it != ag.end(); ++it)
      atom_group.push_back(it->second);
  }

  // identical atoms permutation symmetry groups
  std::vector<std::set<Permutation> > atom_group_perm(atom_group.size());
  for(int g = 0; g < atom_group.size(); ++g) {
    Lapack::SymmetricMatrix group_dist(atom_group[g].size());
    for(int i = 0; i < atom_group[g].size(); ++i)
      for(int j = 0; j < i; ++j)
	group_dist(i, j) = dist(atom_group[g][i], atom_group[g][j]);
    
    atom_group_perm[g] = identical_atoms_permutation_symmetry_group(group_dist, tolerance);
  }

  // all atoms permutation symmetry group
  std::vector<std::vector<Permutation> > perm(1);
  for(int g = 0; g < atom_group.size(); ++g) {
    std::vector<std::vector<Permutation> > new_perm;
    for(int p = 0; p < perm.size(); ++p)
      for(std::set<Permutation>::const_iterator q = atom_group_perm[g].begin(); q != atom_group_perm[g].end(); ++q) {
	btemp = false;
	for(int i = 0; i < atom_group[g].size(); ++i)
	  for(int f = 0; f < g; ++f)
	    for(int j = 0; j < atom_group[f].size(); ++j)
	      if(!are_equal(dist(atom_group[f][j], atom_group[g][i]), dist(atom_group[f][perm[p][f][j]],atom_group[g][(*q)[i]]), tolerance)) {
		btemp = true;
		break;
	      }
	if(btemp)
	  continue;
	new_perm.push_back(perm[p]);
	new_perm.back().push_back(*q);
      }
    perm = new_perm;
  }

  std::set<Permutation> res;
  for(int p = 0; p < perm.size(); ++p) {
    std::vector<int> gp(dist.size());
    for(int g = 0; g < atom_group.size(); ++g)
      for(int a = 0; a < atom_group[g].size(); ++a)
	gp[atom_group[g][a]] = atom_group[g][perm[p][g][a]];
    
    res.insert(Permutation(gp, Permutation::NOCHECK));
  }

  // check that the result is the group
  for(std::set<Permutation>::const_iterator i = res.begin(); i != res.end(); ++i)
    for(std::set<Permutation>::const_iterator j = res.begin(); j != res.end(); ++j)
      if(res.find(*i * *j) == res.end()) {
	std::cerr << funame << "not a group\n";
	throw Error::Logic();
      }

  return res;
}

std::pair<int, int> symmetry_number (const std::vector<Atom>& molecule, double tolerance, int flags) 
{
  const char funame [] = "symmetry_number: ";
  
  static const double min_dist = 0.5;

  int    itemp;
  double dtemp;
  bool   btemp;

  if(tolerance <= 0. || tolerance >= 1.) {
    std::cerr << funame << "tolerance out of range\n";
    throw Error::Range();
  }

  // no atoms
  if(molecule.size() < 2) {
    std::cerr << funame << "no molecule\n";
    throw Error::Range();
  }

  // no dummies
  for(int i = 0; i < molecule.size(); ++i)
    if(molecule[i].number() == AtomBase::DUMMY) {
      std::cerr << funame << "not implemented for dummies\n";
      throw Error::Run();
    }

  // minimal distance check
  for(int i = 0; i < molecule.size(); ++i)
    for(int j = 0; j < i; ++j)
      if((molecule[i] - molecule[j]).vlength() < min_dist) {
	std::cerr << funame << i << "," << j << "-th atom distance is too small\n";
	throw Error::Range();
      }
  
  bool is_plane = false;
  int max_volume_set[3];
  double max_volume = 0.;

  if(molecule.size() < 4)
    is_plane = true;
  else {
    for(int i = 1; i < molecule.size(); ++i)
      for(int j = 1; j < i; ++j)
	for(int k = 1; k < j; ++k) { 
	  dtemp = D3::volume(molecule[i] - molecule[0], molecule[j] - molecule[0], molecule[k] - molecule[0])
	    / (molecule[i] - molecule[0]).vlength() 
	    / (molecule[j] - molecule[0]).vlength() 
	    / (molecule[k] - molecule[0]).vlength();

	  if(dtemp > max_volume) {
	    max_volume_set[0] = i;
	    max_volume_set[1] = j;
	    max_volume_set[2] = k;
	    max_volume = dtemp;
	  }
	  else if(dtemp < -max_volume) {
	    max_volume_set[0] = j;
	    max_volume_set[1] = i;
	    max_volume_set[2] = k;
	    max_volume = -dtemp;
	  }
	}

    if(max_volume < tolerance)
      is_plane = true;
  }

  std::set<Permutation> symm_group = permutation_symmetry_group(molecule, tolerance, flags);

  if(is_plane)
    return std::make_pair<int, int>(symm_group.size(), 1);

  // rotational symmetry number
  int rot = 0;
  for(std::set<Permutation>::const_iterator it = symm_group.begin(); it != symm_group.end(); ++it)
    if(D3::volume(molecule[(*it)[max_volume_set[0]]] - molecule[(*it)[0]], 
		  molecule[(*it)[max_volume_set[1]]] - molecule[(*it)[0]],
		  molecule[(*it)[max_volume_set[2]]] - molecule[(*it)[0]]) > 0.)
      ++rot;

  if(rot != symm_group.size())
    itemp = 1;
  else
    itemp = 2;

  return std::make_pair(rot, itemp);

}

Symmetry::SpaceGroup spatial_symmetry_group (std::vector<Atom>& molecule, double tolerance, int flags) 
{
  const char funame [] = "spatial_symmetry_group: ";
  
  static const double min_dist = 0.5;

  int         itemp;
  double      dtemp;
  bool        btemp;
  std::string stemp;
  D3::Vector  vtemp;

  if(tolerance <= 0. || tolerance >= 1.) {
    std::cerr << funame << "tolerance out of range\n";
    throw Error::Range();
  }

  // no atoms
  if(molecule.size() < 2) {
    std::cerr << funame << "no molecule\n";
    throw Error::Range();
  }

  // no dummies
  for(int i = 0; i < molecule.size(); ++i)
    if(molecule[i].number() == AtomBase::DUMMY) {
      std::cerr << funame << "not implemented for dummies\n";
      throw Error::Run();
    }

  // minimal distance check
  for(int i = 0; i < molecule.size(); ++i)
    for(int j = 0; j < i; ++j)
      if((molecule[i] - molecule[j]).vlength() < min_dist) {
	std::cerr << funame << i << "," << j << "-th atom distance is too small\n";
	throw Error::Range();
      }
  
  // center of mass/charge shift
  vtemp = 0.;
  dtemp = 0.;
  if(flags & IGNORE_ISOTOPE)
    for(int a = 0; a < molecule.size(); ++a) {
      vtemp += molecule[a] * double(molecule[a].number());
      dtemp += double(molecule[a].number());
    }
  else
    for(int a = 0; a < molecule.size(); ++a) {
      vtemp += molecule[a] * molecule[a].mass();
      dtemp += molecule[a].mass();
    }
  vtemp /= dtemp;
  for(int a = 0; a < molecule.size(); ++a)
    molecule[a] -= vtemp;

  // symmetry group
  std::set<Permutation> symm_group = permutation_symmetry_group(molecule, tolerance, flags);
  std::vector<Symmetry::SpaceElement> group_base(1, Symmetry::unity());

  // no symmetry
  if(symm_group.size() == 1)
    return Symmetry::SpaceGroup("C1", group_base, 1);

  D3::Vector nx, ny, nz;

  // is linear?
  bool is_linear = false;
  int max_sector_set [2];
  double max_sector = 0.;

  if(molecule.size() < 3)
    is_linear = true;
  else {
    for(int i = 1; i < molecule.size(); ++i)
      for(int j = 1; j < i; ++j) { 
	vtemp = D3::vprod(molecule[i] - molecule[0], molecule[j] - molecule[0]);
	dtemp = vtemp.vlength() 
	  / (molecule[i] - molecule[0]).vlength() 
	  / (molecule[j] - molecule[0]).vlength();

	if(dtemp > max_sector) {
	  max_sector_set[0] = i;
	  max_sector_set[1] = j;
	  nz = vtemp;
	  max_sector = dtemp;
	}
      }

    if(max_sector < tolerance)
      is_linear = true;
  }
  
  // linear molecule
  if(is_linear) {
    nx = molecule[1] - molecule[0];
    nx.normalize();
    for(int a = 0; a < molecule.size(); ++a) {
      dtemp = vdot(nx, (const D3::Vector&)molecule[a]);
      (D3::Vector&)molecule[a] = 0.;
      molecule[a][0] = dtemp;
    }

    if(symm_group.size() != 2) {
      std::cerr << funame << "linear molecule: symmetry group wrong size: " << symm_group.size() << "\n";
      throw Error::Logic();
    }

    btemp = false;
    for(int a = 0; a < molecule.size(); ++a) {
      int b = (*symm_group.rbegin())[a]; 
      if(b == a) {
	if(btemp) {
	  std::cerr << funame << "linear molecule: more than one central atom\n";
	  throw Error::Logic();
	}
	if(molecule[a][0] > tolerance || molecule[a][0] < -tolerance) {
	  std::cerr << funame << "linear molecule: " << a << "-th central atom: not symmetric \n";
	  throw Error::Logic();
	}
	molecule[a][0] = 0.;
      }
      else {
	dtemp = molecule[b][0] + molecule[a][0];
	if(dtemp > tolerance || dtemp < -tolerance) {
	  std::cerr << funame << "linear molecule: "<< a << "-th and " << b << "-th atoms: not symmetric\n";
	  throw Error::Logic();
	}
	molecule[b][0] = -molecule[a][0];
      }
    }
    
    group_base.push_back(Symmetry::inversion());
    return Symmetry::SpaceGroup("Ci", group_base, 2);
  }// linear molecule

  nz.normalize();

  // is plane?
  bool is_plane = false;
  int max_volume_set [3];
  double max_volume = 0.;

  if(molecule.size() < 4)
    is_plane = true;
  else {
    for(int i = 1; i < molecule.size(); ++i)
      for(int j = 1; j < i; ++j)
	for(int k = 1; k < j; ++k) { 
	  dtemp = D3::volume(molecule[i] - molecule[0], molecule[j] - molecule[0], molecule[k] - molecule[0])
	    / (molecule[i] - molecule[0]).vlength() 
	    / (molecule[j] - molecule[0]).vlength() 
	    / (molecule[k] - molecule[0]).vlength();
	  if(dtemp > max_volume) {
	    max_volume_set[0] = i;
	    max_volume_set[1] = j;
	    max_volume_set[2] = k;
	    max_volume = dtemp;
	  }
	  else if(dtemp < -max_volume) {
	    max_volume_set[0] = j;
	    max_volume_set[1] = i;
	    max_volume_set[2] = k;
	    max_volume = -dtemp;
	  }
	}
    if(max_volume < tolerance)
      is_plane = true;
  }

  std::set<Permutation> rot_group;
  std::ostringstream group_name;

  // plane molecule
  if(is_plane) {
    group_base.push_back(Symmetry::reflection(2));

    // main axis rotation subgroup
    for(std::set<Permutation>::const_iterator it = symm_group.begin(); it != symm_group.end(); ++it)
      if(vdot(D3::vprod(molecule[(*it)[max_sector_set[0]]] - molecule[(*it)[0]], 
			molecule[(*it)[max_sector_set[1]]] - molecule[(*it)[0]]), nz) > 0.)
	rot_group.insert(*it);

    if(rot_group.size() > 1)
      group_base.push_back(Symmetry::rotation(2, (int)rot_group.size()));

    // Cnh symmetry
    if(rot_group.size() == symm_group.size()) {
      nx = molecule[1] - molecule[0];
      ny = D3::vprod(nz, nx);

      D3::Matrix trans(nx, ny);

      for(int a = 0; a < molecule.size(); ++a) {
	molecule[a] = trans * molecule[a];
	molecule[a][2] = 0.;
      }

      if(rot_group.size() > 1)
	group_name << "C" << rot_group.size() << "h";
      else
	group_name << "Ch";

      return Symmetry::SpaceGroup(group_name.str(), group_base, 2 * rot_group.size());
    }// Cnh symmetry
    
    // Dnh symmetry
    if(symm_group.size() != 2 * rot_group.size()) {
      std::cerr << funame << "plane molecule: unknown symmetry group size: " << symm_group.size() << "\n";
      throw Error::General();
    }

    // in-plane rotation
    std::set<Permutation>::const_iterator it;
    for(it = symm_group.begin(); it != symm_group.end(); ++it)
      if(rot_group.find(*it) == rot_group.end())
	break;

    D3::Vector v1 = molecule[1] - molecule[0];
    D3::Vector v2 = molecule[(*it)[1]] - molecule[(*it)[0]];

    if((v1 + v2).vlength() > (v1 - v2).vlength()) {
      nx = v1 + v2;
      ny = D3::vprod(nz, nx);
    }
    else {
      ny = v1 - v2;
      nx = D3::vprod(ny, nz);
    }
     
    D3::Matrix trans(nx, ny);
    for(int a = 0; a < molecule.size(); ++a) {
      molecule[a] = trans * molecule[a];
      molecule[a][2] = 0.;
    }

    group_base.push_back(Symmetry::rotation(0));

    if(rot_group.size() > 1)
      group_name << "D" << rot_group.size() << "h";
    else 
      group_name << "C2v";
    
    return Symmetry::SpaceGroup(group_name.str(), group_base, 4 * rot_group.size());
  }// plane molecule

  Lapack::Matrix mtemp(3);
  for(int i = 0; i < 3; ++i)
    mtemp.row(i) = molecule[max_volume_set[i]] - molecule[0];
  mtemp = mtemp.invert();

  D3::Matrix max_vol_frame;
  for(int i = 0; i < 3; ++i)
    max_vol_frame.row(i) = mtemp.row(i);

  std::map<Permutation, Symmetry::SpaceElement> perm_space_map;

  // rotation subgroup
  for(std::set<Permutation>::const_iterator it = symm_group.begin(); it != symm_group.end(); ++it) {
    D3::Matrix rmatrix;
    for(int i = 0; i < 3; ++i)
      rmatrix.row(i) = molecule[(*it)[max_volume_set[i]]] - molecule[(*it)[0]];

    rmatrix = max_vol_frame * rmatrix;
    rmatrix.orthogonality_check(tolerance);
    rmatrix.orthogonalize();

    if(!rmatrix.inversion())
      rot_group.insert(*it);

    perm_space_map[*it] = Symmetry::SpaceElement(rmatrix);
  }


  if(*rot_group.begin() != Permutation(molecule.size())) {
    std::cerr << funame << "first group element should be identical permutation\n";
    throw Error::Logic();
  }

  itemp = symm_group.size() / rot_group.size();
  if(itemp != 1 && itemp != 2 || symm_group.size() % rot_group.size()) {
    std::cerr << funame << "rotation group order should divide symmetry group order\n";
    throw Error::Logic();
  }

  if(rot_group.size() == 1) {
    

  }

  // cyclic rotational subgroups
  std::map<Permutation, std::set<Permutation> > cg_map;
  for(std::set<Permutation>::const_iterator rit = rot_group.begin(); rit != rot_group.end(); ++rit) {
    if(rit == rot_group.begin())
      continue;

    btemp = false;
    for(std::map<Permutation, std::set<Permutation> >::const_iterator mit = cg_map.begin(); mit != cg_map.end(); ++mit)
      if(mit->second.find(*rit) != mit->second.end()) {
	btemp = true;
	break;
      }

    if(btemp)
      continue;

    std::set<Permutation> new_cg;
    for(Permutation perm = *rit; perm != *rot_group.begin(); perm = perm * *rit)
      new_cg.insert(perm);

    do {
      btemp = false;
      for(std::map<Permutation, std::set<Permutation> >::iterator mit = cg_map.begin(); mit != cg_map.end(); ++mit)
	if(new_cg.find(mit->first) != new_cg.end()) {
	  cg_map.erase(mit);
	  btemp = true;
	  break;
	}
    } while(btemp);

    cg_map[*rit] = new_cg;
  }
  
  //cyclic rotational subgroups dimensions
  std::map<int, int> cg_size;
  for(std::map<Permutation, std::set<Permutation> >::const_iterator mit = cg_map.begin(); mit != cg_map.end(); ++mit)
    ++cg_size[mit->second.size()];

  return Symmetry::SpaceGroup("C1", group_base);
}

