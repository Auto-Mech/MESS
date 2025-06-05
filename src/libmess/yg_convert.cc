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

class InputError {};

int main (int argc, char* argv[])
{
  const char funame [] = "yg_convert: ";

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
  std::string smp_file, out_file;

  KeyGroup YGConvertGroup;

  Key   struc_key("Structure"      );
  Key     smp_key("SmpFile"        );
  Key     out_key("OutFile"        );
  
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

  std::ofstream to(out_file.c_str());
  
  // generalized coordinates
  //
  Dynamic::Coordinates dc;

  for(int i = 0; i < 4; ++i)
    //
    if(!i) {
      //
      dc.write_ang_pos(0)[i] = 1.;
    }
    else
      //
      dc.write_ang_pos(0)[i] = 0.;

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
      
    dist_data[d] = dtemp;
  }

  std::vector<Atom> mf = Structure::fragment(1);

  std::vector<D3::Vector> rel_pos;
  
  while(1) {
    //
    int ray_index;
      
    if(!(from >> ray_index)) {
      //
      std::cerr << "reading " << smp_file << " file done\n";

      break;
    }

    lin.read_line(from);

    for(int i = 0; i < ang_size; ++i)
      //
      lin >> ang_pos[i];

    if(!lin) {
      //
      std::cerr << ray_index << "-th ray sampling: reading angular position failed\n";

      throw InputError();
    }

    // orient second fragment
    //
    for(int i = 0; i < Structure::pos_size(1); ++i)
      //
      dc.write_ang_pos(1)[i] = ang_pos[i + 3];
    
    dc.set_rel_pos(1, rel_pos);

    // energies & output
    //
    for(int d = 0; d < dist_size; ++d) {
      //
      from >> dtemp;

      to << Structure::size() << "\n" << 1 << "  " << 2.2 << "  " << dtemp << "\n";

      for(int a = 0; a < Structure::size(0); ++a)
	//
	to << Structure::fragment(0)[a] << "\n";
      
      for(int a = 0; a < Structure::size(1); ++a) {
	//
	mf[a] = rel_pos[a];

	for(int i = 0; i < 3; ++i)
	  //
	  mf[a][i] += dist_data[d] * ang_pos[i];

	to << mf[a] << "\n";
      }
      to << "\n";
    }
    
    if(!from) {
      //
      std::cerr << ray_index << "-th ray sampling: reading energies failed\n";

      throw InputError();
    }
  }

  return 0;
}
