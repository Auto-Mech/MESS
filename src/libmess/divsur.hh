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

#ifndef DIVSUR_HH
#define DIVSUR_HH

#include "dynamic.hh"
#include "io.hh"

namespace DivSur {

  // Abstract class for geometric objects
  class GeomObject
  {
  public:
    virtual ~GeomObject () {}
  };

  class GeomPlane : public GeomObject, public D3::Plane
  {
  public:
    GeomPlane  (const D3::Vector& n, double d) : D3::Plane(n, d) {}
    ~GeomPlane () {}
  };

  class GeomPivot : public GeomObject, public D3::Vector
  {
  public:
    GeomPivot (const D3::Vector& v) : D3::Vector(v) {}
    ~GeomPivot () {}
  };

  
  // auxiliary class to read the geometrical objects 
  // associated with a fragment
  struct FragData
  {
    std::string name;
    std::map<std::string, ConstSharedPointer<GeomObject> >  data;

    void read (int frag, std::istream&) ;
  };


  /*********************************************************************
   * Abstract elemental surface: the concept of the analytical surface *
   * in the configurational space. It should have three features:      *
   *   1.) to allow randomly orient fragments on its surface           *
   *   2.) to provide a statistical weight of such sampling            *
   *   3.) separate the configurational space into two regions         *
   *********************************************************************/

  class Primitive {
  protected:
    // Normalized generalized reaction coordinate vector and its length
    // Should be redefined for each real surface
    virtual double _set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_vec) const; 

    // randomly generate velocities from the normalized reaction coordinate vector
    static void   _rc2tv (const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, 
			  double temper, Dynamic::Momenta& dm);// thermal
    static double _rc2ev (const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, 
			  double ener, Dynamic::Momenta& dm);// E-resolved
    static double _rc2jv (const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, 
			  double ener, double amom, Dynamic::Momenta& dm); // E,J-resolved
  public:

    // Randomly orient the fragments on the surface
    // Should be redefined for each surface
    virtual void   random_orient   (Dynamic::Coordinates& dv) const =0;

    // Signed distance of the configuration from the surface
    // Should be re-defined for each real surface
    virtual double distance  (const Dynamic::Coordinates& dv) const;

    // output the surface
    // should be redefined for each real surface
    virtual void print (std::ostream&) const = 0;

    // Test if the configuration is inside of the surface
    // Inside corresponds to the positive distance
    bool test (const Dynamic::Coordinates& dv) const;

    // The statistical weight of the configuration on the surface
    double weight (const Dynamic::Coordinates& dv) const;

    // randomly generate velocities for  the given configuration 
    // on the dividing surface. The velocities are directed inside the surface.
    // The output is the statistical weight of the generated configuration
    double t_velocity (double temper,            Dynamic::Vars& dv) const;// thermal
    double e_velocity (double ener,              Dynamic::Vars& dv) const;// E-resolved
    double j_velocity (double ener, double amom, Dynamic::Vars& dv) const;// E,J-resolved

    virtual ~Primitive () {}

  };

  inline bool Primitive::test  (const Dynamic::Coordinates& dv) const 
  {
    if (distance(dv) > 0.) 
      return true; 
    else 
      return false; 
  }

  /************************************************************************
   *        Spherical surface: the distance between two pivot points        *
   *         associated with the corresponding fragments is fixed         *
   ************************************************************************/

  class Sphere : public Primitive {
    D3::Vector _pivot[2];
    double _dist;

    // lab frame pp2pp vector, 1->2
    D3::Vector _lf_pp_12 (const Dynamic::Coordinates& dv) const; 
    double _set_rc (const Dynamic::Coordinates& dv, Dynamic::Momenta& rc_vec) const; 

  public:
    Sphere (ConstSharedPointer<GeomObject>*, double) ;

    double distance  (const Dynamic::Coordinates& dv) const
    { return _dist - _lf_pp_12(dv).vlength(); } 
    void   random_orient   (Dynamic::Coordinates& dv) const;

    ~Sphere () {}

    void print (std::ostream&) const;

    // specific functions
    double max_cm_dist () const 
    { return _dist + _pivot[0].vlength() + _pivot[1].vlength(); }

  };

  inline Sphere::Sphere (ConstSharedPointer<GeomObject>* gop, double d) 
  {
    const char funame [] = "DivSur::Sphere::Sphere (ConstSharedPointer<GeomObject>*, double): "; 

    _dist = d;
    const GeomPivot* pvt;
    for(int frag = 0; frag < 2; ++frag) {
      pvt = dynamic_cast<const GeomPivot*>((const GeomObject*)(gop[frag]));
      if(!pvt) {
	std::cerr << funame << frag << "-th fragment: object is not a pivot point\n";
	throw Error::Form();
      }
      _pivot[frag] = *pvt;
    }
  }


  /********************************************************************************
   *        Plane surface: pivot point of the second fragment belongs to the        *
   *            geometrical plane associated with the first fragment              *
   ********************************************************************************/

  class Plane : public Primitive {
    D3::Plane  _plane;
    D3::Vector _pivot;
    int _plane_frag; // on which fragment plane is

    double _radius; // sampling circle radius

    double _set_rc (const Dynamic::Coordinates& dv, Dynamic::Momenta& rc_vec) const; 

  public:
    Plane (ConstSharedPointer<GeomObject>*) ;
  
    double distance  (const Dynamic::Coordinates& dv) const;
    void   random_orient   (Dynamic::Coordinates& dv) const;

    ~Plane () {}

    void print (std::ostream&) const;

    // specific functions
    void set_circle (double r) {_radius = r + _pivot.vlength() ; }
  };

  inline Plane::Plane (ConstSharedPointer<GeomObject>* gop) 
  {
    const char funame [] = "DivSur::Plane::Plane (ConstSharedPointer<GeomObject>*): "; 

    _radius = -1.;

    const GeomPlane* pln;
    const GeomPivot* pvt;

    for(int frag = 0; frag < 2; ++frag) {
      pln = dynamic_cast<const GeomPlane*>((const GeomObject*)gop[frag]);
      pvt = dynamic_cast<const GeomPivot*>((const GeomObject*)gop[1-frag]);
      if(pln && pvt) {
	_plane = *pln;
	_pivot = *pvt;
	_plane_frag = frag;
	return;
      }
    }
    std::cerr << funame << "no appropriate objects found\n";
    throw Error::Form();
  }

  /********************************************************************************
   *   Set of primitive surfaces. It is a common base for all dividing surfaces   *
   ********************************************************************************/

  class PrimSet
  {
    std::vector<SharedPointer<Primitive> > _sur;  // analytical surfaces

  protected:
    std::map<std::string, int> _fmap;              // find surface by name
    
  public:
    void random_orient (int face, Dynamic::Coordinates& dc)  const { _sur[face]->random_orient(dc); }

    double t_velocity (int face, double temper, Dynamic::Vars& dv) const 
    { return _sur[face]->t_velocity(temper, dv); }
    double e_velocity (int face, double ener, Dynamic::Vars& dv) const 
    { return _sur[face]->e_velocity(ener, dv); }
    double j_velocity (int face, double ener, double amom, Dynamic::Vars& dv) const 
    { return _sur[face]->j_velocity(ener, amom, dv); }

    bool   test     (int face, const Dynamic::Coordinates& dc) const { return _sur[face]->test(dc); } 
    double weight   (int face, const Dynamic::Coordinates& dc) const { return _sur[face]->weight(dc); }
    double distance (int face, const Dynamic::Coordinates& dc) const { return _sur[face]->distance(dc); }

    int size () const { return _sur.size(); }

    void read (std::istream&) ;

    void print (std::ostream&, const std::string&) const;

    virtual ~PrimSet () {}
  };
    
  /********************************************************************************
   *              Basic dividing surface separating configurational               *
   *              space into two regions: the reactants and products              *
   ********************************************************************************/

  class OneSur : public PrimSet, public Dynamic::Condition, public Dynamic::Classifier, public IO::Read
  {
    bool _isinit;

    // product definition
    SharedPointer<Logical::Expr> _prod;

  public:
    OneSur ()                   : _isinit(false) {}
    OneSur (std::istream& from) : _isinit(false) { read(from); }
    virtual ~OneSur () {}

    bool isinit () const { return _isinit; }

    struct SmpRes {
      enum Stat {FACET, INNER, CLOSE, EXCLUDE};
      // does the facet inner part of the configurational 
      // space come into the product
      bool inner;
      Stat stat;
    };
  
    SmpRes facet_test (int face, const Dynamic::Coordinates&) const ;

    // is the configuration in the products state
    bool test (const Dynamic::Coordinates&) const ;

    // which surface did the trajectory crossed
    int  classify (const Dynamic::Coordinates& dv) const ;

    // read surface definition from the stream
    void read (std::istream&) ;

    void print (std::ostream&, const std::string&) const;
  }; // OneSur

  /*************************************************************************************************
   **************** DIVIDING SURFACE BETWEEN MULTIPLE BOUND AND BIMOLECULAR SPECIES ****************
   *************************************************************************************************/

  typedef std::pair<int, int> face_t;

  class MultiSur : private PrimSet, public Dynamic::Classifier, public IO::Read
  {
    bool _isinit;

    // definitions of bound species 
    std::vector<SharedPointer<Logical::Expr> > _species;

  public:
    MultiSur ()                   : _isinit(false) {}
    MultiSur (std::istream& from) : _isinit(false) { read(from); }
    ~MultiSur () {}

    bool isinit () const { return _isinit; }

    // surface random sampling result
    struct SmpRes {
      enum Stat {FACET, INNER, CLOSE, EXCLUDE, FAIL};
      face_t face;
      Stat stat;
    };

    // number of bound and bimolecular species
    int species_size () const { return _species.size() + 1; }

    // which species the configuration belongs to
    int classify (const Dynamic::Coordinates&) const ;

    // does the configuration belong to a given species region
    bool species_test (int spec, const Dynamic::Coordinates&) const ;

    // number of primitives
    int primitive_size () const { return   PrimSet::size(); }

    // find the facet, if any, for a given configuration
    SmpRes facet_test (int sur, const Dynamic::Coordinates&) const ;

    // random configuration on the primitive surface
    void random_orient (int sur, Dynamic::Coordinates& dc)  const { PrimSet::random_orient(sur, dc); }

    // generate random thermal velocity
    double t_velocity (int sur, double temper, Dynamic::Vars& dv) const  
    { return PrimSet::t_velocity(sur, temper, dv); }

    // generate random velocity with given energy 
    double e_velocity (int sur, double ener, Dynamic::Vars& dv) const 
    { return PrimSet::e_velocity(sur, ener, dv); }

    // generate random velocity with given energy and angular momentum in z direction
    double j_velocity (int sur, double ener, double amom, Dynamic::Vars& dv) const 
    { return PrimSet::j_velocity(sur, ener, amom, dv); }

    // primitive test
    bool primitive_test (int sur, const Dynamic::Coordinates& dc) const { return PrimSet::test(sur, dc); }

    // statistical weight 
    double weight (int sur, const Dynamic::Coordinates& dc) const { return PrimSet::weight(sur, dc); }

    // distance to the primitive surface
    double primitive_distance (int sur, const Dynamic::Coordinates& dc) const { return PrimSet::distance(sur, dc); }

    // read the surface from the input stream
    void read (std::istream&) ;

    void print (std::ostream&, const std::string&) const;
  }; // MultiSur

}// DivSur namespace

#endif
