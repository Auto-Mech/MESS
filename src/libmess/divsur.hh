/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2025, Yuri Georgievski <ygeorgi@anl.gov>

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
  //
  class GeomObject
  {
  protected:
    //
    int _type;
    
  public:
    //
    virtual ~GeomObject () {}

    enum {PIVOT, PLANE, CONE};
    
    int type () const { return _type; }
  };

  class GeomPlane : public GeomObject, public D3::Plane
  {
  public:
    //
    GeomPlane  (const D3::Vector& n, double d) : D3::Plane(n, d) { _type = PLANE; }
    
    ~GeomPlane () {}
  };
  
  class GeomPivot : public GeomObject, public D3::Vector
  {
  public:
    //
    GeomPivot (const D3::Vector& v) : D3::Vector(v) { _type = PIVOT; }
    
    ~GeomPivot () {}
  };

  class GeomCone : public GeomObject {
    //
  public:
    //
    GeomCone () : radius(-1.) { _type = CONE; }
    
    D3::Vector head;
    D3::Vector axis;
    double     proj;
    double     radius;
    
    ~GeomCone () {}
  };
  
  // auxiliary class to read the geometrical objects associated with a fragment
  //
  struct FragData : public std::map<std::string, ConstSharedPointer<GeomObject> > {
    //
    void read (int frag, std::istream&);
  };
  
  /*********************************************************************
   * Abstract elemental surface: the concept of the analytical surface *
   * in the configurational space. It should have three features:      *
   *   1.) to allow randomly orient fragments on its surface           *
   *   2.) to provide a statistical weight of such sampling            *
   *   3.) separate the configurational space into two regions         *
   *********************************************************************/
    
  class Primitive {
    //
  protected:
    //
    // Normalized generalized reaction coordinate vector and its length, should be redefined for each real surface
    //
    virtual double _set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_vec) const; 

    // generate random velocities from the normalized reaction coordinate vector
    //
    static void   _rc2tv (const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, double tmpr,              Dynamic::Momenta& dm);// thermal
    static double _rc2ev (const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, double ener,              Dynamic::Momenta& dm);// E-resolved
    static double _rc2jv (const Dynamic::Momenta& rc_vec, const Dynamic::Coordinates& dc, double ener, double amom, Dynamic::Momenta& dm);// E,J-resolved
    
  public:
    //
    // Randomly orient the fragments on the surface, should be redefined for each surface
    //
    virtual void   random_orient   (Dynamic::Coordinates& dv) const =0;

    // signed distance of the configuration from the surface; should be re-defined for each real surface
    //
    virtual double distance  (const Dynamic::Coordinates& dv) const;

    // output the surface; should be redefined for each real surface
    //
    virtual void print (std::ostream&) const = 0;

    // test if the configuration is inside of the surface; inside corresponds to the positive distance
    //
    bool test (const Dynamic::Coordinates& dv) const;

    // statistical weight of the configuration on the surface
    //
    double weight (const Dynamic::Coordinates& dv) const;

    /*
      randomly generate velocities for  the given configuration 
      on the dividing surface. The velocities are directed inside the surface.
      The output is the statistical weight of the generated configuration
    */
    double t_velocity (double temper,            Dynamic::Vars& dv) const;// thermal
    double e_velocity (double ener,              Dynamic::Vars& dv) const;// E-resolved
    double j_velocity (double ener, double amom, Dynamic::Vars& dv) const;// E,J-resolved

    virtual ~Primitive () {}
  };

  inline bool Primitive::test  (const Dynamic::Coordinates& dv) const 
  {
    if(distance(dv) > 0.) {
      //
      return true; 
    }
    else 
      //
      return false; 
  }
    
  /************************************************************************
   *        Spherical surface: the distance between two pivot points        *
   *         associated with the corresponding fragments is fixed         *
   ************************************************************************/

  class Sphere: public Primitive {
    //
    D3::Vector _pivot [2];
    double _dist;

    // lab frame pp2pp vector, 1->2
    //
    D3::Vector _lf_pp_12 (const Dynamic::Coordinates& dv) const; 
    double _set_rc (const Dynamic::Coordinates& dv, Dynamic::Momenta& rc_vec) const; 

  public:
    //
    Sphere (ConstSharedPointer<GeomObject>*, double);
    ~Sphere () {}

    double distance  (const Dynamic::Coordinates& dv) const { return _dist - _lf_pp_12(dv).vlength(); }
    
    void   random_orient   (Dynamic::Coordinates& dv) const;

    void print (std::ostream&) const;

    // specific functions
    //
    double max_cm_dist () const { return _dist + _pivot[0].vlength() + _pivot[1].vlength(); }
  };
    
  inline Sphere::Sphere (ConstSharedPointer<GeomObject>* gop, double d) 
  {
    const char funame [] = "DivSur::Sphere::Sphere: "; 

    _dist = d;
    const GeomPivot* pvt;
    for(int frag = 0; frag < 2; ++frag) {
      //
      pvt = dynamic_cast<const GeomPivot*>((const GeomObject*)(gop[frag]));
      
      if(!pvt) {
	//
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
    //
    D3::Plane  _plane;
    D3::Vector _pivot;
    int _plane_frag; // on which fragment plane is

    double _radius; // sampling circle radius

    double _set_rc (const Dynamic::Coordinates& dv, Dynamic::Momenta& rc_vec) const; 

  public:
    //
    Plane (ConstSharedPointer<GeomObject>*);
    ~Plane () {}
  
    double distance  (const Dynamic::Coordinates& dv) const;
    void   random_orient   (Dynamic::Coordinates& dv) const;
    void   print           (std::ostream&) const;

    // specific functions
    //
    void set_circle (double r) {_radius = r + _pivot.vlength() ; }
  };
    
  inline Plane::Plane (ConstSharedPointer<GeomObject>* gop) 
  {
    const char funame [] = "DivSur::Plane::Plane: "; 

    _radius = -1.;

    const GeomPlane* pln;
    const GeomPivot* pvt;

    for(int frag = 0; frag < 2; ++frag) {
      //
      pln = dynamic_cast<const GeomPlane*>((const GeomObject*)gop[frag]);
      pvt = dynamic_cast<const GeomPivot*>((const GeomObject*)gop[1-frag]);
      
      if(pln && pvt) {
	//
	_plane_frag = frag;
	
	if(Structure::type(_plane_frag) ==  Molecule::MONOATOMIC) {
	  //
	  std::cerr << funame << "monoatomic molecule cannot have cone elements\n";
    
	  throw Error::Logic();
	}

	_plane      = *pln;
	_pivot      = *pvt;

	return;
      }
    }
    
    std::cerr << funame << "no appropriate objects found\n";
    throw Error::Form();
  }

  /********************************************************************************
   *        Plane surface: pivot point of the second fragment belongs to the        *
   *            geometrical plane associated with the first fragment              *
   ********************************************************************************/
    
  class Cone : public Primitive {
    //
    D3::Vector _pivot;     // non-cone fragment pivot point
    D3::Vector _head;      // cone head
    D3::Vector _axis;      // cone axis
    double     _proj;      // cone angle cosine
    double     _radius;    // sampling radius
    
    int        _cone_frag; // on which fragment cone is

    double _set_rc (const Dynamic::Coordinates& dc, Dynamic::Momenta& rc_vec) const; 

  public:
    //
    Cone (ConstSharedPointer<GeomObject>*);
    ~Cone () {}
  
    double distance  (const Dynamic::Coordinates& dv) const;
    void   random_orient   (Dynamic::Coordinates& dv) const;

    void print (std::ostream&) const;
  };
    
  inline Cone::Cone (ConstSharedPointer<GeomObject>* gop) 
  {
    const char funame [] = "DivSur::Cone::Cone: "; 

    const GeomCone*  cp;
    const GeomPivot* pp;

    for(int frag = 0; frag < 2; ++frag) {
      //
      cp = dynamic_cast<const GeomCone*>((const GeomObject*)gop[frag]);
      pp = dynamic_cast<const GeomPivot*>((const GeomObject*)gop[1 - frag]);
      
      if(cp && pp) {
	//
	_cone_frag  = frag;

	if(Structure::type(_cone_frag) ==  Molecule::MONOATOMIC) {
	  //
	  std::cerr << funame << "monoatomic molecule cannot have cone elements\n";
    
	  throw Error::Logic();
	}

	_pivot      = *pp;
	_head       = cp->head;
	_axis       = cp->axis;
	_proj       = cp->proj;
	_radius     = cp->radius;
	
	return;
      }
    }
    std::cerr << funame << "no appropriate objects found\n";
    throw Error::Form();
  }
    
  /*************************************************************************************************
   **************** DIVIDING SURFACE BETWEEN MULTIPLE BOUND AND BIMOLECULAR SPECIES ****************
   *************************************************************************************************/
    
  typedef std::pair<int, int> face_t;

  // surface sampling result
  //
  enum {FACET, INNER, CLOSE, EXCLUDE, FAIL};
  
  struct SmpRes {
    //
    face_t face;
    int stat;
  };

  class MultiSur : public Dynamic::Condition, public Dynamic::Classifier, public IO::Read
  {
    bool _isinit;

    // primitive surfaces definitions
    //
    std::vector<SharedPointer<Primitive> > _prim;
    
    // definitions of bound species
    //
    std::vector<SharedPointer<Logical::Expr> > _species;

  public:
    //
    MultiSur ()                   : _isinit(false) {}
    MultiSur (std::istream& from) : _isinit(false) { read(from); }
    ~MultiSur () {}

    bool isinit () const { return _isinit; }

    // number of bound and bimolecular species
    //
    int species_size () const { return _species.size() + 1; }

    // number of primitives
    //
    int primitive_size () const { return _prim.size(); }

    // which species the configuration belongs to
    //
    int classify (const Dynamic::Coordinates&) const;

    // is the configuration not bimolecular
    //
    bool test (const Dynamic::Coordinates& dc) const { if(classify(dc) != _species.size()) return true; return false; }

    // does the configuration belong to a given species
    //
    bool test (int spec, const Dynamic::Coordinates& dc) const { if(classify(dc) == spec) return true; return false; }

    // find the facet for a given configuration
    //
    SmpRes surface_test (int sur, const Dynamic::Coordinates&) const;

    void random_orient (int sur, Dynamic::Coordinates& dc) const { _prim[sur]->random_orient(dc); }

    double t_velocity (int sur, double tmpr,              Dynamic::Vars& dv) const { return _prim[sur]->t_velocity(tmpr,       dv); }
    double e_velocity (int sur, double ener,              Dynamic::Vars& dv) const { return _prim[sur]->e_velocity(ener,       dv); }
    double j_velocity (int sur, double ener, double amom, Dynamic::Vars& dv) const { return _prim[sur]->j_velocity(ener, amom, dv); }

    double weight    (int sur, const Dynamic::Coordinates& dc) const { return _prim[sur]->weight(dc); }

    bool   prim_test (int sur, const Dynamic::Coordinates& dc) const { return _prim[sur]->test(dc); } 
    double distance  (int sur, const Dynamic::Coordinates& dc) const { return _prim[sur]->distance(dc); }

    // read the surface from the input stream
    //
    void read (std::istream&);

    void print (std::ostream&, const std::string&) const;
  };
    
}// DivSur namespace

#endif
