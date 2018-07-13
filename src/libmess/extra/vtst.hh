

#ifndef CTST_HH
#define CTST_HH

namespace CTST { // classical variational transition state theory

  // an abstract class to represent a number of states density to average over fragment orientations 
  class N_Den {
    double _dist; // distance between the fragments
    double _ener; // energy
  public:
    double distance () const { return _dist; }
    double energy   () const { return _ener; }

    N_Den (double d, double e) : _dist(d), _ener(e) {}
    virtual double operator () const (const Dynamic::Coordinates&) = 0;
  };

  class DenWrap
  {
    ConstSharedPointer<NDen> _den;
	
  public:
    DenWrap () {}
    DenWrap (ConstSharedPointer<N_Den> p) : _den(pf) {}

    void   isinit     () const throw(Error::General);
    double operator() (const Dynamic::Coordinates&) const throw(Error::General);
    ~Wrap () {}
  };

  inline void DenWrap::isinit () const throw(Error::General)
  {
    const char funame [] = "CTST::DenWrap::isinit: ";

    if(!_den) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }
  }

  inline double DenWrap::operator() (const Dynamic::Coordinates& dc) const throw(Error::General)
  {
    isinit();
    return (*_den)(dc);
  }

  // number of states of the fragments at a distance R
  class Ns_Den : N_Den  {

  public:
    Ns_Den (double d, double e) : N_Den(d, e) {}
    double operator () const (const Dynamic::Coordinates&);
  };

  // number of states of the fragments at a distance R resolved over the projection of 
  // the angular momentum on the interfragments axis
  class NsK_Den : public N_Den {
    double _proj; // projection of the angular momentum on the interfragment axis

  public:
    double angular_momentum_projection () const { return _proj; }

    NsK_Den (double d, double e, double k) : N_Den (d, e), _proj(k) {}
    double operator () const (const Dynamic::Coordinates&);
  };

  // total number of states (fragments + orbital) , ej_resolved
  class Nt_Den : N_Den  {
    double _amom; // total angular momentum

  public:
    double angular_momentum () const { return _amom; }

    Nt_Den (double d, double e, double j) : N_Den(d, e), _amom(j) {}
    double operator () const (const Dynamic::Coordinates&);
  };

  // total number of states (fragments + orbital) , ejk_resolved
  class NtK_Den : N_Den  {
    double _amom; // total angular momentum
    double _proj; // total angular momentum projection on the interfragment axis

  public:
    double angular_momentum            () const { return _amom; }
    double angular_momentum_projection () const { return _proj; }
    
    NtK_Den (double d, double e, double j, double k) : N_Den(d, e), _amom(j), _proj(k) {}
    double operator () const (const Dynamic::Coordinates&);
  };


}

#endif
