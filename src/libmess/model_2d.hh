#ifndef MODEL_2D_HH
#define MODEL_2D_HH

#include "model.hh"

namespace Model {

  class Core2 {
    //
    int _mode;
    
    Core2 ();

  protected:
    //
    explicit Core2 (int m);

  public:
    //
    int mode () const { return _mode; }
    
    virtual double ground (double j) const =0;

    virtual void states    (Array<double>&, double de, double j) const =0;

    virtual void convolute (Array<double>&, double de, double j) const =0;
    
    virtual double weight (double t) const =0;
  };
  
  inline Core2::Core2 (int m) : _mode(m) 
  {
    const char funame [] = "Model::Core2::Core2: ";

    if(mode() != NUMBER && mode() != DENSITY && mode() != NOSTATES) {
      //
      std::cerr << funame << "wrong mode\n";
      
      throw Error::Logic();
    }
  }

  class RigidRotor2: public Core2 {
    //
    int _type;

    double _symmetry_factor;
    
    std::vector<double> _rc;

  public:
    //
    enum {LINEAR_ROTOR, SPHERICAL_TOP, PROLATE_TOP, OBLATE_TOP, ASYMMETRIC_TOP};

    static double eps;

    RigidRotor2 (IO::KeyBufferStream&, const std::vector<Atom>&, int);
    
    double ground (double) const;

    void states    (Array<double>&, double de, double j) const;

    void convolute (Array<double>&, double de, double j) const;

    double weight (double t)           const;
  };

}// namespace Model

#endif
