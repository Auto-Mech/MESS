#ifndef MESS_2D_HH
# define MESS_2D_HH

namespace Mess {
  //
  class ReactiveComplex: public Model::ChemGraph {
  };
  
  extern double energy_reference;

  extern double energy_step;

  extern double amom_step;

  extern double amom_max;
  
  extern double temperature;

  // discretized states density/number
  // 
  class PesObject {
    //
  public:
    //
    int size () const;
  
    int amom_size () const;
    
    int ener_size (int a) const;

    int linear_index (int e, int a) const;

    int ener_index (int g) const;

    int amom_index (int g) const;

    int amom_shift (int a) const;

    double states (int) const;

    double weight () const;

    double real_weight () const;

    const std::string& name () const;

  private:
    //
    // ...
  };

  class Well: public PesObject {
    //
  public:
    //
    explicit Well (const Model::Well&);
    
    double collision_frequency () const;

    double kernel(int i, int j) const;
  };

  class Barrier: public PesObject {
    //
  public:
    //
    std::pair<int, int> connect;
  };
    
  extern std::vector<Well> well;

  extern std::vector<Barrier> inner_barrier;

  extern std::vector<Barrier> outer_barrier;

  inline double boltzmann_factor (int e) { return std::exp(energy_step / temperature * e); }
}

#endif
