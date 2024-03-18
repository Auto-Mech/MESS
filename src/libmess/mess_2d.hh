#ifndef MESS_2D_HH
# define MESS_2D_HH

namespace Mess {
  
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
    
    int linear_index (int e, int a) const;

    int ener_index (int g) const;

    int amom_index (int g) const;

    int amom_shift (int a) const;

    int ener_size (int a) const;

    double states (int) const;

    double weight () const;

    double real_weight () const;

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
    
  extern std::vector<Well> well;

  extern std::vector<PesObject> inner_barrier;

  extern std::vector<PesObject> outer_barrier;

  inline double boltzmann_factor (int e) { return std::exp(energy_step / temperature * e); }
}

#endif
