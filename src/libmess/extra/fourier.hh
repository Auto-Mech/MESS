

#ifndef FOURIER_HH
#define FOURIER_HH

// multidimensional fourier expansion
class FourierExpansion {
  int _dim;
  std::map<std::vector<int>, double > _data;
  void _check (const std::vector<int>&) throw(Error::General);

public:

  void add (const std::vector<int>&, double) throw(Error::General);
};

#endif
