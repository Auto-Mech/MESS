

#ifndef THERMO_POLY_HH
#define THERMO_POLY_HH

class ThermoPoly {
  //
  typedef std::vector<double>                     _pos_t;

  typedef std::map<std::multiset<int>, double> _potex_t ;

  _potex_t _global_potex;

  std::vector<double> _frequency;

  void _read_potex (std::ifstream&, _potex_t&);

  _potex_t _local_potex (const _pos_t&, int) const;

  double   _pot (const _pos_t&) const;

public:

  void init (std::istream&);

  std::map<int, double> weight (double temperature) const;
 
};
