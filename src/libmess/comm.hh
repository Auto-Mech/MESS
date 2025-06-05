#ifndef COMM_HH
#define COMM_HH

#include "mess.hh"

namespace Comm {
  //
  enum {STOP_TAG, RUN_TAG, DATA_TAG, FAIL_TAG};
  
  void send_rate_data(const std::map<std::pair<int, int>, double>&);

  void recv_rate_data(int,  std::map<std::pair<int, int>, double>&);
  
  void send_capture_data(const std::map<int, double>&);
  
  void recv_capture_data(int,  std::map<int, double>&);
  
  void send_partition_data(const MasterEquation::Partition&);

  void recv_partition_data(int,  MasterEquation::Partition&);

  void send_log ();

  void recv_log (int);
}

#endif
