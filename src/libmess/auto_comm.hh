#ifndef COMM_HH
#define COMM_HH

#include "new_mess.hh"

namespace Comm {
  //
  enum {STOP_TAG, RUN_TAG, DATA_TAG, FAIL_TAG};
  
  void send_rate_data(const std::list<MasterEquation::RateData>&);

  void recv_rate_data(int,  std::list<MasterEquation::RateData>&);
  
  void send_log ();

  void recv_log (int);
}

#endif
