//#undef INT

#include <mpi.h>

#include "comm.hh"

void Comm::send_rate_data(const std::map<std::pair<int, int>, double>& data)
{
  int    itemp;
  double dtemp;

  itemp = data.size();
  
  MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

  int pair[2];
  
  for(std::map<std::pair<int, int>, double>::const_iterator it = data.begin(); it != data.end(); ++it)  {
    //
    pair[0] = it->first.first;
    
    pair[1] = it->first.second;

    MPI_Send(pair, 2, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

    dtemp = it->second;
    
    MPI_Send(&dtemp, 1, MPI_DOUBLE, 0, DATA_TAG, MPI_COMM_WORLD);
  }
}

void Comm::recv_rate_data(int node,  std::map<std::pair<int, int>, double>& data)
{
  int    itemp;
  double dtemp;

  int n;

  MPI_Status stat;
  
  MPI_Recv(&n, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  int pair[2];

  data.clear();
  
  for(int i = 0; i < n; ++i)  {
    //
    MPI_Recv(pair, 2, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    MPI_Recv(&dtemp, 1, MPI_DOUBLE, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    data[std::make_pair(pair[0], pair[1])] = dtemp;
  }
}
  
void Comm::send_capture_data(const std::map<int, double>& data)
{
  int    itemp;
  double dtemp;

  itemp = data.size();
  
  MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

  for(std::map<int, double>::const_iterator it = data.begin(); it != data.end(); ++it)  {
    //
    itemp = it->first;
    
    MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

    dtemp = it->second;
    
    MPI_Send(&dtemp, 1, MPI_DOUBLE, 0, DATA_TAG, MPI_COMM_WORLD);
  }
}

  
void Comm::recv_capture_data(int node,  std::map<int, double>& data)
{
  int    itemp;
  double dtemp;

  int n;

  MPI_Status stat;
  
  MPI_Recv(&n, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  data.clear();
  
  for(int i = 0; i < n; ++i) {
    //
    MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    MPI_Recv(&dtemp, 1, MPI_DOUBLE, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    data[itemp] = dtemp;
  }
}
  
void Comm::send_partition_data(const MasterEquation::Partition& part)
{
  int    itemp;
  double dtemp;

  itemp = part.size();

  MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

  for(int i = 0; i < part.size(); ++i) {
    //
    itemp = part[i].size();
    
    MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

    Array<int> g(part[i]);
    
    MPI_Send(g, itemp, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);
  }
}

void Comm::recv_partition_data(int node,  MasterEquation::Partition& part)
{
  int    itemp;
  double dtemp;

  int n;

  MPI_Status stat;
  
  MPI_Recv(&n, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  part.clear();

  part.resize(n);
  
  for(int g = 0; g < n; ++g) {
    //
    MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    Array<int> a(itemp);
   
    MPI_Recv(a, a.size(), MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);
    
    for(int i = 0; i < a.size(); ++i)
      //
      part[i].insert(a[i]);
  }
}

void Comm::send_log ()
{
  int itemp;

  itemp = IO::log.str().size();
  
  MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

  MPI_Send(IO::log.str().c_str(), itemp, MPI_CHAR, 0, DATA_TAG, MPI_COMM_WORLD);

  IO::log.clear();

  itemp = IO::aux.str().size();
  
  MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

  MPI_Send(IO::aux.str().c_str(), itemp, MPI_CHAR, 0, DATA_TAG, MPI_COMM_WORLD);

  IO::aux.clear();

  
}

void Comm::recv_log(int node)
{
  int itemp;

  MPI_Status stat;
  
  MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  Array<char> msg(itemp + 1);

  MPI_Recv(msg, itemp, MPI_CHAR, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  msg.back() = 0;

  IO::log << (char*)msg;

  MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  msg.resize(itemp + 1);

  MPI_Recv(msg, itemp, MPI_CHAR, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  msg.back() = 0;

  IO::aux << (char*)msg;
}

  
