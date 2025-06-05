#include "new_comm.hh"

#include <mpi.h>

void Comm::send_rate_data(const std::list<MasterEquation::RateData>& rate_data)
{
  int    itemp;
  double dtemp;

  itemp = rate_data.size();
  
  MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

  std::list<MasterEquation::RateData>::const_iterator rit;
  
  for(rit = rate_data.begin(); rit != rate_data.end(); ++rit) {
    //
    itemp = rit->index_bim_map.size();
    
    MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

    if(itemp) {
      //
      MPI_Send((const int*)rit->index_bim_map, itemp, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

      MPI_Send((const double*)rit->bb_rate, ((RefArr<double>&)(rit->bb_rate)).size(), MPI_DOUBLE, 0, DATA_TAG, MPI_COMM_WORLD);
    }
    
    itemp = rit->well_partition.size();

    MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

    if(rit->well_partition.size()) {
      //
      for(int i = 0; i < rit->well_partition.size(); ++i) {
	//
	Array<int> g(rit->well_partition[i]);

	itemp = g.size();

	MPI_Send(&itemp, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);

	MPI_Send((int*)g, itemp, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);
      }

      itemp = rit->ww_rate.size() * rit->ww_rate.size();
      
      MPI_Send((const double*)rit->ww_rate, itemp, MPI_DOUBLE, 0, DATA_TAG, MPI_COMM_WORLD);

      if(rit->index_bim_map.size()) {
	//
	itemp = rit->ww_rate.size() * rit->bb_rate.size();
	
	MPI_Send((const double*)rit->wb_rate, itemp, MPI_DOUBLE, 0, DATA_TAG, MPI_COMM_WORLD);

	MPI_Send((const double*)rit->bw_rate, itemp, MPI_DOUBLE, 0, DATA_TAG, MPI_COMM_WORLD);
      }
    }
  }
}

void Comm::recv_rate_data(int node, std::list<MasterEquation::RateData>& rate_data_list)
{
  int    itemp;
  double dtemp;

  MPI_Status stat;
  
  rate_data_list.clear();

  int n;
  
  MPI_Recv(&n, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

  for(int i = 0; i < n; ++i) {
    //
    MasterEquation::RateData rd;
    
    MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    if(itemp) {
      //
      rd.index_bim_map.resize(itemp);
    
      rd.bb_rate.resize(itemp);
    
      MPI_Recv((int*)rd.index_bim_map, itemp, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

      MPI_Recv((double*)rd.bb_rate, ((RefArr<double>&)(rd.bb_rate)).size(), MPI_DOUBLE, node, DATA_TAG, MPI_COMM_WORLD, &stat);
    }
    
    MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

    if(itemp) {
      //
      rd.well_partition.resize(itemp);

      rd.ww_rate.resize(itemp);

      for(int j = 0; j < rd.well_partition.size(); ++j) {
	//
	MPI_Recv(&itemp, 1, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

	Array<int> g(itemp);

	MPI_Recv((int*)g, itemp, MPI_INT, node, DATA_TAG, MPI_COMM_WORLD, &stat);

	for(int k = 0; k < g.size(); ++k)
	  //
	  rd.well_partition[j].insert(g[k]);
      }

      itemp =  rd.ww_rate.size() * rd.ww_rate.size();
      
      MPI_Recv((double*)rd.ww_rate, itemp, MPI_DOUBLE, node, DATA_TAG, MPI_COMM_WORLD, &stat);

      if(rd.index_bim_map.size()) {
	//
	rd.wb_rate.resize(rd.ww_rate.size(), rd.bb_rate.size());
      
	rd.bw_rate.resize(rd.bb_rate.size(), rd.ww_rate.size());

	itemp = rd.ww_rate.size() * rd.bb_rate.size();
	
	MPI_Recv((double*)rd.wb_rate, itemp, MPI_DOUBLE, node, DATA_TAG, MPI_COMM_WORLD, &stat);
	
	MPI_Recv((double*)rd.bw_rate, itemp, MPI_DOUBLE, node, DATA_TAG, MPI_COMM_WORLD, &stat);
      }
    }
    
    rate_data_list.push_back(rd);
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

  
