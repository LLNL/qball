#include <mpi.h>
#include <iostream>
using namespace std;
int MPI_Send(void* buf, int count, int datatype, int dest, int tag,
  int comm)
{
  cout << "myMPI_Send: buf=" << buf << " size=" << count << endl;
  return PMPI_Send(buf,count,datatype,dest,tag,comm);
}
