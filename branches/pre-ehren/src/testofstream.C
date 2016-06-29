#include <fstream>
#include <iostream>
#include <cstring>
using namespace std;
#include "Context.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

int main()
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
  
  Context ctxt;
  
  unsigned long long kb = 1024;
  unsigned long long mb = kb * kb;
  unsigned long long gb = mb * kb;
  
  unsigned char* buf = new unsigned char[mb];
  
  memset(buf, 0, mb);
  
  int write_count = four_gb / mb;
  
  cout << " mype=" << ctxt.mype() << endl;
  
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}  
  
  
  
