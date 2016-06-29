
#if LINUX
#define _LARGEFILE_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#elif AIX
#define _LARGE_FILES 1
#define _LARGE_FILE_API 1
#endif


#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
using namespace std;
#include "Context.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
  
  Context ctxt;
  
  unsigned long long kb = 1024;
  unsigned long long mb = kb * kb;
  unsigned long long gb = mb * kb;
  unsigned long long four_gb = 4 * gb;
  unsigned long long ten_gb = 10 * gb;
   cout << "kb:       " << kb
        << "\nmb:       " << mb
        << "\ngb:       " << gb
        << "\nfour_gb:  " << four_gb
        << "\nten_gb:   " << ten_gb
        << endl;
  
  cout << " sizeof(unsigned long long) = " 
       << sizeof(unsigned long long) << endl;
  cout << " sizeof(size_t) = " << sizeof(size_t) << endl;
  cout << " sizeof(string::size_type) = " << sizeof(string::size_type) << endl;
  cout << " sizeof(std::streampos) = " << sizeof(std::streampos) << endl;
  cout << " mype=" << ctxt.mype() << endl;
  
  unsigned long long mysize = mb;
  char* buf = new char[mysize];
  memset(buf, 0, mysize);
  unsigned long long myoffset = mysize * ctxt.mype();
  
  ofstream os("test_file");
  if ( ctxt.oncoutpe() )
  {
    os << "header information" << endl;
  }
  streampos current_offset = os.tellp();
  cout << ctxt.mype() << ": current_offset: " << current_offset << endl;
  os.seekp(myoffset);
  os.write(&buf[0],mysize);
  os.close();
  
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}  
  
  
  
