////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//

#include <config.h>

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
  
  
  
