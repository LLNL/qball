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
// 
// Write large files ( > 2 GB ) on Linux
//

#define _LARGEFILE_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <cstdio>
#include <cstdlib>

#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
  FILE *outfile;
  outfile = fopen( "t.dat", "w" );
  
  const off_t n = atoi(argv[1]);
  const off_t m = 1024*1024/8;
  
  double a[m];
  
  cout << " sizeof(int): " << sizeof(int) << endl;
  cout << " sizeof(size_t): " << sizeof(size_t) << endl;
  cout << " sizeof(long): " << sizeof(long) << endl;
  cout << " sizeof(off_t): " << sizeof(off_t) << endl;
  
  for ( int i = 0; i < n; i++ )
    fwrite(&a[0],sizeof(double),m,outfile);
  
  fclose(outfile);
  
  outfile = fopen("t.dat","r");
  off_t offset = 8*m*n-1000;
  fseeko(outfile,offset,SEEK_SET);
  off_t pos = ftello(outfile);
  
  cout << offset << endl;
  cout << pos << endl;
  fclose(outfile);
}

