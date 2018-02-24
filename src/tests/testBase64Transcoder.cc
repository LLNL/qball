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
////////////////////////////////////////////////////////////////////////////////
//
// testBase64Transcoder.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Base64Transcoder.h"
#include <iostream>
using namespace std;

int main()
{
  const int n = 7;
  double a[n] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  double c[n] = { 0, 0, 0, 0, 0, 0, 0 };
  
  Base64Transcoder xcdr;
  
  int nbytes = n * sizeof(double);
  
  int nchars = xcdr.nchars(nbytes);
  
  cout << " nbytes=" << nbytes << endl;
  cout << " nchars=" << nchars << endl;
  
  char* b = new char[nchars];
  
  xcdr.encode(nbytes,(unsigned char*) &a[0],b);
  
  cout << " b=" << b << endl;
  
  xcdr.decode(nchars,b,(unsigned char*) &c[0]);
  
  for ( int i = 0; i < n; i++ )
    assert(a[i]==c[i]);
    
  cout << " done" << endl;
  
  return 0;
}
