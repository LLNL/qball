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
//
// Test of Intel SSE2 C++ Class library
//
// Compute the real part of a complex scalar product
//
// Compile on Xeon with: icc -O3 -xW thisfile.C
//
// Run with: a.out 100
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>


#include <dvec.h>
#include <complex>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
long long readTSC(void)
{
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
            : "=mr" (ticks.part[0]),
              "=mr" (ticks.part[1])
            : /* no inputs */
            : "eax", "edx");
  return ticks.complete;
}

////////////////////////////////////////////////////////////////////////////////
double f(int n, complex<double> *a, complex<double> *b)
{
  double sum = 0.0;
  for ( int i = 0; i < n; i++ )
  {
    sum += a[i].real()*b[i].real() + a[i].imag()*b[i].imag();
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double ff(int n, complex<double> *a, complex<double> *b)
{
  F64vec2 *aa = (F64vec2*) a;
  F64vec2 *bb = (F64vec2*) b;
  F64vec2 cc(0.0); 
  #pragma unroll(4)
  for ( int i = 0; i < n; i++ )
  {
    cc += aa[i] * bb[i];
  }
  return add_horizontal(cc);
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  int n = atoi(argv[1]);
  
  // allocate 16-byte aligned arrays
  complex<double> *a = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *b = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *aa = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *bb = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  
  // initialize arrays
  for ( int i = 0; i < n; i++ )
  {
    a[i] = (double) i;
    b[i] = complex<double>(0.5*i,i);
    aa[i] = a[i];
    bb[i] = b[i];
  }
  
  long long t0;
  long long tf[10],tff[10];
  
  double sumf,sumff;
  for ( int i = 0; i < 10; i++ )
  {
    t0 = readTSC();
    sumf = f(n,a,b);
    tf[i] = readTSC()-t0;
    
    t0 = readTSC();
    sumff = ff(n,aa,bb);
    tff[i] = readTSC()-t0;
  }
  
  // check correctness
  assert(sumf==sumff);
    
  cout << " cycles per point:" << endl;
  for ( int i = 0; i < 10; i++ )
    cout << "f: "     << ((double)tf[i])/n 
         << "   ff: " << ((double)tff[i])/n << endl;

  return 0;
}
