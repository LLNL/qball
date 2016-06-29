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
#include <dvec.h>
#include <complex>
#include <iostream>
using namespace std;

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

void h(int n, complex<double> *x, complex<double> *y, complex<double>*z)
{
  F64vec2 *xx = (F64vec2*) x;
  F64vec2 *yy = (F64vec2*) y;
  F64vec2 *zz = (F64vec2*) z;
  #pragma unroll(0)
  for ( int i = 0; i < n; i++ )
  {
    F64vec2 ab(xx[i]);
    F64vec2 ba(ab[1],ab[0]);
    F64vec2 cd(yy[i]);
    F64vec2 cc(cd[0]);
    F64vec2 mdd(-cd[1],cd[1]);
    zz[i] = ab * cc + ba * mdd;
  }
}

void g(int n, complex<double> *a, complex<double> *b, complex<double>*c)
{
  #pragma unroll(0)
  for ( int i = 0; i < n; i++ )
  {
    c[i] = a[i] * b[i];
  }
}

void f(int n, double *a, complex<double> *b, complex<double> *c)
{
  // c[i] = a[i] * b[i]
  #pragma unroll(4)
  for ( int i = 0; i < n; i++ )
  {
    c[i] = a[i] * b[i];
  }
}

void ff(int n, double *a, complex<double> *b, complex<double> *c)
{
  // c[i] = a[i] * b[i]
  F64vec2 *bb = (F64vec2*) b;
  F64vec2 *cc = (F64vec2*) c;
  #pragma unroll(4)
  for ( int i = 0; i < n; i++ )
  {
    F64vec2 aa(a[i]); 
    cc[i] = aa * bb[i];
    //_mm_prefetch((char*)&b[i+8],_MM_HINT_T0);
  }
}

int main(int argc, char **argv)
{
  int n = atoi(argv[1]);
  
  // allocate 16-byte aligned arrays
  double *a = (double*) _mm_malloc(n*sizeof(double),16);
  complex<double> *b = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *c = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  double *aa = (double*) _mm_malloc(n*sizeof(double),16);
  complex<double> *bb = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  complex<double> *cc = (complex<double>*)_mm_malloc(n*sizeof(complex<double>),16);
  
  // initialize arrays
  for ( int i = 0; i < n; i++ )
  {
    a[i] = (double) i;
    b[i] = complex<double>(0.5*i,i*i);
    aa[i] = a[i];
    bb[i] = b[i];
  }
  
  long long t0;
  long long tf[10],tff[10];
  
  for ( int i = 0; i < 10; i++ )
  {
    t0 = readTSC();
    f(n,a,b,c);
    tf[i] = readTSC()-t0;
    
    t0 = readTSC();
    ff(n,aa,bb,cc);
    tff[i] = readTSC()-t0;
  }
  
  // check correctness
  for ( int i = 0; i < n; i++ )
  {
    assert(norm(c[i]-cc[i]) == 0.0 );
  }
  
  cout << " cycles per point:" << endl;
  for ( int i = 0; i < 10; i++ )
    cout << "f: " << ((double)tf[i])/n 
         << "   ff: " << ((double)tff[i])/n << endl;

  return 0;
}
