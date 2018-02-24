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
// testLDAFunctional.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

// Test the LDA functional by computing the xc energy of a gaussian
// of width 0.1 a.u. in a cube of side 1.0 a.u.
// With a cube of side 1.0 and 32x32x32 points, 
// The xc energy must be -2.8105 a.u.
// dExc/da must be 0.911682

#include<iostream>
#include "LDAFunctional.h"
#include <cassert>
#include <cmath>
using namespace std;

int main(int argc, char **argv)
{
  // use: testxcf alat np
  if ( argc != 3 )
  {
    cout << " use: testLDAFunctional alat np" << endl;
    return 0;
  }
  assert(argc==3);
  double a = atof(argv[1]);
  double omega = a*a*a;
  int n = atoi(argv[2]);
  int n3 = n*n*n;
  double *rh = new double[n3];
  double *vxc = new double[n3];
  double *exc = new double[n3];
  double excsum = 0.0, dxcsum = 0.0;
  
  double rc = 0.1 * a;
  double pi = 4.0 * atan(1.0);
  double fac = 1.0 / ( pow(pi,1.5) * rc*rc*rc );
  double sum = 0.0;
  
  for ( int i = 0; i < n; i++ )
  {
    double x = ( i * a ) / n - a/2;
    for ( int j = 0; j < n; j++ )
    {
      double y = ( j * a ) / n - a/2;
      for ( int k = 0; k < n; k++ )
      {
        double z = ( k * a ) / n - a/2;
        double r2 = x*x + y*y + z*z;
        int ii = i + n * ( j + n * k );
        rh[ii] = fac * exp( -r2 / (rc*rc) );
        sum += rh[ii];
      }
    }
  }
  sum = sum * omega / n3;   
  // the density should be normalized
  cout << " Integrated density: " << sum << endl;  
      
  LDAFunctional xcf;
  
  int nspin = 1;
  xcf.rho = rh;
  xcf.vxc1 = vxc;
  xcf.exc = exc;
  
  xcf.setxc(n3,nspin);
  
  for ( int i = 0; i < n3; i++ )
    excsum += rh[i] * exc[i];
  for ( int i = 0; i < n3; i++ )
    dxcsum += rh[i] * ( exc[i] - vxc[i] );
  
  cout << " Total LDA xc energy: " << excsum * omega / n3 << endl;
  
  // Note: the energy variation is 3 * dExc/da * delta(a)
  cout << " dExc/da: " << dxcsum * omega / ( n3 * a ) << endl;

  return 0;
}
