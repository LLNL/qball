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
// testXCFunctional.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

// Test an XC functional by computing the xc energy of a gaussian
// of width 0.1 a.u. in a cube of side 1.0 a.u.
// With a cube of side 1.0 and 32x32x32 points, 
// The LDA xc energy must be -2.8105 a.u.
// dExc/da must be 0.911682

#include <iostream>
#include <vector>
#include "LDAFunctional.h"
#include "PBEFunctional.h"
#include "Timer.h"
#include <cassert>
#include <cmath>
using namespace std;

int main(int argc, char **argv)
{
  // use: testxcf alat np
  if ( argc != 3 )
  {
    cout << " use: testXCFunctional alat np" << endl;
    return 0;
  }
  assert(argc==3);
  double a = atof(argv[1]);
  double omega = a*a*a;
  int n = atoi(argv[2]);
  int n3 = n*n*n;
  vector<vector<double> > rh;
  rh.resize(1);
  rh[0].resize(n3);
  
  XCFunctional *xcf_list[2];
  xcf_list[0] = new LDAFunctional(rh);
  xcf_list[1] = new PBEFunctional(rh);
  
  for ( int ixcf = 0; ixcf < 2; ixcf++ )
  {
    Timer tm;
    XCFunctional *xcf = xcf_list[ixcf];
    cout << endl << " Functional name: " << xcf->name() << endl;
  
    double *grad_rho[3];
    if ( xcf->isGGA() )
    {
      grad_rho[0] = xcf->grad_rho[0];
      grad_rho[1] = xcf->grad_rho[1];
      grad_rho[2] = xcf->grad_rho[2];
    }
 
    const double rc = 0.1 * a;
    const double rc2 = rc * rc;
    const double pi = M_PI;
    const double fac = 1.0 / ( pow(pi,1.5) * rc*rc*rc );
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
          rh[0][ii] = fac * exp( -r2 / rc2 );
          sum += rh[0][ii];
 
          if ( xcf->isGGA() )
          {
            grad_rho[0][ii] = - rh[0][ii] * 2.0 * x / rc2;
            grad_rho[1][ii] = - rh[0][ii] * 2.0 * y / rc2;
            grad_rho[2][ii] = - rh[0][ii] * 2.0 * z / rc2;
          }
        }
      }
    }
    sum = sum * omega / n3;
    // the density should be normalized
    cout << " Integrated density: " << sum << endl;

    tm.start();
    xcf->setxc();
    tm.stop();
 
    double *exc = xcf->exc;
    double excsum = 0.0;
    for ( int i = 0; i < n3; i++ )
      excsum += rh[0][i] * exc[i];
 
    cout << " Total xc energy: " << excsum * omega / n3 << endl;

    // if not a GGA, compute dExc/da
    if ( !xcf->isGGA() )
    {
      double *vxc1 = xcf->vxc1;
      double dxcsum = 0.0;
      for ( int i = 0; i < n3; i++ )
        dxcsum += rh[0][i] * ( exc[i] - vxc1[i] );
 
      // Note: the energy variation is 3 * dExc/da * delta(a)
      cout << " dExc/da: " << dxcsum * omega / ( n3 * a ) << endl;
    }
    cout << " " << xcf->name() << " time: " << tm.real() << endl;
  } // ixcf
  return 0;
}
