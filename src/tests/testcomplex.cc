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
#include <complex>
#include <vector>
#include <iomanip>
#include <iostream>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main() {
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {

    cout << setprecision(10);
  
    // create a complex<double> array
    int zsize = 3;
    vector<complex<double> > z(zsize);
    z[0] = complex<double>(0.111111111,-0.444222444);
    z[1] = complex<double>(3.53,2.42);
    z[2] = complex<double>(0.8877889,9.876543);

    for (int i=0; i<zsize; i++) 
      cout << "z[" << i << "] = " << z[i] << endl;

    double* dp = (double*)&z[0];

    for (int i=0; i<2*zsize; i++) 
      cout << "dp[" << i << "] = " << dp[i] << endl;


    const double dtmp = 1.23456789;
    for (int i=0; i<2*zsize; i++)
      dp[i] += dtmp;
    cout << "Added constant " << dtmp << " to dp" << endl;

    for (int i=0; i<zsize; i++) 
      cout << "z[" << i << "] = " << z[i] << endl;



  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}  
  
  
  
