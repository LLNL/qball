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
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <omp.h>
#include <qball/blas.h>
#include <qball/Timer.h>
using namespace std;

#ifdef HAVE_BGQLIBS
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.000001;  // maximum random perturbation to identity matrix

int main(int argc, char **argv)
{
   Timer tm;

   // choose random number seed based on system clock
   srand((unsigned)time(NULL));

   int mype;
   int npes;
#ifdef USE_MPI
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#else
   npes=1;
   mype=0;
#endif

   const int vsize = 10000;
   const int ione = 1;
   vector<complex<double> > avec(vsize);
   vector<complex<double> > bvec(vsize);

   for (int ii=0; ii<vsize; ii++) {
      double drand1 =  rand()*nrandinv*maxrand;
      double drand2 =  rand()*nrandinv*maxrand;
      avec[ii] = complex<double>(drand1,drand2);
   }
   for (int ii=0; ii<vsize; ii++) {
      double drand1 =  rand()*nrandinv*maxrand;
      double drand2 =  rand()*nrandinv*maxrand;
      bvec[ii] = complex<double>(drand1,drand2);
   }

   //test
   complex<double> dotcheck = complex<double>(0.0,0.0);
   for (int ii=0; ii<vsize; ii++)
      dotcheck += conj(avec[ii])*bvec[ii];
   cout << "dotcheck = " << dotcheck << endl;

   tm.start();
#ifdef HPM
   HPM_Start("zvec");
#endif

   complex<double> qv = FC_FUNC(zdotc, ZDOTC)((int*)&vsize,&avec[0],(int*)&ione,&bvec[0],(int*)&ione);

#ifdef HPM
   HPM_Stop("zvec");
#endif
   tm.stop();

   //int nthreads = omp_get_max_threads();
   if (mype == 0)
      cout << "size = " << vsize << ", qv = " << real(qv) << ", zdotc time = " << setprecision(5) << setw(8) << tm.real() << endl;
 
#ifdef USE_MPI
   MPI_Finalize();
#endif
}
