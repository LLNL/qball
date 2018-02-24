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
// test Matrix
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testMatrix input_file
// input_file:
// m_a n_a mb_a nb_a transa
// m_b n_b mb_b nb_b transb
// m_c n_c mb_c nb_c
//

#include <config.h>

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <vector>
#include <map>
#include <math/blas.h>
#include <omp.h>
#include <qball/Timer.h>

#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif

using namespace std;


#ifdef USE_MPI
#include <mpi.h>
#endif

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.000001;  // maximum random perturbation to identity matrix

int main(int argc, char **argv)
{

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

   assert(argc == 5);
   Timer tm;

   int mm = atoi(argv[1]);
   int nn = atoi(argv[2]);
   int kk = atoi(argv[3]);
   int transpose = atoi(argv[4]);  // 1 = transpose, 0 = don't


#ifdef USE_MPI
   //MPI_Bcast(&mm, 1, MPI_INT, 0, MPI_COMM_WORLD);    
   //MPI_Bcast(&nn, 1, MPI_INT, 0, MPI_COMM_WORLD);    
   //MPI_Bcast(&kk, 1, MPI_INT, 0, MPI_COMM_WORLD);    
#endif

   double dzero = 0.0;
   double done = 1.0;
   char cc='t';
   if (!transpose)
      cc = 'n';
   char cn='n';

   vector<double > avec(mm*kk);
   vector<double > bvec(nn*kk);
   vector<double > cvec(mm*nn);
   for (int ii=0; ii<avec.size(); ii++) {
      double drand1 =  rand()*nrandinv*maxrand;
      avec[ii] = drand1;
   }
   for (int ii=0; ii<bvec.size(); ii++) {
      double drand1 =  rand()*nrandinv*maxrand;
      bvec[ii] = drand1;
   }
   for (int ii=0; ii<cvec.size(); ii++) {
      double drand1 =  rand()*nrandinv*maxrand;
      cvec[ii] = drand1;
   }

   const int niter = 100;
   tm.start();
#ifdef HPM
   HPM_Start("dgemm1");
#endif
   for (int iter=0; iter<niter; iter++)
      dgemm(&cc,&cn,&mm,&nn,&kk,&done,&avec[0],&kk,&bvec[0],&kk,&dzero,&cvec[0],&mm);
#ifdef HPM
   HPM_Stop("dgemm1");
#endif
   tm.stop();

   int nthreads = omp_get_max_threads();
   if (mype == 0)
      cout << "M = " << mm << " N = " << nn << " K = " << kk << ", transpose = " << transpose << ", dgemm time = " << setprecision(5) << setw(8) << tm.real()<< " sec, GFlops = " << npes*niter*(2.0e-9*mm*nn*kk) / tm.real() << " on " << npes << " pes, " << nthreads << " threads, niter = " << niter << endl;
 
#ifdef USE_MPI
   MPI_Finalize();
#endif
}
