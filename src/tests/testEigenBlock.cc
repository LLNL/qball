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
// timing for eigensolve
//
// usage:  testEigenBlock nprow npcol m n mb nb

#include <config.h>

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <map>
using namespace std;

#include "Timer.h"

#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Context.h"
#include "Matrix.h"

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.0001;  // maximum random perturbation to identity matrix

int main(int argc, char **argv)
{

   // set up map of timers
   map<string,Timer> tmap;
   tmap["total"].start();

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

   {

      int nprow, npcol, m, n, mb, nb;
      if (argc == 7) {
         nprow = atoi(argv[1]);
         npcol = atoi(argv[2]);
         m = atoi(argv[3]);
         n = atoi(argv[4]);
         mb = atoi(argv[5]);
         nb = atoi(argv[6]);
      }
      else {
         cerr << "Usage:  testEigenBlock nprow npcol m n mb nb" << endl;
#if USE_MPI
         MPI_Abort(MPI_COMM_WORLD,2);
#else
         exit(2);
#endif
      }
      
      Timer tm;
      Context ctxt(nprow,npcol);

      if ( mype == 0 ) {
         cout << " Context " << ctxt.ictxt()
              << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
         cout << " h dimensions " << ctxt.ictxt()
              << ": " << n << " x " << n << " (" << nb << " x " << nb << ")" << endl;
      }
    
      ComplexMatrix h(ctxt,n,n,nb,nb);

      srand48(ctxt.myproc());
      for ( int n = 0; n < h.nloc(); n++ ) {
         complex<double>* p = h.valptr(h.mloc()*n);
         for ( int i = 0; i < h.mloc(); i++ ) {
            double re = drand48();
            double im = drand48();
            p[i] = 0.02 * complex<double>(re,im);
         }
      }
      
      ComplexMatrix s(ctxt,h.n(),h.n(),h.nb(),h.nb());

#ifdef HPM  
      HPM_Start("eigen");
#endif
      valarray<double> w(s.m());
      tmap["heevd"].start();
      h.heevd('l',w,s);
      tmap["heevd"].stop();
      
#ifdef HPM  
      HPM_Stop("eigen");
#endif

      tmap["total"].stop();

      //ewd DEBUG:  print eigenvalues
      if (true)
      {
         if (mype == 0) {
            cout << "Eigenvalues:  " << endl;
            for (int i=0; i<s.m(); i++) {
               cout << "  " << w[i];
               if (i%8 == 0) cout << endl;
            }
            cout << endl;
         }
      }
   
      MPI_Barrier(MPI_COMM_WORLD);
      //ewd DEBUG    

      for ( map<string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
         double time = (*i).second.cpu();
         double tmin = time;
         double tmax = time;
         
         ctxt.dmin(1,1,&tmin,1);
         ctxt.dmax(1,1,&tmax,1);
         if (mype == 0) { 
            
            cout << "  timing "
                 << setw(10) << (*i).first
                 << " : " << setprecision(4) << setw(9) << tmin
                 << " "   << setprecision(4) << setw(9) << tmax << endl;
         }
      }
   }
#ifdef USE_MPI
  MPI_Finalize();
#endif

}
