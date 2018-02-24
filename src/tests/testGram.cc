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

//
// create a square matrix and orthogonalize with Gram-Schmidt
//
// usage:  testGram nprow npcol matrix_size

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <map>
#include <qball/Timer.h>
#include <qball/Context.h>
#include <math/Matrix.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif
using namespace std;

int main(int argc, char **argv)
{
   const bool testOrthog = false;  // explicitly check orthogonalization
   const bool copyToSquareContext = true;  // copy data in place to square context for potrf
   
   // set up map of timers
   map<string,Timer> tmap;

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

      int nprow, npcol, m, n;
      if (argc == 4) {
         nprow = atoi(argv[1]);
         npcol = npes/nprow;
         m = atoi(argv[2]);
         n = atoi(argv[3]);
      }
      else {
         cerr << "Usage:  testGram nprow m n" << endl;
#if USE_MPI
         MPI_Abort(MPI_COMM_WORLD,2);
#else
         exit(2);
#endif
      }
    
      tmap["total"].start();
      tmap["init"].start();
      Context ctxt(nprow,npcol);
      Context ctxtsq(ctxt,npcol,npcol,0,0);

      if ( mype == 0 ) {
         cout << " Context " << ctxt.ictxt()
              << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
      }
    
      int mb = m/nprow + (m%nprow > 0 ? 1 : 0);
      int nb = n/npcol + (n%npcol > 0 ? 1 : 0);
    
      ComplexMatrix c(ctxt,m,n,mb,nb);
      ComplexMatrix s(ctxt,n,n,nb,nb);

      // randomize initial values
      {
         srand48(ctxt.myproc());
         for ( int in = 0; in < nb; in++ ) {
            complex<double>* p = c.valptr(mb*in);
            for ( int im = 0; im < mb; im++ ) {
               double dre = drand48();
               double dim = drand48();
               p[im] = 0.02 * complex<double>(dre,dim);
            }
         }
      }
      tmap["init"].stop();

#ifdef HPM
      HPM_Start("herk");
#endif
      tmap["herk"].start();
      s.herk('l','c',1.0,c,0.0);
      tmap["herk"].stop();
#ifdef HPM
      HPM_Stop("herk");
#endif

#ifdef HPM
      HPM_Start("potrf");
#endif
      tmap["potrf"].start();
      if (copyToSquareContext)
      {
         ComplexMatrix ssq(ctxtsq,n,n,nb,nb);
         if (ssq.active())
         {
            s.copyInPlace(ssq);
            ssq.potrf('l'); // Cholesky decomposition: S = L * L^T
            ssq.copyInPlace(s);
         }
      }
      else
      {
         s.potrf('l'); // Cholesky decomposition: S = L * L^T
      }
      tmap["potrf"].stop();
#ifdef HPM
      HPM_Stop("potrf");
#endif

#ifdef HPM
      HPM_Start("trsm");
#endif
      tmap["trsm"].start();
      c.trsm('r','l','c','n',1.0,s);
      tmap["trsm"].stop();
#ifdef HPM
      HPM_Stop("trsm");
#endif

    
      if (mype == 0) 
         cout << "Done." << endl;
      tmap["total"].stop();

      //ewd DEBUG:  check that orthogonalization was successful
      if (testOrthog)
      {
         if (mype == 0) 
            cout << "Checking orthogonality by computing overlap matrix:  s = " << n << " x " << n << ", nb = " << nb << endl;
         s.gemm('c','n',1.0,c,c,0.0);

         cout << "OVERLAP, mype = " << mype << ", localsize = " << s.localsize() << endl;
       
         if (s.localsize() > 0)
         {
            for ( int in = 0; in < nb; in++ ) {
               complex<double>* sp = s.valptr(nb*in);
               for ( int im = 0; im < nb; im++ ) {
                  if (im == in)
                  {
                     if (abs(real(sp[im])-1.0) > 1.E-6)
                        cout << "ORTHOG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", s val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                  }
                  else
                  {
                     if (abs(real(sp[im])) > 1.E-6 || abs(imag(sp[im]) > 1.E-6))
                        cout << "ORTHOG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", s val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                  }
               }
            }
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      //ewd DEBUG    

    
      for ( map<string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ )
      {
         double time = (*i).second.cpu();
         double tmin = time;
         double tmax = time;
    
         ctxt.dmin(1,1,&tmin,1);
         ctxt.dmax(1,1,&tmax,1);
         if (mype == 0)
         { 
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
