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
// usage:  testPgemm nprow npcol m n

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
#include "Timer.h"
#include "Context.h"
#include "Matrix.h"
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OLD_CTF
#include "cyclopstf.h"
#endif
#ifdef USE_CTF
#include "cyclopstf.hpp"
#endif
#ifdef HAVE_BGQLIBS
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#endif
using namespace std;

int main(int argc, char **argv)
{

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

  /*  
#ifdef USE_CTF
  CTF* myctf_ = new CTF;

  if (mype == 0)
     cout << "Calling myctf_->init..." << endl;

#ifdef HAVE_BGQLIBS
  myctf_->init(MPI_COMM_WORLD,mype,npes,,MACHINE_BGQ);
#else
  myctf_->init(MPI_COMM_WORLD,mype,npes);
#endif
#endif
  */
  
#ifdef USE_OLD_CTF
  {
    int myRank,numPes;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numPes);

    if (mype == 0)
       cout << "Calling CTF_init..." << endl;

    CTF_init(MPI_COMM_WORLD, MACHINE_BGQ, myRank, numPes); 

    if (mype == 0)
       cout << "Calling CTF_init_complex..." << endl;

    CTF_init_complex(MPI_COMM_WORLD, MACHINE_BGQ, myRank, numPes); 
    //CTF_init(MPI_COMM_WORLD, myRank, numPes); 
    //CTF_init_complex(MPI_COMM_WORLD, myRank, numPes); 
  }    
#endif

  {
     int nprow, npcol, m, n, mb, nb;
     if (argc == 6) {
        nprow = atoi(argv[1]);
        //npcol = atoi(argv[2]);
        npcol = npes/nprow;
        m = atoi(argv[2]);
        n = atoi(argv[3]);
        mb = atoi(argv[4]);
        nb = atoi(argv[5]);
     }
     else {
        cerr << "Usage:  testPgemmBlock nprow npcol m n mb nb" << endl;
#if USE_MPI
        MPI_Abort(MPI_COMM_WORLD,2);
#else
        exit(2);
#endif
     }

#ifdef HAVE_BGQLIBS
     Personality_t pers;
     Kernel_GetPersonality(&pers, sizeof(pers));
     const int nTDim = 5;
     vector<int> torusdim(nTDim);
     torusdim[0] = pers.Network_Config.Anodes;
     torusdim[1] = pers.Network_Config.Bnodes;
     torusdim[2] = pers.Network_Config.Cnodes;
     torusdim[3] = pers.Network_Config.Dnodes;
     torusdim[4] = pers.Network_Config.Enodes;
     int nNodes = torusdim[0]*torusdim[1]*torusdim[2]*torusdim[3]*torusdim[4];
     int tasksPerNode = npes / nNodes;
     
     bool torusMult = false;
     int ncol = npes / nprow;
     for (int ii=0; ii<nTDim; ii++)
        if (ncol == torusdim[ii]*tasksPerNode)
           torusMult = true;
     for (int ii=0; ii<nTDim; ii++)
        for (int jj=ii+1; jj<nTDim; jj++)
           if (ncol == torusdim[ii]*torusdim[jj]*tasksPerNode)
              torusMult = true;

     if ( mype == 0 )
     {
        if (torusMult)
           cout << "nrowmax = " << nprow << ", ncol = " << ncol << " is compatible with BG/Q torus " <<
               torusdim[0] << " x " << torusdim[1] << " x " << torusdim[2] << " x " << torusdim[3] << " x " <<
               torusdim[4] << ", tasksPerNode = " << tasksPerNode << endl;
        else
           cout << "<WARNING> nrowmax = " << nprow << ", ncol = " << ncol << " is NOT compatible with BG/Q torus! " <<
               torusdim[0] << " x " << torusdim[1] << " x " << torusdim[2] << " x " << torusdim[3] << " x " <<
               torusdim[4] << ", tasksPerNode = " << tasksPerNode << " </WARNING> " << endl;
     }        
#endif



    if (mype == 0)
       cout << "Initializing matrices..." << endl;

     
     tmap["total"].start();
     tmap["init"].start();
     Context ctxt(nprow,npcol);

     if ( mype == 0 ) {
        cout << " Context " << ctxt.ictxt()
             << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
     }
    
     //int mb = m/nprow + (m%nprow > 0 ? 1 : 0);
     //int nb = n/npcol + (n%npcol > 0 ? 1 : 0);

     if ( mype == 0 )
        cout << " c dimensions " << ctxt.ictxt()
             << ": " << m << " x " << n << " (" << mb << " x " << nb << ")" << endl;

     
     ComplexMatrix c1(ctxt,m,n,mb,nb);
     ComplexMatrix c2(ctxt,m,n,mb,nb);
     ComplexMatrix s1(ctxt,n,n,nb,nb);
     ComplexMatrix s2(ctxt,n,n,nb,nb);

     // randomize initial values:  ewd, this is core dumping on some sizes, make sure local ranges are right
     
     {
        int mloc = c1.mloc();
        int nloc = c1.nloc();

        srand48(ctxt.myproc());
        for ( int in = 0; in < nloc; in++ ) {
           complex<double>* p1 = c1.valptr(mloc*in);
           complex<double>* p2 = c2.valptr(mloc*in);
           for ( int im = 0; im < mloc; im++ ) {
              double dre1 = drand48();
              double dim1 = drand48();
              p1[im] = 0.02 * complex<double>(dre1,dim1);
              double dre2 = drand48();
              double dim2 = drand48();
              p2[im] = 0.02 * complex<double>(dre2,dim2);
           }
        }
     }
     tmap["init"].stop();

     if ( mype == 0 )
        cout << "Initialization complete, calling pzgemm..." << endl;

     MPI_Barrier(MPI_COMM_WORLD);

     tmap["pzgemm-psda1"].start();
     s1.gemm('c','n',1.0,c1,c2,0.0);
     tmap["pzgemm-psda1"].stop();

     MPI_Barrier(MPI_COMM_WORLD);
     
     tmap["pzgemm-psda2"].start();
     c1.gemm('n','n',-1.0,c2,s1,1.0);
     tmap["pzgemm-psda2"].stop();
    
    if (mype == 0) 
      cout << "Done." << endl;
    tmap["total"].stop();


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

#ifdef USE_OLD_CTF
  CTF_exit();
#endif  
#ifdef USE_MPI
  MPI_Finalize();
#endif

}
