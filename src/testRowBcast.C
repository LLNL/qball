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
// testRowBcast.C
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <mpi.h>
#include "Context.h"
using namespace std;

int main(int argc, char **argv)
{
  int mype;
  int npes;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  const int nr = atoi(argv[1]);
  const int nc = npes/nr;
  const int datasize = atoi(argv[2]);
  
  if (mype == 0)
     cout << "Creating main " << nr << " x " << nc << " context." << endl;
  Context* spincontext = new Context(nr,nc);

  Context* my_row_ctxt = 0;
  for ( int irow = 0; irow < spincontext->nprow(); irow++ ) {
     Context* row_ctxt = new Context(*spincontext,1,spincontext->npcol(),irow,0);
     spincontext->barrier();
     if ( irow == spincontext->myrow() )
        my_row_ctxt = row_ctxt;
     else
        delete row_ctxt;
  }

  // fill source vector w. random data
  vector<double> buffer(datasize,0.0);
  vector<double> source(datasize);
  srand48(mype);
  for (int ii=0; ii<source.size(); ii++)
     source[ii] = drand48();
  
  // broadcast data from each column across rows using Context::dbcast, measure time
  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  {
     const int myrow = spincontext->myrow();
     for (int ipcol = 0; ipcol < spincontext->npcol(); ipcol++)
     {
        if (spincontext->mycol() == ipcol)
           spincontext->dbcast_send('r',' ',datasize,1,&source[0],datasize);
        else
           spincontext->dbcast_recv('r',' ',datasize,1,&buffer[0],datasize,myrow,ipcol);
     }
  }
  double t2 = MPI_Wtime();

  {
     double bcastTime = t2-t1;
     double tmin = bcastTime;
     double tmax = bcastTime;
     spincontext->dmin(1,1,&tmin,1);
     spincontext->dmax(1,1,&tmax,1);
     if (mype == 0)
        cout << "ctxt.dbcast timing:  " << spincontext->npcol() << " row bcasts sent/received, " << datasize*sizeof(double)/1024 << " kB each, min = " << tmin << ", max = " << tmax << endl;
  }
  
  // now do the same thing with MPI_Bcast
  MPI_Barrier(MPI_COMM_WORLD);
  double t3 = MPI_Wtime();
  {
     const int myrow = spincontext->myrow();
     for (int ipcol = 0; ipcol < spincontext->npcol(); ipcol++)
     {
        int bcastTask = spincontext->nprow()*ipcol + myrow;
        if (mype == bcastTask)
           buffer = source;
        MPI_Bcast(&buffer[0], datasize, MPI_DOUBLE, 0, my_row_ctxt->comm());
     }
  }
  double t4 = MPI_Wtime();

  {
     double bcastTimeMPI = t4-t3;
     double tmin = bcastTimeMPI;
     double tmax = bcastTimeMPI;
     spincontext->dmin(1,1,&tmin,1);
     spincontext->dmax(1,1,&tmax,1);
     if (mype == 0)
        cout << "MPI_Bcast timing:  " << spincontext->npcol() << " row bcasts sent/received, " << datasize*sizeof(double)/1024 << " kB each, min = " << tmin << ", max = " << tmax << endl;
  }  
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
}
