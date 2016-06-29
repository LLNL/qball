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
// testMPIsubcomm.c
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

#ifdef USE_MPI  
#include <mpi.h>
#endif

int main(int argc, char **argv) {
  int mype;
  int npes;
#ifdef USE_MPI  

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  int nr = atoi(argv[1]);
  int nc = atoi(argv[2]);
  
  vector<int> pmap1;
  MPI_Comm comm1;

  int psize1 = nr*nc;
  pmap1.resize(psize1);
  // column-major order
  int i = 0;
  for ( int ic = 0; ic < nc; ic++ )
    for ( int ir = 0; ir < nr; ir++ ) {
      pmap1[ir+nr*ic] = i;
      i++;
    }
  
  if (mype == 0)
    cout << "Creating main column-major communicator...  psize1 = " << psize1 << endl;

  MPI_Group group_world, subgroup;
  MPI_Comm_group(MPI_COMM_WORLD,&group_world);
  MPI_Group_incl(group_world,psize1,&pmap1[0],&subgroup);
  MPI_Comm_create(MPI_COMM_WORLD,subgroup,&comm1);
  MPI_Group_free(&group_world);
  MPI_Group_free(&subgroup);

  vector<int> pmap2;
  MPI_Comm comm2;

  int psize2 = 0.5*nr*nc;
  pmap2.resize(psize2);
  for (int j=0; j<psize2; j++)
    pmap2[j] = j+psize2;

  if (mype == 256)
    cout << "Creating subcommunicator of pes 256-511, psize2 = " << psize2 << endl;

  MPI_Group group_comm1, subgroup2;
  MPI_Comm_group(comm1,&group_comm1);
  MPI_Group_incl(group_comm1,psize2,&pmap2[0],&subgroup2);
  MPI_Comm_create(comm1,subgroup2,&comm2);
  MPI_Group_free(&group_comm1);
  MPI_Group_free(&subgroup2);

  vector<int> pmap3;
  MPI_Comm comm3;

  comm3 = 0;

  if (mype >= 256) {

    int psize3 = 0.5*psize2;
    pmap3.resize(psize3);
    for (int j=0; j<psize3; j++)
      pmap3[j] = j+psize3;
    
    if (mype == 384)
      cout << "Creating subcommunicator of pes 384-511, psize3 = " << psize3 << endl;

    MPI_Group group_comm2, subgroup3;
    if (mype == 384)
      cout << "MPI_Comm_group call" << endl;
    MPI_Comm_group(comm2,&group_comm2);
    if (mype == 384)
      cout << "MPI_Group_incl call" << endl;
    MPI_Group_incl(group_comm2,psize3,&pmap3[0],&subgroup3);
    if (mype == 384)
    cout << "MPI_Comm_create call" << endl;
    MPI_Comm_create(comm2,subgroup3,&comm3);
    if (mype == 384)
      cout << "MPI_Group_free call" << endl;
    MPI_Group_free(&group_comm2);
    MPI_Group_free(&subgroup3);
    
  }

  MPI_Finalize();

#endif
}
