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
// testtorus.C
//
////////////////////////////////////////////////////////////////////////////////

#include "mpi.h"
#include <iostream>
#include <rts.h>
using namespace std;

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);

  int n = 0, myrank, mysize, rank;
  int x, y, z, t;
  BGLPersonality personality;
  rts_get_personality (&personality, sizeof(personality));

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);
  // PMI_rank2torus(myrank, &x, &y, &z, &t);
  x = personality.xCoord;
  y = personality.yCoord;
  z = personality.zCoord;

  cout << myrank << ": at " << x << " " << y << " " << z << endl;

  MPI_Finalize();
}
