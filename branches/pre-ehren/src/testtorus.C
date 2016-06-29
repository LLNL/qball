////////////////////////////////////////////////////////////////////////////////
//
// testtorus.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testtorus.C,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

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
