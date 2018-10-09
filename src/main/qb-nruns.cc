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
#include <fstream>
#include <mpi.h>
#include <vector>
using namespace std;

#include <qball/qbLink.h>

int main(int argc, char **argv) {

  int mype;
  int npes;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  //if (mype == 0) 
  //  cout << "qb-nruns init, mype = " << mype << ", npes = " << npes << ", argc = " << argc << endl;

  int ninputs = argc-1;
  int npes_sub = npes/ninputs;
  assert(npes_sub > 0);
  if (npes_sub%2 !=0 && mype == 0) 
    cout << "<WARNING> Odd number of pes for each run:  npes_sub = " << npes_sub << " </WARNING>" << endl;
  vector<qbLink*> qb(ninputs);
  vector<string> inpfiles(ninputs);
  vector<string> outfiles(ninputs);
  int firstpe = 0;
  for (int i=0; i<ninputs; i++) {
    int lastpe = firstpe + npes_sub-1;
    string inp(argv[i+1]);
    inpfiles[i] = inp;
    outfiles[i] = inp + ".out";
    qbLink* tmpqb = new qbLink(outfiles[i],firstpe,lastpe);
    qb[i] = tmpqb;
    firstpe += npes_sub;
  }

  for (int i=0; i<ninputs; i++) 
    if (qb[i]->active())
      qb[i]->processInputFile(inpfiles[i]);

  for (int i=0; i<ninputs; i++) 
    if (qb[i] != 0)
      delete qb[i];

  MPI_Finalize();
}






