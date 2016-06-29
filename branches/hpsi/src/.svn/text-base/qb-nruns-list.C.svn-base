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
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include <vector>
using namespace std;

#include "qbLink.h"

int main(int argc, char **argv) {

  int mype;
  int npes;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  //if (mype == 0) 
  //  cout << "qb-nruns init, mype = " << mype << ", npes = " << npes << ", argc = " << argc << endl;

  if (argc != 2) {
    cout << "Usage:  qb-nruns-list [file containing list of input files]\n";
    return 1;
  }

  // open list file containing names of input files
  vector<string> inpfiles;
  int ninputs = -1;
  if (mype == 0) {
    ifstream is;
    is.open(argv[1],ios::in);
    string fname;
    if (is.is_open()) {
      while (!is.eof()) {
        getline(is,fname);
        stringstream ss(fname); // Insert the string into a stream
        string stmp;
        while (ss >> stmp)
          inpfiles.push_back(stmp);
      }
    }
    is.close();
    ninputs = inpfiles.size();
  }

  // broadcast list of input files to all other tasks
  MPI_Bcast(&ninputs, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  if (ninputs <= 0) {
    if (mype == 0)
      cout << "<ERROR> Couldn't read input file list from " << argv[1] << " </ERROR>" << endl;
    return 1;
  }
  inpfiles.resize(ninputs);
  for (int i=0; i<ninputs; i++) {
    int strlen;
    if (mype == 0)
      strlen = inpfiles[i].length();
    MPI_Bcast(&strlen, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    char* buf = new char[strlen+1];
    if (mype == 0) {
      inpfiles[i].copy(buf,string::npos);
      buf[strlen]=0;
      assert(buf[strlen]=='\0');
    }
    MPI_Bcast(buf,strlen+1,MPI_CHAR,0,MPI_COMM_WORLD);
    inpfiles[i] = buf;
    delete [] buf;
  }

  int npes_sub = npes/ninputs;
  assert(npes_sub > 0);
  if (npes_sub%2 !=0 && mype == 0) 
    cout << "<WARNING> Odd number of pes for each run:  npes_sub = " << npes_sub << " </WARNING>" << endl;
  vector<qbLink*> qb(ninputs);
  vector<string> outfiles(ninputs);
  int firstpe = 0;
  for (int i=0; i<ninputs; i++) {
    int lastpe = firstpe + npes_sub-1;
    outfiles[i] = inpfiles[i] + ".out";
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






