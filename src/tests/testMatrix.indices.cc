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
// use: testMatrix input_file [-check] [-ortho]
// input_file:
// nprow npcol
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
using namespace std;

#include <qball/Timer.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <qball/Context.h>
#include <math/Matrix.h>

// mpi_trace functions
//extern "C" {
//  extern  void    trace_start();
//  extern  void    trace_stop();
//}
//

double aa(int i, int j) { return 1.0/(i+1)+2.0/(j+1); }
double bb(int i, int j) { return i-j-3; }

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.01;  // maximum random perturbation to identity matrix

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

  int nprow, npcol;
  
  nprow = 4;
  npcol = npes/nprow;

  int m,n,mb,nb;

  // square matrix
  n = 45;
  nb = n/npcol;
  if (n%npcol != 0) nb++;

  m = n;
  mb = nb;
  
  
  Context ctxt(nprow,npcol);
  DoubleMatrix a(ctxt,m,n,mb,nb);

  // fill matrix
  int mloc = a.mloc();
  for ( int m = 0; m < a.nblocks(); m++ )
    for ( int l = 0; l < a.mblocks(); l++ )
      for ( int y = 0; y < a.nbs(m); y++ )  
        for ( int x = 0; x < a.mbs(l); x++ )
        {
          int i = a.i(l,x);
          int j = a.j(m,y);
          // double aij = a.i(l,x) * 10 + a.j(m,y);
          //double aij = aa(i,j);

          double drand =  rand()*nrandinv*maxrand;
          double aij;
          if (i == j)
            aij = 1.0 + drand;
          else
            aij = drand;

          int iii = x + l*a.mb();
          int jjj = y + m*a.nb();
          int ival = iii + jjj * mloc;

          //ewd DEBUG
          //cout << "mype = " << mype << ", i = " << i << ", j = " << j << ", aij = " << aij << endl;

          a[ival] = aij;
        }

  // print out basic matrix information
  if (mype == 0)
    cout << "Process grid:  " << npes << " pes, " << nprow << " x " << npcol << endl;
  cout << "matrix info, mype = " << mype << ", m x n = " << a.m() << " x " << a.n() << ", mb x nb = " << a.mb() << " x " << a.nb() << ", size = " << a.size() << ", myrow = " << ctxt.myrow() << ", mycol = " << ctxt.mycol() << endl;
  
  // diagonal values for local columns
  int nloc = a.nloc();
  vector<double> locdiag(nloc);
  for (int i=0; i<nloc; i++)
    locdiag[i] = 0.0;

  const int myrow = ctxt.myrow();
  const int mycol = ctxt.mycol();
  for (int i=0; i<n; i++) {
    int iprow = a.pr(i);
    int ipcol = a.pc(i);
    if (ipcol == mycol) {
      if (iprow == myrow) {
        int jloc = a.y(i);
        int index = a.x(i) + mloc*jloc;
        //cout << "mype = " << mype << ", VAL " << i << " = " << a[index] << ", y = " << a.y(i) << endl;
        locdiag[jloc] = a[index];
      }
    }
  }

  // sum locdiag across process rows so that all tasks in a given process column know
  // the diagonal matrix values for their local columns
  ctxt.dsum('c',nloc,1,&locdiag[0],nloc);

  // divide all data by the diagonal
  double* ap = a.valptr();
  for (int i=0; i<nloc; i++) {
    const double ival = locdiag[i];
    if (ival != 0.0) {
      for (int j=0; j<mloc; j++) {
        double old = ap[i*mloc+j];
        ap[i*mloc+j] /= ival;
        cout << "mype = " << mype << ", i = " << i << ", j = " << j << ", old value = " << old << ", new val = " << ap[i*mloc+j] << ", diag for column i = " << ival << endl;
      }
    }
  }
}
