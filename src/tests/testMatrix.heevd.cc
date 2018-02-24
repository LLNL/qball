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
#include <map>
using namespace std;

#include <qball/Timer.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_APC
#include "apc.h"
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

#if USE_APC
  ApcInit();
#endif

  char* infilename = argv[1];
  ifstream infile(infilename);

  assert(argc == 2);

  Timer tm;
  int nprow, npcol;
  int m_a, n_a, mb_a, nb_a;
  int m_b, n_b, mb_b, nb_b;
  int m_c, n_c, mb_c, nb_c;
  char ta, tb;
  if(mype == 0)
  {
    infile >> nprow >> npcol;
    cout<<"nprow="<<nprow<<", npcol="<<npcol<<endl;
    infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
    cout<<"m_a="<<m_a<<", n_a="<<n_a<<endl;
    infile >> m_b >> n_b >> mb_b >> nb_b >> tb;
    cout<<"m_b="<<m_b<<", n_b="<<n_a<<endl;
    infile >> m_c >> n_c >> mb_c >> nb_c;
    cout<<"m_c="<<m_c<<", n_c="<<n_c<<endl;
  }
#ifdef USE_MPI
  MPI_Bcast(&nprow, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&npcol, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&m_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&n_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&mb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&nb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&m_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&n_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&mb_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&nb_b, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&m_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&n_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&mb_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&nb_c, 1, MPI_INT, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&ta, 1, MPI_CHAR, 0, MPI_COMM_WORLD);    
  MPI_Bcast(&tb, 1, MPI_CHAR, 0, MPI_COMM_WORLD);    
#endif
  {  
    if ( ta == 'N' ) ta = 'n';
    if ( tb == 'N' ) tb = 'n';

    Context ctxt(nprow,npcol);

    if ( mype == 0 )
    {
      cout << " Context " << ctxt.ictxt()
           << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
    }

    DoubleMatrix a(ctxt,m_a,n_a,mb_a,nb_a);
    DoubleMatrix b(ctxt,m_b,n_b,mb_b,nb_b);
    DoubleMatrix c(ctxt,m_c,n_c,mb_c,nb_c);

    if ( mype == 0 )
    {
      cout << " m_a x n_a / mb_a x nb_a / ta = "
           << a.m() << "x" << a.n() << " / "
           << a.mb() << "x" << a.nb() << " / " << ta << endl;
      cout << " m_b x n_b / mb_b x nb_b / tb = "
           << b.m() << "x" << b.n() << " / "
           << b.mb() << "x" << b.nb() << " / " << tb << endl;
      cout << " m_c x n_c / mb_c x nb_c      = "
           << c.m() << "x" << c.n() << " / "
           << c.mb() << "x" << c.nb() << endl;
    }

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
            int ival = iii + jjj * a.mloc();
            a[ival] = aij;
          }

    for ( int m = 0; m < b.nblocks(); m++ )
      for ( int l = 0; l < b.mblocks(); l++ )
        for ( int y = 0; y < b.nbs(m); y++ )  
          for ( int x = 0; x < b.mbs(l); x++ )
          {
            int i = b.i(l,x);
            int j = b.j(m,y);
            // double bij = b.i(l,x) * 10 + b.j(m,y);
            double bij = bb(i,j);
            int iii = x + l*b.mb();
            int jjj = y + m*b.nb();
            int ival = iii + jjj * b.mloc();
            b[ival] = bij;
          }

    ComplexMatrix s(ctxt,m_c,n_c,mb_c,nb_c);
    
    for ( int m = 0; m < s.nblocks(); m++ )
       for ( int l = 0; l < s.mblocks(); l++ )
          for ( int y = 0; y < s.nbs(m); y++ )  
             for ( int x = 0; x < s.mbs(l); x++ ) {
                int i = s.i(l,x);
                int j = s.j(m,y);

                double drand1 =  rand()*nrandinv*maxrand;
                double drand2 =  rand()*nrandinv*maxrand;

                double sij;
                if (i == j)
                   sij = 1.0*i + drand1;
                else
                   sij = drand1;
                
                double zij;
                if (i == j)
                   zij = 1.0*i + drand2;
                else 
                   zij = drand2;
              
                int iii = x + l*s.mb();
                int jjj = y + m*s.nb();
                int ival = iii + jjj * s.mloc();
                complex<double> cij = (sij,zij);
                s[ival] = cij;
             }


    assert(m_c == n_c);
    valarray<double> w(s.m());

    if (mype == 0) 
       cout << "Calculating complex eigenvalues on rectangular context..." << endl;

    ComplexMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());

    tmap["heev-rect"].start();
    //s.heev('l',w);
    //s.heev('l',w,z);
    s.heevd('l',w,z);
    tmap["heev-rect"].stop();

    if (mype == 0) {
       cout << "Eigenvalues:  " << endl;
       for (int i=0; i<s.m(); i++) {
          cout << "  " << w[i];
          if (i%8 == 0) cout << endl;
       }
       cout << endl;
    }

    if (mype == 0) 
       cout << "Done." << endl;

    tmap["total"].stop();
    
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
#if USE_APC
  ApcFinalize();
#endif
}
