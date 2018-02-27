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
// create a square matrix and diagonalize with different Scalapack eigensolvers
//
// usage:  testEigenSolvers nprow npcol matrix_size

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

#include <qball/Context.h>
#include <math/matrix.h>

bool runheevr = false;
bool runheevx = false;
bool runheevd = true;
bool runheev = true;

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.0001;  // maximum random perturbation to identity matrix

int main(int argc, char **argv)
{

  // set up map of timers
  map<string,Timer> tmap;

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

  {

    int nprow, npcol, n;
    double xtol;
    if (argc == 4 || argc == 5) {
      nprow = atoi(argv[1]);
      npcol = atoi(argv[2]);
      n = atoi(argv[3]);
      if (argc == 5)
         xtol = atof(argv[4]);
    }
    else {
      cerr << "Usage:  testEigenSolvers nprow npcol matrix_size [tolerance]" << endl;
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
    
    tmap["total"].start();
    tmap["init"].start();
    Context ctxt(nprow,npcol);

    if ( mype == 0 ) {
      cout << " Context " << ctxt.ictxt()
           << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
    }

    bool diagpe = (ctxt.myrow() == ctxt.mycol());
    
    int mb = n/nprow;
    int nb = n/npcol;
    
    ComplexMatrix s(ctxt,n,n,nb,nb);
    DoubleMatrix r(ctxt,n,n,nb,nb);
    
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
            r[ival] = sij;
          }
    tmap["init"].stop();

    if (runheevr)
    {
       if (mype == 0) 
          cout << "Calculating eigenvalues and eigenvectors with heevr..." << endl;
    
       valarray<double> w1(s.m());
       ComplexMatrix z1(s.context(),s.n(),s.n(),s.nb(),s.nb());
       ComplexMatrix c1(s);

       if (mype == 0)
          cout << "matrix dimensions = " << c1.m() << " x " << c1.n() << " (" << c1.mb() << " x " << c1.nb() << " )" << endl;

       tmap["heevr"].start();
       c1.heevr('l',w1,z1);
       tmap["heevr"].stop();
    
       if (mype == 0) {
          cout << "Eigenvalues:  " << endl;
          for (int i=0; i<s.m(); i++) {
             cout << "  " << w1[i];
             if (i%8 == 0) cout << endl;
          }
          cout << endl;

       }      
    }
    
    if (runheevx)
    {
       if (mype == 0 && argc == 5) 
          cout << "Calculating eigenvalues and eigenvectors with heevx (tolerance = " << xtol << "..." << endl;
       else if (mype == 0 && argc == 4) 
          cout << "Calculating eigenvalues and eigenvectors with heevx..." << endl;
    
       valarray<double> w2(s.m());
       ComplexMatrix z2(s.context(),s.n(),s.n(),s.nb(),s.nb());    
       ComplexMatrix c2(s);

       tmap["heevx"].start();
       if (argc == 5)
          c2.heevx('l',w2,z2,xtol);
       else
          c2.heevx('l',w2,z2);
       tmap["heevx"].stop();
    
       if (mype == 0) {
          cout << "Eigenvalues:  " << endl;
          for (int i=0; i<s.m(); i++) {
             cout << "  " << w2[i];
             if (i%8 == 0) cout << endl;
          }
          cout << endl;
       }      
    }

    if (runheevd)
    {
       if (mype == 0) 
          cout << "Calculating eigenvalues and eigenvectors with heevd..." << endl;
    
       valarray<double> w4(s.m());
       ComplexMatrix z4(s.context(),s.n(),s.n(),s.nb(),s.nb());    
       ComplexMatrix c4(s);

       tmap["heevd"].start();
       c4.heevd('l',w4,z4);
       tmap["heevd"].stop();
       
       //ewd test that eigenvectors are orthonormal
       if (true)
       {
          ComplexMatrix zmult(s.context(),s.n(),s.n(),s.nb(),s.nb());    
          zmult.gemm('c','n',1.0,z4,z4,0.0);
       
          if (zmult.localsize() > 0)
          {
             int nb = s.nb();
             for ( int in = 0; in < nb; in++ ) {
                complex<double>* sp = zmult.valptr(nb*in);
                for ( int im = 0; im < nb; im++ ) {
                   if (diagpe && im == in)
                   {
                      if (abs(real(sp[im])-1.0) > 1.E-6)
                         cout << "DIAG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", zmult val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                   }
                   else
                   {
                      if (abs(real(sp[im])) > 1.E-6 || abs(imag(sp[im]) > 1.E-6))
                         cout << "DIAG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", zmult val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                   }
                }
             }
          }
       }
       
       if (mype == 0) {
          cout << "Eigenvalues:  " << endl;
          for (int i=0; i<s.m(); i++) {
             cout << "  " << w4[i];
             if (i%8 == 0) cout << endl;
          }
          cout << endl;
       }
    }

    if (runheev)
    {
       if (mype == 0) 
          cout << "Calculating eigenvalues and eigenvectors with heev..." << endl;
    
       valarray<double> w3(s.m());
       ComplexMatrix z3(s.context(),s.n(),s.n(),s.nb(),s.nb());    
       ComplexMatrix c3(s);

       tmap["heev"].start();
       c3.heev('l',w3,z3);
       tmap["heev"].stop();

       //ewd test that eigenvectors are orthonormal
       if (true)
       {
          ComplexMatrix zmult(s.context(),s.n(),s.n(),s.nb(),s.nb());    
          zmult.gemm('c','n',1.0,z3,z3,0.0);
          
          if (zmult.localsize() > 0)
          {
             int nb = s.nb();
             for ( int in = 0; in < nb; in++ ) {
                complex<double>* sp = zmult.valptr(nb*in);
                for ( int im = 0; im < nb; im++ ) {
                   if (diagpe && im == in)
                   {
                      if (abs(real(sp[im])-1.0) > 1.E-6)
                         cout << "DIAG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", zmult val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                   }
                   else
                   {
                      if (abs(real(sp[im])) > 1.E-6 || abs(imag(sp[im]) > 1.E-6))
                         cout << "DIAG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", zmult val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                   }
                }
             }
          }
       }
       
       if (mype == 0) {
          cout << "Eigenvalues:  " << endl;
          for (int i=0; i<s.m(); i++) {
             cout << "  " << w3[i];
             if (i%8 == 0) cout << endl;
          }
          cout << endl;
       }
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

}
