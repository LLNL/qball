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
// use: testMatrix input_file [-check] 
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
using namespace std;

#include <qball/Timer.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <qball/Context.h>
#include <math/Matrix.h>

double aa(int i, int j) { return 1.0/(i+1)+2.0*i/(j+1); }
double bb(int i, int j) { return i-j-3; }

int main(int argc, char **argv)
{
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

  char* infilename = argv[1];
  ifstream infile(infilename);

  assert(argc == 2 || argc == 3);
  bool tcheck = false;
  if ( argc == 3 )
  {
    if ( !strcmp(argv[2],"-check") )
      tcheck = true;
    else
    {
      cerr << " invalid argv[2]" << endl;
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
  }
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
            double aij = aa(i,j);
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

    tm.start();
    c.gemm(ta,tb,1.0,a,b,0.0);
    tm.stop();


    if ( tcheck )
    {
    cout << " checking results..." << endl;
    for ( int m = 0; m < c.nblocks(); m++ )
      for ( int l = 0; l < c.mblocks(); l++ )
        for ( int y = 0; y < c.nbs(m); y++ )  
          for ( int x = 0; x < c.mbs(l); x++ )
          {
            int i = c.i(l,x);
            int j = c.j(m,y);
            double sum = 0.0;
            int kmax = ( ta == 'n' ) ? a.n() : a.m();
            
            if ( ( ta == 'n' ) && ( tb == 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(i,k) * bb(k,j);
            }
            else if ( ( ta != 'n' ) && ( tb == 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(k,i) * bb(k,j);
            }
            else if ( ( ta == 'n' ) && ( tb != 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(i,k) * bb(j,k);
            }
            else if ( ( ta != 'n' ) && ( tb != 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(k,i) * bb(j,k);
            }
            
            int iii = x + l*c.mb();
            int jjj = y + m*c.nb();
            int ival = iii + jjj * c.mloc();
            if ( fabs( c[ival] - sum ) > 1.e-8 )
            {
              cout << " error at element (" << i << "," << j << ") "
                   << c[ival] << " " << sum << endl;
              exit(1);
            }
          }
          
          
       cout << " results checked" << endl;
    }

    cout << " CPU/Real: " << setw(8) << tm.cpu() 
         << " / " << setw(8) << tm.real();
    if ( tm.real() > 0.0 )
    {
      int kmax = ( ta == 'n' ) ? a.n() : a.m();
      cout << "  MFlops: " 
           << (2.0e-6*m_c*n_c*kmax) / tm.real() << endl;
    }
#if 1    
    double norma=a.nrm2();
    if(mype == 0)cout<<"Norm(a)="<<norma<<endl;
    
    double norm;
#if 0
    if(mype == 0)cout<<"DoubleMatrix::matgather..."<<endl;
    double*  aa=new double[a.m()*a.n()];
    a.matgather(aa, a.m());
    if(mype == 0)cout<<"DoubleMatrix::init..."<<endl;
    b.init(aa, a.m());
    norm=b.nrm2();
    if ( mype == 0 ) cout << "Norm(b)=" << norm << endl;
    if ( fabs(norm-norma)>0.000001 )
       cout << "DoubleMatrix: problem with matgather/init" << endl;
#endif

    if ( c.n() == b.m() && c.m() == b.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::transpose..."<<endl;
      tm.reset();
      tm.start();
      c.transpose(1.0,b,0.0);
      tm.stop();
      if ( mype == 0 ) cout << " transpose time: " << tm.real() << endl;
      norm=c.nrm2();
      if(mype == 0)cout<<"Norm(c)="<<norm<<endl;
    }
    
    if(mype == 0)cout<<"DoubleMatrix::scal..."<<endl;
    c.scal(0.5);
    
    if ( a.m() == b.m() && a.n() == b.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::axpy..."<<endl;
      a.axpy(-2., b);
    }
    
    if ( a.m() == c.m() && a.n() == c.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::operator=..."<<endl;
      c=a;
    }
    
    if(mype == 0)cout<<"DoubleMatrix::nrm2..."<<endl;
    norm=c.nrm2();
    if (mype == 0) cout<<"Norm="<<norm<<endl;
        
#if 1
    a.identity();
    DoubleMatrix a2(a);
    a -= a2;
    norm = a.nrm2();
    if (mype == 0) cout << "Norm(a)=" << norm << endl;
#endif
    
    // Eigenvalues and eigenvectors of c if c is square
    if ( c.m() == c.n() && c.mb() == c.nb() )
    {
      for ( int m = 0; m < c.nblocks(); m++ )
        for ( int l = 0; l < c.mblocks(); l++ )
          for ( int y = 0; y < c.nbs(m); y++ )
            for ( int x = 0; x < c.mbs(l); x++ )
            {
              int i = c.i(l,x);
              int j = c.j(m,y);
              int iii = x + l*c.mb();
              int jjj = y + m*c.nb();
              int ival = iii + jjj * c.mloc();
              if ( i == j )
                c[ival] = i + 1.e-5*drand48();
              else
                c[ival] = 1.e-5*drand48();
            }
      tm.reset();
      tm.start();
      if (mype == 0) cout << "Eigenproblem... ";
      DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
      valarray<double> w(c.m());
      c.syev('l',w,z);
      if (mype == 0) cout << " done" << endl;
      tm.stop();
      if (mype == 0) cout << "Eigenproblem time: " << tm.real() << endl;
    }
    
//  Gram-Schmidt orthogonalization of matrix a
    for ( int m = 0; m < a.nblocks(); m++ )
      for ( int l = 0; l < a.mblocks(); l++ )
        for ( int y = 0; y < a.nbs(m); y++ )  
          for ( int x = 0; x < a.mbs(l); x++ )
          {
            int i = a.i(l,x);
            int j = a.j(m,y);
            double aij = aa(i,j);
            int iii = x + l*a.mb();
            int jjj = y + m*a.nb();
            int ival = iii + jjj * a.mloc();
            if ( i == j ) 
              a[ival] = i + 1.e-6*drand48();
            else
              a[ival] = 1.e-6*drand48();
          }
    DoubleMatrix s(a.context(),a.n(),a.n(),a.nb(),a.nb());
    tm.reset();
    tm.start();
    s.syrk('l','t',2.0,a,0.0);
    if (mype == 0) cout << "Gram syrk time: " << tm.real() << endl;
    
    tm.reset();
    tm.start();
    s.syr('l',-1.0,a,0,'r');
    tm.stop();
    if (mype == 0) cout << "Gram syr time: " << tm.real() << endl;
    
    // Cholesky decomposition
    tm.reset();
    tm.start();
    s.potrf('l'); // Cholesky decomposition: S = L * L^T
    tm.stop();
    if (mype == 0) cout << "Gram Cholesky time: " << tm.real() << endl;
    
    // Triangular solve
    tm.reset();
    tm.start();
    // solve triangular system X * L^T = C
    a.trsm('r','l','t','n',1.0,s);
    tm.stop();
    if (mype == 0) cout << "Gram triangular solve time: " << tm.real() << endl;
#endif
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
}
