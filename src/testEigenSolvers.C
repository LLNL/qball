
// $Id: testEigenSolvers.C,v 1.1 2007/07/24 16:40:14 draeger1 Exp $
//
// create a square matrix and diagonalize with different Scalapack eigensolvers
//
// usage:  testEigenSolvers nprow npcol matrix_size

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

#include "Timer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Context.h"
#include "Matrix.h"

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

  {

    int nprow, npcol, n;
    if (argc == 4) {
      nprow = atoi(argv[1]);
      npcol = atoi(argv[2]);
      n = atoi(argv[3]);
    }
    else {
      cerr << "Usage:  testEigenSolvers nprow npcol matrix_size" << endl;
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
    
    Timer tm;
    Context ctxt(nprow,npcol);

    if ( mype == 0 ) {
      cout << " Context " << ctxt.ictxt()
           << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
    }
    
    int mb = n/nprow;
    int nb = n/npcol;
    
    ComplexMatrix s(ctxt,n,n,nb,nb);
    
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
    /*
    if (mype == 0) {
      cout << "Starting matrix:" << endl;
      for ( int m = 0; m < s.nblocks(); m++ )
        for ( int l = 0; l < s.mblocks(); l++ )
          for ( int y = 0; y < s.nbs(m); y++ )  
            for ( int x = 0; x < s.mbs(l); x++ ) {
              int i = s.i(l,x);
              int j = s.j(m,y);
            
              int iii = x + l*s.mb();
              int jjj = y + m*s.nb();
              int ival = iii + jjj * s.mloc();
              cout << "  " << iii << " " << jjj << "   " << s[ival] << endl;
            }
      cout << endl;      
    }
    */

    valarray<double> w1(s.m());
    valarray<double> w2(s.m());

    if (mype == 0) 
      cout << "Calculating eigenvalues and eigenvectors with heevr..." << endl;
    
    ComplexMatrix z1(s.context(),s.n(),s.n(),s.nb(),s.nb());
    ComplexMatrix z2(s.context(),s.n(),s.n(),s.nb(),s.nb());
    
    ComplexMatrix c1(s);

    if (mype == 0)
      cout << "matrix dimensions = " << c1.m() << " x " << c1.n() << " (" << c1.mb() << " x " << c1.nb() << " )" << endl;
    
    c1.heevr('l',w1,z1);
    
    if (mype == 0) {
      cout << "Eigenvalues:  " << endl;
      for (int i=0; i<s.m(); i++) {
        cout << "  " << w1[i];
        if (i%8 == 0) cout << endl;
      }
      cout << endl;

      /*
      cout << "Eigenvectors:  " << endl;
      for ( int m = 0; m < z1.nblocks(); m++ )
        for ( int l = 0; l < z1.mblocks(); l++ )
          for ( int y = 0; y < z1.nbs(m); y++ )  
            for ( int x = 0; x < z1.mbs(l); x++ ) {
              int i = z1.i(l,x);
              int j = z1.j(m,y);
            
              int iii = x + l*z1.mb();
              int jjj = y + m*z1.nb();
              int ival = iii + jjj * z1.mloc();
              cout << "  " << iii << " " << jjj << "   " << z1[ival] << endl;
            }
      cout << endl;
      */
      
    }
    
    if (mype == 0) 
      cout << "Calculating eigenvalues and eigenvectors with heevd..." << endl;
    
    ComplexMatrix c2(s);
    c2.heevd('l',w2,z2);
    
    if (mype == 0) {
      cout << "Eigenvalues:  " << endl;
      for (int i=0; i<s.m(); i++) {
        cout << "  " << w2[i];
        if (i%8 == 0) cout << endl;
      }
      cout << endl;

      /*
      cout << "Eigenvectors:  " << endl;
      for ( int m = 0; m < z2.nblocks(); m++ )
        for ( int l = 0; l < z2.mblocks(); l++ )
          for ( int y = 0; y < z2.nbs(m); y++ )  
            for ( int x = 0; x < z2.mbs(l); x++ ) {
              int i = z2.i(l,x);
              int j = z2.j(m,y);
            
              int iii = x + l*z2.mb();
              int jjj = y + m*z2.nb();
              int ival = iii + jjj * z2.mloc();
              cout << "  " << iii << " " << jjj << "   " << z2[ival] << endl;
            }
      cout << endl;
      */

      
    }

    /*
    if (mype == 0) 
      cout << "Calculating eigenvalues and eigenvectors with heev..." << endl;
    
    ComplexMatrix c3(s);
    c3.heev('l',w,z);
    
    if (mype == 0) {
      cout << "Eigenvalues:  " << endl;
      for (int i=0; i<s.m(); i++) {
        cout << "  " << w[i];
        if (i%8 == 0) cout << endl;
      }
      cout << endl;
      cout << "Eigenvectors:  " << endl;
      //z.print(cout);
      cout << endl;
    }
    */
    
    if (mype == 0) 
      cout << "Done." << endl;
    tmap["total"].stop();
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif

}
