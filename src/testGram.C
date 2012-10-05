
// $Id: testGram.C,v 1.1 2007/07/24 16:40:14 draeger1 Exp $
//
// create a square matrix and orthogonalize with Gram-Schmidt
//
// usage:  testGram nprow npcol matrix_size

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <map>
#include "Timer.h"
#include "Context.h"
#include "Matrix.h"
#ifdef USE_MPI
#include <mpi.h>
#endif
using namespace std;

int main(int argc, char **argv)
{

  // set up map of timers
  map<string,Timer> tmap;

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

     int nprow, npcol, m, n;
    if (argc == 5) {
      nprow = atoi(argv[1]);
      npcol = atoi(argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
    }
    else {
      cerr << "Usage:  testGram nprow npcol m n" << endl;
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
    
    int mb = m/nprow + (m%nprow > 0 ? 1 : 0);
    int nb = n/npcol + (n%npcol > 0 ? 1 : 0);
    
    ComplexMatrix c(ctxt,m,n,mb,nb);
    ComplexMatrix s(ctxt,n,n,nb,nb);

    // randomize initial values
    {
       srand48(ctxt.myproc());
       for ( int in = 0; in < nb; in++ ) {
          complex<double>* p = c.valptr(mb*in);
          for ( int im = 0; im < mb; im++ ) {
             double dre = drand48();
             double dim = drand48();
             p[im] = 0.02 * complex<double>(dre,dim);
          }
       }
    }
    tmap["init"].stop();

    tmap["herk"].start();
    s.herk('l','c',1.0,c,0.0);
    tmap["herk"].stop();

    tmap["potrf"].start();
    s.potrf('l');
    tmap["potrf"].stop();

    tmap["trsm"].start();
    c.trsm('r','l','c','n',1.0,s);
    tmap["trsm"].stop();

    
    if (mype == 0) 
      cout << "Done." << endl;
    tmap["total"].stop();

    //ewd DEBUG:  check that orthogonalization was successful
    if (true)
    {
       if (mype == 0) 
          cout << "Checking orthogonality by computing overlap matrix:  s = " << n << " x " << n << ", nb = " << nb << endl;
       s.gemm('c','n',1.0,c,c,0.0);

       cout << "OVERLAP, mype = " << mype << ", localsize = " << s.localsize() << endl;
       
       if (s.localsize() > 0)
       {
          for ( int in = 0; in < nb; in++ ) {
             complex<double>* sp = s.valptr(nb*in);
             for ( int im = 0; im < nb; im++ ) {
                if (im == in)
                {
                   if (abs(real(sp[im])-1.0) > 1.E-6)
                      cout << "ORTHOG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", s val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                }
                else
                {
                   if (abs(real(sp[im])) > 1.E-6 || abs(imag(sp[im]) > 1.E-6))
                      cout << "ORTHOG ERROR, mype = " << mype << ", im = " << im << ", in = " << in << ", s val = " << real(sp[im]) << "  " << imag(sp[im]) << endl;
                }
             }
          }
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //ewd DEBUG    

    
    for ( map<string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ )
    {
       double time = (*i).second.cpu();
       double tmin = time;
       double tmax = time;
    
       ctxt.dmin(1,1,&tmin,1);
       ctxt.dmax(1,1,&tmax,1);
       if (mype == 0)
       { 
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
