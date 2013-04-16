// $Id: testPgemm.C,v 1.1 2007/07/24 16:40:14 draeger1 Exp $
//
// usage:  testPgemm nprow npcol m n

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
#ifdef USE_CTF
#include "cyclopstf.h"
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
#endi

#ifdef USE_CTF
  {
    int myRank,numPes;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numPes);
    CTF_init(MPI_COMM_WORLD, MACHINE_BGQ, myRank, numPes); 
    CTF_init_complex(MPI_COMM_WORLD, MACHINE_BGQ, myRank, numPes); 
    //CTF_init(MPI_COMM_WORLD, myRank, numPes); 
    //CTF_init_complex(MPI_COMM_WORLD, myRank, numPes); 
  }    
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
        cerr << "Usage:  testEigenSolvers nprow npcol m n" << endl;
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
    
     ComplexMatrix c1(ctxt,m,n,mb,nb);
     ComplexMatrix c2(ctxt,m,n,mb,nb);
     ComplexMatrix s(ctxt,n,n,nb,nb);

     // randomize initial values
     {
        srand48(ctxt.myproc());
        for ( int in = 0; in < nb; in++ ) {
           complex<double>* p1 = c1.valptr(mb*in);
           complex<double>* p2 = c2.valptr(mb*in);
           for ( int im = 0; im < mb; im++ ) {
              double dre1 = drand48();
              double dim1 = drand48();
              p1[im] = 0.02 * complex<double>(dre1,dim1);
              double dre2 = drand48();
              double dim2 = drand48();
              p2[im] = 0.02 * complex<double>(dre2,dim2);
           }
        }
     }
     tmap["init"].stop();

     tmap["pzgemm"].start();
     s.gemm('c','n',1.0,c1,c2,0.0);
     tmap["pzgemm"].stop();
    
    if (mype == 0) 
      cout << "Done." << endl;
    tmap["total"].stop();


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

#ifdef USE_CTF
  CTF_exit();
#endif  
#ifdef USE_MPI
  MPI_Finalize();
#endif

}
