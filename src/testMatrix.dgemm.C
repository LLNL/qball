// $Id: testMatrix.old2.C,v 1.1 2005/12/05 23:11:16 draeger1 Exp $
//
// test Matrix
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testMatrix input_file
// input_file:
// m_a n_a mb_a nb_a transa
// m_b n_b mb_b nb_b transb
// m_c n_c mb_c nb_c
//

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
#include "blas.h"
#include "omp.h"
#include "Timer.h"

#ifdef BGQ
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif

using namespace std;


#ifdef USE_MPI
#include <mpi.h>
#endif

const double nrandinv = 1./(1.0*RAND_MAX + 1.0);
const double maxrand = 0.000001;  // maximum random perturbation to identity matrix

int main(int argc, char **argv)
{

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

   char* infilename = argv[1];
   ifstream infile(infilename);
   assert(argc == 2);
   Timer tm;
   int m_a, n_a, mb_a, nb_a;
   int m_b, n_b, mb_b, nb_b;
   int m_c, n_c, mb_c, nb_c;
   char ta, tb;
   if(mype == 0)
   {
      infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
      cout<<"m_a="<<m_a<<", n_a="<<n_a<<endl;
      infile >> m_b >> n_b >> mb_b >> nb_b >> tb;
      cout<<"m_b="<<m_b<<", n_b="<<n_a<<endl;
      infile >> m_c >> n_c >> mb_c >> nb_c;
      cout<<"m_c="<<m_c<<", n_c="<<n_c<<endl;
   }
#ifdef USE_MPI
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

   if ( ta == 'N' ) ta = 'n';
   if ( tb == 'N' ) tb = 'n';

   double dzero = 0.0;
   double done = 1.0;
   char ct='t';
   char cn='n';

   int kk = m_a;
   int mm = m_c;
   int nn = n_c;

   vector<double> avec(mm*kk);
   vector<double> bvec(nn*kk);
   vector<double> cvec(mm*nn);
   for (int ii=0; ii<avec.size(); ii++)
      avec[ii] = rand()*nrandinv*maxrand;
   for (int ii=0; ii<bvec.size(); ii++)
      bvec[ii] = rand()*nrandinv*maxrand;
   for (int ii=0; ii<cvec.size(); ii++)
      cvec[ii] = rand()*nrandinv*maxrand;
               
   tm.start();
   HPM_Start("dgemm1");
   dgemm(&ct,&cn,&mm,&nn,&kk,&done,&avec[0],&kk,&bvec[0],&kk,&dzero,&cvec[0],&mm);
   HPM_Stop("dgemm1");
   tm.stop();

   int nthreads = omp_get_max_threads();
   if (mype == 0)
      cout << "M = " << m_c << " N = " << n_c << " K = " << m_a << " dgemm time = " << setprecision(5) << setw(8) << tm.real()<< " sec, GFlops = " << (2.0e-9*m_c*n_c*m_a) / tm.real() << " on " << npes << " pes, " << nthreads << " threads" << endl;
 
#ifdef USE_MPI
   MPI_Finalize();
#endif
}
