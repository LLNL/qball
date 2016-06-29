// $Id: testMatrix.C,v 1.3 2006/01/10 01:15:38 draeger1 Exp $
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

  int nprow, npcol;
  
  nprow = 8;
  npcol = npes/nprow;

  int m,n,mb,nb;

  m = 128;
  n = 2;
  mb = m/nprow;
  if (m%nprow != 0) mb++;
  nb = n/npcol;
  if (n%npcol != 0) nb++;
  
  Context ctxt(nprow,npcol);
  DoubleMatrix a(ctxt,m,n,mb,nb);

  if (mype == 0)
    cout << "Process grid:  " << npes << " pes, " << nprow << " x " << npcol << endl;
  cout << "matrix info, mype = " << mype << ", m x n = " << a.m() << " x " << a.n() << ", mb x nb = " << a.mb() << " x " << a.nb() << ", size = " << a.size() << ", myrow = " << ctxt.myrow() << ", mycol = " << ctxt.mycol() << endl;
  
}
