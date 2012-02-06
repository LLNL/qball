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
#include <vector>
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
