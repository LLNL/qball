#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <mpi.h>
#include <vector>
using namespace std;

#include "qbLink.h"

int main(int argc, char **argv) {

  int mype;
  int npes;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  //if (mype == 0) 
  //  cout << "qb-nruns init, mype = " << mype << ", npes = " << npes << ", argc = " << argc << endl;

  int ninputs = argc-1;
  int npes_sub = npes/ninputs;
  assert(npes_sub > 0);
  if (npes_sub%2 !=0 && mype == 0) 
    cout << "<WARNING> Odd number of pes for each run:  npes_sub = " << npes_sub << " </WARNING>" << endl;
  vector<qbLink*> qb(ninputs);
  vector<string> inpfiles(ninputs);
  vector<string> outfiles(ninputs);
  int firstpe = 0;
  for (int i=0; i<ninputs; i++) {
    int lastpe = firstpe + npes_sub-1;
    string inp(argv[i+1]);
    inpfiles[i] = inp;
    outfiles[i] = inp + ".out";
    qbLink* tmpqb = new qbLink(outfiles[i],firstpe,lastpe);
    qb[i] = tmpqb;
    firstpe += npes_sub;
  }

  for (int i=0; i<ninputs; i++) 
    if (qb[i]->active())
      qb[i]->processInputFile(inpfiles[i]);

  for (int i=0; i<ninputs; i++) 
    if (qb[i] != 0)
      delete qb[i];

  MPI_Finalize();
}






