#include "clooper.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <mpi.h>

// ewd:  These functions are loops that the BG/Q C++ compiler can't SIMDize.
// ewd:  Currently, these are only called when -DBGQ is set, so feel free
// ewd:  to put machine-specific code in here.

void cdLoop(const int size, double complex* v1, double complex* v2, double complex* vout)
{
   for (int ii=0; ii<size; ++ii)
      vout[0] += v1[ii]*v2[ii];
   return;
}

void cdLoop2(const int size, double* v1, double* v2, double* vout)
{
   for (int ii=0; ii<size; ++ii)
   {
      vout[0] += v1[2*ii]*v2[2*ii] - v1[2*ii+1]*v2[2*ii+1];
      vout[1] += v1[2*ii+1]*v2[2*ii] + v1[2*ii]*v2[2*ii+1];
   }
   return;
}
