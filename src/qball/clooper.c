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
   int ii;
   for (ii=0; ii<size; ++ii)
      vout[0] += v1[ii]*v2[ii];
   return;
}

void cdLoop2(const int size, double* v1, double* v2, double* vout)
{
   int ii;
   for (ii=0; ii<size; ++ii)
   {
      vout[0] += v1[2*ii]*v2[2*ii] - v1[2*ii+1]*v2[2*ii+1];
      vout[1] += v1[2*ii+1]*v2[2*ii] + v1[2*ii]*v2[2*ii+1];
   }
   return;
}

void myzdotc(const int size, double complex* v1, double complex* v2, double complex* vout)
{
   int ii;
   for (ii=0; ii<size; ++ii)
      vout[0] += conj(v1[ii])*v2[ii];
   return;
}
