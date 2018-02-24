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
////////////////////////////////////////////////////////////////////////////////
//
// test_fftw.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <qball/Timer.h>

#include <iostream>
#include <complex>
#include <valarray>
using namespace std;
#include <cassert>

#include "fftw.h"

int main(int argc, char**argv)
{
  const int niter = 10;
  const int np = atoi(argv[1]);
  const int nvec = atoi(argv[2]);
  const int ldz = np + 1;

  fftw_plan fwplan, bwplan;

  // resize array zvec holding columns
  valarray<complex<double> > zvec(nvec * ldz);
  
  // initialization of FFT libs

// #define FFTWMEASURE 1
#if FFTWMEASURE
  // FFTWMEASURE
  fwplan = fftw_create_plan(np,FFTW_FORWARD,FFTW_MEASURE|FFTW_IN_PLACE);
  bwplan = fftw_create_plan(np,FFTW_BACKWARD,FFTW_MEASURE|FFTW_IN_PLACE);
#else
  // FFTW_ESTIMATE
  fwplan = fftw_create_plan(np,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
  bwplan = fftw_create_plan(np,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif

  Timer t_fwd,t_bwd;

  for ( int iter = 0; iter < niter; iter++ )
  {
  t_bwd.start();

   /* 
    * void fftw(fftw_plan plan, int howmany,
    *    FFTW_COMPLEX *in, int istride, int idist,
    *    FFTW_COMPLEX *out, int ostride, int odist);
    */
  int ntrans = nvec;
  int inc1 = 1;
  int inc2 = ldz;
  fftw(bwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
  t_bwd.stop();
  t_fwd.start();
  fftw(fwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
  t_fwd.stop();
  }

  fftw_destroy_plan(fwplan);
  fftw_destroy_plan(bwplan);

  cout << " fwd: time per transform (in-place,generic)" 
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_fwd.real()/(niter*nvec) << " microseconds" << endl;

  cout << " bwd: time per transform (in-place,generic)" 
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_bwd.real()/(niter*nvec) << " microseconds" << endl;

#if 1
  // Use out-of-place, specific plan
  
  valarray<complex<double> > zvec_out(zvec.size());
  t_bwd.reset();
  t_fwd.reset();

  fwplan = fftw_create_plan_specific(np,
    FFTW_FORWARD,FFTW_ESTIMATE|FFTW_OUT_OF_PLACE,
    (FFTW_COMPLEX*)&zvec[0],1,(FFTW_COMPLEX*)&zvec_out[0],1);
  bwplan = fftw_create_plan_specific(np,
    FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_OUT_OF_PLACE,
    (FFTW_COMPLEX*)&zvec[0],1,(FFTW_COMPLEX*)&zvec_out[0],1);
    
  for ( int iter = 0; iter < niter; iter++ )
  {

    int ntrans = nvec;
    int inc1 = 1;
    int inc2 = ldz;
    t_bwd.start();
    fftw(bwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                       (FFTW_COMPLEX*)&zvec_out[0],inc1,inc2);
    t_bwd.stop();
  
    t_fwd.start();
    fftw(fwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                       (FFTW_COMPLEX*)&zvec_out[0],inc1,inc2);
    t_fwd.stop();

  }

  fftw_destroy_plan(fwplan);
  fftw_destroy_plan(bwplan);
  
  cout << " fwd: time per transform (out-of-place,specific)" 
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_fwd.real()/(niter*nvec) << " microseconds" << endl;

  cout << " bwd: time per transform (out-of-place,specific)" 
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_bwd.real()/(niter*nvec) << " microseconds" << endl;
#endif
  return 0;
}
