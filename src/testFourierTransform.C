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
//
// testFourierTransform.C
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "Context.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "Timer.h"

int fft_flops(int n)
{
  return 5.0 * n * log((double) n) / log(2.0);
}

int main(int argc, char **argv)
{
  Timer tm;
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that Context objects get destructed before
  // the MPI_Finalize call
  {
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);

  Context ctxt_global;
  int mype = ctxt_global.mype();
  //cout << " Process " << mype << " on " << processor_name << endl;

  D3vector a,b,c,kpoint;
  double ecut;
  if ( argc == 5 )
  {
    a = D3vector(atof(argv[1]),0,0);
    b = D3vector(0,atof(argv[2]),0);
    c = D3vector(0,0,atof(argv[3]));
    ecut = atof(argv[4]);
  }
  else if ( argc == 11 )
  {
    a = D3vector(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    b = D3vector(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    c = D3vector(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    ecut = atof(argv[10]);
  }
  else if ( argc == 14 )
  {
    a = D3vector(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    b = D3vector(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    c = D3vector(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    ecut = atof(argv[10]);
    kpoint = D3vector(atof(argv[11]),atof(argv[12]),atof(argv[13]));
  }
  else
  {
    cout << " use: testFourierTransform a b c ecut(a.u.) [kpoint] " << endl;
  }
  UnitCell cell(a,b,c);
  const double omega = cell.volume();
  
  //cout << " ctxt_global: " << ctxt_global;
  //cout << " ctxt_global.comm(): " << ctxt_global.comm() << endl;
  Context ctxt(ctxt_global.size(),1);
  //cout << " ctxt: " << ctxt;
  //cout << " ctxt.comm(): " << ctxt.comm() << endl;
  Basis basis(ctxt,kpoint);
  basis.resize(cell,cell,ecut);
  
  // transform and interpolate as for wavefunctions

  FourierTransform ft2(basis,2*basis.np(0),2*basis.np(1),2*basis.np(2));
  vector<complex<double> > f2(ft2.np012loc());

  vector<complex<double> > x(basis.localsize());
  vector<complex<double> > x1(basis.localsize());
  vector<complex<double> > x2(basis.localsize());

  double flops = 2*basis.nrod_loc() *      fft_flops(ft2.np2()) +
                 ft2.np1()/2 * ft2.np2() * fft_flops(ft2.np0()) +
                 ft2.np0()   * ft2.np2() * fft_flops(ft2.np1());
  if ( ctxt.oncoutpe() )
  {
    cout << " wfbasis.size() = " << basis.size() << endl;
    cout << " wfbasis.np() = " << basis.np(0) << " " << basis.np(1)
         << " " << basis.np(2) << endl;
    //cout << " flop count: " << flops << endl;
    cout << " wfbasis.nrod_loc(): " << basis.nrod_loc() << endl;
    cout << " zvec.size: " 
         << 2*basis.nrod_loc()*ft2.np2() * sizeof(complex<double>)
         << endl;
  }

  cout.setf(ios::fixed,ios::floatfield);
  cout.setf(ios::right,ios::adjustfield);
  cout << setprecision(6);
  
  const double rc = 1.0;
#if 1
  // Initialize with Fourier coefficients of a normalized gaussian distribution
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    double g2 = basis.g2(i);
    double y = 1.0/omega * exp( -0.25 * g2 * rc*rc );
    x[i] = y;
    x1[i] = y;
    x2[i] = y;
    // x[i] = complex<double>(y,y);
  }
#endif

  /*EWD DEBUG
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.forward(&f2[0],&x[0]);
  tm.stop();
  cout << " fwd1: tm_f_fft:    " << ft2.tm_f_fft.real() << endl;
  cout << " fwd1: tm_f_mpi:    " << ft2.tm_f_mpi.real() << endl;
  cout << " fwd1: tm_f_pack:   " << ft2.tm_f_pack.real() << endl;
  cout << " fwd1: tm_f_unpack: " << ft2.tm_f_unpack.real() << endl;
  cout << " fwd1: tm_f_zero:   " << ft2.tm_f_zero.real() << endl;
  cout << " fwd1: tm_f_map:    " << ft2.tm_f_map.real() << endl;
  cout << " fwd1: tm_f_total:  " << ft2.tm_f_fft.real() +
                                    ft2.tm_f_mpi.real() +
                                    ft2.tm_f_pack.real() +
                                    ft2.tm_f_unpack.real() +
                                    ft2.tm_f_zero.real() +
                                    ft2.tm_f_map.real() << endl;
  cout << " fwd1 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.backward(&x[0],&f2[0]);
  tm.stop();
  cout << " bwd1: tm_b_fft:    " << ft2.tm_b_fft.real() << endl;
  cout << " bwd1: tm_b_mpi:    " << ft2.tm_b_mpi.real() << endl;
  cout << " bwd1: tm_b_pack:   " << ft2.tm_b_pack.real() << endl;
  cout << " bwd1: tm_b_unpack: " << ft2.tm_b_unpack.real() << endl;
  cout << " bwd1: tm_b_zero:   " << ft2.tm_b_zero.real() << endl;
  cout << " bwd1: tm_b_map:    " << ft2.tm_b_map.real() << endl;
  cout << " bwd1: tm_b_total:  " << ft2.tm_b_fft.real() +
                                    ft2.tm_b_mpi.real() +
                                    ft2.tm_b_pack.real() +
                                    ft2.tm_b_unpack.real() +
                                    ft2.tm_b_zero.real() +
                                    ft2.tm_b_map.real() << endl;
  cout << " bwd1 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.forward(&f2[0],&x[0]);
  tm.stop();
  cout << " fwd2: tm_f_fft:    " << ft2.tm_f_fft.real() << endl;
  cout << " fwd2: tm_f_mpi:    " << ft2.tm_f_mpi.real() << endl;
  cout << " fwd2: tm_f_pack:   " << ft2.tm_f_pack.real() << endl;
  cout << " fwd2: tm_f_unpack: " << ft2.tm_f_unpack.real() << endl;
  cout << " fwd2: tm_f_zero:   " << ft2.tm_f_zero.real() << endl;
  cout << " fwd2: tm_f_map:    " << ft2.tm_f_map.real() << endl;
  cout << " fwd2: tm_f_total:  " << ft2.tm_f_fft.real() +
                                    ft2.tm_f_mpi.real() +
                                    ft2.tm_f_pack.real() +
                                    ft2.tm_f_unpack.real() +
                                    ft2.tm_f_zero.real() +
                                    ft2.tm_f_map.real() << endl;

  //cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
  //     << " " << 2*basis.np(2) << " ";
  cout << " fwd2 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.backward(&x[0],&f2[0]);
  tm.stop();
  cout << " bwd2: tm_b_fft:    " << ft2.tm_b_fft.real() << endl;
  cout << " bwd2: tm_b_mpi:    " << ft2.tm_b_mpi.real() << endl;
  cout << " bwd2: tm_b_pack:   " << ft2.tm_b_pack.real() << endl;
  cout << " bwd2: tm_b_unpack: " << ft2.tm_b_unpack.real() << endl;
  cout << " bwd2: tm_b_zero:   " << ft2.tm_b_zero.real() << endl;
  cout << " bwd2: tm_b_map:    " << ft2.tm_b_map.real() << endl;
  cout << " bwd2: tm_b_total:  " << ft2.tm_b_fft.real() +
                                    ft2.tm_b_mpi.real() +
                                    ft2.tm_b_pack.real() +
                                    ft2.tm_b_unpack.real() +
                                    ft2.tm_b_zero.real() +
                                    ft2.tm_b_map.real() << endl;

  //cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
  //     << " " << 2*basis.np(2) << " ";
  cout << " bwd2 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  
  // double transform
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.forward(&f2[0],&x1[0],&x2[0]);
  tm.stop();
  if ( ctxt.oncoutpe() ) {
     cout << " fwd3: tm_f_fft:    " << ft2.tm_f_fft.real() << endl;
     cout << " fwd3: tm_f_mpi:    " << ft2.tm_f_mpi.real() << endl;
     cout << " fwd3: tm_f_pack:   " << ft2.tm_f_pack.real() << endl;
     cout << " fwd3: tm_f_unpack: " << ft2.tm_f_unpack.real() << endl;
     cout << " fwd3: tm_f_zero:   " << ft2.tm_f_zero.real() << endl;
     cout << " fwd3: tm_f_map:    " << ft2.tm_f_map.real() << endl;
     cout << " fwd3: tm_f_total:  " << ft2.tm_f_fft.real() +
         ft2.tm_f_mpi.real() +
         ft2.tm_f_pack.real() +
         ft2.tm_f_unpack.real() +
         ft2.tm_f_zero.real() +
         ft2.tm_f_map.real() << endl;
     cout << " fwd3 time: " << tm.cpu() << " / " << tm.real()
          << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  }
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.backward(&x1[0],&x2[0],&f2[0]);
  tm.stop();
  if ( ctxt.oncoutpe() ) {
     cout << " bwd3: tm_b_fft:    " << ft2.tm_b_fft.real() << endl;
     cout << " bwd3: tm_b_mpi:    " << ft2.tm_b_mpi.real() << endl;
     cout << " bwd3: tm_b_pack:   " << ft2.tm_b_pack.real() << endl;
     cout << " bwd3: tm_b_unpack: " << ft2.tm_b_unpack.real() << endl;
     cout << " bwd3: tm_b_zero:   " << ft2.tm_b_zero.real() << endl;
     cout << " bwd3: tm_b_map:    " << ft2.tm_b_map.real() << endl;
     cout << " bwd3: tm_b_total:  " << ft2.tm_b_fft.real() +
         ft2.tm_b_mpi.real() +
         ft2.tm_b_pack.real() +
         ft2.tm_b_unpack.real() +
         ft2.tm_b_zero.real() +
         ft2.tm_b_map.real() << endl;
     cout << " bwd3 time: " << tm.cpu() << " / " << tm.real()
          << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  }  
#if 1
  //////////////////////////////////////////////////////////////////////////////
  // Integration of a 2-norm normalized plane wave
  //////////////////////////////////////////////////////////////////////////////
  
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    x[i] = 0.0;
  }
  if ( ctxt.myproc() == 0 ) x[1] = 1.0/sqrt(2.0);
  
  ft2.backward(&x[0],&f2[0]);

#if 0
  for ( int i = 0; i < basis.localsize(); i++ )
    cout << basis.kv(3*i) << " " << basis.kv(3*i+1) << " " << basis.kv(3*i+2)
         << "     " << x[i] << endl;
  for ( int i = 0; i < ft.np0(); i++ )
    for ( int j = 0; j < ft.np1(); j++ )
      for ( int k = 0; k < ft.np2_loc(); k++ )
        cout << mype << ": "
             << i << " " << j << " " << k+ft.np2_first() << " "
             << f[ft.index(i,j,k)] << endl;
#endif
             
  // integral of f^2 in r space must be 1.0
  double sum=0.0, tsum = 0.0;
  for ( int i = 0; i < f2.size(); i++ ) 
    tsum += norm(f2[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  
  if ( ctxt.oncoutpe() )
     cout << " sum pw^2: " << sum / ft2.np012() << endl;
  
  //////////////////////////////////////////////////////////////////////////////
  // Integration of a 2-norm normalized gaussian
  //////////////////////////////////////////////////////////////////////////////
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    double g2 = basis.g2(i);
    x[i] = 1.0 / sqrt(omega) * pow(2.0*M_PI*rc*rc,0.75) * 
           exp( -0.25 * g2 * rc*rc );
  }

  // Compute norm in g space
  double gnorm = 0.0;
  for ( int i = 0; i < basis.localsize(); i++ )
    gnorm += 2.0 * norm(x[i]);
  if ( ctxt.oncoutpe() )
    gnorm -= norm(x[0]);
  ctxt.dsum(1,1,&gnorm,1);
  if ( ctxt.oncoutpe() )
     cout << " gaussian gnorm: " << gnorm << endl;
  
  ft2.backward(&x[0],&f2[0]);
  
//   for ( int i = 0; i < basis.localsize(); i++ )
//     cout << basis.kv(3*i) << " " << basis.kv(3*i+1) << " " << basis.kv(3*i+2)
//          << "     " << x[i] << endl;
//   for ( int i = 0; i < ft2.np0(); i++ )
//     for ( int j = 0; j < ft2.np1(); j++ )
//       for ( int k = 0; k < ft2.np2_loc(); k++ )
//         cout << mype << ": "
//              << i << " " << j << " " << k+ft2.np2_first() << " "
//              << f2[ft2.index(i,j,k)] << endl;
             
  // integral of gaussian^2 in r space must be 1.0
  tsum = 0.0;
  for ( int i = 0; i < f2.size(); i++ ) 
    tsum += norm(f2[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  
  if ( ctxt.oncoutpe() )
     cout << " gaussian rnorm: " << sum / ft2.np012() << endl;
#endif
  EWD DEBUG  */
  
  // Define ft from vbasis
  Basis vbasis(ctxt,kpoint);
  vbasis.resize(cell,cell,4.0*ecut);
  if ( ctxt.oncoutpe() )
     cout << " vbasis.np() = " << vbasis.np(0) << " " << vbasis.np(1)
          << " " << vbasis.np(2) << endl;

  /* EWD DEBUG
  FourierTransform vft(basis,vbasis.np(0),vbasis.np(1),vbasis.np(2));
  vector<complex<double> > vf(vft.np012loc());
  vft.backward(&x[0],&vf[0]);
  // integral of gaussian^2 in r space must be 1.0
  tsum = 0.0;
  for ( int i = 0; i < vf.size(); i++ ) 
    tsum += norm(vf[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  
  if ( ctxt.oncoutpe() )
     cout << " gaussian rnorm: " << sum / vft.np012() << endl;
  EWD DEBUG  */

  } // Context scope
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
