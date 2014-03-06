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
// DQBPCGWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "DQBPCGWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Preconditioner.h"
#include "ChargeDensity.h"
#include "AtomSet.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
DQBPCGWavefunctionStepper::DQBPCGWavefunctionStepper(Wavefunction& wf,
                                                     Preconditioner& p,const AtomSet& atoms,
                                                     const ChargeDensity& cd_,vector<vector<double> >& v_r,
                                                     TimerMap& tmap) :
    WavefunctionStepper(wf,tmap), prec_(p), hpsi_(wf,atoms,cd_,v_r), reswf_(wf), hreswf_(wf)
{

  nkp_ = wf_.nkp();
  nspin_ = wf_.nspin();

  cell_moved();
  
}

////////////////////////////////////////////////////////////////////////////////
void DQBPCGWavefunctionStepper::cell_moved()
{
   hpsi_.cell_moved(wf_);
}
////////////////////////////////////////////////////////////////////////////////
void DQBPCGWavefunctionStepper::update(Wavefunction& dwf)
{
   for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
      if (wf_.spinactive(ispin)) {
         for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
            if (wf_.kptactive(ikp)) {
               assert(wf_.sd(ispin,ikp) != 0);
               
               // compute A = V^T H V  and descent direction HV - VA
               if ( wf_.sd(ispin,ikp)->basis().real() ) {
                  // proxy real matrices c, cp
                  DoubleMatrix c_proxy(wf_.sd(ispin,ikp)->c());
                  DoubleMatrix cp_proxy(dwf.sd(ispin,ikp)->c());
                  DoubleMatrix w_proxy(reswf_.sd(ispin,ikp)->c());
                  DoubleMatrix hw_proxy(hreswf_.sd(ispin,ikp)->c());
             
                  DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                                 c_proxy.nb(),c_proxy.nb());
                  DoubleMatrix g(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                                 c_proxy.nb(),c_proxy.nb());

                  const int mloc = wf_.sd(ispin,ikp)->c().mloc();
                  const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
                  const int nloc = wf_.sd(ispin,ikp)->c().nloc();

                  tmap_["dqbpcg_residual"].start();

                  // factor 2.0 in next line: G and -G
                  a.gemm('t','n',2.0,c_proxy,cp_proxy,0.0);
                  // rank-1 update correction
                  a.ger(-1.0,c_proxy,0,cp_proxy,0);
          
                  // w = cp - c * a
                  w_proxy = cp_proxy;
                  w_proxy.gemm('n','n',-1.0,c_proxy,a,1.0);

                  // w now contains the descent direction (HV-VA)
                  tmap_["dqbpcg_residual"].stop();

                  // apply preconditioner to dwf                  
                  tmap_["dqbpcg_prec"].start();
                  const valarray<double>& diag = prec_.diag(ispin,ikp);
                  double* wc = (double*) w_proxy.valptr();

                  for ( int n = 0; n < nloc; n++ ) {
                     double* wcn = &wc[2*mloc*n];
                     for ( int i = 0; i < ngwl; i++ ) {
                        const double fac = diag[i];
                        const double f0 = fac * wcn[2*i];
                        const double f1 = fac * wcn[2*i+1];
                        wcn[2*i] = f0;
                        wcn[2*i+1] = f1;
                     }
                  }
                  // w now contains the preconditioned descent
                  // direction K(HV-VA)
                  tmap_["dqbpcg_prec"].stop();

                  tmap_["dqbpcg_iter0"].start();
                  // w -> (w - psi^t psi dwf)
                  g.gemm('t','n',2.0,c_proxy,cp_proxy,0.0);
                  g.ger(-1.0,c_proxy,0,cp_proxy,0);
                  w_proxy.gemm('n','n',-1.0,c_proxy,g,1.0);
                  
                  // compute H*residual
                  hpsi_.compute(reswf_,hreswf_);


                  tmap_["dqbpcg_iter0"].stop();
                  tmap_["dqbpcg_iterloop"].start();



                  tmap_["dqbpcg_iterloop"].stop();


                  
               }
               else {
                  ComplexMatrix &c_proxy = wf_.sd(ispin,ikp)->c();
                  ComplexMatrix &cp = dwf.sd(ispin,ikp)->c();
                  ComplexMatrix &w_proxy = reswf_.sd(ispin,ikp)->c();
                  ComplexMatrix &hw = hreswf_.sd(ispin,ikp)->c();
                  ComplexMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),c_proxy.nb(),c_proxy.nb());
                  ComplexMatrix g(c_proxy.context(),c_proxy.n(),c_proxy.n(),c_proxy.nb(),c_proxy.nb());
                  const int mloc = wf_.sd(ispin,ikp)->c().mloc();
                  const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
                  const int nloc = wf_.sd(ispin,ikp)->c().nloc();

                  tmap_["dqbpcg_residual"].start();
                  a.gemm('c','n',1.0,c_proxy,cp,0.0);

                  // w = cp - c * a
                  w_proxy.gemm('n','n',-1.0,c_proxy,a,1.0);
                  // dwf.sd->c() now contains the descent direction (HV-VA)
                  tmap_["dqbpcg_residual"].stop();
             
                  // apply preconditioner to dwf             
                  tmap_["dqbpcg_prec"].start();
                  const valarray<double>& diag = prec_.diag(ispin,ikp);
                  double* wc = (double*) w_proxy.valptr();
                  for ( int n = 0; n < nloc; n++ ) {
                     double* wcn = &wc[2*mloc*n];
                     for ( int i = 0; i < ngwl; i++ ) {
                        const double fac = diag[i];
                        const double f0 = fac * wcn[2*i];
                        const double f1 = fac * wcn[2*i+1];
                        wcn[2*i] = f0;
                        wcn[2*i+1] = f1;
                     }
                  }
                  // w now contains the preconditioned descent
                  // direction K(HV-VA)
                  tmap_["dqbpcg_prec"].stop();


                  tmap_["dqbpcg_iter0"].start();                  
                  // dwf -> (dwf - psi^t psi dwf)
                  g.gemm('c','n',1.0,c_proxy,cp,0.0);
                  w_proxy = cp;
                  w_proxy.gemm('n','n',-1.0,c_proxy,g,1.0);
                  
                  // compute H*residual
                  hpsi_.compute(reswf_,hreswf_);





                  
                  tmap_["dqbpcg_iter0"].stop();
                  tmap_["dqbpcg_iterloop"].start();



                  tmap_["dqbpcg_iterloop"].stop();
                  
             
               }
            }
         }
      }
   }
}
