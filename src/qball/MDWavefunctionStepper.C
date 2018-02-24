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
// MDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "MDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
MDWavefunctionStepper::MDWavefunctionStepper(Wavefunction& wf,
  Wavefunction *wfv, double dt, double dt2bye, TimerMap& tmap) :
  wfv_(wfv), dt_(dt), dt2bye_(dt2bye), WavefunctionStepper(wf,tmap)
{
  assert(wfv!=0);
}

////////////////////////////////////////////////////////////////////////////////
void MDWavefunctionStepper::update(Wavefunction& dwf)
{
  // Verlet update of wf using force dwf and wfm stored in *wfv
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    if (wf_.spinactive(ispin))
    {
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      {
        if (wf_.kptactive(ikp))
        {

          tmap_["md_update_wf"].start();
          // Verlet update of wf
          // cp = c + (c - cm) - dt2/m * hpsi
          // This is implemented (for each coefficient) as:
          // cp = 2*c - cm - dt2bye * hpsi
          // cm = c
          // c = cp
          SlaterDet* sd = wf_.sd(ispin,ikp);
          double* cptr = (double*) sd->c().valptr();
          double* cptrm = (double*) wfv_->sd(ispin,ikp)->c().valptr();
          const double* dcptr =
              (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          for ( int n = 0; n < nloc; n++ )
          {
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cm = &cptrm[2*mloc*n];
            const double* dc = &dcptr[2*mloc*n];
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cmtmp = cm[2*i];
              const double cmtmp1 = cm[2*i+1];
              const double dctmp = dc[2*i];
              const double dctmp1 = dc[2*i+1];
              const double cptmp = 2.0*ctmp -  cmtmp -  dt2bye_ * dctmp;
              const double cptmp1= 2.0*ctmp1 - cmtmp1 - dt2bye_ * dctmp1;
              
              c[2*i]    = cptmp;
              c[2*i+1]  = cptmp1;
              cm[2*i]   = ctmp;
              cm[2*i+1] = ctmp1;
            }
          }
          tmap_["md_update_wf"].stop();
        
          //tmap_["riccati"].start();
          //wf_.sd(ispin,ikp)->riccati(*(wfv_->sd(ispin,ikp)),wf_.ultrasoft());
          //tmap_["riccati"].stop();
        }
      }
    }
  }
  ekin_em_ = ekin_ep_;
  ekin_ep_ = ekin_eh();
}

////////////////////////////////////////////////////////////////////////////////
void MDWavefunctionStepper::compute_wfm(Wavefunction& dwf)
{
  // Compute wfm for first MD step using wf, wfv and dwf (= Hpsi)
  // Replace then wfv by wfm

  // Compute cm using c and wavefunction velocity
  // cm = c - dt * v - 0.5 * dt2/m * hpsi
  // replace wfv by wfm
  const double half_dt2bye = 0.5 * dt2bye_;
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    if (wf_.spinactive(ispin))
    {
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      {
        if (wf_.kptactive(ikp))
        {

          SlaterDet* sd = wf_.sd(ispin,ikp);
          
          double* cptr = (double*) sd->c().valptr();
          double* cptrv = (double*) wfv_->sd(ispin,ikp)->c().valptr();
          const double* dcptr =
              (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const vector<double>& occ = sd->occ();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          const bool onrow0 = ( wf_.context().myrow() == 0 );
          for ( int n = 0; n < nloc; n++ )
          {
            const int nglobal = sd->c().j(0,n);
            const double occn = occ[nglobal];
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cv = &cptrv[2*mloc*n];
            const double* dc = &dcptr[2*mloc*n];
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cvtmp = cv[2*i];
              const double cvtmp1 = cv[2*i+1];
              const double dctmp = dc[2*i];
              const double dctmp1 = dc[2*i+1];
              cv[2*i]    = ctmp  - dt_ * cvtmp  - half_dt2bye * dctmp;
              cv[2*i+1]  = ctmp1 - dt_ * cvtmp1 - half_dt2bye * dctmp1;
            }
          }
          //tmap_["riccati"].start();
          //wfv_->sd(ispin,ikp)->riccati(*wf_.sd(ispin,ikp),wf_.ultrasoft());
          //tmap_["riccati"].stop();
        }
      }
    }
  }
  ekin_em_ = 0.0;
  ekin_ep_ = ekin_eh();
  // Note: ekin_ep is a first-order approximation of ekin_e using wf and wfm
  // only. This makes ekin_e consistent with the following update() call
  // Note: *wfv_ now contains wf(t-dt)
}

////////////////////////////////////////////////////////////////////////////////
void MDWavefunctionStepper::compute_wfv(Wavefunction& dwf)
{
  // Compute wfv = (wf - wfm)/dt - 0.5*dtbye*dwf

  assert(dt_!=0.0);
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    if (wf_.spinactive(ispin))
    {
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      {
        if (wf_.kptactive(ikp))
        {

          // compute final velocity wfv
          // v = ( c - cm ) / dt - 0.5 * dt/m * hpsi

          // Note: At this point, *wfv_ contains wf(t-dt)
          
          // hpsi must be orthogonal to the subspace spanned by c
          // compute descent direction H psi - psi (psi^T H psi)
        
          SlaterDet* sd = wf_.sd(ispin,ikp);
          if ( sd->basis().real() )
          {
            // proxy real matrices c, cp
            DoubleMatrix c(sd->c());
            DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
          
            DoubleMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
          
            // factor 2.0 in next line: G and -G
            a.gemm('t','n',2.0,c,cp,0.0);
            // rank-1 update correction
            a.ger(-1.0,c,0,cp,0);
          
            // cp = cp - c * a
            cp.gemm('n','n',-1.0,c,a,1.0);
          }
          else
          {
            ComplexMatrix& c = sd->c();
            ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
            ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
            a.gemm('c','n',1.0,c,cp,0.0);
            // cp = cp - c * a
            cp.gemm('n','n',-1.0,c,a,1.0);
          }
          
          const double dt_inv = 1.0/dt_;
          const double half_dtbye = 0.5 * dt2bye_ / dt_;
          double* cptr = (double*) sd->c().valptr();
          double* cptrv = (double*) wfv_->sd(ispin,ikp)->c().valptr();
          const double* dcptr =
              (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          for ( int n = 0; n < nloc; n++ )
          {
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cv = &cptrv[2*mloc*n];
            const double* dc = &dcptr[2*mloc*n];
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cmtmp = cv[2*i];
              const double cmtmp1 = cv[2*i+1];
              const double dctmp = dc[2*i];
              const double dctmp1 = dc[2*i+1];
              
              cv[2*i]   = ( ctmp  - cmtmp  ) * dt_inv
                  - half_dtbye * dctmp;
              cv[2*i+1] = ( ctmp1 - cmtmp1 ) * dt_inv
                  - half_dtbye * dctmp1;
            }
          }
          // Note: *wfv_ now contains the wavefunction velocity
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double MDWavefunctionStepper::ekin_eh(void)
{
  // compute ekin at time t - 0.5*dt using wf and wfm
  tmap_["ekin_e"].start();
  double ekin_e = 0.0;
  // assume that wfv contains wfm
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    if (wf_.spinactive(ispin))
    {
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      {
        if (wf_.kptactive(ikp))
        {
          const double weight = wf_.weight(ikp);
          SlaterDet* sd = wf_.sd(ispin,ikp);
          double* cptr = (double*) sd->c().valptr();
          double* cptrm = (double*) wfv_->sd(ispin,ikp)->c().valptr();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          // compute electronic kinetic energy at time t-1/2
          const bool onrow0 = ( wf_.context().myrow() == 0 );
          const vector<double>& occ = sd->occ();
          for ( int n = 0; n < nloc; n++ )
          {
            const int nglobal = sd->c().j(0,n);
            const double occn = occ[nglobal];
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cm = &cptrm[2*mloc*n];
            double tmpsum = 0.0;
            if ( sd->basis().real() && onrow0 )
            {
              // correct for double counting of G=0 element
              // i=0 coefficient is real, count only real part
              const double ctmp = c[0];
              const double cmtmp = cm[0];
              tmpsum -= 0.5 * (ctmp - cmtmp)*(ctmp - cmtmp);
            }
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cmtmp = cm[2*i];
              const double cmtmp1 = cm[2*i+1];

              tmpsum += (ctmp -cmtmp )*(ctmp -cmtmp ) +
                  (ctmp1-cmtmp1)*(ctmp1-cmtmp1);
            }
            if ( sd->basis().real() )
              // Note: 2 in next line: from (G,-G)
              ekin_e += weight * ( 2.0 * occn / dt2bye_ ) * tmpsum;
            else
              ekin_e += weight * ( occn / dt2bye_ ) * tmpsum;
          }
        }
      }
    }
  }
  wf_.context().dsum(1,1,&ekin_e,1);
  tmap_["ekin_e"].stop();
  return ekin_e;
}
