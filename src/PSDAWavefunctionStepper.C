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
// PSDAWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "PSDAWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Preconditioner.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PSDAWavefunctionStepper::PSDAWavefunctionStepper(Wavefunction& wf, 
  Preconditioner& p, TimerMap& tmap) : 
  WavefunctionStepper(wf,tmap), prec_(p), wf_last_(wf), dwf_last_(wf) {

  nkp_ = wf_.nkp();
  nspin_ = wf_.nspin();
  extrapolate_.resize(nspin_);
  for (int i=0; i<nspin_; i++)
    extrapolate_[i].resize(nkp_);
  for (int i=0; i<nspin_; i++)
    for (int k=0; k<nkp_; k++)
      extrapolate_[i][k] = false;
}

////////////////////////////////////////////////////////////////////////////////
void PSDAWavefunctionStepper::update(Wavefunction& dwf) {
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);

          //ewd DEBUG
          assert(wf_last_.sd(ispin,ikp) != 0);
          assert(dwf_last_.sd(ispin,ikp) != 0);

          // compute A = V^T H V  and descent direction HV - VA
          if ( wf_.sd(ispin,ikp)->basis().real() ) {
            // proxy real matrices c, cp
            DoubleMatrix c_proxy(wf_.sd(ispin,ikp)->c());
            DoubleMatrix cp_proxy(dwf.sd(ispin,ikp)->c());
          
            DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                           c_proxy.nb(),c_proxy.nb());
            
            tmap_["psda_residual"].start();
            // factor 2.0 in next line: G and -G
            a.gemm('t','n',2.0,c_proxy,cp_proxy,0.0);
            // rank-1 update correction
            a.ger(-1.0,c_proxy,0,cp_proxy,0);
          
            // cp = cp - c * a
            if (wf_.ultrasoft()) {
              // cp = cp - spsi * a
              DoubleMatrix spsi(wf_.sd(ispin,ikp)->spsi());
              cp_proxy.gemm('n','n',-1.0,spsi,a,1.0);
            }
            else {
              // cp = cp - c * a
              cp_proxy.gemm('n','n',-1.0,c_proxy,a,1.0);
            }
            tmap_["psda_residual"].stop();
 
            // dwf.sd->c() now contains the descent direction (HV-VA)
            
            // Apply preconditioner K and store -K(HV-VA) in dwf
	    prec_.apply(*dwf.sd(ispin, ikp), ispin, ikp, /* scale = */ -1.0);
 
            // dwf now contains the preconditioned descent
            // direction -K(HV-VA)
            
	    double* c = (double*) wf_.sd(ispin,ikp)->c().valptr();
            double* c_last = (double*) wf_last_.sd(ispin,ikp)->c().valptr();
            double* dc = (double*) dwf.sd(ispin,ikp)->c().valptr();
            double* dc_last = (double*) dwf_last_.sd(ispin,ikp)->c().valptr();
            const int mloc = wf_.sd(ispin,ikp)->c().mloc();
            const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
            const int nloc = wf_.sd(ispin,ikp)->c().nloc();

            tmap_["psda_update_wf"].start();
            // Anderson extrapolation
            if ( extrapolate_[ispin][ikp] ) {
              double theta = 0.0;
              double a = 0.0, b = 0.0;

              //ewd OMP HERE?
              for ( int lj=0; lj < wf_.sd(ispin,ikp)->c().nblocks(); lj++ )
              {
                 for ( int jj=0; jj < wf_.sd(ispin,ikp)->c().nbs(lj); jj++ )
                 {
                    // global state index
                    const int nglobal = wf_.sd(ispin,ikp)->c().j(lj,jj);
                    const int norig = lj*wf_.sd(ispin,ikp)->c().nb()+jj;
                    const vector<double>& occ = wf_.sd(ispin,ikp)->occ();
                    const double occn = occ[nglobal];

                    //ewd DEBUG
                    //const double occn = 1.0;

                    for ( int i = 0; i < 2*mloc; i++ ) {
                       const double f = dc[i +2*mloc*norig];
                       const double delta_f = f - dc_last[i + 2*mloc*norig];
 
                       // accumulate partial sums of a and b
                       // a = delta_F * F
                  
                       //ewd:  weight sums by occupation to improve convergence in metallic systems
                       a += occn * f * delta_f;
                       b += occn * delta_f * delta_f;
                    }
                 }
              }
              // correct for double counting of asum and bsum on first row
              // factor 2.0: G and -G
              a *= 2.0;
              b *= 2.0;
              int kpc = wf_.mysdctxt(ikp);
              assert(kpc >= 0);
              if ( wf_.sdcontext(ispin,kpc)->myrow() == 0 ) {
                 //ewd:  switching away from local kpoint loops
                 //ewdif ( wf_.sd(ispin,ikp)->context().myrow() == 0 ) {
                 for ( int lj=0; lj < wf_.sd(ispin,ikp)->c().nblocks(); lj++ )
                 {
                    for ( int jj=0; jj < wf_.sd(ispin,ikp)->c().nbs(lj); jj++ )
                    {
                       // global state index
                       const int nglobal = wf_.sd(ispin,ikp)->c().j(lj,jj);
                       const int norig = lj*wf_.sd(ispin,ikp)->c().nb()+jj;
                       const vector<double>& occ = wf_.sd(ispin,ikp)->occ();
                       const double occn = occ[nglobal];
                       
                       const int i = 2*mloc*norig;
                       const double f0 = dc[i];
                       const double f1 = dc[i+1];
                       const double delta_f0 = f0 - dc_last[i];
                       const double delta_f1 = f1 - dc_last[i+1];
                       a -= occn * (f0 * delta_f0 + f1 * delta_f1);
                       b -= occn * (delta_f0 * delta_f0 + delta_f1 * delta_f1);
                    }
                 }
              }
              
              // a and b contain the partial sums of a and b
              double tmpvec[2] = { a, b };
              wf_.sdcontext(ispin,kpc)->dsum(2,1,&tmpvec[0],1);
              //ewdwf_.sd(ispin,ikp)->context().dsum(2,1,&tmpvec[0],1);
              a = tmpvec[0];
              b = tmpvec[1];
              // compute theta = - a / b
              if ( b != 0.0 )
                theta = - a / b;
            
              if ( ikp == 0 && ispin == 0 && wf_.sd(ispin,ikp)->context().oncoutpe() )
                cout << "  <!-- Anderson extrapolation: theta = " << theta;

              if ( theta < -1.0 )
                theta = 0.0;

              theta = min(2.0,theta);
              //ewd try tightening restriction here...
              //theta = min(1.0,theta);

              if ( ikp == 0 && ispin == 0 && wf_.sd(ispin,ikp)->context().oncoutpe() )
                cout <<" ( using " << theta << " )" << " -->"<< endl;
                
              // extrapolation
              for ( int i = 0; i < 2*mloc*nloc; i++ ) {
                // x_bar = x_ + theta * ( x_ - xlast_ ) (store in x_)
                const double x = c[i];
                const double xlast = c_last[i];
                const double xbar = x + theta * ( x - xlast );
 
                // f_bar = f + theta * ( f - flast ) (store in f)
                const double f = dc[i];
                const double flast = dc_last[i];
                const double fbar = f + theta * ( f - flast );
                
                c[i] = xbar + fbar;
                c_last[i] = x;
                dc_last[i] = f;
              }
            }
            else {
              // no extrapolation
              for ( int i = 0; i < 2*mloc*nloc; i++ ) {
                // x_ = x_ + f_
                const double x = c[i];
                const double f = dc[i];
                    
                c[i] = x + f;
                c_last[i] = x;
                dc_last[i] = f;
              }
            }
            extrapolate_[ispin][ikp] = true;
            tmap_["psda_update_wf"].stop();
            
            /*
            // do orthogonalization outside of this routine, so ultrasoft functions
            // can update
            
            enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };
            //const ortho_type ortho = GRAM;
            //const ortho_type ortho = LOWDIN;
            const ortho_type ortho = ORTHO_ALIGN;
 
            switch ( ortho ) {
              case GRAM:
                  tmap_["gram"].start();
                  wf_.sd(ispin,ikp)->gram();
                  tmap_["gram"].stop();
                  break;
                
              case LOWDIN:
                  tmap_["lowdin"].start();
                  wf_.sd(ispin,ikp)->lowdin();
                  tmap_["lowdin"].stop();
                  break;
                
              case ORTHO_ALIGN:
                  tmap_["ortho_align"].start();
                  wf_.sd(ispin,ikp)->ortho_align(*wf_last_.sd(ispin,ikp));
                  tmap_["ortho_align"].stop();
                  break;
                
              case RICCATI:
                  tmap_["riccati"].start();
                  wf_.sd(ispin,ikp)->riccati(*wf_last_.sd(ispin,ikp));
                  tmap_["riccati"].stop();
                  break;
            }
            */
          }
          else {
            ComplexMatrix &c_proxy = wf_.sd(ispin,ikp)->c();
            ComplexMatrix &cp = dwf.sd(ispin,ikp)->c();
            ComplexMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),c_proxy.nb(),c_proxy.nb());

            tmap_["psda_residual"].start();
            a.gemm('c','n',1.0,c_proxy,cp,0.0);
            if (wf_.ultrasoft()) {
              // cp = cp - spsi * a
              ComplexMatrix &spsi = wf_.sd(ispin,ikp)->spsi();
              cp.gemm('n','n',-1.0,spsi,a,1.0);
            }
            else {
              // cp = cp - c * a
              cp.gemm('n','n',-1.0,c_proxy,a,1.0);
            }
            tmap_["psda_residual"].stop();
            
            // dwf.sd->c() now contains the descent direction (HV-VA)
          
            // Apply preconditioner K and store -K(HV-VA) in dwf
            tmap_["psda_prec"].start();
            const valarray<double>& diag = prec_.diag(ispin,ikp);
 
            double* c = (double*) wf_.sd(ispin,ikp)->c().valptr();
            double* c_last = (double*) wf_last_.sd(ispin,ikp)->c().valptr();
            double* dc = (double*) dwf.sd(ispin,ikp)->c().valptr();
            double* dc_last = (double*) dwf_last_.sd(ispin,ikp)->c().valptr();
            const int mloc = wf_.sd(ispin,ikp)->c().mloc();
            const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
            const int nloc = wf_.sd(ispin,ikp)->c().nloc();
            
            // next line: add enhancement factor to descent direction
            // since there is no instability of the Anderson iteration
            // This improves convergence in most cases
            
            //const double psda_enhancement_factor = 2.0;
            for ( int n = 0; n < nloc; n++ ) {
              // note: double mloc length for complex<double> indices
              double* dcn = &dc[2*mloc*n];
              // loop to ngwl only since diag[i] is not defined on [0:mloc-1]
              for ( int i = 0; i < ngwl; i++ ) {
                const double fac = diag[i];
                //const double f0 = -psda_enhancement_factor*fac * dcn[2*i];
                //const double f1 = -psda_enhancement_factor*fac * dcn[2*i+1];
                const double f0 = -fac * dcn[2*i];
                const double f1 = -fac * dcn[2*i+1];
                dcn[2*i] = f0;
                dcn[2*i+1] = f1;
              }
            }
            tmap_["psda_prec"].stop();
          
            // dwf now contains the preconditioned descent
            // direction -K(HV-VA)
            tmap_["psda_update_wf"].start();
            // Anderson extrapolation
            if ( extrapolate_[ispin][ikp] ) {
              double theta = 0.0;
              double a = 0.0, b = 0.0;

              //ewd OMP HERE?
              for ( int lj=0; lj < wf_.sd(ispin,ikp)->c().nblocks(); lj++ )
              {
                 for ( int jj=0; jj < wf_.sd(ispin,ikp)->c().nbs(lj); jj++ )
                 {
                    // global state index
                    const int nglobal = wf_.sd(ispin,ikp)->c().j(lj,jj);
                    const int norig = lj*wf_.sd(ispin,ikp)->c().nb()+jj;
                    const vector<double>& occ = wf_.sd(ispin,ikp)->occ();
                    const double occn = occ[nglobal];


                    //ewd DEBUG
                    //if (wf_.sd(ispin,ikp)->context().myrow() == 0 )
                    //{
                    //   int mype = wf_.sd(ispin,ikp)->context().mype();
                    //   cout << "PSDA, mype = " << mype << ", lj = " << lj << ", jj = " << jj << ", nglobal = " << nglobal << ", occ = " << occ[nglobal] << endl;
                    //}
                    //ewd DEBUG

                    

                    //ewd DEBUG
                    //const double occn = 1.0;

                    for ( int i = 0; i < 2*mloc; i++ ) {
                       const double f = dc[i +2*mloc*norig];
                       const double delta_f = f - dc_last[i + 2*mloc*norig];

                       //ewd:  weight sums by occupation to improve convergence in metallic systems
                       a += occn * f * delta_f;
                       b += occn * delta_f * delta_f;
                    }
                 }
              }
              // a and b contain the partial sums of a and b
              double tmpvec[2] = { a, b };
              int kpc = wf_.mysdctxt(ikp);
              assert(kpc >= 0);
              wf_.sdcontext(ispin,kpc)->dsum(2,1,&tmpvec[0],1);
              //ewdwf_.sd(ispin,ikp)->context().dsum(2,1,&tmpvec[0],1);
              a = tmpvec[0];
              b = tmpvec[1];
              // compute theta = - a / b
              if ( b != 0.0 )
                theta = - a / b;
            
              if ( ikp == 0 && ispin == 0 && wf_.sd(ispin,ikp)->context().oncoutpe() )
                cout << "  <!-- Anderson extrapolation: theta = " << theta;
            
              if ( theta < -1.0 )
                theta = 0.0;
            
              theta = min(2.0,theta);
              //ewd try tightening restriction here...
              //theta = min(1.0,theta);
            
              if ( ikp == 0 && ispin == 0 && wf_.sd(ispin,ikp)->context().oncoutpe() )
                cout <<" ( using " << theta << " )" << " -->"<< endl;
            
              // extrapolation
              for ( int i = 0; i < 2*mloc*nloc; i++ ) {
                // x_bar = x_ + theta * ( x_ - xlast_ ) (store in x_)
                const double x = c[i];
                const double xlast = c_last[i];
                const double xbar = x + theta * ( x - xlast );
                
                // f_bar = f + theta * ( f - flast ) (store in f)
                const double f = dc[i];
                const double flast = dc_last[i];
                const double fbar = f + theta * ( f - flast );
              
                c[i] = xbar + fbar;
                c_last[i] = x;
                dc_last[i] = f;
              }
            }
            else {
              // no extrapolation
              for ( int i = 0; i < 2*mloc*nloc; i++ ) {
                // x_ = x_ + f_
                const double x = c[i];
                const double f = dc[i];
                    
                c[i] = x + f;
                c_last[i] = x;
                dc_last[i] = f;
              }
            }
            extrapolate_[ispin][ikp] = true;
            tmap_["psda_update_wf"].stop();

            /*
            // do orthogonalization outside of this routine, so ultrasoft functions
            // can update   
            
            enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };
            //const ortho_type ortho = GRAM;
            //const ortho_type ortho = LOWDIN;
            const ortho_type ortho = ORTHO_ALIGN;
          
            switch ( ortho ) {
              case GRAM:
                  tmap_["gram"].start();
                  wf_.sd(ispin,ikp)->gram();
                  tmap_["gram"].stop();
                  break;
                
              case LOWDIN:
                  tmap_["lowdin"].start();
                  wf_.sd(ispin,ikp)->lowdin();
                  tmap_["lowdin"].stop();
                  break;
                
              case ORTHO_ALIGN:
                  tmap_["ortho_align"].start();
                  wf_.sd(ispin,ikp)->ortho_align(*wf_last_.sd(ispin,ikp));
                  tmap_["ortho_align"].stop();
                  break;
                  
              case RICCATI:
                  tmap_["riccati"].start();
                  wf_.sd(ispin,ikp)->riccati(*wf_last_.sd(ispin,ikp));
                  tmap_["riccati"].stop();
                  break;
            }

            */
            
          }
        }
      }
    }
  }
}
