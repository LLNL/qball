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
#include "blas.h"
#include <iostream>
using namespace std;

//extern "C" {
//  void dsygvd(int*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, int*, int*);
//  double ddot(const int *, const double *, const int *, const double *, const int *);     
//}


////////////////////////////////////////////////////////////////////////////////
DQBPCGWavefunctionStepper::DQBPCGWavefunctionStepper(Wavefunction& wf,
                                                     Preconditioner& p, const int maxit, const AtomSet& atoms,
                                                     const ChargeDensity& cd_,vector<vector<double> >& v_r,
                                                     TimerMap& tmap) :
    WavefunctionStepper(wf,tmap), prec_(p), hpsi_(wf,atoms,cd_,v_r), reswf_(wf), hreswf_(wf),maxit_(maxit)
{


  nkp_ = wf_.nkp();
  nspin_ = wf_.nspin();


 // eugene add
 cell_moved();

}

////////////////////////////////////////////////////////////////////////////////
void DQBPCGWavefunctionStepper::cell_moved()
{
   hpsi_.cell_moved(wf_);
}
////////////////////////////////////////////////////////////////////////////////
void DQBPCGWavefunctionStepper::update(Wavefunction& hwf)
{ 
   for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
      if (wf_.spinactive(ispin)) {
         for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
            if (wf_.kptactive(ikp)) {
               assert(wf_.sd(ispin,ikp) != 0);
               
               // compute A = V^T H V  and descent direction HV - VA
               if ( wf_.sd(ispin,ikp)->basis().real() ) {

//                  for (int rr_step = 0; rr_step < 1; rr_step++)
//                  {

                     // proxy real matrices for wavefunctions and residuals
                     DoubleMatrix wf_proxy(wf_.sd(ispin,ikp)->c());
                     DoubleMatrix hwf_proxy(hwf.sd(ispin,ikp)->c());
                     DoubleMatrix res_proxy(reswf_.sd(ispin,ikp)->c());
                     DoubleMatrix hres_proxy(hreswf_.sd(ispin,ikp)->c());
                     
                     // the ``conjugate direction''                   
                     DoubleMatrix p(wf_proxy.context(),wf_proxy.m(),wf_proxy.n(),
                                    wf_proxy.mb(),wf_proxy.nb());
                     DoubleMatrix hp(wf_proxy.context(),wf_proxy.m(),wf_proxy.n(),
                                     wf_proxy.mb(),wf_proxy.nb());

                     // ``cross-product'' matrices
                     DoubleMatrix a(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
                                    wf_proxy.nb(),wf_proxy.nb());
                     DoubleMatrix g(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
                                    wf_proxy.nb(),wf_proxy.nb());

                     const int mloc = wf_.sd(ispin,ikp)->c().mloc();
                     const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
                     const int nloc = wf_.sd(ispin,ikp)->c().nloc();
                     
                     // store communication data
                     const Context& ctxt = wf_proxy.context();
                     const int me = wf_proxy.context().myproc(); // eugene: temporary, for printing
                     cout << std::scientific;


                     if (me==0){
                        cout << " m = " << wf_proxy.m() << " n = " << wf_proxy.n() << "mb = " << wf_proxy.mb() << "nb = " << wf_proxy.nb() << endl;
                        cout << " mloc = " << mloc << " ngwl = " << ngwl << " nloc = " << nloc << endl;
                        cout << " nprow = " << ctxt.nprow() << " npcol = " << ctxt.npcol() << endl;
                     }

                     tmap_["dqbpcg_residual"].start();

                     // factor 2.0 in next line: G and -G
                     a.gemm('t','n',2.0,wf_proxy,hwf_proxy,0.0);
                     // rank-1 update correction
                     a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                     // reswf = hwf - wf * a
                     res_proxy = hwf_proxy;
                     res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);
        
                     // reswf now contains the descent direction (HV-VA)
                     double nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                     if (me == 0){ cout << "Init. residual norm est. = " << nrm   << endl; }
                  
                     tmap_["dqbpcg_residual"].stop();


                     // apply preconditioner
                     tmap_["dqbpcg_prec"].start();
                     const valarray<double>& diag = prec_.diag(ispin,ikp);
                     double* wc = (double*) res_proxy.valptr();
            
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

                     // orthogonalize preconditioned residual against wavefunction
                     g.gemm('t','n',2.0,wf_proxy,res_proxy,0.0);
                     g.ger(-1.0,wf_proxy,0,res_proxy,0);
                     res_proxy.gemm('n','n',-1.0,wf_proxy,g,1.0);
 
                     // compute H*residual
                     hpsi_.compute(reswf_,hreswf_);

                     // the 3-by-3 projected matrices  
                     double kmat[9], mmat[9];   

                     // parameters for lapack calls 
                     char jobz = 'V', uplo = 'U';
                     double  ev[3] = {0};
                     int itype = 1; 
                     int info  = 0;

                     // initialize parameters for separate eigensolves  
                     int dimp = 2;  // dimension of the projected problem
                     int ld   = 3;  // leading dimension for kmat and mmat  
                     int lwork  = 1 + 6*ld + 2*ld*ld;
                     int liwork = 3 + 5*ld;
                     double* work   = new double[max(1,lwork)]; 
                     int* iwork  = new int[max(1,liwork)]; 

                     // update wavefunctions by constructing and solving 2x2 problems
                     for ( int j = 0; j < nloc; j++ ) {

                        // get pointer to j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy
                        double* wfcol   =  (double*) wf_proxy.valptr(2*mloc*j);
                        double* hwfcol  =  (double*) hwf_proxy.valptr(2*mloc*j);
                        double* rescol  =  (double*) res_proxy.valptr(2*mloc*j);
                        double* hrescol =  (double*) hres_proxy.valptr(2*mloc*j);
                        double* pcol  =  (double*) p.valptr(2*mloc*j);
                        double* hpcol =  (double*) hp.valptr(2*mloc*j);

                        // form the 2 by 2 lhs projected matrix kmat
                        int ngwl2 = 2*ngwl;
                        int ione = 1;

                        kmat[0]=2.0*ddot(&ngwl2, wfcol, &ione, hwfcol, &ione);
                        if (wf_proxy.context().myrow() == 0){
                           kmat[0] -= wfcol[0]*hwfcol[0];
                        }
                     
                        kmat[1]=2.0*ddot(&ngwl2, rescol, &ione, hwfcol, &ione);
                        if (wf_proxy.context().myrow() == 0){
                           kmat[1] -= rescol[0]*hwfcol[0];
                        }

                        kmat[4]=2.0*ddot(&ngwl2, rescol, &ione, hrescol, &ione);
                        if (wf_proxy.context().myrow() == 0){
                           kmat[4] -= rescol[0]*hrescol[0];
                        }

                        kmat[3] = kmat[1];

                        ctxt.dsum('C', 2, 2, kmat, ld);  
 
                        // form the 2 by 2 rhs projected matrix mmat
                        mmat[0]=2.0*ddot(&ngwl2, wfcol, &ione, wfcol, &ione);
                        if (wf_proxy.context().myrow() == 0){
                           mmat[0] -= wfcol[0]*wfcol[0];
                        }
                     
                        mmat[1]=2.0*ddot(&ngwl2, rescol, &ione, wfcol, &ione);
                        if (wf_proxy.context().myrow() == 0){
                           mmat[1] -= rescol[0]*wfcol[0];
                        }

                        mmat[4]=2.0*ddot(&ngwl2, rescol, &ione, rescol, &ione);
                        if (wf_proxy.context().myrow() == 0){
                           mmat[4] -= rescol[0]*rescol[0];
                        }

                        mmat[3] = mmat[1];

                        ctxt.dsum('C', 2, 2, mmat, ld);  
                     
                        // Solve the 2-by-2 problem
                        dsygvd(&itype, &jobz, &uplo, &dimp, kmat, &ld, mmat, &ld, ev, work, &lwork, iwork, &liwork, &info);
                        assert ( info == 0 );

                        // update the conjugate direction and the wavefunction
                        for ( int i = 0; i < ngwl; i++ ) {
                           pcol[2*i]    = rescol[2*i]*kmat[1]; 
                           pcol[2*i+1]  = rescol[2*i+1]*kmat[1]; 
                           hpcol[2*i]   = hrescol[2*i]*kmat[1]; 
                           hpcol[2*i+1] = hrescol[2*i+1]*kmat[1]; 
                              
                           wfcol[2*i]    = wfcol[2*i]*kmat[0] + pcol[2*i];
                           wfcol[2*i+1]  = wfcol[2*i+1]*kmat[0] + pcol[2*i+1];
                           hwfcol[2*i]   = hwfcol[2*i]*kmat[0] + hpcol[2*i];
                           hwfcol[2*i+1] = hwfcol[2*i+1]*kmat[0] + hpcol[2*i+1];
                        }
                    

                     } // end for j


                     // Now orthonormalize wavefunctions using Cholesky based QR
                     a.gemm('t','n', 2.0, wf_proxy, wf_proxy, 0.0);
                     a.ger(-1.0, wf_proxy,0, wf_proxy,0);
                     a.potrf('U');
                     wf_proxy.trsm('r', 'U', 'n', 'n', 1.0, a);
                     hwf_proxy.trsm('r', 'U', 'n', 'n', 1.0, a);

                     // Compute new residual
                     a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
                     a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                     res_proxy = hwf_proxy;
                     res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);

                     //      nrm = res_proxy.nrm2();
                     nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                     if (me == 0){ cout << "Iteration = 1." << " Residual norm est. = " <<  nrm  << endl; }
 
                     tmap_["dqbpcg_iter0"].stop();


                     tmap_["dqbpcg_iterloop"].start();
               
//// !----             // Start the main loop
                     for (int iter=1; iter<maxit_; iter++){   // loop condition to be changed

                        // apply preconditioner to hwf                  
                        tmap_["dqbpcg_prec"].start();
                        const valarray<double>& diag = prec_.diag(ispin,ikp);
                        double* wc = (double*) res_proxy.valptr();
            
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

                        // orthogonalize preconditioned residual against wavefunction
                        g.gemm('t','n',2.0,wf_proxy,res_proxy,0.0);
                        g.ger(-1.0,wf_proxy,0,res_proxy,0);
                        res_proxy.gemm('n','n',-1.0,wf_proxy,g,1.0);
                  
                        // compute H*residual
                        tmap_["dqbpcg_hpsi"].start();
                        hpsi_.compute(reswf_,hreswf_);
                        tmap_["dqbpcg_hpsi"].stop();

  
                        // orthogonalize conjugate direction  against wavefunction
                        g.gemm('t','n',2.0,wf_proxy,p,0.0);
                        g.ger(-1.0,wf_proxy,0,p,0);
                        p.gemm('n','n',-1.0,wf_proxy,g,1.0);
                        hp.gemm('n','n',-1.0,hwf_proxy,g,1.0);

                        // inner loop
                        // update wavefunctions by constructing and solving 3x3 eigenproblems
                        for ( int j = 0; j < nloc; j++ ) {

                           // get pointer to j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy
                           double* wfcol   =  (double*) wf_proxy.valptr(2*mloc*j);
                           double* hwfcol  =  (double*) hwf_proxy.valptr(2*mloc*j);
                           double* rescol  =  (double*) res_proxy.valptr(2*mloc*j);
                           double* hrescol =  (double*) hres_proxy.valptr(2*mloc*j);
                           double* pcol  =  (double*) p.valptr(2*mloc*j);
                           double* hpcol =  (double*) hp.valptr(2*mloc*j);
         
                          // EV: Normalizing rescol and pcol can be added here if needed (haven't observe any benefit)
                     
                           int ngwl2 = 2*ngwl;
                           int ione = 1;

                           // form the 3 by 3 lhs projected matrix kmat
                           kmat[0]=2.0*ddot(&ngwl2, wfcol, &ione, hwfcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              kmat[0] -= wfcol[0]*hwfcol[0];
                           }
                     
                           kmat[1]=2.0*ddot(&ngwl2, rescol, &ione, hwfcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              kmat[1] -= rescol[0]*hwfcol[0];
                           }

                           kmat[2]=2.0*ddot(&ngwl2, pcol, &ione, hwfcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              kmat[2] -= pcol[0]*hwfcol[0];
                           }

                          
                           kmat[4]=2.0*ddot(&ngwl2, rescol, &ione, hrescol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              kmat[4] -= rescol[0]*hrescol[0];
                           }

                           kmat[5]=2.0*ddot(&ngwl2, pcol, &ione, hrescol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              kmat[5] -= pcol[0]*hrescol[0];
                           }

                           kmat[8]=2.0*ddot(&ngwl2, pcol, &ione, hpcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              kmat[8] -= pcol[0]*hpcol[0];
                           }
                       
                           kmat[3] = kmat[1];
                           kmat[6] = kmat[2];
                           kmat[7] = kmat[5];

                           ctxt.dsum('C', 3, 3, kmat, ld);  

                          // form the 3 by 3 rhs projected matrix mmat
                           mmat[0]=2.0*ddot(&ngwl2, wfcol, &ione, wfcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              mmat[0] -= wfcol[0]*wfcol[0];
                           }
                     
                           mmat[1]=2.0*ddot(&ngwl2, rescol, &ione, wfcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              mmat[1] -= rescol[0]*wfcol[0];
                           }

                           mmat[2]=2.0*ddot(&ngwl2, pcol, &ione, wfcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              mmat[2] -= pcol[0]*wfcol[0];
                           }

                           mmat[4]=2.0*ddot(&ngwl2, rescol, &ione, rescol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              mmat[4] -= rescol[0]*rescol[0];
                           }

                           mmat[5]=2.0*ddot(&ngwl2, pcol, &ione, rescol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              mmat[5] -= pcol[0]*rescol[0];
                           }

                           mmat[8]=2.0*ddot(&ngwl2, pcol, &ione, pcol, &ione);
                           if (wf_proxy.context().myrow() == 0){
                              mmat[8] -= pcol[0]*pcol[0];
                           }

                           mmat[3] = mmat[1];
                           mmat[6] = mmat[2];
                           mmat[7] = mmat[5];

                           ctxt.dsum('C', 3, 3, mmat, ld);  


//if (me == 0){
//   cout << "kmat: " << endl;
//   cout << kmat[0] << " " << kmat[3] << endl;
//   cout << kmat[1] << " " << kmat[4] << endl;
//   cout << "mmat: " << endl;
//   cout << mmat[0] << " " << mmat[3] << endl;
//   cout << mmat[1] << " " << mmat[4] << endl;
//}




                           // Solve the 3-by-3 problem
                           dimp = 3;
                           dsygvd(&itype, &jobz, &uplo, &dimp, kmat, &ld, mmat, &ld, ev, work, &lwork, iwork, &liwork, &info);
                           assert ( info == 0 );
//if (me == 0){cout << "ev = " << ev[0] << "    " << ev[1] << endl;}

                           // update the conjugate direction and the wavefunction
                           for ( int i = 0; i < ngwl; i++ ) {
                              
                              pcol[2*i]    =  kmat[1]*rescol[2*i]  +  kmat[2]*pcol[2*i]; 
                              pcol[2*i+1]  =  kmat[1]*rescol[2*i+1] +  kmat[2]*pcol[2*i+1]; 
                              hpcol[2*i]   =  kmat[1]*hrescol[2*i] +  kmat[2]*hpcol[2*i]; 
                              hpcol[2*i+1] =  kmat[1]*hrescol[2*i+1] + kmat[2]*hpcol[2*i+1];
 
                              wfcol[2*i]    = wfcol[2*i]*kmat[0] + pcol[2*i];
                              wfcol[2*i+1]  = wfcol[2*i+1]*kmat[0] + pcol[2*i+1];
                              hwfcol[2*i]   = hwfcol[2*i]*kmat[0] + hpcol[2*i];
                              hwfcol[2*i+1] = hwfcol[2*i+1]*kmat[0] + hpcol[2*i+1];

                           }

                        } // end for j

                        // Orthogonalize wavefunctions
                        tmap_["dqbpcg_ortho"].start();
                        a.gemm('t','n', 2.0, wf_proxy, wf_proxy, 0.0);
                        a.ger(-1.0, wf_proxy,0, wf_proxy,0);
                        a.potrf('U');
                        wf_proxy.trsm('r', 'U', 'n', 'n', 1.0, a);
                        hwf_proxy.trsm('r', 'U', 'n', 'n', 1.0, a);
                        tmap_["dqbpcg_ortho"].stop();

                        // Compute new residual
                        a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
                        a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                        res_proxy = hwf_proxy;
                        res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);

                        // Print out Frob. norm estimate at this iteration
                        nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                        if (me == 0){ cout << "Iteration = " << iter+1 <<  ". Residual norm est. = " <<  nrm  << endl; }


                     } //end main loop


                     //ewd DEBUG:  calculate eigenvalues
                     if (false)
                     {
                        DoubleMatrix h(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
                                       wf_proxy.nb(),wf_proxy.nb());
                        valarray<double> w(h.m());
                        h.gemm('t','n',2.0,wf_proxy,hwf_proxy,0.0);
                        h.ger(-1.0,wf_proxy,0,hwf_proxy,0);
                        h.syevd('l',w);
                        if (me == 0)
                        {
                           int neig = wf_proxy.n();
                           const double eVolt = 2.0 * 13.6056923;
                           const int neig_line = 8;
                           cout <<    "DEBUG, wfstepper eigenvalues" << endl;
                           for ( int i = 0; i < neig; i++ ) {
                              cout << w[i]*eVolt;
                              if ( i%neig_line == neig_line-1 ) cout << endl;
                           }
                           if ( neig%neig_line != 0 ) cout << endl;
                           
                        }
                     }
                     //ewd DEBUG
                     
                     tmap_["dqbpcg_iterloop"].stop();
                  

//------------- The RR procedure comment begin !!!!! ------------------------------------------
//                     tmap_["dqbpcg_rr"].start();
//
//
//
// 
//                     // The Rayleigh-Ritz procedure
//                     a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
//                     a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
//
//                
////ev                     g.gemm('t','n',2.0, wf_proxy,wf_proxy,0.0);
////ev                     g.ger(-1.0,wf_proxy,0,wf_proxy,0);
////ev
////ev                     // Transform to standard eigenvalue problem
////ev                     g.potrf('U');
////ev                     a.sygst(1,'U', g);
//                  
//
//                     // Eigenvectors and eigenvalues of the reduced problem
//                     DoubleMatrix z(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
//                                    wf_proxy.nb(),wf_proxy.nb());
//                     valarray<double> evalues(z.m());
//
//                     a.syevd('U', evalues, z);
////ev                     z.trsm('l', 'U', 'n', 'n', 1.0, g);
//
////ev                     p.gemm('n','n',1.0,wf_proxy, z,0.0);
////ev                     hp.gemm('n','n',1.0,hwf_proxy, z,0.0);
////ev                     wf_proxy = p;
////ev                     hwf_proxy = hp;
//// use res_proxy and hres_proxy as temporary variables
//                     res_proxy.gemm('n','n',1.0,wf_proxy, z,0.0);
//                     hres_proxy.gemm('n','n',1.0,hwf_proxy, z,0.0);
//                     wf_proxy = res_proxy;
//                     hwf_proxy = hres_proxy;
//
//
//
//
//
////// check orthogonality
///*
//                     g.gemm('t','n',2.0,wf_proxy,res_proxy,0.0);
//                     g.ger(-1.0,wf_proxy,0,res_proxy,0);
//                     double* gc = (double*) g.valptr();
//                     if (me == 0){
//                        for (int i = 0; i<8; i++){ 
//                          for (int j = 0; j<8;j++){
//                             cout << gc[8*i+j] << " ";
//                          }
//                          cout << endl; 
//                        }
//                       cout << endl; 
//                     }
//
//                     g.gemm('t','n',2.0,wf_proxy,wf_proxy,0.0);
//                     g.ger(-1.0,wf_proxy,0,wf_proxy,0);
//                     double* gc = (double*) g.valptr();
//                     if (me == 0){
//                        for (int i = 0; i<8; i++){ 
//                          for (int j = 0; j<8;j++){
//                             cout << gc[8*i+j] << " ";
//                          }
//                          cout << endl; 
//                        }
//                     }
//*/
///////////// end check                      



/////////////
//const double eigenvalue_sum = wf_.dot(hwf);
//if (me == 0){
//   double es = 0;
//   for (int j = 0; j<8; j++){
//     es += evalues[j]; 
//   }
//cout << "Evalues sum = " << es << endl;
//cout << "Evalues sum check = " << eigenvalue_sum << endl;
//
//}
//
////////////

//                     ///////
//                     // Compute new residual
//                     //                 a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
//                     //                 a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
//                     //         
//                     //                 res_proxy = hwf_proxy;
//                     //                 res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);
//                     // New code begin
//                     for ( int j = 0; j < nloc; j++ ) {
//
//                        // get pointer to j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy
//                        double* wf_coeff   =  (double*) wf_proxy.valptr(2*mloc*j);
//                        double* hwf_coeff  =  (double*) hwf_proxy.valptr(2*mloc*j);
//                        double* res_coeff  =  (double*) res_proxy.valptr(2*mloc*j);
//
//                        // Copy the normalized column back to res_proxy and hres_proxy 
//                        for ( int i = 0; i < ngwl; i++ ) {
//             
//                           res_coeff[2*i] =   hwf_coeff[2*i] - wf_coeff[2*i]*evalues[j];
//                           res_coeff[2*i+1] =  hwf_coeff[2*i+1] - wf_coeff[2*i+1]*evalues[j];
//                           
//                        }
//                     }
//                     // New code end
//                
//
//
//
//                     //                    nrm = res_proxy.nrm2();
//                     //                    if (me == 0){ cout << "Final Residual norm = " <<  nrm  << endl; }
//                     //
//                     ////  --- Print out individual residuals
//                     //
//                     int ngwl2 = 2*ngwl;
//                     int ione = 1;
//                     double loc_nrms [8];
//                     double nrms [8];
//                     for (int jj = 0; jj<8; jj++){
//                        double* rescol  =  (double*) res_proxy.valptr(2*mloc*jj);
//                        loc_nrms[jj]=2.0*ddot(&ngwl2, rescol, &ione, rescol, &ione);
////                        if (me == 0){
//                        if (wf_proxy.context().myrow() == 0){
//                           loc_nrms[jj] -= rescol[0]*rescol[0];
//                        }
//                     }
//
//                    MPI_Allreduce(loc_nrms, nrms, 8, MPI_DOUBLE, MPI_SUM, wf_proxy.context().comm() ); 
//                      
// 
//                     if (me == 0){
//                        for (int jj = 0; jj < 8; jj++){ 
//                           cout << "Eigenvalue " << jj+1 << "= " <<  evalues[jj] << ". Res. = " << nrms[jj] << endl;
//                        }
//                     }
//
//
//             
//
//
//                     /////// iend Print out --- needs to be redone
// 
 //                    tmap_["dqbpcg_rr"].stop();
//

                
                     delete [] work;
                     delete [] iwork;
//                     if (me == 0) {cout << "rr_step = " << rr_step << endl;}

//                  } /// end RR STEP

                  
               }
               else {
                  //ev begin comment
                  /*
                    ComplexMatrix &c_pr xy = wf_.sd(ispin,ikp)->c();
                    ComplexMatrix &cp = hwf.sd(ispin,ikp)->c();
                    ComplexMatrix &res_proxy = reswf_.sd(ispin,ikp)->c();
                    ComplexMatrix &hw = hreswf_.sd(ispin,ikp)->c();
                    ComplexMatrix a(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),wf_proxy.nb(),wf_proxy.nb());
                    ComplexMatrix g(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),wf_proxy.nb(),wf_proxy.nb());

                    const int mloc = wf_.sd(ispin,ikp)->c().mloc();
                    const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
                    const int nloc = wf_.sd(ispin,ikp)->c().nloc();


                    tmap_["dqbpcg_residual"].start();
                    a.gemm('c','n',1.0,wf_proxy,cp,0.0);
                    // w = cp - c * a
                    res_proxy = cp;
                    res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);
                    // hwf.sd->c() now contains the descent direction (HV-VA)
                    tmap_["dqbpcg_residual"].stop();
             
                    // apply preconditioner to hwf             
                    tmap_["dqbpcg_prec"].start();
                    const valarray<double>& diag = prec_.diag(ispin,ikp);
                    double* wc = (double*) res_proxy.valptr();
                    //                  const int mloc = res_proxy.mloc();
                    //                  const int nloc = res_proxy.nloc();
                    //                  const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
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
                    // hwf -> (hwf - psi^t psi hwf)
                    g.gemm('c','n',1.0,wf_proxy,cp,0.0);
                    res_proxy.gemm('n','n',-1.0,wf_proxy,g,1.0);
                  
                    // compute H*residual
                    hpsi_.compute(reswf_,hreswf_);





                  
                    tmap_["dqbpcg_iter0"].stop();
                    tmap_["dqbpcg_iterloop"].start();



                    tmap_["dqbpcg_iterloop"].stop();
                  
                    //ev end comment
                    */           
               }
               // eugene: Run the inner loop


            }
         }
      }
   }

   //ewd DEBUG
   //wf_.diag(hwf,true);
   //wf_.printeig();

   const int me = wf_.context().mype();
   const double eigsum = wf_.dot(hwf);
   const double selfsum = wf_.dot(wf_);
   const double hwfsum = hwf.dot(hwf);

   if (me == 0)
      cout << "DEBUG SUMS:  " << eigsum << " " << selfsum << " " << hwfsum << endl;
   //ewd DEBUG
   
}
