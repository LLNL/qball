
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

extern "C" {
//  void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int*);
  void dsygvd_(int*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, int*, int*);
  double ddot_(const int *, const double *, const int *, const double *, const int *);     
}


////////////////////////////////////////////////////////////////////////////////
DQBPCGWavefunctionStepper::DQBPCGWavefunctionStepper(Wavefunction& wf,
                                                     Preconditioner& p,const AtomSet& atoms,
                                                     const ChargeDensity& cd_,vector<vector<double> >& v_r,
                                                     TimerMap& tmap) :
    WavefunctionStepper(wf,tmap), prec_(p), hpsi_(wf,atoms,cd_,v_r), reswf_(wf), hreswf_(wf)
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

for (int rr_step = 0; rr_step < 1; rr_step++){

                  const int maxit = 5;  // set max number of eigensolver iterations 

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


                  // Create single column matrices. 
                  // These are auxiliary objects needed to perform distributed linear algebra 
                  // operations on columns (e.g., scalar products) 
                  DoubleMatrix wfcol(wf_proxy.context(),wf_proxy.m(), 1, wf_proxy.mb(),1);
                  DoubleMatrix rescol(wf_proxy.context(),wf_proxy.m(), 1, wf_proxy.mb(),1);
                  DoubleMatrix pcol(wf_proxy.context(),wf_proxy.m(), 1, wf_proxy.mb(),1);
                  DoubleMatrix hwfcol(wf_proxy.context(),wf_proxy.m(), 1, wf_proxy.mb(),1);
                  DoubleMatrix hrescol(wf_proxy.context(),wf_proxy.m(), 1, wf_proxy.mb(),1);
                  DoubleMatrix hpcol(wf_proxy.context(),wf_proxy.m(), 1, wf_proxy.mb(),1);

                  // ``cross-product'' matrices
                  DoubleMatrix a(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
                                 wf_proxy.nb(),wf_proxy.nb());
                  DoubleMatrix g(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
                                 wf_proxy.nb(),wf_proxy.nb());

                  const int mloc = wf_.sd(ispin,ikp)->c().mloc();
                  const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
                  const int nloc = wf_.sd(ispin,ikp)->c().nloc();
    
                  const int me = wf_proxy.context().myproc(); // eugene: temporary, for printing
                  cout << std::scientific;


if (me==0){
cout << " m = " << wf_proxy.m() << " n = " << wf_proxy.n() << "mb = " << wf_proxy.mb() << "nb = " << wf_proxy.nb() << endl;
cout << " mloc = " << mloc << " ngwl = " << ngwl << " nloc = " << nloc << endl;
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
            
           //       double nrm = res_proxy.nrm2();
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

// eugene: complete the initial iteration


                  // the 3-by-3 projected matrices  
                  double kmat[9], mmat[9];   

                  // parameters for lapack calls 
                  char jobz = 'V', uplo = 'U';
                  double* work;
                  double  ev[3];
                  int itype = 1, dimp, ld, lwork, liwork, info = 0; 
                  int* iwork;

                  // initialize parameters for separate eigensolves  
                  dimp = 2;  // dimension of the projected problem
                  ld   = 3;  // leading dimension for kmat and mmat  
                  lwork  = 1 + 6*ld + 2*ld*ld;
                  liwork = 3 + 5*ld;
                  work   = new double[max(1,lwork)]; 
                  iwork  = new int[max(1,liwork)]; 

                  // Get pointers to values of single column matrices
                  double* wfcol_coeff   = (double*) wfcol.valptr();
                  double* rescol_coeff  = (double*) rescol.valptr();
                  double* pcol_coeff    = (double*) pcol.valptr();
                  double* hwfcol_coeff  = (double*) hwfcol.valptr();
                  double* hrescol_coeff = (double*) hrescol.valptr();
                  double* hpcol_coeff   = (double*) hpcol.valptr();
          
 
                  // initialize the ``conjugate direction'' with the residual
                  p  = res_proxy;
                  hp = hres_proxy; 

                  // update wavefunctions by constructing and solving 2x2 problems
                  for ( int j = 0; j < nloc; j++ ) {

                     // get pointer to j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy
                     double* wf_coeff   =  (double*) wf_proxy.valptr(2*mloc*j);
                     double* hwf_coeff  =  (double*) hwf_proxy.valptr(2*mloc*j);
                     double* res_coeff  =  (double*) res_proxy.valptr(2*mloc*j);
                     double* hres_coeff =  (double*) hres_proxy.valptr(2*mloc*j);

                     // get pointer to j-th column of p and hp 
                     double* p_coeff   =  (double*) p.valptr(2*mloc*j);
                     double* hp_coeff  =  (double*) hp.valptr(2*mloc*j);
                       

                     // copy j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy to 
                     // wfcol, hwfcol, rescol, and hrescol, respectively
                     for ( int i = 0; i < ngwl; i++ ) {

                        wfcol_coeff[2*i]    =  wf_coeff[2*i];
                        wfcol_coeff[2*i+1]  =  wf_coeff[2*i+1];
             
                        hwfcol_coeff[2*i]   =  hwf_coeff[2*i];
                        hwfcol_coeff[2*i+1] =  hwf_coeff[2*i+1];
             
                        rescol_coeff[2*i]   =  res_coeff[2*i];
                        rescol_coeff[2*i+1] =  res_coeff[2*i+1];
            
                        hrescol_coeff[2*i]   =  hres_coeff[2*i];
                        hrescol_coeff[2*i+1] =  hres_coeff[2*i+1];

                     }

//double dd = wfcol.dot(rescol);
//if (me == 0){cout << "dot = " << 2*dd - rescol_coeff[0]*wfcol_coeff[0] << endl;}


                     // form the 2 by 2 lhs projected matrix kmat
                    
                     double mat_loc [9] = {0} ;
                     int ngwl2 = 2*ngwl;
                     int ione = 1;

                     mat_loc[0]=2.0*ddot_(&ngwl2, wfcol_coeff, &ione, hwfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[0] -= wfcol_coeff[0]*hwfcol_coeff[0];
                     }
                     
                     mat_loc[1]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, hwfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[1] -= rescol_coeff[0]*hwfcol_coeff[0];
                     }

                     mat_loc[4]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, hrescol_coeff, &ione);
                     if (me == 0){
                        mat_loc[4] -= rescol_coeff[0]*hrescol_coeff[0];
                     }

                     MPI_Allreduce(mat_loc, kmat, 9, MPI_DOUBLE, MPI_SUM, wf_proxy.context().comm() ); 

                     kmat[3] = kmat[1];

                     mat_loc[0]=2.0*ddot_(&ngwl2, wfcol_coeff, &ione, wfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[0] -= wfcol_coeff[0]*wfcol_coeff[0];
                     }
                     
                     mat_loc[1]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, wfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[1] -= rescol_coeff[0]*wfcol_coeff[0];
                     }

                     mat_loc[4]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, rescol_coeff, &ione);
                     if (me == 0){
                        mat_loc[4] -= rescol_coeff[0]*rescol_coeff[0];
                     }

                     MPI_Allreduce(mat_loc, mmat, 9, MPI_DOUBLE, MPI_SUM, wf_proxy.context().comm() ); 
  
                     mmat[3] = mmat[1];
                    
/*
                     kmat[0]  = 2.0*wfcol.dot(hwfcol);
                     kmat[0] -= wfcol_coeff[0]*hwfcol_coeff[0];      /// !!!! <---- 
  
                     kmat[1]  = 2.0*rescol.dot(hwfcol);
                     kmat[1] -= rescol_coeff[0]*hwfcol_coeff[0];   
                     kmat[3]  = kmat[1];
                     kmat[4]  = 2.0*rescol.dot(hrescol);
                     kmat[4] -= rescol_coeff[0]*hrescol_coeff[0];   


                     // form the 2 by 2 rhs projected matrix mmat
                     mmat[0]  = 2.0*wfcol.dot(wfcol);
                     mmat[0] -= wfcol_coeff[0]*wfcol_coeff[0];   
                     mmat[1]  = 2.0*rescol.dot(wfcol);
                     mmat[1] -= rescol_coeff[0]*wfcol_coeff[0];   
                     mmat[3]  = mmat[1];
                     mmat[4]  = 2.0*rescol.dot(rescol);
                     mmat[4] -= rescol_coeff[0]*rescol_coeff[0];   
*/
                     // Solve the 2-by-2 problem
                     dsygvd_(&itype, &jobz, &uplo, &dimp, kmat, &ld, mmat, &ld, ev, work, &lwork, iwork, &liwork, &info);
                     assert ( info == 0 );

                     // update the conjugate direction and the wavefunction
                     for ( int i = 0; i < ngwl; i++ ) {
                        p_coeff[2*i]  *= kmat[1]; 
                        hp_coeff[2*i] *= kmat[1]; 
                        p_coeff[2*i+1]  *= kmat[1]; 
                        hp_coeff[2*i+1] *= kmat[1]; 
                              
                        wf_coeff[2*i]  = wf_coeff[2*i]*kmat[0] + p_coeff[2*i];
                        hwf_coeff[2*i] = hwf_coeff[2*i]*kmat[0] + hp_coeff[2*i];
                        wf_coeff[2*i+1]  = wf_coeff[2*i+1]*kmat[0] + p_coeff[2*i+1];
                        hwf_coeff[2*i+1] = hwf_coeff[2*i+1]*kmat[0] + hp_coeff[2*i+1];
                     }
                    

                 } // end for j

//cout << " m = " << a.m() << " n = " << a.n() << "mb = " << a.mb() << "nb = " << a.nb() << endl;
//cout << " mloc = " << mloc << " ngwl = " << ngwl << " nloc = " << nloc << endl;
//cout << " m = " << wf_proxy.m() << " n = " << wf_proxy.n() << endl;
//                 double* acoef = (double*) a.valptr();
//                       for (int ii = 0; ii < 64; ii++){
//                            cout << ii << endl;
//                            cout << acoef[ii] << endl;
//                       }

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
               
                 // Start the main loop
                 for (int iter=1; iter<maxit; iter++){   // loop condition to be changed

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
                     hpsi_.compute(reswf_,hreswf_);

  
                     // orthogonalize conjugate direction  against wavefunction
                     g.gemm('t','n',2.0,wf_proxy,p,0.0);
                     g.ger(-1.0,wf_proxy,0,p,0);
                     p.gemm('n','n',-1.0,wf_proxy,g,1.0);
                     hp.gemm('n','n',-1.0,hwf_proxy,g,1.0);

                     // inner loop
                     // update wavefunctions by constructing and solving 2x2 problems
                     for ( int j = 0; j < nloc; j++ ) {

                        // get pointer to j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy
                        double* wf_coeff   =  (double*) wf_proxy.valptr(2*mloc*j);
                        double* hwf_coeff  =  (double*) hwf_proxy.valptr(2*mloc*j);
                        double* res_coeff  =  (double*) res_proxy.valptr(2*mloc*j);
                        double* hres_coeff =  (double*) hres_proxy.valptr(2*mloc*j);

                        // get pointer to j-th column of p and hp 
                        double* p_coeff   =  (double*) p.valptr(2*mloc*j);
                        double* hp_coeff  =  (double*) hp.valptr(2*mloc*j);
                       

                        // copy j-th column of blocks to the single column matrices 
                        for ( int i = 0; i < ngwl; i++ ) {

                           wfcol_coeff[2*i]    =  wf_coeff[2*i];
                           wfcol_coeff[2*i+1]  =  wf_coeff[2*i+1];
             
                           hwfcol_coeff[2*i]   =  hwf_coeff[2*i];
                           hwfcol_coeff[2*i+1] =  hwf_coeff[2*i+1];
             
                           rescol_coeff[2*i]   =  res_coeff[2*i];
                           rescol_coeff[2*i+1] =  res_coeff[2*i+1];
            
                           hrescol_coeff[2*i]   =  hres_coeff[2*i];
                           hrescol_coeff[2*i+1] =  hres_coeff[2*i+1];

                           pcol_coeff[2*i]   =  p_coeff[2*i];
                           pcol_coeff[2*i+1] =  p_coeff[2*i+1];
            
                           hpcol_coeff[2*i]   =  hp_coeff[2*i];
                           hpcol_coeff[2*i+1] =  hp_coeff[2*i+1];
                           
                        }

/*
//// normalize separate res and p columns 
                       //double res_nrm = rescol.nrm2();
//if (me == 0){
//cout << hrescol_coeff[0] << endl;
//cout << hrescol_coeff[1] << endl;
//}

                       double res_nrm = rescol.dot(rescol);
                       res_nrm = sqrt( 2.0*res_nrm - rescol_coeff[0]*rescol_coeff[0] ) ; 
                       rescol *=  1/res_nrm; 
                       hrescol *= 1/res_nrm; 


                      
       //             res_nrm = rescol.nrm2();
       //             res_nrm = sqrt( 2.0*res_nrm*res_nrm - rescol_coeff[0]*rescol_coeff[0] ) ; 
       //             if (me == 0){cout << res_nrm << endl;}                   
 
                        //double p_nrm = rescol.nrm2();
                       // p_nrm = sqrt( 2.0*p_nrm*p_nrm - pcol_coeff[0]*pcol_coeff[0] ) ; 
                      // finish normalization of p vectors....

                       // Copy the normalized column back to res_proxy and hres_proxy 
                       for ( int i = 0; i < ngwl; i++ ) {
             
                          res_coeff[2*i] = rescol_coeff[2*i];
                          res_coeff[2*i+1] = rescol_coeff[2*i+1];  
                          hres_coeff[2*i] =  hrescol_coeff[2*i];
                          hres_coeff[2*i+1] = hrescol_coeff[2*i+1];

            //              pcol_coeff[2*i]   =  p_coeff[2*i];
            //              pcol_coeff[2*i+1] =  p_coeff[2*i+1];
            
            //              hpcol_coeff[2*i]   =  hp_coeff[2*i];
            //              hpcol_coeff[2*i+1] =  hp_coeff[2*i+1];
                           
                       }
//// end normalize separate res and p columns 
*/

//double dd = wfcol.dot(rescol);
//double ddp = wfcol.dot(pcol);
//if (me == 0){cout << "dot = " << 2*dd - rescol_coeff[0]*wfcol_coeff[0] << endl;}
//if (me == 0){cout << "pdot = " << 2*ddp - pcol_coeff[0]*wfcol_coeff[0] << endl;}
                       
                     double mat_loc [9] = {0} ;
                     int ngwl2 = 2*ngwl;
                     int ione = 1;

                     mat_loc[0]=2.0*ddot_(&ngwl2, wfcol_coeff, &ione, hwfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[0] -= wfcol_coeff[0]*hwfcol_coeff[0];
                     }
                     
                     mat_loc[1]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, hwfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[1] -= rescol_coeff[0]*hwfcol_coeff[0];
                     }

                     mat_loc[2]=2.0*ddot_(&ngwl2, pcol_coeff, &ione, hwfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[2] -= pcol_coeff[0]*hwfcol_coeff[0];
                     }



                     mat_loc[4]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, hrescol_coeff, &ione);
                     if (me == 0){
                        mat_loc[4] -= rescol_coeff[0]*hrescol_coeff[0];
                     }

                     mat_loc[5]=2.0*ddot_(&ngwl2, pcol_coeff, &ione, hrescol_coeff, &ione);
                     if (me == 0){
                        mat_loc[5] -= pcol_coeff[0]*hrescol_coeff[0];
                     }

                     mat_loc[8]=2.0*ddot_(&ngwl2, pcol_coeff, &ione, hpcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[5] -= pcol_coeff[0]*hpcol_coeff[0];
                     }


                     MPI_Allreduce(mat_loc, kmat, 9, MPI_DOUBLE, MPI_SUM, wf_proxy.context().comm() ); 
                     kmat[3] = kmat[1];
                     kmat[6]  = kmat[2];
                     kmat[7]  = kmat[5];


                     mat_loc[0]=2.0*ddot_(&ngwl2, wfcol_coeff, &ione, wfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[0] -= wfcol_coeff[0]*wfcol_coeff[0];
                     }
                     
                     mat_loc[1]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, wfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[1] -= rescol_coeff[0]*wfcol_coeff[0];
                     }


                     mat_loc[2]=2.0*ddot_(&ngwl2, pcol_coeff, &ione, wfcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[2] -= pcol_coeff[0]*wfcol_coeff[0];
                     }


                     mat_loc[4]=2.0*ddot_(&ngwl2, rescol_coeff, &ione, rescol_coeff, &ione);
                     if (me == 0){
                        mat_loc[4] -= rescol_coeff[0]*rescol_coeff[0];
                     }


                     mat_loc[5]=2.0*ddot_(&ngwl2, pcol_coeff, &ione, rescol_coeff, &ione);
                     if (me == 0){
                        mat_loc[5] -= pcol_coeff[0]*rescol_coeff[0];
                     }

                     mat_loc[8]=2.0*ddot_(&ngwl2, pcol_coeff, &ione, pcol_coeff, &ione);
                     if (me == 0){
                        mat_loc[5] -= pcol_coeff[0]*pcol_coeff[0];
                     }



                     MPI_Allreduce(mat_loc, mmat, 9, MPI_DOUBLE, MPI_SUM, wf_proxy.context().comm() ); 
                     mmat[3] = mmat[1];
                     mmat[6]  = mmat[2];
                     mmat[7]  = mmat[5];


/*
                       // form the 3 by 3 lhs projected matrix kmat
                       kmat[0]  = 2.0*wfcol.dot(hwfcol);
                       kmat[0] -= wfcol_coeff[0]*hwfcol_coeff[0];   
                       kmat[1]  = 2.0*rescol.dot(hwfcol);
                       kmat[1] -= rescol_coeff[0]*hwfcol_coeff[0];   
                       kmat[2]  = 2.0*pcol.dot(hwfcol);
                       kmat[2] -= pcol_coeff[0]*hwfcol_coeff[0]; 

                       kmat[3]  = kmat[1];
                       kmat[4]  = 2.0*rescol.dot(hrescol);
                       kmat[4] -= rescol_coeff[0]*hrescol_coeff[0];   
                       kmat[5]  = 2.0*pcol.dot(hrescol);
                       kmat[5] -= pcol_coeff[0]*hrescol_coeff[0]; 

                       kmat[6]  = kmat[2];
                       kmat[7]  = kmat[5];
                       kmat[8]  = 2.0*pcol.dot(hpcol);
                       kmat[8] -= pcol_coeff[0]*hpcol_coeff[0]; 


                       // form the 3 by 3 rhs projected matrix mmat
                       mmat[0]  = 2.0*wfcol.dot(wfcol);
                       mmat[0] -= wfcol_coeff[0]*wfcol_coeff[0];   
                       mmat[1]  = 2.0*rescol.dot(wfcol);
                       mmat[1] -= rescol_coeff[0]*wfcol_coeff[0];   
                       mmat[2]  = 2.0*pcol.dot(wfcol);
                       mmat[2] -= pcol_coeff[0]*wfcol_coeff[0]; 


                       mmat[3]  = mmat[1];
                       mmat[4]  = 2.0*rescol.dot(rescol);
                       mmat[4] -= rescol_coeff[0]*rescol_coeff[0];   
                       mmat[5]  = 2.0*pcol.dot(rescol);
                       mmat[5] -= pcol_coeff[0]*rescol_coeff[0]; 


                       mmat[6]  = mmat[2];
                       mmat[7]  = mmat[5];
                       mmat[8]  = 2.0*pcol.dot(pcol);
                       mmat[8] -= pcol_coeff[0]*pcol_coeff[0]; 
*/

                       // Solve the 3-by-3 problem
                       dimp = 3;
                       dsygvd_(&itype, &jobz, &uplo, &dimp, kmat, &ld, mmat, &ld, ev, work, &lwork, iwork, &liwork, &info);
                       assert ( info == 0 );
//if (me==0){cout << ev[0] << " " << ev[1] << " " << ev[2] << endl;}
                       // update the conjugate direction and the wavefunction
                       for ( int i = 0; i < ngwl; i++ ) {
                          p_coeff[2*i]  =   kmat[1]*res_coeff[2*i] +  kmat[2]*p_coeff[2*i]; 
                          hp_coeff[2*i] =   kmat[1]*hres_coeff[2*i] +  kmat[2]*hp_coeff[2*i]; 
                          p_coeff[2*i+1]  =  kmat[1]*res_coeff[2*i+1] +  kmat[2]*p_coeff[2*i+1]; 
                          hp_coeff[2*i+1] =  kmat[1]*hres_coeff[2*i+1] + kmat[2]*hp_coeff[2*i+1]; 
                              
                          wf_coeff[2*i]  = wf_coeff[2*i]*kmat[0] + p_coeff[2*i];
                          hwf_coeff[2*i] = hwf_coeff[2*i]*kmat[0] + hp_coeff[2*i];
                          wf_coeff[2*i+1]  = wf_coeff[2*i+1]*kmat[0] + p_coeff[2*i+1];
                          hwf_coeff[2*i+1] = hwf_coeff[2*i+1]*kmat[0] + hp_coeff[2*i+1];
                        }
                    

                    } // end for j

                    // Orthogonalize wavefunctions
                    a.gemm('t','n', 2.0, wf_proxy, wf_proxy, 0.0);
                    a.ger(-1.0, wf_proxy,0, wf_proxy,0);
                    a.potrf('U');
                    wf_proxy.trsm('r', 'U', 'n', 'n', 1.0, a);
                    hwf_proxy.trsm('r', 'U', 'n', 'n', 1.0, a);
///begin????
//                    p.trsm('r', 'U', 'n', 'n', 1.0, a);
//                    hp.trsm('r', 'U', 'n', 'n', 1.0, a);
///end????

                    // Compute new residual
                    a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
                    a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                    res_proxy = hwf_proxy;
                    res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);

//                    nrm = res_proxy.nrm2();
                    nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                    if (me == 0){ cout << "Iteration = " << iter+1 <<  ". Residual norm est. = " <<  nrm  << endl; }


                 } //end main loop



                 tmap_["dqbpcg_iterloop"].stop();
                  

                 tmap_["dqbpcg_rr"].start();
//   hpsi_.compute(wf_,hwf);
 
                 // The Rayleigh-Ritz procedure
                 a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
                 a.ger(-1.0,wf_proxy,0,hwf_proxy,0);

                
                 g.gemm('t','n',2.0, wf_proxy,wf_proxy,0.0);
                 g.ger(-1.0,wf_proxy,0,wf_proxy,0);

                 // Transform to standard eigenvalue problem
                 g.potrf('U');
                 a.sygst(1,'U', g);
                  

                 // Eigenvectors and eigenvalues of the reduced problem
                 DoubleMatrix z(wf_proxy.context(),wf_proxy.n(),wf_proxy.n(),
                                wf_proxy.nb(),wf_proxy.nb());
                 valarray<double> evalues(z.m());

                 a.syevd('U', evalues, z);
                 z.trsm('l', 'U', 'n', 'n', 1.0, g);

                 p.gemm('n','n',1.0,wf_proxy, z,0.0);
                 hp.gemm('n','n',1.0,hwf_proxy, z,0.0);
                 wf_proxy = p;
                 hwf_proxy = hp;

///////
                    // Compute new residual
//                 a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
//                 a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
//         
//                 res_proxy = hwf_proxy;
//                 res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);
// New code begin
                 for ( int j = 0; j < nloc; j++ ) {

                    // get pointer to j-th column of wf_proxy, hwf_proxy, res_proxy, and hres_proxy
                    double* wf_coeff   =  (double*) wf_proxy.valptr(2*mloc*j);
                    double* hwf_coeff  =  (double*) hwf_proxy.valptr(2*mloc*j);
                    double* res_coeff  =  (double*) res_proxy.valptr(2*mloc*j);

                       // Copy the normalized column back to res_proxy and hres_proxy 
                       for ( int i = 0; i < ngwl; i++ ) {
             
                          res_coeff[2*i] =   hwf_coeff[2*i] - wf_coeff[2*i]*evalues[j];
                          res_coeff[2*i+1] =  hwf_coeff[2*i+1] - wf_coeff[2*i+1]*evalues[j];
                           
                       }
                 }
// New code end
                



//                    nrm = res_proxy.nrm2();
//                    if (me == 0){ cout << "Final Residual norm = " <<  nrm  << endl; }
//
////  --- Print out individual residuals
//
                   double nrms [8];
                   for (int jj = 0; jj<8; jj++){

                     double* res_coeff  =  (double*) res_proxy.valptr(2*mloc*jj);
                     for ( int i = 0; i < ngwl; i++ ) {
                        rescol_coeff[2*i]   =  res_coeff[2*i];
                        rescol_coeff[2*i+1] =  res_coeff[2*i+1];
                     }
                     nrms[jj] = rescol.nrm2();
                   } 

                  if (me == 0){
                     for (int jj = 0; jj < 8; jj++){ 
                       cout << "Eigenvalue " << jj+1 << "= " <<  evalues[jj] << ". Res. = " << nrms[jj] << endl;
                     }
                  }

///////
 
                tmap_["dqbpcg_rr"].stop();


                
                delete [] work;
                delete [] iwork;
if (me == 0) {cout << "rr_step = " << rr_step << endl;}

} /// end RR STEP






                  
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
}
