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
// PPCGWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "PPCGWavefunctionStepper.h"
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
PPCGWavefunctionStepper::PPCGWavefunctionStepper(Wavefunction& wf,
                                                     Preconditioner& p, const int maxit, const AtomSet& atoms,
                                                     const ChargeDensity& cd_,vector<vector<double> >& v_r,
                                                     TimerMap& tmap, const int sbsize, const int qrstep) :
    WavefunctionStepper(wf,tmap), prec_(p), hpsi_(wf,atoms,cd_,v_r), reswf_(wf), hreswf_(wf),maxit_(maxit),
    sbsize_(sbsize), qrstep_(qrstep)               
{


  nkp_ = wf_.nkp();
  nspin_ = wf_.nspin();


 cell_moved();

 assert(sbsize_ > 0);
 assert(qrstep_ > 0);

}

////////////////////////////////////////////////////////////////////////////////
void PPCGWavefunctionStepper::cell_moved()
{
   hpsi_.cell_moved(wf_);
}
////////////////////////////////////////////////////////////////////////////////
void PPCGWavefunctionStepper::update(Wavefunction& hwf)
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
                     int verbose = 0;   // if verbose == 0 then no output
                                        // verbose shoulbe be set to 0 for tests
                                        // ortherwise there is overhead from computing residual norms 

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


                     if (me==0 && verbose){
                        cout << " m = " << wf_proxy.m() << " n = " << wf_proxy.n() << "mb = " << wf_proxy.mb() << "nb = " << wf_proxy.nb() << endl;
                        cout << " mloc = " << mloc << " ngwl = " << ngwl << " nloc = " << nloc << endl;
                        cout << " nprow = " << ctxt.nprow() << " npcol = " << ctxt.npcol() << endl;
                        cout << " sbsize = " << sbsize_ << " qrstep = " << qrstep_ << endl;
                     }

                     tmap_["ppcg_residual"].start();

                     // factor 2.0 in next line: G and -G
                     a.gemm('t','n',2.0,wf_proxy,hwf_proxy,0.0);
                     // rank-1 update correction
                     a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                     // reswf = hwf - wf * a
                     res_proxy = hwf_proxy;
                     res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);
        
                     // reswf now contains the descent direction (HV-VA)
                     if (verbose){
                        double nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                        if (me == 0){ cout << "Init. residual norm est. = " << nrm   << endl; }
                     }                  

                     tmap_["ppcg_residual"].stop();


                     // apply preconditioner
                     tmap_["ppcg_prec"].start();
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
                     tmap_["ppcg_prec"].stop();

                     tmap_["ppcg_iter0"].start();

                     // orthogonalize preconditioned residual against wavefunction
                     g.gemm('t','n',2.0,wf_proxy,res_proxy,0.0);
                     g.ger(-1.0,wf_proxy,0,res_proxy,0);
                     res_proxy.gemm('n','n',-1.0,wf_proxy,g,1.0);
 
                     // compute H*residual
                     tmap_["ppcg_hpsi"].start();
                     hpsi_.compute(reswf_,hreswf_);
                     tmap_["ppcg_hpsi"].stop();

                     // update wavefunctions by constructing and solving 2x2 problems
                     tmap_["ppcg_update"].start();
                     update_wfp0(wf_proxy, hwf_proxy, res_proxy, hres_proxy, p, hp, ngwl, mloc, nloc);
                     tmap_["ppcg_update"].stop();

                     // Now orthonormalize wavefunctions using Cholesky based QR
                     if ( qrstep_ == 1 ) { 
                        tmap_["ppcg_ortho"].start();
                        cholqr(wf_proxy, hwf_proxy);
                        tmap_["ppcg_ortho"].stop();
                     }

                     // Compute new residual
                     a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
                     a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                     res_proxy = hwf_proxy;
                     res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);

                     //      nrm = res_proxy.nrm2();
                     if ( verbose ){
                        double nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                        if (me == 0){ cout << "Iteration = 1." << " Residual norm est. = " <<  nrm  << endl; }
                     }
                     tmap_["ppcg_iter0"].stop();

                     tmap_["ppcg_iterloop"].start();
                     // Start the main loop
                     for (int iter=1; iter<maxit_; iter++){   // loop condition to be changed

                        // apply preconditioner to hwf                  
                        tmap_["ppcg_prec"].start();
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
                        tmap_["ppcg_prec"].stop();

                        // orthogonalize preconditioned residual against wavefunction
                        g.gemm('t','n',2.0,wf_proxy,res_proxy,0.0);
                        g.ger(-1.0,wf_proxy,0,res_proxy,0);
                        res_proxy.gemm('n','n',-1.0,wf_proxy,g,1.0);
                  
                        // compute H*residual
                        tmap_["ppcg_hpsi"].start();
                        hpsi_.compute(reswf_,hreswf_);
                        tmap_["ppcg_hpsi"].stop();

  
                        // orthogonalize conjugate direction  against wavefunction
                        g.gemm('t','n',2.0,wf_proxy,p,0.0);
                        g.ger(-1.0,wf_proxy,0,p,0);
                        p.gemm('n','n',-1.0,wf_proxy,g,1.0);
                        hp.gemm('n','n',-1.0,hwf_proxy,g,1.0);

                        // update wavefunctions by constructing and solving 3x3 eigenproblems
                        tmap_["ppcg_update"].start();
                        update_wfp(wf_proxy, hwf_proxy, res_proxy, hres_proxy, p, hp, ngwl, mloc, nloc);
                        tmap_["ppcg_update"].stop();

                        // Orthogonalize wavefunctions
                        if ( (iter+1)%qrstep_ == 0 ) { 
                           tmap_["ppcg_ortho"].start();
                           cholqr(wf_proxy, hwf_proxy);
                           tmap_["ppcg_ortho"].stop();
                        }

                        // Compute new residual
                        a.gemm('t','n',2.0, wf_proxy,hwf_proxy,0.0);
                        a.ger(-1.0,wf_proxy,0,hwf_proxy,0);
         
                        res_proxy = hwf_proxy;
                        res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);

                        // Print out Frob. norm estimate at this iteration
                        if (verbose) {
                           double nrm = sqrt(2.0*res_proxy.dot(res_proxy));
                           if (me == 0){ cout << "Iteration = " << iter+1 <<  ". Residual norm est. = " <<  nrm  << endl; }
                        }

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
                     
                     tmap_["ppcg_iterloop"].stop();
                  

//------------- The RR procedure comment begin !!!!! ------------------------------------------
//                     tmap_["ppcg_rr"].start();
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
 //                    tmap_["ppcg_rr"].stop();
//

                
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


                    tmap_["ppcg_residual"].start();
                    a.gemm('c','n',1.0,wf_proxy,cp,0.0);
                    // w = cp - c * a
                    res_proxy = cp;
                    res_proxy.gemm('n','n',-1.0,wf_proxy,a,1.0);
                    // hwf.sd->c() now contains the descent direction (HV-VA)
                    tmap_["ppcg_residual"].stop();
             
                    // apply preconditioner to hwf             
                    tmap_["ppcg_prec"].start();
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
                    tmap_["ppcg_prec"].stop();


                    tmap_["ppcg_iter0"].start();                  
                    // hwf -> (hwf - psi^t psi hwf)
                    g.gemm('c','n',1.0,wf_proxy,cp,0.0);
                    res_proxy.gemm('n','n',-1.0,wf_proxy,g,1.0);
                  
                    // compute H*residual
                    hpsi_.compute(reswf_,hreswf_);





                  
                    tmap_["ppcg_iter0"].stop();
                    tmap_["ppcg_iterloop"].start();



                    tmap_["ppcg_iterloop"].stop();
                  
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




////////////////////////////////////////////////////////////////////////////////
// Cholesky decomposition based QR factorization
// Orthonormalize wf and update hwf accordingly
//
void PPCGWavefunctionStepper::cholqr(DoubleMatrix& wf, DoubleMatrix& hwf) 
{
   DoubleMatrix a(wf.context(),wf.n(),wf.n(),wf.nb(),wf.nb());

   a.gemm('t','n', 2.0, wf, wf, 0.0);
   a.ger(-1.0, wf,0, wf,0);
   a.potrf('U');
   wf.trsm('r', 'U', 'n', 'n', 1.0, a);
   hwf.trsm('r', 'U', 'n', 'n', 1.0, a);
}



////////////////////////////////////////////////////////////////////////////////
// update wf, hwf, p, and hp
// inner loop of the PPCG algorithm
//
void PPCGWavefunctionStepper::update_wfp(DoubleMatrix& wf, DoubleMatrix& hwf, 
                                           DoubleMatrix& res,  DoubleMatrix& hres, 
                                           DoubleMatrix&  p, DoubleMatrix& hp, int ngwl, int mloc, int nloc)
{

//   const int me = wf.context().myproc(); // eugene: temporary, for printing
   assert(sbsize_ > 0);
   int nsb  =  nloc / sbsize_ ;          // number of subblocks 
   int sbsize_last = nloc % sbsize_;      // last sublbock contains the remaining columns

   if ( sbsize_last ){
      nsb++;
   }
   else{  
      sbsize_last = sbsize_;
   }

   int* sbsizes = new int[nsb];             // sizes of each subblock  
   for ( int j = 0; j < nsb - 1; j++ ) {
      sbsizes[j] = sbsize_;
   }
   sbsizes[nsb-1] = sbsize_last;

   int ld   = 3*sbsize_;                    // leading dimension of the projected matrices
   int ldsq =  ld*ld;                     
   double* kmat  = new double[ldsq];       // ld-by-ld projected matrices kmat and mmat 
   double* mmat  = new double[ldsq];   

   double* kmat_store  = new double[ldsq*nsb];    // ld-by-(ld*nsb) storage for nsb projected matrices
   double* mmat_store  = new double[ldsq*nsb];
  
   int  ngwl2  = 2*ngwl;
   int  mloc2 = 2*mloc; 

   double dzero = 0, dmone = -1, done = 1, dtwo = 2;
   char trans = 'T', notrans = 'N';

   const Context& ctxt = wf.context();
   int my_proc_row = ctxt.myrow();
 
   // construct and strore nsb small eigenproblems 
   for ( int j = 0; j < nsb ; j ++ ) {

      // get pointers to jth subblocks
      int     bstart  = mloc2*j*sbsize_; 
      double* bwf   =  (double*) wf.valptr(bstart);
      double* bhwf  =  (double*) hwf.valptr(bstart);
      double* bres  =  (double*) res.valptr(bstart);
      double* bhres =  (double*) hres.valptr(bstart);
      double* bp    =  (double*) p.valptr(bstart);
      double* bhp   =  (double*) hp.valptr(bstart);
  
      // EV: Normalizing rescol and pcol can be added here if needed (haven't observe any benefit)

      int  sb = sbsizes[j];

      // Form the 3sbize by 3sbize lhs projected matrix kmat

      // bwf'H*bwf      
      double* dest = kmat; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bwf, &mloc2, bhwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bwf, &mloc2, bhwf, &mloc2, dest, &ld);    
      }

      // bres'H*bwf     
      dest = kmat + sb; 
      dgemm(&trans, &notrans, &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bhwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bhwf, &mloc2, dest, &ld);    
      }

      // bp'H*bwf     
      dest = kmat + 2*sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bp, &mloc2, bhwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bp, &mloc2, bhwf, &mloc2, dest, &ld);    
      }

      // bres'H*bres     
      dest = kmat + sb*ld + sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bhres, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bhres, &mloc2, dest, &ld);    
      }

      // bp'H*bres     
      dest = kmat + sb*ld + 2*sb;
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bp, &mloc2, bhres, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bp, &mloc2, bhres, &mloc2, dest, &ld);    
      }

      // bp'H*bp     
      dest = kmat + 2*sb*ld + 2*sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bp, &mloc2, bhp, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bp, &mloc2, bhp, &mloc2, dest, &ld);    
      }

      // copy the remaining blocks of kmat by symmetry 
      double* sb_src  [3] = { kmat + sb, kmat + 2*sb, kmat + sb*ld + 2*sb };        // pointers to "from" (source) blocks 
      double* sb_dest [3] = { kmat + sb*ld, kmat + 2*sb*ld, kmat + 2*sb*ld + sb }; // pointers to "to" (destination) blocks 
      
      for (int s = 0; s < 3; s++) {
         double* src  = sb_src[s];
         double* dest = sb_dest[s];
         for (int l = 0; l < sbsizes[j]; l++ ){
            for (int i = 0; i < sbsizes[j]; i++ ){
               dest[i*ld] = src[i];
            }
            src  += ld;
            dest += 1;
         } 
      }

      // Form the 3sbize by 3sbize rhs projected matrix mmat

      // bwf'bwf      
      dest = mmat;
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bwf, &mloc2, bwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bwf, &mloc2, bwf, &mloc2, dest, &ld);    
      } 

     // bres'bwf      
      dest = mmat + sb;
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bwf, &mloc2, dest, &ld);    
      }

      // bp'bwf      
      dest = mmat + 2*sb;
      dgemm(&trans, &notrans,  &sb, &sb, &ngwl2, &dtwo, bp, &mloc2, bwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bp, &mloc2, bwf, &mloc2, dest, &ld);    
      }

      // bres'bres     
      dest = mmat + sb*ld + sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bres, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bres, &mloc2, dest, &ld);    
      }

      // bp'bres     
      dest = mmat + sb*ld + 2*sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bp, &mloc2, bres, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bp, &mloc2, bres, &mloc2, dest, &ld);    
      }

      // bp'bp     
      dest = mmat + 2*sb*ld + 2*sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bp, &mloc2, bp, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bp, &mloc2, bp, &mloc2, dest, &ld);    
      }   


      // copy the remaining blocks of mmat by symmetry
      // pointers to 'from' blocks
      sb_src[0]  = mmat + sb;
      sb_src[1]  = mmat + 2*sb;
      sb_src[2]  = mmat + sb*ld + 2*sb;
      // pointers to 'to' blocks
      sb_dest[0] = mmat + sb*ld;
      sb_dest[1] = mmat + 2*sb*ld;
      sb_dest[2] = mmat + 2*sb*ld + sb;

      for (int s = 0; s < 3; s++) {
         double* src  = sb_src[s];
         double* dest = sb_dest[s];
         for (int l = 0; l < sbsizes[j]; l++ ){
            for (int i = 0; i < sbsizes[j]; i++ ){
               dest[i*ld] = src[i];
            }
            src  += ld;
            dest += 1;
         } 
      }

      // store kmat and mmat
      memcpy(kmat_store + j*ldsq, kmat, ldsq*sizeof(double)); 
      memcpy(mmat_store + j*ldsq, mmat, ldsq*sizeof(double)); 

   }
  

   // sum the local parts of kmat and mmat across all processors in the column 
   ctxt.dsum('C', 1, ldsq*nsb, kmat_store, 1);  
   ctxt.dsum('C', 1, ldsq*nsb, mmat_store, 1);  
   
   double* ev     = new double[ld];            // eigenvalues of a projected problems   
   double* buffer = new double[mloc2*sbsize_]; 
   memset(buffer, 0, sizeof(double)*mloc2*sbsize_);   // set the buffer to zero (important!)

   // parameters for blas and lapack calls 
   char jobz  = 'V', uplo = 'U';
   int itype = 1; 
   int info  = 0;
   int lwork  = 1 + 6*ld + 2*ldsq;
   int liwork = 3 + 5*ld;
   double* work   = new double[max(1,lwork)]; 
   int* iwork  = new int[max(1,liwork)]; 


   // solve nsb projected problems and update subblocks

   for ( int j = 0; j < nsb ; j++ ) {

      // get pointers to jth subblocks
      int     bstart  = mloc2*j*sbsize_; 
      double* bwf   =  (double*) wf.valptr(bstart);
      double* bhwf  =  (double*) hwf.valptr(bstart);
      double* bres  =  (double*) res.valptr(bstart);
      double* bhres =  (double*) hres.valptr(bstart);
      double* bp    =  (double*) p.valptr(bstart);
      double* bhp   =  (double*) hp.valptr(bstart);

      int  sb = sbsizes[j];

      memcpy(kmat, kmat_store + j*ldsq, ldsq*sizeof(double)); 
      memcpy(mmat, mmat_store + j*ldsq, ldsq*sizeof(double)); 

      int dimp  = 3*sb;    // size of the current projected problem 
      dsygvd(&itype, &jobz, &uplo, &dimp, kmat, &ld, mmat, &ld, ev, work, &lwork, iwork, &liwork, &info);
      assert ( info == 0 );
  
      // update bwf and bp    
      size_t copy_size = sizeof(double)*mloc2*sb; 
      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bp, &mloc2, kmat + 2*sb, &ld, &dzero, buffer,  &mloc2); 
      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bres, &mloc2, kmat + sb, &ld, &done, buffer,  &mloc2); 
      memcpy(bp, buffer, copy_size); 

      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bwf, &mloc2, kmat, &ld, &done, buffer,  &mloc2); 
      memcpy(bwf, buffer, copy_size); 
     
      // update bhwf and bhp     
      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bhp, &mloc2, kmat + 2*sb, &ld, &dzero, buffer,  &mloc2); 
      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bhres, &mloc2, kmat + sb, &ld, &done, buffer,  &mloc2); 
      memcpy(bhp, buffer, copy_size);  

      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bhwf, &mloc2, kmat, &ld, &done, buffer,  &mloc2); 
      memcpy(bhwf, buffer, copy_size); 

   } // end solve nsb projected problems and update subblocks

   delete [] kmat;
   delete [] kmat_store;
   delete [] mmat;   
   delete [] mmat_store;   
   delete [] buffer;
   delete [] ev;
   delete [] work;
   delete [] iwork;


}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// update wf, hwf, p, and hp
// inner loop of *initial iteration* of the PPCG algorithm
// this code is separated (duplicated) from the update_wfp to avoid
// repeating if statements in the main PPCG loop 
// 
void PPCGWavefunctionStepper::update_wfp0(DoubleMatrix& wf, DoubleMatrix& hwf, 
                                           DoubleMatrix& res,  DoubleMatrix& hres, 
                                           DoubleMatrix&  p, DoubleMatrix& hp, int ngwl, int mloc, int nloc)
{

//   const int me = wf.context().myproc(); // eugene: temporary, for printing
   assert(sbsize_ > 0);
   int nsb  =  nloc / sbsize_ ;          // number of subblocks 
   int sbsize_last = nloc % sbsize_;      // last sublbock contains the remaining columns

   if ( sbsize_last ){
      nsb++;
   }
   else{  
      sbsize_last = sbsize_;
   }

   int* sbsizes = new int[nsb];             // sizes of each subblock  
   for ( int j = 0; j < nsb - 1; j++ ) {
      sbsizes[j] = sbsize_;
   }
   sbsizes[nsb-1] = sbsize_last;

   int ld   = 2*sbsize_;                    // leading dimension of the projected matrices
   int ldsq =  ld*ld;                     
   double* kmat  = new double[ldsq];       // ld-by-ld projected matrices kmat and mmat 
   double* mmat  = new double[ldsq];   

   double* kmat_store  = new double[ldsq*nsb];    // ld-by-(ld*nsb) storage for nsb projected matrices
   double* mmat_store  = new double[ldsq*nsb];
  
   int  ngwl2  = 2*ngwl;
   int  mloc2 = 2*mloc; 

   double dzero = 0, dmone = -1, done = 1, dtwo = 2;
   char trans = 'T', notrans = 'N';

   const Context& ctxt = wf.context();
   int my_proc_row = ctxt.myrow();
 
   // construct and strore nsb small eigenproblems 
   for ( int j = 0; j < nsb ; j ++ ) {

      // get pointers to jth subblocks
      int     bstart  = mloc2*j*sbsize_; 
      double* bwf   =  (double*) wf.valptr(bstart);
      double* bhwf  =  (double*) hwf.valptr(bstart);
      double* bres  =  (double*) res.valptr(bstart);
      double* bhres =  (double*) hres.valptr(bstart);
  
      // EV: Normalizing rescol and pcol can be added here if needed (haven't observe any benefit)

      int  sb = sbsizes[j];

      // Form the 2sbize by 2sbize lhs projected matrix kmat

      // bwf'H*bwf      
      double* dest = kmat; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bwf, &mloc2, bhwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bwf, &mloc2, bhwf, &mloc2, dest, &ld);    
      }

      // bres'H*bwf     
      dest = kmat + sb; 
      dgemm(&trans, &notrans, &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bhwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bhwf, &mloc2, dest, &ld);    
      }

      // bres'H*bres     
      dest = kmat + sb*ld + sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bhres, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bhres, &mloc2, dest, &ld);    
      }

      // copy the remaining blocks of kmat by symmetry 
      double* src  = kmat + sb;
      dest         = kmat + sb*ld;
      for (int l = 0; l < sbsizes[j]; l++ ){
         for (int i = 0; i < sbsizes[j]; i++ ){
            dest[i*ld] = src[i];
         }
         src  += ld;
         dest += 1;
      } 

      // Form the 2sbize by 2sbize rhs projected matrix mmat

      // bwf'bwf      
      dest = mmat;
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bwf, &mloc2, bwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bwf, &mloc2, bwf, &mloc2, dest, &ld);    
      } 

     // bres'bwf      
      dest = mmat + sb;
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bwf, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bwf, &mloc2, dest, &ld);    
      }

      // bres'bres     
      dest = mmat + sb*ld + sb; 
      dgemm(&trans, &notrans,  &sb, &sb,  &ngwl2, &dtwo, bres, &mloc2, bres, &mloc2, &dzero, dest,  &ld); 
      if ( my_proc_row == 0 ) {
         dger(&sb, &sb, &dmone, bres, &mloc2, bres, &mloc2, dest, &ld);    
      }

      // copy the remaining blocks of mmat by symmetry
      src   = mmat + sb;
      dest  = mmat + sb*ld;
      for (int l = 0; l < sbsizes[j]; l++ ){
         for (int i = 0; i < sbsizes[j]; i++ ){
            dest[i*ld] = src[i];
         }
         src  += ld;
         dest += 1;
      } 

      // store kmat and mmat
      memcpy(kmat_store + j*ldsq, kmat, ldsq*sizeof(double)); 
      memcpy(mmat_store + j*ldsq, mmat, ldsq*sizeof(double)); 

   }

   // sum the local parts of kmat and mmat across all processors in the column 
   ctxt.dsum('C', 1, ldsq*nsb, kmat_store, 1);  
   ctxt.dsum('C', 1, ldsq*nsb, mmat_store, 1);  
   
   double* ev     = new double[ld];            // eigenvalues of a projected problems   
   double* buffer = new double[mloc2*sbsize_]; 
   memset(buffer, 0, sizeof(double)*mloc2*sbsize_);   // set the buffer to zero (important!)

   // parameters for blas and lapack calls 
   char jobz  = 'V', uplo = 'U';
   int itype = 1; 
   int info  = 0;
   int lwork  = 1 + 6*ld + 2*ldsq;
   int liwork = 3 + 5*ld;
   double* work   = new double[max(1,lwork)]; 
   int* iwork  = new int[max(1,liwork)]; 


   // solve nsb projected problems and update subblocks

   for ( int j = 0; j < nsb ; j++ ) {

      // get pointers to jth subblocks
      int     bstart  = mloc2*j*sbsize_; 
      double* bwf   =  (double*) wf.valptr(bstart);
      double* bhwf  =  (double*) hwf.valptr(bstart);
      double* bres  =  (double*) res.valptr(bstart);
      double* bhres =  (double*) hres.valptr(bstart);
      double* bp    =  (double*) p.valptr(bstart);
      double* bhp   =  (double*) hp.valptr(bstart);

      int  sb = sbsizes[j];

      memcpy(kmat, kmat_store + j*ldsq, ldsq*sizeof(double)); 
      memcpy(mmat, mmat_store + j*ldsq, ldsq*sizeof(double)); 

      int dimp  = 2*sb;    // size of the current projected problem 
      dsygvd(&itype, &jobz, &uplo, &dimp, kmat, &ld, mmat, &ld, ev, work, &lwork, iwork, &liwork, &info);
      assert ( info == 0 );
  
      // update bwf and bp    
      size_t copy_size = sizeof(double)*mloc2*sb; 
      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bres, &mloc2, kmat + sb, &ld, &dzero, buffer,  &mloc2); 
      memcpy(bp, buffer, copy_size); 

      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bwf, &mloc2, kmat, &ld, &done, buffer,  &mloc2); 
      memcpy(bwf, buffer, copy_size); 
     
      // update bhwf and bhp     
      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bhres, &mloc2, kmat + sb, &ld, &dzero, buffer,  &mloc2); 
      memcpy(bhp, buffer, copy_size);  

      dgemm(&notrans, &notrans,  &ngwl2, &sb,  &sb, &done, bhwf, &mloc2, kmat, &ld, &done, buffer,  &mloc2); 
      memcpy(bhwf, buffer, copy_size); 

   } // end solve nsb projected problems and update subblocks

   delete [] kmat;
   delete [] kmat_store;
   delete [] mmat;   
   delete [] mmat_store;   
   delete [] buffer;
   delete [] ev;
   delete [] work;
   delete [] iwork;


}
////////////////////////////////////////////////////////////////////////////////






