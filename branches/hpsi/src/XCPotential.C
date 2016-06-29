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
// XCPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include "XCPotential.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "blas.h" // daxpy, dcopy
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
XCPotential::XCPotential(ChargeDensity& cd, const string functional_name):
    cd_(cd), ctxt_(cd.vcontext()), vft_(*cd_.vft()), vbasis_(*cd_.vbasis()),
    cd_ecalc_(cd)
{
   tddft_involved_ = false;
   initialize(functional_name);
}
////////////////////////////////////////////////////////////////////////////////
// separate constructor for TDDFT runs
XCPotential::XCPotential(ChargeDensity& cd, const string functional_name, ChargeDensity& cd_ecalc):
  cd_(cd), ctxt_(cd.vcontext()), vft_(*cd_.vft()), vbasis_(*cd_.vbasis()), cd_ecalc_(cd_ecalc)
{
   tddft_involved_ = true;
   initialize(functional_name);
}

////////////////////////////////////////////////////////////////////////////////
void XCPotential::initialize(string functional_name)
{
  if ( functional_name == "LDA" ) {
     if (cd_.nlcc())
        xcf_ = new LDAFunctional(cd_.xcrhor);
     else
        xcf_ = new LDAFunctional(cd_.rhor);
  }
  else if ( functional_name == "PBE" ) {
     if (cd_.nlcc())
        xcf_ = new PBEFunctional(cd_.xcrhor);
     else
        xcf_ = new PBEFunctional(cd_.rhor);
  }
  else if ( functional_name == "PBEsol" ) {
     if (cd_.nlcc())
        xcf_ = new PBESolFunctional(cd_.xcrhor);
     else
        xcf_ = new PBESolFunctional(cd_.rhor);
  }
  else if ( functional_name == "PBErev" ) {
     if (cd_.nlcc())
        xcf_ = new PBERevFunctional(cd_.xcrhor);
     else
        xcf_ = new PBERevFunctional(cd_.rhor);
  }
  else if ( functional_name == "BLYP" ) {
     if (cd_.nlcc())
        xcf_ = new BLYPFunctional(cd_.xcrhor);
     else
        xcf_ = new BLYPFunctional(cd_.rhor);
  }
  else {
    throw XCPotentialException("unknown functional name");
  }
  nspin_ = cd_.rhor.size();
  ngloc_ = vbasis_.localsize();
  np012loc_ = vft_.np012loc();
  
  if ( xcf_->isGGA() ) {
    tmp1.resize(ngloc_);
    if ( nspin_ > 1 )
      tmp2.resize(ngloc_);
    vxctmp.resize(nspin_);
    for ( int ispin = 0; ispin < nspin_; ispin++ )
      vxctmp[ispin].resize(np012loc_);
    tmpr.resize(np012loc_);
  }
}
////////////////////////////////////////////////////////////////////////////////
XCPotential::~XCPotential(void)
{
  delete xcf_;
}

////////////////////////////////////////////////////////////////////////////////
void XCPotential::update(vector<vector<double> >& vr)
{
  // compute exchange-correlation energy and add vxc potential to vr[ispin][ir]
  
  // Input: total electronic density in:
  //   vector<vector<double> >           cd_.rhor[ispin][ir] (real space)
  //   vector<vector<complex<double> > > cd_.rhog[ispin][ig] (Fourier coeffs)
  // The array cd_.rhog is only used if xcf->isGGA() == true
  // to compute the density gradients
  
  // Output: (through member function xcf())
  //
  // exc_, dxc, dxc0_, dxc1_, dxc2_
  //
  // LDA Functional:
  //   exc_, dxc
  //   spin unpolarized: xcf()->exc, xcf()->vxc1
  //   spin polarized:   xcf()->exc, xcf()->vxc1_up, xcf()->vxc1_dn
  //
  // GGA Functional: (through member function xcf())
  //   exc_, dxc, dxc0_, dxc1_, dxc2_
  //   spin unpolarized: xcf()->exc, xcf()->vxc1, xcf()->vxc2
  //   spin polarized:   xcf()->exc_up, xcf()->exc_dn, 
  //                     xcf()->vxc1_up, xcf()->vxc1_dn
  //                     xcf()->vxc2_upup, xcf()->vxc2_dndn, 
  //                     xcf()->vxc2_updn, xcf()->vxc2_dnup
  
   if ( !xcf_->isGGA() )
   {
    // LDA functional
 
    xcf_->setxc();

    exc_ = 0.0;
    const double *const e = xcf_->exc;
    const int size = xcf_->np();
 
    if ( nspin_ == 1 ) {
       // unpolarized
       double* rh;
       if (tddft_involved_)
          rh = (double*)&(cd_ecalc_.rhor[0][0]);
       else
          rh = (double*)xcf_->rho;
      const double *const v = xcf_->vxc1;
      for ( int i = 0; i < size; i++ ) {
        exc_ += rh[i] * e[i];
        vr[0][i] += v[i];
      }
    }
    else {
       // spin polarized
       //const double *const rh_up = xcf_->rho_up;
       //const double *const rh_dn = xcf_->rho_dn;
       double* rh_up;
       double* rh_dn;
       if (tddft_involved_)
       {
          rh_up = &(cd_ecalc_.rhor[0][0]);
          rh_dn = &(cd_ecalc_.rhor[1][0]);
       }
       else
       {
          rh_up = (double*)xcf_->rho_up;
          rh_dn = (double*)xcf_->rho_dn;
       }
       const double *const v_up = xcf_->vxc1_up;
       const double *const v_dn = xcf_->vxc1_dn;
       for ( int i = 0; i < size; i++ ) {                                                         
          exc_ += (rh_up[i] + rh_dn[i]) * e[i];
          vr[0][i] += v_up[i];
          vr[1][i] += v_dn[i];
       }                                                         
    }
    double tsum = exc_ * vbasis_.cell().volume() / vft_.np012();

    //ewd notes:  ctxt_ is a single-column context here
    ctxt_.dsum(1,1,&tsum,1);
    exc_ = tsum;

  }
  else {
    // GGA functional
    exc_ = 0.0;
    int size = xcf_->np();
    
    // compute grad_rho
    const double omega_inv = 1.0 / vbasis_.cell().volume();

    if ( nspin_ == 1 ) {
       const complex<double>* rhogptr = (cd_.nlcc() ? &cd_.xcrhog[0][0] : &cd_.rhog[0][0]);
       for ( int j = 0; j < 3; j++ ) {
        const double *const gxj = vbasis_.gx_ptr(j);
        for ( int ig = 0; ig < ngloc_; ig++ ) {
          /* i*G_j*c(G) */
           //tmp1[ig] = complex<double>(0.0,omega_inv*gxj[ig]) * cd_.rhog[0][ig];
          tmp1[ig] = complex<double>(0.0,omega_inv*gxj[ig]) * rhogptr[ig];
        }
        vft_.backward(&tmp1[0],&tmpr[0]);
        int inc2=2, inc1=1;
        double *grj = xcf_->grad_rho[j];
        dcopy(&np012loc_,(double*)&tmpr[0],&inc2,grj,&inc1);
       }
    }
    else {
      for ( int j = 0; j < 3; j++ ) {
        const double *const gxj = vbasis_.gx_ptr(j);
        //const complex<double>* rhg0 = &cd_.rhog[0][0];
        //const complex<double>* rhg1 = &cd_.rhog[1][0];
        const complex<double>* rhg0 = (cd_.nlcc() ? &cd_.xcrhog[0][0] : &cd_.rhog[0][0]);
        const complex<double>* rhg1 = (cd_.nlcc() ? &cd_.xcrhog[1][0] : &cd_.rhog[1][0]);
        
        for ( int ig = 0; ig < ngloc_; ig++ ) {
          /* i*G_j*c(G) */
          const complex<double> igxj(0.0,omega_inv*gxj[ig]);
          const complex<double> c0 = *rhg0++;
          const complex<double> c1 = *rhg1++;
          tmp1[ig] = igxj * c0; 
          tmp2[ig] = igxj * c1; 
        }
        if (vbasis_.real()) {
          vft_.backward(&tmp1[0],&tmp2[0],&tmpr[0]);
          double *grj_up = xcf_->grad_rho_up[j];
          double *grj_dn = xcf_->grad_rho_dn[j];
          int inc2=2, inc1=1;
          double* p = (double*) &tmpr[0];
          dcopy(&np012loc_,p,  &inc2,grj_up,&inc1);
          dcopy(&np012loc_,p+1,&inc2,grj_dn,&inc1);
        }
        else {
          vft_.backward(&tmp1[0],&tmpr[0]);
          double *grj_up = xcf_->grad_rho_up[j];
          int inc2=2, inc1=1;
          dcopy(&np012loc_,(double*)&tmpr[0],  &inc2,grj_up,&inc1);

          vft_.backward(&tmp2[0],&tmpr[0]);
          double *grj_dn = xcf_->grad_rho_dn[j];
          dcopy(&np012loc_,(double*)&tmpr[0],  &inc2,grj_dn,&inc1);
        }
      } // j
    }
    
    xcf_->setxc();
    
    // compute xc potential
    // take divergence of grad(rho)*vxc2

    // compute components of grad(rho) * vxc2
    if ( nspin_ == 1 )
    {
      for ( int j = 0; j < 3; j++ )
      {
        const double *const gxj = vbasis_.gx_ptr(j);
        const double *const grj = xcf_->grad_rho[j];
        const double *const v2 = xcf_->vxc2;
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          tmpr[ir] = grj[ir] * v2[ir];
        }
        // derivative
        vft_.forward(&tmpr[0],&tmp1[0]);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          // i*G_j*c(G)
          tmp1[ig] *= complex<double>(0.0,gxj[ig]);
        }
        // back to real space
        vft_.backward(&tmp1[0],&tmpr[0]);
        // accumulate div(vxc2*grad_rho) in vxctmp
        double one = 1.0;
        int inc1 = 1, inc2 = 2;
        if ( j == 0 )
        {
          dcopy(&np012loc_,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
        }
        else
        {
          daxpy(&np012loc_,&one,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
        }
      }
    }
    else
    {
      double *v2_upup = xcf_->vxc2_upup;
      double *v2_updn = xcf_->vxc2_updn;
      double *v2_dnup = xcf_->vxc2_dnup;
      double *v2_dndn = xcf_->vxc2_dndn;
      for ( int j = 0; j < 3; j++ )
      {
        const double *gxj = vbasis_.gx_ptr(j);
        const double *grj_up = xcf_->grad_rho_up[j];
        const double *grj_dn = xcf_->grad_rho_dn[j];
        if (vbasis_.real()) {
          for ( int ir = 0; ir < np012loc_; ir++ )
          {
            const double re = v2_upup[ir] * grj_up[ir] + v2_updn[ir] * grj_dn[ir];
            const double im = v2_dnup[ir] * grj_up[ir] + v2_dndn[ir] * grj_dn[ir];
            tmpr[ir] = complex<double>(re,im);
          }
          // derivative
          vft_.forward(&tmpr[0],&tmp1[0],&tmp2[0]);
        }
        else {
          for ( int ir = 0; ir < np012loc_; ir++ )
            tmpr[ir] = complex<double>(v2_upup[ir] * grj_up[ir] + v2_updn[ir] * grj_dn[ir],0.0);
          vft_.forward(&tmpr[0],&tmp1[0]);
          for ( int ir = 0; ir < np012loc_; ir++ )
            tmpr[ir] = complex<double>(v2_dnup[ir] * grj_up[ir] + v2_dndn[ir] * grj_dn[ir],0.0);
          vft_.forward(&tmpr[0],&tmp2[0]);
        }
        
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          // i*G_j*c(G)
          const complex<double> igxj(0.0,gxj[ig]);
          tmp1[ig] *= igxj;
          tmp2[ig] *= igxj;
        }

        if (vbasis_.real()) {
          vft_.backward(&tmp1[0],&tmp2[0],&tmpr[0]);
          // accumulate div(vxc2*grad_rho) in vxctmp
          double one = 1.0;
          int inc1 = 1, inc2 = 2;
          double* p = (double*) &tmpr[0];
          if ( j == 0 )
          {
            dcopy(&np012loc_,p  ,&inc2,&vxctmp[0][0],&inc1);
            dcopy(&np012loc_,p+1,&inc2,&vxctmp[1][0],&inc1);
          }
          else
          {
            daxpy(&np012loc_,&one,p  ,&inc2,&vxctmp[0][0],&inc1);
            daxpy(&np012loc_,&one,p+1,&inc2,&vxctmp[1][0],&inc1);
          }
        }
        else {
          vft_.backward(&tmp1[0],&tmpr[0]);
          // accumulate div(vxc2*grad_rho) in vxctmp
          double one = 1.0;
          int inc1 = 1, inc2 = 2;
          if ( j == 0 )
            dcopy(&np012loc_,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
          else
            daxpy(&np012loc_,&one,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
          vft_.backward(&tmp2[0],&tmpr[0]);
          if ( j == 0 )
            dcopy(&np012loc_,(double*)&tmpr[0],&inc2,&vxctmp[1][0],&inc1);
          else
            daxpy(&np012loc_,&one,(double*)&tmpr[0],&inc2,&vxctmp[1][0],&inc1);
        }
      } // j
    }
    
    // add xc potential to local potential in vr[i]
    // div(vxc2*grad_rho) is stored in vxctmp[ispin][ir]

    double esum=0.0;
    if ( nspin_ == 1 )
    {
      const double *const e = xcf_->exc;
      const double *const v1 = xcf_->vxc1;
      const double *const v2 = xcf_->vxc2;
      //const double *const rh = xcf_->rho;
      double* rh;
      if (tddft_involved_)
         rh = &(cd_ecalc_.rhor[0][0]);
      else
         rh = (double*)xcf_->rho;
      {
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          esum += rh[ir] * e[ir];
          vr[0][ir] += v1[ir] + vxctmp[0][ir];
        }
      }
    }
    else
    {
      const double *const v1_up = xcf_->vxc1_up;
      const double *const v1_dn = xcf_->vxc1_dn;
      const double *const v2_upup = xcf_->vxc2_upup;
      const double *const v2_updn = xcf_->vxc2_updn;
      const double *const v2_dnup = xcf_->vxc2_dnup;
      const double *const v2_dndn = xcf_->vxc2_dndn;
      const double *const eup = xcf_->exc_up;
      const double *const edn = xcf_->exc_dn;
      //const double *const rh_up = xcf_->rho_up;
      //const double *const rh_dn = xcf_->rho_dn;
      double* rh_up;
      double* rh_dn;
      if (tddft_involved_)
      {
         rh_up = &(cd_ecalc_.rhor[0][0]);
         rh_dn = &(cd_ecalc_.rhor[1][0]);
      }
      else
      {
         rh_up = (double*)xcf_->rho_up;
         rh_dn = (double*)xcf_->rho_dn;
      }
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        double r_up = rh_up[ir];
        double r_dn = rh_dn[ir];
        esum += r_up * eup[ir] + r_dn * edn[ir];
        vr[0][ir] += v1_up[ir] + vxctmp[0][ir];
        vr[1][ir] += v1_dn[ir] + vxctmp[1][ir];
      }
    }

    double tsum = esum * vbasis_.cell().volume() / vft_.np012();
    ctxt_.dsum(1,1,&tsum,1);
    exc_ = tsum;

  }
}
////////////////////////////////////////////////////////////////////////////////
// AS: modified version of XCPotential::update(vector<vector<double> >& vr) which
// AS: leaves the potential untouched and only recalculates the energy term
void XCPotential::update_exc(vector<vector<double> >& vr)
{
   assert(tddft_involved_);

   if ( !xcf_->isGGA() )
   {
      // LDA functional
      
      exc_ = 0.0;
      const double *const e = xcf_->exc;
      const int size = xcf_->np();

      if ( nspin_ == 1 )
      {
         // unpolarized
         // AS: was previously const double *const rh = xcf_->rho;
         // AS: however, the energy has to be determined with the new wave function
         const double *const rh = &(cd_ecalc_.rhor[0][0]);
         for ( int i = 0; i < size; i++ )
         {
            exc_ += rh[i] * e[i];
         }
      }
      else
      {
         // spin polarized
         const double *const rh_up = xcf_->rho_up;
         const double *const rh_dn = xcf_->rho_dn;
         for ( int i = 0; i < size; i++ )
         {
            exc_ += (rh_up[i] + rh_dn[i]) * e[i];
         }
      }
      double tsum = exc_ * vbasis_.cell().volume() / vft_.np012();
      ctxt_.dsum(1,1,&tsum,1);
      exc_ = tsum;
   }
   else
   {
      // GGA functional
      exc_ = 0.0;
      int size = xcf_->np();

      // div(vxc2*grad_rho) is stored in vxctmp[ispin][ir]
      
      double esum=0.0;
      if ( nspin_ == 1 )
      {
         const double *const e = xcf_->exc;
         // AS: was previously const double *const rh = xcf_->rho;
         // AS: however, the energy has to be determined with the new wave function
         // AS: untested so far!
         const double *const rh = &(cd_ecalc_.rhor[0][0]);
         {
            for ( int ir = 0; ir < np012loc_; ir++ )
            {
               esum += rh[ir] * e[ir];
            }
         }
      }
      else
      {
         const double *const eup = xcf_->exc_up;
         const double *const edn = xcf_->exc_dn;
         const double *const rh_up = xcf_->rho_up;
         const double *const rh_dn = xcf_->rho_dn;
         for ( int ir = 0; ir < np012loc_; ir++ )
         {
            const double r_up = rh_up[ir];
            const double r_dn = rh_dn[ir];
            esum += r_up * eup[ir] + r_dn * edn[ir];
         }
      }
      double tsum = esum * vbasis_.cell().volume() / vft_.np012();
      ctxt_.dsum(1,1,&tsum,1);
      exc_ = tsum;
   }
}

////////////////////////////////////////////////////////////////////////////////
void XCPotential::compute_stress(valarray<double>& sigma_exc)
{
  // compute exchange-correlation contributions to the stress tensor
  
  if ( !xcf_->isGGA() )
  {
    // LDA functional
 
    dxc_ = 0.0;
    const double *const e = xcf_->exc;
 
    if ( nspin_ == 1 )
    {
      // unpolarized
      const double *const rh = xcf_->rho;
      const double *const v = xcf_->vxc1;
      const int size = xcf_->np();
      for ( int i = 0; i < size; i++ )
      {
        dxc_ += rh[i] * (e[i] - v[i]);
      }
    }
    else
    {
      // spin polarized
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      const double *const v_up = xcf_->vxc1_up;
      const double *const v_dn = xcf_->vxc1_dn;
      const int size = xcf_->np();
      for ( int i = 0; i < size; i++ )                          
      {                                                         
        const double rh = rh_up[i] + rh_dn[i];                  
        dxc_ += rh * e[i] - rh_up[i] * v_up[i] - rh_dn[i] * v_dn[i];  
      }                                                         
    }
    const double fac = 1.0 / vft_.np012();
    double tsum;
    // Next line: factor omega in volume element cancels 1/omega in 
    // definition of sigma_exc
    tsum = - fac * dxc_;
    ctxt_.dsum(1,1,&tsum,1);
    
    // Note: contribution to sigma_exc is a multiple of the identity
    sigma_exc[0] = tsum;
    sigma_exc[1] = tsum;
    sigma_exc[2] = tsum;
    sigma_exc[3] = 0.0;
    sigma_exc[4] = 0.0;
    sigma_exc[5] = 0.0;
  }
  else
  {
    // GGA functional
    
    double dsum=0.0,sum0=0.0,sum1=0.0,sum2=0.0,
           sum3=0.0,sum4=0.0,sum5=0.0;
    if ( nspin_ == 1 )
    {
      const double *const e = xcf_->exc;
      const double *const v1 = xcf_->vxc1;
      const double *const v2 = xcf_->vxc2;
      const double *const rh = xcf_->rho;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        dsum += rh[ir] * ( e[ir] - v1[ir] );
        const double grx = xcf_->grad_rho[0][ir];
        const double gry = xcf_->grad_rho[1][ir];
        const double grz = xcf_->grad_rho[2][ir];
        const double grx2 = grx * grx;
        const double gry2 = gry * gry;
        const double grz2 = grz * grz;
        const double grad2 = grx2 + gry2 + grz2;
        const double v2t = v2[ir];
        sum0 += ( grad2 + grx2 ) * v2t;
        sum1 += ( grad2 + gry2 ) * v2t;
        sum2 += ( grad2 + grz2 ) * v2t;
        sum3 += grx * gry * v2t;
        sum4 += gry * grz * v2t;
        sum5 += grx * grz * v2t;
      }
    }
    else
    {
      const double *const v1_up = xcf_->vxc1_up;
      const double *const v1_dn = xcf_->vxc1_dn;
      const double *const v2_upup = xcf_->vxc2_upup;
      const double *const v2_updn = xcf_->vxc2_updn;
      const double *const v2_dnup = xcf_->vxc2_dnup;
      const double *const v2_dndn = xcf_->vxc2_dndn;
      const double *const eup = xcf_->exc_up;
      const double *const edn = xcf_->exc_dn;
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        const double r_up = rh_up[ir];
        const double r_dn = rh_dn[ir];
        dsum += r_up * ( eup[ir] - v1_up[ir] ) +
                r_dn * ( edn[ir] - v1_dn[ir] );
 
        const double grx_up = xcf_->grad_rho_up[0][ir];
        const double gry_up = xcf_->grad_rho_up[1][ir];
        const double grz_up = xcf_->grad_rho_up[2][ir];
        const double grx2_up = grx_up * grx_up;
        const double gry2_up = gry_up * gry_up;
        const double grz2_up = grz_up * grz_up;
        const double grad2_up = grx2_up + gry2_up + grz2_up;
 
        const double grx_dn = xcf_->grad_rho_dn[0][ir];
        const double gry_dn = xcf_->grad_rho_dn[1][ir];
        const double grz_dn = xcf_->grad_rho_dn[2][ir];
        const double grx2_dn = grx_dn * grx_dn;
        const double gry2_dn = gry_dn * gry_dn;
        const double grz2_dn = grz_dn * grz_dn;
        const double grad2_dn = grx2_dn + gry2_dn + grz2_dn;
 
        const double grad_up_grad_dn = grx_up * grx_dn +
                                       gry_up * gry_dn +
                                       grz_up * grz_dn;

        const double v2_upup_ir = v2_upup[ir];
        const double v2_updn_ir = v2_updn[ir];
        const double v2_dnup_ir = v2_dnup[ir];
        const double v2_dndn_ir = v2_dndn[ir];

        sum0 += v2_upup_ir * ( grad2_up + grx2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + grx_up * grx_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + grx_dn * grx_up ) +
                v2_dndn_ir * ( grad2_dn + grx2_dn );
 
        sum1 += v2_upup_ir * ( grad2_up + gry2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + gry_up * gry_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + gry_dn * gry_up ) +
                v2_dndn_ir * ( grad2_dn + gry2_dn );
 
        sum2 += v2_upup_ir * ( grad2_up + grz2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + grz_up * grz_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + grz_dn * grz_up ) +
                v2_dndn_ir * ( grad2_dn + grz2_dn );
 
        sum3 += v2_upup_ir * grx_up * gry_up +
                v2_updn_ir * grx_up * gry_dn +
                v2_dnup_ir * grx_dn * gry_up +
                v2_dndn_ir * grx_dn * gry_dn;
 
        sum4 += v2_upup_ir * gry_up * grz_up +
                v2_updn_ir * gry_up * grz_dn +
                v2_dnup_ir * gry_dn * grz_up +
                v2_dndn_ir * gry_dn * grz_dn;
 
        sum5 += v2_upup_ir * grx_up * grz_up +
                v2_updn_ir * grx_up * grz_dn +
                v2_dnup_ir * grx_dn * grz_up +
                v2_dndn_ir * grx_dn * grz_dn;
      }
    }
    double fac = 1.0 / vft_.np012();
    double tsum[6];
    // Next line: factor omega in volume element cancels 1/omega in 
    // definition of sigma_exc
    tsum[0] = - fac * ( dsum + sum0 );
    tsum[1] = - fac * ( dsum + sum1 );
    tsum[2] = - fac * ( dsum + sum2 );
    tsum[3] = - fac * sum3;
    tsum[4] = - fac * sum4;
    tsum[5] = - fac * sum5;
    ctxt_.dsum(6,1,&tsum[0],6);
    
    sigma_exc[0] = tsum[0];
    sigma_exc[1] = tsum[1];
    sigma_exc[2] = tsum[2];
    sigma_exc[3] = tsum[3];
    sigma_exc[4] = tsum[4];
    sigma_exc[5] = tsum[5];
  }
}
