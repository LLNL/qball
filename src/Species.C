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
// Species.C:
//
////////////////////////////////////////////////////////////////////////////////

#include "Species.h"
#include "spline.h"
#include "sinft.h"
#include "SphericalIntegration.h"
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
static double simpsn ( int n, double *t )
{
  // extended Simpson's rule, applied to overlapping segments
  // (from Numerical Recipes, first edition, Eq. 4.1.14, p. 122)

  const double c0 =  17.0/48.0, c1 = 59.0/48.0,
               c2 =  43.0/48.0, c3 = 49.0/48.0;
  double sum = c0 * ( t[0] + t[n-1] ) +
               c1 * ( t[1] + t[n-2] ) +
               c2 * ( t[2] + t[n-3] ) +
               c3 * ( t[3] + t[n-4] );
       
  for ( int i = 4; i < n-4; i++ )
  {
    sum += t[i];
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
Species::Species(const Context& ctxt, string name) : ctxt_(ctxt), name_(name),
zval_(-1), mass_(0.0), lmax_(-1), deltar_(0.0), atomic_number_(0), 
llocal_(-1), nquad_(-1), rquad_(0.0),
rcps_(0.0), uri_(""), description_("undefined"), symbol_("Undef"),
initsize_(-1), fix_rcps(false), hubbard_u_(0.0), hubbard_l_(-1),hubbard_alpha_(0.0),
nbeta_(-1),usoft_(false),nlcc_(false)
{}
  
////////////////////////////////////////////////////////////////////////////////
bool Species::initialize(double rcpsval)
{
  // initialize the Species
  rcps_ = rcpsval;

  assert(description_ != "undefined");
  
  const double fpi = 4.0 * M_PI;
  vector<vector<vector<double> > > vnlr;
  
  //const int np = vps_[0].size();
  //ewd store initsize to keep ndft_ from growing too large on subsequent calls
  if (initsize_ == -1) {
    if(!oncv_) {
      initsize_ = vps_[0].size();
    } else {
      initsize_ = projectors_[0][0].size();
    }
  }
  
  const int np = initsize_;

  if (zval_ < 0) throw SpeciesInitException("zval_ < 0");
  if (rcps_ < 0.0) throw SpeciesInitException("rcps_ < 0");
  if (mass_ < 0.0) throw SpeciesInitException("mass_ < 0");
  if (lmax_ < 0 || lmax_ > 3) throw SpeciesInitException("lmax_ <0 or lmax_ >3");

  if (!usoft_ && !oncv_){
    if (vps_.size() < lmax_+1) throw SpeciesInitException("vps_.size < lmax_+1");
  }

  if(!oncv_){
    if (llocal_ < 0 || llocal_ > lmax_) throw SpeciesInitException("llocal_ < 0 || llocal_ > lmax_");
  }
  
  if (!usoft_ && !oncv_) {
    if ( nquad_ == 0 ) // KB potential
    {
      for ( int l = 0; l <= lmax_; l++ )
        if ( l != llocal_ && phi_[l].size() == 0 )
          throw SpeciesInitException("phi_[l] undefined for non-local projector");
    }
  
    if ( nquad_ < 0 ) throw SpeciesInitException("nquad < 0");
    if ( nquad_ > 0 && rquad_ <= 0.0 ) 
      throw SpeciesInitException("semilocal with rquad_ <= 0");
  }
  
  // compute number of non-local projectors nlm_
  nlm_ = 0;
  for ( int l = 0; l <= lmax_; l++ )
  {
    if ( l != llocal_ )
    {
      nlm_ += nchannels_*(2*l + 1);
    }
  }
  
  // compute ndft_: size of radial FFT array
  // ndft_ is the second power of 2 larger than np
  ndft_ = 1;
  while ( ndft_ < np )
  {
    ndft_ *= 2;
    if ( ndft_ > 8192 )
    {
      if (ctxt_.oncoutpe()) 
        cout << "<WARNING> Species::initialize: ndft_ > 8192, np = "
           << np << "</WARNING>" << endl;
      //return false;
    }
  }
  ndft_ *= 2;

  //ewd DEBUG:  try increasing accuracy of vnlg spline by extending potential out to large r
  const double rmax = 40.0;
  while (deltar_*ndft_ < rmax) 
    ndft_ *= 2;
  if (ctxt_.oncoutpe()) 
     cout << "<!-- Species " << name_ << ":  extending grid to rmax = " << rmax 
         << " to increase vnlg resolution (" << ndft_ << " pts) -->" << endl;

  //ewd DEBUG
  //ndft_ *= 2;
  

  rps_.resize(ndft_);
  for ( int i = 0; i < ndft_; i++ )
    rps_[i] = i * deltar_;
    
  if (usoft_) {
    vps_spl_.resize(1);
    vps_[0].resize(ndft_);
    vps_spl_[0].resize(ndft_);
  }
  else if(!oncv_){
    vps_spl_.resize(lmax_+1);
    phi_spl_.resize(lmax_+1);
  
    for ( int l = 0; l <= lmax_; l++ )
    {
      vps_[l].resize(ndft_);
      phi_[l].resize(ndft_);
      vps_spl_[l].resize(ndft_);
      phi_spl_[l].resize(ndft_);
    }
  }
  
  // extend rps and vps_ to full mesh (up to i==ndft_-1)
  
  vector<double> fint(ndft_);
  
  wsg_.resize(lmax_+1);
  for(int ll = 0; ll < lmax_ + 1; ll++) wsg_[ll].resize(nchannels_);
  gspl_.resize(ndft_);
  vlocg_.resize(ndft_);
  vlocg_spl.resize(ndft_);
  
  if (!usoft_) {
    vnlg_.resize(lmax_+1);
    vnlg_spl.resize(lmax_+1);
    for(int ll = 0; ll < lmax_ + 1; ll++){
      vnlg_[ll].resize(nchannels_);
      vnlg_spl[ll].resize(nchannels_);
    }
  }
  
  vector<double> vlocr(ndft_);

  if (usoft_) {
    assert(llocal_ == 0); // SpeciesReader currently assumes this for ultrasoft

    // extend vlocal to full grid
    for ( int i = np; i < ndft_; i++ )
      vps_[0][i] = - zval_ / rps_[i];

    // spline vlocal
    spline(&rps_[0],&vps_[0][0],ndft_,
           SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&vps_spl_[0][0]);

  } else if(oncv_){
    assert(llocal_ == -1);
    
    vloc_.resize(ndft_);
    for (int ip = np; ip < ndft_; ip++ ) vloc_[ip] = - zval_/rps_[ip];

    for(int ll = 0; ll <= lmax_; ll++ ){
      for(int ii = 0; ii < nchannels_; ii++){
	projectors_[ll][ii].resize(ndft_);
	for (int ip = np; ip < ndft_; ip++ ) projectors_[ll][ii][ip] = 0.0;
      }
    }
    
  } else {
    // Extend vps_[l][i] up to ndft_ using -zv/r
    for ( int l = 0; l <= lmax_; l++ )
    {
      for ( int i = np; i < ndft_; i++ )
      {
        vps_[l][i] = - zval_ / rps_[i];
        phi_[l][i] = 0.0;
      }
    }

    // compute spline coefficients of vps_ and phi_
    for ( int l = 0; l <= lmax_; l++ )
    {
      spline(&rps_[0],&vps_[l][0],ndft_,
             SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&vps_spl_[l][0]);
    }
    for ( int l = 0; l <= lmax_; l++ )
    {
      if ( l != llocal_ )
        spline(&rps_[0],&phi_[l][0],ndft_,
               SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&phi_spl_[l][0]);
    }
  }

  if(!usoft_){
    vnlr.resize(lmax_ + 1);
    //  vnlg_ is dimensioned ndft_+1 since it is passed to cosft1
    //  See Numerical Recipes 2nd edition for an explanation.
    for ( int l = 0; l <= lmax_; l++ ){
      vnlr[l].resize(nchannels_);
      for(int ic = 0; ic < nchannels_; ic++) vnlr[l][ic].resize(ndft_);
    }
  }
  
  // local potential: subtract the long range part due to the smeared charge
  // Next line: constant is 2/sqrt(pi)
  // math.h: # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
  if(!oncv_) {
    substract_long_range_part(vps_[llocal_], vlocr);
  } else {
    substract_long_range_part(vloc_, vlocr);
  }
  
  //  Prepare the function vlocr to be used later in the Bessel transforms:
  //
  //  local potential: v(G) = 4 pi \int r^2 vloc(r) sin(Gr)/(Gr) dr
  //  -> store 4 pi r dr vps_(lmax_) in vlocr(i,is)
  //
  //  the Bessel transform is then:
  //  v(G) = 1/G \sum_r sin(Gr) vlocr

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocr[i] *= fpi * rps_[i] * deltar_;
  }
  //  Local potential
  //  Compute Fourier coefficients of the local potential
  //  vlocr[i] contains 4 pi r dr vpsr(lmax_)
  //  v(G) = 4 pi \int r^2 vpsr(r) sin(Gr)/(Gr) dr
  //       = 1/G \sum_r sin(Gr) vlocr
  //
  //  v(G=0) by simpson integration
  //  v(G) = 4 pi \int r^2 vpsr(r) dr
  //       = \sum_r r vlocr
  //
  //  N.B. vlocr[i] contains 4 pi r dr (vps_(lmax_)-v_pseudocharge(r))
  //  Simpson integration up to ndft_ (V is zero beyond that point)
  //  Use fint as temporary array for integration

  for ( int i = 0; i < ndft_; i++ )
  {
    fint[i] = vlocr[i] * rps_[i];
  }
  double v0 = simpsn(ndft_,&fint[0]);

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocg_[i] = vlocr[i];
  }

  sinft(&vlocg_[0],ndft_);

  //  Divide by g's
  gspl_[0] = 0.0;
  vlocg_[0] = v0;
  double fac = M_PI/(ndft_*deltar_);
  for ( int i = 1; i < ndft_; i++ )
  {
    gspl_[i] = i * fac;
    vlocg_[i] /= gspl_[i];
  }

  //ewd DEBUG  
  if (ctxt_.oncoutpe()) 
     cout << "SPECIES.ndft = " << ndft_ << ", np = " << np << ", rmax = " << ndft_*deltar_ << ", gmax = " << gspl_[ndft_-1] << ", hubbard_l = " << hubbard_l_ << endl;

  //  Initialize cubic spline interpolation for local potential Vloc(G)
  //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

  spline(&gspl_[0],&vlocg_[0],ndft_,
         SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&vlocg_spl[0]);

  // Non-local KB projectors
  if ( !usoft_ && nquad_ == 0 && non_local()){

    if(!oncv_){
    
      for ( int l = 0; l <= lmax_; l++ ){
	wsg_[l][0] = 0.0;

	if ( l != llocal_ ){
	  // for KB potentials, compute weights wsg[l]
	  //  Compute weight wsg_[l] by integration on the linear mesh
	  for ( int i = 0; i < ndft_; i++ )
	    {
	      double tmp = phi_[l][i] * rps_[i];
	      fint[i] = ( vps_[l][i] - vps_[llocal_][i] ) * tmp * tmp;
	    }
	  double tmp = simpsn(ndft_,&fint[0]);
	  assert(tmp != 0.0);
	  // Next lines: store 1/<phi|delta_v| phi > in wsg[is][l]
	  wsg_[l][0] = 1.0 / ( deltar_ * tmp );

	  if (ctxt_.oncoutpe()) 
	    cout << "<!-- Kleinman-Bylander normalization term wsg[" << l << "] = " << wsg_[l][0] << " -->" << endl;

	  //   compute non-local projectors:
	  //   w(G) = Ylm(G) i^l 4 pi \int r^2 phi_l(r) j_l(Gr) v_l(r) dr
	  //   -> store 4 pi v_l(r) phi_l(r) dr in vnlr[l][i]
	  //   the Bessel transform is then done by
	  //   l=0: j_0(Gr) = sin(Gr)/(Gr)
	  //        w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
	  //   l=1: j_1(Gr) = sin(Gr)/(Gr)^2 - cos(Gr)/Gr
	  //        w_1(G) = Ylm(G) i^-1 1/G^2 \sum_r sin(Gr) vnlr -
	  //                 Ylm(G) i^-1 1/G   \sum_r cos(Gr) r vnlr
	  //   l=2: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
	  //        w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
	  //                                1/G   \sum_r sin(Gr)*r  vnlr -
	  //                                3/G^2 \sum_r cos(Gr)    vnlr )
	  //   l=3: j_3(Gr) = sin(Gr)*(15/(Gr)^3-6/Gr) - cos(Gr)*(15/(Gr)^2 - 1)
	  //        w_3(G) = Ylm(G) i^-3 (  15/G^3 \sum_r sin(Gr)/r vnlr -
	  //                                 6/G   \sum_r sin(Gr)*r vnlr -
	  //                                15/G^2 \sum_r cos(Gr)   vnlr + 
	  //                                       \sum_r cos(Gr)*r^2 vnlr )
 
	  for ( int i = 0; i < ndft_; i++ )
	    {
	      vnlr[l][0][i] = fpi*deltar_*( vps_[l][i] - vps_[llocal_][i] )*phi_[l][i];
	    }
	}
      }

    } else {

      for ( int l = 0; l < lmax_ + 1; l++ ){
	for (int ic = 0; ic < nchannels_; ic++){
	  wsg_[l][ic] = dij_[l][ic][ic];
	  for (int ip = 0; ip < ndft_; ip++) vnlr[l][ic][ip] = fpi*deltar_*projectors_[l][ic][ip];
	}
      }
      
    }

    //  compute radial Fourier transforms of vnlr
    for ( int l = 0; l <= lmax_; l++ )
    {
      for(int ic = 0; ic < nchannels_; ic++){

	if ( l != llocal_ ){
	  vnlg_[l][ic].resize(ndft_ + 1);
	  vnlg_spl[l][ic].resize(ndft_ + 1);

	  bessel_trans(l, vnlr[l][ic], vnlg_[l][ic]);
	  if ( l == 0 || l == 2 || l == 3){
	    
	    //  Initialize cubic spline interpolation
	    //  Use zero first derivative at G=0 and natural (y"=0) at Gmax
	    spline(&gspl_[0], &vnlg_[l][ic][0], ndft_, SPLINE_FLAT_BC, SPLINE_NATURAL_BC, &vnlg_spl[l][ic][0]);
	    
	  } else if ( l == 1 ){
	    // Initialize spline interpolation
	    // Use natural first derivative at G=0 and natural (y"=0) at Gmax
	    
	    spline(&gspl_[0], &vnlg_[l][ic][0], ndft_, SPLINE_NATURAL_BC, SPLINE_NATURAL_BC, &vnlg_spl[l][ic][0]);
	  }
	} // l != llocal_
	
      }
    } // l
  } // nquad_ == 0 && non_local() (KB projectors)

  // calculate radial Fourier transforms of atomic orbitals for DFT+U
  if (hubbard_l_ > -1)
  {
    assert(!usoft_); //ewd:  DFT+U not implemented for ultrasoft
    assert(!oncv_);  //We don't have orbitals for ONCV
    
    phir_.resize(ndft_+1);
    phig_.resize(ndft_+1);
    phig_spl_.resize(ndft_+1);

    for ( int i = 0; i < ndft_; i++ )
      phir_[i] = fpi * deltar_ * phi_[hubbard_l_][i];

    //ewd:  normalize phir_
    double phinorm = 0.0;
    for ( int i = 0; i < ndft_; i++ )
      phinorm += fpi*deltar_*rps_[i]*rps_[i]*phi_[hubbard_l_][i]*phi_[hubbard_l_][i];
    for ( int i = 0; i < ndft_; i++ )
      phir_[i] /= phinorm;

    if (ctxt_.mype()==0)
      cout << "Species.C:  normalizing orbitals for Hubbard projection by phirnorm = " << phinorm << endl;

    bessel_trans(hubbard_l_,phir_,phig_);

    if ( hubbard_l_ == 0 || hubbard_l_ == 2 || hubbard_l_ == 3 )
    {
      spline(&gspl_[0],&phig_[0],ndft_,SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&phig_spl_[0]);
    }
    else if ( hubbard_l_ == 1 )
    {
      spline(&gspl_[0],&phig_[0],ndft_,
                 SPLINE_NATURAL_BC,SPLINE_NATURAL_BC,&phig_spl_[0]);
    }
  }

  // non-linear core correction
  if (nlcc_)
  {
     // extend to ndft_, pad with zeroes
     rhor_nlcc_.resize(ndft_);
     rhog_nlcc_.resize(ndft_);
     rhog_nlcc_spl_.resize(ndft_);
     for ( int i = np; i < ndft_; i++ )
        rhor_nlcc_[i] = 0.0;
     
     for ( int i = 0; i < ndft_; i++ )
     {
        fint[i] = (rhor_nlcc_[i] * fpi * rps_[i] * deltar_) * rps_[i];
     }
     double rho0 = simpsn(ndft_,&fint[0]);

     for ( int i = 0; i < ndft_; i++ )
     {
        rhog_nlcc_[i] = rhor_nlcc_[i] * fpi * rps_[i] * deltar_;
     }

     sinft(&rhog_nlcc_[0],ndft_);

     //  Divide by g's
     rhog_nlcc_[0] = rho0;
     double fac = M_PI/(ndft_*deltar_);
     for ( int i = 1; i < ndft_; i++ )
     {
        gspl_[i] = i * fac;
        rhog_nlcc_[i] /= gspl_[i];
     }
     spline(&gspl_[0],&rhog_nlcc_[0],ndft_,
            SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&rhog_nlcc_spl_[0]);
  }
  
  // ultrasoft pseudopotentials
  if (usoft_) {
    betalmax_ = 0;
    nbetalm_ = 0;
    for (int b=0; b<nbeta_; b++) {
      nbetalm_ += 2*betal_[b]+1;
      if (betal_[b] > betalmax_)
        betalmax_ = betal_[b];
    }

    betalm_l_.resize(nbetalm_);
    betalm_m_.resize(nbetalm_);
    betaind_.resize(nbetalm_);
    int lmcnt = 0;
    for (int b=0; b<nbeta_; b++) {
      for (int m=-betal_[b]; m<=betal_[b]; m++) {
        betalm_l_[lmcnt] = betal_[b];
        betalm_m_[lmcnt] = m;
        betaind_[lmcnt] = b;
        lmcnt++;
      }
    }
    
    betag_.resize(nbeta_);
    betag_spl_.resize(nbeta_);
    for (int b=0; b<nbeta_; b++) {
      betar_[b].resize(ndft_);
      betag_[b].resize(ndft_);
      betag_spl_[b].resize(ndft_);
    }
    qfung_.resize(nqfun_);
    qfung_spl_.resize(nqfun_);
    for (int q=0; q<nqfun_; q++) {
      qfunr_[q].resize(ndft_);
      int ltotmin = abs(qfunl1_[q]-qfunl2_[q]);
      int ltotmax = qfunl1_[q]+qfunl2_[q];
      int nltot = ltotmax-ltotmin+1;
      qfung_[q].resize(nltot);
      qfung_spl_[q].resize(nltot);
      for (int l=0; l<nltot; l++) {
        qfung_[q][l].resize(ndft_);
        qfung_spl_[q][l].resize(ndft_);
      }
    }

    // extend betar_, qfunr_ to full grid by padding w. zeros
    for (int b=0; b<nbeta_; b++)
      for ( int i = np; i < ndft_; i++ )
        betar_[b][i] = 0.0;
    for (int q=0; q<nqfun_; q++)
      for ( int i = np; i < ndft_; i++ )
        qfunr_[q][i] = 0.0;
    
    vector<double> betaint(ndft_);
    for (int b=0; b<nbeta_; b++) {
      for (int i=0; i<ndft_; i++)
        betaint[i] = betar_[b][i]*fpi*deltar_;
      int betal = betal_[b];
      bessel_trans(betal,betaint,betag_[b]);
      if ( betal == 0 || betal == 2 || betal == 3 )
        spline(&gspl_[0],&betag_[b][0],ndft_,
               SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&betag_spl_[b][0]);
      else if ( betal == 1 )
        spline(&gspl_[0],&betag_[b][0],ndft_,
               SPLINE_NATURAL_BC,SPLINE_NATURAL_BC,&betag_spl_[b][0]);
    }

    // store index of qfun_nm, where n = 0...nbeta-1, m=n..nbeta-1
    qfunind_.resize(nbeta_);
    for (int b=0; b<nbeta_; b++)
      qfunind_[b].resize(nbeta_);
    for (int q=0; q<nqfun_; q++) {
      int b1 = qfunb1_[q];
      int b2 = qfunb2_[q];
      qfunind_[b1][b2] = q;
      qfunind_[b2][b1] = q;
    }
    
    vector<double> qfunint(ndft_);
    for (int q=0; q<nqfun_; q++) {
      qfunint[0] = 0.0;
      for (int i=1; i<ndft_; i++)
        qfunint[i] = qfunr_[q][i]*fpi*deltar_;

      int ltotmin = abs(qfunl1_[q]-qfunl2_[q]);
      int ltotmax = qfunl1_[q]+qfunl2_[q];
      int nltot = ltotmax-ltotmin+1;
      for (int lind=0; lind<nltot; lind++) {
        int qfunl = ltotmin + lind;
        //if rinner defined for this ltot, replace values at r<rinner w. qfcoeff series:
        //   Qnm(r<rinner) = r^ltot SUM qfcoeff[j]*r^2j
        if (rinner_.size() > qfunl) { 
          if (rinner_[qfunl] > 0.0) {
            int rinind = int(rinner_[qfunl]/deltar_) + 1;
            for (int ir=0; ir<rinind; ir++) {
              if (rps_[ir] < rinner_[qfunl]) {
                double rsq = rps_[ir]*rps_[ir];
                double val = qfcoeff_[q][qfunl][0];
                double rpow = rsq;
                int nqf = qfcoeff_[q][qfunl].size();
                for (int j=1; j<nqf; j++) {
                  val += qfcoeff_[q][qfunl][j]*rpow;
                  rpow *= rsq;
                }
                rpow = 1.0;
                for (int l=0; l<qfunl; l++)
                  rpow *= rps_[ir];
                val *= rpow;
                qfunint[ir] = val*fpi*deltar_;
              }
            }
          }
        }

        bessel_trans(qfunl,qfunint,qfung_[q][lind]);
        spline(&gspl_[0],&qfung_[q][lind][0],ndft_,
               SPLINE_NATURAL_BC,SPLINE_FLAT_BC,&qfung_spl_[q][lind][0]);

        /*
        if ( qfunl == 1 )
          spline(&gspl_[0],&qfung_[q][lind][0],ndft_,
               SPLINE_NATURAL_BC,SPLINE_NATURAL_BC,&qfung_spl_[q][lind][0]);
        else 
          spline(&gspl_[0],&qfung_[q][lind][0],ndft_,
               SPLINE_FLAT_BC,SPLINE_NATURAL_BC,&qfung_spl_[q][lind][0]);
        */

#ifdef PRINT_US_SPLINES        
        //ewd DEBUG
        const double gmax = 20.0;
        ostringstream oss1, oss2;
        oss1.width(3);  oss1.fill('0');  oss1 << q;
        oss2.width(1);  oss2.fill('0');  oss2 << qfunl;
        string qfoutfile = "qfung." + symbol_ + "." + oss1.str() + "." + oss2.str() + ".dat";
        ofstream os;
        os.open(qfoutfile.c_str(),ofstream::out);
        for ( int i = 1; i < ndft_-1; i++ )
        {
           double v;
           double g = 0.5*(gspl_[i+1]+gspl_[i]);
           if (g <= gmax) {
              splint(&gspl_[0],&qfung_[q][lind][0],&qfung_spl_[q][lind][0],ndft_,g,&v);
              os << gspl_[i] << "   " << qfung_[q][lind][i] << "   " << v << "    " << g << endl;
           }        
        }        
        os.close();
#endif
      }
    }
    
    // loop through all combinations of ang. momentum and store info 
    // on non-zero combos to avoid oversizing qnmg array
    SphericalIntegration si;

    qnm_lm1_.resize(nbetalm_*nbetalm_);
    qnm_lm2_.resize(nbetalm_*nbetalm_);
    ncgcoeff_.resize(nbetalm_*nbetalm_);
    cgcoeff_.resize(nbetalm_*nbetalm_);
    cgltot_.resize(nbetalm_*nbetalm_);
    cgmtot_.resize(nbetalm_*nbetalm_);

    nqtot_ = 0;
    const double cgeps = 1.E-9;  // ignore terms smaller than this
    for (int lm1 = 0; lm1<nbetalm_; lm1++) {
      for (int lm2 = lm1; lm2<nbetalm_; lm2++) {
        int l1 = betalm_l_[lm1];
        int m1 = betalm_m_[lm1];
        int l2 = betalm_l_[lm2];
        int m2 = betalm_m_[lm2];
        int lmin = abs(l1-l2);
        int lmax = l1+l2;
        int cgcnt = 0;
        for (int ltot=lmin; ltot<=lmax; ltot++) {
          for (int mtot=-ltot; mtot<=ltot; mtot++) {
            double coeff = si.clebsch_gordan_real(l1,m1,l2,m2,ltot,mtot);
            if (abs(coeff) > cgeps) {
              qnm_lm1_[nqtot_] = lm1;
              qnm_lm2_[nqtot_] = lm2;
              cgcoeff_[nqtot_].push_back(coeff);
              cgltot_[nqtot_].push_back(ltot);
              cgmtot_[nqtot_].push_back(mtot);
              cgcnt++;
            }
          }
        }
        if (cgcnt > 0) 
          ncgcoeff_[nqtot_++] = cgcnt;
      }
    }
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////////
void Species::calc_qnmg(Basis* basis_, vector<complex<double> > &qtotg, vector<double> &qaug)
{
  // Q^L_ij(G) was calculated and splined in initialize, so after charge
  // density/wf basis is set we call this function to compute
  //   Q_lilj(G) = SUM_LM c_LM^lilj Y_LM(G) Q^L_ij(G)
  //
  //   where:
  //     li,lj are ang. momentum indices of i,j beta fns (each li corresponds to diff. l,m values)
  //     L is combined total angular momentum of li,lj allowed by sum rule
  //     c_LM^lilj are real Clebsch-Gordan coefficients
  //     Y_LM(G) are real spherical harmonics
  //
  // qaug_ij = INT dr Q_ij(r) (coefficients to overlap matrix) 
  //         = Q_ij(G=0) 
  
  const int ngwl = basis_->localsize();
  const double *kpg   = basis_->kpg_ptr();
  const double *kpg_x = basis_->kpgx_ptr(0);
  const double *kpg_y = basis_->kpgx_ptr(1);
  const double *kpg_z = basis_->kpgx_ptr(2);
  const double omega = basis_->cell().volume();
  const double fpi = 4.0 * M_PI;
  assert(omega != 0.0);
  const double omega_inv = 1./omega;
  assert(nqtot_ > 0);
  qtotg.resize(nqtot_*ngwl);
  for (int i=0; i<nqtot_*ngwl; i++)
    qtotg[i] = complex<double>(0.0,0.0);
  qaug.resize(nbetalm_*nbetalm_);
  for (int i=0; i<nbetalm_*nbetalm_; i++)
    qaug[i] = 0.0;
  
  // calculate qtotg = SUM(ltot,mtot) (-i)^ltot cgcoeff*ylmg(ltot,mtot)*qfung(ltot)
  // mult by structure factor in separate function, so we can update it as atoms move
  for (int qind=0; qind<nqtot_; qind++) {
    int lm1 = qnm_lm1_[qind];
    int lm2 = qnm_lm2_[qind];
    int b1 = betaind_[lm1];
    int b2 = betaind_[lm2];
    int qfuni = qfunind_[b1][b2];
    int ncg = ncgcoeff_[qind];
    for (int icg = 0; icg < ncg; icg++) {
      int ltot = cgltot_[qind][icg];
      int mtot = cgmtot_[qind][icg];
      complex<double> ltotim = complex<double>(1.0,0.0);
      for (int l=0; l<ltot; l++)
        ltotim *= complex<double>(0.0,-1.0);

      double coeff = cgcoeff_[qind][icg];
      if (coeff != 0.0) {
        for (int ig=0; ig<ngwl; ig++) {
          double gnorm = kpg[ig];
          double qfg = qfung(qfuni,ltot,gnorm);
          double ylmg = ylm(ltot,mtot,kpg_x[ig],kpg_y[ig],kpg_z[ig]);
          qtotg[qind*ngwl+ig] += (coeff*ylmg*qfg)*ltotim;
        }
        double ylmzero = ylm(ltot,mtot,0.0,0.0,0.0);
        double qfzero = qfung(qfuni,ltot,0.0);
        qaug[qind] += coeff*ylmzero*qfzero*real(ltotim);
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void Species::qind_to_beta(int qind, int &b1, int &b2) {
  assert(qind<nqtot_);
  assert(qnm_lm1_.size() >= nqtot_);
  assert(qnm_lm2_.size() >= nqtot_);
  int lm1 = qnm_lm1_[qind];
  int lm2 = qnm_lm2_[qind];
  b1 = betaind_[lm1];
  b2 = betaind_[lm2];  
}
////////////////////////////////////////////////////////////////////////////////
void Species::qind_to_betalm(int qind, int &lm1, int &lm2) {
  assert(qind<nqtot_);
  assert(qnm_lm1_.size() >= nqtot_);
  assert(qnm_lm2_.size() >= nqtot_);
  lm1 = qnm_lm1_[qind];
  lm2 = qnm_lm2_[qind];
}
////////////////////////////////////////////////////////////////////////////////
double Species::dzero(int qind) {
  int n,m;
  qind_to_beta(qind,n,m);

  double dz = 0.0;
  int lm1,lm2;
  qind_to_betalm(qind,lm1,lm2);
  int l1 = betalm_l_[lm1];
  int m1 = betalm_m_[lm1];
  int l2 = betalm_l_[lm2];
  int m2 = betalm_m_[lm2];
  if (m1 == m2 && l1 == l2)
    dz = dzero_[n][m];

  return dz;
}
////////////////////////////////////////////////////////////////////////////////
void Species::bessel_trans(int l, vector<double> &fnr, vector<double> &fng) 
{
  // compute integral for radial Fourier transform of fn_l(r) using a Hankel transform:
  //       fn_l(G) = Ylm(G) (-i)^l \int r^2 j_l(Gr) fn_l(r) 4 pi dr
  // input:  fnr = fn_l(r)*4*pi*dr
  // output: fng = \int r^2 j_l(Gr) fnr
  //
  // j_0(Gr) = sin(Gr)/(Gr)
  // j_1(Gr) = sin(Gr)/(Gr)^2 - cos(Gr)/Gr
  // j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
  // j_3(Gr) = sin(Gr)*(15/(Gr)^4 - 6/(Gr)^2) - cos(Gr)*(15/(Gr)^3 - 1/Gr)
  // j_4(Gr) = sin(Gr)*(105/(Gr)^5 - 45/(Gr)^3 + 1/(Gr))
  //                            + cos(Gr)*(-105/(Gr)^4 + 10/(Gr)^2)
  // j_5(Gr) = sin(Gr)*(945/(Gr)^6 - 420/(Gr)^4 + 15/(Gr)^2)
  //                            + cos(Gr)*(-945/(Gr)^5 + 105/(Gr)^3 - 1/(Gr))
  // j_6(Gr) = sin(Gr)*(10395/(Gr)^7 - 4725/(Gr)^5 + 210/(Gr)^3 - 1/Gr)
  //                            + cos(Gr)*(-10395/(Gr)^6 + 1260/(Gr)^4 - 21/(Gr)^2)
  //
  // NOTE:  factor of Ylm(G)*(-i)^l is NOT included in output (fng)
  
  vector<double> fnint_(ndft_);
  double f0 = 0.0;   //  for l=1,2,3 f(G=0) is zero (j_l(Gr) -> 0 as G->0 )

  if ( l == 0 )
  {
    // s projector
    // w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r fnr
    //
    // G=0: Simpson integration up to ndft_
    // Use fnint_ as temporary array for integration

    for ( int i = 0; i < ndft_; i++ )
    {
      fnint_[i] = fnr[i] * rps_[i] * rps_[i];
    }
    f0 = simpsn(ndft_, &fnint_[0]);
    
    for ( int i = 0; i < ndft_; i++ )
    {
      fng[i] = fnr[i] * rps_[i];
    }
 
    sinft(&fng[0],ndft_);
    
    fng[0] = f0;
    // Divide by g
    for ( int i = 1; i < ndft_; i++ )
    {
      fng[i] /= gspl_[i];
    }
  }
  else if ( l == 1 )
  {
    //  p projectors
    //  w_1(G) = Ylm(G) i 1/G^2 \sum_r sin(Gr) fnr -
    //           Ylm(G) i 1/G   \sum_r cos(Gr) r fnr
    //  fnr(i,is,l) contains 4 pi dr phi_l(r) v_l(r)
    
    //  First part: 1/G^2 \sum_r sin(Gr) fnr
    for ( int i = 0; i < ndft_; i++ )
    {
      fng[i] = fnr[i];
    }

    sinft(&fng[0],ndft_);

    // Divide by g**2 and store in fnint_ */
    fnint_[0] = f0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fnint_[i] = fng[i] / ( gspl_[i] * gspl_[i] );
    }

    //  Second part: cosine transform: 1/G   \sum_r cos(Gr) r fnr
    
    for ( int i = 0; i < ndft_; i++ )
    {
      fng[i] = fnr[i] * rps_[i];
    }

    //  N.B. Next line: Initialize also fng[ndft_] to zero
    //  since it is used and modified by cosft1
    //  fng was dimensioned ndft_[is]+1
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);

    // Divide by g and add to fnint_ to get fng
    fng[0] = f0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fng[i] = fnint_[i] - fng[i]/gspl_[i];
    }
  }
  else if ( l == 2 )
  {
    // d projectors
    // d: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
    //    w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  fnr -
    //                            1/G   \sum_r sin(Gr)*r  fnr -
    //                            3/G^2 \sum_r cos(Gr)    fnr )
    // fnr[i] contains 4 pi dr phi_l(r) v_l(r)
    
    // First part: sine transform 3/G^3 \sum_r sin(Gr)/r fnr
    // Note: the integrand is linear near r=0 since fnr(r) ~ r^2
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fng[i] = fnr[i] / rps_[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by 3/G^3 and store in fnint_ */
    fnint_[0] = f0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fnint_[i] = 3.0 * fng[i] / ( gspl_[i] * gspl_[i] * gspl_[i] );
    }
    
    // Second part: sine transform -1/G \sum_r sin(Gr)*r fnr
    for ( int i = 0; i < ndft_; i++ )
    {
      fng[i] = fnr[i] * rps_[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by -1/G and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fnint_[i] += - fng[i] / gspl_[i];
    }
    
    // Third part: cosine transform: -3/G^2 \sum_r cos(Gr) fnr
    
    for ( int i = 0; i < ndft_; i++ )
    {
      fng[i] = fnr[i];
    }
    
    //  N.B. Next line: Initialize also fng[ndft_] to zero
    //  since it is used and modified by cosft1
    //  fng was dimensioned ndft_[is]+1
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -3/G^2 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fnint_[i] += - 3.0 * fng[i] / (gspl_[i] * gspl_[i]);
    }
    
    fng[0] = f0;
    for ( int i = 1; i < ndft_; i++ )
    {
      fng[i] = fnint_[i];
    }
    
  }
  else if ( l == 3 ) {
    // f projectors
    // f: j_3(Gr) = sin(Gr)*(15/(Gr)^4-6/(Gr)^2) - cos(Gr)*(15/(Gr)^3 - 1/(Gr))
    //    w_3(G) = Ylm(G) i^-3 (  15/G^4 \sum_r sin(Gr)/r^2 fnr -
    //                             6/G^2 \sum_r sin(Gr) fnr -
    //                            15/G^3 \sum_r cos(Gr)/r fnr + 
    //                             1/G   \sum_r cos(Gr)*r fnr )
    // fnr[i] contains 4 pi dr phi_l(r) v_l(r)
    
    // First part: sine transform 15/G^4 \sum_r sin(Gr)/r^2 fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i] * rps_[i]);
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by 15/G^4 and store in fnint_ */
    fnint_[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] = 15.0 * fng[i] / ( gspl_[i] * gspl_[i] * gspl_[i] * gspl_[i] );
    }
    
    // Second part: sine transform -6/G^2 \sum_r sin(Gr) fnr
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by -6/G^2 and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -6.0 * fng[i] / ( gspl_[i] * gspl_[i]);
    }
    
    // Third part: cosine transform: -15/G^3 \sum_r cos(Gr)/r fnr
    
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / rps_[i];
    }
    
    //  N.B. Next line: Initialize also fng[ndft_] to zero
    //  since it is used and modified by cosft1
    //  fng was dimensioned ndft_[is]+1
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -15/G^3 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += - 15.0 * fng[i] / (gspl_[i] * gspl_[i] * gspl_[i]);
    }
    
    // Fourth part: cosine transform: 1/G \sum_r cos(Gr)*r fnr
    
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i] * rps_[i];
    }
    
    //  N.B. Next line: Initialize also fng[ndft_] to zero
    //  since it is used and modified by cosft1
    //  fng was dimensioned ndft_[is]+1
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by 1/G and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += fng[i] / gspl_[i];
    }
    
    fng[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnint_[i];
    }
  }
  else if ( l == 4 ) {
     // j_4(Gr) = sin(Gr)*(105/(Gr)^5 - 45/(Gr)^3 + 1/(Gr))
     //                            + cos(Gr)*(-105/(Gr)^4 + 10/(Gr)^2)
    
    // sine transform 105/G^5 \sum_r sin(Gr)/r^3 fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]*rps_[i]);
    }
    sinft(&fng[0],ndft_);
    
    // multiply by 105/G^5 and store in fnint_ */
    fnint_[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] = 105.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // sine transform -45/G^3 \sum_r sin(Gr)/r fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / rps_[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by -45/G^3 and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -45.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // sine transform 1/G \sum_r r*sin(Gr) fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] * rps_[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by 1/G and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += fng[i]/(gspl_[i]);
    }
    
    // cosine transform: -105/G^4 \sum_r cos(Gr)/r^2 fnr
    
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i]/(rps_[i]*rps_[i]);
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -105/G^4 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -105.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // cosine transform: 10/G^2 \sum_r cos(Gr) fnr
    
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i];
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by 10/G^2 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += 10.*fng[i]/(gspl_[i]*gspl_[i]);
    }
    
    fng[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnint_[i];
    }
  }
  else if ( l == 5 ) {
     // j_5(Gr) = sin(Gr)*(945/(Gr)^6 - 420/(Gr)^4 + 15/(Gr)^2)
     //                            + cos(Gr)*(-945/(Gr)^5 + 105/(Gr)^3 - 1/(Gr))
    
    // sine transform 945/G^6 \sum_r sin(Gr)/r^4 fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]*rps_[i]*rps_[i]);
    }
    sinft(&fng[0],ndft_);
    
    // multiply by 945/G^6 and store in fnint_ */
    fnint_[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] = 945.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // sine transform -420/G^4 \sum_r sin(Gr)/r^2 fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]);
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by -420/G^4 and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -420.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // sine transform 15/G^2 \sum_r sin(Gr) fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by 15/G^2 and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += 15.*fng[i]/(gspl_[i]*gspl_[i]);
    }
    
    // cosine transform: -945/G^5 \sum_r cos(Gr)/r^3 fnr
    
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]*rps_[i]);
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -945/G^5 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -945.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // cosine transform: 105/G^3 \sum_r cos(Gr)/r fnr
    
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i]/rps_[i];
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by 105/G^3 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += 105.*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // cosine transform: -1/G \sum_r r*cos(Gr) fnr
    
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i]*rps_[i];
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -1/G and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -1.*fng[i]/(gspl_[i]);
    }
    
    fng[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnint_[i];
    }
  }
  else if ( l == 6 ) {
     // j_6(Gr) = sin(Gr)*(10395/(Gr)^7 - 4725/(Gr)^5 + 210/(Gr)^3 - 1/Gr)
     //                            + cos(Gr)*(-10395/(Gr)^6 + 1260/(Gr)^4 - 21/(Gr)^2)
    
    // sine transform 10395/G^7 \sum_r sin(Gr)/r^5 fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]*rps_[i]*rps_[i]*rps_[i]);
    }
    sinft(&fng[0],ndft_);
    
    // multiply by 10395/G^7 and store in fnint_ */
    fnint_[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] = 10395.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // sine transform -4725/G^5 \sum_r sin(Gr)/r^3 fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]*rps_[i]);
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by -4725/G^5 and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -4725.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // sine transform 210/G^3 \sum_r sin(Gr)/r fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i]/rps_[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by 210/G^3 and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += 210.*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]);
    }

    // sine transform -1/G \sum_r r*sin(Gr) fnr
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i]*rps_[i];
    }
    
    sinft(&fng[0],ndft_);
    
    // multiply by -1/G and accumulate in fnint_ */
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] -= fng[i]/gspl_[i];
    }

    // cosine transform: -10395/G^6 \sum_r cos(Gr)/r^4 fnr
    
    fng[0] = 0.0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnr[i] / (rps_[i]*rps_[i]*rps_[i]*rps_[i]);
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -10395/G^6 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -10395.0*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // cosine transform: 1260/G^4 \sum_r cos(Gr)/r^2 fnr
    
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i]/(rps_[i]*rps_[i]);
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by 1260/G^4 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += 1260.*fng[i]/(gspl_[i]*gspl_[i]*gspl_[i]*gspl_[i]);
    }
    
    // cosine transform: -21/G^2 \sum_r cos(Gr) fnr
    
    for ( int i = 0; i < ndft_; i++ ) {
      fng[i] = fnr[i];
    }
    
    fng[ndft_] = 0.0;
    cosft1(&fng[0],ndft_);
    
    // Multiply by -21/G^2 and add to fnint_
    fnint_[0] += f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fnint_[i] += -21.*fng[i]/(gspl_[i]*gspl_[i]);
    }
    
    fng[0] = f0;
    for ( int i = 1; i < ndft_; i++ ) {
      fng[i] = fnint_[i];
    }
  }


}
////////////////////////////////////////////////////////////////////////////////
double Species::ylm(int l, int m, double gx, double gy, double gz)
{
  const double pi = M_PI;
  const double fpi = 4.0 * pi;

  double f;
  double gnorm = sqrt(gx*gx+gy*gy+gz*gz);
  double gi = 1./gnorm;
  if (gnorm == 0.0) gi = 0.0;
  double gxi = gx*gi;
  double gyi = gy*gi;
  double gzi = gz*gi;
  double gxi2 = gxi*gxi;
  double gyi2 = gyi*gyi;
  double gzi2 = gzi*gzi;
  
  // compute real spherical harmonics 
  if (l == 0) {
    const double s14pi = sqrt(1.0/fpi);
    f = s14pi;
  }
  else if (l == 1) {
    const double s34pi = sqrt(3.0/fpi);  
    if (m == 0) f = s34pi*gzi;
    else if (m == 1) f = s34pi*gxi;
    else if (m == -1) f = s34pi*gyi;
  }
  else if (l == 2) {
    const double s54pi = sqrt(5.0/fpi);
    const double s3 = sqrt(3.0);
    if (m == 0) f = 0.5*s54pi*(3.*gzi2-1.);
    else if (m == 1) f = s3*s54pi*gxi*gzi;
    else if (m == -1) f = s3*s54pi*gyi*gzi;
    else if (m == 2) f = 0.5*s3*s54pi*(gxi2-gyi2);
    else if (m == -2) f = s3*s54pi*gxi*gyi;
  }
  else if (l == 3) {
    const double s74pi = sqrt(7.0/fpi);
    const double s2132pi = sqrt(21.0/(32.*pi));
    const double s3532pi = sqrt(35.0/(32.*pi));
    const double s1054pi = sqrt(105.0/fpi);
    if (m == 0) f = 0.5*s74pi*gzi*(5.*gzi2 - 3.);
    else if (m == 1) f = s2132pi*gxi*(5.*gzi2-1.);
    else if (m == -1) f = s2132pi*gyi*(5.*gzi2-1.);
    else if (m == 2) f = 0.5*s1054pi*gzi*(gxi2-gyi2);
    else if (m == -2) f = s1054pi*gxi*gyi*gzi;
    else if (m == 3) f = s3532pi*gxi*(gxi2-3.*gyi2);
    else if (m == -3) f = s3532pi*gyi*(3.*gxi2-gyi2);
  }
  else if (l == 4) {
    double gxi3 = gxi2*gxi;
    double gyi3 = gyi2*gyi;
    double gzi3 = gzi2*gzi;
    const double s14pi = sqrt(1.0/fpi);
    const double s52pi = sqrt(10.0/fpi);
    const double s54pi = sqrt(5.0/fpi);
    const double s351pi = sqrt(35.0/pi);
    const double s3532pi = sqrt(35.0/(32.*pi));
    if (m == 0) f = 0.375*s14pi*(35.*gzi2*gzi2 - 30.*gzi2 + 3.);
    else if (m == 1) f = 0.75*s52pi*gxi*gzi*(7.*gzi2 - 3.);
    else if (m == -1) f = 0.75*s52pi*gyi*gzi*(7.*gzi2 - 3.);
    else if (m == 2) f = 0.1875*s54pi*(gxi2-gyi2)*(7.*gzi2 - 1.);
    else if (m == -2) f = 0.375*s54pi*gxi*gyi*(7.*gzi2 - 1.);
    else if (m == 3) f = 0.1875*s3532pi*gxi*gzi*(gxi2 - 3.*gyi2);
    else if (m == -3) f = 0.1875*s3532pi*gyi*gzi*(3.*gxi2 - gyi2);
    else if (m == 4) f = 0.1875*s351pi*gxi2*(gxi2-3.*gyi2) - gyi2*(3.*gxi2-gyi2);
    else if (m == -4) f = 0.75*s351pi*gxi*gyi*(gxi2 - gyi2);
  }
  else if (l == 5) {
    double gxi3 = gxi2*gxi;
    double gyi3 = gyi2*gyi;
    double gzi3 = gzi2*gzi;
    const double s11pi = sqrt(11.0/pi);
    const double s77pi = sqrt(77.0/pi);
    const double s1654pi = sqrt(165.0/fpi);
    const double s3852pi = sqrt(385.0/(2.*pi));
    const double s3854pi = sqrt(385.0/fpi);
    const double s11554pi = sqrt(1155.0/fpi);
    if (m == 0) f = 0.0625*s11pi*(63.*gzi3*gzi2-70.*gzi3+15.*gzi);
    else if (m == 1) f = -0.125*s1654pi*gxi*(21.*gzi2*gzi2-14.*gzi2+1.);
    else if (m == -1) f = -0.125*s1654pi*gyi*(21.*gzi2*gzi2-14.*gzi2+1.);
    else if (m == 2) f = 0.125*s11554pi*(gxi2-gyi2)*(3.*gzi3-gzi);
    else if (m == -2) f = 0.5*s11554pi*gxi*gyi*(3.*gzi3-gzi);
    else if (m == 3) f = -0.0625*s3852pi*(gxi3-3.*gxi*gyi2)*(9.*gzi2-1.);
    else if (m == -3) f = -0.0625*s3852pi*(3.*gxi2*gyi-gyi3)*(9.*gzi2-1.);
    else if (m == 4) f = 0.375*s3854pi*gzi*(gxi2*gxi2-6.*gxi2*gyi2+gyi2*gyi2);
    else if (m == -4) f = 1.5*s3854pi*gzi*(gxi3*gyi-gxi*gyi3);
    else if (m == 5) f = -0.1875*s77pi*(gxi3*gxi2-10.*gxi3*gyi2+gxi*gyi2*gyi2);
    else if (m == -5) f = -0.1875*s77pi*(5.*gxi2*gxi2*gyi-10.*gxi2*gyi3+gyi3*gyi2);
  }
  else if (l == 6) {
    double gxi3 = gxi2*gxi;
    double gyi3 = gyi2*gyi;
    double gzi3 = gzi2*gzi;
    const double s13pi = sqrt(13./pi);
    const double s2734pi = sqrt(273./fpi);
    const double s10012pi = sqrt(1001.*0.5/pi);
    const double s13652pi = sqrt(1365.*0.5/pi);
    const double s30032pi = sqrt(3003.*0.5/pi);
    const double s914pi = sqrt(91./fpi);
    if (m == 0) f = 0.03125*s13pi*(231.*gzi3*gzi3-315.*gzi2*gzi2+105.*gzi2-5.);
    else if (m == 1) f = -0.125*s2734pi*gxi*(33.*gzi3*gzi2-30.*gzi3+5.*gzi);
    else if (m == -1) f = -0.125*s2734pi*gyi*(33.*gzi3*gzi2-30.*gzi3+5.*gzi);
    else if (m == 2) f = 0.015625*s13652pi*(gxi2-gyi2)*(33.*gzi2*gzi2-18.*gzi2+1.);
    else if (m == -2) f = 0.0625*s13652pi*gxi*gyi*(33.*gzi2*gzi2-18.*gzi2+1.);
    else if (m == 3) f = -0.0625*s13652pi*(gxi3-3.*gxi*gyi2)*(11.*gzi3-3.*gzi);
    else if (m == -3) f = -0.0625*s13652pi*(3.*gxi2*gyi-gyi3)*(11.*gzi3-3.*gzi);
    else if (m == 4) f = 0.1875*s914pi*(gxi2*gxi2-6.*gxi2*gyi2+gyi2*gyi2)*(11.*gzi2-1.);
    else if (m == -4) f = 0.375*s914pi*(gxi3*gyi-gxi*gyi3)*(11.*gzi2-1.);
    else if (m == 5) f = -0.1875*s10012pi*(gxi3*gxi2-10.*gxi3*gyi2+5.*gxi*gyi2*gyi2)*gzi;
    else if (m == -5) f = -0.1875*s10012pi*(5.*gxi2*gxi2*gyi-10.*gxi2*gyi3+gyi3*gyi2)*gzi;
    else if (m == 6) f = 0.03125*s30032pi*(gxi3*gxi3-15.*gxi2*gxi2*gyi2+15.*gxi2*gyi2*gyi2-gyi3*gyi3);
    else if (m == -6) f = 0.03125*s30032pi*(6.*gxi3*gxi2*gyi-20.*gxi3*gyi3+6.*gxi*gyi3*gyi2);
  }
  return f;
}
////////////////////////////////////////////////////////////////////////////////
void Species::vpsr(int l, double r, double &v)
{
  if ( l > lmax_ || r > rps_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(&rps_[0],&vps_[l][0],&vps_spl_[l][0],ndft_,r,&v);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::dvpsr(int l, double r, double &v, double &dv)
{
  if ( l > lmax_ || r > rps_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(&rps_[0],&vps_[l][0],&vps_spl_[l][0],ndft_,r,&v,&dv);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::vlocg(double g, double &v)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(&gspl_[0],&vlocg_[0],&vlocg_spl[0],ndft_,g,&v);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::dvlocg(double g, double &v, double &dv)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(&gspl_[0],&vlocg_[0],&vlocg_spl[0],ndft_,g,&v,&dv);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::vnlg(int l, int ic, double g, double &v)
{
  assert ( l >= 0 && l <= lmax_ );
  if ( l == llocal_ || g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(&gspl_[0], &vnlg_[l][ic][0], &vnlg_spl[l][ic][0], ndft_, g, &v);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::dvnlg(int l, int ic, double g, double &v, double &dv)
{
  assert ( l >= 0 && l <= lmax_ );
  if ( l == llocal_ || g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(&gspl_[0], &vnlg_[l][ic][0], &vnlg_spl[l][ic][0], ndft_, g, &v, &dv);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::phig(double g, double &v)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(&gspl_[0],&phig_[0],&phig_spl_[0],ndft_,g,&v);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::dphig(double g, double &v, double &dv)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(&gspl_[0],&phig_[0],&phig_spl_[0],ndft_,g,&v,&dv);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::betag(int b, double g, double &v)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(&gspl_[0],&betag_[b][0],&betag_spl_[b][0],ndft_,g,&v);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Species::dbetag(int b, double g, double &v, double &dv)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(&gspl_[0],&betag_[b][0],&betag_spl_[b][0],ndft_,g,&v,&dv);
  }
}
////////////////////////////////////////////////////////////////////////////////
double Species::qfung(int q, int ltot, double g)
{
  //need to subtract abs(qnm_l1-qnm_l2) from ltot to get function index
  int l = ltot - abs(qfunl1_[q]-qfunl2_[q]);  // q  must be 0...nqfun, NOT nqtot
  assert(l >= 0);

  //ewd DEBUG
  if (q >= nqtot_)
    cout << "ERROR:  Species.qfung, q = " << q << ", nqtot_ = " << nqtot_ << ", ltot = " << ltot << endl;
  
  double v;
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(&gspl_[0],&qfung_[q][l][0],&qfung_spl_[q][l][0],ndft_,g,&v);
  }
  return v;
}

////////////////////////////////////////////////////////////////////////////////
double Species::rhog_nlcc(double g)
{
  double rho;
  if ( g > gspl_[ndft_-1] )
  {
    rho = 0.0;
  }
  else
  {
    splint(&gspl_[0],&rhog_nlcc_[0],&rhog_nlcc_spl_[0],ndft_,g,&rho);
  }
  return rho;
}

////////////////////////////////////////////////////////////////////////////////
double Species::rhopsg( double g )
{
  double arg = 0.25 * rcps_ * rcps_ * g * g;
  return -zval_ * exp( -arg );
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator << ( ostream &os, Species &s )
{
  // XML output of species
  // If the uri is known, use href to refer to it
  if ( s.uri() != "" )
  {
    os <<"<species name=\"" << s.name() 
       << "\" href=\"" << s.uri() << "\"/>" << endl;
  }
  else
  {
    os <<"<species name=\"" << s.name() << "\">" << endl;
    os << "<description>" << s.description() << "</description>" << endl;
    os << "<symbol>" << s.symbol() << "</symbol>" << endl;
    os << "<atomic_number>" << s.atomic_number() << "</atomic_number>" << endl;
    os << "<mass>" << s.mass() << "</mass>" << endl;
    os << "<norm_conserving_pseudopotential>" << endl;
    os << "<valence_charge>" << s.zval() << "</valence_charge>" << endl;
    os << "<lmax>" << s.lmax() << "</lmax>" << endl;
    os << "<llocal>" << s.llocal() << "</llocal>" << endl;
    os << "<nquad>" << s.nquad() << "</nquad>" << endl;
    os << "<rquad>" << s.rquad() << "</rquad>" << endl;
    os << "<mesh_spacing>" << s.deltar() << "</mesh_spacing>" << endl;
    os.setf(ios::fixed,ios::floatfield);
    os << setprecision(6);
    for ( int l = 0; l <= s.lmax(); l++ )
    {
      const int size = s.vps()[l].size();
      os << "<projector l=\"" << l << "\" size=\"" << size 
         << "\">" << endl;
      os << "<radial_potential>\n";
      for ( int i = 0; i < size; i++ )
        os << s.vps()[l][i] << "\n";
      os << "</radial_potential>\n";
      if ( l < s.phi().size() && s.phi()[l].size() == size )
      {
        os << "<radial_function>\n";
        for ( int i = 0; i < size; i++ )
          os << s.phi()[l][i] << "\n";
        os << "</radial_function>\n";
      }
      os << "</projector>" << endl;
    }
    os << "</norm_conserving_pseudopotential>" << endl;
    os << "</species>" << endl;
  }
  
  return os;
}

////////////////////////////////////////////////////////////////////////////////
void Species::printsys(ostream& os) const {
  os.setf(ios::left,ios::adjustfield);
  if (fix_rcps)
    os << "species " << name() << " " << uri() << " " << rcps() << endl;
  else
    os << "species " << name() << " " << uri() << endl;
}

////////////////////////////////////////////////////////////////////////////////
void Species::info(ostream &os)
{
  os.setf(ios::left,ios::adjustfield);

  os << " name_ = " << name() << endl;
  os << " description_ = " << description() << endl;
  os << " uri_ = " << uri() << endl;
  os << " symbol_ = " << symbol() << endl;
  os << " atomic_number_ = " << atomic_number() << endl;
  // describe type of potential
  if (usoft_) {
    os << " Vanderbilt ultrasoft potential" << endl;
  }
  else {
    if ( nquad() == 0 )
    {
      if ( lmax() == 0 )
        os << " local potential" << endl;
      else
        os << " Kleinman-Bylander potential" << endl;
    }
    else
    {
      os << " Semi-local potential with " << nquad()
         << " quadrature points in [0.0, " << rquad() << "]" << endl;
      os << " local (within 1.e-6) beyond r = " << rcut_loc(1.e-6) << endl;
      os << " local (within 1.e-5) beyond r = " << rcut_loc(1.e-5) << endl;
      os << " local (within 1.e-4) beyond r = " << rcut_loc(1.e-4) << endl;
      os << " local (within 1.e-3) beyond r = " << rcut_loc(1.e-3) << endl;
    }
  }
  
  os << " valence charge = " << zval()
     << " / ionic mass_ = " << mass()
     << " (amu)" << endl;
  os << " lmax_ =   " << lmax() << endl;
  os << " llocal_ = " << llocal() << endl;
  os << " rcps_ =   " << rcps() << endl; 
}

////////////////////////////////////////////////////////////////////////////////
double Species::rcut_loc(double epsilon)
{
  // find radius at which the largest deviation delta_vnl(r) is < epsilon
  double delta = 0.0;
  int i = ndft_-1;
  while ( ( delta < epsilon ) && i > 0 )
  {
    i--;
    for ( int l = 0; l <= lmax_; l++ )
    {
      if ( l != llocal_ )
      {
        // compute deviation vps_[l][i]
        double dv = fabs( vps_[l][i] - vps_[llocal_][i] );
        delta = dv > delta ? dv : delta;
      }
    }
  }
  // adjust i so that delta_v[i] < epsilon
  if ( i < ndft_-1 ) i++;
  
  return rps_[i];
}
////////////////////////////////////////////////////////////////////////////////
void Species::set_hubbard_u(double uval, int lval)
{
  hubbard_u_ = uval;
  hubbard_l_ = lval;

  if (lval > lmax_) 
    if (ctxt_.oncoutpe()) 
      cout << "<WARNING> Species::set_hubbard_u: hubbard_l > lmax for species " << name_ << "! </WARNING>" << endl;
    
  if (lval == llocal_)
  {
    if (ctxt_.oncoutpe()) 
      cout << "<ERROR> Species::set_hubbard_u: hubbard_l = llocal for species " << name_ << "! This is not currently supported.</ERROR>" << endl;
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Species::substract_long_range_part(const vector<double> & vloc, vector<double> & vloc_sr) const {
  // local potential: subtract the long range part due to the smeared charge
  // Next line: constant is 2/sqrt(pi)
  // math.h: # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
  
  vloc_sr[0] = vloc[0] + (zval_/rcps_)*M_2_SQRTPI;
  for (int ip = 1; ip < ndft_; ip++){
    vloc_sr[ip] = vloc[ip] + (zval_/rps_[ip]) * erf( rps_[ip]/rcps_ );
  }
  
}
