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
// EnergyFunctional.C
//
////////////////////////////////////////////////////////////////////////////////

#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "StructureFactor.h"
#include "XCPotential.h"
#include "NonLocalPotential.h"
#include "ConfinementPotential.h"
#include "EnthalpyFunctional.h"
#include "HubbardPotential.h"

#include "Timer.h"
#include "blas.h"

#include <iostream>
#include <iomanip>
#include <algorithm> // fill()
using namespace std;

////////////////////////////////////////////////////////////////////////////////
EnergyFunctional::EnergyFunctional(const Sample& s, const Wavefunction& wf, ChargeDensity& cd)
    : s_(s), wf_(wf), cd_(cd) {
  const AtomSet& atoms = s_.atoms;
  
  const bool compute_stress = ( s_.ctrl.stress == "ON" );  // if stress off, don't store dtwnl

  sigma_ekin.resize(6);
  sigma_econf.resize(6);
  sigma_eps.resize(6);
  sigma_ehart.resize(6);
  sigma_exc.resize(6);
  sigma_enl.resize(6);
  sigma_esr.resize(6);
  sigma.resize(6);
  
  vbasis_ = cd_.vbasis();
  //cout << vbasis_->context().mype() << ": vbasis_->context() = " 
  //     << vbasis_->context() << endl;
  
  // define FT's on vbasis contexts
  
  vft = cd_.vft();
  int np0v = vft->np0();
  int np1v = vft->np1();
  int np2v = vft->np2();
  
  v_r.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    v_r[ispin].resize(vft->np012loc());
  }
  tmp_r.resize(vft->np012loc());

  if ( s_.ctxt_.oncoutpe() ) {
    cout << "  <!-- EnergyFunctional: charge density basis: " << vbasis_->size() << " plane waves, ngloc = " << vbasis_->localsize() << " -->" << endl;
    cout << "  <!-- EnergyFunctional: np0v,np1v,np2v: " << np0v << " "
         << np1v << " " << np2v << " -->" << endl;
    cout << "  <!-- EnergyFunctional: vft->np012(): "
         << vft->np012() << " -->" << endl;
  }
  
  const int ngloc = vbasis_->localsize();
  //cout << " EnergyFunctional: ngloc: " << ngloc << endl;
  
  nsp_ = atoms.nsp();
  
  vps.resize(nsp_);
  dvps.resize(nsp_);
  rhops.resize(nsp_);
  
  zv_.resize(nsp_);
  rcps_.resize(nsp_);
  na_.resize(nsp_);
  namax_ = 0;
  
  for ( int is = 0; is < nsp_; is++ )
  {
    vps[is].resize(ngloc);
    dvps[is].resize(ngloc);
    rhops[is].resize(ngloc);
    if ( atoms.na(is) > namax_ ) namax_ = atoms.na(is);
  }

  // AS: compute total electronic density used for setting up the Hamiltonian
  if (s_.ctrl.tddft_involved)
  {
     hamil_cd_ = new ChargeDensity(s_,*s_.hamil_wf);
     hamil_rhoelg.resize(ngloc);
     hamil_rhogt.resize(ngloc);
     
     (*hamil_cd_).update_density();
     update_hamiltonian();

     // AS: the charge density based on hamil_wf has to be used
     xcp = new XCPotential((*hamil_cd_),s_.ctrl.xc,cd_);
  }
  else
  {
     xcp = new XCPotential(cd_,s_.ctrl.xc);
  }

  nlp.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    nlp[ispin].resize(wf_.nkptloc());
  
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for (int kloc=0; kloc<wf_.nkptloc(); kloc++) 
      nlp[ispin][kloc] = 0;

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    if (wf_.spinactive(ispin)) 
      for (int kloc=0; kloc<wf_.nkptloc(); kloc++) 
         nlp[ispin][kloc] = new NonLocalPotential((AtomSet&)s_.atoms, *wf_.sdloc(ispin,kloc),compute_stress);

  if (s_.ctrl.extra_memory >= 5) // use extra memory in large, huge mode
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
      if (wf_.spinactive(ispin)) 
        for (int kloc=0; kloc<wf_.nkptloc(); kloc++) 
          nlp[ispin][kloc]->use_highmem();
  
  if (s_.ctrl.dft_plus_u) 
    hubp_ = new HubbardPotential((AtomSet&)s_.atoms, s_.symmetries, wf_);
  
  vion_local_g.resize(ngloc);
  dvion_local_g.resize(ngloc);
  vlocal_g.resize(ngloc);
  vtemp.resize(ngloc);
  rhoelg.resize(ngloc);
  rhogt.resize(ngloc);
  rhopst.resize(ngloc);

  veff_g.resize(wf_.nspin());
  if (s_.ctrl.ultrasoft) {
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
      veff_g[ispin].resize(ngloc);
  }
  vxc_g.resize(wf_.nspin());
  if (s_.ctrl.nlcc || s_.ctrl.ultrasoft) {
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
      vxc_g[ispin].resize(ngloc);
  }
  
  tau0.resize(nsp_);
  taum.resize(nsp_);
  fion_esr.resize(nsp_);
  ftmp.resize(3*namax_);
  
  eself_ = 0.0;

  // set Ewald parameters to ensure convergence of sum and avoid finite-size errors
  s_.atoms.set_rcps(wf_.ecut());

  for ( int is = 0; is < nsp_; is++ ) {
    Species *s = atoms.species_list[is];
    
    const int na = atoms.na(is);
    tau0[is].resize(3*na);
    taum[is].resize(3*na);
    fion_esr[is].resize(3*na);
    
    eself_ += na * s->eself();
    na_[is] = na;
    
    zv_[is] = s->zval();
    rcps_[is] = s->rcps();
  }
  

  // FT for interpolation of wavefunctions on the fine grid
  ft.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    ft[ispin].resize(wf_.nkp());

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ ) 
      ft[ispin][ikp] = cd_.ft(ispin,ikp);

  // Confinement potentials
  cfp.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    cfp[ispin].resize(wf_.nkp());

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      cfp[ispin][ikp] = 0;

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          cfp[ispin][ikp] = new ConfinementPotential(s_.ctrl.ecuts,s_.ctrl.facs,
                                            s_.ctrl.sigmas,wf_.sd(ispin,ikp)->basis());
        }
      }
    }
  }
  
  sf.init(tau0,*vbasis_);
  
  cell_moved(compute_stress);  //ewd:  compute_stress = false, don't store dtwnl
  
  atoms_moved();

  epvf = 0;
  if (s_.ctrl.enthalpy_pressure != 0.0) 
    epvf = new EnthalpyFunctional(cd_,s_.ctrl.enthalpy_pressure,s_.ctrl.enthalpy_threshold);

  eharris_ = 0.0;
  
}

////////////////////////////////////////////////////////////////////////////////
EnergyFunctional::~EnergyFunctional(void) {
  delete xcp;
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    if (wf_.spinactive(ispin)) {
      int nlpsize = nlp[ispin].size();
      for (int i=0; i<nlpsize; i++) 
        delete nlp[ispin][i];
    }
  
  if (s_.ctrl.dft_plus_u)
    delete hubp_;
  
  if (epvf != 0)
    delete epvf;

  if (s_.ctrl.tddft_involved)
     delete hamil_cd_;
  hamil_cd_ = 0;
}
////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::update_hamiltonian(void)
// ewd:  updates hamil_rhoelg from hamil_cd_
{
   const Wavefunction& wf = s_.wf;
   const int ngloc = vbasis_->localsize();
   const UnitCell& cell = wf.cell();
   const double omega = cell.volume();
   const double omega_inv = 1.0 / omega;

   if ( wf.nspin() == 1 )
   {
      for ( int ig = 0; ig < ngloc; ig++ )
      {
         hamil_rhoelg[ig] = omega_inv * (*hamil_cd_).rhog[0][ig];
      }
   }
   else
   {
      for ( int ig = 0; ig < ngloc; ig++ )
      {
         hamil_rhoelg[ig] = omega_inv * ( (*hamil_cd_).rhog[0][ig] + (*hamil_cd_).rhog[1][ig] );
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::print_timing() {
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    uint64_t count = (*i).second.counts();
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.mype()==0 ) {
       cout << left << setw(34) << "<timing where=\"energy_functional\""
            << setw(8) << " name=\""
            << setw(15) << (*i).first << "\""
            << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
            << " max=\"" << setprecision(3) << setw(9) << tmax << "\""
            << " count=\"" << setw(9) << count << "\"/>"
            << endl;
    }
  }
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
     if (wf_.spinactive(ispin)) {
        int nlpsize = nlp[ispin].size();
        for (int i=0; i<nlpsize; i++) 
           nlp[ispin][i]->print_timing();
     }

}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::update_vhxc(void) {
  // called when the charge density has changed
  // update Hartree and xc potentials using the charge density cd_
  // compute Hartree and xc energies

  const UnitCell& cell = wf_.cell();
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const double *const g2i = vbasis_->g2i_ptr();
  const double fpi = 4.0 * M_PI;
  const int ngloc = vbasis_->localsize();
  double tsum[2];
  tsum[0] = 0.0;
  tsum[1] = 0.0;
  
  // compute total electronic density: rhoelg = rho_up + rho_dn
  if ( wf_.nspin() == 1 ) {
    for ( int ig = 0; ig < ngloc; ig++ ) {
      rhoelg[ig] = omega_inv * cd_.rhog[0][ig];
    }
  }
  else {
    for ( int ig = 0; ig < ngloc; ig++ ) {
      rhoelg[ig] = omega_inv * ( cd_.rhog[0][ig] + cd_.rhog[1][ig] );
    }
  }
  
  // update XC energy and potential
  tmap["exc"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for (int i=0; i<vft->np012loc(); i++)
      v_r[ispin][i] = 0.0;
  
  //fill(v_r[ispin].begin(),v_r[ispin].end(),0.0);

  xcp->update(v_r);
  exc_ = xcp->exc();
  tmap["exc"].stop();


  /*
  //ewd DEBUG: calculate integral of (vxc(r)+vhart(r))*rhor(r) to compare w. PWSCF
  double tvxc = 0.0;
  for (int i=0; i<vft->np012loc(); i++) 
    tvxc += -1.*v_r[0][i]*cd_.rhor[0][i];
  double tfac = omega/(double)vft->np012loc();
  tvxc *= tfac;

  double thart = 0.0;
  double teh = 0.0;
  for ( int ig = 0; ig < ngloc; ig++ ) {
    vlocal_g[ig] = fpi * rhoelg[ig] * g2i[ig];
    teh += norm(rhoelg[ig])*g2i[ig];
  }
  teh *= omega*fpi*0.5;
  cout << "EF.EHART: ehart = " << teh << endl;
  
  vft->backward(&vlocal_g[0],&tmp_r[0]);
  for (int i=0; i<vft->np012loc(); i++) 
    thart += -1.*real(tmp_r[i])*cd_.rhor[0][i];
  thart *= tfac;
  cout << "EF.DEBAND: " << thart+tvxc << ", thart = " << thart << ", tvxc = " << tvxc << endl;
  */
  //ewd DEBUG

  //ewd DEBUG
  //for (int i=0; i<vft->np012loc(); i++) 
  //  cout << "EF.RHOR mype = " << s_.ctxt_.mype() << ", ir = " << i << ", rhor = " << cd_.rhor[0][i] << endl;
  //for ( int ig = 0; ig < ngloc; ig++ ) 
  //  cout << "EF.RHOG mype = " << s_.ctxt_.mype() << ", ig = " << ig << ", rhog = " << cd_.rhog[0][ig] << endl;

  
  
  // we need xc potential in reciprocal space for ultrasoft
  if (s_.ctrl.ultrasoft || s_.ctrl.nlcc) {
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
      vector<complex<double> > vrtmp(vft->np012loc());
      for (int i=0; i<vft->np012loc(); i++) 
        vrtmp[i] = complex<double>(v_r[ispin][i], 0.0);

      vft->forward(&vrtmp[0],&vxc_g[ispin][0]);    
    }
  }
  
  // update electronic enthalpy energy and potential
  epv_ = 0.0;
  if (s_.ctrl.enthalpy_pressure != 0.0) {
    epvf->update(v_r);
    epv_ = epvf->epv();
    const double gpa = 29421.0120;
    if (s_.ctxt_.oncoutpe()) {
      cout << "<!-- Enthalpy functional:  pressure = " << gpa*s_.ctrl.enthalpy_pressure << " GPa, threshold = " << s_.ctrl.enthalpy_threshold << " -->" << endl; 
      cout << setprecision(8);
      cout << "<electronic_volume> " << epvf->evol() << " </electronic_volume>" << endl;
    }
  }
  
  // compute local potential energy: 
  // integral of el. charge times ionic local pot.
  int len=2*ngloc,inc1=1;
  if (vbasis_->real()) { 
    tsum[0] = 2.0 * ddot(&len,(double*)&rhoelg[0],&inc1,
                         (double*)&vion_local_g[0],&inc1);
    // remove double counting for G=0
    if ( vbasis_->context().myrow() == 0 ) {
      tsum[0] -= real(conj(rhoelg[0])*vion_local_g[0]);
    }
  }
  else {
    tsum[0] = ddot(&len,(double*)&rhoelg[0],&inc1,
                         (double*)&vion_local_g[0],&inc1);
  }
  tsum[0] *= omega; // tsum[0] contains eps
  
  // Hartree energy
  ehart_ = 0.0;

  if ( s_.ctrl.esm_bc == "" )
  {
     double ehsum = 0.0;
     if (s_.ctrl.tddft_involved)
     {
        // AS: the individual terms of the sum are complex in general
        complex<double> ehsum_cplx = complex<double>(0.0,0.0);
        for ( int ig = 0; ig < ngloc; ig++ )
        {
           const complex<double> tmp = rhoelg[ig] + rhopst[ig];
           const complex<double> hamil_tmp = hamil_rhoelg[ig] + rhopst[ig];
           ehsum_cplx += tmp * conj(hamil_tmp) * complex<double>(g2i[ig],0.0);
           // AS: the next line is needed for correct stress and/or forces
           rhogt[ig] = tmp;
           hamil_rhogt[ig] = hamil_tmp;
        }
        ehsum = real(ehsum_cplx);     
     }
     else
     {
        for ( int ig = 0; ig < ngloc; ig++ ) {
           const complex<double> tmp = rhoelg[ig] + rhopst[ig];
           ehsum += norm(tmp) * g2i[ig];
           rhogt[ig] = tmp;
        }
     }
     
     // factor 1/2 from definition of Ehart cancels with half sum over G
     // Note: rhogt[ig] includes a factor 1/Omega
     // Factor omega in next line yields prefactor 4 pi / omega in
     
     double vfact = vbasis_->real() ? 1.0 : 0.5;
     tsum[1] = vfact * omega * fpi * ehsum;
     // tsum[1] contains ehart
  
     vbasis_->context().dsum(2,1,&tsum[0],2);
     eps_   = tsum[0];
     ehart_ = tsum[1];
  
     // compute vlocal_g = vion_local_g + vhart_g
     // where vhart_g = 4 * pi * (rhoelg + rhopst) * g2i  
     if (s_.ctrl.tddft_involved)  // AS: the charge density based on hamil_wf has to be used
        for ( int ig = 0; ig < ngloc; ig++ )
           vlocal_g[ig] = vion_local_g[ig] + fpi * hamil_rhogt[ig] * g2i[ig];
     else
        for ( int ig = 0; ig < ngloc; ig++ )
           vlocal_g[ig] = vion_local_g[ig] + fpi * rhogt[ig] * g2i[ig];

  }
  else
  {
    // else replace with ESM hartree energy and potential //
     
     assert(!s_.ctrl.tddft_involved);  // I haven't checked that this works with TDDFT

     const string esm_bc = s_.ctrl.esm_bc; // boundary conditions (bc1, bc2, or bc3) ( default "" )
     const double esm_w = s_.ctrl.esm_w;   // position offset of ESM region ( default 0.0 )
     const int esm_nfit = 4;               // number of fitting parameters/points for smoothing
                                          // shouldn't need to change; can probably hard-code
     const Basis& basis = *vbasis_;
     const int np2v = vft->np2();
     const double tpi = 2.0 * M_PI;
     const complex<double> ci( 0.0, 1.0 );
     const complex<double> c0( 0.0, 0.0 );
     const double L = cell.amat(8);
     const double z0 = L/2.0;
     const double z1 = z0 + fabs(esm_w);
     const int np2half = np2v/2;
     vector<complex<double> > vhart_g( ngloc, c0 );

     for ( int ig = 0; ig < ngloc; ig++ )
     {
        rhogt[ig] = rhoelg[ig] + rhopst[ig];
     }

#if _OPENMP
#pragma omp parallel for
#endif     
     for ( int irod = 0; irod < basis.nrod_loc(); irod++ )
     {
        vector<complex<double> > vg2 ( np2v, c0 );
        vector<complex<double> > vg2b( np2v, c0 );

        int k0 = basis.rod_h(irod);
        int k1 = basis.rod_k(irod);
//
// if g_parallel != 0:
//
        if ( ( k0 != 0 ) || ( k1 != 0 ) ) 
        {
           complex<double> tmp = c0;
           complex<double> tmp1 = c0;
           complex<double> tmp2 = c0;
           const double t0 = k0 * cell.bmat(0) + k1 * cell.bmat(1);
           const double t1 = k0 * cell.bmat(3) + k1 * cell.bmat(4);
           double gp = sqrt( t0*t0 + t1*t1 );
           for ( int il = 0; il < basis.rod_size(irod); il++ )
           {
              int ig = basis.rod_first(irod) + il;
              int k2 = basis.rod_lmin(irod) + il;  // gz index
              int iz = ( k2 < 0 ) ? k2 + np2v : k2; // z index for vg2
              double kn = double(k2) * tpi/L;
              double cc = cos( kn * z0 );
              double ss = sin( kn * z0 );
              vg2[iz] = (fpi * rhogt[ig]) / (gp*gp + kn*kn);
              // should we be using rhoelg or rhogt here? need to check...
              if ( esm_bc == "bc1" ) 
              {
                 complex<double> ci( 2.0, 3.0 );
                 tmp1 += rhogt[ig] * (cc + ci*ss) / (gp - ci*kn);
                 tmp2 += rhogt[ig] * (cc - ci*ss) / (gp + ci*kn);
              }
              else if ( esm_bc == "bc2" ) 
              {
                 tmp = ((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.0*gp); 
                 tmp1 += rhogt[ig] * (cc + ci*ss) / (gp*gp + kn*kn) * tmp;
                 tmp = ((gp-ci*kn)*exp(gp*(z1-z0))+(gp+ci*kn)*exp(-gp*(z1-z0)))/(2.0*gp); 
                 tmp2 += rhogt[ig] * (cc - ci*ss) / (gp*gp + kn*kn) * tmp;
              }
              else if ( esm_bc == "bc3" ) 
              {
                 tmp = ((gp+ci*kn)*exp(gp*(z1-z0))+(gp-ci*kn)*exp(-gp*(z1-z0)))/(2.0*gp); 
                 tmp1 += rhogt[ig] * (cc + ci*ss) / (gp*gp + kn*kn) * tmp;
                 tmp = (gp - ci*kn) / gp;
                 tmp2 += rhogt[ig] * (cc - ci*ss) / (gp*gp + kn*kn) * tmp;
              }
           }

           vft->backward_1z( vg2, vg2b ); 

           //ewd DEBUG
           //cout << "DEBUG1, mype = " << s_.ctxt_.mype() << ", vg2b[0] = " << vg2b[0] << ", vg2b[1] = " << vg2b[1] << ", vg2b[2] = " << vg2b[2] << endl;

           
           for ( int iz = 0; iz < np2v; iz++ )
           {
              double z;
              if ( iz <= np2half ) 
              {
                 z = double(iz) / double(np2v) * L;
              }
              else
              {
                 z = double(iz - np2v) / double(np2v) * L;
              }
              if ( esm_bc == "bc1" )
              {
                 vg2b[iz] += -tpi/gp*(exp(gp*(z-z0))*tmp1+exp(-gp*(z+z0))*tmp2);
              }
              else if ( esm_bc == "bc2" )
              {
                 vg2b[iz] += -fpi*(exp(gp*(z-z1))-exp(-gp*(z+3.0*z1)))*tmp1 
                     /(1.0-exp(-4.0*gp*z1)) 
                     +fpi*(exp(gp*(z-3.0*z1))-exp(-gp*(z+z1)))*tmp2 
                     /(1.0-exp(-4.0*gp*z1));
              }
              else if ( esm_bc == "bc3" )
              {
                 vg2b[iz]+= -fpi*exp(gp*(z-z1))*tmp1 
                     +tpi*(exp(gp*(z-z0-2.0*z1))-exp(-gp*(z+z0)))*tmp2;
              }
           }
           
           vft->forward_1z( vg2b, vg2 );
           
           //ewd DEBUG
           //cout << "DEBUG2, mype = " << s_.ctxt_.mype() << ", vg2[0] = " << vg2[0] << ", vg2[1] = " << vg2[1] << ", vg2[2] = " << vg2[2] << endl;

           for ( int il = 0; il < basis.rod_size(irod); il++ )
           {
              int ig = basis.rod_first(irod) + il;
              int k2 = basis.rod_lmin(irod) + il;
              int iz = ( k2 < 0 ) ? k2 + np2v : k2;
              vhart_g[ig] = vg2[iz] * 2.0;
           }
        } 
        //
        // else g_parallel = 0:
        //
        else // if ( (k0 == 0) && (k1 == 0) )
        {
           complex<double> rhog0 = c0; // charge density for gamma
           complex<double> tmp1 = c0;
           complex<double> tmp2 = c0;
           complex<double> tmp3 = c0;
           complex<double> tmp4 = c0;
           complex<double> f1 = c0;
           complex<double> f2 = c0;
           complex<double> f3 = c0;
           complex<double> f4 = c0;
           int nz_l = np2half + esm_nfit;
           int nz_r = np2half - esm_nfit;
           double z_l = double(nz_l) * L / double(np2v) - L;
           double z_r = double(nz_r) * L / double(np2v);
           assert( basis.rod_lmin(irod) == 0 );
           
           for ( int il = -basis.rod_size(irod)+1; il < basis.rod_size(irod); il++ )
           {
              complex<double> rhog = c0;
              if ( il < 0 )
              {
                 int ig = basis.rod_first(irod) - il;
                 rhog = conj( rhogt[ig] );
              }
              else
              {
                 int ig = basis.rod_first(irod) + il;
                 rhog = rhogt[ig];
              }
              int k2 = basis.rod_lmin(irod) + il;
              int iz = ( k2 < 0 ) ? k2 + np2v : k2;
              
              if ( k2 == 0 )
              {
                 rhog0 = rhog;
                 if ( esm_bc == "bc1" )
                 {
                    vg2[iz] = -tpi * z0*z0 * rhog;
                 }
                 else if ( esm_bc == "bc2" )
                 {
                    vg2[iz] = tpi * (2.0*z1 - z0) * z0 * rhog;
                 }
                 else if ( esm_bc == "bc3" )
                 {
                    vg2[iz] = tpi * (4.0*z1 - z0) * z0 * rhog;
                 }
              }
              else
              {
                 double kn = double(k2) * tpi/L;
                 double cc = cos( kn * z0 );
                 double ss = sin( kn * z0 );
                 if ( esm_bc == "bc1" )
                 {
                    tmp1 += rhog * ci * (cc + ci*ss) / kn;
                    tmp2 += rhog * ci * (cc - ci*ss) / kn;
                    tmp3 += rhog * cc / (kn*kn);
                 }
                 else if ( esm_bc == "bc2" )
                 {
                    tmp1 += rhog * (cc + ci*ss) / (kn*kn);
                    tmp2 += rhog * (cc - ci*ss) / (kn*kn);
                    tmp3 += rhog * ci * cc / kn;
                    tmp4 += rhog * ss / kn;
                 }
                 else if ( esm_bc == "bc3" )
                 {
                    tmp1 += rhog * (cc + ci*ss) / (kn*kn);
                    tmp2 += rhog * (cc - ci*ss) / kn;
                    tmp3 += rhog * (cc + ci*ss) / kn;
                 }
                 vg2[iz] = fpi * rhog / (kn*kn);
                 // for smoothing
                 double c_r = cos( kn * z_r );
                 double s_r = sin( kn * z_r );
                 double c_l = cos( kn * z_l );
                 double s_l = sin( kn * z_l );
                 f1 += fpi * rhog * (c_r + ci*s_r) / (kn*kn);
                 f2 += fpi * rhog * (c_l + ci*s_l) / (kn*kn);
                 f3 += fpi * ci * rhog * (c_r + ci*s_r) / kn;
                 f4 += fpi * ci * rhog * (c_l + ci*s_l) / kn;
              }
           }

           vft->backward_1z( vg2, vg2b );

           //ewd DEBUG
           //cout << "DEBUG3, mype = " << s_.ctxt_.mype() << ", vg2b[0] = " << vg2b[0] << ", vg2b[1] = " << vg2b[1] << ", vg2b[2] = " << vg2b[2] << endl;

           for ( int iz = 0; iz < np2v; iz++ )
           {
              double z = 0.0;
              if ( iz <= np2half )
              {
                 z = double(iz) / double(np2v) * L;
              }
              else
              {
                 z = double(iz - np2v) / double(np2v) * L;
              }
              if ( esm_bc == "bc1" )
              {
                 vg2b[iz] += (-tpi * z*z * rhog0) - (tpi * (z-z0) * tmp1) 
                     - (tpi * (z+z0) * tmp2) - (fpi * tmp3);
              }
              else if ( esm_bc == "bc2" )
              {
                 vg2b[iz] += (-tpi * z*z * rhog0) - (tpi * (z+z1) * tmp1 / z1) 
                     + (tpi * (z-z1) * tmp2 / z1) 
                     - (fpi * z * (z1-z0) / z1 * tmp3) 
                     + (fpi * (z1-z0) * tmp4);
              }
              else if ( esm_bc == "bc3" )
              {
                 vg2b[iz] += (-tpi * (z*z + 2.0*z*z0) * rhog0) - (fpi * tmp1) 
                     - (fpi * ci * (z-z0) * tmp2) - (fpi * ci * (z1-z0) * tmp3);
              }
           }

           // for smoothing
           if ( esm_bc == "bc1" )
           {
              f1 += (-tpi * z_r*z_r * rhog0) - (tpi * (z_r-z0) * tmp1) 
                  - (tpi * (z_r+z0) * tmp2) - (fpi * tmp3) - (tpi * z0*z0 * rhog0);
              f2 += (-tpi * z_l*z_l * rhog0) - (tpi * (z_l-z0) * tmp1)
                  - (tpi * (z_l+z0) * tmp2) - (fpi * tmp3) - (tpi * z0*z0 * rhog0);
              f3 += (-tpi * tmp1) - (tpi * tmp2) - (fpi * z_r * rhog0);
              f4 += (-tpi * tmp1) - (tpi * tmp2) - (fpi * z_l * rhog0);
           }
           else if ( esm_bc == "bc2" )
           {
              f1 += (-tpi * z_r*z_r * rhog0) - (tpi * (z_r+z1) * tmp1 / z1) 
                  + (tpi * (z_r-z1) * tmp2 / z1) - (fpi * z_r * (z1-z0) / z1 * tmp3) 
                  + (fpi * (z1-z0) * tmp4) + (tpi * (2.0*z1 - z0) * z0 * rhog0);
              f2 += (-tpi * z_l*z_l * rhog0) - (tpi * (z_l+z1) * tmp1 / z1) 
                  + (tpi * (z_l-z1) * tmp2 / z1) - (fpi * z_l * (z1-z0) / z1 * tmp3)
                  + (fpi * (z1-z0) * tmp4) + (tpi * (2.0*z1 - z0) * z0 * rhog0);
              f3 += (-fpi * z_r * rhog0) - (tpi * tmp1 / z1) + (tpi * tmp2 / z1)
                  - (fpi * (z1-z0) / z1 * tmp3);
              f4 += (-fpi * z_l * rhog0) - (tpi * tmp1 / z1) + (tpi * tmp2 / z1)
                  - (fpi * (z1-z0) / z1 * tmp3);
           }
           else if ( esm_bc == "bc3" )
           {
              f1 += (-tpi * (z_r*z_r + 2.0*z_r*z0) * rhog0) - (fpi * tmp1)
                  - (fpi * ci * (z_r-z1) * tmp2) - (fpi * ci * (z1-z0) * tmp3)
                  + (tpi * (4.0*z1 - z0) * z0 * rhog0);
              f2 += (-tpi * (z_l*z_l + 2.0*z_l*z0) * rhog0) - (fpi * tmp1) 
                  - (fpi * ci * (z_l-z1) * tmp2) - (fpi * ci * (z1-z0) * tmp3)
                  + (tpi * (4.0*z1 - z0) * z0 * rhog0);
              f3 += (-tpi * (2.0*z_r + 2.0*z0) * rhog0) - (fpi * ci * tmp2);
              f4 += (-tpi * (2.0*z_l + 2.0*z0) * rhog0) - (fpi * ci * tmp2);
           }
           // for smoothing
           // factor 2 will be multiplied later (at vhart_g <= vg2)
           z_l += L;
           double z_l2 = z_l * z_l;
           double z_r2 = z_r * z_r;
           double z_t3 = z_l - z_r;
           z_t3 = z_t3 * z_t3 * z_t3;
           complex<double> a0 = (f1 * z_l2 * (z_l-3.0*z_r) + z_r * (f3 * z_l2 * (-z_l+z_r)
                                                                    + z_r * (f2 * (3.0*z_l-z_r) + f4 * z_l * (-z_l+z_r)))) / z_t3;
           complex<double> a1 = (f3 * z_l2*z_l + z_l * (6.0*f1 - 6.0*f2 + (f3 + 2.0*f4) * z_l) 
                                 * z_r  - (2.0*f3 + f4) * z_l * z_r2 -f4 * z_r2*z_r)/ z_t3;
           complex<double> a2 = (-3.0 * f1 * (z_l+z_r) + 3.0 * f2 * (z_l+z_r) - (z_l-z_r) 
                                 * (2.0*f3*z_l + f4*z_l + f3*z_r + 2.0*f4*z_r)) / z_t3;
           complex<double> a3 = (2.0*f1 - 2.0*f2 + (f3+f4) * (z_l-z_r)) / z_t3;

           for ( int iz = nz_r; iz <= nz_l; iz++ )
           {
              double z = double(iz) / double(np2v) * L;
              vg2b[iz] = a0 + a1*z + a2*z*z + a3*z*z*z;
           }

           vft->forward_1z( vg2b, vg2 );

           //ewd DEBUG
           //cout << "DEBUG4, mype = " << s_.ctxt_.mype() << ", vg2[0] = " << vg2[0] << ", vg2[1] = " << vg2[1] << ", vg2[2] = " << vg2[2] << endl;

           for ( int il = 0; il < basis.rod_size(irod); il++ )
           {
              int ig = basis.rod_first(irod) + il;
              int k2 = basis.rod_lmin(irod) + il;
              int iz = ( k2 < 0 ) ? k2 + np2v : k2;
              vhart_g[ig] = vg2[iz] * 2.0;
           }
        }
     } // irod loop

     // ESM Hartree energy term

     double eh = 0.0;
     for ( int ig = 0; ig < ngloc; ig++ )
     {
        //
        //  Add Hartree potential into local potential
        //
        vlocal_g[ig] = vion_local_g[ig] + vhart_g[ig];
        //
        //  Add contribution to Hartree energy
        //
        if ( vbasis_->context().mype() == 0 && ig == 0 ) // gamma only correction for G=0
        {
           eh += real( vhart_g[ig] * conj( rhogt[ig] ) ) * 0.5;
        }
        else
        {
           eh += real( vhart_g[ig] * conj( rhogt[ig] ) );
        }
     }

     // factor 1/2 from definition of Ehart cancels with half sum over G
     tsum[1] = eh * omega;
     vbasis_->context().dsum(2,1,&tsum[0],2);
     eps_   = tsum[0];
     ehart_ = tsum[1];

     //ewd DEBUG
     // print vhart for debug
     if (false)
     {
        for ( int irod = 0; irod < basis.nrod_loc(); irod++ )
        {
           int k0 = basis.rod_h(irod);
           int k1 = basis.rod_k(irod);
           if( k0 != 0 || k1 != 0 ) continue;
           
           // g_parallel = 0
           vector<complex<double> > vg2 ( np2v, c0 );
           vector<complex<double> > vg2b( np2v, c0 );
           
           assert( basis.rod_lmin(irod) == 0 );
           
           for ( int il = -basis.rod_size(irod)+1; il < basis.rod_size(irod); il++ )
           {
              complex<double> v = 0.0;
              if ( il < 0 )
              {
                 int ig = basis.rod_first(irod) - il;
                 v = conj( vhart_g[ig] );
              }
              else
              {
                 int ig = basis.rod_first(irod) + il;
                 v = vhart_g[ig];
              }
              int k2 = basis.rod_lmin(irod) + il;
              int iz = ( k2 < 0 ) ? k2 + np2v : k2;
              vg2[iz] = v;
           }
           
           vft->backward_1z( vg2, vg2b );
           
           for( int k2 = np2half+1-np2v; k2 <= np2half; ++k2 ) {
              double z = double(k2) / double(np2v) * L;
              int iz = ( k2 < 0 ) ? k2 + np2v : k2;
              cout << "ESM_debug " << z << " " << real(vg2b[iz]) << endl;
           }
           break;
        } // print vhart for debug
     }

     
  }  // end if ESM

  // v_eff = vxc_g + vion_local_g + vhart_g
  if (s_.ctrl.ultrasoft)
     for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
        for ( int ig = 0; ig < ngloc; ig++ )
           veff_g[ispin][ig] = vxc_g[ispin][ig] + vlocal_g[ig];

        
  // FT to tmpr_r
  vft->backward(&vlocal_g[0],&tmp_r[0]);
 
  // add local potential in tmp_r to v_r[ispin][i]
  // v_r contains the xc potential
  const int size = tmp_r.size();
  if ( wf_.nspin() == 1 ) {
     for ( int i = 0; i < size; i++ ) {
        v_r[0][i] += real(tmp_r[i]);
     }
  }
  else {
     for ( int i = 0; i < size; i++ ) {
        const double vloc = real(tmp_r[i]);
        v_r[0][i] += vloc;
        v_r[1][i] += vloc;
     }
  }
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::update_harris(void) {
   // update terms needed to compute Harris-Foulkes estimate of the energy
   // called when the charge density is calculated, before mixing
   //
   // NOTE: this routine fills v_r, rhoelg with unmixed values.  Call update_vhxc after mixing.
   
   tmap["harris"].start();

   const UnitCell& cell = wf_.cell();
   const double omega = cell.volume();
   const double omega_inv = 1.0 / omega;
   const double *const g2i = vbasis_->g2i_ptr();
   const double fpi = 4.0 * M_PI;
   const int ngloc = vbasis_->localsize();
   double tsum[2];
  tsum[0] = 0.0;
  tsum[1] = 0.0;
   
   // compute total electronic density: rhoelg = rho_up + rho_dn
   if ( wf_.nspin() == 1 ) {
      for ( int ig = 0; ig < ngloc; ig++ ) {
         rhoelg[ig] = omega_inv * cd_.rhog[0][ig];
      }
   }
   else {
      for ( int ig = 0; ig < ngloc; ig++ ) {
         rhoelg[ig] = omega_inv * ( cd_.rhog[0][ig] + cd_.rhog[1][ig] );
      }
   }
  
   // update XC energy and potential
  xcp->update(v_r);
  eharris_ = xcp->exc();

  // compute local potential energy: 
  // integral of el. charge times ionic local pot.
  int len=2*ngloc,inc1=1;
  if (vbasis_->real()) { 
    tsum[0] = 2.0 * ddot(&len,(double*)&rhoelg[0],&inc1,
                         (double*)&vion_local_g[0],&inc1);
    // remove double counting for G=0
    if ( vbasis_->context().myrow() == 0 ) {
      tsum[0] -= real(conj(rhoelg[0])*vion_local_g[0]);
    }
  }
  else {
    tsum[0] = ddot(&len,(double*)&rhoelg[0],&inc1,
                         (double*)&vion_local_g[0],&inc1);
  }
  tsum[0] *= omega;
  
  // Hartree energy
  double ehsum = 0.0;
  for ( int ig = 0; ig < ngloc; ig++ ) {
    const complex<double> tmp = rhoelg[ig] + rhopst[ig];
    ehsum += norm(tmp) * g2i[ig];
  }
  double vfact = vbasis_->real() ? 1.0 : 0.5;
  tsum[1] = vfact * omega * fpi * ehsum;
  
  vbasis_->context().dsum(2,1,&tsum[0],2);


  //ewd DEBUG
  if (false && vbasis_->context().mype() == 0)
     cout << "EF.UPDATE_HARRIS, exc = " << eharris_ << ", ehsum = " << tsum[1] << ", eps = " << tsum[0] << endl;

  eharris_ += tsum[0] + tsum[1];
  
  tmap["harris"].stop();

}

////////////////////////////////////////////////////////////////////////////////
// AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the potential
// untouched and only recalculates the energy terms
void EnergyFunctional::update_exc_ehart_eps(void)
{
  const Wavefunction& wf = s_.wf;
  const UnitCell& cell = wf.cell();
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const double *const g2i = vbasis_->g2i_ptr();
  const double fpi = 4.0 * M_PI;
  const int ngloc = vbasis_->localsize();
  double tsum[2];

  // compute total electronic density: rhoelg = rho_up + rho_dn
  if ( wf.nspin() == 1 )
  {
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhoelg[ig] = omega_inv * cd_.rhog[0][ig];
    }
  }
  else
  {
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhoelg[ig] = omega_inv * ( cd_.rhog[0][ig] + cd_.rhog[1][ig] );
    }
  }

  // update XC energy and potential
  tmap["exc"].start();
  xcp->update_exc(v_r);
  exc_ = xcp->exc();
  tmap["exc"].stop();

  // compute local potential energy:
  // integral of el. charge times ionic local pot.
  // AS: the current charge density has to be used
  int len=2*ngloc,inc1=1;
  tsum[0] = 2.0 * ddot(&len,(double*)&rhoelg[0],&inc1,
         (double*)&vion_local_g[0],&inc1);
  // remove double counting for G=0
  if ( vbasis_->context().myrow() == 0 )
  {
    tsum[0] -= real(conj(rhoelg[0])*vion_local_g[0]);
  }
  tsum[0] *= omega; // tsum[0] contains eps

  // Hartree energy
  // AS: both charge densities have to be used
  ehart_ = 0.0;
  double ehsum = 0.0;
  complex<double> ehsum_cplx = complex<double>(0.0,0.0);
  for ( int ig = 0; ig < ngloc; ig++ )
  {
    const complex<double> tmp = rhoelg[ig] + rhopst[ig];
    const complex<double> hamil_tmp = hamil_rhoelg[ig] + rhopst[ig];
    ehsum_cplx += tmp * conj(hamil_tmp) * g2i[ig];
  }
  ehsum = real(ehsum_cplx);

  // factor 1/2 from definition of Ehart cancels with half sum over G
  // Note: rhogt[ig] includes a factor 1/Omega
  // Factor omega in next line yields prefactor 4 pi / omega in
  tsum[1] = omega * fpi * ehsum;
  // tsum[1] contains ehart

  vbasis_->context().dsum(2,1,&tsum[0],2);
  eps_   = tsum[0];
  ehart_ = tsum[1];
}

////////////////////////////////////////////////////////////////////////////////
double EnergyFunctional::energy(bool compute_hpsi, Wavefunction& dwf,
              bool compute_forces, vector<vector<double> >& fion,
                                bool compute_stress, valarray<double>& sigma) {
  const bool debug_stress = compute_stress && 
    s_.ctrl.debug.find("STRESS") != string::npos;
  const double fpi = 4.0 * M_PI;

  const UnitCell& cell(wf_.cell());
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const int ngloc = vbasis_->localsize();
  const double *const g2i = vbasis_->g2i_ptr();
  const bool use_confinement = s_.ctrl.ecuts > 0.0;
  
  if ( compute_hpsi ) {
    for ( int ispin = 0; ispin < dwf.nspin(); ispin++ ) {
      if (dwf.spinactive(ispin)) {
        for ( int ikp=0; ikp<dwf.nkp(); ikp++) {
          if (dwf.kptactive(ikp)) {
            assert(dwf.sd(ispin,ikp) != 0);
            dwf.sd(ispin,ikp)->c().clear();
          }
        }        
      }
    }
  }

  // kinetic energy
  tmap["ekin"].start();
  
  // compute ekin, confinement energy, stress from ekin and econf
  // ekin = sum_G |c_G|^2  G^2
  // econf = sum_G |c_G|^2 fstress[G]
  // stress_ekin_ij = (1/Omega) sum_G |c_G|^2 * 2 * G_i * G_j
  // stress_econf_ij = (1/Omega) sum_G |c_G|^2 * dfstress[G] * G_i * G_j
  ekin_ = 0.0;
  econf_ = 0.0;
  sigma_ekin = 0.0;
  sigma_econf = 0.0;
  valarray<double> sum(0.0,14), tsum(0.0,14);

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          const double weight = wf_.weight(ikp);
          const SlaterDet& sd = *(wf_.sd(ispin,ikp));
          const Basis& wfbasis = wf_.sd(ispin,ikp)->basis();

          // factor fac in next lines: 2.0 for G and -G (if basis is real) and
          // 0.5 from 1/(2m)
          // note: if basis is real, the factor of 2.0 for G=0 need not be
          // corrected since G^2 = 0
          // Note: the calculation of fac in next line is valid only for nkp=1
          // If k!=0, kpg2(0) !=0 and the ig=0 coefficient must be dealt with 
          // separately
          const double fac = wfbasis.real() ? 1.0 : 0.5;
          const ComplexMatrix& c = sd.c();
          const Context& sdctxt = sd.context();
              
          // compute psi2sum(G) = fac * sum_G occ(n) psi2(n,G)
          const int ngwloc = wfbasis.localsize();
          valarray<double> psi2sum(ngwloc);
          const complex<double>* p = c.cvalptr();
          const int mloc = c.mloc();
          const int nloc = c.nloc();
          const double * const occ = sd.occ_ptr();
          const double *const kpg2  = wfbasis.kpg2_ptr();
          const double *const kpg_x = wfbasis.kpgx_ptr(0);
          const double *const kpg_y = wfbasis.kpgx_ptr(1);
          const double *const kpg_z = wfbasis.kpgx_ptr(2);
          tsum = 0.0;

          //ewd:  this pragma causes incorrect results
          //#pragma omp parallel for
          for ( int ig = 0; ig < ngwloc; ig++ ) {
            double tmpsum = 0.0;
            for ( int lj=0; lj < c.nblocks(); lj++ )
            {
               for ( int jj=0; jj < c.nbs(lj); jj++ )
               {
                  // global state index
                  const int nglobal = c.j(lj,jj);
                  const int norig = lj*c.nb()+jj;
                  const double psi2 = norm(p[ig+norig*mloc]);
                  tmpsum += fac * occ[nglobal] * psi2;
               }
            }
            psi2sum[ig] = tmpsum;
        
            // accumulate contributions to ekin,econf,sigma_ekin,sigma_econf in tsum
            const double psi2s = psi2sum[ig];
            // tsum[0]: ekin partial sum
            tsum[0] += psi2s * kpg2[ig];
              
            if ( compute_stress ) {
              const double tkpgx = kpg_x[ig];
              const double tkpgy = kpg_y[ig];
              const double tkpgz = kpg_z[ig];
 
              const double fac_ekin = 2.0 * psi2s;
                
              tsum[1]  += fac_ekin * tkpgx * tkpgx;
              tsum[2]  += fac_ekin * tkpgy * tkpgy;
              tsum[3]  += fac_ekin * tkpgz * tkpgz;
              tsum[4]  += fac_ekin * tkpgx * tkpgy;
              tsum[5]  += fac_ekin * tkpgy * tkpgz;
              tsum[6]  += fac_ekin * tkpgx * tkpgz;
              
            }
            // tsum[0-6] contains the contributions to
            // ekin, sigma_ekin, from vector ig
          } // ig
          
          if ( use_confinement ) {
            const valarray<double>& fstress = cfp[ispin][ikp]->fstress();
            const valarray<double>& dfstress = cfp[ispin][ikp]->dfstress();
            for ( int ig = 0; ig < ngwloc; ig++ ) {
              const double psi2s = psi2sum[ig];
              // tsum[7]: econf partial sum
              tsum[7] += psi2s * fstress[ig];
              
              if ( compute_stress ) {
                const double tkpgx = kpg_x[ig];
                const double tkpgy = kpg_y[ig];
                const double tkpgz = kpg_z[ig];
              
                const double fac_econf = psi2s * dfstress[ig];
                tsum[8]  += fac_econf * tkpgx * tkpgx;
                tsum[9]  += fac_econf * tkpgy * tkpgy;
                tsum[10] += fac_econf * tkpgz * tkpgz;
                tsum[11] += fac_econf * tkpgx * tkpgy;
                tsum[12] += fac_econf * tkpgy * tkpgz;
                tsum[13] += fac_econf * tkpgx * tkpgz;
              }
              // tsum[7-13] contains the contributions to
              // econf,sigma_econf from vector ig
            } // ig
          }
          sum += weight * tsum;
        }
      }
    }
  }
  assert(wf_.weightsum() > 0.0);
  sum /= wf_.weightsum();

  // sum contains the contributions to ekin, etc.. from this task
  wf_.wfcontext()->dsum(14,1,&sum[0],14);
 
  ekin_  = sum[0];
  econf_ = sum[7];
  if (compute_stress) {
    sigma_ekin[0] = sum[1];
    sigma_ekin[1] = sum[2];
    sigma_ekin[2] = sum[3];
    sigma_ekin[3] = sum[4];
    sigma_ekin[4] = sum[5];
    sigma_ekin[5] = sum[6];
 
    sigma_econf[0] = sum[8];
    sigma_econf[1] = sum[9];
    sigma_econf[2] = sum[10];
    sigma_econf[3] = sum[11];
    sigma_econf[4] = sum[12];
    sigma_econf[5] = sum[13];
  
    sigma_ekin *= omega_inv;
    sigma_econf *= omega_inv;

    // symmetrize sigma_ekin and sigma_econf
    if ( s_.symmetries.nsym() > 0) {
      const int nsym_ = s_.symmetries.nsym();
      valarray<double> sigma_ekin_xtal(6);
      valarray<double> sigma_ekin_xtal_sym(6);
      // convert sigma to crystal coordinates to match symmetry operators
      wf_.cell().cart_to_crystal(&sigma_ekin[0],&sigma_ekin_xtal[0]);

      // symmetrize sigma_ekin_xtal
      valarray<double> sigma_sum(6);
      for (int i=0; i<6; i++)
        sigma_sum[i] = sigma_ekin_xtal[i];  // identity operation
      for ( int isym = 0; isym < nsym_; isym++) {
        s_.symmetries.symlist[isym]->applyToTensor(&sigma_ekin_xtal[0],&sigma_ekin_xtal_sym[0]);
        for (int i=0; i<6; i++)
          sigma_sum[i] += sigma_ekin_xtal_sym[i];
      }
      double nsyminv = 1./((double)nsym_+1.);
      for (int i=0; i<6; i++)
        sigma_sum[i] *= nsyminv;
    
      // convert sigma back to Cartesian coordinates
      wf_.cell().crystal_to_cart(&sigma_sum[0],&sigma_ekin[0]);

      valarray<double> sigma_econf_xtal(6);
      valarray<double> sigma_econf_xtal_sym(6);
      wf_.cell().cart_to_crystal(&sigma_econf[0],&sigma_econf_xtal[0]);
      // symmetrize sigma_econf_xtal
      for (int i=0; i<6; i++)
        sigma_sum[i] = sigma_econf_xtal[i];  // identity operation
      for ( int isym = 0; isym < nsym_; isym++) {
        s_.symmetries.symlist[isym]->applyToTensor(&sigma_econf_xtal[0],&sigma_econf_xtal_sym[0]);
        for (int i=0; i<6; i++)
          sigma_sum[i] += sigma_econf_xtal_sym[i];
      }
      for (int i=0; i<6; i++)
        sigma_sum[i] *= nsyminv;
    
      // convert sigma back to Cartesian coordinates
      wf_.cell().crystal_to_cart(&sigma_sum[0],&sigma_econf[0]);
    }
  }  

  tmap["ekin"].stop();
  
  // Stress from Eps
  sigma_eps = 0.0;
  if ( compute_stress ) {
    tsum = 0.0;
    const double *const gi  = vbasis_->gi_ptr();
    const double *const g_x = vbasis_->gx_ptr(0);
    const double *const g_y = vbasis_->gx_ptr(1);
    const double *const g_z = vbasis_->gx_ptr(2);
    for ( int ig = 0; ig < ngloc; ig++ ) {
      // factor of 2 in next line: G and -G
      // note: gi[0] == 0.0 in next line (no division by zero)
      const double fac = 2.0 * gi[ig] * 
        real( conj(rhoelg[ig]) * dvion_local_g[ig] );
      
      const double tgx = g_x[ig];
      const double tgy = g_y[ig];
      const double tgz = g_z[ig];
      
      tsum[0] += fac * tgx * tgx;
      tsum[1] += fac * tgy * tgy;
      tsum[2] += fac * tgz * tgz;
      tsum[3] += fac * tgx * tgy;
      tsum[4] += fac * tgy * tgz;
      tsum[5] += fac * tgx * tgz;
    }
    vbasis_->context().dsum(6,1,&tsum[0],6);

    sigma_eps[0] = eps_ * omega_inv + tsum[0];
    sigma_eps[1] = eps_ * omega_inv + tsum[1];
    sigma_eps[2] = eps_ * omega_inv + tsum[2];
    sigma_eps[3] = tsum[3];
    sigma_eps[4] = tsum[4];
    sigma_eps[5] = tsum[5];
  }
  
  // Stress from Hartree energy
  if ( compute_stress ) {
    tsum = 0.0;
    const double *const g_x = vbasis_->gx_ptr(0);
    const double *const g_y = vbasis_->gx_ptr(1);
    const double *const g_z = vbasis_->gx_ptr(2);
  
    for ( int ig = 0; ig < ngloc; ig++ ) {
      const double temp = norm(rhogt[ig]) * g2i[ig] * g2i[ig];
      const double tgx = g_x[ig];
      const double tgy = g_y[ig];
      const double tgz = g_z[ig];

      tsum[0] += temp * tgx * tgx;
      tsum[1] += temp * tgy * tgy;
      tsum[2] += temp * tgz * tgz;
      tsum[3] += temp * tgx * tgy;
      tsum[4] += temp * tgy * tgz;
      tsum[5] += temp * tgx * tgz;
    }

    for ( int is = 0; is < nsp_; is++ ) {
      // Factor 2 in next line: 2*real(x) = ( x + h.c. )
      // Factor 0.5 in next line: definition of sigma
      const double fac2 = 2.0 * 0.25 * rcps_[is]*rcps_[is];
      complex<double> *s = &sf.sfac[is][0];
      for ( int ig = 0; ig < ngloc; ig++ ) {
        const complex<double> sg = s[ig];
        const complex<double> rg = rhogt[ig];
        // next line: keep only real part
        const double temp = fac2 *
        ( rg.real() * sg.real() +
          rg.imag() * sg.imag() )
        * rhops[is][ig] * g2i[ig];
        
        const double tgx = g_x[ig];
        const double tgy = g_y[ig];
        const double tgz = g_z[ig];

        tsum[0] += temp * tgx * tgx;
        tsum[1] += temp * tgy * tgy;
        tsum[2] += temp * tgz * tgz;
        tsum[3] += temp * tgx * tgy;
        tsum[4] += temp * tgy * tgz;
        tsum[5] += temp * tgx * tgz;
      }
    }
    vbasis_->context().dsum(6,1,&tsum[0],6);
    // Factor in next line:
    //  factor 2 from G and -G
    //  factor fpi from definition of sigma
    //  no factor 1/Omega^2 (already included in rhogt[ig] above)
    sigma_ehart[0] = ehart_ * omega_inv - 2.0 * fpi * tsum[0];
    sigma_ehart[1] = ehart_ * omega_inv - 2.0 * fpi * tsum[1];
    sigma_ehart[2] = ehart_ * omega_inv - 2.0 * fpi * tsum[2];
    sigma_ehart[3] = - 2.0 * fpi * tsum[3];
    sigma_ehart[4] = - 2.0 * fpi * tsum[4];
    sigma_ehart[5] = - 2.0 * fpi * tsum[5];
  } // compute_stress
  
  // Stress from exchange-correlation
  if ( compute_stress ) {
    xcp->compute_stress(sigma_exc);
  }
  
  // zero ionic forces
  if (compute_forces) {
    for ( int is = 0; is < fion.size(); is++ )
      for ( int ia = 0; ia < fion[is].size(); ia++ )
        fion[is][ia] = 0.0;
  }

  vector<vector<double> > fion_nl;
  if (compute_forces) {
    fion_nl.resize(nsp_);
    for (int is=0; is<nsp_; is++) 
      fion_nl[is].resize(3*na_[is]);
  }
  
  // Non local energy
  tmap["nonlocal"].start();
  // calculate nonlocal energy, averaged over all kpoints on local sdcontext
  double enlsum[2] = {0.0, 0.0};

  // temporary array for averaging non-local stress over k-points
  valarray<double> tsigma_enl;
  tsigma_enl.resize(6);
  if (compute_stress)
    for (int i=0; i<6; i++) 
      sigma_enl[i] = 0.0;

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for (int kloc=0; kloc<wf_.nkptloc(); kloc++) {

        //ewd DEBUG:  calculate v_r*rhor to compare against PWSCF deband
        //double deband = 0.0;
        //for (int i=0; i<vft->np012loc(); i++)
        // deband += v_r[ispin][i]*cd_.rhor[ispin][i];
        //double fac = omega/(double)vft->np012loc();
        //deband *= fac;
        //cout << "EF.DEBAND deband = " << deband << "    , fac = " << fac << endl;

        double wt = dwf.weight(dwf.kptloc(kloc));
        enlsum[0] += wt*nlp[ispin][kloc]->energy(compute_hpsi,*dwf.sdloc(ispin,kloc),
                                              compute_forces, fion_nl, compute_stress,
                                              tsigma_enl,veff_g[ispin]);
        enlsum[1] += wt;
        if (compute_forces) {
          for (int is=0; is<nsp_; is++) {
            for (int ia=0; ia<na_[is]; ia++) {
              fion[is][3*ia+0] += wt*fion_nl[is][3*ia];
              fion[is][3*ia+1] += wt*fion_nl[is][3*ia+1];
              fion[is][3*ia+2] += wt*fion_nl[is][3*ia+2];
            }
          }
        }
        if (compute_stress)
          for (int i=0; i<6; i++) 
            sigma_enl[i] += wt*tsigma_enl[i];
      }
    }
  }
  dwf.wfcontext()->dsum('r',2,1,&enlsum[0],2);     // weighted average over all kpoints
  assert(enlsum[1] != 0.0);
  if (wf_.nspin() == 2)
    enlsum[1] *= 0.5;
  enl_ = enlsum[0]/enlsum[1];

  if (compute_forces) {
    for (int is=0; is<nsp_; is++) 
      dwf.wfcontext()->dsum('r',3*na_[is],1,&fion[is][0],3*na_[is]);
    double wtsuminv = 1.0/enlsum[1];
    for (int is=0; is<nsp_; is++) {
      for (int ia=0; ia<na_[is]; ia++) {
        fion[is][3*ia+0] *= wtsuminv;
        fion[is][3*ia+1] *= wtsuminv;
        fion[is][3*ia+2] *= wtsuminv;
      }
    }
  }
  if (compute_stress) {
    dwf.wfcontext()->dsum('r',6,1,&sigma_enl[0],6);
    double wtsuminv = 1.0/enlsum[1];
    for (int i=0; i<6; i++) 
      sigma_enl[i] *= wtsuminv;

    // symmetrize sigma_enl
    if ( s_.symmetries.nsym() > 0) {
      const int nsym_ = s_.symmetries.nsym();

      valarray<double> sigma_enl_xtal(6);
      valarray<double> sigma_enl_xtal_sym(6);
      // convert sigma to crystal coordinates to match symmetry operators
      wf_.cell().cart_to_crystal(&sigma_enl[0],&sigma_enl_xtal[0]);

      // symmetrize sigma_enl_xtal
      valarray<double> sigma_sum(6);
      for (int i=0; i<6; i++)
        sigma_sum[i] = sigma_enl_xtal[i];  // identity operation
      for ( int isym = 0; isym < nsym_; isym++) {
        s_.symmetries.symlist[isym]->applyToTensor(&sigma_enl_xtal[0],&sigma_enl_xtal_sym[0]);
        for (int i=0; i<6; i++)
          sigma_sum[i] += sigma_enl_xtal_sym[i];
      }
      double nsyminv = 1./((double)nsym_+1.);
      for (int i=0; i<6; i++)
        sigma_sum[i] *= nsyminv;

      // convert sigma back to Cartesian coordinates
      wf_.cell().crystal_to_cart(&sigma_sum[0],&sigma_enl[0]);
    }
  }
  tmap["nonlocal"].stop();

  // DFT+U contribution
  ehub_ = 0.0;
  if (s_.ctrl.dft_plus_u) {
    tmap["hubbard"].start();
    
    // calculate Hubbard contribution to energy, averaged over all kpoints on local sdcontext
    ehub_ = hubp_->energy(compute_hpsi,dwf);

    /*
    double ehubsum[2] = {0.0, 0.0};

    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
      if (wf_.spinactive(ispin)) {
        for (int kloc=0; kloc<wf_.nkptloc(); kloc++) {
          double wt = dwf.weight(dwf.kptloc(kloc));
          ehubsum[0] += wt*hubp[ispin][kloc]->energy(compute_hpsi,*dwf.sdloc(ispin,kloc),ispin,dwf.kpoint(dwf.kptloc(kloc)));
          ehubsum[1] += wt;
        }
      }
    }
    dwf.wfcontext()->dsum('r',2,1,&ehubsum[0],2);     // weighted average over all kpoints
    assert(ehubsum[1] != 0.0);
    if (wf_.nspin() == 2)
      ehubsum[1] *= 0.5;
    ehub_ = ehubsum[0]/ehubsum[1];
    */
    tmap["hubbard"].stop();
  }
  
  ecoul_ = ehart_ + esr_ - eself_;

  ets_ = 0.0;
  if ( s_.ctrl.smearing_width > 0.0 ) {
    const double wf_entropy = wf_.entropy();
    // used when smearing was in Kelvin
    //const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    const double boltz = 0.5; // convert from Ry to Ha
    ets_ = - wf_entropy * s_.ctrl.smearing_width * boltz;
  }
  etotal_ = ekin_ + econf_ + eps_ + enl_ + ecoul_ + exc_ + ets_ + epv_ + ehub_;

  //ewd DEBUG
  //print(cout);
  
  //ewd DEBUG
  /*
  cout << "DEBUG:  ehart_ = " << ehart_ << endl;
  cout << "DEBUG:  PWSCF ewald_ehart = " << pwscf_ewald_ << endl;
  pwscf_ewald_ += esr_ - eself_;
  cout << "DEBUG:  PWSCF ewald_tot = " << pwscf_ewald_ << endl;
  */

  // AS: write the different terms to the total energy to the output
  // AS: helpful especially for testing
  if ( (s_.ctxt_.oncoutpe()) && (s_.ctrl.non_selfc_energy) )
  {
     cout << " AS: EKIN_ " << ekin_ << endl;
     cout << " AS: ECONF_ " << econf_ << endl;
     cout << " AS: EPS_ " << eps_ << endl;
     cout << " AS: ENL_ " << enl_ << endl;
     cout << " AS: EHART_ " << ehart_ << endl;
     cout << " AS: ESR_ " << esr_ << endl;
     cout << " AS: ESELF_ " << eself_ << endl;
     cout << " AS: EXC_ " << exc_ << endl;
     cout << " AS: ETS_ " << ets_ << endl;
     cout << " AS: EEXF_ " << eexf_ << endl;
  }
  
  if ( compute_hpsi ) {
    tmap["hpsi"].start();
    //assert(wf_.nspin()==1);
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
      if (wf_.spinactive(ispin)) {
        for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
          if (wf_.kptactive(ikp)) {
            assert(wf_.sd(ispin,ikp) != 0);
            const SlaterDet& sd = *(wf_.sd(ispin,ikp));
            SlaterDet& sdp = *(dwf.sd(ispin,ikp));
            const ComplexMatrix& c = sd.c();
            const Basis& wfbasis = sd.basis();

            ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
            const int mloc = cp.mloc();
            const double* kpg2 = wfbasis.kpg2_ptr();
            const int ngwloc = wfbasis.localsize();
          
            // Laplacian
            if ( use_confinement ) {
              for ( int n = 0; n < sd.nstloc(); n++ ) {
                assert(cfp[ispin][ikp]!=0); // cfp must be non-zero if this ikp active
                const valarray<double>& fstress = cfp[ispin][ikp]->fstress();
                for ( int ig = 0; ig < ngwloc; ig++ ) {
                  cp[ig+mloc*n] += 0.5 * ( kpg2[ig] + fstress[ig] ) * 
                      c[ig+mloc*n];
                }
              }
            }
            else {
               //ewd:  OMP HERE?
               // AS: e_n are currently only hard-coded and renormalization is disabled by default
               // double as_renorm[4] = {-28.70987,-28.70987,-28.70987,-3.13510};
               //#pragma omp parallel for
               for ( int n = 0; n < sd.nstloc(); n++ ) {
                for ( int ig = 0; ig < ngwloc; ig++ ) {
                  cp[ig+mloc*n] += 0.5 * kpg2[ig] * c[ig+mloc*n];

                  // AS: energy_renorm allows applying a "shifted Hamiltonian" (H-e_n), e_n being the shift
                  // AS: for each individual state n, to the wave function during the time propgation
                  // AS: by subtracting the phase: cp -= 1.0 * as_renorm[n] / 27.2116 * c[ig+mloc*n];
                  // AS: this should prevent the time propagation from blowing up
                  // if (energy_renorm) {
                  // AS: renormalize each state individually
                  // cp[ig+mloc*n] -= 1.0 * as_renorm[n] / 27.2116 * c[ig+mloc*n];
                  // AS: renormalize all states with the same value
                  // cp[ig+mloc*n] -= 1.0 * (-0.8204477) * c[ig+mloc*n];   

                }
              }
            }
            sd.rs_mul_add(*ft[ispin][ikp], &v_r[ispin][0], sdp);
          }
        }
      }
    }
    tmap["hpsi"].stop();
  } // if compute_hpsi

  if ( compute_forces ) {
    const int* idx = vbasis_->idx_ptr();
    const double* gx0 = vbasis_->gx_ptr(0);
    const double* gx1 = vbasis_->gx_ptr(1);
    const double* gx2 = vbasis_->gx_ptr(2);
    vector<complex<double> > trhog;
    for ( int is = 0; is < nsp_; is++ ) {
       for ( int ig = 0; ig < ngloc; ig++ ) {
          double tmp = fpi * rhops[is][ig] * g2i[ig];
          vtemp[ig] =  tmp * conj(rhogt[ig]) + vps[is][ig] * conj(rhoelg[ig]);
       }
       Species *s = s_.atoms.species_list[is];
       if (s->nlcc())
       {
          cd_.nlcc_forceden(is,trhog);
          const double spinfac = 1./(double)wf_.nspin();
          for ( int ig = 0; ig < ngloc; ig++ )
             for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
                vtemp[ig] += spinfac*omega_inv*trhog[ig]*conj(vxc_g[ispin][ig]);
       }
       
       memset((void*)&ftmp[0],0,3*namax_*sizeof(double));
       // loop over atoms of species is
       for ( int ia = 0; ia < na_[is]; ia++ ) {
          double sum0=0.0,sum1=0.0,sum2=0.0;
          double *c0 = sf.cos0_ptr(is,ia);
          double *c1 = sf.cos1_ptr(is,ia);
          double *c2 = sf.cos2_ptr(is,ia);
          double *s0 = sf.sin0_ptr(is,ia);
          double *s1 = sf.sin1_ptr(is,ia);
          double *s2 = sf.sin2_ptr(is,ia);
          //#pragma omp parallel for
          for ( int ig = 0; ig < ngloc; ig++ ) {
             // compute Exp[igr] in 3D as a product of 1D contributions
             // complex<double> teigr = ei0[kv[3*ig+0]] *
             //                         ei1[kv[3*ig+1]] *
             //                         ei2[kv[3*ig+2]];
             const int iii = ig+ig+ig;
             const int kx = idx[iii];
             const int ky = idx[iii+1];
             const int kz = idx[iii+2];
          
             const double cos_a = c0[kx];
             const double cos_b = c1[ky];
             const double cos_c = c2[kz];
             
             const double sin_a = s0[kx];
             const double sin_b = s1[ky];
             const double sin_c = s2[kz];
 
             // Next line: exp(-i*gr) =
             // (cos_a - I sin_a)*(cos_b - I sin_b)*(cos_c - I sin_c)
             double teigr_re = 
                 cos_a*cos_b*cos_c - sin_a*sin_b*cos_c -
                 sin_a*cos_b*sin_c - cos_a*sin_b*sin_c;
             double teigr_im = 
                 sin_a*sin_b*sin_c - sin_a*cos_b*cos_c -
                 cos_a*sin_b*cos_c - cos_a*cos_b*sin_c;
             
             /* fion is real */
             double tmp = teigr_re * vtemp[ig].imag() + 
                 teigr_im * vtemp[ig].real();

             sum0 += tmp * gx0[ig];
             sum1 += tmp * gx1[ig];
             sum2 += tmp * gx2[ig];
             
          }
          ftmp[3*ia]   = sum0;
          ftmp[3*ia+1] = sum1;
          ftmp[3*ia+2] = sum2;
       }
       int len = 3*na_[is];
       vbasis_->context().dsum(len,1,&ftmp[0],len);

       double vfac = vbasis_->real() ? 2.0*omega : omega;
       for ( int i = 0; i < 3*na_[is]; i++ )
          fion[is][i] += fion_esr[is][i] - vfac * ftmp[i];
    }
  }
  
  if (compute_stress)
    sigma = sigma_ekin + sigma_econf + sigma_eps + sigma_enl +
          sigma_ehart + sigma_exc + sigma_esr;
  if ( debug_stress && s_.ctxt_.oncoutpe() ) {
  //if ( compute_stress && s_.ctxt_.oncoutpe() ) {
    // ewd:  add more significant figures to conversion
    const double gpa = 29421.0120;
    //const double gpa = 29421.5;
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout << setprecision(8);
    cout << " <stress_tensor unit=\"GPa\">\n"
         << "   <debugsigma_ekin_xx> " << setw(12) << gpa*sigma_ekin[0] << " </debugsigma_ekin_xx>\n"
         << "   <debugsigma_ekin_yy> " << setw(12) << gpa*sigma_ekin[1] << " </debugsigma_ekin_yy>\n"
         << "   <debugsigma_ekin_zz> " << setw(12) << gpa*sigma_ekin[2] << " </debugsigma_ekin_zz>\n"
         << "   <debugsigma_ekin_xy> " << setw(12) << gpa*sigma_ekin[3] << " </debugsigma_ekin_xy>\n"
         << "   <debugsigma_ekin_yz> " << setw(12) << gpa*sigma_ekin[4] << " </debugsigma_ekin_yz>\n"
         << "   <debugsigma_ekin_xz> " << setw(12) << gpa*sigma_ekin[5] << " </debugsigma_ekin_xz>\n"
         << endl
         << "   <debugsigma_econf_xx> " << setw(12) << gpa*sigma_econf[0] << " </debugsigma_econf_xx>\n"
         << "   <debugsigma_econf_yy> " << setw(12) << gpa*sigma_econf[1] << " </debugsigma_econf_yy>\n"
         << "   <debugsigma_econf_zz> " << setw(12) << gpa*sigma_econf[2] << " </debugsigma_econf_zz>\n"
         << "   <debugsigma_econf_xy> " << setw(12) << gpa*sigma_econf[3] << " </debugsigma_econf_xy>\n"
         << "   <debugsigma_econf_yz> " << setw(12) << gpa*sigma_econf[4] << " </debugsigma_econf_yz>\n"
         << "   <debugsigma_econf_xz> " << setw(12) << gpa*sigma_econf[5] << " </debugsigma_econf_xz>\n"
         << endl
         << "   <debugsigma_eps_xx> " << setw(12) << gpa*sigma_eps[0] << " </debugsigma_eps_xx>\n"
         << "   <debugsigma_eps_yy> " << setw(12) << gpa*sigma_eps[1] << " </debugsigma_eps_yy>\n"
         << "   <debugsigma_eps_zz> " << setw(12) << gpa*sigma_eps[2] << " </debugsigma_eps_zz>\n"
         << "   <debugsigma_eps_xy> " << setw(12) << gpa*sigma_eps[3] << " </debugsigma_eps_xy>\n"
         << "   <debugsigma_eps_yz> " << setw(12) << gpa*sigma_eps[4] << " </debugsigma_eps_yz>\n"
         << "   <debugsigma_eps_xz> " << setw(12) << gpa*sigma_eps[5] << " </debugsigma_eps_xz>\n"
         << endl
         << "   <debugsigma_enl_xx> " << setw(12) << gpa*sigma_enl[0] << " </debugsigma_enl_xx>\n"
         << "   <debugsigma_enl_yy> " << setw(12) << gpa*sigma_enl[1] << " </debugsigma_enl_yy>\n"
         << "   <debugsigma_enl_zz> " << setw(12) << gpa*sigma_enl[2] << " </debugsigma_enl_zz>\n"
         << "   <debugsigma_enl_xy> " << setw(12) << gpa*sigma_enl[3] << " </debugsigma_enl_xy>\n"
         << "   <debugsigma_enl_yz> " << setw(12) << gpa*sigma_enl[4] << " </debugsigma_enl_yz>\n"
         << "   <debugsigma_enl_xz> " << setw(12) << gpa*sigma_enl[5] << " </debugsigma_enl_xz>\n"
         << endl
         << "   <debugsigma_ehart_xx> " << setw(12) << gpa*sigma_ehart[0] << " </debugsigma_ehart_xx>\n"
         << "   <debugsigma_ehart_yy> " << setw(12) << gpa*sigma_ehart[1] << " </debugsigma_ehart_yy>\n"
         << "   <debugsigma_ehart_zz> " << setw(12) << gpa*sigma_ehart[2] << " </debugsigma_ehart_zz>\n"
         << "   <debugsigma_ehart_xy> " << setw(12) << gpa*sigma_ehart[3] << " </debugsigma_ehart_xy>\n"
         << "   <debugsigma_ehart_yz> " << setw(12) << gpa*sigma_ehart[4] << " </debugsigma_ehart_yz>\n"
         << "   <debugsigma_ehart_xz> " << setw(12) << gpa*sigma_ehart[5] << " </debugsigma_ehart_xz>\n"
         << endl
         << "   <debugsigma_exc_xx> " << setw(12) << gpa*sigma_exc[0] << " </debugsigma_exc_xx>\n"
         << "   <debugsigma_exc_yy> " << setw(12) << gpa*sigma_exc[1] << " </debugsigma_exc_yy>\n"
         << "   <debugsigma_exc_zz> " << setw(12) << gpa*sigma_exc[2] << " </debugsigma_exc_zz>\n"
         << "   <debugsigma_exc_xy> " << setw(12) << gpa*sigma_exc[3] << " </debugsigma_exc_xy>\n"
         << "   <debugsigma_exc_yz> " << setw(12) << gpa*sigma_exc[4] << " </debugsigma_exc_yz>\n"
         << "   <debugsigma_exc_xz> " << setw(12) << gpa*sigma_exc[5] << " </debugsigma_exc_xz>\n"
         << endl
         << "   <debugsigma_esr_xx> " << setw(12) << gpa*sigma_esr[0] << " </debugsigma_esr_xx>\n"
         << "   <debugsigma_esr_yy> " << setw(12) << gpa*sigma_esr[1] << " </debugsigma_esr_yy>\n"
         << "   <debugsigma_esr_zz> " << setw(12) << gpa*sigma_esr[2] << " </debugsigma_esr_zz>\n"
         << "   <debugsigma_esr_xy> " << setw(12) << gpa*sigma_esr[3] << " </debugsigma_esr_xy>\n"
         << "   <debugsigma_esr_yz> " << setw(12) << gpa*sigma_esr[4] << " </debugsigma_esr_yz>\n"
         << "   <debugsigma_esr_xz> " << setw(12) << gpa*sigma_esr[5] << " </debugsigma_esr_xz>\n"
         << endl
         << "   <debugsigma_eks_xx> " << setw(12) << gpa*sigma[0] << " </debugsigma_eks_xx>\n"
         << "   <debugsigma_eks_yy> " << setw(12) << gpa*sigma[1] << " </debugsigma_eks_yy>\n"
         << "   <debugsigma_eks_zz> " << setw(12) << gpa*sigma[2] << " </debugsigma_eks_zz>\n"
         << "   <debugsigma_eks_xy> " << setw(12) << gpa*sigma[3] << " </debugsigma_eks_xy>\n"
         << "   <debugsigma_eks_yz> " << setw(12) << gpa*sigma[4] << " </debugsigma_eks_yz>\n"
         << "   <debugsigma_eks_xz> " << setw(12) << gpa*sigma[5] << " </debugsigma_eks_xz>\n"
         << " </stress_tensor>" << endl;
  }
  
  return etotal_;
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::atoms_moved(void)
{
  const AtomSet& atoms = s_.atoms;
  int ngloc = vbasis_->localsize();

  // fill tau0, taum with values in atom_list
  atoms.get_positions(tau0,true);
  sf.update(tau0,*vbasis_);
  
  // compute Fourier coefficients of the local potential
  memset( (void*)&vion_local_g[0], 0, 2*ngloc*sizeof(double) );
  memset( (void*)&dvion_local_g[0], 0, 2*ngloc*sizeof(double) );
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );

  for ( int is = 0; is < atoms.nsp(); is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const complex<double> sg = s[ig];
      rhopst[ig] += sg * rhops[is][ig];
      vion_local_g[ig] += sg * vps[is][ig];
      dvion_local_g[ig] += sg * dvps[is][ig];
    }
  }
  
  // compute esr: pseudocharge repulsion energy
  const UnitCell& cell = wf_.cell();
  const double omega_inv = 1.0 / cell.volume();
  
  esr_  = 0.0;
  sigma_esr = 0.0;
  for ( int is = 0; is < nsp_; is++ )
    for ( int i = 0; i < fion_esr[is].size(); i++ )
      fion_esr[is][i] = 0.0;

  for ( int k = 0; k < nsp_; k++ ) {
    for ( int j = k; j < nsp_; j++ ) {
      double rckj = sqrt ( rcps_[k]*rcps_[k]+rcps_[j]*rcps_[j] );
      int lax = na_[k];
      if ( k == j ) lax--;

      // choose number of images in real-space Ewald sum such that all points
      // at least (e.g.) five sigmas away are included
      // (2*newald+1)*cell.min_wsdist() >= 5*sigma
      const int nsigmas = 5;
      const int newald = (int)( (0.5*nsigmas*rckj/cell.min_wsdist()) + 0.5);
      if ( s_.ctxt_.oncoutpe() ) 
        cout << "<!-- EnergyFunctional:  number of images in real-space ewald sum = " << newald << " for rckj = " << rckj << " -->" << endl;
      
      for ( int l = 0; l < lax; l++ ) {
        int inf = 0;
        if ( k == j ) inf = l+1;
        for ( int m = inf; m < na_[j]; m++ ) {

          for (int ucx=-newald; ucx<=newald; ucx++) {
            for (int ucy=-newald; ucy<=newald; ucy++) {
              for (int ucz=-newald; ucz<=newald; ucz++) {
                
                double xlm = tau0[k][3*l+0] - tau0[j][3*m+0];
                double ylm = tau0[k][3*l+1] - tau0[j][3*m+1];
                double zlm = tau0[k][3*l+2] - tau0[j][3*m+2];
                D3vector vlm(xlm,ylm,zlm);
                cell.fold_in_ws(vlm);

                vlm = vlm + ucx*cell.a(0) + ucy*cell.a(1) + ucz*cell.a(2);
                xlm = vlm.x;
                ylm = vlm.y;
                zlm = vlm.z;
                double rlm = sqrt(xlm*xlm + ylm*ylm + zlm*zlm);
                double arg = rlm / rckj;
                double esr_lm = zv_[k] * zv_[j] * erfc(arg) / rlm;
                esr_ += esr_lm;

                double desr_erfc = 2.0 * zv_[k]*zv_[j]*exp(-arg*arg)/rckj/sqrt(M_PI);
                
                // desrdr = (1/r) d Esr / dr
                double desrdr = -(esr_lm+desr_erfc) / ( rlm*rlm );
                fion_esr[k][3*l+0] -= desrdr * xlm;
                fion_esr[j][3*m+0] += desrdr * xlm;
                fion_esr[k][3*l+1] -= desrdr * ylm;
                fion_esr[j][3*m+1] += desrdr * ylm;
                fion_esr[k][3*l+2] -= desrdr * zlm;
                fion_esr[j][3*m+2] += desrdr * zlm;

                sigma_esr[0] += desrdr * xlm * xlm;
                sigma_esr[1] += desrdr * ylm * ylm;
                sigma_esr[2] += desrdr * zlm * zlm;
                sigma_esr[3] += desrdr * xlm * ylm;
                sigma_esr[4] += desrdr * ylm * zlm;
                sigma_esr[5] += desrdr * xlm * zlm;

              }
            }
          }

        }
      }
    }
  }
  sigma_esr *= - omega_inv;

  // update ultrasoft potentials
  if (s_.ctrl.ultrasoft)
    for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) 
      if (wf_.spinactive(ispin)) 
        for (int k=0; k<nlp[ispin].size(); k++)
          nlp[ispin][k]->update_usfns(vbasis_);

}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::cell_moved(const bool compute_stress) {
  const UnitCell& cell = wf_.cell();
  // resize vbasis_
  vbasis_->resize(cell,wf_.refcell(),vbasis_->ecut());
  
  const int ngloc = vbasis_->localsize();
  const double omega = cell.volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;  
  
  const AtomSet& atoms = s_.atoms;

  for ( int is = 0; is < nsp_; is++ ) {
    Species *s = atoms.species_list[is];
    const double * const g = vbasis_->g_ptr();
    double v,dv;  
    for ( int ig = 0; ig < ngloc; ig++ ) {
      rhops[is][ig] = s->rhopsg(g[ig]) * omega_inv;
      s->dvlocg(g[ig],v,dv);
      vps[is][ig] =  v * omega_inv;
      dvps[is][ig] =  dv * omega_inv;
    }    
  }
  
  // Update confinement potentials
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) 
    if (wf_.spinactive(ispin)) 
      for ( int kp = 0; kp < wf_.nkp(); kp++ ) 
        if ( cfp[ispin][kp] != 0 ) 
          cfp[ispin][kp]->update();
  
  // update non-local potential
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) 
    if (wf_.spinactive(ispin)) 
      for (int k=0; k<nlp[ispin].size(); k++) {
        nlp[ispin][k]->update_twnl(compute_stress);
        if (s_.ctrl.ultrasoft)
          nlp[ispin][k]->update_usfns(vbasis_);
      }
  
  // update Hubbard potential
  if (s_.ctrl.dft_plus_u)
    hubp_->update_phiylm();
}

////////////////////////////////////////////////////////////////////////////////
double EnergyFunctional::casino_ewald(void)
{
   //ewd:  calculate Ewald energy term
   //ewd:  for debugging and CASINO output
   const int ngloc = vbasis_->localsize();
   const double fpi = 4.0 * M_PI;
   const double vfact = vbasis_->real() ? 1.0 : 0.5;
   const double omega = wf_.cell().volume();
   const double *const g2i = vbasis_->g2i_ptr();
   double tsum;
   double ewald = 0.0;
   for ( int ig = 0; ig < ngloc; ig++ )
      ewald += norm(rhopst[ig]) * g2i[ig];
   tsum = vfact * omega * fpi * ewald;
   vbasis_->context().dsum(1,1,&tsum,1);
   ewald = tsum + esr_ - eself_;


   //ewd DEBUG
   if ( false && s_.ctxt_.mype()==0 )
      cout << "EF.EWALD:  tsum = " << tsum << ", esr = " << esr_ << ", eself = " << eself_ << endl;
   
   return ewald;
}  
////////////////////////////////////////////////////////////////////////////////
double EnergyFunctional::casino_vloc(void)
{
   //ewd:  calculate Ewald energy term
   //ewd:  for debugging and CASINO output
   const int ngloc = vbasis_->localsize();
   const double fpi = 4.0 * M_PI;
   const double vfact = vbasis_->real() ? 1.0 : 0.5;
   const double omega = wf_.cell().volume();
   const double *const g2i = vbasis_->g2i_ptr();
   double tsum;
   double casino_vloc = 0.0;
   for ( int ig = 0; ig < ngloc; ig++ )
      casino_vloc += real(vion_local_g[ig] * conj(rhoelg[ig]));
   tsum = casino_vloc;
   vbasis_->context().dsum(1,1,&tsum,1);
   casino_vloc = tsum;

   //ewd DEBUG
   if ( false && s_.ctxt_.mype()==0 )
      cout << "EF.CASINO_VLOC:  tsum = " << tsum << endl;
   
   return casino_vloc;
}  
////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::print(ostream& os) const
{
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setprecision(8);
  os << "  <ekin>   " << setw(15) << ekin() << " </ekin>\n"
     << "  <econf>  " << setw(15) << econf() << " </econf>\n"
     << "  <eps>    " << setw(15) << eps() << " </eps>\n"
     << "  <enl>    " << setw(15) << enl() << " </enl>\n"
     << "  <ecoul>  " << setw(15) << ecoul() << " </ecoul>\n"
     << "  <exc>    " << setw(15) << exc() << " </exc>\n"
     << "  <esr>    " << setw(15) << esr() << " </esr>\n"
     << "  <eself>  " << setw(15) << eself() << " </eself>\n"
     << "  <ets>    " << setw(15) << ets() << " </ets>\n";
  if (s_.ctrl.enthalpy_pressure != 0.0)
    os << "  <epv>    " << setw(15) << epv() << " </epv>\n";
  if (s_.ctrl.dft_plus_u) 
    os << "  <ehub>   " << setw(15) << ehub() << " </ehub>\n";
  os << "  <etotal> " << setw(15) << etotal() << " </etotal>\n"
     << flush;
}
  
////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::print_memory(ostream&os, double& totsum, double& locsum) const
{
  int kmult = wf_.nspin()*wf_.nkp();
  int kmultloc = wf_.nkptloc();
  nlp[0][0]->print_memory(os,kmult,kmultloc,totsum,locsum);
}
  
////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const EnergyFunctional& e )
{ 
  e.print(os); 
  return os;
}
