////////////////////////////////////////////////////////////////////////////////
//
// ChargeDensity.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.C,v 1.33 2010/02/10 22:50:31 draeger1 Exp $

#include "ChargeDensity.h"
#include "Sample.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "FourierTransform.h"
#include "SlaterDet.h"
#include "SymmetrySet.h"
#include "AtomSet.h"
#include "Species.h"
#include "PrintMem.h"
#include "StructureFactor.h"
#include "blas.h"

#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ChargeDensity::ChargeDensity(Sample& s) : wf_(s.wf),atoms_(s.atoms),ctxt_(s.wf.context()) {
  ultrasoft_ = s.ctrl.ultrasoft;
  nlcc_ = s.ctrl.nlcc;
  highmem_ = false;
  if (s.ctrl.extra_memory >= 9)
    highmem_ = true;
  D3vector tmpkpoint(0.0,0.0,0.0);
  vcontext_ = wf_.sdloc(0)->basis().context();
  vbasis_ = new Basis(vcontext_, tmpkpoint, ultrasoft_);
  if (s.ctrl.ecutden > 0.0)
     vbasis_->resize(wf_.cell(),wf_.refcell(),s.ctrl.ecutden);
  else
     vbasis_->resize(wf_.cell(),wf_.refcell(),4.0*wf_.ecut());
  
  // define vft_, FT on vbasis context for transforming the density  
  int np0v = vbasis_->np(0);
  int np1v = vbasis_->np(1);
  int np2v = vbasis_->np(2);
  vft_ = new FourierTransform(*vbasis_,np0v,np1v,np2v);
  rhor.resize(wf_.nspin());
  rhog.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    rhor[ispin].resize(vft_->np012loc());
    rhog[ispin].resize(vbasis_->localsize());
  }
  rhotmp.resize(vft_->np012loc());

  if (nlcc_)
  {
     xcrhor.resize(wf_.nspin());
     xcrhog.resize(wf_.nspin());
     for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
        xcrhor[ispin].resize(vft_->np012loc());
        xcrhog[ispin].resize(vbasis_->localsize());
     }
  }
  
  // if nsym>0, identify symmetry-equivalent points for later averaging across procs
  nsym_ = s.symmetries.nsym();
  if (nsym_ > 0) {
    int np012loc = vft_->np012loc();
    symindexloc.resize(np012loc);
    symmultloc.resize(np012loc);
    for (int i=0; i<np012loc; i++) {
      symindexloc[i] = -1;
      symmultloc[i] = 0;
    }
    nsymgrp_ = 0;
    int np2first = vft_->np2_first(vcontext_.myrow());
    int np2last = np2first + vft_->np2_loc(vcontext_.myrow());

    for (int sy=0; sy<nsym_; sy++) {
      if (s.symmetries.symlist[sy]->setGrid(np0v,np1v,np2v)) {
        if ( vcontext_.oncoutpe() )
          cout << "<WARNING> symmetry " << sy << " does not match grid! </WARNING>" << endl;
      }
    }
    if ( vcontext_.oncoutpe() )
      cout << "<!-- nsym = " << nsym_ << " -->" << endl;

    vector<vector<vector<int > > > doneit;
    doneit.resize(np0v);
    for (int i=0; i<np0v; i++) 
      doneit[i].resize(np1v);
    for (int i=0; i<np0v; i++) 
      for (int j=0; j<np1v; j++) 
        doneit[i][j].resize(np2v);
    
    for (int i=0; i<np0v; i++) 
      for (int j=0; j<np1v; j++) 
        for (int k=0; k<np2v; k++) 
          doneit[i][j][k] = 0;

    for (int i=0; i<np0v; i++) {
      for (int j=0; j<np1v; j++) {
        for (int k=0; k<np2v; k++) {
          if (doneit[i][j][k] == 0) {
            doneit[i][j][k] = 1;

            // check if k is on local proc
            if (k >= np2first && k < np2last) {
              int locindex = (k-np2first)*np1v*np0v + j*np0v + i;
              symindexloc[locindex]=nsymgrp_;
              symmultloc[locindex]++;
            }

            for (int sy=0; sy<nsym_; sy++) {

              // map i,j,k to local index, set symindexloc[that] = nsymgrp_
              // symmultloc[that]++

              int istar,jstar,kstar;
              s.symmetries.symlist[sy]->applyToGridPoint(i,j,k,istar,jstar,kstar);
              doneit[istar][jstar][kstar] = 1;
              // check if kstar is on local proc
              if (kstar >= np2first && kstar < np2last) {
                int locindex = (kstar-np2first)*np1v*np0v + jstar*np0v + istar;
                symindexloc[locindex]=nsymgrp_;
                symmultloc[locindex]++;
              }
            }
            nsymgrp_++;
            
          }
        }
      }
    }

    if (vcontext_.oncoutpe()) 
      cout << "<!-- Finished mapping symmetric grid points:  nsymgrp = " << nsymgrp_ << ", np012 = " << np0v*np1v*np2v << ", vbasis.size = " << vbasis_->size() << ", nsym = " << nsym_ << " -->" << endl;

  }

  // FT for interpolation of wavefunctions on the fine grid
  ft_.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    ft_[ispin].resize(wf_.nkp());

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
    for (int kp=0; kp<wf_.nkp(); kp++)
      ft_[ispin][kp] = 0;

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          ft_[ispin][ikp] =
              new FourierTransform(wf_.sd(ispin,ikp)->basis(),np0v,np1v,np2v);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
ChargeDensity::~ChargeDensity(void) {
  delete vbasis_;
  delete vft_;
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) 
    for ( int ikp = 0; ikp < ft_[ispin].size(); ikp++ )
      delete ft_[ispin][ikp];
  //ewd clean up output -- add back in with verbose or timing option?
  /*
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    //if ( ctxt_.myproc()==0 ) {
    if ( ctxt_.mype()==0 ) {
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
  */
}

////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_density() {
  assert(rhor.size() == wf_.nspin());
  const double omega = vbasis_->cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;
  vector<complex<double> > rhogus;
  
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    assert(rhor[ispin].size() == vft_->np012loc() );
    assert(rhotmp.size() == vft_->np012loc() );
    double ncnorm;
    
    for ( int i = 0; i < vft_->np012loc(); i++ )
      rhor[ispin][i] = 0.0;

    // compute rhor contribution from wf orbitals
    if ( wf_.spinactive(ispin)) {
      tmap["charge_compute"].start();
      // on each sdcontext, rhor is filled with the weighted contributions to 
      // the total density, for each state and kpoint
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          double wt = wf_.weight(ikp)/wf_.weightsum();
          // add contribution for this kpoint to rhor
          wf_.sd(ispin,ikp)->compute_density(*ft_[ispin][ikp],wt, &rhor[ispin][0]);
        }
      }
      tmap["charge_compute"].stop();
    }

    // sum across rows of wfcontext to calculate total density on both spincontexts
    tmap["charge_rowsum"].start();
    wf_.wfcontext()->dsum('r',vft_->np012loc(),1,&rhor[ispin][0],vft_->np012loc());
    tmap["charge_rowsum"].stop();

    // add ultrasoft contribution to charge density in reciprocal space
    //ewd note:  this may need to be tweaked to work with spin/k-points
    if (ultrasoft_) { 
      const int nsp = atoms_.nsp();
      const int ngwl = vbasis_->localsize();
      const int mloc = vbasis_->maxlocalsize();
      vector<vector<double> > tau;
      atoms_.get_positions(tau,true);
      
      rhogus.resize(ngwl);
      for (int ig=0; ig<ngwl; ig++)
        rhogus[ig] = complex<double>(0.0,0.0);
      
      // calc summat = SUM occ_i <beta|phi_i> <phi_i|beta>
      if ( wf_.spinactive(ispin)) {
        for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
          if (wf_.kptactive(ikp)) {
            assert(wf_.sd(ispin,ikp) != 0);
            double wt = wf_.weight(ikp)/wf_.weightsum();
            SlaterDet* sdp = wf_.sd(ispin,ikp);
            const Basis& basis_ = sdp->basis();
            const vector<double>& occ = sdp->occ();
            //ewd OPT
            //sdp->calc_betapsi();
            //ewd NO-OPT: remove all OPT comments for now
            sdp->calc_betapsi();
            for (int is=0; is<nsp; is++) {
              const ComplexMatrix* bpsi = sdp->betapsi(is);
              const complex<double>* bp = bpsi->cvalptr();  
              const int nbase = sdp->context().mycol() * sdp->c().nb();
              Species *s = atoms_.species_list[is];
              if (s->ultrasoft()) { 
                int naloc_t = atoms_.usloc_atind_t[is].size();
                int na = atoms_.na(is);
                int nqtot = s->nqtot();
                int nbeta = s->nbeta();
                int nbetalm = s->nbetalm();
                int nstloc = bpsi->nloc();
                int bp_mloc = bpsi->mloc();
                vector<complex<double> > summat(nqtot*na);
                for (int i=0; i<nqtot*na; i++)
                  summat[i] = complex<double>(0.0,0.0);
               
                if (bpsi->size() > 0) {
                  for (int ibl = 0; ibl < naloc_t; ibl++) {
                    int ia = atoms_.usloc_atind_t[is][ibl]; // ia = absolute atom index of betapsi local ind
                    for (int qind=0; qind < nqtot; qind++) {
                      int lm1,lm2;
                      s->qind_to_betalm(qind,lm1,lm2);
                      double mult = 2.0;          // off-diag. terms
                      if (lm1 == lm2) mult = 1.0;
                      for (int n=0; n<nstloc; n++) {
                        const double sdocc = occ[n + nbase];
                        int ind1 = bp_mloc*n + ibl*nbetalm + lm1;
                        int ind2 = bp_mloc*n + ibl*nbetalm + lm2;
                        summat[ia*nqtot + qind] += omega_inv*wt*mult*sdocc*(conj(bp[ind1])*bp[ind2]);
                      }
                    }
                  }
                }
                int dsize = 2*na*nqtot;
                double* psummat = (double*)&summat[0];
                // sum over process rows = sum over all states
                sdp->context().dsum('r',dsize,1,&psummat[0],dsize);
                // sum over process columns = summat for all atoms now on all pes
                sdp->context().dsum('c',dsize,1,&psummat[0],dsize);
                
                if (highmem_) {
                  // multiply summat w. qnmg, then sum across rows to get ultrasoft
                  // contribution to charge density
                  complex<double> zone = complex<double>(1.0,0.0);
                  int one = 1;
                  char cn='n';
                  int sumlocsize = nqtot*atoms_.usloc_nat[is]; // want to multiply chunk of summat that matches local qnmg 
                  int sumind0 = nqtot*atoms_.usloc_atind[is][0];
                  if (sumlocsize > 0 && ngwl > 0)
                    zgemm(&cn,&cn,(int*)&ngwl,&one,&sumlocsize,&zone,&qnmg_[is][0],
                          (int*)&ngwl,&summat[sumind0],&sumlocsize,&zone,&rhogus[0],
                          (int*)&ngwl);
                }
                else {
                  // multiply qnm(G), structure factor and summat
                  const double *const gx = vbasis_->gx_ptr(0);
                  const double *const gy = vbasis_->gx_ptr(1);
                  const double *const gz = vbasis_->gx_ptr(2);
                  int naloc = atoms_.usloc_nat[is];
                  int ialoc0 = atoms_.usloc_atind[is][0];
                  for (int ig=0; ig<ngwl; ig++) {
                    for (int ialoc=0; ialoc<naloc; ialoc++) {
                      int ia = ialoc0 + ialoc;
                      const double arg = tau[is][3*ia]*gx[ig] + tau[is][3*ia+1]*gy[ig] + tau[is][3*ia+2]*gz[ig];
                      double sgr = sin(arg);
                      double cgr = cos(arg);
                      complex<double> rhosum = complex<double>(0.0,0.0);
                      for (int qind=0; qind<nqtot; qind++)
                        rhosum += qnmg_[is][qind*ngwl+ig]*summat[nqtot*ia+qind];
                      rhogus[ig] += omega_inv*rhosum*complex<double>(cgr,-sgr);
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      // sum rhogus across rows to average ultrasoft contribution over all atoms, kpts, spin
      int sumsize = 2*ngwl;
      double* rp = (double*)&rhogus[0];
      wf_.wfcontext()->dsum('r',sumsize,1,&rp[0],sumsize);

      // transform ultrasoft density contribution to real space
      vft_->backward(&rhogus[0],&rhotmp[0]);

      double uscharge_ = 0.0;
      for ( int i = 0; i < rhor[ispin].size(); i++ )
        uscharge_ += rhotmp[i].real();
      uscharge_ *= omega/(double)vft_->np012();
      vcontext_.dsum('c',1,1,&uscharge_,1);
      if ( wf_.spincontext(ispin)->myproc() == 0 ) 
        cout << "<!-- ultrasoft contribution to total charge = " << uscharge_ << " -->" << endl;

      // add ultrasoft contribution to rhor
      const int rhor_size = rhor[ispin].size();
      double *const prhor = &rhor[ispin][0];
      for ( int i = 0; i < rhor_size; i++ ) 
        prhor[i] += rhotmp[i].real();
    }

    tmap["charge_sym"].start();
    if (nsym_ > 0) {
      int np012loc = vft_->np012loc();

      int subset = 16777216; // if total density exceeds 128 megs, break sum into parts
      double nsyminv = 1./((double)nsym_+1.);
      if (nsymgrp_ <= subset) {
        double symrho[nsymgrp_];
        for (int i=0; i<nsymgrp_; i++)
          symrho[i] = 0.0;
      
        const double *const prhor = &rhor[ispin][0];
        for (int i=0; i<np012loc; i++) 
          symrho[symindexloc[i]] += nsyminv*symmultloc[i]*prhor[i];

        //wf_.spincontext(ispin)->dsum('c',nsymgrp_,1,&symrho[0],nsymgrp_);
        vcontext_.dsum('c',nsymgrp_,1,&symrho[0],nsymgrp_);
        for (int i=0; i<np012loc; i++)
          rhor[ispin][i] = symrho[symindexloc[i]];
      }
      else {
        int nsubsets = nsymgrp_/subset;
        if (nsymgrp_%subset != 0) nsubsets++;
        double symrho[subset];
        for (int n=0; n<nsubsets; n++) {
          int min = n*subset;
          int max = min+subset;
          if (max > nsymgrp_) max = nsymgrp_;
          int size = max-min+1;

          for (int i=0; i<subset; i++)
            symrho[i] = 0.0;
          const double *const prhor = &rhor[ispin][0];
          for (int i=0; i<np012loc; i++)
            if (symindexloc[i] >= min && symindexloc[i] < max)
              symrho[symindexloc[i]-min] += nsyminv*symmultloc[i]*prhor[i];
            
          vcontext_.dsum('c',size,1,&symrho[0],size);
          for (int i=0; i<np012loc; i++)
            if (symindexloc[i] >= min && symindexloc[i] < max)
              rhor[ispin][i] = symrho[symindexloc[i]-min];
        }
      }
    }
    tmap["charge_sym"].stop();
    
    // check integral of charge density
    // compute Fourier coefficients of the charge density
    double sum = 0.0;
    const int rhor_size = rhor[ispin].size();
    const double *const prhor = &rhor[ispin][0];
    tmap["charge_integral"].start();

    for ( int i = 0; i < rhor_size; i++ ) {
      const double prh = prhor[i];
      sum += prh;
      rhotmp[i] = complex<double>(omega * prh, 0.0);
    }
    sum *= omega / vft_->np012();
    
    if ( wf_.spinactive(ispin)) {
      wf_.spincontext(ispin)->dsum('c',1,1,&sum,1);
      if ( wf_.spincontext(ispin)->myproc() == 0 ) {
        cout.setf(ios::fixed,ios::floatfield);
        cout.setf(ios::right,ios::adjustfield);
        cout << "  <!-- total_electronic_charge: " << setprecision(8) << sum 
             << ", spin = " << ispin << " -->" << endl;
        // print out warning if total charge is not an integer value
        const double chgthresh = 1.E-7;
        int chgint = (int) (sum+chgthresh);
        double nonint = abs( sum - (double) chgint);
        if (nonint > chgthresh) 
          cout << "<WARNING> Total electronic charge has non-integer value!! </WARNING>" << endl;
      }
    }
    tmap["charge_vft"].start();
    vft_->forward(&rhotmp[0],&rhog[ispin][0]);
    tmap["charge_vft"].stop();
  
  }
  if (nlcc_)
     add_nlccden();
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_rhor(void) {
  // recalculate rhor from rhog
  assert(rhor.size() == wf_.nspin());
  const double omega = vbasis_->cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    assert(rhor[ispin].size() == vft_->np012loc() );
    assert(rhotmp.size() == vft_->np012loc() );
      
    vft_->backward(&rhog[ispin][0],&rhotmp[0]);

    const int rhor_size = rhor[ispin].size();
    double *const prhor = &rhor[ispin][0];
    for ( int i = 0; i < rhor_size; i++ ) {
      prhor[i] = rhotmp[i].real() * omega_inv;
    }
  }
  if (nlcc_)
     add_nlccden();
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::reshape_rhor(const Context& oldvctxt, const Context& newvctxt) {
  // when nrowmax is changed, sd objects are rebuilt, often with new vcontexts, 
  // thus rhor is no longer distributed correctly.  This routine assumes that rhor
  // contains the correct values, redistributes them on the new context and then
  // updates rhog from rhor.
  // 
  // Note that this routine is currently incomplete -- although
  // vbasis, vft, rhotmp, rhor and rhog objects are updated, the ft_
  // array and symmetry arrays are not.

  //ewd SPIN
  assert(wf_.nspin() == 1);
  
  const Context* wfctxt = wf_.spincontext(0);
  vcontext_ = newvctxt;
  D3vector tmpkpoint(0.0,0.0,0.0);
  Basis* oldvbasis_ = new Basis(oldvctxt, tmpkpoint);
  double ecut = vbasis_->ecut();
  oldvbasis_->resize(wf_.cell(),wf_.refcell(),ecut);
  int oldnp0v = oldvbasis_->np(0);
  int oldnp1v = oldvbasis_->np(1);
  int oldnp2v = oldvbasis_->np(2);
  FourierTransform* oldvft_ = new FourierTransform(*oldvbasis_,oldnp0v,oldnp1v,oldnp2v);

  if (vbasis_ != 0) 
    delete vbasis_;
  vbasis_ = new Basis(vcontext_, tmpkpoint);
  vbasis_->resize(wf_.cell(),wf_.refcell(),ecut);
  int np0v = vbasis_->np(0);
  int np1v = vbasis_->np(1);
  int np2v = vbasis_->np(2);
  if (vft_ != 0) 
    delete vft_;
  vft_ = new FourierTransform(*vbasis_,np0v,np1v,np2v);

  // collect rhor on first process of each column of old context
  int oldrow = oldvctxt.myrow();
  int oldcol = oldvctxt.mycol();
  vector<double> rhortmp;
  if (oldrow == 0)
    rhortmp.resize(oldvft_->np012());

  for ( int i = 0; i < oldvctxt.nprow(); i++ ) {    
    if ( i == oldrow ) {
      int size = oldvft_->np012loc();
      oldvctxt.isend(1,1,&size,1,0,oldcol);
      if (size > 0) {
        vector<double> rhorsend(size);
        for (int j = 0; j < size; j++)
          rhorsend[j] = rhor[0][j];
        oldvctxt.dsend(size,1,&rhorsend[0],1,0,oldcol);
      }
    }
  }
  if ( oldrow == 0) {
    int offset = 0;
    for ( int i = 0; i < oldvctxt.nprow(); i++ ) {
      int size = 0;
      oldvctxt.irecv(1,1,&size,1,i,oldcol);
      if (size > 0) {
        oldvctxt.drecv(size,1,&rhortmp[offset],1,i,oldcol);
        offset += size;
      }
    }
    assert(offset == oldvft_->np012());
  }

  // rhor is now complete on pe0:  redistribute on new context
  int newrow = wfctxt->myrow();
  int newcol = wfctxt->mycol();
  int newsize = vft_->np012loc();
  if (newrow == 0)
    rhortmp.resize(vft_->np012());

  // send full density to first pe of columns on new context
  if (wfctxt->onpe0()) {
    for ( int i = 1; i < wfctxt->npcol(); i++ ) {
      int fullsize = vft_->np012();
      wfctxt->isend(1,1,&fullsize,1,0,i);
      wfctxt->dsend(fullsize,1,&rhortmp[0],1,0,i);
    }
  }      
  if (newrow == 0 && newcol > 0) {
    int tmpsize;
    wfctxt->irecv(1,1,&tmpsize,1,0,0);
    assert(tmpsize == vft_->np012());
    wfctxt->drecv(tmpsize,1,&rhortmp[0],1,0,0);
  }

  // now distribute density across rows
  wfctxt->isend(1,1,&newsize,1,0,newcol);
  if (newrow == 0) {
    int offset = 0;
    for ( int i = 0; i < wfctxt->nprow(); i++ ) {
      int tmpsize;
      wfctxt->irecv(1,1,&tmpsize,1,i,newcol);
      vector<double> rhorsend(tmpsize);
      for (int j = 0; j < tmpsize; j++)
        rhorsend[j] = rhortmp[offset+j];
      wfctxt->isend(1,1,&tmpsize,1,i,newcol);
      wfctxt->dsend(tmpsize,1,&rhorsend[0],1,i,newcol);
      offset += tmpsize;
    }
  }

  for ( int i = 0; i < wfctxt->nprow(); i++ ) {
    if (i == newrow) { 
      vector<double> rhorrecv(newsize);
      rhor[0].resize(newsize);
      int tmp2;
      wfctxt->irecv(1,1,&tmp2,1,0,newcol);
      assert(tmp2==newsize);
      wfctxt->drecv(newsize,1,&rhorrecv[0],1,0,newcol);
      for (int j = 0; j < newsize; j++)
        rhor[0][j] = rhorrecv[j];
    }
  }

  // recalculate rhog from rhor
  rhotmp.resize(newsize);
  int rhogsize = vbasis_->localsize();
  rhog[0].resize(rhogsize);
  const double omega = wf_.cell().volume();
  const int rhor_size = rhor[0].size();
  for ( int i = 0; i < rhor_size; i++ ) {
    rhotmp[i] = complex<double>(omega * rhor[0][i], 0.0);
  }
  vector<complex<double> > rhogtmp(rhogsize);
  vft_->forward(&rhotmp[0],&rhogtmp[0]);
  complex<double> *rhogp = &rhog[0][0];
  for (int j = 0; j < rhogsize; j++)
    rhogp[j] = rhogtmp[j];

  delete oldvbasis_;
  delete oldvft_;

  if (nlcc_)
     update_nlcc();

}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_usfns() {
  // fill Q_nm(G), multiply structure factor term

  const int nsp = atoms_.nsp();
  const int ngwl = vbasis_->localsize();
  const double omega = vbasis_->cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  if (highmem_) { 
    vector<vector<double> > tau;
    atoms_.get_positions(tau,true);

    qnmg_.resize(nsp);
    for (int is=0; is<nsp; is++) {
      Species *s = atoms_.species_list[is];
      if (s->ultrasoft()) { 
        int naloc = atoms_.usloc_nat[is];
        int nqtot = s->nqtot();
        qnmg_[is].resize(nqtot*naloc*ngwl);
        for (int ibl = 0; ibl < naloc; ibl++) {
          int ia = atoms_.usloc_atind[is][ibl];

          // calculate structure factor for Q_nm(G)
          vector<double> cgr(ngwl); 
          vector<double> sgr(ngwl); 
          const double *const kpgx = vbasis_->kpgx_ptr(0);
          const double *const kpgy = vbasis_->kpgx_ptr(1);
          const double *const kpgz = vbasis_->kpgx_ptr(2);
          for ( int ig = 0; ig < ngwl; ig++ ) {
            const double arg = tau[is][3*ia]*kpgx[ig] + tau[is][3*ia+1]*kpgy[ig] + tau[is][3*ia+2]*kpgz[ig];
            sgr[ig] = sin(arg);
            cgr[ig] = cos(arg);
          }

          vector<double> qaug;
          vector<complex<double> > qnm;
          s->calc_qnmg(vbasis_,qnm,qaug);

          for (int qind=0; qind < nqtot; qind++) {
            const int ind0 = ibl*nqtot*ngwl + qind*ngwl;
            for ( int ig = 0; ig < ngwl; ig++ ) 
              // qnmg needs inverse volume factor from FFT
              qnmg_[is][ind0+ig] = omega_inv*qnm[qind*ngwl+ig]*complex<double>(cgr[ig],-sgr[ig]);
          }
        }
      }
    }
  }
  else {
    qnmg_.resize(nsp);
    for (int is=0; is<nsp; is++) {
      Species *s = atoms_.species_list[is];
      if (s->ultrasoft()) { 
        int nqtot = s->nqtot();
        qnmg_[is].resize(nqtot*ngwl);
        vector<double> qaug;
        s->calc_qnmg(vbasis_,qnmg_[is],qaug);
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_nlcc() {
  // sum up nlcc core density*sfact(ia)

  const double omega = vbasis_->cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  vector<vector<double> > tau;
  atoms_.get_positions(tau,true);
  StructureFactor sf;
  sf.init(tau,*vbasis_);
  sf.update(tau,*vbasis_);

  rhornlcc_.resize(wf_.nspin());
  rhognlcc.resize(wf_.nspin());
  const int ngwl = vbasis_->localsize();
  const double* gptr = vbasis_->g_ptr();
  for (int ispin=0; ispin<wf_.nspin(); ispin++)
  {
     if ( wf_.spinactive(ispin)) 
     {
        rhognlcc[ispin].resize(ngwl);
        memset( (void*)&rhognlcc[ispin][0], 0, 2*ngwl*sizeof(double) );
        for (int is=0; is<atoms_.nsp(); is++) {
           Species *s = atoms_.species_list[is];
           complex<double> *sfac = &sf.sfac[is][0];
           if (s->nlcc()) { 
              for ( int ig = 0; ig < ngwl; ig++ ) {
                 double gval = gptr[ig];
                 double rhogcore = s->rhog_nlcc(gval);
                 // need factor of omega_inv for rhog transform, canceled by sfac
                 rhognlcc[ispin][ig] += rhogcore*sfac[ig];
              }
           }
        }
      
        vft_->backward(&rhognlcc[ispin][0],&rhotmp[0]);

        const int rhor_size = rhor[ispin].size();
        rhornlcc_[ispin].resize(rhor_size);
        for ( int i = 0; i < rhor_size; i++ )
           // ewd: should we follow PWSCF's example of forcing positive FFT(rhognlcc)?
           //rhornlcc_[ispin][i] = fabs(rhotmp[i].real())*omega_inv;
           rhornlcc_[ispin][i] = rhotmp[i].real()*omega_inv;
     }
  }

  add_nlccden();
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::add_nlccden()
{
   if (rhognlcc.size() != rhog.size())
      update_nlcc();
   else if (rhornlcc_.size() != rhor.size())
      update_nlcc();
   else if (rhognlcc[0].size() != rhog[0].size())
      update_nlcc();
   else if (rhornlcc_[0].size() != rhor[0].size())
      update_nlcc();
   
   const int ngwl = vbasis_->localsize();
   for (int ispin=0; ispin<wf_.nspin(); ispin++)
   {
      if ( wf_.spinactive(ispin)) 
      {
         for ( int ig = 0; ig < ngwl; ig++ )
            xcrhog[ispin][ig] = rhog[ispin][ig] + rhognlcc[ispin][ig];
         const int rhor_size = rhor[ispin].size();
         for ( int i = 0; i < rhor_size; i++ )
            xcrhor[ispin][i] = rhor[ispin][i] + rhornlcc_[ispin][i];
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::calc_drhogus(vector<vector<complex<double> > > &drhogus)
{
  // compute derivative of rhog with respect to ion coordinates
  if (ultrasoft_) { 
    const int nsp = atoms_.nsp();
    const int ngwl = vbasis_->localsize();
    const int mloc = vbasis_->maxlocalsize();
    const double omega = vbasis_->cell().volume();
    assert(omega != 0.0);
    const double omega_inv = 1.0 / omega;
    vector<vector<double> > tau;
    atoms_.get_positions(tau,true);

    drhogus.resize(3);
    for (int j=0; j<3; j++)
      drhogus[j].resize(mloc);
    
    for (int j=0; j<3; j++)
      for (int ig=0; ig<mloc; ig++)
        drhogus[j][ig] = complex<double>(0.0,0.0);

    // calc summat = SUM occ_i <beta|phi_i> <phi_i|beta>
    for (int ispin=0; ispin<wf_.nspin(); ispin++) {
      if ( wf_.spinactive(ispin)) {
        for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
          if (wf_.kptactive(ikp)) {
            assert(wf_.sd(ispin,ikp) != 0);
            double wt = wf_.weight(ikp)/wf_.weightsum();
            SlaterDet* sdp = wf_.sd(ispin,ikp);
            const Basis& basis_ = sdp->basis();
            const vector<double>& occ = sdp->occ();
            //ewd OPT
            //sdp->calc_betapsi();
            //ewd NO-OPT
            sdp->calc_betapsi();
            for (int is=0; is<nsp; is++) {
              const ComplexMatrix* bpsi = sdp->betapsi(is);
              const complex<double>* bp = bpsi->cvalptr();  
              const int nbase = sdp->context().mycol() * sdp->c().nb();
              Species *s = atoms_.species_list[is];
              if (s->ultrasoft()) { 
                int naloc_t = atoms_.usloc_atind_t[is].size();
                int na = atoms_.na(is);
                int nqtot = s->nqtot();
                int nbeta = s->nbeta();
                int nbetalm = s->nbetalm();
                int nstloc = bpsi->nloc();
                int bp_mloc = bpsi->mloc();

                for (int j=0; j<3; j++) {          // force component 
                  vector<complex<double> > summat(nqtot*na);
                  vector<complex<double> > dsummat(nqtot*na);
                  for (int i=0; i<nqtot*na; i++) {
                    summat[i] = complex<double>(0.0,0.0);
                    dsummat[i] = complex<double>(0.0,0.0);
                  }
              
                  sdp->calc_dbetapsi(j);
                  if (bpsi->size() > 0) {
                    const ComplexMatrix* dbetapsi = sdp->dbetapsi(is);
                    const complex<double>* dbpj = dbetapsi->cvalptr();
                    for (int ibl = 0; ibl < naloc_t; ibl++) {
                      int ia = atoms_.usloc_atind_t[is][ibl]; // ia = absolute atom index of betapsi local ind
                      for (int qind=0; qind < nqtot; qind++) {
                        int lm1,lm2;
                        s->qind_to_betalm(qind,lm1,lm2);
                        double mult = 2.0;          // off-diag. terms
                        if (lm1 == lm2) mult = 1.0;
                        for (int n=0; n<nstloc; n++) {
                          const double sdocc = occ[n + nbase];
                          int ind1 = bp_mloc*n + ibl*nbetalm + lm1;
                          int ind2 = bp_mloc*n + ibl*nbetalm + lm2;
                          summat[ia*nqtot + qind] += omega_inv*wt*mult*sdocc*(conj(bp[ind1])*bp[ind2]);
                          dsummat[ia*nqtot + qind] += omega_inv*wt*mult*sdocc*(conj(dbpj[ind1])*bp[ind2] + conj(bp[ind1])*dbpj[ind2]);
                        }
                      }
                    }
                  }
                  int dsize = 2*na*nqtot;
                  double* psummat = (double*)&summat[0];
                  double* pdsummat = (double*)&dsummat[0];
                  // sum over process rows = sum over all states
                  sdp->context().dsum('r',dsize,1,&psummat[0],dsize);
                  sdp->context().dsum('r',dsize,1,&pdsummat[0],dsize);
                  // sum over process columns = summat for all atoms now on all pes
                  sdp->context().dsum('c',dsize,1,&psummat[0],dsize);
                  sdp->context().dsum('c',dsize,1,&pdsummat[0],dsize);
                  if (highmem_) {
                    complex<double> zone = complex<double>(1.0,0.0);
                    int one = 1;
                    char cn='n';
                    // want to multiply chunk of summat that matches local qnmg 
                    int sumlocsize = nqtot*atoms_.usloc_nat[is]; 
                    int sumind0 = nqtot*atoms_.usloc_atind[is][0];
                    if (sumlocsize > 0 && ngwl > 0)
                      zgemm(&cn,&cn,(int*)&ngwl,&one,&sumlocsize,&zone,&qnmg_[is][0],
                            (int*)&ngwl,&dsummat[sumind0],&sumlocsize,&zone,&drhogus[j][0],
                            (int*)&ngwl);

                    // add dqnm(G)/dR
                    int naloc = atoms_.usloc_nat[is];
                    int ialoc0 = atoms_.usloc_atind[is][0];
                    const double *const gxj = vbasis_->gx_ptr(j);
                    for (int ig=0; ig<ngwl; ig++) {
                      complex<double> dsf = complex<double>(0.0,-gxj[ig]);
                      for (int ialoc=0; ialoc<naloc; ialoc++) {
                        int ia = ialoc0 + ialoc;
                        complex<double> rhosum = complex<double>(0.0,0.0);
                        for (int qind=0; qind<nqtot; qind++)
                          drhogus[j][ig] += dsf*qnmg_[is][qind*ngwl+ig]*summat[nqtot*ia+qind];
                      }
                    }
                  }
                  else {
                    // multiply qnm(G), structure factor and dsummat, add dqnm(G)/dR
                    const double *const gx = vbasis_->gx_ptr(0);
                    const double *const gy = vbasis_->gx_ptr(1);
                    const double *const gz = vbasis_->gx_ptr(2);
                    const double *const gxj = vbasis_->gx_ptr(j);
                    int naloc = atoms_.usloc_nat[is];
                    int ialoc0 = atoms_.usloc_atind[is][0];
                    for (int ig=0; ig<ngwl; ig++) {
                      complex<double> dsf = complex<double>(0.0,-gxj[ig]);
                      for (int ialoc=0; ialoc<naloc; ialoc++) {
                        int ia = ialoc0 + ialoc;
                        const double arg = tau[is][3*ia]*gx[ig] + tau[is][3*ia+1]*gy[ig] + tau[is][3*ia+2]*gz[ig];
                        double sgr = sin(arg);
                        double cgr = cos(arg);
                        complex<double> rhosum = complex<double>(0.0,0.0);
                        for (int qind=0; qind<nqtot; qind++)
                          rhosum += qnmg_[is][qind*ngwl+ig]*dsummat[nqtot*ia+qind] + dsf*qnmg_[is][qind*ngwl+ig]*summat[nqtot*ia+qind];
                        drhogus[j][ig] += omega_inv*rhosum*complex<double>(cgr,-sgr);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // sum drhogus across rows to average ultrasoft contribution over all atoms, kpts, spin
    for (int j=0; j<3; j++) {
      int sumsize = 2*mloc;
      double* rp = (double*)&drhogus[j][0];
      wf_.wfcontext()->dsum('r',sumsize,1,&rp[0],sumsize);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::print_memory(ostream& os, double& totsum, double& locsum) const
{
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setprecision(3);
  PrintMem pm;

  const int ngw = vbasis_->size();
  const int ngwl = vbasis_->localsize();
  const int np012 = vft_->np012();
  const int np012loc = vft_->np012loc();
  const int nspin = wf_.nspin();
  
  double rhor_size = (double)(wf_.wfcontext()->npcol()*np012*nspin*sizeof(double));
  double rhor_locsize = (double)(np012loc*sizeof(double));
  double rhog_size = (double)(wf_.wfcontext()->npcol()*ngw*nspin*sizeof(complex<double>));
  double rhog_locsize = (double)(ngwl*sizeof(complex<double>));

  double qnmg_size = 0.0;
  double qnmg_locsize = 0.0;
  if (ultrasoft_) {
    const int nsp = atoms_.nsp();
    for (int is=0; is<nsp; is++) {
      Species *s = atoms_.species_list[is];
      if (s->ultrasoft()) { 
        int na = atoms_.na(is);
        int naloc = atoms_.usloc_atind[is].size();
        int nqtot = s->nqtot();
        if (highmem_) {
          qnmg_size += (double)(nspin*nqtot*na*ngw*sizeof(complex<double>));
          qnmg_locsize += (double)(nqtot*naloc*ngwl*sizeof(complex<double>));
        }
        else {
          qnmg_size += (double)(nspin*nqtot*ngw*sizeof(complex<double>));
          qnmg_locsize += (double)(nqtot*ngwl*sizeof(complex<double>));
        }
      }
    }
  }

  totsum += rhor_size + rhog_size;
  locsum += rhor_locsize + rhog_locsize;
  if (ultrasoft_) {
    totsum += qnmg_size;
    locsum += qnmg_locsize;
  }
  if (nlcc_) {
     totsum += rhor_size + rhog_size;
     locsum += rhor_locsize + rhog_locsize;
  }
  
  string rhor_unit = pm.memunit(rhor_size);
  string rhor_locunit = pm.memunit(rhor_locsize);
  os << "<!-- memory cd.rhor     :  " << setw(7) << rhor_size << rhor_unit << "  (" << rhor_locsize << rhor_locunit << " local) -->" << endl;
  string rhog_unit = pm.memunit(rhog_size);
  string rhog_locunit = pm.memunit(rhog_locsize);
  os << "<!-- memory cd.rhog     :  " << setw(7) << rhog_size << rhog_unit << "  (" << rhog_locsize << rhog_locunit << " local) -->" << endl;

  if (nlcc_) {
     os << "<!-- memory cd.nlccrhor :  " << setw(7) << rhor_size << rhor_unit << "  (" << rhor_locsize << rhor_locunit << " local) -->" << endl;
     os << "<!-- memory cd.nlccrhog :  " << setw(7) << rhog_size << rhog_unit << "  (" << rhog_locsize << rhog_locunit << " local) -->" << endl;
  }
  
  if (ultrasoft_) {
    string qnmg_unit = pm.memunit(qnmg_size);
    string qnmg_locunit = pm.memunit(qnmg_locsize);
    os << "<!-- memory cd.qnmg     :  " << setw(7) << qnmg_size << qnmg_unit << "  (" << qnmg_locsize << qnmg_locunit << " local) -->" << endl;
  }

}
