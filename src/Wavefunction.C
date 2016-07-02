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
// Wavefunction.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Wavefunction.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Basis.h"
#include "Context.h"
#include "jacobi.h"
#include "profile.h"
#include <vector>
#include <iomanip>
#include <sstream>
#if USE_CSTDIO_LFS
#include <cstdio>
#endif
#include "fstream"
#ifdef BGQ
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Wavefunction::Wavefunction(const Context& ctxt) : ctxt_(ctxt),nel_(0),nempty_(0),
                                                  nspin_(1), deltaspin_(0),
                                                  ecut_(0.0), nrowmax_(32),
                                                  nparallelkpts_(0), mu_(0.0),
                                                  deltacharge_(0.0)
{
  // create a default wavefunction: one k point, k=0
  kpoint_.resize(1);
  kpoint_[0] = D3vector(0,0,0);
  nkptloc_ = 1;
  kpt_added_ = false;
  kptloc_.resize(nkptloc_);
  kptloc_[0] = 0;
  spinloc_ = 0;
  weight_.resize(1);
  weight_[0] = 1.0;
  weightsum_ = 1.0;
  compute_nst();
  hasdata_ = false;
  mbset_ = -1;
  nbset_ = -1;
  mblks_ = 1;
  nblks_ = 1;
  //ewdallocate();
  ultrasoft_ = false;
  force_complex_wf_ = false;
  wf_phase_real_ = false;
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::Wavefunction(const Wavefunction& wf) : ctxt_(wf.ctxt_), 
nel_(wf.nel_), nempty_(wf.nempty_), nspin_(wf.nspin_), 
deltaspin_(wf.deltaspin_), nrowmax_(wf.nrowmax_), 
nparallelkpts_(wf.nparallelkpts_), kpt_added_(wf.kpt_added_),
nkptloc_(wf.nkptloc_), spinloc_(wf.spinloc_),
cell_(wf.cell_), refcell_(wf.refcell_), 
ecut_(wf.ecut_), weightsum_(wf.weightsum_), 
ultrasoft_(wf.ultrasoft_), force_complex_wf_(wf.force_complex_wf_),
wf_phase_real_(wf.wf_phase_real_),mbset_(wf.mbset_),nbset_(wf.nbset_),
mblks_(wf.mblks_),nblks_(wf.nblks_),mu_(wf.mu_),deltacharge_(wf.deltacharge_)
{
  // Create a Wavefunction using the dimensions of the argument
  compute_nst();
  
  // Next lines: do special allocation of contexts to ensure that 
  // contexts are same as those of wf
  
  // create sd contexts and SlaterDets
  assert(ctxt_.active());
  
  int nkppar = wf.sdcontextsize(0);
  spincontext_.resize(nspin_);
  sdcontext_.resize(nspin_);
  sd_.resize(nspin_);

  sdcontextsq_.resize(nspin_);
  for (int ispin=0; ispin<nspin_; ispin++) {
    sdcontext_[ispin].resize(nkppar);
    sdcontextsq_[ispin].resize(nkppar);
    sd_[ispin].resize(wf.nkp());
  }
  
  kpoint_.resize(wf.nkp());
  weight_.resize(wf.nkp());
  for (int ispin=0; ispin<nspin_; ispin++) 
    for (int kp=0; kp<wf.nkp(); kp++) 
      sd_[ispin][kp] = 0;
  
  for (int kp=0; kp<wf.nkp(); kp++) {
    kpoint_[kp] = wf.kpoint(kp);
    weight_[kp] = wf.weight(kp);
  }
  kptloc_.resize(nkptloc_);
  for ( int kloc=0; kloc<wf.nkptloc(); kloc++)
    kptloc_[kloc] = wf.kptloc(kloc);
  mysdctxt_.resize(wf.nkp());

  // store number of local k-points on each sdcontext on all procs
  nkptloc_list_.resize(nkppar);
  for (int k=0; k<nkppar; k++)
    nkptloc_list_[k] = wf.nkptloc(k);

  wfcontext_ = new Context(*wf.wfcontext());
  for ( int ispin = 0; ispin < nspin(); ispin++ ) {
    spincontext_[ispin] = 0;
    if (wf.spinactive(ispin)) {
      spincontext_[ispin] = new Context(*wf.spincontext(ispin));
      for ( int ikp = 0; ikp < wf.sdcontextsize(ispin); ikp++ ) {
        sdcontext_[ispin][ikp] = 0;
        if (wf.sdcontext(ispin,ikp) != 0 ) {
          if (wf.sdcontext(ispin,ikp)->active() ) {
            sdcontext_[ispin][ikp] = new Context(*wf.sdcontext(ispin,ikp));
            sdcontextsq_[ispin][ikp] = new Context(*wf.sdcontextsq(ispin,ikp));


            /*
            Context* my_col_ctxt = 0;
            for ( int icol = 0; icol < sdcontext_[ispin][ikp]->npcol(); icol++ ) {
              Context* col_ctxt = new Context(*sdcontext_[ispin][ikp],sdcontext_[ispin][ikp]->nprow(),1,0,icol);
              sdcontext_[ispin][ikp]->barrier();
              if ( icol == sdcontext_[ispin][ikp]->mycol() ) {
                my_col_ctxt = col_ctxt;
              }

              else
                delete col_ctxt;
            }
            */


            //Don't make new column contexts from scratch, construct them from those in previous SlaterDet
            for ( int kloc=0; kloc<wf.nkptloc(); kloc++) {
               int kp = wf.kptloc(kloc);         // global index of local kpoint
               if ( wf.sd(ispin,kp) != 0 ) {
                  Context* col_ctxt = new Context(wf.sd(ispin,kp)->col_ctxt());
                  sd_[ispin][kp] = new SlaterDet(*sdcontext_[ispin][ikp],*col_ctxt,
                                                 *sdcontextsq_[ispin][ikp],kpoint_[kp],
                                                 ultrasoft_,force_complex_wf_);
                  mysdctxt_[kp] = ikp;
                  if (wf.sd(ispin,kp)->highmem())
                     sd_[ispin][kp]->set_highmem();

               }
            }
          }
        }
      }
    }
  }  

  kptproc0_.resize(nspin_);
  nkplocproc0_.resize(nspin_);
  for (int ispin=0; ispin<nspin_; ispin++) {
    kptproc0_[ispin].resize(nkppar);
    nkplocproc0_[ispin].resize(nkppar);
  }
  for (int ispin=0; ispin<nspin_; ispin++) {
    for (int kloc=0; kloc<nkppar; kloc++) {
      kptproc0_[ispin][kloc] = wf.kptproc0(ispin,kloc);
      nkplocproc0_[ispin][kloc] = wf.nkplocproc0(ispin,kloc);
    }
  }
  resize(cell_,refcell_,ecut_);
  reset();

  hasdata_ = true;   // wf has been allocated

}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::~Wavefunction() {
  deallocate();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print_timing(void) {
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ ) {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    uint64_t count = (*i).second.counts();
    //if ( ctxt_.myproc()==0 ) {
    if ( ctxt_.mype()==0 ) {
       cout << left << setw(34) << "<timing where=\"wavefunction\""
            << setw(8) << " name=\""
            << setw(15) << (*i).first << "\""
            << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
            << " max=\"" << setprecision(3) << setw(9) << tmax << "\""
            << " count=\"" << setw(9) << count << "\"/>"
            << endl;
    }
  }

  if ( spincontext_[0] != 0 ) 
     if (sdcontext_[0][0] != 0 ) 
        if (sdcontext_[0][0]->active() ) 
           sd_[0][0]->print_timing();

}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::allocate(void) {

   tmap["allocate"].start();

   // create sd contexts and SlaterDets
  const int nkpts = nkp();
  assert(ctxt_.active());
  
  spincontext_.resize(nspin_);
  sdcontext_.resize(nspin_);
  sdcontextsq_.resize(nspin_);
  sd_.resize(nspin_);

  // determine dimensions of wfcontext
  assert(nrowmax_>0);
  const int size = ctxt_.size();
  int npr = nrowmax_;
  if (npr > size) npr = size;
  while ( size%npr != 0 ) npr--;
  // npr now divides size
  int npc = size/npr;

  wfcontext_ = new Context(ctxt_.comm(),npr,npc);
  int npcs = npc/nspin_;
  int scol1 = npcs;
  if (npcs == 0) npcs = 1;

  // determine dimensions of sdcontext
  int nkppar = nparallelkpts_;
  if (nparallelkpts_ == 0 || nkppar > nkpts) nkppar = nkpts; 
  while (npcs%nkppar != 0) nkppar--;
  assert(nkppar > 0 && nkppar <= nkpts);
  int npck = npcs/nkppar;

  for (int ispin=0; ispin<nspin_; ispin++) {
    sdcontext_[ispin].resize(nkppar);
    sdcontextsq_[ispin].resize(nkppar);
    sd_[ispin].resize(nkpts);
  }
  for (int ispin=0; ispin<nspin_; ispin++) {
    spincontext_[ispin] = 0;
    for (int k=0; k<nkppar; k++)
    {
      sdcontext_[ispin][k] = 0;
      sdcontextsq_[ispin][k] = 0;
    }
    for (int k=0; k<nkpts; k++)
      sd_[ispin][k] = 0;
  }
  kptproc0_.resize(nspin_);
  nkplocproc0_.resize(nspin_);
  for (int ispin=0; ispin<nspin_; ispin++) {
    kptproc0_[ispin].resize(nkppar);
    nkplocproc0_[ispin].resize(nkppar);
  }
  spincontext_[0] = new Context(*wfcontext_,npr,npcs,0,0);
  if (nspin_ == 2) 
    spincontext_[1] = new Context(*wfcontext_,npr,npcs,0,scol1);

  for (int ispin=0; ispin<nspin_; ispin++) {
    if (spinactive(ispin)) {
      spinloc_ = ispin;
      spincontext_[ispin]->set_coutpe(0);
    }
  }
  mysdctxt_.resize(nkpts);
  for (int k=0; k<nkpts; k++)
    mysdctxt_[k] = -1;
  
  // create list of local k-points for each sdcontext to be used to create
  // SlaterDet objects
  int nkpleftover = nkpts%nkppar;
  nkptloc_ = nkpts/nkppar;
  kptloc_.resize(nkptloc_);

  // store number of local k-points on each sdcontext on all procs
  nkptloc_list_.resize(nkppar);
  for (int k=0; k<nkppar; k++)
    nkptloc_list_[k] = nkptloc_;
  for (int k=0; k<nkpleftover; k++)
    nkptloc_list_[k]++;
  
  for (int ispin=0; ispin<nspin_; ispin++) {
    if (nkpts == 1) {
      if ( spinactive(ispin) ) {
         //if ( ctxt_.oncoutpe() )
         if ( spincontext_[ispin]->myproc() == 0)
            cout << "<!-- Creating SlaterDet context " << spincontext_[ispin]->nprow() << "x" << spincontext_[ispin]->npcol() << " from spincontext, ispin = " << ispin << " -->" << endl;
        sdcontext_[ispin][0] = new Context(*spincontext_[ispin],spincontext_[ispin]->nprow(),spincontext_[ispin]->npcol(),0,0);

        //ewd:  right now, we're assuming default block size, i.e. square matrices live only on npcol x npcol tasks.  May want to generalize by computing subset of tasks square matrices live on from block size.

        sdcontextsq_[ispin][0] = new Context(*sdcontext_[ispin][0],sdcontext_[ispin][0]->npcol(),sdcontext_[ispin][0]->npcol(),0,0);
        
        Context* my_col_ctxt = 0;
        for ( int icol = 0; icol < sdcontext_[ispin][0]->npcol(); icol++ ) {
          Context* col_ctxt = new Context(*sdcontext_[ispin][0],sdcontext_[ispin][0]->nprow(),1,0,icol);
          sdcontext_[ispin][0]->barrier();
          if ( icol == sdcontext_[ispin][0]->mycol() )
            my_col_ctxt = col_ctxt;
          else
            delete col_ctxt;
        }

        sd_[ispin][0] = new SlaterDet(*sdcontext_[ispin][0],*my_col_ctxt,*sdcontextsq_[ispin][0],
                                      kpoint_[0],ultrasoft_,force_complex_wf_);
        mysdctxt_[0] = 0;
      }
    }
    else {
      // divide each spincontext into subcontexts over kpoints
      if ( spinactive(ispin) ) {
        int kpcnt = 0;
        for (int k=0; k<nkppar; k++) {
          int kcol0 = k*npck;
          int tnpck = (npck > 0) ? npck : 1; 
          if ( ctxt_.oncoutpe() )
            cout << "<!-- Creating subcontext " << npr << "x" << tnpck << " at row 0, col " << kcol0 << " of spincontext " << npr << "x" << npc << " -->" << endl;
          Context* subctxt_ = new Context(*spincontext_[ispin],npr,tnpck,0,kcol0);
          if (subctxt_->active()) {
            sdcontext_[ispin][k] = subctxt_;
            sdcontextsq_[ispin][k] = new Context(*sdcontext_[ispin][k],sdcontext_[ispin][k]->npcol(),sdcontext_[ispin][k]->npcol(),0,0);
            //sdcontextsq_[ispin][k] = new Context(*spincontext_[ispin],tnpck,tnpck,0,kcol0);
            // index of first kpoint for this context
            int kp0 = k*nkptloc_ + (nkpleftover > k ? k : nkpleftover);

            // adjust nkptloc_ when kpoints not equally distributed, i.e. nkpleftover > 0
            if (nkpleftover > k)
              nkptloc_++;

            Context* my_col_ctxt = 0;
            for ( int icol = 0; icol < sdcontext_[ispin][k]->npcol(); icol++ ) {
              Context* col_ctxt = new Context(*sdcontext_[ispin][k],sdcontext_[ispin][k]->nprow(),1,0,icol);
              sdcontext_[ispin][k]->barrier();
              if ( icol == sdcontext_[ispin][k]->mycol() )
                my_col_ctxt = col_ctxt;
              else
                delete col_ctxt;
            }

            // fill kptloc array, create SlaterDet objects for each local kpoint
            int localcnt = 0;
            kptloc_.resize(nkptloc_);
            for (int kp=kp0; kp<kp0+nkptloc_; kp++) {
              if (sdcontext_[ispin][k]->myproc() == 0) {
                int tmpcoutpe = sdcontext_[ispin][k]->mype();
                sdcontext_[ispin][k]->set_coutpe(tmpcoutpe);
              }              
              if (sdcontext_[ispin][k]->oncoutpe())
                cout << "<!-- creating SlaterDet for kp = " << kp << " -->" << endl;

              sd_[ispin][kp] = new SlaterDet(*sdcontext_[ispin][k],*my_col_ctxt,*sdcontextsq_[ispin][k],
                                             kpoint_[kp],ultrasoft_,force_complex_wf_);
              mysdctxt_[kp] = k;

              kptloc_[localcnt++] = kp;
            }
          }
          else {
            delete subctxt_;
          }
        }
      }
    }
  }

  for (int ispin=0; ispin<nspin_; ispin++) {
    int tnkploc_;
    int spe0 = ispin*npr*npcs;
    if (nspin_ == 2)
      if (spinactive(0) && spinactive(1))
        spe0 = 0;
    kptproc0_[ispin][0] = spe0;
    nkplocproc0_[ispin][0] = nkptloc_;
    for (int k=1; k<nkppar; k++) {
      int kp0 = k*npck*npr;
      kptproc0_[ispin][k] = spe0+kp0;
#if USE_MPI
      MPI_Barrier(ctxt_.comm());
      if ( ctxt_.oncoutpe() ) {
        MPI_Status status;
        MPI_Recv(&tnkploc_,1,MPI_INT,kp0,kp0,ctxt_.comm(),&status);
        nkplocproc0_[ispin][k] = tnkploc_;
      }
      else if ( ctxt_.mype() == kp0 ) {
        // send nkploc_ to pe0
        tnkploc_ = nkptloc_;
        MPI_Send(&tnkploc_,1,MPI_INT,ctxt_.coutpe(),kp0,ctxt_.comm());
      }
#endif
    }
  }

  // once we've allocated, assume wavefunction has data from here on out
  if (!hasdata_)
     hasdata_ = true;

  resize(cell_,refcell_,ecut_);
  reset();

  //ewd DEBUG:  do we need this?
  update_occ(0.0,0);
  if ( ctxt_.oncoutpe() )
     cout << "<!-- Updated occupation of wf -->" << endl;
  //if ( ctxt_.oncoutpe() )
    //cout << "<!-- Occupation NOT updated during wf allocate, testing... -->" << endl;

  tmap["allocate"].stop();

}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_hasdata(bool hasd) {
  if (!hasdata_ && hasd)
    allocate();
  hasdata_ = hasd;
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::deallocate(void) {
   hasdata_ = false;
   for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if ( spincontext_[ispin] != 0 ) {
         for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
            if (sdcontext_[ispin][ikp] != 0 ) {
               if (sdcontext_[ispin][ikp]->active() ) {
                  for ( int kloc=0; kloc<nkptloc_; kloc++) {
                     int kp = kptloc_[kloc];
                     if ( sd_[ispin][kp] != 0 )
                        delete sd_[ispin][kp];
                  }
               }
               delete sdcontext_[ispin][ikp];
               delete sdcontextsq_[ispin][ikp];
            }
         }
         delete spincontext_[ispin];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::clear(void) {
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
      if (sdcontext_[ispin][ikp] != 0 ) {
        if (sdcontext_[ispin][ikp]->active() ) {
          for ( int kloc=0; kloc<nkptloc_; kloc++) {
            int kp = kptloc_[kloc];
            if ( sd_[ispin][kp] != 0 )
              sd(ispin,kp)->c().clear();
          }
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_ultrasoft(bool us) {
  ultrasoft_ = us;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_local_block(int mb, int nb) {
   mbset_ = mb;
   nbset_ = nb;
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nblocks(int mblks, int nblks) {
   mblks_ = mblks;
   nblks_ = nblks;
}
////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nkp(void) const { return kpoint_.size(); } 

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nel() const { return nel_; } // total number of electrons

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nst() const { 
  if ( nspin_ == 1 )
    return nst_[0];
  else
    return nst_[0]+nst_[1];
}

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nst(int ispin) const {
  assert(ispin >= 0 && ispin < 2);
  return nst_[ispin];
}

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nempty() const { return nempty_; } // number of empty states

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nspin() const { return nspin_; } // number of spin channels

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::deltacharge() const { return deltacharge_; }

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::mu() const { return mu_; }

////////////////////////////////////////////////////////////////////////////////
bool Wavefunction::kptactive(int k) const {
  // loop through local k-points, see if any have absolute index k
  bool match = false;
  for ( int kloc=0; kloc<nkptloc_; kloc++) {
    int kp = kptloc_[kloc];
    if (kp == k)
      match = true;
  }
  return match;
}
////////////////////////////////////////////////////////////////////////////////
bool Wavefunction::spinactive(int ispin) const {
  bool active = false;
  if (spincontext_[ispin] != 0)
    if (spincontext_[ispin]->active())
      active = true;
  return active;
}
////////////////////////////////////////////////////////////////////////////////
double Wavefunction::entropy(void) const {
  // first value is local contribution to weighted average, second is local contrib to norm
  double entropysum[2] = {0.0, 0.0}; 
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          entropysum[0] += weight_[ikp]*(sd(ispin,ikp)->entropy(nspin_));
          entropysum[1] += weight_[ikp];
        }
      }
    }  
  }  
  // average entropy sum over k-points and spin
  wfcontext()->dsum('r',2,1,&entropysum[0],2);
  assert(entropysum[1] != 0.0);
  entropysum[0] /= entropysum[1];
  return entropysum[0];
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::resize(const UnitCell& cell, const UnitCell& refcell, 
  double ecut) {

   try {
    // resize all SlaterDets using cell, refcell, ecut and nst_[ispin]
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
          if (kptactive(ikp)) {
            assert(sd_[ispin][ikp] != 0);
            sd_[ispin][ikp]->set_local_block(mbset_,nbset_);
            sd_[ispin][ikp]->set_nblocks(mblks_,nblks_);
            sd_[ispin][ikp]->resize(cell,refcell,ecut,nst_[ispin]);
          }
        }
      }
    }
    cell_ = cell;
    refcell_ = refcell;
    ecut_ = ecut;
  }
  catch ( const SlaterDetException& sdex ) {
    cout << "<ERROR> Wavefunction: SlaterDetException during resize: </ERROR>" << endl
         << sdex.msg << endl;
    // no resize took place
    return;
  }
  catch ( bad_alloc )
  {
    cout << "<ERROR> Wavefunction: insufficient memory for resize operation </ERROR>" << endl;
    return;
  }
  
}
  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::reset(void) {
  // reset all SlaterDets
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->reset();
        }
      }
    }
  }
}  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_highmem(void) {
  // toggle highmem flag on SlaterDets
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->set_highmem();
        }
      }
    }
  }
}  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::compute_nst(void) {
  // recompute nst from nel_, deltaspin_, nempty_

  nst_.resize(nspin_);
  if ( nspin_ == 1 )
  {
    nst_[0] = ( nel_ + 1 ) / 2 + nempty_;
  }
  else
  {
    // nspin == 2
    nst_[0] = ( nel_ + 1 ) / 2 + deltaspin_ + nempty_;
    nst_[1] = nel_ / 2 - deltaspin_ + nempty_;
  }
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nel(int nel) {
  if ( nel == nel_ ) return;
  if ( nel < 0 )
  {
    cout << "<ERROR> Wavefunction::set_nel: nel < 0 </ERROR>" << endl;
    return;
  }
  
  nel_ = nel;
  compute_nst();
  if (hasdata_)
    resize(cell_,refcell_,ecut_);
}
  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nempty(int nempty) {
  if ( nempty == nempty_ ) return;
  if ( nempty < 0 )
  {
    cout << "<ERROR> Wavefunction::set_nempty: negative value </ERROR>" << endl;
    return;
  }
  nempty_ = nempty;
  compute_nst();
  if (hasdata_)
    resize(cell_,refcell_,ecut_);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nspin(int nspin) {
  assert(nspin==1 || nspin==2);
  if ( nspin == nspin_ ) return;
  
  if (hasdata_) {
    deallocate();
    cout << "<!-- Wavefunction::set_nspin: " << nspin << " deallocate done -->" << endl;
  }

  nspin_ = nspin;
  compute_nst();

  if (hasdata_) {
    allocate();
    cout << "<!-- Wavefunction::set_nspin: " << nspin << " allocate done -->" << endl;
    resize(cell_,refcell_,ecut_);
    reset();
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_deltaspin(int deltaspin) {
  deltaspin_ = deltaspin;
  
  // force spin polarized calculation, even for deltaspin == 0
  if (nspin_ == 1) {
    set_nspin(2);
  }
  // spin polarization changed -- we don't need to reallocate, just update nst and SlaterDets
  else {
    compute_nst();
    if (hasdata_) {
      resize(cell_,refcell_,ecut_);
      reset();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_deltacharge(double deltacharge)
{
   if ( deltacharge == deltacharge_ ) return;
   deltacharge_ = deltacharge;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nrowmax(int n) {
  assert(n>0);
  if ( n > ctxt_.size() ) {
    cout << "<WARNING> Wavefunction::set_nrowmax: nrowmax > ctxt_.size(). </WARNING>" << endl;
    n = ctxt_.size();
  }

  // set nrowmax to be compatible with torus dimensions on BG/Q
#ifdef USE_CTF
#ifdef BGQ
  Personality_t pers;
  Kernel_GetPersonality(&pers, sizeof(pers));
  const int nTDim = 5;
  vector<int> torusdim(nTDim);
  torusdim[0] = pers.Network_Config.Anodes;
  torusdim[1] = pers.Network_Config.Bnodes;
  torusdim[2] = pers.Network_Config.Cnodes;
  torusdim[3] = pers.Network_Config.Dnodes;
  torusdim[4] = pers.Network_Config.Enodes;

  bool torusMult = false;
  int ncol = ctxt_.size() / n;
  for (int ii=0; ii<nTDim; ii++)
     if (ncol == torusdim[ii])
        torusMult = true;
  for (int ii=0; ii<nTDim; ii++)
     for (int jj=ii+1; jj<nTDim; jj++)
        if (ncol == torusdim[ii]*torusdim[jj])
           torusMult = true;

  if ( ctxt_.oncoutpe() )
  {
     if (torusMult)
        cout << "Wavefunction::set_nrowmax:  nrowmax = " << n << " is compatible with BG/Q torus " <<
            torusdim[0] << " x " << torusdim[1] << " x " << torusdim[2] << " x " << torusdim[3] << " x " <<
            torusdim[4] << endl;
     else
        cout << "<WARNING> Wavefunction::set_nrowmax:  nrowmax = " << n << " is NOT compatible with BG/Q torus! " <<
            torusdim[0] << " x " << torusdim[1] << " x " << torusdim[2] << " x " << torusdim[3] << " x " <<
            torusdim[4] << " </WARNING> " << endl;
  }        
#endif
#endif
  
  if (nrowmax_ != n) {
     nrowmax_ = n;
     
    if (hasdata_)
       deallocate(); 
    //reshape();  // single col_ctxt not being reshaped correctly, just deallocate
  }
}
  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nparallelkpts(int n) {
  if (nparallelkpts_ != n) {
    nparallelkpts_ = n;
    if (hasdata_) 
       deallocate();  // single col_ctxt not being reshaped correctly, just deallocate
       //reshape();
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::reshape(void) {

  // create sd contexts and SlaterDets
  const int nkpts = nkp();
  assert(ctxt_.active());
  spincontext_.resize(nspin_);
  sdcontext_.resize(nspin_);
  sdcontextsq_.resize(nspin_);
  sd_.resize(nspin_);

  // determine dimensions of spincontext
  assert(nrowmax_>0);
  const int size = ctxt_.size();
  int npr = nrowmax_;
  while ( size%npr != 0 ) npr--;
  // npr now divides size
  int npc = size/npr;

  if (wfcontext_ != 0) delete wfcontext_;
  wfcontext_ = new Context(ctxt_.comm(),npr,npc);
  int npcs = npc/nspin_;
  int scol1 = npcs;
  if (npcs == 0) npcs = 1;

  // determine dimensions of sdcontext
  int nkppar = nparallelkpts_;
  if (nparallelkpts_ == 0 || nkppar > nkpts) nkppar = nkpts; 
  while (npcs%nkppar != 0) nkppar--;
  assert(nkppar > 0 && nkppar <= nkpts);
  int npck = npcs/nkppar;
  
  bool spin_changed = false;
  if (spincontext_[0] != 0)
    if (spincontext_[0]->size() != npr*npcs)
      spin_changed = true;

  for (int ispin=0; ispin<nspin_; ispin++) {
    sdcontext_[ispin].resize(nkppar);
    sdcontextsq_[ispin].resize(nkppar);
    sd_[ispin].resize(nkpts);
  }
  for (int ispin=0; ispin<nspin_; ispin++) 
    if (spincontext_[ispin] != 0)
      delete spincontext_[ispin];

  kptproc0_.resize(nspin_);
  nkplocproc0_.resize(nspin_);
  for (int ispin=0; ispin<nspin_; ispin++) {
    kptproc0_[ispin].resize(nkppar);
    nkplocproc0_[ispin].resize(nkppar);
  }

  spincontext_[0] = new Context(*wfcontext_,npr,npcs,0,0);
  if (nspin_ == 2) 
    spincontext_[1] = new Context(*wfcontext_,npr,npcs,0,scol1);

  for (int ispin=0; ispin<nspin_; ispin++) {
    if (spinactive(ispin)) {
      spinloc_ = ispin;
      spincontext_[ispin]->set_coutpe(0);
    }
  }
  
  // create list of local k-points for each sdcontext to be used to create
  // SlaterDet objects
  int nkpleftover = nkpts%nkppar;
  nkptloc_ = nkpts/nkppar;
  kptloc_.resize(nkptloc_);

  // store number of local k-points on each sdcontext on all procs
  nkptloc_list_.resize(nkppar);
  for (int k=0; k<nkppar; k++)
    nkptloc_list_[k] = nkptloc_;
  for (int k=0; k<nkpleftover; k++)
    nkptloc_list_[k]++;

  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if ( spinactive(ispin) && !spin_changed) {
      if (nkp() == 1) {
        Context* oldctxt = sdcontext_[ispin][0];
        Context* newctxt = new Context(*spincontext_[ispin],spincontext_[ispin]->nprow(),spincontext_[ispin]->npcol());
        Context* newctxtsq = new Context(*newctxt,newctxt->npcol(),newctxt->npcol(),0,0);

        if (oldctxt->size() == newctxt->size()) {
          Context* my_col_ctxt = 0;
          for ( int icol = 0; icol < newctxt->npcol(); icol++ ) {
            Context* col_ctxt = new Context(*newctxt,newctxt->nprow(),1,0,icol);
            newctxt->barrier();
            if ( icol == newctxt->mycol() )
              my_col_ctxt = col_ctxt;
            else
              delete col_ctxt;
          }
          sd_[ispin][0]->reshape(*newctxt,*my_col_ctxt,*newctxtsq,true);
        }
      }
      else {
        int kpcnt = 0;
        for (int k=0; k<nkppar; k++) {
          int kcol0 = k*npck;
          if ( ctxt_.oncoutpe() )
            cout << "<!-- Creating subcontext " << npr << "x" << npck << " at row 0, col " <<
              kcol0 << " of spincontext " << npr << "x" << npc << " -->" << endl;
          Context* subctxt_ = new Context(*spincontext_[ispin],npr,npck,0,kcol0);
          Context* subctxtsq_ = new Context(*subctxt_,subctxt_->npcol(),subctxt_->npcol(),0,0);
          if (subctxt_->active()) {
            Context* my_col_ctxt = 0;
            for ( int icol = 0; icol < subctxt_->npcol(); icol++ ) {
              Context* col_ctxt = new Context(*subctxt_,subctxt_->nprow(),1,0,icol); 
              if ( icol == subctxt_->mycol() ) 
                my_col_ctxt = col_ctxt; 
              else 
                delete col_ctxt; 
            } 
            Context* oldctxt = 0;
            if (nkptloc_ > 1)
              oldctxt = sdcontext_[ispin][k];
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int ikp = kptloc_[kloc];
              bool setnewctxt = (kloc == nkptloc_-1) ? true : false;
              sd_[ispin][ikp]->reshape(*subctxt_,*my_col_ctxt,*subctxtsq_,setnewctxt);
            }
          }
          else {
            sdcontext_[ispin][k] = 0;
            delete subctxt_;
          }
        }
      }
    }
    else if (spin_changed) {
      if (nkp() == 1) {
        SlaterDet* tmpsd = 0;
        Context* newctxt = 0;
        Context* newctxtsq = 0;
        if (spinactive(ispin)) {
          newctxt = new Context(*spincontext_[ispin],spincontext_[ispin]->nprow(),spincontext_[ispin]->npcol());
          newctxtsq = new Context(*newctxt,newctxt->npcol(),newctxt->npcol(),0,0);
          
          Context* my_col_ctxt = 0;
          for ( int icol = 0; icol < newctxt->npcol(); icol++ ) {
            Context* col_ctxt = new Context(*newctxt,newctxt->nprow(),1,0,icol);
            newctxt->barrier();
            if ( icol == newctxt->mycol() )
              my_col_ctxt = col_ctxt;
            else
              delete col_ctxt;
          }
          tmpsd = new SlaterDet(*newctxt,*my_col_ctxt,*newctxtsq,kpoint_[0],ultrasoft_,force_complex_wf_);
          tmpsd->resize(cell_,refcell_,ecut_,nst_[ispin]);
        }
        sd_[ispin][0]->copyTo(tmpsd);

        delete sd_[ispin][0];
        delete sdcontext_[ispin][0];
        delete sdcontextsq_[ispin][0];
        sdcontext_[ispin][0] = newctxt;
        sdcontextsq_[ispin][0] = newctxtsq;
        sd_[ispin][0] = tmpsd;
      }
      else {
        int kpcnt = 0;
        for (int k=0; k<nkppar; k++) {
          int kcol0 = k*npck;
          if ( ctxt_.oncoutpe() )
            cout << "<!-- Creating subcontext " << npr << "x" << npck << " at row 0, col " <<
              kcol0 << " of spincontext " << npr << "x" << npc << " -->" << endl;


          for ( int kloc=0; kloc<nkptloc_; kloc++) {
            int ikp = kptloc_[kloc];

            SlaterDet* tmpsd = 0;
            Context* subctxt_ = 0;
            Context* subctxtsq_ = 0;

            if (spinactive(ispin)) {
              subctxt_ = new Context(*spincontext_[ispin],npr,npck,0,kcol0);
              subctxtsq_ = new Context(*subctxt_,subctxt_->npcol(),subctxt_->npcol(),0,0);
              if (subctxt_->active()) {
                Context* my_col_ctxt = 0;
                for ( int icol = 0; icol < subctxt_->npcol(); icol++ ) {
                  Context* col_ctxt = new Context(*subctxt_,subctxt_->nprow(),1,0,icol); 
                  if ( icol == subctxt_->mycol() ) 
                    my_col_ctxt = col_ctxt; 
                  else 
                    delete col_ctxt; 
                }

                tmpsd = new SlaterDet(*subctxt_,*my_col_ctxt,*subctxtsq_,kpoint_[ikp],
                                      ultrasoft_,force_complex_wf_);
              }
            }
            sd_[ispin][ikp]->copyTo(tmpsd);

            delete sd_[ispin][ikp];
            delete sdcontext_[ispin][ikp];
            delete sdcontextsq_[ispin][ikp];
            sdcontext_[ispin][ikp] = subctxt_;
            sdcontextsq_[ispin][ikp] = subctxtsq_;
            sd_[ispin][ikp] = tmpsd;
          }
        }
      }
    }
    
  }
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::add_kpoint(D3vector kpoint, double weight) {

  // add_kpoint adds new k-point to kpoint_ vector, except on first call when
  // default gamma point is replaced

  if ( ctxt_.oncoutpe() )
    cout << "Adding kpoint = " << kpoint << ", weight = " << weight << endl;

  if (kpt_added_) {
    for ( int i = 0; i < kpoint_.size(); i++ ) {
      if ( kpoint == kpoint_[i] ) {
        if ( ctxt_.oncoutpe() )
          cout << "<WARNING> Wavefunction::add_kpoint: warning: kpoint already defined </WARNING>" 
            << endl;
      }
    }
  }
  
  if (hasdata_)
    deallocate();
  
  if (kpt_added_) {
    kpoint_.push_back(kpoint);
    weight_.push_back(weight);
    weightsum_ = 0.0;
    for ( int i = 0; i < weight_.size(); i++ )
      weightsum_ += weight_[i];
  }
  else {
    kpt_added_ = true;
    assert(kpoint_.size()==1);
    assert(weight_.size()==1);
    kpoint_[0] = kpoint;
    weight_[0] = weightsum_ = weight;
  }
  assert(weightsum_ > 0.0);

  if (hasdata_) {
    if ( ctxt_.oncoutpe() )
      cout << "<!-- Wavefunction::add_kpoint: " << nkp() << " k-points defined, calling allocate -->" << endl;
    allocate();
    if ( ctxt_.oncoutpe() )
      cout << " Wavefunction::add_kpoint: " << kpoint << " allocate done" << endl;
    resize(cell_,refcell_,ecut_);
    reset();

    /* ewd:  I don't think this is necessary
    update_occ(0.0,0);
    if ( ctxt_.oncoutpe() )
      cout << "Updated occupation of wf" << endl;
    */

  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::del_kpoint(D3vector kpoint) {
  cout << "<ERROR> Wavefunction::del_kpoint: not implemented </ERROR>" << endl;
  assert(false);

  // remember to recalculate weightsum_ after removing kpoint
  weightsum_ = 0.0;
  for ( int i = 0; i < weight_.size(); i++ )
    weightsum_ += weight_[i];
  assert(weightsum_ > 0.0);

}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize(double amplitude, bool highmem) {
  if (!hasdata_) {
    hasdata_ = true;

    //ewd DEBUG
    if (ctxt_.mype() == 0)
      cout << "<!-- Randomize_wf:  allocating wavefunction... -->" << endl;

    allocate();
  }
  if (highmem)
    set_highmem();

  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->randomize(amplitude);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize_us(double amplitude, AtomSet& as, bool highmem) {
  if (!hasdata_) {
    hasdata_ = true;
    if (ctxt_.mype() == 0)
      cout << "<!-- Randomize_wf_us:  allocating wavefunction... -->" << endl;
    allocate();
  }

  if (highmem)
    set_highmem();

  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->randomize_us(amplitude,as);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize_real(double amplitude)
{
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->randomize_real(amplitude);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// AS: shift state n_state by the vector (shift_x, shift_y, shift_z)
void Wavefunction::shift_wf(double shift_x,double shift_y,double shift_z, int n_state)
{
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->shift_wf(shift_x,shift_y,shift_z,n_state);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// AS: change phase of the wave function to make it real for Gamma only
void Wavefunction::phase_wf_real(void)
{
  // AS: DEBUG
   if ( ctxt_.oncoutpe() )
      cout << " AS: changing phase of WF to make it real" << endl;

  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->phase_wf_real();
      sd_[ispin][ikp]->gram();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::rescale(double factor) {
  if (!hasdata_) {
    if ( ctxt_.oncoutpe() )
      cout << "<ERROR> Wavefunction.rescale called before allocate()! </ERROR>" << endl;
    assert(false);
  }
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->rescale(factor);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::promote_occ(double occ_change, int origin_level, int destination_level, int ispin)
{
   assert(ispin < nspin_);
   {
      for ( int ikp = 0; ikp < nkp(); ikp++ )
      {
         sd_[ispin][ikp]->promote_occ(occ_change, origin_level, destination_level);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::update_occ(double temp, int ngauss) {
  // update occupation numbers using eigenvalues in SlaterDet
  vector<double> totalcharge(nspin_);

  if ( temp == 0.0 ) {
    // zero temperature
    if ( nspin_ == 1 ) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          sd_[0][ikp]->update_occ(nel_,nspin_);
        }
      }
    }
    else if ( nspin_ == 2 ) {
      const int nocc_up = (nel_+1)/2+deltaspin_;
      const int nocc_dn = nel_/2 - deltaspin_;
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp) && spinactive(0)) {  
          sd_[0][ikp]->update_occ(nocc_up,nspin_);
        }
        if (kptactive(ikp) && spinactive(1)) {  
          sd_[1][ikp]->update_occ(nocc_dn,nspin_);
        }
      }
    }
    else {
      // incorrect value of nspin_
      assert(false);
    }
  }
  else {
    // finite temperature
    const double eVolt = 0.036749023; // 1 eV in Hartree
    //const int maxiter = 1000;
    const int maxiter = 10000;
 
    // loop to find value of mu
    double mu[2] = {0.0,0.0};
    //double dmu[2] = {2.0 * eVolt,2.0*eVolt};
    double dmu[2] = {10.0 * eVolt,10.0*eVolt};
    enum direction { up, down };
    direction dir[2] = {up, up};

    double rhosum[2] = {0.0,0.0};
    double rhonorm[2] = {0.0,0.0};
    double kptweight[2] = {0.0,0.0};
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
          if (kptactive(ikp)) {
            assert(sd_[ispin][ikp] != 0);
            sd_[ispin][ikp]->update_occ(nspin_,mu[ispin],temp,ngauss);
            rhosum[ispin] += weight_[ikp]*sd(ispin,ikp)->total_charge();
            rhonorm[ispin] += weight_[ikp];
          }
        }
        spincontext_[ispin]->dsum('r',1,1,&rhosum[ispin],1);
        spincontext_[ispin]->dsum('r',1,1,&rhonorm[ispin],1);
        assert(rhonorm[ispin] != 0.0);
        kptweight[ispin] = rhonorm[ispin];
      }
    }

    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {

        if (nspin_ == 1)
          totalcharge[ispin] = (double) nel_ + deltacharge_;
        else if ( nspin_ == 2 ) {
          if (ispin == 0)
             totalcharge[ispin] = (nel_+1)/2+deltaspin_ + deltacharge_;  // ewd:  need a factor of 0.5 in front of deltacharge??
          else if (ispin == 1)
            totalcharge[ispin] = nel_/2 - deltaspin_ + deltacharge_;
        }
        
        int niter = 0;
        double wtcharge = totalcharge[ispin]*kptweight[ispin];
        while ( niter < maxiter && fabs(rhosum[ispin] - wtcharge) > 1.e-8*kptweight[ispin] ) {
          niter++;
          if ( rhosum[ispin] < wtcharge ) {
            if ( dir[ispin] == down ) dmu[ispin] /= 2.0;
            mu[ispin] += dmu[ispin];
            dir[ispin] = up;
          }
          else {
            if ( dir[ispin] == up ) dmu[ispin] /= 2.0;
            mu[ispin] -= dmu[ispin];
            dir[ispin] = down;
          }

          rhosum[ispin] = 0.0;
          rhonorm[ispin] = 0.0;
      
          for ( int ikp=0; ikp<nkp(); ikp++) {
            if (kptactive(ikp)) {
              assert(sd_[ispin][ikp] != 0);
              sd_[ispin][ikp]->update_occ(nspin_,mu[ispin],temp,ngauss);
              rhosum[ispin] += weight_[ikp]*sd(ispin,ikp)->total_charge();
              rhonorm[ispin] += weight_[ikp];
            }
          }
          spincontext_[ispin]->dsum('r',1,1,&rhosum[ispin],1);
          spincontext_[ispin]->dsum('r',1,1,&rhonorm[ispin],1);
          assert(rhonorm[ispin] != 0.0);
          kptweight[ispin] = rhonorm[ispin];
        }
 
        if ( niter == maxiter ) {
          double tmp = fabs(rhosum[ispin] - wtcharge);
          cout << "<ERROR> Wavefunction::update_occ: mu did not converge in "
               << maxiter << " iterations : mu = " << mu[ispin] << ", dmu = " << dmu[ispin] << ", rhosum[ispin] = " << rhosum[ispin] << ", totalcharge = " << totalcharge[ispin] << ", wtcharge = " << wtcharge << ", diff = " << tmp << ", ispin = " << ispin << " </ERROR>" << endl;
          //ewd DEBUG DEBUG DEBUG
          //ctxt_.abort(1);
        }

        mu_ = mu[ispin];   //ewd:  should we store mu separately for each spin channel or not?

      }
    }

    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        if ( spincontext_[ispin]->myproc() == 0 ) {
          //!! print on one process only
          if (nspin_ == 1) 
            cout << " <!-- Wavefunction::update_occ: sum = "
                 << rhosum[0]/rhonorm[0] << " -->" << endl;
          else
            cout << " <!-- Wavefunction::update_occ: ispin = " << ispin << ", sum = "
                 << rhosum[ispin]/rhonorm[ispin] << " -->" << endl;


          string smearname = "Methfessel-Paxton";
          if (ngauss == -1) 
            smearname = "Fermi";
          else if (ngauss == 0) 
            smearname = "Gaussian";
      
          cout << " <!-- Wavefunction::update_occ: using " << smearname << " smearing, mu = "
               << setprecision(4) << mu[ispin] / eVolt << " eV" << ", ispin = " << ispin << " -->" << endl;
        }
      }
    }



  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::init_usfns(AtomSet* atoms) {
  if (ultrasoft_) {
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
          if (kptactive(ikp)) {
            assert(sd_[ispin][ikp] != 0);
            sd_[ispin][ikp]->init_usfns(atoms);
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::update_usfns() {
  if (ultrasoft_) {
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
          if (kptactive(ikp)) {
            assert(sd_[ispin][ikp] != 0);
            sd_[ispin][ikp]->update_usfns();
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::calc_spsi() {
  if (ultrasoft_) {
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
          if (kptactive(ikp)) {
            assert(sd_[ispin][ikp] != 0);
            sd_[ispin][ikp]->calc_spsi();
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::gram(void) {
  QB_Pstart(12,gram);
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->gram();
        }
      }
    }
  }
  QB_Pstop(gram);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::riccati(Wavefunction& wf) {
  assert(wf.context() == ctxt_);
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->riccati(*wf.sd_[ispin][ikp]);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::align(Wavefunction& wf) {
  assert(wf.context() == ctxt_);
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->align(*wf.sd_[ispin][ikp]);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::dot(const Wavefunction& wf) const {
  assert(wf.context() == ctxt_);
  double sum[2] = {0.0,0.0};
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          //sum[0] += weight_[ikp]*sd_[ispin][ikp]->dot(*wf.sd_[ispin][ikp]);

          double tmp = sd_[ispin][ikp]->dot(*wf.sd_[ispin][ikp]);
          sum[0] += weight_[ikp]*tmp;
          sum[1] += weight_[ikp];
        }
      }
    }
  }
  // dot product is an average over spin and k-points
  wfcontext()->dsum('r',2,1,&sum[0],2);
  assert(sum[1] != 0.0);
  sum[0] /= sum[1];
  return sum[0];
}

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::sdot(const Wavefunction& wf) const {
  // for ultrasoft, want to dot wf with S*wf

  assert(wf.context() == ctxt_);
  double sum[2] = {0.0,0.0};
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          //sum[0] += weight_[ikp]*sd_[ispin][ikp]->dot(*wf.sd_[ispin][ikp]);

          double tmp = sd_[ispin][ikp]->sdot(*wf.sd_[ispin][ikp]);
          sum[0] += weight_[ikp]*tmp;
          sum[1] += weight_[ikp];
        }
      }
    }
  }
  // dot product is an average over spin and k-points
  wfcontext()->dsum('r',2,1,&sum[0],2);
  assert(sum[1] != 0.0);
  sum[0] /= sum[1];
  return sum[0];
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::diag(Wavefunction& dwf, bool eigvec) {

   //bool copyToSquareContext = true;
   bool copyToSquareContext = false;

   QB_Pstart(13,diag);
  // subspace diagonalization of <*this | dwf>
  // if eigvec==true, eigenvectors are computed and stored in *this, dwf is 
  // overwritten
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];
              if ( sd_[ispin][kp] != 0 ) {

                // compute eigenvalues
                if ( sd(ispin,kp)->basis().real() ) {
                  // proxy real matrices c, cp
                  DoubleMatrix c(sd(ispin,kp)->c());
                  DoubleMatrix sc(sd(ispin,kp)->spsi());
                  DoubleMatrix cp(dwf.sd(ispin,kp)->c());
 
                  DoubleMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
                  valarray<double> w(h.m());

                  if (ultrasoft_)
                  {  
                     assert(false);  // ultrasoft needs complex data
                  }
                  else            // diagonalize psi*hpsi
                  {
                     tmap["diag-gemm1"].start();
                     // factor 2.0 in next line: G and -G
                     h.gemm('t','n',2.0,c,cp,0.0);
                     // rank-1 update correction
                     h.ger(-1.0,c,0,cp,0);
                     tmap["diag-gemm1"].stop();
                  }                  
                  
                  if (copyToSquareContext)
                  {

                     //if (sdcontext_[ispin][ikp]->oncoutpe()) 
                     //  cout << "<!-- Wavefunction::diag: using in-place data move to speed up eigensolve. -->" << endl; 
                  
                    DoubleMatrix hsq(*sdcontextsq_[ispin][ikp],c.n(),c.n(),c.nb(),c.nb());
                    if (hsq.active())
                       h.copyInPlace(hsq);

                    if ( eigvec ) {
                       DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
                       z.clear();
                       tmap["diag-syevd"].start();
                       if (hsq.active())
                       {
                          DoubleMatrix zsq(*sdcontextsq_[ispin][ikp],c.n(),c.n(),c.nb(),c.nb());
                          hsq.syevd('l',w,zsq);
                          // copy eigenvectors back to default context
                          zsq.copyInPlace(z);
                       }
                       tmap["diag-syevd"].stop();

                       // add barrier to keep timing clear
                       c.context().barrier();
                       
                       cp = c;
                       tmap["diag-gemm2"].start();
                       c.gemm('n','n',1.0,cp,z,0.0);
                       tmap["diag-gemm2"].stop();
                    }
                    else {
                       tmap["diag-syevd"].start();
                       if (hsq.active())
                          hsq.syevd('l',w);  // need to copy w to other tasks?
                       tmap["diag-syevd"].stop();
                    }

                    // all tasks need eigenvalues to calculate occupation
                    MPI_Bcast(&w[0], c.n(), MPI_DOUBLE, 0, c.context().comm());

                  }
                  else {

#if 0
                    //ewd: test jacobi eigensolver
                    vector<double> w(h.m());
                    DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
                    const int maxsweep = 30;
                    int nsweep = jacobi(maxsweep,1.e-6,h,z,w);
                    if ( eigvec )
                    {
                      cp = c;
                      c.gemm('n','n',1.0,cp,z,0.0);
                    }
#else                    
                    if ( eigvec ) {
                      DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());

                      //ewd:  for now, always use syevd.  We need to decide
                      // how to toggle between syev and syevd.
                      tmap["diag-syevd"].start();
                      h.syevd('l',w,z);
                      tmap["diag-syevd"].stop();
                      cp = c;
                      tmap["diag-gemm2"].start();
                      c.gemm('n','n',1.0,cp,z,0.0);
                      tmap["diag-gemm2"].stop();
                    }
                    else {
                      tmap["diag-syevd"].start();
                      h.syevd('l',w);
                      tmap["diag-syevd"].stop();
                      //h.syev('l',w);
                    }
#endif
                  }
                  // set eigenvalues in SlaterDet
                  sd(ispin,kp)->set_eig(w);
                }
                else {
                  ComplexMatrix& c = sd(ispin,kp)->c();
                  ComplexMatrix& cp = dwf.sd(ispin,kp)->c();
                  ComplexMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
                  tmap["diag-gemm1"].start();
                  h.gemm('c','n',1.0,c,cp,0.0);
                  tmap["diag-gemm1"].stop();
                  valarray<double> w(h.m());

                  if (copyToSquareContext)
                  {
                     //if (sdcontext_[ispin][ikp]->oncoutpe()) 
                     //  cout << "<!-- Wavefunction::diag: using in-place data move to speed up eigensolve. -->" << endl; 
                    ComplexMatrix hsq(*sdcontextsq_[ispin][ikp],c.n(),c.n(),c.nb(),c.nb());
                    if (hsq.active())
                       h.copyInPlace(hsq);

                    if ( eigvec ) {
                       ComplexMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
                       z.clear();
                       tmap["diag-heevd"].start();
                       if (hsq.active())
                       {
                          ComplexMatrix zsq(*sdcontextsq_[ispin][ikp],c.n(),c.n(),c.nb(),c.nb());
                          hsq.heevd('l',w,zsq);
                          // copy eigenvectors back to default context
                          zsq.copyInPlace(z);
                       }
                       tmap["diag-heevd"].stop();

                       // add barrier to keep timing clear
                       c.context().barrier();
                       
                       cp = c;
                       tmap["diag-gemm2"].start();
                       c.gemm('n','n',1.0,cp,z,0.0);
                       tmap["diag-gemm2"].stop();
                    }
                    else {
                       tmap["diag-heevd"].start();
                       if (hsq.active())
                       {
                          hsq.heevd('l',w);
                       }
                       tmap["diag-heevd"].start();
                    }

                    // all tasks need eigenvalues to calculate occupation
                    MPI_Bcast(&w[0], c.n(), MPI_DOUBLE, 0, c.context().comm());

                  }
                  else {   // no context reshaping
                    if ( eigvec ) {
                      ComplexMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
                      //h.heevx('l',w,z);
                      tmap["diag-heevd"].start();
                      h.heevd('l',w,z);
                      tmap["diag-heevd"].stop();
                      //h.heev('l',w,z);
                        
                      cp = c;
                      tmap["diag-gemm2"].start();
                      c.gemm('n','n',1.0,cp,z,0.0);
                      tmap["diag-gemm2"].stop();
                      
                    }
                    else {
                       tmap["diag-heev"].start();
                       h.heev('l',w);
                       tmap["diag-heev"].stop();
                    }
                  }
                  // set eigenvalues in SlaterDet
                  sd(ispin,kp)->set_eig(w);      // need to copy w to !ctxtsq_.active() tasks?
                }
              }
            }
          }
        }
      }
    }
  }
  QB_Pstop(diag);
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::extrap_real(const double dt, const AtomSet& as) {
  // transform states to real space and move them according to distance
  // from neighboring atoms using weighting function. 
  // (algorithm idea suggested by Miguel Morales, adapted from Umrigar's QMC 
  // "space warp" technique)

  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];         // global index of local kpoint
              if ( sd_[ispin][kp] != 0 ) {
                int nstloc = sd_[ispin][kp]->nstloc();
                int nb = sd_[ispin][kp]->c().nb();
                
                const Basis& basis = sd_[ispin][kp]->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
                const int myproc = basis.context().myproc();
                const int np0 = ft.np0();
                const int np1 = ft.np1();
                const int np2 = ft.np2();
                const int np2first = ft.np2_first(myproc);
                const int np012loc = ft.np012loc();
                const int np012 = ft.np012();

                for ( int n = 0; n < nstloc; n++ ) {
                  int mloc = sd_[ispin][kp]->c().mloc();
                  vector<complex<double> > wftmp(np012loc);
                  vector<complex<double> > wftmp2(np012loc);
                  vector<complex<double> > twfr(np012);

                  ComplexMatrix& c = sd_[ispin][kp]->c();
                  ft.backward(c.cvalptr(mloc*n),&wftmp[0]);

                  // copy full state to all procs in column
                  for (int i=0; i<np012; i++)
                    twfr[i] = 0.0;
                  int offset = np0*np1*np2first;
                  for (int i=0; i<np012loc; i++)
                    twfr[i+offset] = wftmp[i];

                  double *a = (double*) &twfr[0];
                  int size = 2*np012;
                  basis.context().dsum('c',size,1,&a[0],size);

                  // loop over local grid points, atoms in AtomSet to 
                  // calculate distance weighting fn
                  vector<vector<double> > rion;
                  vector<vector<double> > vion;
                  rion.resize(as.nsp());
                  vion.resize(as.nsp());
                  for ( int is = 0; is < rion.size(); is++ ) {
                    rion[is].resize(3*as.na(is));
                    vion[is].resize(3*as.na(is));
                  }
                  as.get_positions(rion,true);
                  as.get_velocities(vion,true);
                
                  vector<vector<vector<double> > > exwt;
                  exwt.resize(as.nsp());
                  for ( int is = 0; is < exwt.size(); is++ ) {
                    exwt[is].resize(as.na(is));
                    for ( int ia = 0; ia < exwt[is].size(); ia++ )
                      exwt[is][ia].resize(np012loc);
                  }
                  vector<double> exnorm;
                  exnorm.resize(np012loc);
                  for (int i=0; i<np012loc; i++)
                    exnorm[i] = 0.0;
                  
                  int ip0 = 0;
                  int ip1 = 0;
                  int ip2 = np2first;
                  D3vector a0 = cell_.a(0);
                  D3vector a1 = cell_.a(1);
                  D3vector a2 = cell_.a(2);
                  double a0len = length(a0);
                  double a1len = length(a1);
                  double a2len = length(a2);
                  double a0delta = a0len/(double)np0;
                  double a1delta = a1len/(double)np1;
                  double a2delta = a2len/(double)np2;
                  for (int i = 0; i < np012loc; i++) {
                    D3vector re = ((double)ip0/(double)np0)*a0 + ((double)ip1/(double)np1)*a1 + ((double)ip2/(double)np2)*a2;
                    for ( int is = 0; is < as.nsp(); is++ ) {
                      for ( int ia = 0; ia < as.na(is); ia++ ) {
                        // use (r-R)^-4 as distance weighting function
                        D3vector ri(rion[is][3*ia],rion[is][3*ia+1],rion[is][3*ia+2]);
                        D3vector rie = re - ri;
                        cell_.fold_in_ws(rie);
                        const double val = norm(rie);
                        double twt = 1./(val*val);
                        twt = min(twt,1.E+20);
                        exwt[is][ia][i] = twt;
                        exnorm[i] += twt;
                      }
                    }
                    ip0++;
                    if (ip0 >= np0) {
                      ip0 = 0;
                      ip1++;
                      if (ip1 >= np1) { 
                        ip1 = 0;
                        ip2++;
                      }
                    }
                  }

                  // now use weight functions to calculate new grid points for 
                  // real-space wf values
                  ip0 = 0;
                  ip1 = 0;
                  ip2 = np2first;
                  for (int i = 0; i < np012loc; i++) {
                    D3vector re = ((double)ip0/(double)np0)*a0 + ((double)ip1/(double)np1)*a1 + ((double)ip2/(double)np2)*a2;
                    D3vector renew = re;
                    for ( int is = 0; is < as.nsp(); is++ ) {
                      for ( int ia = 0; ia < as.na(is); ia++ ) {
                        // atomwt should be between 0 and dt
                        double atomwt = dt*exwt[is][ia][i]/exnorm[i]; 
                        D3vector vi(vion[is][3*ia],vion[is][3*ia+1],vion[is][3*ia+2]);
                        // wf at re moves to re+dt*v, move value 
                        // at -dt*v to re
                        renew = renew - atomwt*vi;
                      }
                    }
                    int jp0 = ip0;
                    int jp1 = ip1;
                    int jp2 = ip2;
                    int jp0b,jp1b,jp2b;
                    jp0b = jp0+1;
                    if (jp0b >= np0) jp0b -= np0;
                    jp1b = jp1+1;
                    if (jp1b >= np1) jp1b -= np1;
                    jp2b = jp2+1;
                    if (jp2b >= np2) jp2b -= np2;
                    D3vector re_delta = renew - re;
                    // we assume re_delta < grid spacing
                    if (re_delta.x < 0.0) {
                      jp0--;
                      if (jp0 < 0) jp0 += np0;
                      re_delta.x = a0delta+re_delta.x;
                      jp0b = jp0+1;
                      if (jp0b >= np0) jp0b -= np0;
                    }
                    if (re_delta.y < 0.0) {
                      jp1--;
                      if (jp1 < 0) jp1 += np1;
                      re_delta.y = a1delta+re_delta.y;
                      jp1b = jp1+1;
                      if (jp1b >= np1) jp1b -= np1;
                    }
                    if (re_delta.z < 0.0) {
                      jp2--;
                      if (jp2 < 0) jp2 += np2;
                      re_delta.z = a2delta+re_delta.z;
                      jp2b = jp2+1;
                      if (jp2b >= np2) jp2b -= np2;
                    }
                    int jindex = jp0 + np0*jp1 + np0*np1*jp2;
                    //int jindloc = jp0 + np0*jp1 + np0*np1*(jp2-np2first);
                    
                    // simple linear interpolation
                    double t = re_delta.x/a0delta;
                    double u = re_delta.y/a1delta;
                    double v = re_delta.z/a2delta;

                    int jinterp[8];
                    jinterp[0] = jp0 + np0*jp1 + np0*np1*jp2;    // i,j,k      
                    jinterp[1] = jp0b + np0*jp1 + np0*np1*jp2;   // i+1,j,k    
                    jinterp[2] = jp0 + np0*jp1b + np0*np1*jp2;   // i,j+1,k    
                    jinterp[3] = jp0 + np0*jp1 + np0*np1*jp2b;   // i,j,k+1    
                    jinterp[4] = jp0b + np0*jp1b + np0*np1*jp2;  // i+1,j+1,k  
                    jinterp[5] = jp0b + np0*jp1 + np0*np1*jp2b;  // i+1,j,k+1  
                    jinterp[6] = jp0 + np0*jp1b + np0*np1*jp2b;  // i,j+1,k+1  
                    jinterp[7] = jp0b + np0*jp1b + np0*np1*jp2b; // i+1,j+1,k+1
                    
                    complex<double> wfi = (1.-t)*(1.-u)*(1.-v)*twfr[jinterp[0]] + t*(1.-u)*(1.-v)*twfr[jinterp[1]] + (1.-t)*u*(1.-v)*twfr[jinterp[2]] + (1.-t)*(1.-u)*v*twfr[jinterp[3]] + t*u*(1.-v)*twfr[jinterp[4]] + t*(1.-u)*v*twfr[jinterp[5]] + (1.-t)*u*v*twfr[jinterp[6]] + t*u*v*twfr[jinterp[7]];
                    
                    //twfrsum[jindex] += wftmp[i];
                    wftmp2[i] = wfi;
                    
                    ip0++;
                    if (ip0 >= np0) {
                      ip0 = 0;
                      ip1++;
                      if (ip1 >= np1) { 
                        ip1 = 0;
                        ip2++;
                      }
                    }
                  }

                  for (int i = 0; i < np012loc; i++)
                    wftmp[i] = wftmp2[i];
                  
                  ft.forward(&wftmp[0],c.valptr(mloc*n));

                } // for n < nstloc
              }
            }
          }
        }
      }
    }
  }




}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print(ostream& os, string encoding, string tag) const {

  //ewd SPIN
  assert(nspin_==1);

  if ( ctxt_.oncoutpe() ) {
    os << "<" << tag << " ecut=\"" << ecut_ << "\""
       << " nspin=\"" << nspin_ << "\""
       << " nel=\"" << nel_ << "\""
       << " nempty=\"" << nempty_ << "\">" << endl;
    os << setprecision(10);
    os << "<domain a=\""
       << cell_.a(0) << "\"\n        b=\""
       << cell_.a(1) << "\"\n        c=\""
       << cell_.a(2) << "\"/>" << endl;
    os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
       <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
       <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
  }
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          sd_[ispin][ikp]->print(os,encoding);
        }
      }
    }
  }
  
  if ( ctxt_.oncoutpe() )
    os << "</" << tag << ">" << endl;
}


////////////////////////////////////////////////////////////////////////////////
void Wavefunction::printeig(void) {

  const int neig_line = 8;

  int nkpoints = nkp();
  int nkppar = sdcontext_[0].size();
  
  for ( int ispin = 0; ispin < nspin_ ; ispin++) {
    int kcnt = 0;
    for ( int ikp = 0; ikp < nkppar; ikp++ ) {
      int kpe0 = kptproc0_[ispin][ikp];

      if (ctxt_.mype() == kpe0) { // this process has the data we need to print
        if (ctxt_.oncoutpe()) { // data is already on coutpe, print it
          for ( int kloc=0; kloc<nkptloc_; kloc++) {
            int kp = kptloc_[kloc];
            const double eVolt = 2.0 * 13.6056923;
            ostringstream oss3;
            oss3.width(1);  oss3.fill('0');  oss3 << ispin;
            D3vector tkpvec = sd(ispin,kp)->kpoint();
            int neig = sd(ispin,kp)->nst();
            cout <<    "  <eigenvalues spin=\"" << oss3.str() << "\" kpoint=\"" << setprecision(8) << 
                sd(ispin,kp)->kpoint() << "\" n=\"" << neig << "\">" << endl;
            for ( int i = 0; i < neig; i++ ) {
              cout << setw(12) << setprecision(5) << sd(ispin,kp)->eig(i)*eVolt;
              if ( i%neig_line == neig_line-1 ) cout << endl;
            }
            if ( neig%neig_line != 0 ) cout << endl;
            cout << "  </eigenvalues>" << endl;
          }
        }
        else { // send data to coutpe
          for ( int kloc=0; kloc<nkptloc_; kloc++) {
            int kp = kptloc_[kloc];
            D3vector tkpvec = sd(ispin,kp)->kpoint();
            double tkpoint[3] = {tkpvec.x, tkpvec.y, tkpvec.z};
            MPI_Send(&tkpoint,3,MPI_DOUBLE,ctxt_.coutpe(),kpe0,ctxt_.comm());
            int neig = sd(ispin,kp)->nst();
            MPI_Send(&neig,1,MPI_INT,ctxt_.coutpe(),kpe0,ctxt_.comm());
            double teig[neig];
            for (int n=0; n<neig; n++)
              teig[n] = sd(ispin,kp)->eig(n);
            MPI_Send(&teig,neig,MPI_DOUBLE,ctxt_.coutpe(),kpe0,ctxt_.comm());
          }
        }
      }
      if (ctxt_.oncoutpe() && ctxt_.mype() != kpe0) { // receive data from kpe0, print it
        for ( int kloc=0; kloc<nkplocproc0_[ispin][ikp]; kloc++) {
          int kp = kptloc_[kloc];
          double tkpoint[3];
          MPI_Status status;
          MPI_Recv(&tkpoint,3,MPI_DOUBLE,kpe0,kpe0,ctxt_.comm(),&status);
          int neig;
          MPI_Recv(&neig,1,MPI_INT,kpe0,kpe0,ctxt_.comm(),&status);
          double teig[neig];
          MPI_Recv(&teig,neig,MPI_DOUBLE,kpe0,kpe0,ctxt_.comm(),&status);
          const double eVolt = 2.0 * 13.6056923;
          ostringstream oss3;
          oss3.width(1);  oss3.fill('0');  oss3 << ispin;
          cout <<    "  <eigenvalues spin=\"" << oss3.str() << "\" kpoint=\"" << setprecision(8) << tkpoint[0] << " " << tkpoint[1] << " " << tkpoint[2]
               << "\" n=\"" << neig << "\">" << endl;
          for ( int i = 0; i < neig; i++ ) {
            cout << setw(12) << setprecision(5) << teig[i]*eVolt;
            if ( i%neig_line == neig_line-1 ) cout << endl;
          }
          if ( neig%neig_line != 0 ) cout << endl;
          cout << "  </eigenvalues>" << endl;
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::printocc(void) {
  
  const int nocc_line = 10;

  cout.setf(ios::right,ios::adjustfield);
  cout.setf(ios::fixed,ios::floatfield);

  int nkpoints = nkp();
  int nkppar = sdcontext_[0].size();
  
  for ( int ispin = 0; ispin < nspin_ ; ispin++) {
    int kcnt = 0;
    for ( int ikp = 0; ikp < nkppar; ikp++ ) {
      int kpe0 = kptproc0_[ispin][ikp];

      if (ctxt_.mype() == kpe0) { // this process has the data we need to print
        if (ctxt_.oncoutpe()) { // data is already on coutpe, print it
          for ( int kloc=0; kloc<nkptloc_; kloc++) {
            int kp = kptloc_[kloc];
            const double eVolt = 2.0 * 13.6056923;
            ostringstream oss3;
            oss3.width(1);  oss3.fill('0');  oss3 << ispin;
            D3vector tkpvec = sd(ispin,kp)->kpoint();
            int nocc = sd(ispin,kp)->nst();
            cout <<    "  <occupation spin=\"" << oss3.str() << "\" kpoint=\"" << setprecision(8) << 
                sd(ispin,kp)->kpoint() << "\" n=\"" << nocc << "\">" << endl;
            for ( int i = 0; i < nocc; i++ ) {
              cout << setw(8) << setprecision(4) << sd(ispin,kp)->occ(i);
              if ( i%nocc_line == nocc_line-1 ) cout << endl;
            }
            if ( nocc%nocc_line != 0 ) cout << endl;
            cout << "  </occupation>" << endl;
          }
        }
        else { // send data to coutpe
          for ( int kloc=0; kloc<nkptloc_; kloc++) {
            int kp = kptloc_[kloc];
            D3vector tkpvec = sd(ispin,kp)->kpoint();
            double tkpoint[3] = {tkpvec.x, tkpvec.y, tkpvec.z};
            MPI_Send(&tkpoint,3,MPI_DOUBLE,ctxt_.coutpe(),kpe0,ctxt_.comm());
            int nocc = sd(ispin,kp)->nst();
            MPI_Send(&nocc,1,MPI_INT,ctxt_.coutpe(),kpe0,ctxt_.comm());
            double tocc[nocc];
            for (int n=0; n<nocc; n++)
              tocc[n] = sd(ispin,kp)->occ(n);
            MPI_Send(&tocc,nocc,MPI_DOUBLE,ctxt_.coutpe(),kpe0,ctxt_.comm());
          }
        }
      }
      if (ctxt_.oncoutpe() && ctxt_.mype() != kpe0) { // receive data from kpe0, print it
        for ( int kloc=0; kloc<nkplocproc0_[ispin][ikp]; kloc++) {
          int kp = kptloc_[kloc];
          double tkpoint[3];
          MPI_Status status;
          MPI_Recv(&tkpoint,3,MPI_DOUBLE,kpe0,kpe0,ctxt_.comm(),&status);
          int nocc;
          MPI_Recv(&nocc,1,MPI_INT,kpe0,kpe0,ctxt_.comm(),&status);
          double tocc[nocc];
          MPI_Recv(&tocc,nocc,MPI_DOUBLE,kpe0,kpe0,ctxt_.comm(),&status);
          ostringstream oss3;
          oss3.width(1);  oss3.fill('0');  oss3 << ispin;
          cout <<    "  <occupation spin=\"" << oss3.str() << "\" kpoint=\"" << setprecision(8) << tkpoint[0] << " " << tkpoint[1] << " " << tkpoint[2]
               << "\" n=\"" << nocc << "\">" << endl;
          for ( int i = 0; i < nocc; i++ ) {
            cout << setw(8) << setprecision(4) << tocc[i];
            if ( i%nocc_line == nocc_line-1 ) cout << endl;
          }
          if ( nocc%nocc_line != 0 ) cout << endl;
          cout << "  </occupation>" << endl;
        }
      }
      cout << setprecision(8);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write(SharedFilePtr& sfp, string encoding, string tag) const {
  sfp.sync();

  if ( ctxt_.oncoutpe() ) {
    ostringstream os;
    os << "<" << tag << " ecut=\"" << ecut_ << "\""
       << " nspin=\"" << nspin_ << "\""
       << " delta_spin=\"" << deltaspin_ << "\""
       << " nel=\"" << nel_ << "\""
       << " nempty=\"" << nempty_ << "\">" << endl;
    os << setprecision(10);
    os << "<domain a=\""
       << cell_.a(0) << "\"\n        b=\""
       << cell_.a(1) << "\"\n        c=\""
       << cell_.a(2) << "\"/>" << endl;
    if ( refcell_.volume() != 0.0 )
      {
        os << "<reference_domain a=\""
           << refcell_.a(0) << "\"\n        b=\""
           << refcell_.a(1) << "\"\n        c=\""
           << refcell_.a(2) << "\"/>" << endl;
      }
    os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
       <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
       <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
    string str(os.str());
    int len = str.size();
    MPI_Status status;
    int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)str.c_str(),
                                len,MPI_CHAR,&status);
    if ( err != 0 )
      cout << " Wavefunction::write: error in MPI_File_write" << endl;
    sfp.advance(len);
  }


  for ( int ikp = 0; ikp < nkp(); ikp++ ) {
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin) && kptactive(ikp)) {
        assert(sd_[ispin][ikp] != 0);
        sd_[ispin][ikp]->write(sfp,encoding,weight_[ikp],ispin,nspin_);
      }
      sfp.sync();
    }
  }
  sfp.sync();

  /*
  //ewd:  this worked before spin was added
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        for ( int kloc=0; kloc<nkptloc_list_[ikp]; kloc++) {
          if (sdcontext_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() ) {
            int kp = kptloc_[kloc];
            if ( sd_[ispin][kp] != 0 ) 
              sd_[ispin][kp]->write(sfp,encoding,weight_[kp],ispin,nspin_);
          }
          sfp.sync();
        }
      }
    }
  }
  sfp.sync();
  */

  
  /*
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
    for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
      if (sdcontext_[ispin][ikp] != 0 ) {
        if (sdcontext_[ispin][ikp]->oncoutpe() ) {
          ostringstream os;
          os << "</" << tag << ">" << endl;
          string str(os.str());
          int len = str.size();
          MPI_Status status;
          int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)str.c_str(),
                                  len,MPI_CHAR,&status);
          if ( err != 0 )
            cout << " Wavefunction::write: error in MPI_File_write" << endl;
          sfp.advance(len);
        }
      }
    }
    }
  }
  */

  if ( ctxt_.oncoutpe() ) {
      ostringstream os;
      os << "</" << tag << ">" << endl;
      string str(os.str());
      int len = str.size();
      MPI_Status status;
      int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)str.c_str(),
                                  len,MPI_CHAR,&status);
      if ( err != 0 )
        cout << " Wavefunction::write: error in MPI_File_write" << endl;
      sfp.advance(len);
    }
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write_dump(string filebase) {

  // make unique filename for each process by appending process number to filename

  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  ostringstream oss;
  oss.width(6);  oss.fill('0');  oss << mype;
  string mypefile = filebase + oss.str(); 
  ofstream os;

  bool ifempty = (nempty_ > 0);

  os.open(mypefile.c_str(),ofstream::binary);
  //fopenFILE* PEFILE = fopen(mypefile.c_str(),"wb");
  
  // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
  os.write(mypefile.c_str(),sizeof(char)*mypefile.length());
  os.flush();
#endif
  
  int wkploc = nkptloc_;
  int kp0 = kptloc_[0];         // global index of local kpoint
  int wnst = sd(0,kp0)->nst();
  //os.write((char*)&wkploc,sizeof(int));
  //os.write((char*)&wnst,sizeof(int));

  // write out wave function
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);
          int mloc = sd_[ispin][ikp]->c().mloc();
          int nloc = sd_[ispin][ikp]->c().nloc();
          int ngwloc = sd_[ispin][ikp]->basis().localsize();
          const complex<double>* p = sd_[ispin][ikp]->c().cvalptr();

          // implement single large write for less I/O node contention on BG/Q
          //for ( int n = 0; n < nloc; n++ )
          //   os.write((char*)&p[n*mloc],sizeof(complex<double>)*ngwloc);

          os.write((char*)&p[0],sizeof(complex<double>)*nloc*mloc);
          os.flush();
          
          //fopen for ( int n = 0; n < nloc; n++ )
          //fopen    fwrite(&p[n*mloc],sizeof(complex<double>),ngwloc,PEFILE);
          
          // if there are empty states, save occupation and eigenvalues
          if (ifempty) {
            int nst = sd_[ispin][ikp]->nst();
            const double* peig = sd_[ispin][ikp]->eig_ptr();
            const double* pocc = sd_[ispin][ikp]->occ_ptr();
            os.write((char*)&peig[0],sizeof(double)*nst);
            os.write((char*)&pocc[0],sizeof(double)*nst);
            os.flush();
            //fopen fwrite(&peig[0],sizeof(double),nst,PEFILE);
            //fopen fwrite(&pocc[0],sizeof(double),nst,PEFILE);
          }
        }
      }
    }
  }  
  os.close();
  //fopen fclose(PEFILE);

}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write_fast(string filebase) {

   // write_fast uses C fwrite calls to dump checkpoint data to a subset of
   // files (ideal for machines like BG/Q where the number of I/O nodes is << npes)

   int mype, npes;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  MPI_Comm_size(MPI_COMM_WORLD,&npes);
#else
  mype = 0;
  npes = 1;
#endif

  int nFiles = npes;
  int filenum = mype;
  int writerTask = mype;
  int nTasksPerFile = 1;
#ifdef BGQ
  nTasksPerFile = 32;
  nFiles = npes/nTasksPerFile;
#endif

  if (mype == 0)
     cout << "Wavefunction::write_fast:  writing data from " << npes << " tasks to " << nFiles << " files." << endl;

  // if nFiles != npes, need to calculate which tasks write to which file
  if (nFiles != npes)
  {
     nTasksPerFile = ( npes%nFiles == 0 ? npes/nFiles : npes/nFiles + 1);
     assert(nTasksPerFile > 1);
     filenum = mype/nTasksPerFile;
     writerTask = filenum*nTasksPerFile;
  }
  
  // all tasks send their data to their writer task, who writes it to file

  ostringstream oss;
  oss.width(6);  oss.fill('0');  oss << filenum;
  string mypefile = filebase + oss.str(); 
  ofstream os;
  bool ifempty = (nempty_ > 0);

  if (mype == writerTask)  // open file
     os.open(mypefile.c_str(),ofstream::binary);
     //fopenFILE* PEFILE = fopen(mypefile.c_str(),"wb");
    
  // write out wave function
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp=0; ikp<nkp(); ikp++) {
        if (kptactive(ikp)) {
          assert(sd_[ispin][ikp] != 0);

          if (mype == writerTask)  // write local data
          {
          int mloc = sd_[ispin][ikp]->c().mloc();
          int nloc = sd_[ispin][ikp]->c().nloc();
          const complex<double>* p = sd_[ispin][ikp]->c().cvalptr();

          os.write((char*)&p[0],sizeof(complex<double>)*nloc*mloc);
          os.flush();
          //fopen for ( int n = 0; n < nloc; n++ )
          //fopen    fwrite(&p[n*mloc],sizeof(complex<double>),ngwloc,PEFILE);
          }
          
          for (int jj=1; jj<nTasksPerFile; jj++)
          {
             int dataTask = jj + writerTask;
             if (mype == dataTask)
             {
                int mloc = sd_[ispin][ikp]->c().mloc();
                int nloc = sd_[ispin][ikp]->c().nloc();
                int nComplex = mloc*nloc;
                MPI_Send(&nComplex,1,MPI_INT,writerTask,dataTask,MPI_COMM_WORLD);
                complex<double>* p = sd_[ispin][ikp]->c().valptr();
                MPI_Send(&p[0],nComplex,MPI_DOUBLE_COMPLEX,writerTask,dataTask,MPI_COMM_WORLD);
             }
             if (mype == writerTask)
             {
                int nComplex;
                MPI_Recv(&nComplex,1,MPI_INT,dataTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                vector<complex<double> > data;
                data.resize(nComplex);
                MPI_Recv(&data[0],nComplex,MPI_DOUBLE_COMPLEX,dataTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                os.write((char*)&data[0],sizeof(complex<double>)*nComplex);
                os.flush();                
             }
          }
          
          // if there are empty states, save occupation and eigenvalues
          if (ifempty && mype == writerTask) {
             int nst = sd_[ispin][ikp]->nst();
             const double* peig = sd_[ispin][ikp]->eig_ptr();
             const double* pocc = sd_[ispin][ikp]->occ_ptr();
             os.write((char*)&peig[0],sizeof(double)*nst);
             os.write((char*)&pocc[0],sizeof(double)*nst);
             os.flush();
          }
        }
      }
    }
  }  
  
  if (mype == writerTask)  // close file
     os.close();
  //fopen fclose(PEFILE);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write_states(string filebase, string format) {

   int mype, npes;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  MPI_Comm_size(MPI_COMM_WORLD,&npes);
#else
  mype = 0;
  npes = 1;
#endif

  bool ifempty = (nempty_ > 0);
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {
            Context* tctxt = sdcontext_[ispin][ikp];
            int nprow = tctxt->nprow();
            int npcol = tctxt->npcol();
            int prow = tctxt->myrow();
            int pcol = tctxt->mycol();
            int writerTask = mype - (mype%nprow);
            assert(writerTask >= 0);
            if (mype == writerTask)
               assert(prow == 0);
            
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];         // global index of local kpoint
              if ( sd_[ispin][kp] != 0 ) {
                int nstloc = sd_[ispin][kp]->nstloc();
                int nb = sd_[ispin][kp]->c().nb();

                const Basis& basis = sd_[ispin][kp]->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));

                for ( int lj=0; lj < sd_[ispin][kp]->c().nblocks(); lj++ )
                {
                   for ( int jj=0; jj < sd_[ispin][kp]->c().nbs(lj); jj++ )
                   {
                      // global state index
                      const int nglobal = sd_[ispin][kp]->c().j(lj,jj);
                      const int norig = lj*sd_[ispin][kp]->c().nb()+jj;

                      ofstream os;
                      if (mype == writerTask) {
                         // write out wavefunction for this state and k-point
                         ostringstream oss1,oss2,oss3;
                         oss1.width(5);  oss1.fill('0');  oss1 << nglobal;
                         oss2.width(4);  oss2.fill('0');  oss2 << kp;
                         oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                         string statefile;
                         if (nspin_ == 1) 
                            statefile = filebase + "k" + oss2.str() + "n" + oss1.str(); 
                         else
                            statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + "n" + oss1.str(); 
                         
                         if (format == "molmol" || format == "text") 
                            statefile = statefile + ".iso";
                         else if (format == "gopenmol") 
                            statefile = statefile + ".plt";
                         
                         if (format == "binary") 
                            os.open(statefile.c_str(),ofstream::binary);
                         else {
                            os.open(statefile.c_str(),ofstream::out);
                            os.setf(ios::scientific,ios::floatfield);
                            os << setprecision(8);
                         }

                         // headers for visualization formats
                         if (format == "molmol" || format == "text") {
                            D3vector a0 = cell_.a(0);
                            D3vector a1 = cell_.a(1);
                            D3vector a2 = cell_.a(2);
                            if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                                a2.x != 0.0 || a2.y != 0.0 )
                               os << "Error writing header:  Molmol isosurface requires rectangular box!" << endl;
                            double dx = a0.x/(double)ft.np0();
                            double dy = a1.y/(double)ft.np1();
                            double dz = a2.z/(double)ft.np2();
                            //    origin    npoints       grid spacing
                            os << "0.0 " << ft.np0() << " " << dx << endl;
                            os << "0.0 " << ft.np1() << " " << dy << endl;
                            os << "0.0 " << ft.np2() << " " << dz << endl;
                         }
                         else if (format == "gopenmol") {
                            D3vector a0 = cell_.a(0);
                            D3vector a1 = cell_.a(1);
                            D3vector a2 = cell_.a(2);
                            if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                                a2.x != 0.0 || a2.y != 0.0 )
                               os << "Error writing header:  gOpenMol isosurface requires rectangular box!" << endl;
                            os << "3 200" << endl;  //ewd copying this from another utility, not sure what it means
                            os << ft.np0() << " " << ft.np1() << " " << ft.np2() << endl;
                            os << "0.0 " << a2.z << endl;
                            os << "0.0 " << a1.y << endl;
                            os << "0.0 " << a0.x << endl;
                         }
                      }
                   
                      int mloc = sd_[ispin][kp]->c().mloc();
                      vector<complex<double> > wftmp(ft.np012loc());
                  
                      ComplexMatrix& c = sd_[ispin][kp]->c();
                      ft.backward(c.cvalptr(mloc*norig),&wftmp[0]);

                      // write out wave function
                      if (mype == writerTask)  // write local data
                      {
                         int size = ft.np012loc();
                         os.write((char*)&wftmp[0],sizeof(complex<double>)*size);
                         os.flush();
                      }
          
                      for (int jj=1; jj<nprow; jj++)
                      {
                         int dataTask = jj + writerTask;
                         if (mype == writerTask)
                         {
                            int nComplex;
                            MPI_Recv(&nComplex,1,MPI_INT,dataTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                            vector<complex<double> > data;
                            data.resize(nComplex);
                            MPI_Recv(&data[0],nComplex,MPI_DOUBLE_COMPLEX,dataTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                            os.write((char*)&data[0],sizeof(complex<double>)*nComplex);
                            os.flush();                
                         }
                         if (mype == dataTask)
                         {
                            int size = ft.np012loc();
                            MPI_Send(&size,1,MPI_INT,writerTask,dataTask,MPI_COMM_WORLD);
                            MPI_Send(&wftmp[0],size,MPI_DOUBLE_COMPLEX,writerTask,dataTask,MPI_COMM_WORLD);
                         }
                      }
                      if (mype == writerTask)  // close file
                         os.close();
                   }
                }

                // if there are empty states, save occupation and eigenvalues
                if (ifempty && prow == 0 && pcol == 0 && format == "binary") {
                   // write out occupation for this state and k-point
                   ofstream os;
                   ostringstream oss2,oss3;
                   oss2.width(4);  oss2.fill('0');  oss2 << kp;
                   oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                   string statefile;
                   if (nspin_ == 1) 
                      statefile = filebase + "k" + oss2.str() + ".occ"; 
                   else
                      statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + ".occ";
                   os.open(statefile.c_str(),ofstream::binary);
                   
                   int nst = sd_[ispin][kp]->nst();
                   const double* peig = sd_[ispin][kp]->eig_ptr();
                   const double* pocc = sd_[ispin][kp]->occ_ptr();
                   os.write((char*)&peig[0],sizeof(double)*nst);
                   os.write((char*)&pocc[0],sizeof(double)*nst);
                   os.flush();
                   os.close();
                }
              }
            }
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::read_dump(string filebase) {

  if (!hasdata_) {
    hasdata_ = true;
    allocate();
  }

  // get unique filename for each process by appending process number to filename

  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  ostringstream oss;
  oss.width(6);  oss.fill('0');  oss << mype;
  string mypefile = filebase + oss.str(); 
  ifstream is;

  bool ifempty = (nempty_ > 0);

  is.open(mypefile.c_str(),ofstream::binary);

  if (is.is_open()) {

     // hack to make checkpointing work with BlueGene compilers
#ifdef BGQ
     int len = mypefile.length();
     char* tmpfilename = new char[256];
     is.read(tmpfilename,sizeof(char)*mypefile.length());
#endif
    int wkploc = nkptloc_;
    int kp0 = kptloc_[0];         // global index of local kpoint
    int wnst = sd(0,kp0)->nst();
    //is.read((char*)&wkploc,sizeof(int));
    //is.read((char*)&wnst,sizeof(int));

    // read in wave function
    for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
          if (kptactive(ikp)) {
            assert(sd_[ispin][ikp] != 0);
            int mloc = sd_[ispin][ikp]->c().mloc();
            int nloc = sd_[ispin][ikp]->c().nloc();
            int ngwloc = sd_[ispin][ikp]->basis().localsize();
            const complex<double>* p = sd_[ispin][ikp]->c().cvalptr();
            for ( int n = 0; n < nloc; n++ )
              is.read((char*)&p[n*mloc],sizeof(complex<double>)*ngwloc);
            
            // if there are empty states, load occupation and eigenvalues
            if (ifempty) {
              int nst = sd_[ispin][ikp]->nst();
              const double* peig = sd_[ispin][ikp]->eig_ptr();
              const double* pocc = sd_[ispin][ikp]->occ_ptr();
              is.read((char*)&peig[0],sizeof(double)*nst);
              is.read((char*)&pocc[0],sizeof(double)*nst);
            }
          }
        }
      }
    }  
  }
  else {
     if ( ctxt_.oncoutpe())
        cout << "<!-- LoadCmd: " << filebase << " checkpoint files not found, skipping load. -->" << endl;
  }
  is.close();
  
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::read_fast(string filebase) {

  if (!hasdata_) {
    hasdata_ = true;
    allocate();
  }
   // read_fast reads checkpoint data written using write_fast, distributes

   int mype, npes;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  MPI_Comm_size(MPI_COMM_WORLD,&npes);
#else
  mype = 0;
  npes = 1;
#endif

  int nFiles = npes;
  int filenum = mype;
  int readerTask = mype;
  int nTasksPerFile = 1;
#ifdef BGQ
  nTasksPerFile = 32;
  nFiles = npes/nTasksPerFile;
#endif

  if (mype == 0)
     cout << "Wavefunction::read_fast:  reading data for " << npes << " tasks from " << nFiles << " files." << endl;

  // if nFiles != npes, need to calculate which tasks wrote to which file
  if (nFiles != npes)
  {
     nTasksPerFile = ( npes%nFiles == 0 ? npes/nFiles : npes/nFiles + 1);
     assert(nTasksPerFile > 1);
     filenum = mype/nTasksPerFile;
     readerTask = filenum*nTasksPerFile;
  }
  
  ostringstream oss;
  oss.width(6);  oss.fill('0');  oss << filenum;
  string mypefile = filebase + oss.str(); 
  ifstream is;
  bool ifempty = (nempty_ > 0);

  if (mype == readerTask)
     is.open(mypefile.c_str(),ofstream::binary);

  // read in wave function
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
     if (spinactive(ispin)) {
        for ( int ikp=0; ikp<nkp(); ikp++) {
           if (kptactive(ikp)) {
              assert(sd_[ispin][ikp] != 0);
              
              int fileFound = -1;
              if (mype == readerTask)  // read local data
              {
                 int mloc = sd_[ispin][ikp]->c().mloc();
                 int nloc = sd_[ispin][ikp]->c().nloc();
                 const complex<double>* p = sd_[ispin][ikp]->c().cvalptr();
                 if (is.is_open())
                 {
                    fileFound = 1;
                    is.read((char*)&p[0],sizeof(complex<double>)*nloc*mloc);
                 }
                 else
                 {
                    fileFound = 0;
                    if ( ctxt_.oncoutpe())
                       cout << "<!-- LoadCmd: " << filebase << " checkpoint files not found, skipping load. -->" << endl;
                 }                       
              }                       
               
              for (int jj=1; jj<nTasksPerFile; jj++)
              {
                 int dataTask = jj + readerTask;
                 if (mype == dataTask)
                 {
                    MPI_Recv(&fileFound,1,MPI_INT,readerTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    if (fileFound == 1)
                    {
                       int mloc = sd_[ispin][ikp]->c().mloc();
                       int nloc = sd_[ispin][ikp]->c().nloc();
                       int nComplex = mloc*nloc;
                       MPI_Send(&nComplex,1,MPI_INT,readerTask,dataTask,MPI_COMM_WORLD);
                       complex<double>* p = sd_[ispin][ikp]->c().valptr();
                       MPI_Recv(&p[0],nComplex,MPI_DOUBLE_COMPLEX,readerTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    }
                 }
                 if (mype == readerTask)
                 {
                    MPI_Send(&fileFound,1,MPI_INT,dataTask,dataTask,MPI_COMM_WORLD);
                    if (fileFound == 1)
                    {
                       int nComplex;
                       MPI_Recv(&nComplex,1,MPI_INT,dataTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                       vector<complex<double> > data;
                       data.resize(nComplex);
                       if (is.is_open())
                          is.read((char*)&data[0],sizeof(complex<double>)*nComplex);
                       MPI_Send(&data[0],nComplex,MPI_DOUBLE_COMPLEX,dataTask,dataTask,MPI_COMM_WORLD);
                    }
                 }
              }
              if (fileFound == 0)
                 return;
              
              // if there are empty states, load occupation and eigenvalues
              if (ifempty)
              {
                 int nst = sd_[ispin][ikp]->nst();
                 double* peig = (double*)sd_[ispin][ikp]->eig_ptr();
                 double* pocc = (double*)sd_[ispin][ikp]->occ_ptr();
                 if (mype == readerTask) {
                    is.read((char*)&peig[0],sizeof(double)*nst);
                    is.read((char*)&pocc[0],sizeof(double)*nst);
                    for (int jj=1; jj<nTasksPerFile; jj++)
                    {
                       int dataTask = jj + readerTask;
                       MPI_Send(&peig[0],nst,MPI_DOUBLE,dataTask,dataTask,MPI_COMM_WORLD);
                       MPI_Send(&pocc[0],nst,MPI_DOUBLE,dataTask,dataTask,MPI_COMM_WORLD);
                    }
                 }
                 else
                 {
                    MPI_Recv(&peig[0],nst,MPI_DOUBLE,readerTask,mype,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    MPI_Recv(&pocc[0],nst,MPI_DOUBLE,readerTask,mype,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                 }
              }
           }
        }
     }  
  }
  
  if (mype == readerTask)
     is.close();
  
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::read_states(string filebase) {
   
   if (!hasdata_) {
      hasdata_ = true;
      allocate();
   }
   int mype, npes;
#if USE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&mype);
   MPI_Comm_size(MPI_COMM_WORLD,&npes);
#else
   mype = 0;
   npes = 1;
#endif
   
   bool ifempty = (nempty_ > 0);
   if (!hasdata_) {
      hasdata_ = true;
      allocate();
   }

   int checkPointFound = 0;
   for ( int ispin = 0; ispin < nspin_; ispin++ ) {
      if (spinactive(ispin)) {
         for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
            if (sdcontext_[ispin][ikp] != 0 ) {
               if (sdcontext_[ispin][ikp]->active() ) {
                  Context* tctxt = sdcontext_[ispin][ikp];
                  int nprow = tctxt->nprow();
                  int npcol = tctxt->npcol();
                  int prow = tctxt->myrow();
                  int pcol = tctxt->mycol();
                  int readerTask = mype - (mype%nprow);
                  assert(readerTask >= 0);
                  if (mype == readerTask)
                     assert(prow == 0);

                  for ( int kloc=0; kloc<nkptloc_; kloc++) {
                     int kp = kptloc_[kloc];         // global index of local kpoint
                     if ( sd_[ispin][kp] != 0 ) {
                        int nstloc = sd_[ispin][kp]->nstloc();
                        int nb = sd_[ispin][kp]->c().nb();

                        const Basis& basis = sd_[ispin][kp]->basis();
                        FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
                        
                        for ( int n = 0; n < nstloc; n++ ) {
                           // global n index
                           const int nn = sd(ispin,kp)->c().j(0,n);
                           
                           vector<complex<double> > wftmp(ft.np012loc());

                           ifstream is;
                           int fileFound = 0;
                           if (mype == readerTask) {
                              // read in wavefunction for this state and k-point
                              ostringstream oss1,oss2,oss3;
                              oss1.width(5);  oss1.fill('0');  oss1 << nn;
                              oss2.width(4);  oss2.fill('0');  oss2 << kp;
                              oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                              string statefile;
                              if (nspin_ == 1) 
                                 statefile = filebase + "k" + oss2.str() + "n" + oss1.str(); 
                              else
                                 statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + "n" + oss1.str(); 
                              is.open(statefile.c_str(),ofstream::binary);
                              if (is.is_open()) {
                                 // read local data
                                 int size = ft.np012loc();
                                 is.read((char*)&wftmp[0],sizeof(complex<double>)*size);
                                 fileFound = 1;
                              }
                              else {
                                 fileFound = -1;
                                 checkPointFound = -1;
                                 if ( ctxt_.oncoutpe())
                                    cout << "<!-- LoadCmd: " << filebase << " checkpoint files not found, skipping load. -->" << endl;
                              }
                           }

                           for (int jj=1; jj<nprow; jj++)
                           {
                              int dataTask = jj + readerTask;
                              if (mype == readerTask)
                              {
                                 MPI_Send(&fileFound,1,MPI_INT,dataTask,dataTask,MPI_COMM_WORLD);
                                 if (fileFound == 1)
                                 {
                                    int nComplex;
                                    MPI_Recv(&nComplex,1,MPI_INT,dataTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                    vector<complex<double> > data;
                                    data.resize(nComplex);
                                    is.read((char*)&data[0],sizeof(complex<double>)*nComplex);
                                    MPI_Send(&data[0],nComplex,MPI_DOUBLE_COMPLEX,dataTask,dataTask,MPI_COMM_WORLD);
                                 }
                              }
                              if (mype == dataTask)
                              {
                                 MPI_Recv(&fileFound,1,MPI_INT,readerTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                 if (fileFound == 1)
                                 {
                                    int size = ft.np012loc();
                                    MPI_Send(&size,1,MPI_INT,readerTask,dataTask,MPI_COMM_WORLD);
                                    MPI_Recv(&wftmp[0],size,MPI_DOUBLE_COMPLEX,readerTask,dataTask,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                                 }
                              }
                           }

                           // this causes a hang when nstloc = 0 on some tasks
                           //if (fileFound <= 0)
                           //   return;

                           if (fileFound == 1)
                           {
                              if (mype == readerTask)
                                 is.close();
                  
                              //ewd when k=0, force complex part to zero
                              //if (basis.real()) 
                              //   for ( int i = 0; i < ft.np012loc(); i++ )
                              //      wftmp[i] = complex<double>(real(wftmp[i]),0.0);
                           
                              ComplexMatrix& c = sd_[ispin][kp]->c();
                              int mloc = sd_[ispin][kp]->c().mloc();
                              ft.forward(&wftmp[0],c.valptr(mloc*n));
                           }
                        } // for n < nstloc

                        // sync up file found info across all tasks
                        tctxt->imin(1,1,&checkPointFound,1);
                        if (checkPointFound < 0)
                           return;
                        
                        // if there are empty states, load occupation and eigenvalues
                        if (ifempty)
                        {
                           ifstream is;
                           int nProcs = nprow*npcol;
                           int procZero = mype - pcol*nprow - prow; // want the first task in this spin, kpoint context
                           if (mype == procZero) {
                              // get occupation for this state and k-point
                              ostringstream oss2,oss3;
                              oss2.width(4);  oss2.fill('0');  oss2 << kp;
                              oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                              string statefile;
                              if (nspin_ == 1) 
                                 statefile = filebase + "k" + oss2.str() + ".occ"; 
                              else
                                 statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + ".occ";
                              is.open(statefile.c_str(),ofstream::binary);
                           }
                           
                           int nst = sd_[ispin][ikp]->nst();
                           double* peig = (double*)sd_[ispin][ikp]->eig_ptr();
                           double* pocc = (double*)sd_[ispin][ikp]->occ_ptr();
                           if (mype == procZero) {
                              is.read((char*)&peig[0],sizeof(double)*nst);
                              is.read((char*)&pocc[0],sizeof(double)*nst);
                              is.close();
                              for (int jj=1; jj<nProcs; jj++)
                              {
                                 int dataTask = jj + procZero;
                                 MPI_Send(&peig[0],nst,MPI_DOUBLE,dataTask,dataTask,MPI_COMM_WORLD);
                                 MPI_Send(&pocc[0],nst,MPI_DOUBLE,dataTask,dataTask,MPI_COMM_WORLD);
                              }
                           }
                           else
                           {
                              MPI_Recv(&peig[0],nst,MPI_DOUBLE,procZero,mype,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                              MPI_Recv(&pocc[0],nst,MPI_DOUBLE,procZero,mype,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
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
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write_states_old(string filebase, string format) {
  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {
            Context* tctxt = sdcontext_[ispin][ikp];
            int prow = tctxt->myrow();
            int pcol = tctxt->mycol();
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];         // global index of local kpoint
              if ( sd_[ispin][kp] != 0 ) {
                int nstloc = sd_[ispin][kp]->nstloc();
                int nb = sd_[ispin][kp]->c().nb();

                const Basis& basis = sd_[ispin][kp]->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));

                for ( int n = 0; n < nstloc; n++ ) {
                  // global n index
                  const int nn = pcol*nb + n;

                  ofstream os;
                  if (tctxt->myrow() == 0) {
                    // write out wavefunction for this state and k-point
                    ostringstream oss1,oss2,oss3;
                    oss1.width(5);  oss1.fill('0');  oss1 << nn;
                    oss2.width(4);  oss2.fill('0');  oss2 << kp;
                    oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                    string statefile;
                    if (nspin_ == 1) 
                      statefile = filebase + "k" + oss2.str() + "n" + oss1.str(); 
                    else
                      statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + "n" + oss1.str(); 
                    
                    if (format == "molmol" || format == "text") 
                      statefile = statefile + ".iso";
                    else if (format == "gopenmol") 
                      statefile = statefile + ".plt";

                    if (format == "binary") 
                      os.open(statefile.c_str(),ofstream::binary);
                    else {
                      os.open(statefile.c_str(),ofstream::out);
                      os.setf(ios::scientific,ios::floatfield);
                      os << setprecision(8);
                    }

                    // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                    if (format == "binary") {
                       os.write(statefile.c_str(),sizeof(char)*statefile.length());
                       os.flush();
                    }
#endif

                    // headers for visualization formats
                    if (format == "molmol" || format == "text") {
                      D3vector a0 = cell_.a(0);
                      D3vector a1 = cell_.a(1);
                      D3vector a2 = cell_.a(2);
                      if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                          a2.x != 0.0 || a2.y != 0.0 )
                        os << "Error writing header:  Molmol isosurface requires rectangular box!" << endl;
                      double dx = a0.x/(double)ft.np0();
                      double dy = a1.y/(double)ft.np1();
                      double dz = a2.z/(double)ft.np2();
                      //    origin    npoints       grid spacing
                      os << "0.0 " << ft.np0() << " " << dx << endl;
                      os << "0.0 " << ft.np1() << " " << dy << endl;
                      os << "0.0 " << ft.np2() << " " << dz << endl;
                    }
                    else if (format == "gopenmol") {
                      D3vector a0 = cell_.a(0);
                      D3vector a1 = cell_.a(1);
                      D3vector a2 = cell_.a(2);
                      if( a0.y != 0.0 || a0.z != 0.0 || a1.x != 0.0 || a1.z != 0.0 ||
                          a2.x != 0.0 || a2.y != 0.0 )
                        os << "Error writing header:  gOpenMol isosurface requires rectangular box!" << endl;
                      os << "3 200" << endl;  //ewd copying this from another utility, not sure what it means
                      os << ft.np0() << " " << ft.np1() << " " << ft.np2() << endl;
                      os << "0.0 " << a2.z << endl;
                      os << "0.0 " << a1.y << endl;
                      os << "0.0 " << a0.x << endl;
                    }
                  }     
                  int mloc = sd_[ispin][kp]->c().mloc();
                  vector<complex<double> > wftmp(ft.np012loc());
                  vector<double> wftmpr(2*ft.np012loc());
                  
                  ComplexMatrix& c = sd_[ispin][kp]->c();
                  ft.backward(c.cvalptr(mloc*n),&wftmp[0]);

                  // copy to double array
                  double *a = (double*) &wftmp[0];
                  for ( int i = 0; i < 2*ft.np012loc(); i++ )
                    wftmpr[i] = a[i];
                  
                  // send data to first proc in context column
                  for ( int i = 0; i < tctxt->nprow(); i++ ) {
                    if ( i == prow ) {
                      int size = 2*ft.np012loc();
                      tctxt->isend(1,1,&size,1,0,pcol);
                      tctxt->dsend(size,1,&wftmpr[0],1,0,pcol);
                    }
                  }
                  // receive data, print to file
                  if (tctxt->myrow() == 0) {
                    for ( int i = 0; i < tctxt->nprow(); i++ ) {
                      int size = 0;
                      tctxt->irecv(1,1,&size,1,i,pcol);
                      tctxt->drecv(size,1,&wftmpr[0],1,i,pcol);
                      if (format == "binary")
                      {
                        os.write((char*)&wftmpr[0],sizeof(double)*size);
                        os.flush();
                      }
                      else if (format == "molmol" || format == "gopenmol")
                      {
                        // write out |wf|^2 on grid, with x varying fastest
                        ostringstream oss;
                        oss.setf(ios::scientific,ios::floatfield);
                        oss << setprecision(5);
                        for (int j=0; j<size; j+=2)
                          oss << wftmpr[j]*wftmpr[j]+wftmpr[j+1]*wftmpr[j+1] << endl;
                        string tos = oss.str();
                        os.write(tos.c_str(),tos.length());
                        os.flush();
                      }
                      else if (format == "text")
                      {
                        // write out wf on grid, with x varying fastest
                        ostringstream oss;
                        oss.setf(ios::scientific,ios::floatfield);
                        oss << setprecision(5);
                        for (int j=0; j<size; j+=2)
                          os << wftmpr[j] << "  " << wftmpr[j+1] << "  " << wftmpr[j]*wftmpr[j]+wftmpr[j+1]*wftmpr[j+1] << endl;
                        string tos = oss.str();
                        os.write(tos.c_str(),tos.length());
                        os.flush();
                      }
                    }
                    os.close();
                  }
                  
                } // for n < nstloc
              
                // if there are empty states, save occupation and eigenvalues
                if (nempty_ > 0 || format == "binary") {
                  ofstream os;
                  if (tctxt->myrow() == 0) {
                    // write out occupation for this state and k-point
                    ostringstream oss2,oss3;
                    oss2.width(4);  oss2.fill('0');  oss2 << kp;
                    oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                    string statefile;
                    if (nspin_ == 1) 
                      statefile = filebase + "k" + oss2.str() + ".occ"; 
                    else
                      statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + ".occ";
                    os.open(statefile.c_str(),ofstream::binary);
                  
                    int nst = sd_[ispin][kp]->nst();
                    const double* peig = sd_[ispin][kp]->eig_ptr();
                    const double* pocc = sd_[ispin][kp]->occ_ptr();
                    // hack to make checkpointing work w. BlueGene compilers
#ifdef BGQ
                    os.write(statefile.c_str(),sizeof(char)*statefile.length());
#endif
                    os.write((char*)&peig[0],sizeof(double)*nst);
                    os.write((char*)&pocc[0],sizeof(double)*nst);
                    os.flush();
                    os.close();
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

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::read_states_old(string filebase) {

  if (!hasdata_) {
    hasdata_ = true;
    allocate();
  }

  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {
            
            Context* tctxt = sdcontext_[ispin][ikp];
            int prow = tctxt->myrow();
            int pcol = tctxt->mycol();
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];         // global index of local kpoint
              if ( sd_[ispin][kp] != 0 ) {
                int nstloc = sd_[ispin][kp]->nstloc();
                int nb = sd_[ispin][kp]->c().nb();

                const Basis& basis = sd_[ispin][kp]->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));

                for ( int n = 0; n < nstloc; n++ ) {
                  // global n index
                  const int nn = pcol*nb + n;

                  int mloc = sd_[ispin][kp]->c().mloc();
                  vector<double> wftmpr(2*ft.np012loc());

                  ifstream is;
                  if (tctxt->myrow() == 0) {
                    // read in wavefunction for this state and k-point
                    ostringstream oss1,oss2,oss3;
                    oss1.width(5);  oss1.fill('0');  oss1 << nn;
                    oss2.width(4);  oss2.fill('0');  oss2 << kp;
                    oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                    string statefile;
                    if (nspin_ == 1) 
                      statefile = filebase + "k" + oss2.str() + "n" + oss1.str(); 
                    else
                      statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + "n" + oss1.str(); 
                    is.open(statefile.c_str(),ofstream::binary);
                    if (is.is_open()) {

                      // hack to make checkpointing work with BlueGene compilers
#ifdef BGQ
                      int len = statefile.length();
                      char* tmpfilename = new char[256];
                      is.read(tmpfilename,sizeof(char)*statefile.length());
#endif
                       
                      // open files and reads data
                      for ( int i = 0; i < tctxt->nprow(); i++ ) {
                        int size = 2*ft.np2_loc(i)*ft.np0()*ft.np1();
                        is.read((char*)&wftmpr[0],sizeof(double)*size);
                        tctxt->isend(1,1,&size,1,i,pcol);
                        if (size > 0)
                          tctxt->dsend(size,1,&wftmpr[0],1,i,pcol);
                      }
                    }
                    else { // file not found
                      if ( ctxt_.oncoutpe() && n == 0 && ikp == 0 && kloc == 0)
                         cout << "<!-- LoadCmd: " << filebase << " checkpoint files not found, skipping load. -->" << endl;

                      for ( int i = 0; i < tctxt->nprow(); i++ ) {
                        int size = -1;
                        tctxt->isend(1,1,&size,1,i,pcol);
                      }
                    }
                  }

                  // receive data, fill arrays
                  for ( int i = 0; i < tctxt->nprow(); i++ ) {
                    if ( i == prow ) {
                      int size = 0;
                      tctxt->irecv(1,1,&size,1,0,pcol);
                      if (size > 0)
                        tctxt->drecv(size,1,&wftmpr[0],1,0,pcol);

                      if (size >= 0) {
                        assert(size==2*ft.np012loc());
                        // copy to complex array
                        vector<complex<double> > wftmp(ft.np012loc());      

                        //for ( int i = 0; i < ft.np012loc(); i++ )
                        //  wftmp[i] = complex<double>(wftmpr[2*i],wftmpr[2*i+1]);
                        //ewd when k=0, force complex part to zero
                        if (basis.real()) 
                          for ( int i = 0; i < ft.np012loc(); i++ )
                            wftmp[i] = complex<double>(wftmpr[2*i],0.0);
                        else 
                          for ( int i = 0; i < ft.np012loc(); i++ )
                            wftmp[i] = complex<double>(wftmpr[2*i],wftmpr[2*i+1]);
                        ComplexMatrix& c = sd_[ispin][kp]->c();
                        ft.forward(&wftmp[0],c.valptr(mloc*n));
                      }
                    }
                  }

                  if (tctxt->myrow() == 0) 
                    is.close();

                } // for n < nstloc

                
                // if there are empty states, read occupation and eigenvalues
                if (nempty_ > 0) {
                  ifstream is;
                  int nst = sd_[ispin][kp]->nst();
                  vector<double> eigtmp(nst);
                  vector<double> occtmp(nst);

                  if (tctxt->myrow() == 0) {
                    // read occupation for this state and k-point
                    ostringstream oss2,oss3;
                    oss2.width(4);  oss2.fill('0');  oss2 << kp;
                    oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                    string statefile;
                    if (nspin_ == 1) 
                      statefile = filebase + "k" + oss2.str() + ".occ"; 
                    else
                      statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + ".occ";
                    is.open(statefile.c_str(),ofstream::binary);

                    if (is.is_open()) {
                      // hack to make checkpointing work with BlueGene compilers
#ifdef BGQ
                      int len = statefile.length();
                      char* tmpfilename = new char[256];
                      is.read(tmpfilename,sizeof(char)*statefile.length());
#endif

                      is.read((char*)&eigtmp[0],sizeof(double)*nst);
                      is.read((char*)&occtmp[0],sizeof(double)*nst);
                      is.close();

                      for ( int i = 0; i < tctxt->nprow(); i++ ) {
                        tctxt->isend(1,1,&nst,1,i,pcol);
                        tctxt->dsend(nst,1,&occtmp[0],1,i,pcol);
                        tctxt->dsend(nst,1,&eigtmp[0],1,i,pcol);
                      }
                    }
                    else { // file not found
                      if ( ctxt_.oncoutpe())
                        cout << "<!-- LoadCmd:  occupation file not found. -->" << endl;
                      for ( int i = 0; i < tctxt->nprow(); i++ ) {
                        int size = -1;
                        tctxt->isend(1,1,&size,1,i,pcol);
                      }
                    }
                  }

                  for ( int i = 0; i < tctxt->nprow(); i++ ) {
                    if ( i == prow ) {
                      int size = 0;
                      tctxt->irecv(1,1,&size,1,0,pcol);

                      if (size > 0) {
                        assert(size==nst);
                        tctxt->drecv(size,1,&occtmp[0],1,0,pcol);
                        sd_[ispin][kp]->set_occ(occtmp);
                        tctxt->drecv(size,1,&eigtmp[0],1,0,pcol);
                        sd_[ispin][kp]->set_eig(eigtmp);
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
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write_mditer(string filebase, int mditer) {
  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  // write out mditer
  if (mype == 0)
  {
     string mditerfile = filebase + ".mditer";
     cout << "<!-- Wavefunction::write_states, writing mditer " << mditer << " to file " << mditerfile << " -->" << endl;
     ofstream osmd;
     osmd.open(mditerfile.c_str());
     osmd << mditer << endl;
     osmd.flush();
     osmd.close();
  }
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::read_mditer(string filebase, int& mditer) {
  int mype;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
#else
  mype = 0;
#endif

  // check for mditer file
  int mdtmp = -1;
  if (mype == 0)
  {
     string mditerfile = filebase + ".mditer";
     ifstream ismd;
     ismd.open(mditerfile.c_str());
     if (ismd.is_open()) {
        //ismd.read((char*)&mdtmp,sizeof(int));
        ismd >> mdtmp;
        ismd.close();
     }
  }
#if USE_MPI
  MPI_Bcast(&mdtmp, 1, MPI_INT, 0, MPI_COMM_WORLD);     
#endif
  
  if (mdtmp > 0 && mdtmp < 999999999)
     mditer = mdtmp;
}
  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::info(ostream& os, string tag)
{
  if ( ctxt_.oncoutpe() )
  {
    os << "<" << tag << " ecut=\"" << ecut_ << "\""
       << " nspin=\"" << nspin_ << "\""
       << " nel=\"" << nel_ << "\""
       << " nempty=\"" << nempty_ << "\">" << endl;
    os.setf(ios::fixed,ios::floatfield);
    os << "<cell a0=\""
       << setprecision(6) << cell_.a(0) << "\"\n      a1=\""
       << cell_.a(1) << "\"\n      a2=\""
       << cell_.a(2) << "\"/>" << endl;
    os << "<reciprocal_lattice b0=\""
       << setprecision(6) << cell_.b(0) << "\"\n      b1=\""
       << cell_.b(1) << "\"\n      b2=\""
       << cell_.b(2) << "\"/>" << endl;
    os << "<refcell a0=\""
       << refcell_.a(0) << "\"\n         a1=\""
       << refcell_.a(1) << "\"\n         a2=\""
       << refcell_.a(2) << "\"/>" << endl;
  }

  //ewd:  don't print sd->info for all k-points, as it can make the output file huge 
  // and the dimensions are usually the same for all k-points.  Instead, just print
  // info for first k-point
  /*
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
    for ( int ikp=0; ikp<nkp(); ikp++) {
    if (kptactive(ikp)) {
    assert(sd_[ispin][ikp] != 0);
        sd_[ispin][ikp]->info(os);
      }
    }
    }
  }
  */

  if ( ctxt_.myproc() == 0)
    sd_[0][0]->info(os);
  
  if ( ctxt_.oncoutpe() )
    os << "</" << tag << ">" << endl;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print_casino(ostream& os, int ispin, int kk) const {

   bool onproc0 = false;
   bool sdctxt_active = false;
   int kpc = -1;
   if (spinactive(ispin)) {

      kpc = mysdctxt(kk); // index of sdcontext for this kpoint   
      if (kpc >= 0)
         if (sdcontext_[ispin][kpc] != 0 ) 
            if (sdcontext_[ispin][kpc]->active() )
            {
               sdctxt_active = true;
               if (sdcontext_[ispin][kpc]->myproc() == 0 )
                  onproc0 = true;
            }

      if (sdctxt_active)
         if (sdcontext_[ispin][kpc]->mycol() == 0)
            sd_[ispin][kk]->basis().print_casino(os);
   
   }
   // pes not involved in basis print need to wait before sending data to proc 0
   if (sdctxt_active)
      sdcontext_[ispin][kpc]->barrier();
   
   if ( onproc0 ) {
      os << "WAVE FUNCTION" << endl;
      os << "-------------" << endl;
      os << "Number of k-points" << endl;
      os << nkp() << endl;
   }

   // wrapper loop to limit data going to proc 0 and filling up MPI buffers
   for (int nn=0; nn < nst_[ispin]; nn++)
   {
      if (sdctxt_active)
      {
         assert(sd_[ispin][kk] != 0);
         const vector<double>& eig = sd(ispin,kk)->eig();
         const int nloc = sd(ispin,kk)->c().nloc();
         for (int nl=0; nl < nloc; nl++) {
            const int nglobal = sd(ispin,kk)->c().j(0,nl);
            if (nn == nglobal) {
               // send all data to proc 0 for printing
               if (sdcontext_[ispin][kpc]->myrow() == 0)
               {
                  double eig_send = eig[nn];
                  sdcontext_[ispin][kpc]->dsend(1,1,&eig_send,1,0,0);
               }
               for ( int i = 0; i < sdcontext_[ispin][kpc]->nprow(); i++ )
               {
                  if ( i == sdcontext_[ispin][kpc]->myrow() )
                  {
                     int mloc = sd_[ispin][kk]->c().mloc();
                     int nloc = sd_[ispin][kk]->c().nloc();
                     int ngwloc = sd_[ispin][kk]->basis().localsize();
                     double* p = (double*) sd_[ispin][kk]->c().cvalptr();
                     int size = 2*ngwloc;
                     sdcontext_[ispin][kpc]->isend(1,1,&ngwloc,1,0,0);
                     sdcontext_[ispin][kpc]->dsend(size,1,&p[2*nl*mloc],1,0,0);
                  }
               }
            }
         }
      }
         
      if (onproc0)
      {
         if (nn == 0) {
            os << "k-point # ; # of bands (up spin/down spin) ; k-point coords (au)" << endl;
            os << " " << kk+1 << " " << nst_[ispin] << " " << ispin << " " << kpoint_[kk] << endl;
         }
         bool kk_real = ( kpoint_[kk] == D3vector(0.0,0.0,0.0) );
            
         // which process column is state nn data coming from?
         int nb = sd_[ispin][kk]->c().nb();
         int npcol = sdcontext_[ispin][kpc]->npcol();
         int colind = nn/nb;
         int col_recv = colind%npcol;
            
         double eig_recv;
         sdcontext_[ispin][kpc]->drecv(1,1,&eig_recv,1,0,col_recv);
         os << "Band, spin, eigenvalue (au)" << endl;
         os << "   " << nn+1 << " " << ispin+1 << " " << eig_recv << endl;
         os << "Eigenvector coefficients" << endl;
            
         for ( int i = 0; i < sdcontext_[ispin][kpc]->nprow(); i++ ) {
            int size = 0;
            sdcontext_[ispin][kpc]->irecv(1,1,&size,1,i,col_recv);
            int size_recv = 2*size;
            vector<double> readtmp(size_recv);
            sdcontext_[ispin][kpc]->drecv(size_recv,1,&readtmp[0],1,i,col_recv);
               
            int ngwread = size;
            vector<complex<double> > ctmp(ngwread);
            for (int j=0; j<ngwread; j++) 
               ctmp[j] = complex<double>(readtmp[2*j],readtmp[2*j+1]);
            if (kk_real) {
               // print out G-vectors for full complex basis
               if (i == 0) {
                  os << ctmp[0] << endl;
                  for (int j=1; j<ngwread; j++) {
                     os << ctmp[j] << endl;
                     os << conj(ctmp[j]) << endl;
                  }
               }
               else {
                  for (int j=0; j<ngwread; j++) {
                     os << ctmp[j] << endl;
                     os << conj(ctmp[j]) << endl;
                  }
               }
            }
            else {
               for (int j=0; j<ngwread; j++) {
                  os << ctmp[j] << endl;
               }
            }
         }
      }
      if (sdctxt_active)
         sdcontext_[ispin][kpc]->barrier();
   }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print_vmd(string filebase, const AtomSet& as) const {
   int maxnst = nst_[0];
   if (nspin_ > 1)
      if (nst_[1] > maxnst) maxnst = nst_[1];
   
   for (int n=0; n<maxnst; n++)
      print_vmd(filebase, as, n);
   
}
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print_vmd(string filebase, const AtomSet& as, const int statenum) const {
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {

            Context* tctxt = sdcontext_[ispin][ikp];
            int prow = tctxt->myrow();
            int pcol = tctxt->mycol();
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];         // global index of local kpoint
              if ( sd_[ispin][kp] != 0 ) {
                int nstloc = sd_[ispin][kp]->nstloc();
                int nb = sd_[ispin][kp]->c().nb();

                const Basis& basis = sd_[ispin][kp]->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
                
                D3vector a0 = cell_.a(0);
                D3vector a1 = cell_.a(1);
                D3vector a2 = cell_.a(2);
                const int np0 = ft.np0();
                const int np1 = ft.np1();
                const int np2 = ft.np2();
                D3vector dft0 = a0/(double)np0;
                D3vector dft1 = a1/(double)np1;
                D3vector dft2 = a2/(double)np2;

                for ( int n = 0; n < nstloc; n++ ) {
                  // global n index
                  const int nn = pcol*nb + n;
                  if (nn == statenum || statenum < 0)  // statenum < 0 means print all states
                  {
                  
                     ofstream os;
                     os.setf(ios::scientific,ios::floatfield);
                     os << setprecision(8);
                     if (tctxt->myrow() == 0) {
                        // write out wavefunction for this state and k-point
                        ostringstream oss1,oss2,oss3;
                        oss1.width(5);  oss1.fill('0');  oss1 << nn;
                        oss2.width(4);  oss2.fill('0');  oss2 << kp;
                        oss3.width(1);  oss3.fill('0');  oss3 << ispin;
                        string statefile;
                        if (nspin_ == 1) 
                           statefile = filebase + "k" + oss2.str() + "n" + oss1.str(); 
                        else
                           statefile = filebase + "s" + oss3.str() + "k" + oss2.str() + "n" + oss1.str(); 
                        statefile = statefile + ".cube";
                        os.open(statefile.c_str(),ofstream::out);
                        
                        // write out VMD CUBE format header
                        os << "Qbox wavefunction in VMD CUBE format" << endl;
                        os << "  state " << nn << ", k-point " << kp << " = " << kpoint_[kp] << endl;

                        // get atom positions
                        vector<vector<double> > rion;
                        rion.resize(as.nsp());
                        int natoms_total = 0;
                        for ( int is = 0; is < as.nsp(); is++ ) {
                           rion[is].resize(3*as.na(is));
                           natoms_total += as.na(is);
                        }
                        as.get_positions(rion,true);
                        D3vector origin(0.0,0.0,0.0);
                        os << natoms_total << " " << origin << endl;
                        
                        // print FFT grid info
                        os << np0 << " " << dft0 << endl;
                        os << np1 << " " << dft1 << endl;
                        os << np2 << " " << dft2 << endl;

                        // print atom coordinates
                        for ( int is = 0; is < as.nsp(); is++ ) {
                           const int atnum = as.atomic_number(is);
                           double atnumd = (double)atnum;
                           for ( int ia = 0; ia < as.na(is); ia++ ) 
                              os << atnum << " " << atnumd << " " << rion[is][3*ia] << " " << rion[is][3*ia+1] << " " << rion[is][3*ia+2] << endl;
                        }
                     }
                
                     // print isosurface:  values in six columns with z fast
                     int mloc = sd_[ispin][kp]->c().mloc();
                     vector<complex<double> > wftmp(ft.np012loc());
                     vector<double> wftmpr(2*ft.np012loc());
                  
                     ComplexMatrix& c = sd_[ispin][kp]->c();
                     ft.backward(c.cvalptr(mloc*n),&wftmp[0]);
                  
                     // copy |wf|^2 to double array for communication
                     double *a = (double*) &wftmp[0];
                     for ( int i = 0; i < ft.np012loc(); i++ )
                        wftmpr[i] = a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1];
                
                     // send data to first proc in context column
                     for ( int i = 0; i < tctxt->nprow(); i++ ) {
                        if ( i == prow ) {
                           int size = ft.np012loc();
                           tctxt->isend(1,1,&size,1,0,pcol);
                           tctxt->dsend(size,1,&wftmpr[0],1,0,pcol);
                        }
                     }
                     // receive data, store for reordering on output
                     if (tctxt->myrow() == 0) {
                        vector<double> wftmprecv(ft.np012());
                        int recvoffset = 0;
                        for ( int i = 0; i < tctxt->nprow(); i++ ) {
                           int size = 0;
                           tctxt->irecv(1,1,&size,1,i,pcol);
                           tctxt->drecv(size,1,&wftmprecv[recvoffset],1,i,pcol);
                           recvoffset += size;
                        }
                        
                        // write wf data to file
                        int cnt = 0;
                        for (int ii = 0; ii < np0; ii++) {
                           ostringstream oss;
                           oss.setf(ios::scientific,ios::floatfield);
                           oss << setprecision(5);
                           for (int jj = 0; jj < np1; jj++) {
                              for (int kk = 0; kk < np2; kk++) {
                                 int index = ii + jj*np0 + kk*np0*np1;
                                 oss << wftmprecv[index] << " ";
                                 cnt++;
                                 if (cnt >= 6) {
                                    cnt = 0;
                                    oss << endl;
                                 }
                              }
                           }
                           string tos = oss.str();
                           //os << tos.c_str();
                           os.write(tos.c_str(),tos.length());
                        }
                        os.close();

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
}

////////////////////////////////////////////////////////////////////////////////
// AS: is true when wave function has to be forced to complex also for kpoint == (0,0,0)
bool Wavefunction::force_complex_set(void) const { return force_complex_wf_; }

////////////////////////////////////////////////////////////////////////////////
// AS: enable or disable forcing of complex wave functions
void Wavefunction::force_complex(bool new_force_complex_wf)
{
   //ewd
   //deallocate();

  force_complex_wf_ = new_force_complex_wf;

  /* ewd
  kpoint_.resize(1);
  kpoint_[0] = D3vector(0,0,0);
  weight_.resize(1);
  weight_[0] = 1.0;
  compute_nst();
  //create_contexts();
  allocate();
  */  
}

////////////////////////////////////////////////////////////////////////////////
// AS: is true when the wave function is made real for Gamma only
bool Wavefunction::phase_real_set(void) const { return wf_phase_real_; }

////////////////////////////////////////////////////////////////////////////////
// AS: change phase of the wave function to make it real for Gamma only
void Wavefunction::phase_real(bool new_wf_phase_real)
{
  wf_phase_real_ = new_wf_phase_real;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, Wavefunction& wf)
{
  wf.print(os,"text","wavefunction");
  return os;
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction& Wavefunction::operator=(const Wavefunction& wf)
{
  if ( this == &wf ) return *this;
  assert(ctxt_ == wf.ctxt_);
  assert(nel_ == wf.nel_);
  assert(nempty_== wf.nempty_);
  assert(nspin_ == wf.nspin_);
  assert(nrowmax_ == wf.nrowmax_);
  assert(nparallelkpts_ == wf.nparallelkpts_);
  assert(deltaspin_ == wf.deltaspin_);
  assert(cell_ == wf.cell_);
  assert(refcell_ == wf.refcell_);
  assert(ecut_ == wf.ecut_);
  
  weight_ = wf.weight_;
  kpoint_ = wf.kpoint_;
  ultrasoft_ = wf.ultrasoft_;
  force_complex_wf_ = wf.force_complex_wf_;
  wf_phase_real_ = wf.wf_phase_real_;
  
  nkptloc_ = wf.nkptloc_;
  kptloc_ = wf.kptloc_;
  mysdctxt_ = wf.mysdctxt_;
  for ( int ispin = 0; ispin < nspin_; ispin++ ) {
    if (spinactive(ispin)) {
      for ( int ikp = 0; ikp < sdcontext_[ispin].size(); ikp++ ) {
        if (sdcontext_[ispin][ikp] != 0 ) {
          if (sdcontext_[ispin][ikp]->active() ) {
            for ( int kloc=0; kloc<nkptloc_; kloc++) {
              int kp = kptloc_[kloc];
              if ( sd_[ispin][kp] != 0 ) 
                *sd_[ispin][kp] = *wf.sd_[ispin][kp];
            }
          }
        }
      }
    }
  }
  return *this;
}
