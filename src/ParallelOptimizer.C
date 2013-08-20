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
///////////////////////////////////////////////////////////////////////////////
//
// ParallelOptimizer.C
//
////////////////////////////////////////////////////////////////////////////////

#include "ParallelOptimizer.h"
#include "EnergyFunctional.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "PSDAWavefunctionStepper.h"
#include "JDWavefunctionStepper.h"
#include "Preconditioner.h"

#include <iostream>
#include <iomanip>
using namespace std;

#ifdef USE_APC
#include "apc.h"
#endif

////////////////////////////////////////////////////////////////////////////////
ParallelOptimizer::ParallelOptimizer(Sample& s) : s_(s) {

  fion.resize(s_.atoms.nsp()+s_.atoms.nsp_mm());
  for ( int is = 0; is < fion.size(); is++ )
    fion[is].resize(3*s_.atoms.na(is));
    
  sigma_eks.resize(6);
  sigma_kin.resize(6);
  sigma_ext.resize(6);
  sigma.resize(6);
}

////////////////////////////////////////////////////////////////////////////////
ParallelOptimizer::~ParallelOptimizer() {
}

////////////////////////////////////////////////////////////////////////////////
void ParallelOptimizer::optimize(int niter, int nitscf, int nite) {

   if ( s_.ctxt_.oncoutpe() ) 
      cout << "  <!-- ParallelOptimizer started. -->" << endl;
   
   const int maxnrowmax_ = 4096;

  niter_ = niter;
  nitscf_ = nitscf;
  nite_ = nite;

  int npes;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
#else
  npes = 1;
#endif
   if ( s_.ctxt_.oncoutpe() ) 
      cout << "  <!-- ParallelOptimizer called with npes = " << npes << ". -->" << endl;

  Timer tm_partot;
  tm_partot.start();

  int nkpts = s_.wf.nkp();
  int nspin = s_.wf.nspin();

  // factorize npes to get set of candidate parameters for nrowmax, nparallelkpts
  list<int> npe_factors = factorize(npes);
  list<int> npkp_factors = factorize(npes,nkpts); // will return only 1 if nkpts==1

  int nparkp = npkp_factors.back();
  int nrowmax = npes/(nspin*nparkp);

  // make context rectangular for initial guess
  if (nrowmax%2 == 0 && nrowmax > 4) 
    nrowmax /= 2;   
  if (nrowmax%2 == 0 && nrowmax > 16) 
    nrowmax /= 2;   
  //ewd: just added
  if (nrowmax%2 == 0 && nrowmax > 32)
    nrowmax /= 2;   

  // if nkpts > 1, optimize nparallelkpts
  double kpttime = 1.E+20;
  if (nkpts > 1) {     
    kpttime = runtime(nrowmax,nparkp,nspin,s_.ctrl.reshape_context,true);
    if ( s_.ctxt_.oncoutpe() ) 
      cout << "  <!-- ParallelOptimizer: nparallelkpts = " << nparkp << ", nrowmax = " << nrowmax << ", runtime = " << kpttime << " -->" << endl;

    if (npkp_factors.size() > 1) { // choose nparallelkpts from common factors 
      bool improve = true;
      list<int>::reverse_iterator lit = npkp_factors.rbegin();
      assert(*lit == nparkp);
      lit++;
      while (improve && lit != npkp_factors.rend()) { 
        int testnpkp = *lit;
        int npcol = npes/(nspin*testnpkp*nrowmax);
        int nprow = nrowmax;
        if (npcol > nrowmax) {
          nprow = npcol;
        }
        double testtime = runtime(nrowmax,testnpkp,nspin,s_.ctrl.reshape_context,true);
        if ( s_.ctxt_.oncoutpe() ) 
          cout << "  <!-- ParallelOptimizer: nparallelkpts = " << testnpkp << ", nrowmax = " << nprow    << ", runtime = " << testtime << " -->" << endl;
        if (testtime < kpttime) {
          kpttime = testtime;
          nparkp = testnpkp;
          nrowmax = nprow;
        }
        else {
          improve = false;
        }
        lit++;
      }
    }
    else {  // nkpts does not have common factors with npes, sweep through npe_factors
      bool improve = true;
      bool stophere = false;
      list<int>::iterator lit = npe_factors.begin();
      lit++;
      while (!stophere && lit != npe_factors.end() && *lit <= nkpts) {
        int testnpkp = *lit;
        int nprow = nrowmax;
        if (nprow > npes/(nspin*testnpkp)) 
          nprow = npes/(nspin*testnpkp);
        int npcol = npes/(nspin*testnpkp*nrowmax);
        if (npcol > nrowmax)
          nprow = npcol;
        double testtime = runtime(nprow,testnpkp,nspin,s_.ctrl.reshape_context,true);
        if ( s_.ctxt_.oncoutpe() ) 
          cout << "  <!-- ParallelOptimizer: nparallelkpts = " << testnpkp << ", nrowmax = " << nprow << ", runtime = " << testtime << " -->" << endl;
        if (testtime < kpttime) {
          kpttime = testtime;
          nparkp = testnpkp;
          nrowmax = nprow;
          improve = true;
        }
        else {
          if (improve)
            improve = false;
          else 
            stophere = true; // stop after two consecutive values don't improve timing
        }
        lit++;
      }
    }
  } // if nkpts > 1

  // once nparallelkpts is chosen, optimize nrowmax 
  s_.wf.set_nparallelkpts(nparkp);
  int npesleft = npes/nparkp;
  list<int> nrow_factors = factorize(npesleft);

  double nrowtime = kpttime;  // best timing so far
  bool improve = true;
  bool reshape = false;
  list<int>::reverse_iterator lit = nrow_factors.rbegin();
  while (improve && lit != nrow_factors.rend()) { 
    int testnrow = *lit;
    int npcol = npesleft/testnrow;
    if (npcol > testnrow) {
      improve = false;
    }
    else if (testnrow > maxnrowmax_) {
      // do nothing, try next nrowmax
    }
    else {
      bool testreshape = false;
      double testtime = runtime(testnrow,1,nspin,false,true);
      if ( s_.ctxt_.oncoutpe() ) 
        cout << "  <!-- ParallelOptimizer: nparallelkpts = " << nparkp << ", nrowmax = " << testnrow << ", reshape = false, runtime = " << testtime << " -->" << endl;

      // try turning on context reshaping, see if time improves
      //ewd:  disable context reshaping for now
      if (false && testnrow > 2*npcol) {
        double reshapetime = runtime(testnrow,1,nspin,true,true);
        if ( s_.ctxt_.oncoutpe() ) 
          cout << "  <!-- ParallelOptimizer: nparallelkpts = " << nparkp << ", nrowmax = " << testnrow << ", reshape = true, runtime = " << reshapetime << " -->" << endl;
      
        if (reshapetime < testtime)
        {
           testtime = reshapetime;
           testreshape = true;
        }
      }
      
      if (testtime < nrowtime) {
        if (testreshape)
          reshape = true;
        nrowtime = testtime;
        nrowmax = testnrow;
      }
      else {
        improve = false;
      }
    }
    lit++;
  }

  string reshapetxt = "OFF";
  if (reshape)
    reshapetxt = "ON";
  if ( s_.ctxt_.oncoutpe() ) 
    cout << "<!-- ParallelOptimizer:  best runtime = " << nrowtime << ", nrowmax = " << nrowmax << ", nparallelkpts = " << nparkp << ", reshape_context = " << reshapetxt << " -->" << endl;
  s_.wf.set_nrowmax(nrowmax);
  s_.wf.set_nparallelkpts(nparkp);

  // randomize wavefunction
  bool highmem = false;
  if (s_.ctrl.extra_memory >= 3)
    highmem = true;
  if (s_.ctrl.ultrasoft)
     s_.wf.randomize_us(0.02,s_.atoms,highmem);     // this will force allocate() and resize()
  else
     s_.wf.randomize(0.02,highmem);              // this will force allocate() and resize()
  s_.ctrl.reshape_context = reshape;
  s_.wf.set_reshape_context(reshape);


  tm_partot.stop();
  // print iteration time
  double time = tm_partot.real();
  double tmin = time;
  double tmax = time;
  s_.ctxt_.dmin(1,1,&tmin,1);
  s_.ctxt_.dmax(1,1,&tmax,1);
  if ( s_.ctxt_.oncoutpe() ) {
    cout << "  <!-- timing "
         << setw(15) << "paropt_total"
         << " : " << setprecision(3) << setw(9) << tmin
         << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
  }

}
////////////////////////////////////////////////////////////////////////////////
double ParallelOptimizer::runtime(int nrowmax, int npark, int nspin, bool reshape, bool print_timing) {

   if (s_.ctxt_.oncoutpe())
      cout << "<!--ParallelOptimizer::runtime called w. nparallelkpts = " << npark << ", nrowmax = " << nrowmax << ", nspin = " << nspin << ", reshape = " << reshape << "-->" << endl;

   const int nempty = s_.wf.nempty();
   const bool compute_eigvec = nempty > 0 || s_.ctrl.wf_diag == "T";

  AtomSet& atoms = s_.atoms;
  const UnitCell& cell = s_.wf.cell();
  const UnitCell& refcell = s_.wf.refcell();
  const double omega = cell.volume();
  const double dt = s_.ctrl.dt;

  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;
  
  const bool atoms_move = ( niter_ > 0 && atoms_dyn != "LOCKED" );
  const bool compute_hpsi = ( wf_dyn != "LOCKED" );
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool use_confinement = ( s_.ctrl.ecuts > 0.0 );
  const bool use_preconditioner = wf_dyn == "PSD" || wf_dyn == "PSDA";  
  const bool ultrasoft = s_.ctrl.ultrasoft;
  const bool usdiag = (ultrasoft && atoms_move);
  const bool nlcc = s_.ctrl.nlcc;

  const string charge_mixing = s_.ctrl.charge_mixing;

  //Wavefunction wf(s_.ctxt_);
  Wavefunction& wf = s_.wf;
  wf.set_nrowmax(nrowmax);
  wf.set_nparallelkpts(npark);
  bool highmem = false;
  if (s_.ctrl.extra_memory >= 3)
    highmem = true;
  if (ultrasoft)
    wf.set_ultrasoft(ultrasoft);

  if (ultrasoft)
     wf.randomize_us(0.02,atoms,highmem);  // this will force allocate() and resize()
  else
     wf.randomize(0.02,highmem);              // this will force allocate() and resize()
  wf.set_reshape_context(reshape);
  
  if ( nempty > 0 )
    wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);

  // simulation objects for this parallel distribution
  ChargeDensity cd_(s_);
  EnergyFunctional ef_(s_,wf,cd_);
  Wavefunction dwf(wf);
  dwf.set_reshape_context(reshape);
  cd_.set_nlcc(nlcc);

  // use extra memory for SlaterDets if memory variable = normal, large or huge
  if (s_.ctrl.extra_memory >= 3) 
    wf.set_highmem();

  Preconditioner *preconditioner = 0;
  if ( use_preconditioner ) {
    preconditioner = new Preconditioner(s_,wf,ef_);
  }
  
  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "SD" )
    wf_stepper = new SDWavefunctionStepper(wf,1.0,tmap);
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(wf,*preconditioner,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(wf,*preconditioner,tmap);  
  else if ( wf_dyn == "JD" )
    wf_stepper = new JDWavefunctionStepper(wf,*preconditioner,ef_,tmap);  

  // if ultrasoft, calculate position-dependent functions
  double time_cd_usfns = 0.0;
  if (ultrasoft) {
    Timer tm_cdus;
    tm_cdus.start();
    cd_.update_usfns();
    tm_cdus.stop();
    time_cd_usfns = tm_cdus.real();
    s_.ctxt_.dmax(1,1,&time_cd_usfns,1);
  }

  // if non-linear core correction defined, calculate position-dependent density
  if (nlcc)
     cd_.update_nlcc();
    
  double time_cd_update = 0.0;
  {
    Timer tm_cd;
    tm_cd.start();
    cd_.update_density();
    tm_cd.stop();
    time_cd_update = tm_cd.real();
    s_.ctxt_.dmax(1,1,&time_cd_update,1);
  }
  
  double time_cd_rhor = 0.0;
  {
    Timer tm_cd;
    tm_cd.start();
    cd_.update_rhor();
    tm_cd.stop();
    time_cd_rhor = tm_cd.real();
    s_.ctxt_.dmax(1,1,&time_cd_rhor,1);
  }
  
  double time_ef_nonscf = 0.0;
  {
    Timer tm_ef;
    tm_ef.start();
    double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
    tm_ef.stop();
    time_ef_nonscf = tm_ef.real();
    s_.ctxt_.dmax(1,1,&time_ef_nonscf,1);
  }

  double time_ef_ionic = 0.0;
  if ( compute_forces || compute_stress ) {
    Timer tm_ef;
    tm_ef.start();
    double energy = ef_.energy(false,dwf,compute_forces,fion,
                               compute_stress,sigma_eks);
    tm_ef.stop();
    time_ef_ionic = tm_ef.real();
    s_.ctxt_.dmax(1,1,&time_ef_ionic,1);
  }

  double time_ef_vhxc = 0.0;
  {
    Timer tm_ef;
    tm_ef.start();
    ef_.update_vhxc();
    tm_ef.stop();
    time_ef_vhxc = tm_ef.real();
    s_.ctxt_.dmax(1,1,&time_ef_vhxc,1);
  }

  double time_sd_usfns = 0.0;
  if (ultrasoft) {
    D3vector gamma(0.0,0.0,0.0);
    int nkptloc = wf.nkptloc();
    int kp = wf.kptloc(0);         // global index of local kpoint
    if (nkptloc > 1 && wf.kpoint(kp) == gamma) 
      kp = wf.kptloc(1);  // use more expensive non-gamma calculation for timing
    Timer tm_sdus;
    tm_sdus.start();
    wf.sd(0,kp)->update_usfns();
    tm_sdus.stop();
    time_sd_usfns = nkptloc*tm_sdus.real();
    s_.ctxt_.dmax(1,1,&time_sd_usfns,1);
  }

  double time_wf_gram = 0.0;
  {
    D3vector gamma(0.0,0.0,0.0);
    int nkptloc = wf.nkptloc();
    int kp = wf.kptloc(0);         // global index of local kpoint
    if (nkptloc > 1 && wf.kpoint(kp) == gamma) 
      kp = wf.kptloc(1);  // use more expensive non-gamma calculation for timing
    Timer tm_gram;
    tm_gram.start();
    wf.sd(0,kp)->gram();
    tm_gram.stop();
    time_sd_usfns = nkptloc*tm_gram.real();
    s_.ctxt_.dmax(1,1,&time_wf_gram,1);
  }

  double time_sd_betapsi = 0.0;
  if (ultrasoft) {
    D3vector gamma(0.0,0.0,0.0);
    int nkptloc = wf.nkptloc();
    int kp = wf.kptloc(0);         // global index of local kpoint
    if (nkptloc > 1 && wf.kpoint(kp) == gamma) 
      kp = wf.kptloc(1);  // use more expensive non-gamma calculation for timing
    Timer tm_bpsi;
    tm_bpsi.start();
    wf.sd(0,kp)->calc_betapsi();
    tm_bpsi.stop();
    time_sd_betapsi = nkptloc*tm_bpsi.real();
    s_.ctxt_.dmax(1,1,&time_sd_betapsi,1);
  }

  double time_wf_dot = 0.0;
  {
    Timer tm_wf;
    tm_wf.start();
    double eigsum = wf.dot(dwf);
    tm_wf.stop();
    time_wf_dot = tm_wf.real();
    s_.ctxt_.dmax(1,1,&time_wf_dot,1);
  }

  double time_wf_align = 0.0;
  if ( compute_forces && compute_eigvec ) {
    Timer tm_wf;
    tm_wf.start();
    dwf.align(wf);
    tm_wf.stop();
    time_wf_align = tm_wf.real();
    s_.ctxt_.dmax(1,1,&time_wf_align,1);
  }

  double time_wf_ortho_align = 0.0;
  if ( compute_forces) {
    D3vector gamma(0.0,0.0,0.0);
    Timer tm_wf;
    tm_wf.start();
    int nkptloc = wf.nkptloc();
    int kp = wf.kptloc(0);         // global index of local kpoint
    if (nkptloc > 1 && wf.kpoint(kp) == gamma) 
      kp = wf.kptloc(1);  // use more expensive non-gamma calculation for timing
    dwf.sd(0,kp)->ortho_align(*dwf.sd(0,kp));
    tm_wf.stop();
    time_wf_ortho_align = nkptloc*tm_wf.real();
    s_.ctxt_.dmax(1,1,&time_wf_ortho_align,1);
  }

  double time_wf_diag = 0.0;
  if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" ) {
    Timer tm_wf;
    tm_wf.start();
    dwf.set_ultrasoft(s_.ctrl.ultrasoft);
    dwf.diag(dwf,compute_eigvec);
    tm_wf.stop();
    time_wf_diag = tm_wf.real();
    s_.ctxt_.dmax(1,1,&time_wf_diag,1);
  }

  double time_wf_update = 0.0;
  {
    Timer tm_wf;
    tm_wf.start();
    wf_stepper->update(dwf);
    tm_wf.stop();
    time_wf_update = tm_wf.real();
    s_.ctxt_.dmax(1,1,&time_wf_update,1);
  }


  if ( s_.ctxt_.oncoutpe() && print_timing ) {
    cout << "    <!-- ParallelOptimizer:  EnergyFunctional timing: ionic =  " << time_ef_ionic << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  EnergyFunctional timing: nonscf = " << time_ef_nonscf << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  EnergyFunctional timing: vhxc = " << time_ef_vhxc << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  ChargeDensity timing: update = " << time_cd_update << " -->" << endl;
    if (ultrasoft)
       cout << "    <!-- ParallelOptimizer:  ChargeDensity timing: usfns = " << time_cd_usfns << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  ChargeDensity timing: rhor = " << time_cd_rhor << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  Wavefunction timing: diag = " << time_wf_diag << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  Wavefunction timing: dot = " << time_wf_dot << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  Wavefunction timing: align = " << time_wf_align << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  Wavefunction timing: gram = " << time_wf_gram << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  Wavefunction timing: ortho_align = " << time_wf_ortho_align << " -->" << endl;
    cout << "    <!-- ParallelOptimizer:  Wavefunction timing: update = " << time_wf_update << " -->" << endl;
    if (ultrasoft) {
       cout << "    <!-- ParallelOptimizer:  SlaterDet timing: update_usfns = " << time_sd_usfns << " -->" << endl;
       cout << "    <!-- ParallelOptimizer:  SlaterDet timing: calc_betapsi = " << time_sd_betapsi << " -->" << endl;
    }       
  }


  double nonscf_tot = time_ef_nonscf + time_wf_dot + time_wf_update + time_sd_usfns + time_wf_gram + time_sd_betapsi;
  double scf_tot = time_cd_update + time_ef_vhxc + nite_*nonscf_tot;
  if (nite_ > 1) 
    scf_tot += time_cd_rhor;
  if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" ) 
    scf_tot += time_ef_nonscf + time_wf_diag;

  double ionic_tot = nitscf_*scf_tot + time_cd_update + time_ef_vhxc + time_ef_ionic + time_cd_usfns;
  if ( compute_forces && compute_eigvec )
    ionic_tot += time_wf_align;
  if ( compute_forces )
    ionic_tot += time_wf_ortho_align;

  // delete stepper
  if ( wf_stepper != 0 ) delete wf_stepper;
  
  // delete preconditioner
  if ( use_preconditioner ) delete preconditioner;

  int tmpiter = niter_ > 0 ? niter_ : 1;
  return tmpiter*ionic_tot;

}


////////////////////////////////////////////////////////////////////////////////
list<int> ParallelOptimizer::factorize(int n) {
  // find factors of n

  list<int> factors;
  int maxfacts = (int) sqrt((double)n);
  for (int i=1; i<=maxfacts; i++) {
    if (n%i == 0) {
      factors.push_back(i);
      int rem = n/i;
      if (rem != i) 
        factors.push_back(rem);
    }
  }
  factors.sort();
  return factors;
}

////////////////////////////////////////////////////////////////////////////////
list<int> ParallelOptimizer::factorize(int n1, int n2) {
  // find common factors of n1, n2

  list<int> factors;
  int n = n1 > n2 ? n1 : n2;
  int maxfacts = (int) sqrt((double)n);
  for (int i=1; i<=n; i++) 
    if (n1%i == 0 && n2%i == 0) 
      factors.push_back(i);

  factors.sort();
  factors.unique();

  return factors;
}

