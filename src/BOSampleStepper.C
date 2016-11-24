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
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "BOSampleStepper.h"
#include "EnergyFunctional.h"
#include "Extrapolator.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "JDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "PSDAWavefunctionStepper.h"
#include "SDIonicStepper.h"
#include "SDAIonicStepper.h"
#include "CGIonicStepper.h"
#include "MDIonicStepper.h"
#include "BMDIonicStepper.h"
#include "SDCellStepper.h"
#include "FCPStepper.h"
#include "Preconditioner.h"
#include "AndersonMixer.h"
#include "MLWFTransform.h"
#include "SimpleConvergenceDetector.h"
#include "Hugoniostat.h"
#include "PrintMem.h"
#include "FourierTransform.h"
#include "profile.h"
#include <fstream>
#include <sys/stat.h>
#ifdef USE_APC
#include "apc.h"
#endif
#include <iostream>
#include <iomanip>
#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
extern "C" void summary_start(void);
extern "C" void summary_stop(void);
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, int nitscf, int nite) :
  SampleStepper(s), cd_(s), ef_(s,s.wf,cd_),dwf(s.wf), wfv(s.wfv), nitscf_(nitscf),
  nite_(nite), initial_atomic_density(false) {}

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::~BOSampleStepper()
{
}
////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::initialize_density(void)
{
  // initialize cd_ with a sum of atomic densities

  double atom_radius[] =
      {
          0.00,  0.80,  0.59,  3.16,  2.12,  // null H  He  Li  Be
          1.64,  1.27,  1.06,  0.91,  0.79,  //  B   C   N   O   F
          0.72,  3.59,  2.74,  2.23,  2.10,  // Ne  Na  Mg  Al  Si
          1.85,  1.66,  1.49,  1.34,  4.59,  //  P   S  Cl  Ar   K
          3.67,  3.48,  3.33,  3.23,  3.14,  // Ca  Sc  Ti   V  Cr
          3.04,  2.95,  2.87,  2.82,  2.74,  // Mn  Fe  Co  Ni  Cu
          2.68,  2.57,  2.36,  2.15,  1.95,  // Zn  Ga  Ge  As  Se
          1.78,  1.66,  5.01,  4.14,  4.01,  // Br  Kr  Rb  Sr   Y
          3.89,  3.74,  3.59,  3.46,  3.36,  // Zr  Nb  Mo  Tc  Ru
          3.27,  3.19,  3.12,  3.04,  2.95,  // Rh  Pd  Ag  Cd  In
          2.74,  2.51,  2.32,  2.17,  2.04,  // Sn  Sb  Te   I  Xe
          5.63,  4.78,  0.00,  0.00,  4.67,  // Cs  Ba  La  Ce  Pr
          3.89,  3.87,  4.50,  4.37,  4.40,  // Nd  Pm  Sm  Eu  Gd
          4.25,  4.31,  0.00,  4.27,  4.20,  // Tb  Dy  Ho  Er  Tm
          4.20,  4.10,  3.93,  3.78,  3.65,  // Yb  Lu  Hf  Ta   W
          4.25,  4.31,  0.00,  4.27,  4.20,  // Tb  Dy  Ho  Er  Tm
          4.20,  4.10,  3.93,  3.78,  3.65,  // Yb  Lu  Hf  Ta   W
          3.55,  3.50,  3.40,  3.34,  3.29,  // Re  Os  Ir  Pt  Au
          3.23,  2.95,  2.91,  2.70,  2.55,  // Hg  Tl  Pb  Bi  Po
          4.00,  2.27,  4.00,  4.00,  4.00,  // At  Rn  Fr  Ra  Ac
          4.00,  4.00,  4.00,  4.00,  4.00,  // Th  Pa   U  Np  Pu
          4.00,  4.00,  4.00,  4.00,  4.00,  // Am  Cm  Bk  Cf  Es
          4.00,  4.00,  4.00,  4.00          // Fm  Md  No  Lr
      };
  

  if (!s_.wf.hasdata())
    s_.wf.set_hasdata(true);

  const AtomSet& atoms = s_.atoms;
  const Basis* const vbasis = cd_.vbasis();
  const int ngloc = vbasis->localsize();
  vector<vector<complex<double> > > rhops;
  const int nsp = atoms.nsp();
  rhops.resize(nsp);
  vector<complex<double> > rhopst(ngloc);
  const double * const g2 = vbasis->g2_ptr();
  
  for ( int is = 0; is < nsp; is++ )
  {
    rhops[is].resize(ngloc);
    Species *s = atoms.species_list[is];
    const int zval = s->zval();
    const int atomic_number = s->atomic_number();
    assert(atomic_number < sizeof(atom_radius)/sizeof(double));
    // scaling factor 2.0 in next line is empirically adjusted
    double rc = 2.0 * atom_radius[atomic_number];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      double arg = 0.25 * rc * rc * g2[ig];
      rhops[is][ig] = zval * exp( -arg );
    }
  }
  
  vector<vector<double> > tau0;
  tau0.resize(nsp);
  for ( int is = 0; is < nsp; is++ )
    tau0.resize(3*atoms.na(is));
  atoms.get_positions(tau0);
  StructureFactor sf;
  sf.init(tau0,*vbasis);
  sf.update(tau0,*vbasis);
  
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < nsp; is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const complex<double> sg = s[ig];
      rhopst[ig] += sg * rhops[is][ig];
    }
  }

  for (int ispin = 0; ispin < s_.wf.nspin(); ispin++) {
    //if (s_.wf.spinactive(ispin)) {
      cd_.rhog[ispin].resize(ngloc);
      for (int i = 0; i < ngloc; i++) {
        cd_.rhog[ispin][i] = rhopst[i];
      }
      //}
  }
  initial_atomic_density = true;
}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::step(int niter)
{
  const bool onpe0 = s_.ctxt_.onpe0();

  const bool anderson_charge_mixing = ( s_.ctrl.charge_mix_ndim > 0 );
  
  // determine whether eigenvectors must be computed
  // eigenvectors are computed if explicitly requested with wf_diag==T
  // or if the SlaterDet has fractionally occupied states
  //const bool fractional_occ = (s_.wf.nspin() * s_.wf.nel() != 2 * s_.wf.nst());
  int occtest = (2 * s_.wf.nst()) - s_.wf.nspin() * s_.wf.nel();
  const bool fractional_occ = (occtest != 0 && occtest != 1);
  const bool compute_eigvec = fractional_occ || s_.ctrl.wf_diag == "T";
  const bool compute_mlwf = s_.ctrl.wf_diag == "MLWF";
  const bool compute_mlwfc = s_.ctrl.wf_diag == "MLWFC";
  enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };

  if (fractional_occ && onpe0)
       cout << "<!-- BOSampleStepper:  fractional occupation detected. -->" << endl;
  else if (onpe0)
       cout << "<!-- BOSampleStepper:  fractional occupation not detected. -->" << endl;
   
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  const int nspin = wf.nspin();

  const UnitCell& cell = wf.cell();

  const double dt = s_.ctrl.dt;

  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;

  //ewd check for case where PSDA used (incorrectly) with nite = 1 and empty states
  if (wf_dyn == "PSDA" && nite_ == 0 && compute_eigvec) {
    if ( onpe0 ) {
      cout << "<ERROR> BOSampleStepper:  PSDA unstable with empty states and nite = 0. </ERROR>" << endl;
      cout << "<ERROR> BOSampleStepper:  Increase nite or use wf_dyn = PSD. </ERROR>" << endl;
    }
    return;
  }
  
  const bool extrapolate_wf = (atoms_dyn == "MD" && s_.ctrl.wf_extrap != "OFF");
  Extrapolator extrapolator;
  Wavefunction* wfmm;
  if ( extrapolate_wf && ( s_.ctrl.wf_extrap == "ASP" || s_.ctrl.wf_extrap == "NTC" ) ) 
    wfmm = new Wavefunction(wf);

  // Next lines: special value of niter = 0: GS calculation only
  const bool atoms_move = ( niter > 0 && atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool cell_moves = ( niter > 0 && compute_stress &&
                            cell_dyn != "LOCKED" );
  // GS-only calculation:
  const bool gs_only = !atoms_move && !cell_moves;
  const bool use_confinement = ( s_.ctrl.ecuts > 0.0 );

  const string charge_mixing = s_.ctrl.charge_mixing;

  const bool ultrasoft = s_.ctrl.ultrasoft;
  const bool usdiag = (ultrasoft && atoms_move);
  const bool nlcc = s_.ctrl.nlcc;
  cd_.set_nlcc(nlcc);
  
  //ewd check that MLWF not being used with ultrasoft (not yet implemented)
  if (ultrasoft && (compute_mlwf || compute_mlwfc)) {
    if ( onpe0 ) 
      cout << "<ERROR> BOSampleStepper:  Maximally-localized Wannier Functions not yet implemented with ultrasoft. </ERROR>" << endl;
    return;
  }
  if (onpe0 && ultrasoft && compute_stress)
      cout << "<WARNING> BOSampleStepper:  stress not yet implemented with ultrasoft!  Results WILL be wrong. </WARNING>" << endl;
  if (onpe0 && nlcc && compute_stress)
      cout << "<WARNING> BOSampleStepper:  stress not yet implemented with non-linear core corrections!  Results WILL be wrong. </WARNING>" << endl;
  
  // use extra memory for SlaterDets if memory variable = normal, large or huge
  if (s_.ctrl.extra_memory >= 3) 
    wf.set_highmem();
  
  if (!gs_only && s_.ctrl.dft_plus_u && onpe0)
    cout << "<WARNING> Forces and stress not currently implemented with DFT+U! </WARNING>" << endl;
  
  Timer tm_iter;

  const bool use_preconditioner = wf_dyn == "PSD" || wf_dyn == "PSDA" || wf_dyn == "JD";
  Preconditioner *preconditioner = 0;
  if ( use_preconditioner )
  {
    // create a preconditioner using the information about wf in s_.wf
    // and the information about the hessian in df
    preconditioner = new Preconditioner(s_,s_.wf,ef_);
  }

  // initialize occupation
  wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);

  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "SD" )
  {
    const double emass = s_.ctrl.emass;
    double dt2bye = (emass == 0.0) ? 0.5 / wf.ecut() : dt*dt/emass;

    // divide dt2bye by facs coefficient if stress == ON
    const double facs = 2.0;
    if ( s_.ctrl.stress == "ON" )
    {
      dt2bye /= facs;
    }
    wf_stepper = new SDWavefunctionStepper(wf,dt2bye,tmap);
  }
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(wf,*preconditioner,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(wf,*preconditioner,tmap);  
  else if ( wf_dyn == "JD" )
    wf_stepper = new JDWavefunctionStepper(wf,*preconditioner,ef_,tmap);  
  // wf_stepper == 0 indicates that wf_dyn == LOCKED

  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "SDA" )
    ionic_stepper = new SDAIonicStepper(s_);
  else if ( atoms_dyn == "CG" )
    ionic_stepper = new CGIonicStepper(s_);
  else if ( atoms_dyn == "MD" || atoms_dyn == "IMPULSIVE")
    ionic_stepper = new MDIonicStepper(s_);
  else if ( atoms_dyn == "BMD" )
    ionic_stepper = new BMDIonicStepper(s_);

  if ( ionic_stepper )
    ionic_stepper->setup_constraints();

  double fmu = 0.0;
  FCPStepper* fcp_stepper = 0;
  if ( s_.ctrl.fcp_thermostat == "SCALING" )
     fcp_stepper = new FCPStepper(s_);

  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);

  // Allocate wavefunction velocity if not available
  if ( atoms_move && extrapolate_wf )
  {
    if ( s_.wfv == 0 )
    {
      s_.wfv = new Wavefunction(wf);
      s_.wfv->clear();
    }
  }

  MLWFTransform* mlwft=0;

  if ( compute_mlwf || compute_mlwfc )
  {
    // MLWF can be computed at the gamma point only
    // There must be a single k-point, and it must be gamma
    if ( wf.nkp() > 1 || ( wf.nkp()==1 && wf.kpoint(0) != D3vector(0,0,0) ) )
    {
      if ( onpe0 )
      {
        cout << "<ERROR> BOSampleStepper::step: MLWF can be computed at k=0 only </ERROR>"
             << endl;
        cout << "<ERROR> BOSampleStepper::step: cannot run </ERROR>" << endl;
      }
      return;
    }
    if (wf.nspin() > 1 && s_.ctxt_.oncoutpe()) 
      cout << "<ERROR> nspin > 1!  MLWF doesn't currently work with spin-polarized systems </ERROR>" << endl;
    assert(wf.nspin()==1);
    mlwft = new MLWFTransform(*wf.sd(0,0));
  }

  //ewd:  experimental Hugoniostat for Kyle and Sebastien
  Hugoniostat* hugstat = 0;
  if (s_.ctrl.hugoniostat == "ON") {
    hugstat = new Hugoniostat(s_.ctrl.hug_etot,s_.ctrl.hug_volume,s_.ctrl.hug_pressure,ionic_stepper->temp(),s_.ctxt_.oncoutpe());
    hugstat->set_deltatemp(s_.ctrl.hug_deltatemp);
    hugstat->set_updatefreq(s_.ctrl.hug_freq);
  }

  // Charge mixing variables
  vector<vector<complex<double> > > rhog_in;
  vector<vector<complex<double> > > drhog;
  vector<vector<complex<double> > > rhobar;
  vector<vector<complex<double> > > drhobar;

  rhog_in.resize(nspin);
  drhog.resize(nspin);
  rhobar.resize(nspin);
  drhobar.resize(nspin);
  for (int ispin = 0; ispin < nspin; ispin++) {
    rhog_in[ispin].resize(cd_.rhog[ispin].size());
    drhog[ispin].resize(cd_.rhog[ispin].size());
    rhobar[ispin].resize(cd_.rhog[ispin].size());
    drhobar[ispin].resize(cd_.rhog[ispin].size());
  }

  vector<double> wkerker(rhog_in[0].size());
  vector<double> wls(rhog_in[0].size());

  const int anderson_ndim = s_.ctrl.charge_mix_ndim;

  vector<AndersonMixer*> mixer;
  mixer.resize(nspin);
  for (int ispin = 0; ispin < nspin; ispin++)
    //if (s_.wf.spinactive(ispin))
      mixer[ispin] = new AndersonMixer(2*rhog_in[ispin].size(),anderson_ndim,&cd_.vcontext());
      
  // compute Kerker preconditioning
  // real space Kerker cutoff in a.u.
  const double rc_Kerker = s_.ctrl.charge_mix_rcut;
  const double *const g2 = cd_.vbasis()->g2_ptr();

  // define q1 cutoff for row weighting of LS charge mixing
  // Use rc1 = 3 a.u. default cutoff
  double rc1 = 3.0;
  // check if override from the debug variable
  // use: set debug RC1 <value>
  if ( s_.ctrl.debug.find("RC1") != string::npos )
  {
    istringstream is(s_.ctrl.debug);
    string s;
    is >> s >> rc1;
    cout << " override rc1 value: rc1 = " << rc1 << endl;
    assert(rc1 >= 0.0);
  }

  if ( rc1 != 0.0 )
  {
    const double q1 = 2.0 * M_PI / rc1;
    for ( int i=0; i < wls.size(); i++ )
    {
      if ( g2[i] != 0.0 )
        wls[i] = sqrt(g2[i] / ( g2[i] + q1*q1 ));
      else
        wls[i] = 1.0;
    }
  }
  else
  {
    for ( int i=0; i < wls.size(); i++ )
      wls[i] = 1.0;
  }

  if ( rc_Kerker > 0.0 )
  {
    const double q0_kerker = 2 * M_PI / rc_Kerker;
    const double q0_kerker2 = q0_kerker * q0_kerker;
    for ( int i=0; i < wkerker.size(); i++ )
      wkerker[i] = g2[i] / ( g2[i] + q0_kerker2 );
  }
  else
  {
    for ( int i=0; i < wkerker.size(); i++ )
      wkerker[i] = 1.0;
  }

  // if ultrasoft, calculate position-dependent functions
  if (ultrasoft) {
     tmap["init-usfns"].start();
     cd_.update_usfns();
     wf.update_usfns();
     tmap["init-usfns"].stop();
  }
  // if non-linear core correction defined, calculate position-dependent density
  if (nlcc)
     cd_.update_nlcc();
  
  double tmpthresh = (atoms_move ? s_.ctrl.threshold_force : 0.0 );
  SimpleConvergenceDetector conv_force(s_.ctrl.threshold_force_nsteps, tmpthresh);
  tmpthresh = (compute_stress ? s_.ctrl.threshold_stress : 0.0 );
  SimpleConvergenceDetector conv_stress(s_.ctrl.threshold_stress_nsteps, tmpthresh);

  // calculate time available to avoid exceeding run_timer
  double tbase, tleft;
  bool testtimer = true;
  if (s_.ctrl.run_timer > 0.0) {
    tbase = MPI_Wtime();
    tleft = s_.ctrl.run_timer - (tbase - s_.ctrl.time_init);
  }
  
  // Next line: special case of niter=0: compute GS only        
  for ( int iter = 0; iter < max(niter,1); iter++ )
  {

    // check timing
    if (s_.ctrl.run_timer > 0.0 && niter > 1 && iter > 1 && testtimer) {
      double tnow = MPI_Wtime();
      double sofar = tnow - tbase;
      double tleft = s_.ctrl.run_timer - ( tnow - s_.ctrl.time_init); 
      double titer = sofar/iter; // avg. time per iteration
      double tmaxiter;
      MPI_Allreduce(&titer, &tmaxiter, 1, MPI_DOUBLE, MPI_MAX, s_.ctxt_.comm());
      double maxtime = tmaxiter*(niter-iter);
      if (maxtime > tleft) { // we'll exceed timer, lower niter
        s_.ctrl.timer_hit = true;
        int newiter = (int)(tleft/tmaxiter);
        int tmpniter = iter + newiter - 1;
        if (tmpniter != niter) {
          niter = tmpniter;
          if ( s_.ctxt_.oncoutpe() )
            cout << "<!-- estimated ionic iteration time = " << setprecision(3) << tmaxiter << " sec will exceed run_timer, lowering niter to " << niter << " -->" << endl;
        }
        else {
          testtimer = false;
          if ( s_.ctxt_.oncoutpe() )
            cout << "<!-- estimated ionic iteration time = " << setprecision(3) << tmaxiter << " sec is stable, using niter = " << niter << " -->" << endl;
        }
      }
    }
    
    // ionic iteration
    // test convergence of forces and/or stress
    bool ionic_converge = true;
    if (atoms_move && !conv_force.isConverged())
      ionic_converge = false;
    if (compute_stress && !conv_stress.isConverged())
      ionic_converge = false;
    if (gs_only)
      ionic_converge = false;

    if (ionic_converge) {
      if ( s_.ctxt_.oncoutpe() ) {
        cout.setf(ios::scientific,ios::floatfield);
        cout << "  <!-- BOSampleStepper: ionic convergence reached:  -->" << endl;
        if (atoms_move)
          cout << "  <!-- BOSampleStepper:  ionic convergence, maximum forces varied by less than " << conv_force.threshold() << " a.u. over last " << conv_force.nsteps() << " steps -->" << endl;
        if (compute_stress)
          cout << "  <!-- BOSampleStepper:  ionic convergence, maximum stress varied by less than " << conv_stress.threshold() << " GPa over last " << conv_stress.nsteps() << " steps -->" << endl;
      }
      iter = niter;
    }
    else {

      // check for change in basis size which can cause errors in charge mixing between iterations
      for (int ispin = 0; ispin < nspin; ispin++) {
        //if (s_.wf.spinactive(ispin)) {
          if (rhog_in[ispin].size() != cd_.rhog[ispin].size()) {
            if ( s_.ctxt_.oncoutpe() ) 
              cout << "<ERROR> ChargeDensity basis mismatch!  A ref_cell may be needed. </ERROR>" << endl;
            s_.ctxt_.abort(1);
          }
          //}
      }

      if (iter > 0) {
        // calculate max force and stress values, add to convergence detectors
        if (atoms_move) {
          double maxforce = 0.0;
          for ( int is = 0; is < atoms.atom_list.size(); is++ ) 
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) 
              for (int q=0; q<3; q++) 
                if (fabs(fion[is][3*ia+q]) > maxforce) maxforce = fabs(fion[is][3*ia+q]);

          conv_force.addValue(maxforce);
        }
        if (compute_stress) {
          // ewd:  add more significant figures to conversion
          const double gpa = 29421.0120;
          //const double gpa = 29421.5;
          double maxstress = 0.0;
          for (int q=0; q<6; q++) 
            if (fabs(sigma_eks[q]) > maxstress) maxstress = fabs(sigma_eks[q]);

          conv_stress.addValue(gpa*maxstress);
        }
      }
      

      tm_iter.start();
#ifdef USE_APC
      ApcStart(1);
#endif
      
      if ( onpe0 )
        cout << "<iteration count=\"" << iter+1 << "\">\n";

      if ( ionic_stepper )
        atoms.sync();

      // compute energy and ionic forces using existing wavefunction
      
      if ( !gs_only )
      {
        tmap["charge"].start();
        if ( initial_atomic_density ) {
          cd_.update_rhor();
          if (ultrasoft)
            //ewd DEBUG:  need to renormalize ultrasoft density so tot_charge is correct?
            assert(false);
        }
        else
          cd_.update_density();
        tmap["charge"].stop();

        ef_.update_vhxc();
        const bool compute_forces = true;

        // if symmetry is used, need to calculate set of symmetry-equivalent atoms for 
        // force averaging
        if ( s_.symmetries.nsym() > 0) {
          atoms.findSymmetricAtoms(s_.symmetries);
        }

        // need eigenvalues to compute forces w. ultrasoft
        //ewd12-21-11 if (ultrasoft) { 
        //ewd10-5-12if (ultrasoft && iter == 0) { 
        if (ultrasoft) { 
          tmap["gram"].start();
          s_.wf.gram();
          tmap["gram"].stop();
          ef_.energy(true,dwf,false,fion,false,sigma_eks);
          tmap["diag"].start();
          s_.wf.diag(dwf,true);
          tmap["diag"].stop();
          tmap["usfns"].start();
          s_.wf.update_usfns();
          tmap["usfns"].stop();
          s_.wf.printeig();
        }
        
        double energy =
            ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

        // average forces over symmetric atoms
        if ( compute_forces && s_.symmetries.nsym() > 0) {
          const int nsym_ = s_.symmetries.nsym();
          for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
              // start with identity symmetry operation
              D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
              for ( int isym = 0; isym < nsym_; isym++) {
                int ja = atoms.symatomid(is,ia,isym);
                D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
                //fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
                //ewd need to convert force to crystal coordinates, as symmetry matrix
                // is in crystal coords
                D3vector ftmp_xtal = cell.cart_to_crystal(ftmp);
                D3vector fsym_xtal = s_.symmetries.symlist[isym]->applyToVector(ftmp_xtal,false);
                fsym = fsym + cell.crystal_to_cart(fsym_xtal);
              }
              fion[is][3*ia] = fsym.x/(double)(nsym_+1);
              fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
              fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
            }
          }       
        }

        if (compute_forces && atoms.add_fion_ext()) {
          for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
              D3vector ftmp = atoms.get_fion_ext(is,ia);
              fion[is][3*ia] += ftmp.x;
              fion[is][3*ia+1] += ftmp.y;
              fion[is][3*ia+2] += ftmp.z;
            }
          }
        }

        fmu = - s_.wf.mu() + s_.ctrl.fcp_mu;
        
        if ( onpe0 )
        {
          cout.setf(ios::fixed,ios::floatfield);
          cout.setf(ios::right,ios::adjustfield);
          cout << "  <ekin>   " << setprecision(8)
               << setw(15) << ef_.ekin() << " </ekin>\n";
          if ( use_confinement )
            cout << "  <econf>  " << setw(15) << ef_.econf() << " </econf>\n";
          cout << "  <eps>    " << setw(15) << ef_.eps() << " </eps>\n"
               << "  <enl>    " << setw(15) << ef_.enl() << " </enl>\n"
               << "  <ecoul>  " << setw(15) << ef_.ecoul() << " </ecoul>\n"
               << "  <exc>    " << setw(15) << ef_.exc() << " </exc>\n"
               << "  <esr>    " << setw(15) << ef_.esr() << " </esr>\n"
               << "  <eself>  " << setw(15) << ef_.eself() << " </eself>\n"
               << "  <ets>    " << setw(15) << ef_.ets() << " </ets>\n"
               << "  <etotal> " << setw(15) << ef_.etotal() << " </etotal>\n";
          if ( compute_stress )
          {
            const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
            const double enthalpy = ef_.etotal() + pext * cell.volume();
            cout << "  <pv>     " << setw(15) << pext * cell.volume()
                 << " </pv>" << endl;
            cout << "  <enthalpy> " << setw(15) << enthalpy << " </enthalpy>\n"
                 << flush;
          }
        }
        
        // Add empirical forces (if any) to ions
        double epot_mmion = 0.0;
        if (atoms.empirical_list.size() > 0) {
          const int nspqm_ = atoms.atom_list.size();
          for (int k=0; k<atoms.empirical_list.size(); k++) {
          
            int is1 = atoms.empirical_list[k]->is1;
            int is2 = atoms.empirical_list[k]->is2;

            double prefac = 1.0;
            if (is1 == is2) prefac = 0.5;
            
            vector<Atom*>::iterator pa1; 
            vector<Atom*>::iterator pa2; 
            if (is1 >= nspqm_) 
              pa1 = atoms.mmatom_list[is1-nspqm_].begin();
            else 
              pa1 = atoms.atom_list[is1].begin();
            for (int ia=0; ia<atoms.na(is1); ia++) {
              if (is2 >= nspqm_)
                pa2 = atoms.mmatom_list[is2-nspqm_].begin();
              else 
                pa2 = atoms.atom_list[is2].begin();
              for (int ja=0; ja<atoms.na(is2); ja++) {
                D3vector r12 = (*pa1)->position() - (*pa2)->position();
                cell.fold_in_ws(r12);
                
                double r12val = length(r12);
                if (r12val > 0.0) 
                  epot_mmion += prefac*atoms.empirical_list[k]->pot(r12val);
              
                if ( compute_forces) {
                  D3vector emp_force = prefac*atoms.empirical_list[k]->force(r12);
                  fion[is1][3*ia] += emp_force.x;
                  fion[is1][3*ia+1] += emp_force.y;
                  fion[is1][3*ia+2] += emp_force.z;
                  fion[is2][3*ja] -= emp_force.x;
                  fion[is2][3*ja+1] -= emp_force.y;
                  fion[is2][3*ja+2] -= emp_force.z;
                }
                pa2++;
              }
              pa1++;
            }
          }
        }
      
        if ( iter > 0 && ionic_stepper )
        {
          ionic_stepper->compute_v(energy,fion);
          if ( fcp_stepper )
             fcp_stepper->compute_v(fmu);
        }
        // at this point, positions r0, velocities v0 and forces fion are
        // consistent
        double ekin_ion = 0.0, temp_ion = 0.0;
        if ( ionic_stepper )
        {
          ekin_ion = ionic_stepper->ekin();
          temp_ion = ionic_stepper->temp();
        }
        
        // if hugoniostat is turned on, add current total energy, pressure and volume
        if (s_.ctrl.hugoniostat == "ON") { 
          double pext = 0.0;
          if ( compute_stress )
            pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
          hugstat->addValues(ef_.etotal(),cell.volume(),pext,ionic_stepper->temp());
          if (hugstat->updatenow) {
            double currthtemp = s_.ctrl.th_temp;
            hugstat->updateTemp(currthtemp);
            s_.ctrl.th_temp = currthtemp;
          }
        }

        // print positions, velocities and forces at time t0
        if ( onpe0 )
        {
          cout << "<atomset>" << endl;
          cout << atoms.cell();
          for ( int is = 0; is < atoms.atom_list.size(); is++ )
          {
            int i = 0;
            for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
            {
              Atom* pa = atoms.atom_list[is][ia];
              cout << "  <atom name=\"" << pa->name() << "\""
                   << " species=\"" << pa->species()
                   << "\">\n"
                   << "    <position> " << pa->position() << " </position>\n"
                   << "    <velocity> " << pa->velocity() << " </velocity>\n"
                   << "    <force> "
                   << fion[is][i] << " "
                   << fion[is][i+1] << " "
                   << fion[is][i+2]
                   << " </force>\n";
              cout << "  </atom>" << endl;
              i += 3;
            }
          }
          // MM atoms
          for ( int is = 0; is < atoms.mmatom_list.size(); is++ ) {
            int i = 0;
            const int offset = atoms.atom_list.size();
            for ( int ia = 0; ia < atoms.mmatom_list[is].size(); ia++ ) {
              Atom* pa = atoms.mmatom_list[is][ia];
              cout << "  <mmatom name=\"" << pa->name() << "\""
                   << " species=\"" << pa->species()
                   << "\">\n"
                   << "    <position> " 
                   << ionic_stepper->r0(is+offset,i) << " "
                   << ionic_stepper->r0(is+offset,i+1) << " " 
                   << ionic_stepper->r0(is+offset,i+2) << " </position>\n"
                   << "    <velocity> " 
                   << ionic_stepper->v0(is+offset,i) << " "
                   << ionic_stepper->v0(is+offset,i+1) << " " 
                   << ionic_stepper->v0(is+offset,i+2) << " </velocity>\n"
                   << "    <force> " 
                   << fion[is+offset][i] << " "
                   << fion[is+offset][i+1] << " " 
                   << fion[is+offset][i+2]
                   << " </force>\n  </mmatom>" << endl;
              
              i += 3;
            }
          }
          cout << "</atomset>" << endl;
          cout << "  <econst> " << energy+ekin_ion << " </econst>\n";
          if (atoms.mmatom_list.size() > 0) 
            cout << "  <epot_mmion> " << epot_mmion << " </epot_mmion>\n";
          cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
          cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";
        }
        
        if ( atoms_move )
        {
           if ( s_.constraints.size() > 0 )
          {
            s_.constraints.compute_forces(ionic_stepper->r0(), fion);
            if ( onpe0 )
            {
              s_.constraints.list_constraints(cout);
            }
          }
          // move atoms to new position: r0 <- r0 + v0*dt + dt2/m * fion
          ionic_stepper->compute_r(energy,fion);
          ef_.atoms_moved();
          if ( fcp_stepper )
             fcp_stepper->compute_r(fmu);
          if (ultrasoft) {
             tmap["usfns"].start();
             cd_.update_usfns();
             wf.update_usfns();
             tmap["usfns"].stop();
          }
          tmap["charge"].start();
          if (nlcc)
             cd_.update_nlcc();
          tmap["charge"].stop();
        }

        if ( compute_stress )
        {
          compute_sigma();
          print_stress();

          if ( cell_moves )
          {
            cell_stepper->compute_new_cell(sigma);

            // Update cell
            cell_stepper->update_cell();

            ef_.cell_moved(compute_stress);
            ef_.atoms_moved(); // modifications of the cell also move ions
            if (ultrasoft) {
              tmap["usfns"].start();
              cd_.update_usfns();
              wf.update_usfns();
              tmap["usfns"].stop();
            }
            tmap["charge"].start();
            if (nlcc)
               cd_.update_nlcc();
            tmap["charge"].stop();

            if ( use_preconditioner )
              preconditioner->update();
          }
        }
      } // if !gs_only

      // Recalculate ground state wavefunctions
      
      // wavefunction extrapolation
      if ( atoms_move && extrapolate_wf ) {
	extrapolator.extrapolate_wavefunction(s_.ctrl.wf_extrap, s_.wf, s_.wfv, wfmm, ultrasoft, nspin, iter, dt, s_.ctxt_);
      }

      // do nitscf self-consistent iterations, each with nite electronic steps
      if ( wf_stepper != 0 )
      {
        wf_stepper->preprocess();
        if ( anderson_charge_mixing )
          for (int ispin = 0; ispin < nspin; ispin++)
            //if (s_.wf.spinactive(ispin))
              mixer[ispin]->restart();

        SimpleConvergenceDetector conv_scf(s_.ctrl.threshold_scf_nsteps, s_.ctrl.threshold_scf);
        bool convflag = false;
        double etot_harris;

#ifdef USE_MPIP
        MPI_Pcontrol(1);
#endif        
#ifdef HPM  
  HPM_Start("scfloop");
  summary_start();
#endif
#ifdef TAU
  QB_Pstart(14,scfloop);
#endif
        // SCF LOOP
        tmap["scfloop"].start();
        for ( int itscf = 0; itscf < nitscf_; itscf++ )
        {

          // check timing
          if (niter <= 1 && s_.ctrl.run_timer > 0.0 && nitscf_ > 1 && itscf > 1 && testtimer) {
            double tnow = MPI_Wtime();
            double sofar = tnow - tbase;
            double tleft = s_.ctrl.run_timer - ( tnow - s_.ctrl.time_init); 
            double tscf = sofar/itscf; // avg. time per iteration
            double tmaxiter;
            MPI_Allreduce(&tscf, &tmaxiter, 1, MPI_DOUBLE, MPI_MAX, s_.ctxt_.comm());
            double maxtime = tmaxiter*(nitscf_-itscf);
            if (maxtime > tleft) { // we'll exceed timer, lower nitscf_
              s_.ctrl.timer_hit = true;
              int newiter = (int)(tleft/tmaxiter);
              int tmpnitscf = itscf + newiter - 1;
              if (tmpnitscf != nitscf_) {
                nitscf_ = tmpnitscf;
                if ( s_.ctxt_.oncoutpe() )
                  cout << "<!-- estimated scf iteration time = " << setprecision(3) << tmaxiter << " sec will exceed run_timer, changing nscf to " << nitscf_ << " -->" << endl;
              }
              else {
                testtimer = false;
                if ( s_.ctxt_.oncoutpe() )
                  cout << "<!-- estimated scf iteration time = " << setprecision(3) << tmaxiter << " sec is stable, using nscf = " << nitscf_ << " -->" << endl;
              }
            }
          }
          
          // check convergence
          if (conv_scf.isConverged()) {
            if ( s_.ctxt_.oncoutpe() ) {
              cout.setf(ios::scientific,ios::floatfield);
              if (fractional_occ)
              {
                 cout << "  <!-- BOSampleStepper: scf convergence at itscf = " << itscf << ", Harris-Foulkes energy varied by less than " << setprecision(2) 
                   << conv_scf.threshold() << " a.u. over " << conv_scf.nsteps() 
                   << " scf steps. -->" << endl;
              cout.flush();
              }
              else
              {
                 cout << "  <!-- BOSampleStepper: scf convergence at itscf = " << itscf << ", scf energy varied by less than " << setprecision(2) 
                   << conv_scf.threshold() << " a.u. over " << conv_scf.nsteps() 
                   << " scf steps. -->" << endl;
              cout.flush();
              }
            }
            itscf = nitscf_;
            convflag = true;
          }          
          // continue itscf loop
          else {
            if (itscf > 0) {
               if (ultrasoft || !fractional_occ)
                  conv_scf.addValue(ef_.etotal()); // Harris-Foulkes estimate not implemented for ultrasoft yet (need enl_)
               else
                  conv_scf.addValue(etot_harris);
            }            
            if ( nite_ > 1 && onpe0 )
              cout << "  <!-- BOSampleStepper: start scf iteration -->" << endl;

            // compute new density in cd_.rhog
            tmap["charge"].start();
            if ( itscf==0 && initial_atomic_density ) {
              cd_.update_rhor();
              if (ultrasoft)
                assert(false);  // ewd DEBUG: need to implement this for ultrasoft
            }
            else
               // QB_Pstart("charge");
               cd_.update_density();
               // QB_Pstop("charge");
            tmap["charge"].stop();

            if (fractional_occ)
            {
               tmap["scf_ef"].start();
               ef_.update_harris();
               tmap["scf_ef"].stop();
            }
            
            if ( charge_mixing != "off" && nite_ > 0) {
              if ( itscf == 0) {
                //ewd:  read rhog_in from checkpoint if possible
                for ( int ispin = 0; ispin < nspin; ispin++ ) {
                  int rhogflag = 1;

                  if (cell_dyn == "LOCKED" && (atoms_dyn == "LOCKED" || dt == 0.0) && s_.rhog_last.size() == rhog_in.size()) {
                     if (s_.rhog_last[ispin].size() == rhog_in[ispin].size()) {
                        rhogflag = 0;
                        for ( int i=0; i < rhog_in[ispin].size(); i++ ) 
                           rhog_in[ispin][i] = s_.rhog_last[ispin][i];
                     }
                  }
                  
                  if (rhogflag) {
                    for ( int i=0; i < rhog_in[ispin].size(); i++ ) 
                      rhog_in[ispin][i] = cd_.rhog[ispin][i];
                  }
                }
              }
              // itscf > 0
              else {

                for ( int ispin = 0; ispin < nspin; ispin++ ) {
                  // compute correction drhog
                  for ( int i=0; i < rhog_in[ispin].size(); i++ )
                  {
                    drhog[ispin][i] = (cd_.rhog[ispin][i] - rhog_in[ispin][i]);
                  }

                  const double alpha = s_.ctrl.charge_mix_coeff;
                  // Anderson acceleration
                  if ( anderson_charge_mixing )
                  {
                    // row weighting of LS calculation
                    for ( int i=0; i < drhog[ispin].size(); i++ )
                      drhog[ispin][i] /= wls[i];
                    
                    //ewd:  try running this on every task
                    mixer[ispin]->update((double*)&rhog_in[ispin][0],(double*)&drhog[ispin][0],
                                         (double*)&rhobar[ispin][0],(double*)&drhobar[ispin][0]);

                    for ( int i=0; i < drhobar[ispin].size(); i++ )
                      drhobar[ispin][i] *= wls[i];
                    
                    for ( int i=0; i < rhog_in[ispin].size(); i++ )
                      rhog_in[ispin][i] = rhobar[ispin][i] + alpha * drhobar[ispin][i] * wkerker[i];
                  }
                  else {
                    for ( int i=0; i < rhog_in[ispin].size(); i++ )
                      rhog_in[ispin][i] += alpha * drhog[ispin][i] * wkerker[i];
                  }
                  
                  // Apply correction
                  for ( int i=0; i < rhog_in[ispin].size(); i++ )
                    cd_.rhog[ispin][i] = rhog_in[ispin][i];
                }
                cd_.update_rhor();
              }              
            }

            //QB_Pstart(update_vhxc);
            tmap["scf_ef"].start();
            ef_.update_vhxc();
            tmap["scf_ef"].stop();
            //QB_Pstop(update_vhxc);

            // reset stepper only if multiple non-selfconsistent steps
            if ( nite_ > 0 ) wf_stepper->preprocess();
            int nitemin_ = ( nite_ > 0 ? nite_ : 1);
            for ( int ite = 0; ite < nitemin_; ite++ )
            {
               //QB_Pstart(energy+hamiltonian_update);
               tmap["scf_ef"].start();
               double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
               tmap["scf_ef"].stop();
               //QB_Pstop(energy+hamiltonian_update);

               // compute the sum of eigenvalues (with fixed weight)
               // to measure convergence of the subspace update
               // compute trace of the Hamiltonian matrix Y^T H Y
               // scalar product of Y and (HY): tr Y^T (HY) = sum_ij Y_ij (HY)_ij
               // Note: since the hamiltonian is hermitian and dwf=H*wf
               // the dot product in the following line is real

               if (ultrasoft) { 
                  const double us_eigenvalue_sum = s_.wf.sdot(dwf);
                  if ( onpe0 )
                     cout  << "  <eigenvalue_sum> "
                           << us_eigenvalue_sum << " </eigenvalue_sum>" << endl;
               }
               else {
                  const double eigenvalue_sum = s_.wf.dot(dwf);
                  if ( onpe0 )
                     cout << "  <eigenvalue_sum> "
                          << eigenvalue_sum << " </eigenvalue_sum>" << endl;
               }
               //QB_Pstart(pzgemm used to update wf)               
               wf_stepper->update(dwf); 
               //QB_Pstop(pzgemm used to update wf)
               
              if (ultrasoft)
                  wf.update_usfns();

               // update ultrasoft functions if needed, call gram
               for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ ) {
                  if (s_.wf.spinactive(ispin)) {
                     for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ ) {
                        if (s_.wf.kptactive(ikp)) {
                           assert(s_.wf.sd(ispin,ikp) != 0);
                           //if (ultrasoft) { 
                           //  tmap["usfns"].start();
                           //  s_.wf.sd(ispin,ikp)->update_usfns(); // calculate betapsi, spsi
                           //  tmap["usfns"].stop();
                           //}
                           tmap["gram"].start();
                           s_.wf.sd(ispin,ikp)->gram();
                           tmap["gram"].stop();
                           
                           //if (ultrasoft) { 
                           //  tmap["usfns"].start();
                           //  s_.wf.sd(ispin,ikp)->calc_betapsi(); // calculate betapsi
                           //  tmap["usfns"].stop();
                           //}
                        }
                     }
                  }
               }
               
               if ( onpe0 )
               {
                  cout.setf(ios::fixed,ios::floatfield);
                  cout.setf(ios::right,ios::adjustfield);
                  cout << "  <etotal_int scf_iter=\"" << itscf << "\"> " << setw(15) << setprecision(8)
                       << energy << " </etotal_int>\n";
                  if ( compute_stress )
                  {
                     const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
                     const double enthalpy = energy + pext * cell.volume();
                     cout << "  <enthalpy_int> " << setw(15)
                          << enthalpy << " </enthalpy_int>\n"
                          << flush;
                  }
               }
               if (ite == 0)
               {
                  if (fractional_occ)
                  {
                     etot_harris = ef_.etotal_harris();
                     if (onpe0)
                        cout << "  <eharris> " << setw(15) << setprecision(8) << etot_harris << " </eharris>\n";

                     //ewd DEBUG
                     if (false && onpe0)
                     {
                        cout << "EHARRIS.DEBUG:  " << ef_.ekin() << "  " << ef_.econf() << "  " << ef_.enl() << "  " << ef_.ets() << "  " << ef_.epv() << "  " << ef_.ehub() << "  " << ef_.eharris() << "  " << ef_.esr() << "  " << ef_.eself() << endl;
                     }
                     //ewd DEBUG

                  }
               }
            } // for ite

            // subspace diagonalization
            if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" || usdiag)
            {
               tmap["scf_ef"].start();
               ef_.energy(true,dwf,false,fion,false,sigma_eks);
               tmap["scf_ef"].stop();
               tmap["diag"].start();
               s_.wf.diag(dwf,compute_eigvec);
               tmap["diag"].stop();

               // update ultrasoft functions w. new eigenvectors
               if (compute_eigvec && ultrasoft)
                  s_.wf.update_usfns();

               if (itscf%s_.ctrl.iprint == 0)
                  s_.wf.printeig();
            }

            // update occupation numbers if fractionally occupied states
            if ( fractional_occ )
            {
              wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);
              if (itscf%s_.ctrl.iprint == 0)
                s_.wf.printocc();
              const double wf_entropy = wf.entropy();
              if ( onpe0 )
              {
                cout << "  <!-- Wavefunction entropy: " << wf_entropy << " -->" << endl;
                //const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
                cout << "  <!-- Entropy contribution to free energy: "
                     << - wf_entropy * s_.ctrl.smearing_width * 2.0 << " -->" << endl;
              }
            }
          
            if ( nite_ > 1 && onpe0 )
              cout << "  <!-- BOSampleStepper: end scf iteration -->" << endl;
          }
        } // for itscf
        tmap["scfloop"].stop();

         //ewd DEBUG
        if (false)
        {
            double eewald = ef_.casino_ewald();
            double evloc = ef_.casino_vloc();
            double omega = wf.cell().volume();
            if ( onpe0 )
               cout << "CASINO.debug, ewald = " << eewald << ", eloc = " << evloc << ", volume = " << wf.cell().volume() << ", ehart = " << ef_.ehart() << endl;
         }
#ifdef TAU  
  QB_Pstop(scfloop);
#endif
#ifdef HPM
  summary_stop();
  HPM_Stop("scfloop");
#endif
#ifdef USE_MPIP
        MPI_Pcontrol(0);
#endif        

        if (!convflag && s_.ctrl.threshold_scf > 0.0) 
          if ( s_.ctxt_.oncoutpe() ) 
            cout << "<WARNING> Ionic iteration finished without triggering scf convergence threshold, consider increasing nscf = " << nitscf_ << " </WARNING>" << endl;
                      
        if ( compute_mlwf || compute_mlwfc )
        {
           tmap["mlwf"].start();
           SlaterDet& sd = *(wf.sd(0,0));
          mlwft->compute_transform();

          if ( compute_mlwf )
            mlwft->apply_transform(sd);
          
          if ( onpe0 )
          {
            cout << " <mlwf_set size=\"" << sd.nst() << "\">" << endl;
            for ( int i = 0; i < sd.nst(); i++ )
            {
              D3vector ctr = mlwft->center(i);
              double sp = mlwft->spread(i);
              cout.setf(ios::fixed, ios::floatfield);
              cout.setf(ios::right, ios::adjustfield);
              cout << "   <mlwf center=\"" << setprecision(6)
                   << setw(12) << ctr.x
                   << setw(12) << ctr.y
                   << setw(12) << ctr.z
                   << " \" spread=\" " << sp << " \"/>"
                   << endl;
            }
            cout << " </mlwf_set>" << endl;
            D3vector edipole = mlwft->dipole();
            cout << " <electronic_dipole> " << edipole
                 << " </electronic_dipole>" << endl;
            D3vector idipole = atoms.dipole();
            cout << " <ionic_dipole> " << idipole
                 << " </ionic_dipole>" << endl;
            cout << " <total_dipole> " << idipole + edipole
                 << " </total_dipole>" << endl;
            cout << " <total_dipole_length> " << length(idipole + edipole)
                 << " </total_dipole_length>" << endl;
          }
          tmap["mlwf"].stop();
        }


        // If GS calculation only, print energy and atomset at end of iterations


        bool fastend = false;

        if ( gs_only && !fastend)
        {

           tmap["postscf"].start();
           // need eigenvalues to compute forces w. ultrasoft
           if (ultrasoft) { 
              ef_.energy(true,dwf,false,fion,false,sigma_eks);
              //tmap["diag"].start();
              s_.wf.diag(dwf,compute_eigvec);
              //s_.wf.diag(dwf,true);  // ewd:  why true?  Why did I do this??
              //tmap["diag"].stop();
              
              // update ultrasoft functions w. new eigenvectors
              if (compute_eigvec && ultrasoft) { 
                 //tmap["usfns"].start();
                 s_.wf.update_usfns();
                 //tmap["usfns"].stop();
              }
              s_.wf.printeig();
           }
           
           //tmap["charge"].start();
           if ( !initial_atomic_density )
              cd_.update_density();
           //tmap["charge"].stop();

          ef_.update_vhxc();
          const bool compute_forces = true;
          ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

          if (s_.atoms.add_fion_ext()) {
            for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
              for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
                D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
                fion[is][3*ia] += ftmp.x;
                fion[is][3*ia+1] += ftmp.y;
                fion[is][3*ia+2] += ftmp.z;
              }
            }
          }

          if ( onpe0 )
          {
            cout << ef_;
            cout << "<atomset>" << endl;
            cout << atoms.cell();
            for ( int is = 0; is < atoms.atom_list.size(); is++ )
            {
              int i = 0;
              for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
              {
                Atom* pa = atoms.atom_list[is][ia];
                cout << "  <atom name=\"" << pa->name() << "\""
                     << " species=\"" << pa->species()
                     << "\">\n"
                     << "    <position> " << pa->position() << " </position>\n"
                     << "    <velocity> " << pa->velocity() << " </velocity>\n"
                     << "    <force> "
                     << fion[is][i] << " "
                     << fion[is][i+1] << " "
                     << fion[is][i+2]
                     << " </force>\n";
                cout << "  </atom>" << endl;
                i += 3;
              }
            }
            cout << "</atomset>" << endl;
            if ( compute_stress )
            {
              compute_sigma();
              print_stress();
            }
          }
          tmap["postscf"].stop();
        }
        else if (gs_only && fastend) {
           if ( onpe0 )
              cout << ef_;
        }
        wf_stepper->postprocess();
      }
      else
      {
         tmap["postscf"].start();
         // wf_stepper == 0, wf_dyn == LOCKED
         // evaluate and print energy
         //tmap["charge"].start();
         cd_.update_density();
         //tmap["charge"].stop();
         ef_.update_vhxc();
         ef_.energy(true,dwf,false,fion,false,sigma_eks);
         if ( onpe0 )
         {
            cout << ef_;
         }
         tmap["postscf"].stop();
      }

#ifdef USE_APC
      ApcStop(1);
#endif
      // print iteration time
      double time = tm_iter.real();
      double tmin = time;
      double tmax = time;
      s_.ctxt_.dmin(1,1,&tmin,1);
      s_.ctxt_.dmax(1,1,&tmax,1);
      if ( onpe0 )
      {
       cout << left << setw(34) << "<timing where=\"run\""
           << setw(24) << " name=\" iteration\""
             << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
             << " max=\"" << setprecision(3) << setw(9) << tmax << "\""
             << " count=\"" << setw(9) << 1 << "\"/>"
             << endl;
        cout << "</iteration>" << endl;
      }

      s_.ctrl.mditer++;

      // if savedenfreq variable set, save density in VMD Cube format
      if (s_.ctrl.savedenfreq > 0)
      {
         if (s_.ctrl.mditer%s_.ctrl.savedenfreq == 0 || s_.ctrl.mditer == 1)
         {
            //string filebase = "density.";
            string filebase = s_.ctrl.savedenfilebase;
            ostringstream oss;
            oss.width(7);  oss.fill('0');  oss << s_.ctrl.mditer;
            string denfilename = filebase + "." + oss.str() + ".cube";
            string format = "binary";

            const Context* wfctxt = s_.wf.spincontext(0);
            FourierTransform* ft_ = cd_.vft();
            if (wfctxt->mycol() == 0) {
               vector<double> rhortmp(ft_->np012loc());
               for (int j = 0; j < ft_->np012loc(); j++)
                  rhortmp[j] = cd_.rhor[0][j];
    
               for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                  if ( i == wfctxt->myrow() ) {
                     int size = ft_->np012loc();
                     wfctxt->isend(1,1,&size,1,0,0);
                     wfctxt->dsend(size,1,&rhortmp[0],1,0,0);
                  }
               }
               if ( wfctxt->oncoutpe() ) {
                  ofstream os;
                  os.open(denfilename.c_str(),ofstream::out);    // text output
                  os.setf(ios::fixed,ios::floatfield);
                  os << setprecision(8);
                  vector<double> tmprecv(ft_->np012());
                  int recvoffset = 0;

                  D3vector a0 = s_.wf.cell().a(0);
                  D3vector a1 = s_.wf.cell().a(1);
                  D3vector a2 = s_.wf.cell().a(2);
                  const int np0 = ft_->np0();
                  const int np1 = ft_->np1();
                  const int np2 = ft_->np2();
                  D3vector dft0 = a0/(double)np0;
                  D3vector dft1 = a1/(double)np1;
                  D3vector dft2 = a2/(double)np2;
                  
                  for ( int i = 0; i < wfctxt->nprow(); i++ ) {
                     int size = 0;
                     wfctxt->irecv(1,1,&size,1,i,0);
                     wfctxt->drecv(size,1,&tmprecv[recvoffset],1,i,0);
                     recvoffset += size;

                     if (i==0) {
                        // write out VMD CUBE format header
                        os << "Qbox wavefunction in VMD CUBE format" << endl;
                        os << "  electron density" << endl;

                        // get atom positions
                        AtomSet& as = s_.atoms;
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
                  }

                  // write density data to file
                  int cnt = 0;
                  for (int ii = 0; ii < np0; ii++) {
                     ostringstream oss;
                     oss.setf(ios::scientific,ios::floatfield);
                     oss << setprecision(5);
                     for (int jj = 0; jj < np1; jj++) {
                        for (int kk = 0; kk < np2; kk++) {
                           int index = ii + jj*np0 + kk*np0*np1;
                           oss << tmprecv[index] << " ";
                           cnt++;
                           if (cnt >= 6) {
                              cnt = 0;
                              oss << endl;
                           }
                        }
                     }
                     string tos = oss.str();
                     os.write(tos.c_str(),tos.length());
                  }
                  os.close();
               }
            }
         }
      }

      // if savewffreq variable set, save |wf|^2 (one or all) in VMD Cube format
      if (s_.ctrl.savewffreq > 0)
      {
         if (s_.ctrl.mditer%s_.ctrl.savewffreq == 0 || s_.ctrl.mditer == 1)
         {
            string filebase = s_.ctrl.savewffilebase;
            ostringstream oss;
            oss.width(7);  oss.fill('0');  oss << s_.ctrl.mditer;
            string wffilename = filebase + "." + oss.str() + ".cube";
            s_.wf.print_vmd(wffilename,atoms,s_.ctrl.savewfstate);            
         }
      }

      // if savefreq variable set, checkpoint
      if (s_.ctrl.savefreq > 0)
      {
         if (s_.ctrl.savefreq == 1 || (s_.ctrl.mditer > 0 && s_.ctrl.mditer%s_.ctrl.savefreq == 0) )
         {
            tmap["mdsave"].start();

            // compute ionic forces at last position to update velocities
            // consistently with last position
            cd_.update_density();

            ef_.update_vhxc();
            const bool compute_forces = true;
            double energy =
                ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

            // average forces over symmetric atoms
            if ( compute_forces && s_.symmetries.nsym() > 0) {
               int nsym_ = s_.symmetries.nsym();
               for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
                  for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
                     // start with identity symmetry operation
                     D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
                     for ( int isym = 0; isym < nsym_; isym++) {
                        int ja = s_.atoms.symatomid(is,ia,isym);
                        D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
                        fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
                     }
                     fion[is][3*ia] = fsym.x/(double)(nsym_+1);
                     fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
                     fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
                  }
               }       
            }

            if (compute_forces && s_.atoms.add_fion_ext()) {
               for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
                  for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
                     D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
                     fion[is][3*ia] += ftmp.x;
                     fion[is][3*ia+1] += ftmp.y;
                     fion[is][3*ia+2] += ftmp.z;
                  }
               }
            }
      
            ionic_stepper->compute_v(energy,fion);
            // positions r0 and velocities v0 are consistent
            if ( fcp_stepper )
               fcp_stepper->compute_v(fmu);

            // create output directory if it doesn't exist
            string dirbase = "md.";
            string filebase = "mdchk";
            ostringstream oss;
            oss.width(7);  oss.fill('0');  oss << s_.ctrl.mditer;
            string dirstr = dirbase + oss.str();
            string format = "binary";
            if ( onpe0 )
            {
               int mode = 0775;
               struct stat statbuf;
               int rc = stat(dirstr.c_str(), &statbuf);
               if (rc == -1)
               {
                  cout << "Creating directory: " << dirstr << endl;
                  rc = mkdir(dirstr.c_str(), mode);
                  rc = stat(dirstr.c_str(), &statbuf);
               }
               if (rc != 0 || !(statbuf.st_mode))
               {
                  cout << "<ERROR> Can't stat directory " << dirstr << " </ERROR> " << endl;
                  MPI_Abort(MPI_COMM_WORLD,2);
               }
            }
            string filestr = dirstr + "/" + filebase;
            s_.wf.write_states(filestr,format);
            s_.wf.write_mditer(filestr,s_.ctrl.mditer);

            // write .sys file
            if ( onpe0 )
            {
               string sysfilename = dirstr + "/" + "mdsave.sys";
               ofstream os;
               os.open(sysfilename.c_str(),ofstream::out);

               // cell info
               string cmd("set cell ");
               s_.wf.cell().printsys(os,cmd);
               
               // ref cell info, if necessary
               if ( s_.wf.refcell().volume() != 0.0 ) {
                  string refcmd("set ref_cell ");
                  s_.wf.refcell().printsys(os,refcmd);
               }

               // species info
               const int nspqm_ = s_.atoms.nsp();
               for (int i=0; i<nspqm_; i++)
                  s_.atoms.species_list[i]->printsys(os);
               
               const int nspmm_ = s_.atoms.nsp_mm();
               for (int i=0; i<nspmm_; i++)
                  s_.atoms.mmspecies_list[i]->printsys(os);
               
               // atom coordinates and info
               s_.atoms.printsys(os);
               os.close();
            }
            tmap["mdsave"].stop();
         }
      }

      if ( atoms_move )
        s_.constraints.update_constraints(dt);

    } // else (not converged)
  } // for iter
  // print memory usage of main data objects
  if (s_.ctxt_.oncoutpe()) {
    double memtot = 0.0;
    double loctot = 0.0;

    int kmult = s_.wf.nspin()*s_.wf.nkp();
    int kmultloc = s_.wf.nkptloc();
    ef_.print_memory(cout,memtot,loctot);
    s_.wf.sd(0,0)->print_memory(cout,kmult,kmultloc,memtot,loctot);

    cd_.print_memory(cout,memtot,loctot);
    
    PrintMem pm;
    string totunit = pm.memunit(memtot);
    string locunit = pm.memunit(loctot);
    cout << "<!-- memory total       :  " << setw(7) << memtot << totunit << "  (" << loctot << locunit << " local) -->" << endl;
  }
  

#ifdef PRINTALL
  //ewd print timing
  cd_.print_timing();
  ef_.print_timing();
  s_.wf.print_timing();

  // print timer map
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
     double time = (*i).second.real();
     double tmin = time;
     double tmax = time;
     s_.ctxt_.dmin(1,1,&tmin,1);
     s_.ctxt_.dmax(1,1,&tmax,1);
     uint64_t count = (*i).second.counts();
     if ( s_.ctxt_.mype()==0 )
     {
        cout << left << setw(34) << "<timing where=\"run\""
             << setw(8) << " name=\""
             << setw(15) << (*i).first << "\""
             << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
             << " max=\"" << setprecision(3) << setw(9) << tmax << "\""
             << " count=\"" << setw(9) << count << "\"/>"
             << endl;
     }
  }
#endif  
  
  
  //ewd save last charge density 
  if (cell_dyn == "LOCKED") {
    s_.rhog_last.resize(s_.wf.nspin());
    for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
    {
      if (s_.wf.spinactive(ispin))
      {
        s_.rhog_last[ispin].resize(rhog_in[ispin].size());
        for ( int i=0; i < rhog_in[ispin].size(); i++ )
          s_.rhog_last[ispin][i] = rhog_in[ispin][i];
      }
    }
  }
  
  if ( atoms_move )
  {

    // need eigenvalues to compute forces w. ultrasoft
    if (ultrasoft) { 
       ef_.energy(true,dwf,false,fion,false,sigma_eks);
      tmap["diag"].start();
      //s_.wf.diag(dwf,compute_eigvec);
      s_.wf.diag(dwf,true);
      tmap["diag"].stop();
      s_.wf.printeig();
    }

    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();

    ef_.update_vhxc();
    const bool compute_forces = true;
    double energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

    // average forces over symmetric atoms
    if ( compute_forces && s_.symmetries.nsym() > 0) {
      int nsym_ = s_.symmetries.nsym();
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          // start with identity symmetry operation
          D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
          for ( int isym = 0; isym < nsym_; isym++) {
            int ja = s_.atoms.symatomid(is,ia,isym);
            D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
            fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
          }
          fion[is][3*ia] = fsym.x/(double)(nsym_+1);
          fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
          fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
        }
      }       
    }

    if (compute_forces && s_.atoms.add_fion_ext()) {
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
          fion[is][3*ia] += ftmp.x;
          fion[is][3*ia+1] += ftmp.y;
          fion[is][3*ia+2] += ftmp.z;
        }
      }
    }
      
    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
    if ( fcp_stepper )
       fcp_stepper->compute_v(fmu);
  }

  if ( atoms_move && extrapolate_wf )
  {
    // compute wavefunction velocity after last iteration
    // s_.wfv contains the previous wavefunction

    s_.wfv->align(s_.wf);

    if (ultrasoft) {
      for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ ) {
        if (s_.wf.spinactive(ispin)) {
          for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ ) {
            if (s_.wf.kptactive(ikp)) {
              assert(s_.wf.sd(ispin,ikp) != 0);
              tmap["usfns"].start();
              s_.wf.sd(ispin,ikp)->update_usfns();
              tmap["usfns"].stop();
            }
          }
        }
      }
    }

    for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
    {
      if (s_.wf.spinactive(ispin))
      {
        for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
        {
          if (s_.wf.kptactive(ikp))
          {
            assert(s_.wf.sd(ispin,ikp) != 0);
            double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
            double* cm = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
            const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
            const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
            const int len = 2*mloc*nloc;
            double dt_inv;
            if (dt == 0.0) dt_inv = 0.0;
            else dt_inv = 1.0 / dt;
            if (s_.ctrl.wf_extrap == "NTC")
            {
              double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xmm = cmm[i];
                cm[i] = dt_inv * ( x - xmm );
              }
              tmap["gram"].start();
              s_.wf.sd(ispin,ikp)->gram();
              tmap["gram"].stop();
            }
            else if (s_.ctrl.wf_extrap == "SIMPLE" || s_.ctrl.wf_extrap == "ASP" )
            {
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xm = cm[i];
                cm[i] = dt_inv * ( x - xm );
              }
              tmap["gram"].start();
              s_.wf.sd(ispin,ikp)->gram();
              tmap["gram"].stop();
            }
          }
        }
      }
    }

    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();

    ef_.update_vhxc();
    const bool compute_forces = true;
    double energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

    // average forces over symmetric atoms
    if ( compute_forces && s_.symmetries.nsym() > 0) {
      int nsym_ = s_.symmetries.nsym();
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          // start with identity symmetry operation
          D3vector fsym(fion[is][3*ia],fion[is][3*ia+1],fion[is][3*ia+2]);
          for ( int isym = 0; isym < nsym_; isym++) {
            int ja = s_.atoms.symatomid(is,ia,isym);
            D3vector ftmp(fion[is][3*ja],fion[is][3*ja+1],fion[is][3*ja+2]);
            fsym = fsym + s_.symmetries.symlist[isym]->applyToVector(ftmp,false);
          }
          fion[is][3*ia] = fsym.x/(double)(nsym_+1);
          fion[is][3*ia+1] = fsym.y/(double)(nsym_+1);
          fion[is][3*ia+2] = fsym.z/(double)(nsym_+1);
        }
      }       
    }

    if (compute_forces && s_.atoms.add_fion_ext()) {
      for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
          D3vector ftmp = s_.atoms.get_fion_ext(is,ia);
          fion[is][3*ia] += ftmp.x;
          fion[is][3*ia+1] += ftmp.y;
          fion[is][3*ia+2] += ftmp.z;
        }
      }
    }
      
    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
    if ( fcp_stepper )
       fcp_stepper->compute_v(fmu);
  }
  else
  {
    // delete wavefunction velocities
    if ( s_.wfv != 0 )
      delete s_.wfv;
    s_.wfv = 0;
  }

  delete mlwft;

  // delete steppers
  delete wf_stepper;
  delete ionic_stepper;
  delete cell_stepper;

  // delete preconditioner
  if ( use_preconditioner ) delete preconditioner;
  if ( s_.ctrl.wf_extrap == "NTC" || s_.ctrl.wf_extrap == "NTC" ) 
     delete wfmm;

  // delete Hugoniostat
  if (hugstat != 0) delete hugstat;

  for (int ispin = 0; ispin < nspin; ispin++) 
    if (s_.wf.spinactive(ispin))
      delete mixer[ispin];

  initial_atomic_density = false;
}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::get_forces(vector<vector<double> >& f) const {
  AtomSet& atoms = s_.atoms;
  if (f.size() < atoms.atom_list.size()) 
    f.resize(atoms.atom_list.size());

  for ( int is = 0; is < atoms.atom_list.size(); is++ ) {
    if (f[is].size() != 3*atoms.atom_list[is].size())
      f[is].resize(3*atoms.atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ ) {
      f[is][i] = fion[is][i];
      f[is][i+1] = fion[is][i+1];
      f[is][i+2] = fion[is][i+2];
      i+=3;
    }
  }

  // get forces on MM atoms, if any
  if (f.size() < atoms.atom_list.size()+atoms.mmatom_list.size())
    f.resize(atoms.atom_list.size() + atoms.mmatom_list.size());
  const int offset = atoms.atom_list.size();
  for ( int is = 0; is < atoms.mmatom_list.size(); is++ ) {
    if (f[is+offset].size() != 3*atoms.mmatom_list[is].size())
      f[is+offset].resize(3*atoms.atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atoms.mmatom_list[is].size(); ia++ ) {
      D3vector t = atoms.mmatom_list[is][ia]->position();
      f[is+offset][i] = fion[is+offset][i];
      f[is+offset][i+1] = fion[is+offset][i+1];
      f[is+offset][i+2] = fion[is+offset][i+2];
      i+=3;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double BOSampleStepper::get_energy(string ename) {
  if (ename == "ekin") { return ef_.ekin(); }
  else if (ename == "econf") { return ef_.econf(); }
  else if (ename == "eps") { return ef_.eps(); }
  else if (ename == "enl") { return ef_.enl(); }
  else if (ename == "ehart") { return ef_.ehart(); }
  else if (ename == "ecoul") { return ef_.ecoul(); }
  else if (ename == "exc") { return ef_.exc(); }
  else if (ename == "esr") { return ef_.esr(); }
  else if (ename == "eself") { return ef_.eself(); }
  else if (ename == "ets") { return ef_.ets(); }
  else if (ename == "etotal") { return ef_.etotal(); }
  else { return ef_.etotal(); }
}

////////////////////////////////////////////////////////////////////////////////
valarray<double> BOSampleStepper::get_stress(string sname) {
  if (sname == "total") { return sigma; }
  else if (sname == "ext") { return sigma_ext; }
  else if (sname == "eks") { return sigma_eks; }
  else if (sname == "kin") { return sigma_kin; }
  else { 
    valarray<double> nullv(1);
    nullv[0] = 0.0;
    return nullv; 
  }
}
