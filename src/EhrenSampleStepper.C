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
// EhrenSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "EhrenSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SOTDWavefunctionStepper.h"
#include "SORKTDWavefunctionStepper.h"
#include "FORKTDWavefunctionStepper.h"
#include "TDEULERWavefunctionStepper.h"
#include "SDIonicStepper.h"
#include "SDAIonicStepper.h"
#include "CGIonicStepper.h"
#include "MDIonicStepper.h"
#include "BMDIonicStepper.h"
#include "SDCellStepper.h"
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
#include <deque>
#ifdef HPM
#include <bgpm/include/bgpm.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
EhrenSampleStepper::EhrenSampleStepper(Sample& s, int nitscf, int nite) :
  SampleStepper(s), cd_(s), ef_(s,s.wf,cd_),dwf(s.wf), wfv(s.wfv), nitscf_(nitscf),
  nite_(nite)
{
   const string wf_dyn = s_.ctrl.wf_dyn;
   tddft_involved_ = s_.ctrl.tddft_involved;

   Wavefunction* wf0 = new Wavefunction(dwf);
   Wavefunction* wf1 = new Wavefunction(s.wf);
   Wavefunction* wf2 = new Wavefunction(s.wf);
   Wavefunction* wf3 = new Wavefunction(s.wf);
   Wavefunction* wf4 = new Wavefunction(s.wf);
   
   wfdeque.push_back ( wf0 );
   wfdeque.push_back ( wf1 );
   wfdeque.push_back ( wf2 );
   wfdeque.push_back ( wf3 );
   wfdeque.push_back ( wf4 );

}
////////////////////////////////////////////////////////////////////////////////
EhrenSampleStepper::~EhrenSampleStepper()
{
   // delete any created Wavefunctions
   for (deque<Wavefunction*>::const_iterator it = wfdeque.begin(); it != wfdeque.end(); ++it)
      delete *it;
}
////////////////////////////////////////////////////////////////////////////////
void EhrenSampleStepper::step(int niter)
{
  const bool oncoutpe = s_.ctxt_.oncoutpe();
  if (!tddft_involved_)
  {
     if ( oncoutpe )
        cout << "<ERROR> EhrenSampleStepper can not be used when tddft_involved = false!</ERROR>" << endl;
      return;
   }

  // determine whether eigenvectors must be computed
  // eigenvectors are computed if explicitly requested with wf_diag==T
  // or if the SlaterDet has fractionally occupied states
  int occtest = (2 * s_.wf.nst()) - s_.wf.nspin() * s_.wf.nel();
  const bool fractional_occ = (occtest != 0 && occtest != 1);
  const bool compute_eigvec = fractional_occ || s_.ctrl.wf_diag == "T";
  const bool compute_mlwf = s_.ctrl.wf_diag == "MLWF";
  const bool compute_mlwfc = s_.ctrl.wf_diag == "MLWFC";
  enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };

  if (fractional_occ && oncoutpe)
       cout << "<!-- EhrenSampleStepper:  fractional occupation detected. -->" << endl;
  else if (oncoutpe)
       cout << "<!-- EhrenSampleStepper:  fractional occupation not detected. -->" << endl;
   
  if (s_.ctrl.reshape_context)
     s_.wf.set_reshape_context(s_.ctrl.reshape_context);
  
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  const int nspin = wf.nspin();
  const UnitCell& cell = wf.cell();
  const double dt = s_.ctrl.dt;
  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;

  // AS: flag used below to control initial step for SOTD
  bool wfv_is_new = false;
  const bool extrapolate_wf = ( atoms_dyn == "MD" ) && ( nite_ == 1 ) && (!tddft_involved_);

  const bool atoms_move = ( niter > 0 && atoms_dyn != "LOCKED" );
  const bool compute_forces = true;
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool cell_moves = ( niter > 0 && compute_stress &&
                            cell_dyn != "LOCKED" );
  const bool use_confinement = ( s_.ctrl.ecuts > 0.0 );

  const bool ultrasoft = s_.ctrl.ultrasoft;
  const bool usdiag = (ultrasoft && atoms_move);
  const bool nlcc = s_.ctrl.nlcc;
  cd_.set_nlcc(nlcc);
  
  //ewd check that MLWF not being used with ultrasoft (not yet implemented)
  if (ultrasoft && (compute_mlwf || compute_mlwfc)) {
    if ( oncoutpe ) 
      cout << "<ERROR> EhrenSampleStepper:  Maximally-localized Wannier Functions not yet implemented with ultrasoft. </ERROR>" << endl;
    return;
  }
  if (oncoutpe && ultrasoft && compute_stress)
      cout << "<WARNING> EhrenSampleStepper:  stress not yet implemented with ultrasoft!  Results WILL be wrong. </WARNING>" << endl;
  if (oncoutpe && nlcc && compute_stress)
      cout << "<WARNING> EhrenSampleStepper:  stress not yet implemented with non-linear core corrections!  Results WILL be wrong. </WARNING>" << endl;
  
  // use extra memory for SlaterDets if memory variable = normal, large or huge
  if (s_.ctrl.extra_memory >= 3) 
    wf.set_highmem();
  
  if (s_.ctrl.dft_plus_u && oncoutpe)
    cout << "<WARNING> Forces and stress not currently implemented with DFT+U! </WARNING>" << endl;
  
  Timer tm_iter;

  //EWD TDDFT DIFF
  // initialize occupation
  //wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);

  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "TDEULER" )
     wf_stepper = new TDEULERWavefunctionStepper(wf,s_.ctrl.tddt,tmap);
  else if ( wf_dyn == "SOTD" )
  {
     // AS: prepare the array to store the |psi(t-tddt)> for the second-order propagation
     if ( s_.wfv == 0 )
     {
        s_.wfv = new Wavefunction(wf);
        s_.wfv->clear();
        wfv_is_new = true;
     }

     // AS: initially prepare the deque with the respective wave functions
     wfdeque.push_back ( s_.wfv );
     *wfdeque[0]=*(s_.wfv);

     wfdeque.push_back ( s_.wfv );
     *wfdeque[1]=wf;

     wf_stepper = new SOTDWavefunctionStepper(wf,s_.ctrl.tddt,tmap,&wfdeque);
  }
  else if ( wf_dyn == "SORKTD" )
     wf_stepper = new SORKTDWavefunctionStepper(wf,s_.ctrl.tddt,tmap,&wfdeque);
  else if ( wf_dyn == "FORKTD" )
     wf_stepper = new FORKTDWavefunctionStepper(wf,s_.ctrl.tddt,tmap,&wfdeque);
  else
  {
     if ( oncoutpe )
        cout << "<ERROR> EhrenSampleStepper:  undefined wf_dyn = " << wf_dyn << "! </ERROR>" << endl;
     return;
  }
  
  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "SDA" )
    ionic_stepper = new SDAIonicStepper(s_);
  else if ( atoms_dyn == "CG" )
    ionic_stepper = new CGIonicStepper(s_);
  // AS: IMPULSIVE does MD with constant velocities
  else if ( atoms_dyn == "MD" || atoms_dyn == "IMPULSIVE" )
    ionic_stepper = new MDIonicStepper(s_);
  else if ( atoms_dyn == "BMD" )
    ionic_stepper = new BMDIonicStepper(s_);

  if ( ionic_stepper )
    ionic_stepper->setup_constraints();

  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);

  MLWFTransform* mlwft=0;

  if ( compute_mlwf || compute_mlwfc )
  {
    // MLWF can be computed at the gamma point only
    // There must be a single k-point, and it must be gamma
    if ( wf.nkp() > 1 || ( wf.nkp()==1 && wf.kpoint(0) != D3vector(0,0,0) ) )
    {
      if ( oncoutpe )
      {
        cout << "<ERROR> EhrenSampleStepper::step: MLWF can be computed at k=0 only </ERROR>"
             << endl;
        cout << "<ERROR> EhrenSampleStepper::step: cannot run </ERROR>" << endl;
      }
      return;
    }
    if (wf.nspin() > 1 && s_.ctxt_.oncoutpe()) 
      cout << "<ERROR> nspin > 1!  MLWF doesn't currently work with spin-polarized systems </ERROR>" << endl;
    assert(wf.nspin()==1);
    mlwft = new MLWFTransform(*wf.sd(0,0));
  }

  // if ultrasoft, calculate position-dependent functions
  if (ultrasoft) {
     tmap["usfns"].start();
     cd_.update_usfns();
     wf.update_usfns();
     tmap["usfns"].stop();
  }
  // if non-linear core correction defined, calculate position-dependent density
  if (nlcc)
     cd_.update_nlcc();
  
  // calculate time available to avoid exceeding run_timer
  double tbase, tleft;
  bool testtimer = true;
  if (s_.ctrl.run_timer > 0.0) {
    tbase = MPI_Wtime();
    tleft = s_.ctrl.run_timer - (tbase - s_.ctrl.time_init);
  }

  //EWD:  check if hamil_wf needs to be constructed
  if (s_.hamil_wf == 0)
  {
     s_.hamil_wf = new Wavefunction(s_.wf);
     (*s_.hamil_wf) = s_.wf;
     (*s_.hamil_wf).update_occ(0.0,0);
  }
  

#ifdef HPM  
       HPM_Start("iterloop");
#endif
#ifdef TAU
       QB_Pstart(14,scfloop);
#endif
  tmap["total_niter"].start();
  for ( int iter = 0; iter < niter; iter++ )
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
    
    tm_iter.start();
#ifdef USE_APC
    ApcStart(1);
#endif
      
    if ( oncoutpe )
       cout << "<iteration count=\"" << iter+1 << "\">\n";

    if ( ionic_stepper )
       atoms.sync();

    // compute energy and ionic forces using existing wavefunction
      
    tmap["charge"].start();
    cd_.update_density();
    ( ef_.hamil_cd() )->update_density();
    tmap["charge"].stop();

    tmap["efn"].start();
    ef_.update_hamiltonian();
    ef_.update_vhxc();
    tmap["efn"].stop();

    // if symmetry is used, need to calculate set of symmetry-equivalent atoms for 
    // force averaging
    if ( s_.symmetries.nsym() > 0) {
       atoms.findSymmetricAtoms(s_.symmetries);
    }

    // need eigenvalues to compute forces w. ultrasoft
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
        
    tmap["efn"].start();
    double energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
    tmap["efn"].stop();

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

    if ( oncoutpe && iter%s_.ctrl.iprint == 0)
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
        
      
    if ( iter > 0 && ionic_stepper )
    {
       ionic_stepper->compute_v(energy,fion);
    }
    // at this point, positions r0, velocities v0 and forces fion are
    // consistent
    double ekin_ion = 0.0, temp_ion = 0.0;
    if ( ionic_stepper )
    {
       ekin_ion = ionic_stepper->ekin();
       temp_ion = ionic_stepper->temp();
    }

    // print positions, velocities and forces at time t0
    if ( oncoutpe && iter%s_.ctrl.iprint == 0)
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
       cout << "</atomset>" << endl;
       cout << "  <econst> " << energy+ekin_ion << " </econst>\n";
       cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
       cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";
    }
        
    if ( atoms_move )
    {
       if ( s_.constraints.size() > 0 )
       {
          s_.constraints.compute_forces(ionic_stepper->r0(), fion);
          if ( oncoutpe )
          {
             s_.constraints.list_constraints(cout);
          }
       }
       // move atoms to new position: r0 <- r0 + v0*dt + dt2/m * fion
       ionic_stepper->compute_r(energy,fion);
       tmap["efn"].start();
       ef_.atoms_moved();
       tmap["efn"].stop();
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

          tmap["efn"].start();
          ef_.cell_moved();
          ef_.atoms_moved(); // modifications of the cell also move ions
          tmap["efn"].stop();
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
    }

    wf_stepper->preprocess();

    // AS: update the Hamiltonian, the potential, and the energy before propagation
    // starts and/or after the mixing
    tmap["efn"].start();
    ef_.update_hamiltonian();
    ef_.update_vhxc();
    //ef_.update_exc_ehart_eps();
    tmap["efn"].stop();

    // AS: the first step is an EULER one in the case of the second-order propagation
    if ( ( wf_dyn == "SOTD" ) && wfv_is_new ) {
       // AS: |psi(t-tddt)> is the wave function in the very beginning
       *wfdeque[0]=wf;

       // AS: |psi(t)> is calculated by doing one TDEULER step using the EULER stepper
       WavefunctionStepper* wf_init_stepper = new TDEULERWavefunctionStepper(wf,s_.ctrl.tddt,tmap);
               
       tmap["efn"].start();
       cout << wf_dyn << " initialization energy: " << ef_.energy(true,dwf,false,fion,false,sigma_eks) << endl;
       tmap["efn"].stop();
       cout << wf_dyn << " initialization expectation value: " << s_.wf.dot(dwf) << endl;

       wf_init_stepper->update(dwf);
       // AS: now we have |psi(t)>
       *wfdeque[1]=wf;
       (*s_.hamil_wf)=s_.wf;

       // AS: the wave functions used in the Hamiltonian are NOT updated here

       //EWD: has the wf changed since we computed cd above??

       tmap["charge"].start();
       cd_.update_density();
       ( ef_.hamil_cd() )->update_density();
       tmap["charge"].stop();

       // AS: update the Hamiltonian after the first EULER step
       tmap["efn"].start();
       ef_.update_hamiltonian();
       ef_.update_vhxc();
       tmap["efn"].stop();
       // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the potential
       // untouched and only recalculates the energy terms
       //ef_.update_exc_ehart_eps();   
       // AS: we are ready to go with the SOTD
       if ( wf_dyn == "SOTD" ) wfv_is_new=false;
    }
             
    if (wf_dyn=="FORKTD") {
       // AS: initialize wfdeque[0] with 0, it will be k_1
       *wfdeque[0]=dwf;
       // AS: initialize wfdeque[1] with 0, it will be k_2
       *wfdeque[1]=s_.wf;
       // AS: initialize wfdeque[2] with 0, it will be k_3
       *wfdeque[2]=s_.wf;
       // AS: initialize wfdeque[3] with 0, it will be k_4
       *wfdeque[3]=s_.wf;
       // AS: initialize wfdeque[4] with wf, to keep a copy
       *wfdeque[4]=s_.wf;
    }
             
    tmap["efn"].start();
    energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
    tmap["efn"].stop();
    // compute the sum of eigenvalues (with fixed weight)
    // to measure convergence of the subspace update
    // compute trace of the Hamiltonian matrix Y^T H Y
    // scalar product of Y and (HY): tr Y^T (HY) = sum_ij Y_ij (HY)_ij
    // Note: since the hamiltonian is hermitian and dwf=H*wf
    // the dot product in the following line is real

    if (false)
    {
       if (ultrasoft) { 
          const double us_eigenvalue_sum = s_.wf.sdot(dwf);
          if ( oncoutpe )
             cout  << "  <eigenvalue_sum> "
                   << us_eigenvalue_sum << " </eigenvalue_sum>" << endl;
       }
       else {
          const double eigenvalue_sum = s_.wf.dot(dwf);
          if ( oncoutpe )
             cout << "  <eigenvalue_sum> "
                  << eigenvalue_sum << " </eigenvalue_sum>" << endl;
       }
    }
    
    // AS: calculation and output of < psi(t) | H(t_0) | psi(t) > , i.e., the expectation values
    // AS: in the case the wave functions are propagated in time
    // AS: code taken from Wavefunction::diag(Wavefunction& dwf, bool eigvec)
    if ( s_.ctrl.wf_diag == "EIGVAL" )
    {
       if ( oncoutpe ) cout << "<" << wf_dyn << " expectation set>" << endl;

       // AS: this is the new version, not copying code, but copying the wave function instead
       Wavefunction* to_diag_wf1 = new Wavefunction(s_.wf);
       (*to_diag_wf1) = (s_.wf);
       tmap["diag"].start();
       (*to_diag_wf1).diag(dwf,false);
       tmap["diag"].stop();

       if ( oncoutpe )
       {
          for ( int ispin = 0; ispin < (*to_diag_wf1).nspin(); ispin++ )
          {
             for ( int ikp = 0; ikp < (*to_diag_wf1).nkp(); ikp++ )
             {
                const int nst = (*to_diag_wf1).sd(ispin,ikp)->nst();
                const double eVolt = 2.0 * 13.6058;
                cout <<    "  <" << wf_dyn << " expectation values spin=\"" << ispin
                     << "\" kpoint=\"" << s_.wf.sd(ispin,ikp)->kpoint()
                     << "\" n=\"" << nst << "\">" << endl;
                cout << "  ";
                for ( int i = 0; i < nst; i++ )
                {
                   cout << setw(15) << setprecision(8)
                        << (*to_diag_wf1).sd(ispin,ikp)->eig(i)*eVolt;
                   if ( i%5 == 4 ) cout << endl;
                }
                if ( nst%5 != 0 ) cout << endl;
                cout << "  </" << wf_dyn << " expectation values>" << endl;
             }
          }
       }

       if ( oncoutpe ) cout << "</" << wf_dyn << " expectation set>" << endl;
       delete(to_diag_wf1);
    }
             
    // AS: calculation and output of < psi(t) | psi(t) > , i.e., the orthonormalization
    // AS: code adopted from SlaterDet::gram(void)
    if ( s_.ctrl.wf_diag == "EIGVAL" )
    {
       for ( int ispin = 0; ispin < (wf).nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < (wf).nkp(); ikp++ )
          {
             ComplexMatrix ortho(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                                                                                                      
             ortho.gemm('c','n',1.0,(wf).sd(ispin,ikp)->c(),(wf).sd(ispin,ikp)->c(),0.0);
             if ( oncoutpe )
             {
                cout << "ortho: " << endl;
             }
             //ewd:  this is not going to print correctly, as tasks will not print data in order
             cout << ortho;   
          }
       }
    }

    if (wf_dyn=="SORKTD") {
       wfdeque.push_back( s_.wfv );
       // AS: initialize wfdeque[0] with 0, it will be k_1
       (*wfdeque[0]).clear();

       wfdeque.push_back( s_.wfv );
       // AS: initialize wfdeque[1] with 0, it will be k_2
       (*wfdeque[1]).clear();

       wfdeque.push_back( s_.wfv );
       // AS: initialize wfdeque[2] with wf, to keep a copy
       (*wfdeque[2])=s_.wf;

       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             // AS: here we turn wfdeque[0] into k_1
             ((*wfdeque[0]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
             // AS: here we turn wf into |psi(t)+0.5*k_1>
             ((s_.wf).sd(ispin,ikp)->c()).axpy(0.5, (*wfdeque[0]).sd(ispin,ikp)->c() );
          }
       }
                   
       // AS: for the correct output of the nonscf energy
       if (s_.ctrl.non_selfc_energy)
       {
          // AS: the wave functions used in the Hamiltonian are NOT updated here
          tmap["charge"].start();
          cd_.update_density();
          tmap["charge"].stop();
          // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
          // potential untouched and only recalculates the energy terms
          //ef_.update_exc_ehart_eps();
       }
                   
       // AS: update the Hamiltonian with psi(t)+0.5*k_1
       (*s_.hamil_wf)=s_.wf;
       tmap["charge"].start();
       ( ef_.hamil_cd() )->update_density();
       tmap["charge"].stop();
       tmap["efn"].start();
       ef_.update_hamiltonian();
       ef_.update_vhxc();
       // AS: apply the Hamiltonian to |psi(t)+0.5*k_1>
       ef_.energy(true,dwf,false,fion,false,sigma_eks);
       tmap["efn"].stop();

       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             // AS: here we turn wfdeque[1] into k_2
             ((*wfdeque[1]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
          }
       }
    }

    if (wf_dyn=="FORKTD") {
       tmap["forktd_tot"].start();

       // AS: initialize wfdeque[0] with 0, it will be k_1
       (*wfdeque[0]).clear();
       // AS: initialize wfdeque[1] with 0, it will be k_2
       (*wfdeque[1]).clear();
       // AS: initialize wfdeque[2] with 0, it will be k_3
       (*wfdeque[2]).clear();
       // AS: initialize wfdeque[3] with 0, it will be k_4
       (*wfdeque[3]).clear();
       // AS: initialize wfdeque[4] with wf, to keep a copy
       (*wfdeque[4])=s_.wf;

       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             // AS: here we turn wfdeque[0] into k_1
             ((*wfdeque[0]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
             // AS: here we turn wf into |psi(t)+0.5*k_1>
             ((s_.wf).sd(ispin,ikp)->c()).axpy(0.5, (*wfdeque[0]).sd(ispin,ikp)->c() );
          }
       }

       // AS: for the correct output of the nonscf energy
       if (s_.ctrl.non_selfc_energy)
       {
          // AS: the wave functions used in the Hamiltonian are NOT updated here
          tmap["charge"].start();
          cd_.update_density();
          tmap["charge"].stop();

          // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
          // potential untouched and only recalculates the energy terms
          //ef_.update_exc_ehart_eps();
       }

       // AS: set the Hamiltonian accordingly
       (*s_.hamil_wf)=s_.wf;

       tmap["charge"].start();
       ( ef_.hamil_cd() )->update_density();
       tmap["charge"].stop();

       tmap["efn"].start();
       ef_.update_hamiltonian();
       ef_.update_vhxc();
       // AS: apply the Hamiltonian to |psi(t)+0.5*k_1>
       ef_.energy(true,dwf,false,fion,false,sigma_eks);
       tmap["efn"].stop();

       // AS: setting wf back to |psi(t)>
       s_.wf=(*wfdeque[4]);

       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             // AS: here we turn wfdeque[1] into k_2
             ((*wfdeque[1]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
             // AS: here we turn wf into |psi(t)+0.5*k_2>
             ((s_.wf).sd(ispin,ikp)->c()).axpy(0.5, (*wfdeque[1]).sd(ispin,ikp)->c() );
          }
       }
       
       // AS: for the correct output of the energy
       if (s_.ctrl.non_selfc_energy)
       {
          // AS: the wave functions used in the Hamiltonian are NOT updated here
          tmap["charge"].start();
          cd_.update_density();
          tmap["charge"].stop();

          // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
          // potential untouched and only recalculates the energy terms
          //ef_.update_exc_ehart_eps();
       }

       // AS: set the Hamiltonian accordingly
       (*s_.hamil_wf)=s_.wf;
       
       tmap["charge"].start();
       ( ef_.hamil_cd() )->update_density();
       tmap["charge"].stop();

       tmap["efn"].start();
       ef_.update_hamiltonian();
       ef_.update_vhxc();
       // AS: apply the Hamiltonian to |psi(t)+0.5*k_2>
       ef_.energy(true,dwf,false,fion,false,sigma_eks);
       tmap["efn"].stop();

       // AS: setting wf back to |psi(t)>
       s_.wf=(*wfdeque[4]);

       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             // AS: here we turn wfdeque[2] into k_3
             ((*wfdeque[2]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
             // AS: here we turn wf into |psi(t)+k_3>
             ((s_.wf).sd(ispin,ikp)->c()).axpy(1.0, (*wfdeque[2]).sd(ispin,ikp)->c() );
          }
       }
       
       // AS: for the correct output of the energy
       if (s_.ctrl.non_selfc_energy)
       {
          // AS: the wave functions used in the Hamiltonian are NOT updated here
          tmap["charge"].start();
          cd_.update_density();
          tmap["charge"].stop();
          // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
          // potential untouched and only recalculates the energy terms
          //ef_.update_exc_ehart_eps();
       }

       // AS: set the Hamiltonian accordingly
       (*s_.hamil_wf)=s_.wf;

       tmap["charge"].start();
       ( ef_.hamil_cd() )->update_density();
       tmap["charge"].stop();

       tmap["efn"].start();
       ef_.update_hamiltonian();
       ef_.update_vhxc();
       // AS: apply the Hamiltonian to |psi(t)+k_3>
       ef_.energy(true,dwf,false,fion,false,sigma_eks);
       tmap["efn"].stop();

       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             // AS: here we turn wfdeque[3] into k_4
             ((*wfdeque[3]).sd(ispin,ikp)->c()).axpy(-1.0*s_.ctrl.tddt*(complex<double>(0,1)), dwf.sd(ispin,ikp)->c() );
          }
       }
       tmap["forktd_tot"].stop();
    }

    // AS: energy renormalization is currently disabled
    // output of the renormalized energy
    // double wf_dyn_energy;
    // double wf_dyn_eigenvalue_sum;
    // if ( tddft_involved ) {
    // wf_dyn_energy = ef_.energy(true,dwf,false,fion,false,sigma_eks,true);
    // wf_dyn_eigenvalue_sum = real(s_.wf.dot(dwf));
    // if ( oncoutpe )
    // {
    //   cout << wf_dyn << " energy: " << wf_dyn_energy << endl;
    //   cout << wf_dyn << " expectation value: " << wf_dyn_eigenvalue_sum << endl;
    // }
    // }    
                
    wf_stepper->update(dwf);



    // update ultrasoft functions if needed, call gram
    //ewd:  take this out?
    /*
    if (ultrasoft)
       wf.update_usfns();
                
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
    */
    
    
    // AS: change the phase of the wave function if the respective variable is set
    if ( s_.wf.phase_real_set() )
    {
       tmap["phase"].start();
       s_.wf.phase_wf_real();
       tmap["phase"].stop();
    }
    
    if ( oncoutpe )
    {
       cout.setf(ios::fixed,ios::floatfield);
       cout.setf(ios::right,ios::adjustfield);
       cout << "  <etotal_int> " << setw(15) << setprecision(8)
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
    
    if (s_.ctrl.iprint > 0 && iter%s_.ctrl.iprint == 0)
       s_.wf.printocc(); 
    
    // AS: for the correct output of the energy
    if ( s_.ctrl.non_selfc_energy && (wf_dyn!="SORKTD") && (wf_dyn!="FORKTD") )
    {
       // AS: the wave functions used in the Hamiltonian are NOT updated here
       tmap["charge"].start();
       cd_.update_density();
       tmap["charge"].stop();

       // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the
       // potential untouched and only
       // AS: recalculates the energy terms
       tmap["efn"].start();
       ef_.update_exc_ehart_eps();
       tmap["efn"].stop();
    }

    // update occupation numbers if fractionally occupied states
    //EWD:  we shouldn't be doing this for TDDFT, right?
    if ( false && fractional_occ )
    {
       wf.update_occ(s_.ctrl.smearing_width,s_.ctrl.smearing_ngauss);
       s_.wf.printocc();
       const double wf_entropy = wf.entropy();
       if ( oncoutpe )
       {
          cout << "  <!-- Wavefunction entropy: " << wf_entropy << " -->" << endl;
          //const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
          cout << "  <!-- Entropy contribution to free energy: "
               << - wf_entropy * s_.ctrl.smearing_width * 2.0 << " -->" << endl;
       }
    }
          
    // AS: after the first set of non-selfconsistent steps the Hamiltonian is no longer fixed
    // AS: this makes the Hamiltonian time dependent
    (*s_.hamil_wf)=s_.wf;
    tmap["charge"].start();
    ( ef_.hamil_cd() )->update_density();
    tmap["charge"].stop();

    // AS: keep the previous wave function
    if ( wf_dyn == "SOTD" ) *(s_.wfv)=*wfdeque[0];

    // AS: non-adiabatic overlap for succeeding steps of the Born-Oppenheimer trajectory is calculated
    if ( ( atoms_move ) && ( s_.ctrl.na_overlap_min >= 0 ) )
    {
       // AS: after the first step previous_wf is simply set to the current wf
       if ( s_.previous_wf == 0 )
       {
          s_.previous_wf = new Wavefunction(wf);
          (*s_.previous_wf) = s_.wf;
       }
           
       // AS: ATTENTION: currently the min/max indices are not used yet! Overlap is calculated for all the states!
       for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
       {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
             if ( wf.force_complex_set() ) {
                ComplexMatrix overlap1(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                ComplexMatrix overlap2(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                    
                // AS: The following two lines should give something like the d_jk
                // overlap1.gemm('c','n',-1.0*(complex<double>(0,1))/2.0/s_.ctrl.dt,(*s_.previous_wf).sd(ispin,ikp)->c(),wf.sd(ispin,ikp)->c(),0.0);
                // overlap1.gemm('c','n',1.0*(complex<double>(0,1))/2.0/s_.ctrl.dt,wf.sd(ispin,ikp)->c(),(*s_.previous_wf).sd(ispin,ikp)->c(),1.0);
                    
                // AS: The following two lines are the same stuff without the complex i
                overlap1.gemm('c','n',-1.0/2.0/s_.ctrl.dt,(*s_.previous_wf).sd(ispin,ikp)->c(),wf.sd(ispin,ikp)->c(),0.0);
                overlap2.gemm('c','n',1.0/2.0/s_.ctrl.dt,wf.sd(ispin,ikp)->c(),(*s_.previous_wf).sd(ispin,ikp)->c(),0.0);

                // AS: this is the normalization
                // overlap1.gemm('c','n',1.0,wf.sd(ispin,ikp)->c(),wf.sd(ispin,ikp)->c(),0.0);
                    
                cout << "OVERLAP complex: " << endl;
                cout << "overlap1" << overlap1 << endl;
                cout << "overlap2" << overlap2 << endl;
                overlap1 += overlap2;
                cout << "overlap1+overlap2" << overlap1 << endl;
             }
             else
             {
                DoubleMatrix overlap1(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());
                DoubleMatrix overlap2(wf.sd(ispin,ikp)->context(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).nb(),(wf.sd(ispin,ikp)->c()).nb());

                // DoubleMatrix proxies
                DoubleMatrix previous_wf_proxy( (*s_.previous_wf).sd(ispin,ikp)->c() );
                DoubleMatrix wf_proxy( wf.sd(ispin,ikp)->c() );

                // AS: The following two lines are something like the d_jk without the complex i
                // AS: 2*(-1.0/2.0/s_.ctrl.dt) is the factor
                overlap1.gemm('c','n',-1.0/s_.ctrl.dt,previous_wf_proxy,wf_proxy,0.0);
                // AS: -1.0*(-1.0/2.0/s_.ctrl.dt) is the factor
                overlap1.ger(1.0/2.0/s_.ctrl.dt,previous_wf_proxy,0,wf_proxy,0);
                overlap2.gemm('c','n',1.0/s_.ctrl.dt,wf_proxy,previous_wf_proxy,0.0);
                overlap2.ger(-1.0/2.0/s_.ctrl.dt,wf_proxy,0,previous_wf_proxy,0);

                // AS: the following three lines reproduce the correct normalization
                // overlap1.gemm('c','n',2.0,wf_proxy,wf_proxy,0.0);
                // rank-1 update using first row of sdc_proxy() and c_proxy
                // overlap1.ger(-1.0,wf_proxy,0,wf_proxy,0);

                cout << "OVERLAP real: " << endl;
                cout << "overlap1" << overlap1 << endl;
                cout << "overlap2" << overlap2 << endl;
                overlap1 += overlap2;
                cout << "overlap1+overlap2" << overlap1 << endl;
             }
          }
       }

       // AS: keep a copy of the wave function to calculate the overlap during the following iteration
       delete s_.previous_wf;
       s_.previous_wf = new Wavefunction(s_.wf);
       (*s_.previous_wf) = s_.wf;
    }

    if ( compute_mlwf || compute_mlwfc )
    {
       SlaterDet& sd = *(wf.sd(0,0));
       mlwft->compute_transform();

       if ( compute_mlwf )
          mlwft->apply_transform(sd);
          
          if ( oncoutpe )
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
    }

    wf_stepper->postprocess();
 
#ifdef USE_APC
    ApcStop(1);
#endif
    // print iteration time
    double time = tm_iter.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( oncoutpe )
    {
       cout << left << setw(34) << "<timing where=\"run\""
            << setw(24) << " name=\" iteration\""
            << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
            << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
            << endl;
       cout << "</iteration>" << endl;
    }

    s_.ctrl.mditer++;

    // if savedenfreq variable set, save density in VMD Cube format
    if (s_.ctrl.savedenfreq > 0)
    {
       if (s_.ctrl.mditer%s_.ctrl.savedenfreq == 0 || s_.ctrl.mditer == 1)
       {
          string filebase = "density.";
          ostringstream oss;
          oss.width(7);  oss.fill('0');  oss << s_.ctrl.mditer;
          string denfilename = filebase + oss.str() + ".cube";
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

    // if savefreq variable set, checkpoint
    if (s_.ctrl.savefreq > 0)
    {
       if (s_.ctrl.savefreq == 1 || (s_.ctrl.mditer > 0 && s_.ctrl.mditer%s_.ctrl.savefreq == 0) )
       {
          // create output directory if it doesn't exist
          string dirbase = "md.";
          string filebase = "mdchk";
          ostringstream oss;
          oss.width(7);  oss.fill('0');  oss << s_.ctrl.mditer;
          string dirstr = dirbase + oss.str();
          string format = "binary";
          if ( s_.ctxt_.mype()==0 )
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
          if ( s_.ctxt_.mype()==0 )
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
          
          if (s_.ctrl.tddft_involved)
          {
             // write s_.hamil_wf
             string hamwffile = filestr + "hamwf";
             if ( s_.ctxt_.mype()==0 )
                cout << "<!-- MDSaveCmd:  wf write finished, writing hamil_wf to " << hamwffile << "... -->" << endl;
             s_.hamil_wf->write_states(hamwffile,format);
          }
          
       }
    }

    
    if ( atoms_move )
       s_.constraints.update_constraints(dt);

  } // for iter

  tmap["total_niter"].stop();
#ifdef TAU  
  QB_Pstop(scfloop);
#endif
#ifdef HPM
  HPM_Stop("iterloop");
#endif

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
  // print timer map
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
     double time = (*i).second.real();
     double tmin = time;
     double tmax = time;
     s_.ctxt_.dmin(1,1,&tmin,1);
     s_.ctxt_.dmax(1,1,&tmax,1);
     if ( s_.ctxt_.mype()==0 )
     {
       cout << left << setw(34) << "<timing where=\"run\""
           << setw(8) << " name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
     }
  }
#endif  
  
  if ( atoms_move )
  {

    // need eigenvalues to compute forces w. ultrasoft
    if (ultrasoft) { 
      ef_.energy(true,dwf,false,fion,false,sigma_eks);
      tmap["post_diag"].start();
      //s_.wf.diag(dwf,compute_eigvec);
      s_.wf.diag(dwf,true);
      tmap["post_diag"].stop();
      s_.wf.printeig();
    }

    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["post_charge"].start();
    cd_.update_density();
    tmap["post_charge"].stop();

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
  }

  // compute ionic forces at last position to update velocities
  // consistently with last position
  tmap["post_charge"].start();
  cd_.update_density();
  tmap["post_charge"].stop();

  ef_.update_vhxc();
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

  if (wf_dyn != "SOTD")
  {
     // delete wavefunction velocities
     // AS: only if we don't need them anymore; for SOTD we do need them later on
     if ( s_.wfv != 0 )
        delete s_.wfv;
     s_.wfv = 0;
  }

  delete mlwft;

  // delete steppers
  delete wf_stepper;
  delete ionic_stepper;
  delete cell_stepper;

}

////////////////////////////////////////////////////////////////////////////////
void EhrenSampleStepper::get_forces(vector<vector<double> >& f) const {
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
double EhrenSampleStepper::get_energy(string ename) {
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
valarray<double> EhrenSampleStepper::get_stress(string sname) {
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
