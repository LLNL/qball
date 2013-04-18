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
// Wavefunction.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "D3vector.h"
#include "UnitCell.h"
#include "AtomSet.h"
#include "SharedFilePtr.h"
#include <vector>
#if USE_CSTDIO_LFS
#include <cstdio>
#endif
using namespace std;

class SlaterDet;
class Context;

class Wavefunction {
  private:

  const Context& ctxt_;

  int nel_;           // number of electrons
  int nempty_;        // number of empty states
  int nspin_;         // number of spins (1 or 2)
  int deltaspin_;     // number of spin excitations
  
  int nrowmax_;        // maximum number of rows of a spincontext
  int nparallelkpts_;  // maximum number of kpoints to run in parallel
  int nkptloc_;        // number of kpoints of local processor
  vector<int> nkptloc_list_;  // full list of nkpts on each subcontext
  vector<int> kptloc_; // index of kpoint in kpoint_ array on local processor
  vector<int> mysdctxt_;  // index of sdcontext corresponding to global k-point index kp
  bool kpt_added_;     // if user adds kpoint with Kpoint command, replace default gamma
  vector<vector<int> > kptproc0_;    // pe number of proc 0 for spin i, sdcontext j
  vector<vector<int> > nkplocproc0_; // number of local kpoints for spin i, sdcontext j
  int spinloc_;        // index of local spin
  bool ultrasoft_;     // use overlap matrix S for diagonalization, orthogonalization
  bool force_complex_wf_; // AS: force complex basis [even if kpoint == (0,0,0)]
                          // AS: necessary for time propagation of wave functions
  bool wf_phase_real_; // AS: change phase of the wave function to make it real for Gamma only
  
  UnitCell cell_ ;    // unit cell
  UnitCell refcell_ ; // reference cell
  double   ecut_ ;    // energy cutoff
  
  vector<double>      weight_;  // weight[ikp]
  vector<D3vector>    kpoint_;  // kpoint[ikp]
  double              weightsum_;
  
  vector<int> nst_;                       // nst_[ispin]
  Context* wfcontext_;                    // contains spincontext(s)
  vector<Context*> spincontext_;          // spincontext[ispin]
  vector<vector<Context*> > sdcontext_;   // sdcontext_[ispin][ikp]
  vector<vector<SlaterDet*> > sd_;        // sd[ispin][ikp]
  vector<vector<Context*> > sdcontextsq_;   // sdcontextsq_[ispin][ikp]
  bool reshape_context_;
  
  bool hasdata_;   // wait to allocate until a load, randomize or run command
  void allocate(); // create contexts and allocate SlaterDet's 
  void deallocate();
  void compute_nst();
  void resize(); // resize SlaterDets if ecut,cell,refcell,or nst have changed
  void reshape(); // reshape SlaterDets onto new parallel distribution while
                  // preserving existing data
  
  public:
  
  Wavefunction(const Context& ctxt);
  Wavefunction(const Wavefunction& wf);
  ~Wavefunction();
  Wavefunction& operator=(const Wavefunction& wf);
  
  const Context& context(void) const { return ctxt_; }
  const UnitCell& cell(void) const { return cell_; }
  const UnitCell& refcell(void) const { return refcell_; }
  const D3vector kpoint(int ikp) const { return kpoint_[ikp]; }
  double weight(int ikp) const { return weight_[ikp]; }
  double weightsum(void) const { return weightsum_; }
  double ecut(void) const { return ecut_; }
  SlaterDet* sd(int ispin, int ikp) const { return sd_[ispin][ikp]; }
  const Context* wfcontext(void) const { return wfcontext_; }
  const Context* sdcontext(int ispin, int ikp) const 
    { return sdcontext_[ispin][ikp]; }
  const Context* sdcontextsq(int ispin, int ikp) { return sdcontextsq_[ispin][ikp]; }
  const Context* spincontext(int ispin) const 
    { return spincontext_[ispin]; }
  int nkp(void) const;            // number of k points
  int nel(void) const;            // total number of electrons
  int nst(int ispin) const;       // number of states of spin ispin
  int nst(void) const;            // number of states
  int nempty(void) const;         // number of empty states
  int nspin(void) const;          // number of spins
  int deltaspin(void) const { return deltaspin_; } // number of spin excitations
  int nrowmax(void) const { return nrowmax_; }
  int nparallelkpts(void) const { return nparallelkpts_; }
  int sdcontextsize(int ispin) const { return sdcontext_[ispin].size(); }
  int nkptloc(void) const { return nkptloc_; }
  int nkptloc(int sdnum) const { return nkptloc_list_[sdnum]; }
  int mysdctxt(int kp) const { return mysdctxt_[kp]; }
  int kptloc(int k) const { return kptloc_[k]; }
  bool kptactive(int k) const;
  bool spinactive(int ispin) const;
  int spinloc(void) const { return spinloc_; }
  SlaterDet* sdloc(int k) const { return sd_[spinloc_][kptloc_[k]]; }
  SlaterDet* sdloc(int ispin, int k) const { return sd_[ispin][kptloc_[k]]; }
  int kptproc0(int ispin, int kloc) const { return kptproc0_[ispin][kloc]; }
  int nkplocproc0(int ispin, int kloc) const { return nkplocproc0_[ispin][kloc]; }
  double spin(void) const;        // total spin
  
  bool hasdata(void) { return hasdata_; }
  void set_hasdata(bool hasd);
  void resize(const UnitCell& cell, const UnitCell& refcell, double ecut);
  void resize(double ec) { resize(cell_,refcell_,ec); }
  void set_ecut(double ecut) { ecut_ = ecut; }
  void set_cell(const UnitCell& cell) { cell_ = cell; }
  void set_refcell(const UnitCell& refcell) { refcell_ = refcell; }
  void reset(void); // initialize with lowest plane waves
  void clear(void); // initialize with zero
  void set_nel(int nel);
  void set_nempty(int nempty);
  void set_nspin(int nspin);
  void set_deltaspin(int deltaspin);
  void set_nrowmax(int n);
  void set_nparallelkpts(int n);
  void add_kpoint(D3vector kpoint, double weight);
  void del_kpoint(D3vector kpoint);
  void set_reshape_context(bool reshape);
  void set_ultrasoft(bool us);
  bool ultrasoft(void) { return ultrasoft_; }
  void init_usfns(AtomSet* atoms);
  void update_usfns();
  void calc_spsi();
  void set_highmem(void);
  
  void randomize(double amplitude, bool highmem);
  void randomize_us(double amplitude, AtomSet& as, bool highmem);
  void randomize_real(double amplitude);
  // AS: shift state n_state by the vector (shift_x, shift_y, shift_z)
  void shift_wf(double shift_x,double shift_y,double shift_z,int n_state);
  // AS: change phase of the wave function to make it real for Gamma only
  void phase_wf_real(void);

  void rescale(double factor);
  
  void update_occ(double smearingwidth, int ngauss);
  double entropy(void) const; // dimensionless entropy
  void gram(void);
  void riccati(Wavefunction& wf);
  void align(Wavefunction& wf);
  void diag(Wavefunction& dwf, bool eigvec);
  void extrap_real(const double dt, const AtomSet& as);
  
  double dot(const Wavefunction& wf) const;
  double sdot(const Wavefunction& wf) const;
  
  void print(ostream& os, string encoding, string tag) const;
  void print_casino(ostream& os) const;
  void print_vmd(string filebase, const AtomSet& as) const;
  void printeig(void);
  void printocc(void);
  void write(SharedFilePtr& fh, string encoding, string tag) const;
  void write_dump(string filebase);
  void write_fast(string filebase);
  void write_states(string filebase, string format);
  void write_states_old(string filebase, string format);
  void write_mditer(string filebase, int mditer);
  void read_dump(string filebase);
  void read_fast(string filebase);
  void read_states(string filebase);
  void read_states_old(string filebase);
  void read_mditer(string filebase, int& mditer);
  void info(ostream& os, string tag);
  // AS: is true when wave function has to be forced to complex also for kpoint == (0,0,0)
  bool force_complex_set(void) const;
  // AS: enable or disable forcing of complex wave functions
  void force_complex(bool new_force_complex_wf);
  // AS: is true when the wave function is made real for Gamma only
  bool phase_real_set(void) const;
  // AS: change phase of the wave function to make it real for Gamma only
  void phase_real(bool new_wf_phase_real);
};
ostream& operator << ( ostream& os, Wavefunction& wf );
#endif
