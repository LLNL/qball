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
// EnergyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef ENERGYFUNCTIONAL_H
#define ENERGYFUNCTIONAL_H

#include <complex>
#include <vector>
#include <valarray>
#include <map>
#include <string>
#include "ChargeDensity.h"
#include "StructureFactor.h"
#include "SelfConsistentPotential.h"
#include "Timer.h"
using namespace std;

class Sample;
class Basis;
class AtomSet;
class Wavefunction;
class UnitCell;
class FourierTransform;
class XCPotential;
class NonLocalPotential;
class ConfinementPotential;
class EnthalpyFunctional;
class HubbardPotential;

typedef map<string,Timer> TimerMap;

class EnergyFunctional
{
  private:

   friend class SelfConsistentPotential;
   
  const Sample& s_;
  const Wavefunction& wf_;
  ChargeDensity& cd_;
  Basis* vbasis_;
  FourierTransform *vft;
  vector<vector<FourierTransform*> > ft;
  StructureFactor sf;
  XCPotential* xcp;
  EnthalpyFunctional* epvf;
  vector<vector<NonLocalPotential*> > nlp;
  vector<vector<ConfinementPotential*> > cfp; // cfp[ispin][ikp]
  HubbardPotential* hubp_;
  ChargeDensity *hamil_cd_;   // AS: specifies the density used for setting up the Hamiltonian
  vector<complex<double> > hamil_rhoelg, hamil_rhogt;   // AS: specifies the density used for setting up the Hamiltonian
  
  vector<vector<double> > vps, dvps, rhops;
  vector<complex<double> > tmp_r, vion_local_g, dvion_local_g, vlocal_g,
      rhopst, rhogt, rhoelg, vtemp;
  vector<double> ftmp;
  
  vector<vector<double> > tau0, taum, fion_esr;
  vector<double> zv_, rcps_;
  vector<int> na_;
  int namax_;
  int nsp_;
  double ekin_, econf_, eps_, enl_, ehub_, ehart_, 
      ecoul_, exc_, esr_, eself_, ets_, epv_, eexf_, etotal_;
  double eharris_;  // terms for Harris-Foulkes estimate for convergence detection
  
  valarray<double> sigma_ekin,sigma_econf,sigma_eps,sigma_ehart,sigma_exc,
    sigma_enl, sigma_esr, sigma;

  public:

  vector<vector<double> > v_r;
  vector<vector<complex<double> > > vxc_g;
  vector<vector<complex<double> > > veff_g;
  mutable TimerMap tmap;
  
  ChargeDensity *hamil_cd() { return hamil_cd_; } ;   // AS: specifies the density used for setting up the Hamiltonian

  double energy(Wavefunction& wf, bool compute_hpsi, Wavefunction& dwf,
    bool compute_forces, vector<vector<double> >& fion,
                bool compute_stress, valarray<double>& sigma);
  
  double etotal(void) const { return etotal_; }
  double ekin(void) const { return ekin_; }
  double econf(void) const { return econf_; }
  double eps(void) const { return eps_; }
  double enl(void) const { return enl_; }
  double ehart(void) const { return ehart_; }
  double ecoul(void) const { return ecoul_; }
  double exc(void) const { return exc_; }
  double esr(void) const { return esr_; }
  double eself(void) const { return eself_; }
  double ets(void) const { return ets_; }
  double epv(void) const { return epv_; }
  double eexf(void) const { return eexf_; }
  double ehub(void) const { return ehub_; }
  double eharris(void) const { return eharris_; }
  double etotal_harris(void) const { return (ekin_ + econf_ + enl_ + ets_ + epv_ + ehub_ + eharris_ + esr_ - eself_); };  
  double casino_ewald(void);
  double casino_vloc(void);
  const ConfinementPotential *confpot(int ispin, int ikp) const { return cfp[ispin][ikp]; }
  
  void update_vhxc(void);
  void update_harris(void);

  // AS: modified version of EnergyFunctional::update_vhxc(void) which leaves the potential untouched and only
  // AS: recalculates the energy terms
  void update_exc_ehart_eps(void);
  // AS: update the Hamiltonian in the case of TDDFT (i.e. selfconsistent) time propagation
  void update_hamiltonian(void);
  
  void atoms_moved(void);
  void cell_moved(const bool compute_stress);
  
  void print(ostream& os) const;
  void print_memory(ostream&os, double& totsum, double& locsum) const;
  void print_timing();

  SelfConsistentPotential get_self_consistent_potential() const
  {
     return SelfConsistentPotential(*this);
  }

  void set_self_consistent_potential(const SelfConsistentPotential & potential)
  {
     v_r = potential.v_r;
     hamil_rhoelg = potential.hamil_rhoelg;
     rhoelg = potential.rhoelg;
     eps_ = potential.eps_;
     ehart_ = potential.ehart_;
     exc_ = potential.exc_;
     esr_ = potential.esr_;
  }
  
  EnergyFunctional(const Sample& s, const Wavefunction& wf, ChargeDensity& cd);
  ~EnergyFunctional();
};
ostream& operator << ( ostream& os, const EnergyFunctional& e );
#endif
