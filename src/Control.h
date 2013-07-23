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
// Control.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <vector>

struct Control
{
  // control variables
  string debug, timing;
  string wf_dyn, atoms_dyn; // dynamics string flags 
  int nite;
  double emass;       // electron mass
    
  string smearing;        // method to use for smearing (fermi or gaussian)
  double smearing_width;  // smearing width (default = 0.0, i.e. no smearing)
  int smearing_ngauss;    // order of Methfessel-Paxton expansion
  double ecutprec;

  string wf_diag;
  string wf_extrap;
  // AS: control the calculation of the total energy during non-selfconsistent electronic steps
  bool non_selfc_energy; 

  string tcp;
  double tcp_rcut;
  double tcp_sigma;
  
  double gms_mix; // mixing factor for generalized minimum spread functions
  
  string thermostat;
  double th_temp,th_time, th_width; // thermostat control
  string center_of_mass;  // can be "fixed", "free", or "unset"
  string hugoniostat;    // "ON" or "OFF"
  double hug_etot;       // reference energy
  double hug_volume;     // reference volume
  double hug_pressure;   // reference pressure
  double hug_deltatemp;  // size of initial temperature step
  int hug_freq;          // how often (ionic steps) to update temperature
  
  string stress;
  string cell_dyn;
  string cell_lock;
  double cell_mass;
  int cell_stepfreq;  // only move cell every cell_stepfreq ionic steps
  double ecuts,sigmas,facs; // confinement energy parameters
  double ext_stress[6]; // external stress tensor: xx,yy,zz,xy,yz,xz
  
  string xc;
  string spin;
  int delta_spin;
  int nkpoints; // number of kpoints; wait to allocate Wavefunction until all are input

  double dt;
  double tddt; // AS: time step for the wave function propagation
  int na_overlap_min; // AS: minimum band index for the calculation of non-adiabatic overlaps
  int na_overlap_max; // AS: maximum band index for the calculation of non-adiabatic overlaps
  int iprint;
  int timeout;
  double threshold_scf;      // energy threshold for electronic iteration loops
  int threshold_scf_nsteps;  // number of scf steps over which threshold is evaluated
  double threshold_force,threshold_stress;  // thresholds for ionic iteration loops
  int threshold_force_nsteps,threshold_stress_nsteps;  // number of ionic steps over which threshold is evaluated
  
  string charge_mixing;
  double charge_mix_coeff;
  double charge_mix_rcut;
  int charge_mix_ndim;
    
  double enthalpy_pressure;
  double enthalpy_threshold;

  bool reshape_context;

  double run_timer;  // total maximum run time in seconds
  double time_init;
  bool timer_hit;
  bool timer_mdsavecmd;
  bool timer_savecmd;
  bool timer_savesyscmd;
    
  bool tddft_involved;
  bool dft_plus_u;
  bool ultrasoft;
  bool nlcc;         // non-linear core correction
  double ecutden;

  int extra_memory;  // guides use of extra memory to speed computation
  int mditer; // store global iteration count to help with checkpointing

  int savefreq;     // if > 0, checkpoint within iteration loop
  int savedenfreq;  // if > 0, checkpoint within iteration loop
};
#endif
