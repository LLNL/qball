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
// RunCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "RunCmd.h"
#include<iostream>
using namespace std;
#include <qball/BOSampleStepper.h>
#include <qball/CPSampleStepper.h>
#include <qball/EhrenSampleStepper.h>

#include<ctime>
#include<cassert>

int RunCmd::action(int argc, char **argv)
{

  if ( argc < 2 || argc > 5)
  {
    if ( ui->oncoutpe() )
    {
      cout << " use: run [-atomic_density] niter" << endl;
      cout << "      run [-atomic_density] niter nitscf" << endl;
      cout << "      run [-atomic_density] niter nitscf nite" << endl;
      return 1;
    }
  }

  if (s->ctrl.timer_hit) {
    if ( ui->oncoutpe() )
      cout << " <!-- RunCmd: run_timer exceeded in previous run, skipping current run command. -->" << endl;
    return 0;
  }
  
  if ( s->wf.nst() == 0 )
  {
    if ( ui->oncoutpe() )
      cout << " <!-- RunCmd: no states, cannot run -->" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->oncoutpe() )
      cout << " <!-- RunCmd: ecut = 0.0, cannot run -->" << endl;
    return 1;
  }
  
  SampleStepper* stepper;
  
  int iarg = 1;
  bool atomic_density = false;
  if ( !strcmp(argv[iarg],"-atomic_density") )
  {
    atomic_density = true;
    iarg++;
    argc--;
  }

  int niter = atoi(argv[iarg]);
  int nite = 0;
  int nitscf = 1;
  if ( argc == 3 )
  {
    // run niter nitscf
    nitscf = atoi(argv[iarg+1]);
  }
  else if ( argc == 4 )
  {
    // run niter nitscf nite
    nitscf = atoi(argv[iarg+1]);
    nite = atoi(argv[iarg+2]);
  }

  if (!s->wf.hasdata()) {
    if ( !atomic_density )
    {
       if (ui->oncoutpe())
          cout << "<INFO> Wavefunction has been neither loaded nor initialized! Calling randomize_wf... </INFO>" << endl;
       //ewd:  add this to avoid users running without initializing wf
       //ewd:  note:  this can cause problems if users set ecut before calling .sys file
       double amp = 0.02;
       bool highmem = false;
       if (s->ctrl.extra_memory >= 3)
          highmem = true;
       if (s->ctrl.ultrasoft)
          s->wf.randomize_us(amp,s->atoms,highmem);
       else
          s->wf.randomize(amp,highmem);
    }
  }

  // check if hamil_wf needs to be constructed
  if (s->ctrl.tddft_involved && s->hamil_wf == 0)
  {
     if (ui->oncoutpe())
        cout << "<INFO> s->hamil_wf has not been initialized, calling constructor. </INFO>" << endl;
     s->hamil_wf = new Wavefunction(s->wf);
     (*s->hamil_wf) = s->wf;
     (*s->hamil_wf).update_occ(0.0,0);
  }
  
  if ( s->ctrl.wf_dyn == "MD" )
    stepper = new CPSampleStepper(*s);
  else if (s->ctrl.wf_dyn == "TDEULER" || s->ctrl.wf_dyn == "SOTD" || s->ctrl.wf_dyn == "SORKTD" || s->ctrl.wf_dyn == "FORKTD" || s->ctrl.wf_dyn == "ETRS" || s->ctrl.wf_dyn == "AETRS")
    stepper = new EhrenSampleStepper(*s,nitscf,nite);
  else
    stepper = new BOSampleStepper(*s,nitscf,nite);
  
  assert(stepper!=0);
  
  if (ui->oncoutpe() )
    cout << "<run niter_ionic=\"" << niter << "\" niter_scf=\"" << nitscf << "\" niter_nonscf=\"" << nite << "\">" << endl;

  if ( atomic_density )
    stepper->initialize_density();

  s->wf.info(cout,"wavefunction");
  stepper->step(niter);
  if (ui->oncoutpe() )
    cout << "</run>" << endl;
  
  delete stepper;
  
  return 0;
}
