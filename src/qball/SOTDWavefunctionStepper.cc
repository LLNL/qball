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
// SOTDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SOTDWavefunctionStepper.C,v 1.8 2011-05-03 15:56:19 schleife Exp $

#include <config.h>

#include "SOTDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
using namespace std;

// AS: this class provides the SOTD wave function stepper
// AS: SOTD enables second-order wave-function propagation according to |psi(t+tddt)> = |psi(t-tddt)> - 2*i*tddt*|H psi(t)>

////////////////////////////////////////////////////////////////////////////////
SOTDWavefunctionStepper::SOTDWavefunctionStepper(Wavefunction& wf, double tddt,
 TimerMap& tmap, deque<Wavefunction*> *wfdeque) : tddt_(tddt), WavefunctionStepper(wf,tmap), wfdeque_(wfdeque)
{ }

////////////////////////////////////////////////////////////////////////////////
void SOTDWavefunctionStepper::update(Wavefunction& dwf)
{
  // AS: create some space for the new wave function
  (*wfdeque_).push_back( &wf_ );

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["sd_update_wf"].start();

      (*(*wfdeque_)[2]).sd(ispin,ikp)->c() = (*(*wfdeque_)[0]).sd(ispin,ikp)->c();
      ((*(*wfdeque_)[2]).sd(ispin,ikp)->c()).axpy(-2.0*tddt_*(complex<double>(0,1)),dwf.sd(ispin,ikp)->c());

      tmap_["sd_update_wf"].stop();
    }
  }

  // AS: trash the oldest wave function
  (*wfdeque_).pop_front();

  // AS: copy back to wf, because EnergyFunctional is applied based on wf
  wf_=*(*wfdeque_)[1];
}
