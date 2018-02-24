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
// TDEULERWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDEULERWavefunctionStepper.C,v 1.8 2011-03-31 15:56:19 schleife Exp $

#include <config.h>

#include "TDEULERWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

// AS: this class provides the TDEULER wave function stepper

////////////////////////////////////////////////////////////////////////////////
TDEULERWavefunctionStepper::TDEULERWavefunctionStepper(Wavefunction& wf, double tddt,
  TimerMap& tmap) :
  tddt_(tddt), WavefunctionStepper(wf,tmap)
{}

////////////////////////////////////////////////////////////////////////////////
void TDEULERWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["sd_update_wf"].start();
      // AS: TDEULER enables first-order wave-function propagation according to |psi(t+tddt)> = |psi(t)> - i*tddt*|H psi(t)>
      wf_.sd(ispin,ikp)->c().axpy(-tddt_*(complex<double>(0,1)),dwf.sd(ispin,ikp)->c());
      // AS: if there was an energy renormalization to prevent time propagation from blowing up
      // AS: it could be undone here; this is just an example:
      // wf_.sd(ispin,ikp)->c() *= exp( complex<double>(0,1) * tddt_ * (2.263) );
      // AS: for the Euler integrator it does not change anything, possibly we need it later
      tmap_["sd_update_wf"].stop();
    }
  }
}
