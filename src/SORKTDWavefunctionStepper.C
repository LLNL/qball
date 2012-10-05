////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SORKTDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SORKTDWavefunctionStepper.C,v 1.8 2011-06-02 15:56:19 schleife Exp $

#include "SORKTDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
using namespace std;

// AS: this class provides the SORKTD wave function stepper
// AS: SORKTD enables fourth-order wave-function propagation of:
// AS: d |psi(t)> / dt = - (i / h) H[psi(t)] |psi(t)>
// AS: according to:
// AS: k_1 = - dt (i / h) H[psi(t)] |psi(t)>
// AS: k_2 = - dt (i / h) H[psi(t)+0.5*k_1] |psi(t)+0.5*k_1>
// AS: |psi(t+dt)> = |psi(t)> + k_2
// AS:             = |psi(t)> - dt (i / h) H[psi(t) - 0.5 * dt (i / h) H[psi(t)] |psi(t)>] |psi(t) - 0.5 * dt (i / h) H[psi(t)] |psi(t)>

////////////////////////////////////////////////////////////////////////////////
SORKTDWavefunctionStepper::SORKTDWavefunctionStepper(Wavefunction& wf, double tddt,
 TimerMap& tmap, deque<Wavefunction*> *wfdeque) : tddt_(tddt), WavefunctionStepper(wf,tmap), wfdeque_(wfdeque)
{ }

////////////////////////////////////////////////////////////////////////////////
void SORKTDWavefunctionStepper::update(Wavefunction& dwf)
{
  wf_ = (*(*wfdeque_)[2]);

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["sd_update_wf"].start();

      (wf_.sd(ispin,ikp)->c()).axpy( 1.0, (*(*wfdeque_)[1]).sd(ispin,ikp)->c() );

      tmap_["sd_update_wf"].stop();
    }
  }

  // AS: trash the old deque
  (*wfdeque_).pop_front();
  (*wfdeque_).pop_front();
  (*wfdeque_).pop_front();
}
