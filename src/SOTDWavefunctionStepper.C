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
// SOTDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SOTDWavefunctionStepper.C,v 1.8 2011-05-03 15:56:19 schleife Exp $

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
