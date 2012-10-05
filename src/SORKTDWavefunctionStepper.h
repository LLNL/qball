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
// SORKTDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SORKTDWavefunctionStepper.h,v 1.5 2011-06-02 15:56:19 schleife Exp $

#ifndef SORKTDWAVEFUNCTIONSTEPPER_H
#define SORKTDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

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

class SORKTDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  protected:

  deque<Wavefunction*> *wfdeque_ ;

  public:

  void update(Wavefunction& dwf);

  SORKTDWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, deque<Wavefunction*> *wfdeque);
  ~SORKTDWavefunctionStepper() {};
};
#endif
