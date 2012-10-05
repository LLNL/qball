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
// SOTDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SOTDWavefunctionStepper.h,v 1.5 2011-05-03 15:56:19 schleife Exp $

#ifndef SOTDWAVEFUNCTIONSTEPPER_H
#define SOTDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

#include <deque>
using namespace std;

// AS: this class provides the SOTD wave function stepper
// AS: SOTD enables second-order wave-function propagation according to |psi(t+tddt)> = |psi(t-tddt)> - 2*i*tddt*|H psi(t)>

class SOTDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  protected:

  deque<Wavefunction*> *wfdeque_ ;
  // AS: the following block worked without dequeue
  // Wavefunction *prev_wf_;

  public:

  void update(Wavefunction& dwf);

  SOTDWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, deque<Wavefunction*> *wfdeque);
  // AS: the following block worked without dequeue
  //SOTDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap, Wavefunction* prev_wf);
  ~SOTDWavefunctionStepper() {};
};
#endif
