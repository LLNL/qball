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
// TDEULERWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDEULERWavefunctionStepper.h,v 1.5 2011-03-31 15:56:19 schleife Exp $

#ifndef TDEULERWAVEFUNCTIONSTEPPER_H
#define TDEULERWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

// AS: this class provides the TDEULER wave function stepper
// AS: TDEULER enables first-order wave-function propagation according to |psi(t+tddt)> = |psi(t)> - i*tddt*|H psi(t)>

class TDEULERWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  public:

  void update(Wavefunction& dwf);

  TDEULERWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap);
  ~TDEULERWavefunctionStepper() {};
};
#endif
