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
// JDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: JDWavefunctionStepper.h,v 1.2 2009/09/08 05:36:49 fgygi Exp $

#ifndef JDWAVEFUNCTIONSTEPPER_H
#define JDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
class Preconditioner;
class EnergyFunctional;

class JDWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  EnergyFunctional& ef_;
  Wavefunction wft_, dwft_;

  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void) {}

  JDWavefunctionStepper(Wavefunction& wf, Preconditioner& p,
                        EnergyFunctional& ef, TimerMap& tmap);
  ~JDWavefunctionStepper() {};
};
#endif
