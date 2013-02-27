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
