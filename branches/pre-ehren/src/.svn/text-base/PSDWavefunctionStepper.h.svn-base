////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.h,v 1.2 2009/03/25 22:30:34 draeger1 Exp $

#ifndef PSDWAVEFUNCTIONSTEPPER_H
#define PSDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
class Preconditioner;

class PSDWavefunctionStepper : public WavefunctionStepper
{
  private:
  
  Preconditioner& prec_;

  public:

  void update(Wavefunction& dwf);
  
  PSDWavefunctionStepper(Wavefunction& wf, Preconditioner& p, TimerMap& tmap);
  ~PSDWavefunctionStepper() {};
};
#endif
