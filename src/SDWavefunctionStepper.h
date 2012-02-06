////////////////////////////////////////////////////////////////////////////////
//
// SDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.h,v 1.2 2009/03/25 22:30:34 draeger1 Exp $

#ifndef SDWAVEFUNCTIONSTEPPER_H
#define SDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class SDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double alpha_;

  public:

  void update(Wavefunction& dwf);

  SDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap);
  ~SDWavefunctionStepper() {};
};
#endif
