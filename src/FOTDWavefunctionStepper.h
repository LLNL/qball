////////////////////////////////////////////////////////////////////////////////
//
// FOTDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FOTDWavefunctionStepper.h,v 1.5 2011-06-02 15:56:19 schleife Exp $

#ifndef FOTDWAVEFUNCTIONSTEPPER_H
#define FOTDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

#include <deque>
using namespace std;

// AS: this class provides the FOTD wave function stepper
// AS: FOTD enables fourth-order wave-function propagation according to
// AS: |psi(t+2tddt)> = 8|psi(t+tddt)> + 12*i*tddt*|H psi(t)> - 8|psi(t-tddt)> + |psi(t-2tddt)> 

class FOTDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  protected:

  deque<Wavefunction> *wfdeque_ ;
  // AS: the following block worked without dequeue
  // Wavefunction *next_wf_;
  // Wavefunction *wf1_;
  // Wavefunction *wf2_;

  public:

  void update(Wavefunction& dwf);

  FOTDWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, deque<Wavefunction> *wfdeque);
  // AS: the following block worked without dequeue
  // FOTDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap, Wavefunction* next_wf, Wavefunction* wf1, Wavefunction* wf2);
  ~FOTDWavefunctionStepper() {};
};
#endif
