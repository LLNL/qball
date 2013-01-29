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
