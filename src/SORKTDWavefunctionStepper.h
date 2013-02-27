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
