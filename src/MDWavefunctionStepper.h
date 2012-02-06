////////////////////////////////////////////////////////////////////////////////
//
// MDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDWavefunctionStepper.h,v 1.3 2009/03/27 00:53:24 draeger1 Exp $

#ifndef MDWAVEFUNCTIONSTEPPER_H
#define MDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class MDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double dt_;
  double dt2bye_;
  Wavefunction *wfv_;

  double ekin_ep_, ekin_em_;
  double ekin_eh(void);

  public:

  void update(Wavefunction& dwf);
  void compute_wfm(Wavefunction& dwf);
  void compute_wfv(Wavefunction& dwf);
  double ekin(void) const { return 0.5*(ekin_ep_ + ekin_em_); }

  MDWavefunctionStepper(Wavefunction& wf, Wavefunction* wfv,
    double dt, double dt2bye, TimerMap& tmap);
  ~MDWavefunctionStepper() {};
};
#endif
