////////////////////////////////////////////////////////////////////////////////
//
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDAWavefunctionStepper.h,v 1.5 2009/09/25 23:18:11 draeger1 Exp $

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
class Preconditioner;

class PSDAWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  Wavefunction wf_last_, dwf_last_;

  // Anderson acceleration flag
  int nkp_;
  int nspin_;
  vector<vector<bool> > extrapolate_;
  
  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void)
  {
    for (int i=0; i<nspin_; i++)
      for (int k=0; k<nkp_; k++)
        extrapolate_[i][k] = false;
  }

  PSDAWavefunctionStepper(Wavefunction& wf, Preconditioner& p, TimerMap& tmap);
  ~PSDAWavefunctionStepper() {};
};
#endif
