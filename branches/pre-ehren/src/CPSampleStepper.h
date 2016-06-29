////////////////////////////////////////////////////////////////////////////////
//
// CPSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CPSampleStepper.h,v 1.3 2008/07/16 00:14:18 draeger1 Exp $

#ifndef CPSAMPLESTEPPER_H
#define CPSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "ChargeDensity.h"
#include "Sample.h"
#include "Wavefunction.h"
class MDWavefunctionStepper;
class MDIonicStepper;
using namespace std;

class CPSampleStepper : public SampleStepper
{
  private:
  
  ChargeDensity cd_;
  EnergyFunctional ef_;
  Wavefunction dwf;
  Wavefunction* wfv;
  
  MDWavefunctionStepper* mdwf_stepper;
  MDIonicStepper* mdionic_stepper;

  // Do not allow construction of CPSampleStepper unrelated to a Sample
  CPSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);
  void get_forces(vector<vector<double> > &f) const;
  double get_energy(string ename);
  valarray<double> get_stress(string sname);

  CPSampleStepper(Sample& s);
  ~CPSampleStepper();
};
#endif
