////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.h,v 1.4 2009/10/07 18:53:22 draeger1 Exp $

#ifndef BOSAMPLESTEPPER_H
#define BOSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
class WavefunctionStepper;
class IonicStepper;
using namespace std;

class BOSampleStepper : public SampleStepper
{
  private:
  
  Wavefunction dwf;
  Wavefunction* wfv;
  int nitscf_;
  int nite_;
  ChargeDensity cd_;
  EnergyFunctional ef_;
  
  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  bool initial_atomic_density;
  
  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);
  void initialize_density(void);
  
  void get_forces(vector<vector<double> > &f) const;
  double get_energy(string ename);
  valarray<double> get_stress(string sname);

  BOSampleStepper(Sample& s, int nitscf, int nite);
  ~BOSampleStepper();
};
#endif
