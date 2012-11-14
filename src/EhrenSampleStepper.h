////////////////////////////////////////////////////////////////////////////////
//
// EhrenSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EhrenSampleStepper.h,v 1.4 2009/10/07 18:53:22 draeger1 Exp $

#ifndef EHRENSAMPLESTEPPER_H
#define EHRENSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
#include "PlotCmd.h"
#include <deque>
class WavefunctionStepper;
class IonicStepper;
using namespace std;

class EhrenSampleStepper : public SampleStepper
{
  private:
  
  Wavefunction dwf;
  Wavefunction* wfv;
  deque<Wavefunction*> wfdeque;
  
  int nitscf_;
  int nite_;
  ChargeDensity cd_;
  EnergyFunctional ef_;
  
  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  bool tddft_involved_;
  
  // Do not allow construction of EhrenSampleStepper unrelated to a Sample
  EhrenSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);
  
  void get_forces(vector<vector<double> > &f) const;
  double get_energy(string ename);
  valarray<double> get_stress(string sname);

  EhrenSampleStepper(Sample& s, int nitscf, int nite);
  ~EhrenSampleStepper();
};
#endif
