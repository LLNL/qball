////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.h,v 1.4 2009/10/07 18:53:22 draeger1 Exp $

#ifndef SAMPLESTEPPER_H
#define SAMPLESTEPPER_H

#include "Sample.h"
#include "Timer.h"
#include <map>
#include <string>
#include <valarray>
using namespace std;

typedef map<string,Timer> TimerMap;

class SampleStepper
{
  protected:
  
  Sample& s_;
  AtomSet& atoms_;
  int                       nsp_;
  int                       ndofs_;
  vector<int>               na_;      // number of atoms per species na_[nsp_]
  vector<double>            pmass_;   // pmass_[nsp_]
  
  vector<vector<double> > fion;
  valarray<double> sigma_eks, sigma_kin, sigma_ext, sigma;

  // Do not allow construction of SampleStepper unrelated to a Sample
  SampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  virtual void step(int niter) = 0;
  virtual void get_forces(vector<vector<double> > &f) const = 0;
  virtual double get_energy(string ename) = 0;
  virtual valarray<double> get_stress(string ename) = 0;
  void print_stress(void);
  void compute_sigma(void); // compute kinetic contribution to stress
  virtual void initialize_density() {}
  
  SampleStepper(Sample& s);
  virtual ~SampleStepper(void);
};
#endif
