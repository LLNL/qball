////////////////////////////////////////////////////////////////////////////////
//
// EnthalpyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnthalpyFunctional.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef ENTHALPYFUNCTIONAL_H
#define ENTHALPYFUNCTIONAL_H

#include "ChargeDensity.h"
#include <vector>
#include <valarray>
#include <complex>
using namespace std;

class Basis;
class FourierTransform;

class EnthalpyFunctional {
  private:
  
  const Context& ctxt_;  
  const ChargeDensity& cd_;
  
  double epv_;
  double evol_;
  int np012loc_;
  double pressure_;
  double threshold_;
  double sigma_;
  int nsigma_;
  double stepdx_;
  vector<double> stepfn_;

  FourierTransform& vft_;
  Basis& vbasis_;

  public:

  EnthalpyFunctional(const ChargeDensity& cd, double pressure, double threshold);
  ~EnthalpyFunctional();
  void update(vector<vector<double> >& vr);
  double epv(void) { return epv_; }
  double evol(void) { return evol_; }
};

#endif
