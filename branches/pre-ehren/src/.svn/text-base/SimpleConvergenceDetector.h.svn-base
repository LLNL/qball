////////////////////////////////////////////////////////////////////////////////
//
// SimpleConvergenceDetector.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SimpleConvergenceDetector.h,v 1.2 2007/10/24 19:05:23 draeger1 Exp $

#ifndef SIMPLECONVERGENCEDETECTOR_H
#define SIMPLECONVERGENCEDETECTOR_H

#include <list>
#include "ConvergenceDetector.h"
using namespace std;

class SimpleConvergenceDetector : public ConvergenceDetector {
  private:

  int nsteps_;
  list<double> recent_history_; // store last nsteps values
  double lastval_;

  public:

  SimpleConvergenceDetector(int nsteps, double threshold);

  double threshold(void) { return threshold_; };
  int nsteps(void) {return nsteps_; };

  void addValue(double val);
  bool isConverged(void);
  ~SimpleConvergenceDetector() {}
    
};
#endif
