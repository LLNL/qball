////////////////////////////////////////////////////////////////////////////////
//
// SimpleConvergenceDetector.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SimpleConvergenceDetector.C,v 1.2 2007/10/24 19:05:23 draeger1 Exp $

#include <cmath>
#include "SimpleConvergenceDetector.h"
#include "ConvergenceDetector.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
SimpleConvergenceDetector::SimpleConvergenceDetector(int nsteps, double threshold) : 
  nsteps_(nsteps) {
  threshold_ = threshold;
}

////////////////////////////////////////////////////////////////////////////////
void SimpleConvergenceDetector::addValue(double val) {

  lastval_ = val;
  if (recent_history_.size() >= nsteps_) 
    recent_history_.pop_front();

  recent_history_.push_back(val);
}

////////////////////////////////////////////////////////////////////////////////
bool SimpleConvergenceDetector::isConverged() {
  bool converge = false;
  if (recent_history_.size() >= nsteps_ && threshold_ > 0.0) {
    double max = lastval_;
    double min = lastval_;
    for (list<double>::iterator li = recent_history_.begin(); li != recent_history_.end(); li++) {
      if ((*li) > max) max = (*li);
      if ((*li) < min) min = (*li);
    }
    if (fabs(max-min) < threshold_) 
      converge = true;
  }
  return converge;
}

