////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SimpleConvergenceDetector.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

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

