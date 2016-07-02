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
// EnthalpyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

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
