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
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////

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
