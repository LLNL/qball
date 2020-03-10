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
// EhrenSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef EHRENSAMPLESTEPPER_H
#define EHRENSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
#include <ui/PlotCmd.h>
#include <deque>
#include "CurrentDensity.h"

class WavefunctionStepper;
class IonicStepper;
using namespace std;

class EhrenSampleStepper : public SampleStepper
{
  private:
  
  Wavefunction dwf;
  Wavefunction* wfv;
  deque<Wavefunction*> wfdeque;
  VectorPotential* tempvp_;
  
  int nitscf_;
  int nite_;
  ChargeDensity cd_;
  CurrentDensity currd_;
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

// Local Variables:
// mode: c++
// End:
