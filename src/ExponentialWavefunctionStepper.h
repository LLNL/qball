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
// ExponentialWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ExponentialWavefunctionStepper.h,v 1.5 2011-06-02 15:56:19 schleife Exp $

#ifndef EXPONENTIALWAVEFUNCTIONSTEPPER_H
#define EXPONENTIALWAVEFUNCTIONSTEPPER_H

#include "EnergyFunctional.h"
#include "SelfConsistentPotential.h"
#include "Wavefunction.h"
#include "WavefunctionStepper.h"

#include <deque>
using namespace std;

class ExponentialWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;
  int order_;
  int stored_iter_;
  bool approximated_;
  bool merge_exp_;
  std::vector<SelfConsistentPotential> potential_;
  Wavefunction expwf_;
  Wavefunction wfhalf_;
  Wavefunction newwf_; 

  protected:

  EnergyFunctional & ef_;
  Sample & s_;
  void exponential(int num_exp, double dt1, double dt2, Wavefunction * dwf = 0);

  public:
  void preupdate();
  void update(Wavefunction& dwf);

  ExponentialWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s, bool approximated);
  ~ExponentialWavefunctionStepper() {};
};
#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
