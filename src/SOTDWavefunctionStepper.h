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
// SOTDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SOTDWavefunctionStepper.h,v 1.5 2011-05-03 15:56:19 schleife Exp $

#include <config.h>

#ifndef SOTDWAVEFUNCTIONSTEPPER_H
#define SOTDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

#include <deque>
using namespace std;

// AS: this class provides the SOTD wave function stepper
// AS: SOTD enables second-order wave-function propagation according to |psi(t+tddt)> = |psi(t-tddt)> - 2*i*tddt*|H psi(t)>

class SOTDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  protected:

  deque<Wavefunction*> *wfdeque_ ;
  // AS: the following block worked without dequeue
  // Wavefunction *prev_wf_;

  public:

  void update(Wavefunction& dwf);

  SOTDWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, deque<Wavefunction*> *wfdeque);
  // AS: the following block worked without dequeue
  //SOTDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap, Wavefunction* prev_wf);
  ~SOTDWavefunctionStepper() {};
};
#endif

// Local Variables:
// mode: c++
// End:
