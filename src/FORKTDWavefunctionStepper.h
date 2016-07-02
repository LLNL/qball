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
// FORKTDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FORKTDWavefunctionStepper.h,v 1.5 2011-06-02 15:56:19 schleife Exp $

#include <config.h>

#ifndef FORKTDWAVEFUNCTIONSTEPPER_H
#define FORKTDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

#include <deque>
using namespace std;

// AS: this class provides the FORKTD wave function stepper
// AS: FORKTD enables fourth-order wave-function propagation of:
// AS: d |psi(t)> / dt = - (i / h) H[psi(t)] |psi(t)>
// AS: according to:
// AS: k_1 = - dt (i / h) H[psi(t)] |psi(t)>
// AS: k_2 = - dt (i / h) H[psi(t)+0.5*k_1] |psi(t)+0.5*k_1>
// AS: k_3 = - dt (i / h) H[psi(t)+0.5*k_2] |psi(t)+0.5*k_2>
// AS: k_4 = - dt (i / h) H[psi(t)+k_3] |psi(t)+k_3>
// AS: |psi(t+dt)> = |psi(t)> + 1/6 * k_1 + 1/3 * k_2 + 1/3 * k_3 + 1/6 * k_4

class FORKTDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  protected:

  deque<Wavefunction*> *wfdeque_ ;

  public:

  void update(Wavefunction& dwf);

  FORKTDWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, deque<Wavefunction*> *wfdeque);
  ~FORKTDWavefunctionStepper() {};
};
#endif
