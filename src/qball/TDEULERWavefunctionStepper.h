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
// TDEULERWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDEULERWavefunctionStepper.h,v 1.5 2011-03-31 15:56:19 schleife Exp $

#include <config.h>

#ifndef TDEULERWAVEFUNCTIONSTEPPER_H
#define TDEULERWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

// AS: this class provides the TDEULER wave function stepper
// AS: TDEULER enables first-order wave-function propagation according to |psi(t+tddt)> = |psi(t)> - i*tddt*|H psi(t)>

class TDEULERWavefunctionStepper : public WavefunctionStepper
{
  private:

  double tddt_;

  public:

  void update(Wavefunction& dwf);

  TDEULERWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap);
  ~TDEULERWavefunctionStepper() {};
};
#endif

// Local Variables:
// mode: c++
// End:
