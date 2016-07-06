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
// JDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: JDWavefunctionStepper.h,v 1.2 2009/09/08 05:36:49 fgygi Exp $

#include <config.h>

#ifndef JDWAVEFUNCTIONSTEPPER_H
#define JDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
class Preconditioner;
class EnergyFunctional;

class JDWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  EnergyFunctional& ef_;
  Wavefunction wft_, dwft_;

  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void) {}

  JDWavefunctionStepper(Wavefunction& wf, Preconditioner& p,
                        EnergyFunctional& ef, TimerMap& tmap);
  ~JDWavefunctionStepper() {};
};
#endif
