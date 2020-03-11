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
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H
#include "Timer.h"
#include <map>
#include <string>

typedef std::map<std::string,Timer> TimerMap;
class Wavefunction;

class WavefunctionStepper
{
  private:
  
  protected:
  Wavefunction& wf_;
  TimerMap& tmap_;
  
  public:

  virtual void preupdate(void) {}
  //virtual void update(Wavefunction& dwf) = 0;
  virtual void update(Wavefunction& dwf,int time = 0.0) = 0;
  virtual void preprocess(void) {}
  virtual void postprocess(void) {}

  WavefunctionStepper(Wavefunction& wf, TimerMap& tmap) : 
  wf_(wf), tmap_(tmap)
  {}
  virtual ~WavefunctionStepper() {}
};
#endif

// Local Variables:
// mode: c++
// End:
