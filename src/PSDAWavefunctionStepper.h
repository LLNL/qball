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
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
class Preconditioner;

class PSDAWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  Wavefunction wf_last_, dwf_last_;

  // Anderson acceleration flag
  int nkp_;
  int nspin_;
  vector<vector<bool> > extrapolate_;
  
  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void)
  {
    for (int i=0; i<nspin_; i++)
      for (int k=0; k<nkp_; k++)
        extrapolate_[i][k] = false;
  }

  PSDAWavefunctionStepper(Wavefunction& wf, Preconditioner& p, TimerMap& tmap);
  ~PSDAWavefunctionStepper() {};
};
#endif

// Local Variables:
// mode: c++
// End:
