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
// PPCGWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PPCGWAVEFUNCTIONSTEPPER_H
#define PPCGWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
#include "HPsi.h"
class Preconditioner;
class AtomSet;
class ChargeDensity;

class PPCGWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  HPsi hpsi_;
  Wavefunction reswf_,hreswf_;
  int maxit_;
  
  int nkp_;
  int nspin_;

// E.V. new members begin
  int  sbsize_;  
  int  qrstep_;     
  void cholqr(DoubleMatrix&, DoubleMatrix&);
  void update_wfp0(DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, int, int, int);     
  void update_wfp(DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, DoubleMatrix&, int, int, int);     
// E.V. new members end
  
  public:

  void cell_moved(void);
  void update(Wavefunction& dwf);
  virtual void preprocess(void)
  {
  }

  PPCGWavefunctionStepper(Wavefunction& wf, Preconditioner& p, const int maxit,
                            const AtomSet& atoms, const ChargeDensity& cd_,
                            vector<vector<double> >& v_r, TimerMap& tmap, 
                            int sbsize, int qrstep);
  ~PPCGWavefunctionStepper() {};
};
#endif
