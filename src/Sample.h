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
// Sample.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SAMPLE_H
#define SAMPLE_H

#include "AtomSet.h"
#include "Wavefunction.h"
#include "Control.h"
#include "ConstraintSet.h"
#include "SymmetrySet.h"
#include <vector>
#include <complex>

class Context;

class Sample {
  private:
  
  public:
  
  const Context& ctxt_;

  AtomSet atoms;
  ConstraintSet constraints;
  Wavefunction wf;
  // AS: hamil_wf is a pointer which by default points to wf.
  // AS: in order to construct an Hamiltonian from different wave functions (charge densities)
  // AS: than wf, make hamil_wf point to these wave functions instead
  Wavefunction* hamil_wf;
  // AS: keep a copy of the wave function during Born-Oppenheimer MD when non-adiabatic overlaps are to be calculated
  Wavefunction* previous_wf;
  Wavefunction* wfv; // wavefunction velocity
  Control ctrl;
  SymmetrySet symmetries;
  vector<vector<complex<double> > > rhog_last; // previous charge density (to avoid discontinuity in restart)

  Sample(const Context& ctxt) : ctxt_(ctxt), atoms(ctxt), wf(ctxt), wfv(0), hamil_wf(0),
      symmetries(ctxt), constraints(ctxt) { ctrl.sigmas = 0.5; ctrl.facs = 2.0; }
  ~Sample(void) { delete wfv; }
  void reset(void)
  {
    atoms.reset();
    constraints.reset();
    //extforces.reset();
    wf.reset();
    delete wfv;
  }

};

#endif
