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
// RMMDIISWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "RMMDIISWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Preconditioner.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
RMMDIISWavefunctionStepper::RMMDIISWavefunctionStepper(Wavefunction& wf, 
  Preconditioner& p, TimerMap& tmap) : 
  WavefunctionStepper(wf,tmap), prec_(p)
{}

void RMMDIISWavefunctionStepper::preprocess(void) {
  iter_ = 0;
  std::cout << "etotal PREPROCESS " << iter_ << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
void RMMDIISWavefunctionStepper::update(Wavefunction& dwf) {
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          assert(wf_.sd(ispin,ikp)->basis().real());
	  assert(! wf_.ultrasoft());

	  std::cout << "etotal update " << iter_ << std::endl;
	  
          tmap_["rmmdiis_residual"].start();

	  // proxy real matrices c, cp
	  DoubleMatrix c(wf_.sd(ispin,ikp)->c());
	  DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
          
	  DoubleMatrix a(c.context(), c.n(), c.n(), c.nb(), c.nb());
	  
	  // factor 2.0 in next line: G and -G
	  a.gemm('t', 'n', 2.0, c, cp, 0.0);
	  
	  // rank-1 update correction
	  a.ger(-1.0, c, 0, cp, 0);

	  // cp = cp - c * a
	  cp.gemm('n','n',-1.0,c,a,1.0);
	
	  tmap_["rmmdiis_residual"].stop();
	  
          // dwf.sd->c() now contains the descent direction (HV-VA)

	  prec_.apply(*dwf.sd(ispin, ikp), ispin, ikp);

	  wf_.sd(ispin,ikp)->c().axpy(-1.0, dwf.sd(ispin, ikp)->c());

        }
      }
    }
  }

  iter_++;
}

void RMMDIISWavefunctionStepper::postprocess(void) {
  std::cout << "etotal PostPROCESS " << iter_ << std::endl;
}
