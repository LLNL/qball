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
// Preconditioner.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "Preconditioner.h"
#include "EnergyFunctional.h"
#include "Wavefunction.h"
#include "Sample.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ConfinementPotential.h"

////////////////////////////////////////////////////////////////////////////////
Preconditioner::Preconditioner(const Sample& s, const Wavefunction& wf, const EnergyFunctional& ef) : 
    s_(s), wf_(wf), ef_(ef)
{
  update();
}

////////////////////////////////////////////////////////////////////////////////
void Preconditioner::update(void)
{
  // reinitialize preconditioner
  bool use_confinement = s_.ctrl.ecuts > 0.0;
  // If ecutprec is zero, use ecut
  const double ecutpr = s_.ctrl.ecutprec> 0.0 ? s_.ctrl.ecutprec : wf_.ecut();
  
  diag_.resize(wf_.nspin());
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      diag_[ispin].resize(wf_.nkp());
      for ( int ikp=0; ikp<wf_.nkp(); ikp++) {
        if (wf_.kptactive(ikp)) {
          assert(wf_.sd(ispin,ikp) != 0);
          // Only resize and initialize diag_ if ikp is active on this task
          const Basis& basis = wf_.sd(ispin,ikp)->basis();
          const int ngwloc = basis.localsize();
          diag_[ispin][ikp].resize(ngwloc);
          const double *kpg2_ptr = basis.kpg2_ptr();
          
          if ( use_confinement ) {
            const valarray<double>& fstress = ef_.confpot(ispin,ikp)->fstress();
            for ( int ig = 0; ig < ngwloc; ig++ ) {
              double e = 0.5 * ( kpg2_ptr[ig] + fstress[ig] );
              diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;
            }
          }
          else {
            for ( int ig = 0; ig < ngwloc; ig++ ) {
              double e = 0.5 * kpg2_ptr[ig];
              diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;
            }
          }
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
void Preconditioner::apply(SlaterDet & sd, int ispin, int ikp){
          
  const valarray<double>& precdiag = diag(ispin, ikp);
  
  double* dcoeff =  (double*) sd.c().cvalptr();
  const int mloc = sd.c().mloc();
  const int ngwl = sd.basis().localsize();
  const int nloc = sd.c().nloc();
  
  for ( int n = 0; n < nloc; n++ ) {
    // note: double mloc length for complex<double> indices
    double* dc = &dcoeff[2*mloc*n];
    
    // loop to ngwl only since diag[i] is not defined on [0:mloc-1]
    for ( int i = 0; i < ngwl; i++ ) {
      double fac = precdiag[i];
      dc[2*i]     *= fac;
      dc[2*i + 1] *= fac;
    }
  }
 
}
