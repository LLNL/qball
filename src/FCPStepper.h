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
// FCPStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FCPStepper.h,v 1.13 2010-04-16 22:41:55 fgygi Exp $

#include <config.h>

//
// FCPStepper is used in the following way
//
// input: r0,v0
//
// compute energy e0(r0) and forces f0(r0)
// for ( k=0; k<n; k++ )
// {
//   // r0,v0,e0,f0 known
//   stepper->compute_r(e0,f0)
//   {
//     computes rp using r0, v0 and f0
//     restores constraints on rp using rp and r0
//     updates rm<-r0, r0<-rp and update atomset positions
//   }
//
//   compute f0(r0)
//
//   stepper->compute_v(e0,f0)
//   {
//     computes v0 using r0,rm,f0
//     restores constraints on v0 using r0, v0
//     modifies velocities using the thermostat (rescaling)
//     updates atomset velocities
//   }
// }
// r0,v0,f0 consistent at this point
//

#ifndef FCPSTEPPER_H
#define FCPSTEPPER_H

#include "Sample.h"

class FCPStepper
{
  private:

  double th_temp_;
  double th_time_;
  double th_width_;
  double pmass_;
  double v0_;
  double r0_;
  double rm_;
  double eta_;
  double ekin_;
  Sample& s_;
  double dt_;
  void compute_ekin(void);

  public:

  FCPStepper(Sample& s) : s_(s), dt_(s.ctrl.dt)
  {
    th_temp_  = s.ctrl.fcp_th_temp;
    th_time_  = s.ctrl.fcp_th_time;
    th_width_ = s.ctrl.fcp_th_width;
    if( th_temp_  < 0.0 ) th_temp_  = s.ctrl.th_temp;
    if( th_time_  < 0.0 ) th_time_  = s.ctrl.th_time;
    if( th_width_ < 0.0 ) th_width_ = s.ctrl.th_width;
    pmass_ = s.ctrl.fcp_pmass;
    v0_    = 0.0;
    r0_    = s.wf.deltacharge();
    rm_    = 0.0;
    eta_   = 0.0;
    ekin_  = 0.0;
    compute_ekin();
  }

  void compute_r(double f0);
  void compute_v(double f0);
  double temp(void) const
  {
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    return 2.0 * ( ekin_ / boltz );
  }
};

#endif
