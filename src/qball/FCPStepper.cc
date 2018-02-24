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
// FCPStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FCPStepper.C,v 1.23 2010-04-16 22:41:55 fgygi Exp $

#include <config.h>

#include "FCPStepper.h"
#include <stdlib.h>
using namespace std;


////////////////////////////////////////////////////////////////////////////////
void FCPStepper::compute_r(double f0)
{
  // f0 contains forces at r0
  // Compute new positions rp using the velocity Verlet algorithm
  // enforce constraints for rp
  // update rm <- r0, r0 <- rp, and update atomset

  // compute rp
  const double dt2by2m = dt_ * dt_ / ( 2.0 * pmass_ );
  double rp = r0_ + v0_ * dt_ + dt2by2m * f0;
  rm_ = r0_;
  r0_ = rp;
  s_.wf.set_deltacharge( r0_ );
}

////////////////////////////////////////////////////////////////////////////////
void FCPStepper::compute_v(double f0)
{
  // compute velocities v0_ using r0_, rm_ and f0(r0)
  // enforce constraints for vp
  // adjust velocities with the thermostat
  // compute ekin

  assert(dt_ != 0.0);
  compute_ekin();

  const double dtby2m = dt_ / ( 2.0 * pmass_ );
  const double vhalf = ( r0_ - rm_ ) / dt_;
  v0_ = vhalf + dtby2m * f0;

  eta_ = tanh ( ( temp() - th_temp_ ) / th_width_ ) / th_time_;
  if ( s_.ctxt_.onpe0() )
  {
    cout << "  fcp thermostat: temp=" << temp() << endl;
    cout << "  fcp thermostat: tref=" << th_temp_ << endl;
    cout << "  fcp thermostat: eta=" << eta_ << endl;
  }

  const double fac = 1.0 - eta_ * fabs(dt_);
  v0_ *= fac;
}

////////////////////////////////////////////////////////////////////////////////
void FCPStepper::compute_ekin(void)
{
  const double v = v0_;
  ekin_ = 0.5 * pmass_ * v * v;
}
