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
// BMDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "BMDIonicStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void BMDIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  // f0 contains forces at r0
  // e0 contains energy at r0
  // Compute new positions rp using the velocity Verlet algorithm
  // enforce constraints for rp
  // update rm <- r0, r0 <- rp, and update atomset

  if ( e0 > em_ )
  {
    // energy has increased, reset r0 to rm and make a steepest descent step
    r0_ = rm_;
    for ( int is = 0; is < r0_.size(); is++ )
    {
      const double dt2by2m = dt_ * dt_ / ( 2.0 * pmass_[is] );
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        rp_[is][i] = rm_[is][i] + dt2by2m * fm_[is][i];
      }
    }
    constraints_.enforce_r(r0_,rp_);
    r0_ = rp_;
    e0_ = em_;
    atoms_.set_positions(r0_);
  }
  else
  {
    // compute rp
    for ( int is = 0; is < r0_.size(); is++ )
    {
      const double dt2by2m = dt_ * dt_ / ( 2.0 * pmass_[is] );
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        rp_[is][i] = r0_[is][i] + v0_[is][i] * dt_ + dt2by2m * f0[is][i];
      }
    }
    constraints_.enforce_r(r0_,rp_);
    rm_ = r0_;
    fm_ = f0;
    em_ = e0_;
    r0_ = rp_;
    atoms_.set_positions(r0_);
  }
}

////////////////////////////////////////////////////////////////////////////////
void BMDIonicStepper::compute_v(double e0, const vector<vector< double> >& f0)
{
  // compute velocities v0_ using r0_, rm_ and f0(r0)
  // enforce constraints for vp
  // adjust velocities with the thermostat

  e0_ = e0;
  assert(dt_ != 0.0);
  for ( int is = 0; is < v0_.size(); is++ )
  {
    assert(pmass_[is] > 0.0);
    const double dtby2m = dt_ / ( 2.0 * pmass_[is] );
    for ( int i = 0; i < v0_[is].size(); i++ )
    {
      const double vhalf = ( r0_[is][i] - rm_[is][i] ) / dt_;
      v0_[is][i] = vhalf + dtby2m * f0[is][i];
    }
  }

  // check if energy increased
  if ( e0_ > em_ )
  {
    // reset the velocity to zero
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        v0_[is][i] = 0.0;
      }
    }
  }
  else
  {
    // accelerate
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        v0_[is][i] *= 1.1;
      }
    }
  }
  compute_ekin();
}

////////////////////////////////////////////////////////////////////////////////
void BMDIonicStepper::compute_ekin(void)
{
  ekin_ = 0.0;
  for ( int is = 0; is < v0_.size(); is++ )
  {
    for ( int i = 0; i < v0_[is].size(); i++ )
    {
      const double v = v0_[is][i];
      ekin_ += 0.5 * pmass_[is] * v * v;
    }
  }
}
