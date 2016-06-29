////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// BMDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BMDIonicStepper.C,v 1.1 2010/01/16 01:26:35 draeger1 Exp $

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
