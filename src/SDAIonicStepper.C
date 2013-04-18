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
// SDAIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "SDAIonicStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SDAIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  double fp0;
  bool wolfe1, wolfe2;

  // check if the largest component of the force f0 is larger than max_force
  const double max_force = 0.1;
  double largest_force = 0.0;
  for ( int is = 0; is < r0_.size(); is++ )
  {
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      largest_force = max(largest_force,fabs(f0[is][i]));
    }
  }

  if ( largest_force > max_force )
  {
    if ( s_.ctxt_.oncoutpe() )
      cout << "  SDAIonicStepper: force exceeds limit, taking SD step " << endl;
    // take a steepest descent step with limited displacement and exit
    const double alpha_sd = max_force/largest_force;
    // SD step
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        rp_[is][i] = r0_[is][i] + alpha_sd * f0[is][i];
      }
    }
    constraints_.enforce_r(r0_,rp_);
    rm_ = r0_;
    r0_ = rp_;
    atoms_.set_positions(r0_);
    // reset the SDA algorithm
    first_step_ = true;
    return;
  }

  // SDA algorithm

  if ( !first_step_ )
  {
    wolfe1 = e0 < ec_ + fpc_ * sigma1_ * alpha_;
    // fp0 = -proj(f0,pc)
    fp0 = 0.0;
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        fp0 -= f0[is][i] * pc_[is][i];
      }
    }
    wolfe2 = fabs(fp0) < sigma2_ * fabs(fpc_);
    if ( s_.ctxt_.oncoutpe() )
    {
      cout << "  SDAIonicStepper: fpc = " << fpc_ << endl;
      cout << "  SDAIonicStepper: fp0 = " << fp0 << endl;
      cout << "  SDAIonicStepper: ec = " << ec_ << " e0 = " << e0 <<  endl;
      cout << "  SDAIonicStepper: ec_ + fpc_ * sigma1_ * alpha_ ="
           << ec_ + fpc_ * sigma1_ * alpha_ << endl;
      cout << "  SDAIonicStepper: wolfe1/wolfe2 = "
           << wolfe1 << "/" << wolfe2 << endl;
    }
  }

  if ( first_step_ || (wolfe1 && wolfe2) )
  {
    // set new descent direction
    // pc = f0
    fc_ = f0;
    pc_ = fc_;
    // fpc = d_e / d_alpha in direction pc
    fpc_ = 0.0;
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        fpc_ -= fc_[is][i] * pc_[is][i];
      }
    }
    ec_ = e0;
    rc_ = r0_;
    fp0 = fpc_;
    // reset line minimizer
    linmin_.reset();
  }

  alpha_ = linmin_.newalpha(alpha_,e0,fp0);

  if ( s_.ctxt_.oncoutpe() )
    cout << "  SDAIonicStepper: alpha = " << alpha_ << endl;

  // rp = rc + alpha * pc
  for ( int is = 0; is < r0_.size(); is++ )
  {
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = rc_[is][i] + alpha_ * pc_[is][i];
    }
  }

  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);

  first_step_ = false;
}

