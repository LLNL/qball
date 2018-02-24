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
// SDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "SDIonicStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SDIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  // Steepest descent step
  for ( int is = 0; is < r0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = r0_[is][i] + dt2bym * f0[is][i];
    }
  }
  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
}
