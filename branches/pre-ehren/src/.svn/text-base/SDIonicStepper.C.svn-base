////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.C,v 1.5 2009/03/27 00:53:24 draeger1 Exp $

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
