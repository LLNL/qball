////////////////////////////////////////////////////////////////////////////////
//
// SDAIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDAIonicStepper.C,v 1.5 2009/03/25 22:30:34 draeger1 Exp $

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
    if ( s_.ctxt_.onpe0() )
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
    if ( s_.ctxt_.onpe0() )
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

  if ( s_.ctxt_.onpe0() )
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

