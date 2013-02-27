////////////////////////////////////////////////////////////////////////////////
//
// BMDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BMDIonicStepper.h,v 1.1 2010/01/16 01:26:35 draeger1 Exp $

//
// IonicStepper is used in the following way
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

#ifndef BMDIONICSTEPPER_H
#define BMDIONICSTEPPER_H

#include "IonicStepper.h"

class BMDIonicStepper : public IonicStepper
{
  private:

  double e0_,em_;
  double ekin_;
  std::vector<std::vector< double> > fm_;
  void compute_ekin(void);

  public:

  BMDIonicStepper(Sample& s) : IonicStepper(s)
  {
    e0_ = 0.0;
    em_ = 0.0;
    ekin_ = 0.0;
    atoms_.get_positions(r0_);
    atoms_.get_velocities(v0_);
    fm_.resize(r0_.size());
    for ( int is = 0; is < fm_.size(); is++ )
      fm_[is].resize(r0_[is].size());
    compute_ekin();
  }

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0);
  double ekin(void) const { return ekin_; }
  double temp(void) const
  {
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    if ( ndofs_ > 0 )
      return 2.0 * ( ekin_ / boltz ) / ndofs_;
    else
      return 0.0;
  }
};

#endif
