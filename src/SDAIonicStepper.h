////////////////////////////////////////////////////////////////////////////////
//
// SDAIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDAIonicStepper.h,v 1.3 2009/03/27 00:53:24 draeger1 Exp $

#ifndef SDAIONICSTEPPER_H
#define SDAIONICSTEPPER_H

#include "IonicStepper.h"
#include "LineMinimizer.h"
#include <vector>

class SDAIonicStepper : public IonicStepper
{
  private:

  bool first_step_;
  std::vector<std::vector< double> > rc_;
  std::vector<std::vector< double> > pc_;
  std::vector<std::vector< double> > fc_;
  double ec_, fpc_;
  double alpha_, sigma1_, sigma2_;
  LineMinimizer linmin_;

  public:

  SDAIonicStepper(Sample& s) : IonicStepper(s), first_step_(true),
  sigma1_(0.1), sigma2_(0.5) { linmin_.set_sigma1(sigma1_); }

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
};

#endif
