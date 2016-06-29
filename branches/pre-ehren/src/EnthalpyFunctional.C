////////////////////////////////////////////////////////////////////////////////
//
// EnthalpyFunctional.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnthalpyFunctional.C,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#include "EnthalpyFunctional.h"
#include "Basis.h"
#include "FourierTransform.h"
#include <cassert>
#include <cmath>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
EnthalpyFunctional::EnthalpyFunctional(const ChargeDensity& cd, double pressure, double threshold):
cd_(cd), ctxt_(cd.vcontext()), vft_(*cd_.vft()), vbasis_(*cd_.vbasis()),
pressure_(pressure), threshold_(threshold) {
  np012loc_ = vft_.np012loc();
  sigma_ = threshold_/3.0;

  // smeared step function between rho=threshold-nsigma*sigma and rho=threshold+nsigma*sigma
  nsigma_ = 6;
  const int npts = 300;
  stepfn_.resize(npts);
  stepdx_ = 2.*nsigma_*sigma_/npts;
  for (int i=0; i<npts; i++) {
    //double x = i*stepdx_;
    //stepfn_[i] = erf(x/(sigma_*sqrt(2.)));

    // form used in PWSCF
    double x = i*stepdx_ - nsigma_*sigma_;
    stepfn_[i] = 0.5*(1.0+erf(x/(sigma_*sqrt(2.))));
  }

}

////////////////////////////////////////////////////////////////////////////////
EnthalpyFunctional::~EnthalpyFunctional(void) {

}

////////////////////////////////////////////////////////////////////////////////
void EnthalpyFunctional::update(vector<vector<double> >& vr) {

  const int nspin_ = cd_.rhor.size();
  const int size = vr[0].size();

  evol_ = 0.0;

  const int npts = stepfn_.size();

  const double fac1 = pressure_/(sigma_*sqrt(2.*3.141592653589793));
  const double fac2 = 0.5/(sigma_*sigma_);

  if (nspin_ == 1) { 
    for ( int i = 0; i < size; i++ ) {
      double tmprho = cd_.rhor[0][i];
      if (tmprho > threshold_-nsigma_*sigma_) {
        //int rhoind = (int)(tmprho/stepdx_);
        int rhoind = (int)((tmprho - threshold_ + nsigma_*sigma_)/stepdx_);
        if (rhoind < npts-1) { 
          //evol_ += stepfn_[rhoind];
          // use linear interpolation between data points:
          //   stepfn[i] + ((stepfn[i+1]-stepfn[i])/dx)*[(tmprho-(thresh-nsigma*sigma))-rho[at pt i])]
          //evol_ += stepfn_[rhoind];
          evol_ += stepfn_[rhoind] + (stepfn_[rhoind+1]-stepfn_[rhoind])*(tmprho-threshold_+nsigma_*sigma_-(double)rhoind*stepdx_)/stepdx_;
        }
        else {
          evol_ += 1.0;
        }
      }
      double delrho = tmprho-threshold_;
      vr[0][i] += fac1*exp(-fac2*delrho*delrho);
    }
  }
  else {
    cout << "<ERROR> EnthalpyFunctional does not support nspin > 1! </ERROR>" << endl;
    assert(false);
  }

  // sum volume contributions over rows
  double tsum = evol_ * vbasis_.cell().volume() / vft_.np012();
  ctxt_.dsum(1,1,&tsum,1);
  evol_ = tsum;

  epv_ = evol_*pressure_;
}
