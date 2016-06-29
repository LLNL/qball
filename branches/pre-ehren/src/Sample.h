////////////////////////////////////////////////////////////////////////////////
//
// Sample.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Sample.h,v 1.6 2009/09/25 23:18:11 draeger1 Exp $

#ifndef SAMPLE_H
#define SAMPLE_H

#include "AtomSet.h"
#include "Wavefunction.h"
#include "Control.h"
#include "ConstraintSet.h"
#include "SymmetrySet.h"
#include <vector>
#include <complex>

class Context;

class Sample {
  private:
  
  public:
  
  const Context& ctxt_;

  AtomSet atoms;
  ConstraintSet constraints;
  Wavefunction wf;
  Wavefunction* wfv; // wavefunction velocity
  Control ctrl;
  SymmetrySet symmetries;
  vector<vector<complex<double> > > rhog_last; // previous charge density (to avoid discontinuity in restart)

  Sample(const Context& ctxt) : ctxt_(ctxt), atoms(ctxt), wf(ctxt), wfv(0), symmetries(ctxt),
      constraints(ctxt) { ctrl.sigmas = 0.5; ctrl.facs = 2.0; }
  ~Sample(void) { delete wfv; }

};

#endif
