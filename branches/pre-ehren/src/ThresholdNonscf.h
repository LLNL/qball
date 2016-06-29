////////////////////////////////////////////////////////////////////////////////
//
// ThresholdNonscf.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ThresholdNonscf.h,v 1.4 2008/04/07 22:00:37 draeger1 Exp $

#ifndef THRESHOLDNONSCF_H
#define THRESHOLDNONSCF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThresholdNonscf : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "threshold_nonscf"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> threshold_nonscf takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> threshold_nonscf must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.threshold_nonscf = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.threshold_nonscf;
     return st.str();
  }

  ThresholdNonscf(Sample *sample) : s(sample) { s->ctrl.threshold_nonscf = 0.0; };
};
#endif
