////////////////////////////////////////////////////////////////////////////////
//
// ThresholdForce.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ThresholdForce.h,v 1.4 2007/05/29 19:44:13 draeger1 Exp $

#ifndef THRESHOLDFORCE_H
#define THRESHOLDFORCE_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThresholdForce : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "threshold_force"; };

  int set ( int argc, char **argv ) {
    if ( argc < 2 || argc > 3 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> threshold_force takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> threshold_force must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.threshold_force = v;
    if (argc == 3) {
      int n = atoi(argv[2]);
      if ( n < 2 ) {
        if ( ui->oncoutpe() )
          cout << " <ERROR> nsteps must be greater than 1 </ERROR>" << endl;
        return 1;
      }
      s->ctrl.threshold_force_nsteps = n;
    }
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.threshold_force << "  " << s->ctrl.threshold_force_nsteps;
     return st.str();
  }

  ThresholdForce(Sample *sample) : s(sample) { 
    s->ctrl.threshold_force = 0.0; 
    s->ctrl.threshold_force_nsteps = 2; 
  };
};
#endif
