////////////////////////////////////////////////////////////////////////////////
//
// ThresholdStress.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ThresholdStress.h,v 1.4 2007/05/29 19:44:13 draeger1 Exp $

#ifndef THRESHOLDSTRESS_H
#define THRESHOLDSTRESS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThresholdStress : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "threshold_stress"; };

  int set ( int argc, char **argv ) {
    if ( argc < 2 || argc > 3 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> threshold_stress takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> threshold_stress must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.threshold_stress = v;
    if (argc == 3) {
      int n = atoi(argv[2]);
      if ( n < 2 ) {
        if ( ui->oncoutpe() )
          cout << " <ERROR> nsteps must be greater than 1 </ERROR>" << endl;
        return 1;
      }
      s->ctrl.threshold_stress_nsteps = n;
    }
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.threshold_stress << "  " << s->ctrl.threshold_stress_nsteps;
     return st.str();
  }

  ThresholdStress(Sample *sample) : s(sample) { 
    s->ctrl.threshold_stress = 0.0; 
    s->ctrl.threshold_stress_nsteps = 2; 
  };
};
#endif
