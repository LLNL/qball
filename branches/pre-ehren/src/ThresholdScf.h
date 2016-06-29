////////////////////////////////////////////////////////////////////////////////
//
// ThresholdScf.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ThresholdScf.h,v 1.6 2008/04/07 22:00:37 draeger1 Exp $

#ifndef THRESHOLDSCF_H
#define THRESHOLDSCF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThresholdScf : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "threshold_scf"; };

  int set ( int argc, char **argv ) {
    if ( argc < 2 || argc > 3 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> threshold_scf takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> threshold_scf must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.threshold_scf = v;
    if (argc == 3) {
      int n = atoi(argv[2]);
      if ( n < 2 ) {
        if ( ui->oncoutpe() )
          cout << " <ERROR> nsteps must be greater than 1 </ERROR>" << endl;
        return 1;
      }
      s->ctrl.threshold_scf_nsteps = n;
    }
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.threshold_scf << "  " << s->ctrl.threshold_scf_nsteps;
     return st.str();
  }

  ThresholdScf(Sample *sample) : s(sample) { 
    s->ctrl.threshold_scf = 0.0; 
    s->ctrl.threshold_scf_nsteps = 2; 
  };
};
#endif
