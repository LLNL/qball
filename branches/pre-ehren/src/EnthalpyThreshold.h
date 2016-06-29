////////////////////////////////////////////////////////////////////////////////
//
// EnthalpyThreshold.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnthalpyThreshold.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef ENTHALPYTHRESHOLD_H
#define ENTHALPYTHRESHOLD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class EnthalpyThreshold : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "enthalpy_threshold"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> enthalpy_threshold only takes one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    s->ctrl.enthalpy_threshold = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.enthalpy_threshold;
     return st.str();
  }

  EnthalpyThreshold(Sample *sample) : s(sample) { s->ctrl.enthalpy_threshold = 1.E-4; };
};
#endif
