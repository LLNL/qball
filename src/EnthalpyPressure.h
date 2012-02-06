////////////////////////////////////////////////////////////////////////////////
//
// EnthalpyPressure.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnthalpyPressure.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef ENTHALPYPRESSURE_H
#define ENTHALPYPRESSURE_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class EnthalpyPressure : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "enthalpy_pressure"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> enthalpy_pressure only takes one value </ERROR>" << endl;
      return 1;
    }
    const double gpa = 29421.0120;
    double v = atof(argv[1])/gpa;   // convert input value from GPa to a.u.
    s->ctrl.enthalpy_pressure = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.enthalpy_pressure;
     return st.str();
  }

  EnthalpyPressure(Sample *sample) : s(sample) { s->ctrl.enthalpy_pressure = 0.0; };
};
#endif
