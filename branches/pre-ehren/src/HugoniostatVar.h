////////////////////////////////////////////////////////////////////////////////
//
// Hugoniostat.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HugoniostatVar.h,v 1.1 2008/05/05 19:37:43 draeger1 Exp $

#ifndef HUGONIOSTAT_H
#define HUGONIOSTAT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class HugoniostatVar : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "hugoniostat"; };

  int set ( int argc, char **argv ) {
    if ( argc != 4 ) {
      if ( ui->oncoutpe() )
      cout << "<ERROR> hugoniostat requires three input values:  etot, volume, pressure </ERROR>" << endl;
      return 1;
    }

    double eref = atof(argv[1]);
    double vref = atof(argv[2]);
    double pref = atof(argv[3]);
    s->ctrl.hugoniostat = "ON";
    s->ctrl.hug_etot = eref;
    s->ctrl.hug_volume = vref;
    s->ctrl.hug_pressure = pref;

    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.hugoniostat;
     return st.str();
  }

  HugoniostatVar(Sample *sample) : s(sample) { s->ctrl.hugoniostat = "OFF"; };
};
#endif
