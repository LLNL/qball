////////////////////////////////////////////////////////////////////////////////
//
// Thermostat.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Thermostat.h,v 1.5 2008/04/07 22:00:37 draeger1 Exp $

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Thermostat : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "thermostat"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << "<ERROR> thermostat takes only one value </ERROR>" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "SCALING" || v == "ANDERSEN" || v == "LOWE" || v == "OFF" ) ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> Thermostat must be SCALING or ANDERSEN or LOWE or OFF </ERROR>"
             << endl;
      return 1;
    }

    s->ctrl.thermostat = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.thermostat;
     return st.str();
  }

  Thermostat(Sample *sample) : s(sample) { s->ctrl.thermostat = "OFF"; };
};
#endif
