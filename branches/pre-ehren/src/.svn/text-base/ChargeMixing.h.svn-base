////////////////////////////////////////////////////////////////////////////////
//
// ChargeMixing.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeMixing.h,v 1.5 2009/04/22 16:37:50 draeger1 Exp $

#ifndef CHARGEMIXING_H
#define CHARGEMIXING_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ChargeMixing : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "charge_mixing"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> charge_mixing takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "simple" || v == "SIMPLE" || v == "Simple" || v == "anderson" || v == "ANDERSON" || v == "Anderson" || v == "off" || v == "OFF")  ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> charge_mixing must be simple, anderson, or off </ERROR>" << endl;
      return 1;
    }
    if (v == "SIMPLE" || v == "Simple") v = "simple";
    if (v == "ANDERSON" || v == "Anderson") v = "anderson";
    if (v == "OFF") v = "off";

    s->ctrl.charge_mixing = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.charge_mixing;
     return st.str();
  }

  ChargeMixing(Sample *sample) : s(sample) { s->ctrl.charge_mixing = "anderson"; }
};
#endif
