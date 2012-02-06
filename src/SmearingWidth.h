////////////////////////////////////////////////////////////////////////////////
//
// SmearingWidth.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SmearingWidth.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef SMEARINGWIDTH_H
#define SMEARINGWIDTH_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class SmearingWidth : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "smearing_width"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> smearing_width takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> smearing_width must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.smearing_width = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.smearing_width;
     return st.str();
  }

  SmearingWidth(Sample *sample) : s(sample) { s->ctrl.smearing_width = 0.0; }
};
#endif
