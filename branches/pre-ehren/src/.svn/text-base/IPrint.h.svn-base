////////////////////////////////////////////////////////////////////////////////
//
// IPrint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: IPrint.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef IPRINT_H
#define IPRINT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// Iprint variable is used to control the output of eigenvalues and occupation
// numbers (which can be substantial for many states and/or k-points)

class IPrint : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "iprint"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> iprint takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v < 0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> iprint must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.iprint = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.iprint;
     return st.str();
  }

  IPrint(Sample *sample) : s(sample) { s->ctrl.iprint = 1; }
};
#endif
