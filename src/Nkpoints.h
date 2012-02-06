////////////////////////////////////////////////////////////////////////////////
//
// Nkpoints.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Nkpoints.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef NKPOINTS_H
#define NKPOINTS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// Nkpoints variable is used solely to avoid Wavefunction deallocate and allocate calls 
// after each Kpoint is defined, which can be quite expensive on large platforms like BG/L

class Nkpoints : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "nkpoints"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> nkpoints takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v <= 0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> nkpoints must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.nkpoints = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.nkpoints;
     return st.str();
  }

  Nkpoints(Sample *sample) : s(sample) { s->ctrl.nkpoints = 1; }
};
#endif
