////////////////////////////////////////////////////////////////////////////////
//
// ReshapeContext.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ReshapeContext.h,v 1.4 2008/04/07 22:00:37 draeger1 Exp $

#ifndef RESHAPECONTEXT_H
#define RESHAPECONTEXT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ReshapeContext : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "reshape_context"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> reshape_context takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( v == "ON" || v == "on" || v == "yes" || v == "YES" || v == "T" || v == "true") {
      s->ctrl.reshape_context = true;
    }
    else {
      s->ctrl.reshape_context = false;
    }
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.reshape_context;
     return st.str();
  }

  ReshapeContext(Sample *sample) : s(sample) { s->ctrl.reshape_context = false; }
};
#endif
