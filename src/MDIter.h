////////////////////////////////////////////////////////////////////////////////
//
// MDIter.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDIter.h,v 1.5 2009/04/22 16:37:50 draeger1 Exp $

#ifndef MDITER_H
#define MDITER_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class MDIter : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "mditer"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> mditer takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v < 0 )
    {
       if ( ui->oncoutpe() )
          cout << " <ERROR> mditer must be non-negative </ERROR>" << endl;
       return 1;
    }

    s->ctrl.mditer = v;
    
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.extra_memory;
     return st.str();
  }

  MDIter(Sample *sample) : s(sample) { s->ctrl.mditer = 0; }
};
#endif
