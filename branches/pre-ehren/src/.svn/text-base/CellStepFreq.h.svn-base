////////////////////////////////////////////////////////////////////////////////
//
// CellStepFreq.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellStepFreq.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef CELLSTEPFREQ_H
#define CELLSTEPFREQ_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class CellStepFreq : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "cell_stepfreq"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> cell_stepfreq takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v <= 0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> cell_stepfreq must be positive </ERROR>" << endl;
      return 1;
    }

    s->ctrl.cell_stepfreq = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_stepfreq;
     return st.str();
  }

  CellStepFreq(Sample *sample) : s(sample) { s->ctrl.cell_stepfreq = 1; }
};
#endif
