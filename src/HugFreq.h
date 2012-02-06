////////////////////////////////////////////////////////////////////////////////
//
// HugFreq.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HugFreq.h,v 1.1 2008/05/05 19:37:43 draeger1 Exp $

#ifndef HUGFREQ_H
#define HUGFREQ_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class HugFreq : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "hug_freq"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> hug_freq takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int f = atoi(argv[1]);
    if (f < 1) { 
      if ( ui->oncoutpe() )
        cout << " <ERROR> hug_freq must be positive and non-zero </ERROR>" << endl;
      return 1;
    }

    s->ctrl.hug_freq = f;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.hug_freq;
     return st.str();
  }

  HugFreq(Sample *sample) : s(sample) { s->ctrl.hug_freq = 5; }
};
#endif
