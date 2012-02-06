////////////////////////////////////////////////////////////////////////////////
//
// Memory.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Memory.h,v 1.5 2009/04/22 16:37:50 draeger1 Exp $

#ifndef MEMORY_H
#define MEMORY_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Memory : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "memory"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> memory takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "huge" || v == "HUGE" || v == "large" || v == "LARGE" || v == "normal" || v == "NORMAL" || v == "compact" || v == "COMPACT") ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> memory must be huge, large or compact </ERROR>" << endl;
      return 1;
    }
    if (v == "HUGE" || v == "huge") s->ctrl.extra_memory = 9;
    if (v == "LARGE" || v == "large") s->ctrl.extra_memory = 5;
    if (v == "NORMAL" || v == "normal") s->ctrl.extra_memory = 3;
    if (v == "COMPACT" || v == "compact") s->ctrl.extra_memory = 0;

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

  Memory(Sample *sample) : s(sample) { s->ctrl.extra_memory = 3; }
};
#endif
