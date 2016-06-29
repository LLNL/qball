////////////////////////////////////////////////////////////////////////////////
//
// HugDeltaTemp.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HugDeltaTemp.h,v 1.1 2008/05/05 19:37:43 draeger1 Exp $

#ifndef HUGDELTATEMP_H
#define HUGDELTATEMP_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class HugDeltaTemp : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "hug_deltatemp"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> hug_deltatemp takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    s->ctrl.hug_deltatemp = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.hug_deltatemp;
     return st.str();
  }

  HugDeltaTemp(Sample *sample) : s(sample) { s->ctrl.hug_deltatemp = 10.0; }
};
#endif
