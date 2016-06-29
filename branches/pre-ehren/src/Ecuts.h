////////////////////////////////////////////////////////////////////////////////
//
// Ecuts.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Ecuts.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef ECUTS_H
#define ECUTS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Ecuts : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "ecuts"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> ecuts takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> ecuts must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.ecuts = 0.5 * v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->ctrl.ecuts;
     return st.str();
  }

  Ecuts(Sample *sample) : s(sample) { s->ctrl.ecuts = 0.0; }
};
#endif
