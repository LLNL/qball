////////////////////////////////////////////////////////////////////////////////
//
// Ecutden.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Ecutden.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef ECUTDEN_H
#define ECUTDEN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Ecutden : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "ecutden"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> ecutden takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> ecutden must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.ecutden = 0.5*v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->ctrl.ecutden;
     return st.str();
  }

  Ecutden(Sample *sample) : s(sample) { s->ctrl.ecutden = 0.0; }
};
#endif
