////////////////////////////////////////////////////////////////////////////////
//
// Ecutprec.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Ecutprec.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef ECUTPREC_H
#define ECUTPREC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Ecutprec : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "ecutprec"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> ecutprec takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> ecutprec must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.ecutprec = 0.5*v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->ctrl.ecutprec;
     return st.str();
  }

  Ecutprec(Sample *sample) : s(sample) { s->ctrl.ecutprec = 0.0; }
};
#endif
