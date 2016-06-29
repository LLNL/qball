////////////////////////////////////////////////////////////////////////////////
//
// Nempty.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Nempty.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef NEMPTY_H
#define NEMPTY_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Nempty : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "nempty"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> nempty takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v < 0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> nempty must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->wf.set_nempty(v);
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.nempty();
     return st.str();
  }

  Nempty(Sample *sample) : s(sample) {};
};
#endif
