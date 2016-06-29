////////////////////////////////////////////////////////////////////////////////
//
// Stress.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Stress.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef STRESS_H
#define STRESS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Stress : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "stress"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> stress takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> stress must be ON or OFF </ERROR>" << endl;
      return 1;
    }

    s->ctrl.stress = v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.stress;
     return st.str();
  }

  Stress(Sample *sample) : s(sample) { s->ctrl.stress = "OFF"; }
};
#endif
