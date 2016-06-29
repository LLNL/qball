////////////////////////////////////////////////////////////////////////////////
//
// Debug.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Debug.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef DEBUG_H
#define DEBUG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Debug : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "debug"; };

  int set ( int argc, char **argv )
  {
    const string allowed_values("OFF STRESS");
    if ( argc < 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> debug must be in " << allowed_values << endl;
      return 1;
    }
    
    string v;
    for ( int iarg = 1; iarg < argc; iarg++ )
    {
      string vt = argv[iarg];
      if ( allowed_values.find(vt) == string::npos )
      {
        // not a valid argument
        if ( ui->oncoutpe() )
        {
          cout << vt << " is not a valid debug option  </ERROR>" << endl;
          cout << " <ERROR> valid debug options are: " << allowed_values << endl;
        }
        return 1;
      }
        
      v += " " + vt;
    }
    
    // check if OFF is found in the string v
    if ( v.find("OFF") != string::npos )
      v = "OFF";

    s->ctrl.debug = v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.debug;
     return st.str();
  }

  Debug(Sample *sample) : s(sample) { s->ctrl.debug = "OFF"; }
};
#endif
