////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// Debug.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

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
