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
// RunTimer.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef RUNTIMER_H
#define RUNTIMER_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class RunTimer : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "run_timer"; };

  int set ( int argc, char **argv )
  {
    if (! (argc == 2 || argc == 3))
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> run_timer syntax error </ERROR>" << endl;
      return 1;
    }
    
    if (argc == 2) {
      int len = strlen(argv[1]);
      char c = argv[1][len-1];
      if (c == 's') {
        char tmp[len-1];
        for (int i=0; i<len-1; i++)
          tmp[i] = argv[1][i];
        s->ctrl.run_timer = atof(tmp);
      }
      else if (c == 'm') {
        char tmp[len-1];
        for (int i=0; i<len-1; i++)
          tmp[i] = argv[1][i];
        s->ctrl.run_timer = 60.*atof(tmp);
      }
      else if (c == 'h') {
        char tmp[len-1];
        for (int i=0; i<len-1; i++)
          tmp[i] = argv[1][i];
        s->ctrl.run_timer = 3600.*atof(tmp);
      }
      else if (c == 'd') {
        char tmp[len-1];
        for (int i=0; i<len-1; i++)
          tmp[i] = argv[1][i];
        s->ctrl.run_timer = 24.*3600.*atof(tmp);
      }
      else {
        s->ctrl.run_timer = atof(argv[1]);
      }
    }
    else if (argc == 3) {
      string t = argv[2];
      if (t == "d" || t == "day" || t == "days" || t == "D" || t == "Days")
        s->ctrl.run_timer = 24.*3600.*atof(argv[1]);
      else if (t == "h" || t == "hour" || t == "hours" || t == "H" || t == "Hours")
        s->ctrl.run_timer = 3600.*atof(argv[1]);
      else if (t == "m" || t == "min" || t == "mins" || t == "minutes" || t == "Min" || t == "Minutes" || t == "minute")
        s->ctrl.run_timer = 60.*atof(argv[1]);
      else
        s->ctrl.run_timer = atof(argv[1]);
    }      
    if ( ui->oncoutpe() )
      cout << " <!--RunTimer:  run_timer set to " << s->ctrl.run_timer << " seconds. -->" << endl;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.run_timer << " sec";
     return st.str();
  }

  RunTimer(Sample *sample) : s(sample) { s->ctrl.run_timer = 0.0; };
};
#endif
