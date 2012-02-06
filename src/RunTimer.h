////////////////////////////////////////////////////////////////////////////////
//
// RunTimer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunTimer.h,v 1.2 2010/08/26 17:44:16 draeger1 Exp $

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

  char *name ( void ) const { return "run_timer"; };

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
