////////////////////////////////////////////////////////////////////////////////
//
// Spin.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Spin.h,v 1.1 2009/09/25 23:18:11 draeger1 Exp $

#ifndef SPIN_H
#define SPIN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Spin : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "spin"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> spin takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> spin must be non-negative </ERROR>" << endl;
      return 1;
    }

    int deltaspin = (int)v;
    s->wf.set_deltaspin(deltaspin);
    s->ctrl.delta_spin = deltaspin;

    int nup = s->wf.nst(0);
    int ndown = s->wf.nst(1);
    double testdelta = (double)(nup-ndown)/2.0;
    
    if ( ui->oncoutpe() ) {
      cout << "<nel_up> " << nup << " </nel_up>" << endl;
      cout << "<nel_down> " << ndown << " </nel_down>" << endl;
      cout << "<total_spin> " << testdelta << " </total_spin>" << endl;
    }
    if (testdelta != v) {
      if ( ui->oncoutpe() ) {
        int nel = s->wf.nel();
        cout << " <ERROR> spin does not match nel = " << nel << "! </ERROR>" << endl;
      }
      return 1;
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.deltaspin();
     return st.str();
  }

  Spin(Sample *sample) : s(sample) { s->ctrl.delta_spin = 0; };
};
#endif
