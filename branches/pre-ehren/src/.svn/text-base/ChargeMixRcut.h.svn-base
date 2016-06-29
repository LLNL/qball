////////////////////////////////////////////////////////////////////////////////
//
// ChargeMixRcut.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeMixRcut.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef CHARGEMIXRCUT_H
#define CHARGEMIXRCUT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ChargeMixRcut : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "charge_mix_rcut"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> charge_mix_rcut takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> charge_mix_rcut must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.charge_mix_rcut = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.charge_mix_rcut;
     return st.str();
  }

  ChargeMixRcut(Sample *sample) : s(sample) { s->ctrl.charge_mix_rcut = 10.0; };
};
#endif
