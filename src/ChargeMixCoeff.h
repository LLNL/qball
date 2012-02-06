////////////////////////////////////////////////////////////////////////////////
//
// ChargeMixCoeff.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeMixCoeff.h,v 1.4 2008/09/19 21:24:50 draeger1 Exp $

#ifndef CHARGEMIXCOEFF_H
#define CHARGEMIXCOEFF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ChargeMixCoeff : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "charge_mix_coeff"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> charge_mix_coeff takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> charge_mix_coeff must be non-negative </ERROR>" << endl;
      return 1;
    }
    else if ( v > 1.0 ) {
      if ( ui->oncoutpe() ) {
        cout << " <WARNING> charge_mix_coeff should be less than or equal to one. </WARNING>" << endl;
        cout << " <WARNING> Setting charge_mix_coeff = 1.0 </WARNING>" << endl;
      }
      v = 1.0;
    }
    s->ctrl.charge_mix_coeff = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.charge_mix_coeff;
     return st.str();
  }

  ChargeMixCoeff(Sample *sample) : s(sample) { s->ctrl.charge_mix_coeff = 0.5; };
};
#endif
