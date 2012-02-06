////////////////////////////////////////////////////////////////////////////////
//
// CellMass.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellMass.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef CELLMASS_H
#define CELLMASS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class CellMass : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "cell_mass"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> cell_mass takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v <= 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> cell_mass must be positive </ERROR>" << endl;
      return 1;
    }

    s->ctrl.cell_mass = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_mass;
     return st.str();
  }

  CellMass(Sample *sample) : s(sample) { s->ctrl.cell_mass = 10000.0; }
};
#endif
