////////////////////////////////////////////////////////////////////////////////
//
// CellDyn.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellDyn.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef CELLDYN_H
#define CELLDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class CellDyn : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "cell_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> cell_dyn takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "LOCKED" || v == "SD" ) )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> cell_dyn must be in [LOCKED,SD] </ERROR>" << endl;
      return 1;
    }

    s->ctrl.cell_dyn = v;
        
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_dyn;
     return st.str();
  }

  CellDyn(Sample *sample) : s(sample) { s->ctrl.cell_dyn = "LOCKED"; };
};
#endif
