////////////////////////////////////////////////////////////////////////////////
//
// WfDiag.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WfDiag.h,v 1.4 2009/09/29 00:06:50 draeger1 Exp $

#ifndef WFDIAG_H
#define WFDIAG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class WfDiag : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "wf_diag"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> wf_diag takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "T" || v == "F" || v == "EIGVAL" || v == "MLWF" || v == "MLWFC") )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> wf_diag must be T, F, EIGVAL, MLWF, or MLWFC </ERROR>" << endl;
      return 1;
    }

    s->ctrl.wf_diag = v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.wf_diag;
     return st.str();
  }

  WfDiag(Sample *sample) : s(sample) { s->ctrl.wf_diag = "F"; };
};
#endif
