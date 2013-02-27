////////////////////////////////////////////////////////////////////////////////
//
// WF_Phase_Real.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WF_Phase_Real.h,v 1.0 2011-03-09 15:53:20 schleife Exp $

#ifndef WF_PHASE_REAL_H
#define WF_PHASE_REAL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements a status variable to force complex basis [even if kpoint == (0,0,0)]
// AS: necessary for time propagation of wave functions

class WF_Phase_RealVar : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "phase_real_wf"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " phase_real_wf takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " phase_real_wf must be ON or OFF" << endl;
      return 1;
    }

    if ( v == "ON" )
    {
      s->wf.phase_real( true ) ;
    }

    if ( v == "OFF" )
    {
      s->wf.phase_real( false ) ;
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) ;
     if ( s->wf.phase_real_set() ) st << "ON" ; else st << "OFF" ;
     return st.str();
  }

  WF_Phase_RealVar(Sample *sample) : s(sample) { s->wf.phase_real( false ); }
};
#endif
