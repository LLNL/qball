////////////////////////////////////////////////////////////////////////////////
//
// Force_Complex_WF.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Force_Complex_WF.h,v 1.0 2011-03-09 15:53:20 schleife Exp $

#ifndef FORCE_COMPLEX_WF_H
#define FORCE_COMPLEX_WF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements a status variable to force complex basis [even if kpoint == (0,0,0)]
// AS: necessary for time propagation of wave functions

class Force_Complex_WF : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "force_complex_wf"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " force_complex_wf takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " force_complex_wf must be ON or OFF" << endl;
      return 1;
    }

    if ( v == "ON" )
    {

       s->wf.force_complex( true ) ;

       //ewd:  change how this is done so we don't have to set variable in a given order
       /*
         if (s->wf.nel() > 0) {
           if ( ui->oncoutpe() )
             cout << " force_complex_wf can only be set before atoms are specified!" << endl;
           return 1;
         }
         else
         { 
           s->wf.force_complex( true ) ;
         }
       */
       
    }

    if ( v == "OFF" )
    {
      s->wf.force_complex( false ) ;
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
     if ( s->wf.force_complex_set() ) st << "ON" ; else st << "OFF" ;
     return st.str();
  }

  Force_Complex_WF(Sample *sample) : s(sample) { s->wf.force_complex( false ); }
};
#endif
