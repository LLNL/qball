////////////////////////////////////////////////////////////////////////////////
//
// ExtStress.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ExtStress.h,v 1.4 2008/04/07 22:00:37 draeger1 Exp $

#ifndef EXTSTRESS_H
#define EXTSTRESS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ExtStress : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "ext_stress"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 7 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> ext_stress must be specified as s_xx,s_yy,s_zz,s_xy,s_yz,s_xz" 
           << endl;
      return 1;
    }
    
    for ( int i = 0; i < 6; i++ )
      s->ctrl.ext_stress[i] = atof(argv[i+1]);
      
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     for ( int i = 0; i < 6; i++ )
       st << setw(10) << s->ctrl.ext_stress[i];
     return st.str();
  }

  ExtStress(Sample *sample) : s(sample) {
    s->ctrl.ext_stress[0]=0.0;
    s->ctrl.ext_stress[1]=0.0;
    s->ctrl.ext_stress[2]=0.0;
    s->ctrl.ext_stress[3]=0.0;
    s->ctrl.ext_stress[4]=0.0;
    s->ctrl.ext_stress[5]=0.0;
  };
};
#endif
