////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// ExtStress.h
//
////////////////////////////////////////////////////////////////////////////////

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
