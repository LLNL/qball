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
// Xc.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef XC_H
#define XC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Xc : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "xc"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> xc takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "LDA" || v == "PBE" || v == "PBEsol" || v == "PBErev" || v == "BLYP" ) )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> xc must be LDA, PBE, PBEsol, PBErev or BLYP </ERROR>" << endl;
      return 1;
    }

    s->ctrl.xc= v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.xc;
     return st.str();
  }

  Xc(Sample *sample) : s(sample) { s->ctrl.xc = "LDA"; };
};
#endif
