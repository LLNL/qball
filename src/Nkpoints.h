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
// Nkpoints.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef NKPOINTS_H
#define NKPOINTS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// Nkpoints variable is used solely to avoid Wavefunction deallocate and allocate calls 
// after each Kpoint is defined, which can be quite expensive on large platforms like BG/L

class Nkpoints : public Var {
  Sample *s;

  public:

  char const*name ( void ) const { return "nkpoints"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> nkpoints takes only one value </ERROR>" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v <= 0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> nkpoints must be non-negative </ERROR>" << endl;
      return 1;
    }

    s->ctrl.nkpoints = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.nkpoints;
     return st.str();
  }

  Nkpoints(Sample *sample) : s(sample) { s->ctrl.nkpoints = 1; }
};
#endif
