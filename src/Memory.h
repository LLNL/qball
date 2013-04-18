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
// Memory.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MEMORY_H
#define MEMORY_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Memory : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "memory"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> memory takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "huge" || v == "HUGE" || v == "large" || v == "LARGE" || v == "normal" || v == "NORMAL" || v == "compact" || v == "COMPACT") ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> memory must be huge, large or compact </ERROR>" << endl;
      return 1;
    }
    if (v == "HUGE" || v == "huge") s->ctrl.extra_memory = 9;
    if (v == "LARGE" || v == "large") s->ctrl.extra_memory = 5;
    if (v == "NORMAL" || v == "normal") s->ctrl.extra_memory = 3;
    if (v == "COMPACT" || v == "compact") s->ctrl.extra_memory = 0;

    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.extra_memory;
     return st.str();
  }

  Memory(Sample *sample) : s(sample) { s->ctrl.extra_memory = 3; }
};
#endif
