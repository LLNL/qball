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
// HugDeltaTemp.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef HUGDELTATEMP_H
#define HUGDELTATEMP_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class HugDeltaTemp : public Var {
  Sample *s;

  public:

  char const*name ( void ) const { return "hug_deltatemp"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> hug_deltatemp takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    s->ctrl.hug_deltatemp = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.hug_deltatemp;
     return st.str();
  }

  HugDeltaTemp(Sample *sample) : s(sample) { s->ctrl.hug_deltatemp = 10.0; }
};
#endif
