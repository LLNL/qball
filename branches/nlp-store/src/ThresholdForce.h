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
// ThresholdForce.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef THRESHOLDFORCE_H
#define THRESHOLDFORCE_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThresholdForce : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "threshold_force"; };

  int set ( int argc, char **argv ) {
    if ( argc < 2 || argc > 3 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> threshold_force takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> threshold_force must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.threshold_force = v;
    if (argc == 3) {
      int n = atoi(argv[2]);
      if ( n < 2 ) {
        if ( ui->oncoutpe() )
          cout << " <ERROR> nsteps must be greater than 1 </ERROR>" << endl;
        return 1;
      }
      s->ctrl.threshold_force_nsteps = n;
    }
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.threshold_force << "  " << s->ctrl.threshold_force_nsteps;
     return st.str();
  }

  ThresholdForce(Sample *sample) : s(sample) { 
    s->ctrl.threshold_force = 0.0; 
    s->ctrl.threshold_force_nsteps = 2; 
  };
};
#endif
