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
// ThresholdScf.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef THRESHOLDSCF_H
#define THRESHOLDSCF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThresholdScf : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "threshold_scf"; };

  int set ( int argc, char **argv ) {
    if ( argc < 2 || argc > 3 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> threshold_scf takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> threshold_scf must be non-negative </ERROR>" << endl;
      return 1;
    }
    s->ctrl.threshold_scf = v;
    if (argc == 3) {
      int n = atoi(argv[2]);
      if ( n < 2 ) {
        if ( ui->oncoutpe() )
          cout << " <ERROR> nsteps must be greater than 1 </ERROR>" << endl;
        return 1;
      }
      s->ctrl.threshold_scf_nsteps = n;
    }
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.threshold_scf << "  " << s->ctrl.threshold_scf_nsteps;
     return st.str();
  }

  ThresholdScf(Sample *sample) : s(sample) { 
    s->ctrl.threshold_scf = 0.0; 
    s->ctrl.threshold_scf_nsteps = 2; 
  };
};
#endif
