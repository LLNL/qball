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
// Smearing.h
//
////////////////////////////////////////////////////////////////////////////////
//
// set method for smearing occupation numbers, currently either Fermi-Dirac or Gaussian
// (Methfessel and Paxton, PRB 40, 3616 (1989)).  In the latter case, the number of terms 
// in the expansion can be given as an optional second argument.  

#include <config.h>

#ifndef SMEARING_H
#define SMEARING_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Smearing : public Var {
  Sample *s;

  public:

  char const*name ( void ) const { return "smearing"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 && argc != 3) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> smearing takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "fermi" || v == "FERMI" || v == "Fermi" || v == "gaussian" || v == "GAUSSIAN" || v == "Gaussian")  ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> smearing must be fermi or gaussian </ERROR>" << endl;
      return 1;
    }
    if (v == "FERMI" || v == "Fermi") v = "fermi";
    if (v == "GAUSSIAN" || v == "Gaussian") v = "gaussian";

    s->ctrl.smearing_ngauss = 0;
    if (argc == 3 && v == "gaussian") {
      int n = atoi(argv[2]);
      s->ctrl.smearing_ngauss = n;
    }
    if (v == "fermi")
      s->ctrl.smearing_ngauss = -1;

    s->ctrl.smearing = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.smearing;
     return st.str();
  }

  Smearing(Sample *sample) : s(sample) { 
    s->ctrl.smearing = "fermi"; 
    s->ctrl.smearing_ngauss = -1;    // default: Fermi smearing
  }
};
#endif
