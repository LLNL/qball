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
// CenterOfMassc.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef CENTEROFMASS_H
#define CENTEROFMASS_H

#include <string>
#include <iostream>
#include <cstdlib>

#include "Sample.h"

// CenterOfMass variable determines whether or not to subtract the
// center of mass velocity from particle velocities.  The default
// value is "unset", allowing us to turn it on as a default with
// atoms_dyn = MD while still giving the user the chance to turn it
// off by setting this variable to free before setting atoms_dyn MD.

class CenterOfMass : public Var {
  Sample *s;

  public:

  char const*name ( void ) const { return "center_of_mass"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> center_of_mass takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];

    if (v == "FIXED" || v=="fixed") { 
      v = "fixed";
    }
    else if (v == "FREE" || v=="free") { 
      v = "free";
    }
    else {
      if ( ui->oncoutpe() )
        cout << " <ERROR> center_of_mass must be set to either free or fixed </ERROR>" << endl;
      return 1;
    }

    s->ctrl.center_of_mass = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.center_of_mass;
     return st.str();
  }

  CenterOfMass(Sample *sample) : s(sample) { s->ctrl.center_of_mass = "unset"; }
};
#endif
