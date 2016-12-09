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
// Dt.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef DT_H
#define DT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Dt : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "dt"; };

  int set ( int argc, char **argv ) {
    string unit_name;

    if ( argc == 3 ){
      unit_name = argv[2];
    } else if ( argc == 2 ) {
      unit_name = "atomictime";
      ui->warning("Missing units for the 'dt' variable. Assuming 'atomictime'.");
    } else {
      ui->error("The variable 'dt' requires two arguments: the value and the units.");
      return 1;
    }

    cout << "NAME " << unit_name << " - " << argv[2] << endl;
    
    Unit unit(Dimensions::time, unit_name);

    if(!unit.exists()) {
      ui->error("Unknown time unit '" + unit_name + "'.");
      return 1;
    }
    
    double value = unit.to_atomic(atof(argv[1]));
    
    if ( value < 0.0 ) {
      ui->error("The variable 'dt' must be non-negative");
      return 1;
    }

    s->ctrl.dt = value;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.dt;
     return st.str();
  }

  Dt(Sample *sample) : s(sample) { s->ctrl.dt = 3.0; }
};
#endif
