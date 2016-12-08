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
// Ecut.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef ECUT_H
#define ECUT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "Unit.h"

class Ecut : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "ecut"; };

  int set ( int argc, char **argv ) {

    string unit_name;
    
    if ( argc == 2 ) {
      unit_name = "rydberg";
      if ( ui->oncoutpe() ) cout << " <WARNING> units missing for variable ecut, assuming rydberg </WARNING>" << endl; 
    } else if ( argc != 3 ) {
      if ( ui->oncoutpe() ) cout << " <ERROR> ecut takes only one value followed by its units </ERROR>" << endl; 
      return 1;
    } else {
      unit_name = argv[2];
    }
    
    Unit unit = Unit::Energy(unit_name);

    if(!unit.exists()) {
      if ( ui->oncoutpe() ) cout << " <ERROR> unknown energy unit '" << unit_name << "' </ERROR>" << endl; 
      return 1;
    }
    
    double value = unit.to_atomic(atof(argv[1]));

    if ( value < 0.0 ) {
      if ( ui->oncoutpe() ) cout << " <ERROR> ecut must be non-negative </ERROR>" << endl;
      return 1;
    }
    
    if ( s->wf.ecut() == value ) return 0;

    if (s->wf.hasdata()) {
      s->wf.resize(value);
      if ( s->wfv != 0 ) {
        s->wfv->resize(value);
        s->wfv->clear();
      }
    } else {
      s->wf.set_ecut(value);
      if ( s->wfv != 0 ) s->wfv->set_ecut(value);
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->wf.ecut();
     return st.str();
  }

  Ecut(Sample *sample) : s(sample) {};
};
#endif
