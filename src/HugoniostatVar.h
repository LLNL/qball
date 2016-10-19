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
// Hugoniostat.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef HUGONIOSTAT_H
#define HUGONIOSTAT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class HugoniostatVar : public Var {
  Sample *s;

  public:

  char const*name ( void ) const { return "hugoniostat"; };

  int set ( int argc, char **argv ) {
    if ( argc != 4 ) {
      if ( ui->oncoutpe() )
      cout << "<ERROR> hugoniostat requires three input values:  etot, volume, pressure </ERROR>" << endl;
      return 1;
    }

    double eref = atof(argv[1]);
    double vref = atof(argv[2]);
    double pref = atof(argv[3]);
    s->ctrl.hugoniostat = "ON";
    s->ctrl.hug_etot = eref;
    s->ctrl.hug_volume = vref;
    s->ctrl.hug_pressure = pref;

    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.hugoniostat;
     return st.str();
  }

  HugoniostatVar(Sample *sample) : s(sample) { s->ctrl.hugoniostat = "OFF"; };
};
#endif
