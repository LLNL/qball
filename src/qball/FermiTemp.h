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
// FermiTemp.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Deprecated variable:  included to provide backward compatibility with older Qbox
// input files.  Fermi_temp inputs smearing in Kelvin, converted to Ry to be compatible
// with smearing and smearing_width variables.

#include <config.h>

#ifndef FERMITEMP_H
#define FERMITEMP_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class FermiTemp : public Var {
  Sample *s;

  public:

  char const*name ( void ) const { return "fermi_temp"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2) {
      if ( ui->oncoutpe() )
      cout << "<ERROR> fermi_temp takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);

    const double boltz = 1.0 / ( 11605.0 * 13.6058 ); // convert from Kelvin to Ry

    ui->warning(string("\nThe fermi_temp variable is deprecated. Automatically converting input\n") +
		string("to the following commands instead:\n\n") +
		string("  set smearing fermi\n") +
		string("  set smearing_width ") + argv[1] + string(" kelvin\n"));

    s->ctrl.smearing = "fermi";
    s->ctrl.smearing_width = v*boltz;
    s->ctrl.smearing_ngauss = -1;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.smearing_width;
     return st.str();
  }

  FermiTemp(Sample *sample) : s(sample) { }
};
#endif

// Local Variables:
// mode: c++
// End:
