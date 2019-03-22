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
// LaserEnvelope.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef LASERENVELOPE_H
#define LASERENVELOPE_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include <qball/Sample.h>

class LaserEnvelope : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "laser_envelope"; };
  
  int set ( int argc, char **argv ) {
    if ( argc == 2 ) {
      s->ctrl.envelope_type = argv[1];

      if ( s->ctrl.envelope_type != "constant" ){ 
        ui->error(argv[1]);
        ui->error("laser_envelope must be [type] [center of envelope] [width of envelope]. type can be: constant or gaussian. ");
        return 1;
      }
    }

    else if ( argc == 4 ) { 
    s->ctrl.envelope_type = argv[1];
    s->ctrl.envelope_center = atof(argv[2]);
    s->ctrl.envelope_width = atof(argv[3]);
    }

    else {
      ui->error("laser_envelope must be [type] [center of envelope] [width of envelope]. type can be: constant or gaussian. ");
      return 1;
    }

    return 0;
  }
  
  string print (void) const
  {
    ostringstream st;
    st.setf(ios::left,ios::adjustfield);
    st << setw(10) << name() << " = ";
    st.setf(ios::right,ios::adjustfield);
    st << setw(10) << s->ctrl.envelope_type;
    st << " " << setw(10) <<  s->ctrl.envelope_center << " " << setw(10) << s->ctrl.envelope_width;
    return st.str();
  }
  
  LaserEnvelope(Sample *sample) : s(sample) { s->ctrl.envelope_type = "constant"; s->ctrl.envelope_center = 0.0 ; s->ctrl.envelope_width = 0.0;}

};

#endif

// Local Variables:
// mode: c++
// End:
