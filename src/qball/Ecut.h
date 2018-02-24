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

class Ecut : public StandardVar {
  Sample *s;

  public:

  int set ( int argc, char **argv ) {

    double value;
    
    int error = parse(argc, argv, value, StandardVar::non_negative);
    if(error != 0) return error;

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

  Ecut(Sample *sample) : StandardVar("ecut", Dimensions::energy, "rydberg"), s(sample) {};
};
#endif

// Local Variables:
// mode: c++
// End:
