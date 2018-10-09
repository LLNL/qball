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
// SmearingWidth.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SMEARINGWIDTH_H
#define SMEARINGWIDTH_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>
#include "StandardVar.h"

class SmearingWidth : public StandardVar {
  Sample *s;

  public:

  int set ( int argc, char **argv ) {
    return parse(argc, argv,  s->ctrl.smearing_width, StandardVar::non_negative);
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.smearing_width;
     return st.str();
  }

  SmearingWidth(Sample *sample):
    StandardVar("smearing_width", Dimensions::energy, "hartree"),
    s(sample) {
    s->ctrl.smearing_width = 0.0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
