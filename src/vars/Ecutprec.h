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
// Ecutprec.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef ECUTPREC_H
#define ECUTPREC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>

class Ecutprec : public StandardVar {
  Sample *s;

  public:

  int set ( int argc, char **argv ) {
    return parse(argc, argv, s->ctrl.ecutprec, StandardVar::non_negative);
  }
  
  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->ctrl.ecutprec;
     return st.str();
  }

  Ecutprec(Sample *sample) :
    StandardVar("ecutprec", Dimensions::energy, "rydberg"),
    s(sample) {
    s->ctrl.ecutprec = 0.0;
  }

};
#endif

// Local Variables:
// mode: c++
// End:
