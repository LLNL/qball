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
// TDDt.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDDt.h,v 1.5 2011-04-04 15:56:18 schleife Exp $

#include <config.h>

#ifndef TDDT_H
#define TDDT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>
#include "StandardVar.h"

// AS: this class implements tddt as status variable
// AS: which specifies the dt for the wave function time propagation

class TDDt : public StandardVar {
  Sample *s;

  public:

  int set ( int argc, char **argv ) {
    return parse(argc, argv, s->ctrl.tddt, StandardVar::non_zero);
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.tddt;
     return st.str();
  }

  TDDt(Sample *sample):
    StandardVar("TD_dt", Dimensions::time, "atomictime"),
    s(sample) {
    s->ctrl.tddt = 0.01;
  }

};

#endif

// Local Variables:
// mode: c++
// End:
