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
// FcpThWidth.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FcpThWidth.h,v 1.4 2008-09-08 15:56:19 fgygi Exp $

#include <config.h>

#ifndef FCPTHWIDTH_H
#define FCPTHWIDTH_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>

class FcpThWidth : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "fcp_th_width"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " fcp_th_width takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v <= 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " fcp_th_width must be positive" << endl;
      return 1;
    }

    s->ctrl.fcp_th_width = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.fcp_th_width;
     return st.str();
  }

  // Use same value of th_width if not specified (implied as negative value).
  FcpThWidth(Sample *sample) : s(sample) { s->ctrl.fcp_th_width = -1.0; }
};
#endif

// Local Variables:
// mode: c++
// End:
