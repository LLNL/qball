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
// Spin.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SPIN_H
#define SPIN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Spin : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "spin"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> spin takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> spin must be non-negative </ERROR>" << endl;
      return 1;
    }

    int deltaspin = (int)v;
    s->wf.set_deltaspin(deltaspin);
    s->ctrl.delta_spin = deltaspin;

    int nup = s->wf.nst(0) - s->wf.nempty();
    int ndown = s->wf.nst(1) - s->wf.nempty();
    double testdelta = (double)(nup-ndown)/2.0;
    
    if ( ui->oncoutpe() ) {
      cout << "<nel_up> " << nup << " </nel_up>" << endl;
      cout << "<nel_down> " << ndown << " </nel_down>" << endl;
      cout << "<total_spin> " << testdelta << " </total_spin>" << endl;
    }
    if (testdelta != v) {
      if ( ui->oncoutpe() ) {
        int nel = s->wf.nel();
        cout << " <ERROR> spin does not match nel = " << nel << "! </ERROR>" << endl;
      }
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
     st << setw(10) << s->wf.deltaspin();
     return st.str();
  }

  Spin(Sample *sample) : s(sample) { s->ctrl.delta_spin = 0; };
};
#endif
