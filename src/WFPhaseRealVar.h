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
// WF_Phase_Real.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WF_Phase_Real.h,v 1.0 2011-03-09 15:53:20 schleife Exp $

#include <config.h>

#ifndef WF_PHASE_REAL_H
#define WF_PHASE_REAL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements a status variable to force complex basis [even if kpoint == (0,0,0)]
// AS: necessary for time propagation of wave functions

class WF_Phase_RealVar : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "phase_real_wf"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " phase_real_wf takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " phase_real_wf must be ON or OFF" << endl;
      return 1;
    }

    if ( v == "ON" )
    {
      s->wf.phase_real( true ) ;
    }

    if ( v == "OFF" )
    {
      s->wf.phase_real( false ) ;
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) ;
     if ( s->wf.phase_real_set() ) st << "ON" ; else st << "OFF" ;
     return st.str();
  }

  WF_Phase_RealVar(Sample *sample) : s(sample) { s->wf.phase_real( false ); }
};
#endif
