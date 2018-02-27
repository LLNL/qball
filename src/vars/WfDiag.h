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
// WfDiag.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef WFDIAG_H
#define WFDIAG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class WfDiag : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "wf_diag"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> wf_diag takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "T" || v == "F" || v == "EIGVAL" || v == "MLWF" || v == "MLWFC") )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> wf_diag must be T, F, EIGVAL, MLWF, or MLWFC </ERROR>" << endl;
      return 1;
    }

    s->ctrl.wf_diag = v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.wf_diag;
     return st.str();
  }

  WfDiag(Sample *sample) : s(sample) { s->ctrl.wf_diag = "F"; };
};
#endif

// Local Variables:
// mode: c++
// End:
