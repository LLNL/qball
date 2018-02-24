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
// SaveDenFreq.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SAVEDENFREQ_H
#define SAVEDENFREQ_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class SaveDenFreq : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "savedenfreq"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 && argc != 3)
    {
      if ( ui->oncoutpe() )
      cout << " savedenfreq takes only one or two values" << endl;
      return 1;
    }
    int v = atoi(argv[1]);
    s->ctrl.savedenfreq = v;
    if (argc == 3)
       s->ctrl.savedenfilebase = argv[2];

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << name() << ":  ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.savedenfreq;
     return st.str();
  }

  SaveDenFreq(Sample *sample) : s(sample) {
     s->ctrl.savedenfreq = -1 ;
     s->ctrl.savedenfilebase = "density";
  }
};
#endif

// Local Variables:
// mode: c++
// End:
