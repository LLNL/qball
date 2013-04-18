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
// Print_Density_Every.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Print_Density_Every.h,v 1.5 2011-05-25 15:56:18 schleife Exp $

#ifndef PRINT_DENSITY_EVERY
#define PRINT_DENSITY_EVERY

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements a function to print the density every N number of MD steps

class Print_Density_Every : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "print_density_every"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " print_density_every takes only one value" << endl;
      return 1;
    }

    int v = atoi(argv[1]);

    s->ctrl.print_density_every = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << name() << ":  ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.print_density_every;
     return st.str();
  }

  Print_Density_Every(Sample *sample) : s(sample) { s->ctrl.print_density_every = -1 ; }
};
#endif
