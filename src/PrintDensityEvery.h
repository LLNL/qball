////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
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
