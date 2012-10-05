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
// TDDt.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDDt.h,v 1.5 2011-04-04 15:56:18 schleife Exp $

#ifndef TDDT_H
#define TDDT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements tddt as status variable
// AS: which specifies the dt for the wave function time propagation

class TDDt : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "TD_dt"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " TD_dt takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v == 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " TD_dt must be non-zero" << endl;
      return 1;
    }

    s->ctrl.tddt = v;
    return 0;
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

  TDDt(Sample *sample) : s(sample) { s->ctrl.tddt = 0.01; }
};
#endif
