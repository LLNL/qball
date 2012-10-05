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
// Non_Selfconsistent_Energy_Output.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Force_Complex_WF.h,v 1.0 2011-05-27 11:53:20 schleife Exp $

#ifndef NON_SELFCONSISTENT_ENERGY_OUTPUT_H
#define NON_SELFCONSISTENT_ENERGY_OUTPUT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements a status variable to control the calculation of the total energy
// AS: also during non-selfconsistent electronic steps
// AS: turning this variable to ON increases the computational cost by a factor of two

class Non_Selfconsistent_Energy_Output : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "non_selfc_energy"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " non_selfc_energy takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->oncoutpe() )
        cout << " non_selfc_energy must be ON or OFF" << endl;
      return 1;
    }

    if ( v == "ON" )
    {
      s->ctrl.non_selfc_energy = true;
    }

    if ( v == "OFF" )
    {
      s->ctrl.non_selfc_energy = false;
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
     if ( s->ctrl.non_selfc_energy ) st << "ON" ; else st << "OFF" ;
     return st.str();
  }

  Non_Selfconsistent_Energy_Output(Sample *sample) : s(sample) { s->ctrl.non_selfc_energy=false; }
};
#endif
