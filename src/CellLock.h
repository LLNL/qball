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
// CellLock.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef CELLLOCK_H
#define CELLLOCK_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class CellLock : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "cell_lock"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> cell_lock takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "OFF" || v == "A" || v == "B" || v == "C" ||
            v == "AB" || v == "AC" || v == "BC" || v == "ABC" ||
            v == "S"  || v == "AS" || v == "BS" || v == "CS" ||
            v == "ABS" || v == "ACS" || v == "BCS" || v == "R" ||
            v == "V"  || v == "AV" || v == "BV" || v == "CV" ||
            v == "ABV"  || v == "ACV" || v == "BCV" ||
            v == "SV"  || v == "ASV" || v == "BSV" || v == "CSV" ||
            v == "ABSV"  || v == "ACSV" || v == "BCSV"))
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> cell_lock must be in "
             << "[OFF,A,B,C,AB,AC,BC,S,AS,BS,CS,ABS,ACS,BCS,R,AV,BV,CV,ABV,ACV,BCV,SV,ASV,BSV,CSV,ABSV,ACSV,BCSV]" << endl;
      return 1;
    }

    s->ctrl.cell_lock = v;
        
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_lock;
     return st.str();
  }

  CellLock(Sample *sample) : s(sample) { s->ctrl.cell_lock = "OFF"; }
};
#endif
