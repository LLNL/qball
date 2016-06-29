////////////////////////////////////////////////////////////////////////////////
//
// CellLock.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellLock.h,v 1.5 2008/04/07 22:00:37 draeger1 Exp $

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
