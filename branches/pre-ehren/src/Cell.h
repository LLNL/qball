////////////////////////////////////////////////////////////////////////////////
//
// Cell.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Cell.h,v 1.5 2009/12/18 18:40:13 draeger1 Exp $

#ifndef CELL_H
#define CELL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Cell : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "cell"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 10 )
    {
      if ( ui->oncoutpe() )
      cout << "<ERROR> cell must be specified with 3 vectors (9 values) </ERROR>" << endl;
      return 1;
    }
    
    D3vector a0(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector a1(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector a2(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a0,a1,a2);
    
    if ( cell.volume() < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << "<ERROR> cell volume must be positive </ERROR>" << endl;
      return 1;
    }

    if (s->wf.hasdata()) {
      s->wf.resize(cell,s->wf.refcell(),s->wf.ecut());
      if ( s->wfv != 0 )
      {
        s->wfv->resize(cell,s->wf.refcell(),s->wf.ecut());
        s->wfv->clear();
      }
    }
    else {
      s->wf.set_cell(cell);
      if ( s->wfv != 0 )
        s->wfv->set_cell(cell);
    }

    s->atoms.set_cell(a0,a1,a2);
    
    if ( ui->oncoutpe() )
    {
      cout << "  <unitcell>\n"
           << s->wf.cell() 
           << "  </unitcell>" << endl;
    }
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.cell();
     return st.str();
  }

  Cell(Sample *sample) : s(sample) {};
};
#endif
