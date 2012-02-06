////////////////////////////////////////////////////////////////////////////////
//
// RefCell.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RefCell.h,v 1.4 2009/12/18 18:40:13 draeger1 Exp $

#ifndef REFCELL_H
#define REFCELL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class RefCell : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "ref_cell"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 10 )
    {
      if ( ui->oncoutpe() )
      cout << "<ERROR> ref_cell must be specified with 3 vectors (9 values) </ERROR>" << endl;
      return 1;
    }
    
    D3vector a0(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector a1(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector a2(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell ref_cell(a0,a1,a2);
    
    if ( ref_cell.volume() < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << "<ERROR> ref_cell volume must be positive </ERROR>" << endl;
      return 1;
    }

    if (s->wf.hasdata()) {
      s->wf.resize(s->wf.cell(), ref_cell,s->wf.ecut());
      if ( s->wfv != 0 )
      {
        s->wfv->resize(s->wf.cell(),s->wf.refcell(),s->wf.ecut());
        s->wfv->clear();
      }
    }
    else {
      s->wf.set_refcell(ref_cell);
      if ( s->wfv != 0 )
        s->wfv->set_refcell(ref_cell);
    }
    
    if ( ui->oncoutpe() )
    {
      cout << "  <reference_unit_cell>\n"
           << s->wf.refcell() 
           << "  </reference_unit_cell>" << endl;
    }
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.refcell();
     return st.str();
  }

  RefCell(Sample *sample) : s(sample) {}
};
#endif
