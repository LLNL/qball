////////////////////////////////////////////////////////////////////////////////
//
// Ecut.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Ecut.h,v 1.4 2009/12/18 18:40:13 draeger1 Exp $

#ifndef ECUT_H
#define ECUT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Ecut : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "ecut"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> ecut takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> ecut must be non-negative </ERROR>" << endl;
      return 1;
    }
    
    if ( s->wf.ecut() == 0.5 * v )
      return 0;

    if (s->wf.hasdata()) {
      s->wf.resize(0.5*v);
      if ( s->wfv != 0 ) 
      {
        s->wfv->resize(0.5*v);
        s->wfv->clear();
      }
    }
    else {
      s->wf.set_ecut(0.5*v);
      if ( s->wfv != 0 )
        s->wfv->set_ecut(0.5*v);

      /*
      //ewd:  add this to avoid users running without initializing wf
      //ewd:  note:  this can cause problems if users set ecut before calling .sys file
      double amp = 0.02;
      bool highmem = false;
      if (s->ctrl.extra_memory >= 3)
         highmem = true;
      if (s->ctrl.ultrasoft)
         s->wf.randomize_us(amp,s->atoms,highmem);
      else
         s->wf.randomize(amp,highmem);
      */
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->wf.ecut();
     return st.str();
  }

  Ecut(Sample *sample) : s(sample) {};
};
#endif
