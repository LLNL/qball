////////////////////////////////////////////////////////////////////////////////
//
// NA_overlaps.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NA_overlaps.h,v 1.5 2011-05-25 15:56:18 schleife Exp $

#ifndef NA_OVERLAPS_H
#define NA_OVERLAPS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

// AS: this class implements na_overlap_min and na_overlap_max as status variables
// AS: they specify the minimum and maximum band index to calculate (non-adiabatic) overlaps for

class NA_overlaps : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "NA_overlaps"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 3 )
    {
      if ( ui->oncoutpe() )
      cout << " NA_overlaps takes only two values" << endl;
      return 1;
    }

    int v = atoi(argv[1]);
    int w = atoi(argv[2]);
    if ( w <= v )
    {
      if ( ui->oncoutpe() )
        cout << " second band index must be larger than first one" << endl;
      return 1;
    }

    if ( ( w < 0 ) or (v < 0) )
    {
      if ( ui->oncoutpe() )
        cout << " both band indices must be larger than or equal to zero" << endl;
      return 1;
    }

    s->ctrl.na_overlap_min = v;
    s->ctrl.na_overlap_max = w;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << name() << ": min = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.na_overlap_min << " ; max = " << s->ctrl.na_overlap_max;
     return st.str();
  }

  NA_overlaps(Sample *sample) : s(sample) { s->ctrl.na_overlap_min = -1 ; s->ctrl.na_overlap_max = -1 ; }
};
#endif
