////////////////////////////////////////////////////////////////////////////////
//
// FermiTemp.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FermiTemp.h,v 1.6 2008/04/07 22:00:37 draeger1 Exp $
//
// Deprecated variable:  included to provide backward compatibility with older Qbox
// input files.  Fermi_temp inputs smearing in Kelvin, converted to Ry to be compatible
// with smearing and smearing_width variables.

#ifndef FERMITEMP_H
#define FERMITEMP_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class FermiTemp : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "fermi_temp"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2) {
      if ( ui->oncoutpe() )
      cout << "<ERROR> fermi_temp takes only one value </ERROR>" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);

    const double boltz = 1.0 / ( 11605.0 * 13.6058 ); // convert from Kelvin to Ry

    if ( ui->oncoutpe() ) {
      cout << "<!-- Note:  fermi_temp variable is deprecated.  Automatically converting input" << " -->" << endl;
      cout << "<!-- to the following commands instead:" << " -->" << endl;
      cout << "<!--   set smearing fermi" << " -->" << endl;
      cout << "<!--   set smearing_width " << v*boltz << " -->" << endl;
    }

    s->ctrl.smearing = "fermi";
    s->ctrl.smearing_width = v*boltz;
    s->ctrl.smearing_ngauss = -1;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.smearing_width;
     return st.str();
  }

  FermiTemp(Sample *sample) : s(sample) { }
};
#endif
