////////////////////////////////////////////////////////////////////////////////
//
// Smearing.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Smearing.h,v 1.6 2008/04/07 22:00:37 draeger1 Exp $
//
// set method for smearing occupation numbers, currently either Fermi-Dirac or Gaussian
// (Methfessel and Paxton, PRB 40, 3616 (1989)).  In the latter case, the number of terms 
// in the expansion can be given as an optional second argument.  

#ifndef SMEARING_H
#define SMEARING_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Smearing : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "smearing"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 && argc != 3) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> smearing takes only one or two values </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "fermi" || v == "FERMI" || v == "Fermi" || v == "gaussian" || v == "GAUSSIAN" || v == "Gaussian")  ) {
      if ( ui->oncoutpe() )
        cout << " <ERROR> smearing must be fermi or gaussian </ERROR>" << endl;
      return 1;
    }
    if (v == "FERMI" || v == "Fermi") v = "fermi";
    if (v == "GAUSSIAN" || v == "Gaussian") v = "gaussian";

    s->ctrl.smearing_ngauss = 0;
    if (argc == 3 && v == "gaussian") {
      int n = atoi(argv[2]);
      s->ctrl.smearing_ngauss = n;
    }
    if (v == "fermi")
      s->ctrl.smearing_ngauss = -1;

    s->ctrl.smearing = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.smearing;
     return st.str();
  }

  Smearing(Sample *sample) : s(sample) { 
    s->ctrl.smearing = "fermi"; 
    s->ctrl.smearing_ngauss = -1;    // default: Fermi smearing
  }
};
#endif
