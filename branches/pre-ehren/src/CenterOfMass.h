////////////////////////////////////////////////////////////////////////////////
//
// CenterOfMassc.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CenterOfMass.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef CENTEROFMASS_H
#define CENTEROFMASS_H

#include <string>
#include <iostream>
#include <cstdlib>

#include "Sample.h"

// CenterOfMass variable determines whether or not to subtract the
// center of mass velocity from particle velocities.  The default
// value is "unset", allowing us to turn it on as a default with
// atoms_dyn = MD while still giving the user the chance to turn it
// off by setting this variable to free before setting atoms_dyn MD.

class CenterOfMass : public Var {
  Sample *s;

  public:

  char *name ( void ) const { return "center_of_mass"; };

  int set ( int argc, char **argv ) {
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
      cout << " <ERROR> center_of_mass takes only one value </ERROR>" << endl;
      return 1;
    }
    
    string v = argv[1];

    if (v == "FIXED" || v=="fixed") { 
      v = "fixed";
    }
    else if (v == "FREE" || v=="free") { 
      v = "free";
    }
    else {
      if ( ui->oncoutpe() )
        cout << " <ERROR> center_of_mass must be set to either free or fixed </ERROR>" << endl;
      return 1;
    }

    s->ctrl.center_of_mass = v;
    return 0;
  }

  string print (void) const {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.center_of_mass;
     return st.str();
  }

  CenterOfMass(Sample *sample) : s(sample) { s->ctrl.center_of_mass = "unset"; }
};
#endif
