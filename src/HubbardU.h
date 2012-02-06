////////////////////////////////////////////////////////////////////////////////
//
// HubbardU.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HubbardU.h,v 1.2 2010/08/30 20:55:43 draeger1 Exp $

#ifndef HUBBARDU_H
#define HUBBARDU_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class HubbardU : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "hubbard_u"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> use:  set hubbard_u [species] [U (eV)] [orbital l local] </ERROR>" << endl;
      return 1;
    }

    s->ctrl.dft_plus_u = true;
    
    string species = argv[1];
    double uval = atof(argv[2]);
    int lval = atoi(argv[3]);

    if ( s->atoms.findSpecies(species)) {
      int isp = s->atoms.isp(species);
      s->atoms.species_list[isp]->set_hubbard_u(uval,lval);
      if ( ui->oncoutpe() )
        cout << "<!-- Hubbard U = " << uval << " Hartree for species " << species << " -->" << endl;
    }
    else {
      cout << " <ERROR> set hubbard_u:  species " << species << " not found! </ERROR>" << endl;
      return 1;
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.dft_plus_u;
     return st.str();
  }

  HubbardU(Sample *sample) : s(sample) { s->ctrl.dft_plus_u = false; };
};
#endif
