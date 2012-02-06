////////////////////////////////////////////////////////////////////////////////
//
// MMSpeciesCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MMSpeciesCmd.h,v 1.4 2008/04/07 22:00:37 draeger1 Exp $

#ifndef MMSPECIESCMD_H
#define MMSPECIESCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class MMSpeciesCmd : public Cmd {
  public:

  Sample *s;

  MMSpeciesCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "mmspecies"; }

  char *help_msg(void) const {
    return 
    "\n mmspecies\n\n"
    " syntax: mmspecies name atom_mass\n\n"
    "   The mmspecies command defines an empirical species name.\n\n";
  }

  int action(int argc, char **argv) {

    if ( argc != 3) {
      if ( ui->oncoutpe() )
        cout << "<!-- use: mmspecies name atom_mass -->" << endl;
      return 1;
    }

    double mass = atof(argv[2]);
    MMSpecies* sp = new MMSpecies(s->ctxt_,argv[1],mass);

    if (!s->atoms.addMMSpecies(sp,argv[1])) {
      if ( ui->oncoutpe() )
        cout << "<ERROR> MMSpeciesCmd: could not add species " << argv[1] << "</ERROR>" << endl;
      return 1;
    }
  
    return 0;
  }

};
#endif
