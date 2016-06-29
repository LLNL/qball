////////////////////////////////////////////////////////////////////////////////
//
// UnlockCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: UnlockCmd.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef UNLOCKCMD_H
#define UNLOCKCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"
#include "Species.h"
#include "MMSpecies.h"
#include "Atom.h"

class UnlockCmd : public Cmd {
  public:

  Sample *s;

  UnlockCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "unlock"; }

  char *help_msg(void) const {
    return 
    "\n unlock\n\n"
    " syntax: unlock [atom name|species name]\n\n"
    "   The unlock command unlocks either a single atom or all atoms within a species.\n\n";
  }

  int action(int argc, char **argv) {
    string name;
    // atom must be defined with only one argument
    if ( argc != 2 ) {
      if ( ui->oncoutpe() )
        cout << "<!-- use: unlock [atom name|species name] -->" << endl;
      return 1;
    }
  
    name = argv[1];
    if ( s->atoms.findSpecies(name)) {
      int isp = s->atoms.isp(name);
      for ( int ia = 0; ia < s->atoms.atom_list[isp].size(); ia++ ) {
        Atom* pa = s->atoms.atom_list[isp][ia];
        pa->unlock_atom();
        if ( ui->oncoutpe() )
          cout << "<!-- Atom " << pa->name() << " unlocked. -->" << endl;
      }
    }
    else if (s->atoms.findMMSpecies(name) ) {
      int isp = s->atoms.isp_mm(name);
      for ( int ia = 0; ia < s->atoms.atom_list[isp].size(); ia++ ) {
        Atom* pa = s->atoms.mmatom_list[isp][ia];
        pa->unlock_atom();
        if ( ui->oncoutpe() )
          cout << "<!-- Atom " << pa->name() << " unlocked. -->" << endl;
      }
    }
    else if (s->atoms.findAtom(name) ) {
      Atom *a = s->atoms.findAtom(name);
      a->unlock_atom();
      if ( ui->oncoutpe() )
        cout << "<!-- Atom " << a->name() << " unlocked. -->" << endl;
    }
    else if (s->atoms.findMMAtom(name) ) {
      Atom *a = s->atoms.findMMAtom(name);
      a->unlock_atom();
      if ( ui->oncoutpe() )
        cout << "<!-- Atom " << a->name() << " unlocked. -->" << endl;
    }
    else {
      if ( ui->oncoutpe() )
        cout << "<ERROR>UnlockCmd:  " << name << " not found!</ERROR>" << endl;
      return 1;
    }
    return 0;
  }

};
#endif
