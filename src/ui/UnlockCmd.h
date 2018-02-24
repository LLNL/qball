////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// UnlockCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef UNLOCKCMD_H
#define UNLOCKCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include <qball/Sample.h>
#include <qball/Species.h>
#include "MMSpecies.h"
#include <qball/Atom.h>

class UnlockCmd : public Cmd {
  public:

  Sample *s;

  UnlockCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "unlock"; }

  char const*help_msg(void) const {
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

// Local Variables:
// mode: c++
// End:
