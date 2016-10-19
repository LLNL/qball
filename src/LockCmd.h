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
// LockCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef LOCKCMD_H
#define LOCKCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"
#include "Species.h"
#include "MMSpecies.h"
#include "Atom.h"

class LockCmd : public Cmd {
  public:

  Sample *s;

  LockCmd(Sample *sample) : s(sample) { };

  char const*name(void) const { return "lock"; }

  char const*help_msg(void) const {
    return 
    "\n lock\n\n"
    " syntax: lock [-rescale|-norescale] [atom name|species name]\n\n"
    "   The lock command fixes either a single atom or all atoms within a species.\n\n";
  }


  int action(int argc, char **argv) {
    string name;
    // atom must be defined with only one or two arguments
    if ( argc < 2 || argc > 3) {
      if ( ui->oncoutpe() )
        cout << "<!-- use: lock [-rescale|-norescale] [atom name|species name] -->" << endl;
      return 1;
    }
  
    name = argv[1];
    bool rescale = true;
    if (name == "-rescale") {
      rescale = true;
      name = argv[2];
    }
    else if (name == "-norescale") {
      rescale = false;
      name = argv[2];
    }

    if ( s->atoms.findSpecies(name)) {
      int isp = s->atoms.isp(name);
      for ( int ia = 0; ia < s->atoms.atom_list[isp].size(); ia++ ) {
        Atom* pa = s->atoms.atom_list[isp][ia];
        pa->lock_atom();
        pa->set_rescale(rescale);
        if ( ui->oncoutpe() )
          cout << "<!-- Atom " << pa->name() << " locked. -->" << endl;
      }
    }
    else if (s->atoms.findMMSpecies(name) ) {
      int isp = s->atoms.isp_mm(name);
      for ( int ia = 0; ia < s->atoms.atom_list[isp].size(); ia++ ) {
        Atom* pa = s->atoms.mmatom_list[isp][ia];
        pa->lock_atom();
        pa->set_rescale(rescale);
        if ( ui->oncoutpe() )
          cout << "<!-- Atom " << pa->name() << " locked. -->" << endl;
      }
    }
    else if (s->atoms.findAtom(name) ) {
      Atom *a = s->atoms.findAtom(name);
      a->lock_atom();
      a->set_rescale(rescale);
      if ( ui->oncoutpe() )
        cout << "<!-- Atom " << a->name() << " locked. -->" << endl;
    }
    else if (s->atoms.findMMAtom(name) ) {
      Atom *a = s->atoms.findMMAtom(name);
      a->lock_atom();
      a->set_rescale(rescale);
      if ( ui->oncoutpe() )
        cout << "<!-- Atom " << a->name() << " locked. -->" << endl;
    }
    else {
      if ( ui->oncoutpe() )
        cout << "<ERROR>LockCmd:  " << name << " not found!</ERROR>" << endl;
      return 1;
    }
    return 0;
  }

};
#endif
