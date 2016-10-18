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
// MMAtomCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef MMATOMCMD_H
#define MMATOMCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class MMAtomCmd : public Cmd {
  public:

  Sample *s;

  MMAtomCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "mmatom"; }
  char const*help_msg(void) const {
    return 
    "\n mmatom\n\n"
    " syntax: mmatom name species x y z [vx vy vz]\n\n"
    "   The mmatom command defines a new classical atom and adds it to the atom\n"
    "   list.  The name can be any character string, the species must be the name\n"
    "   of a species defined with the mmspecies command.\n\n";
  }

  int action(int argc, char **argv) {
    string name;
    string species;
    D3vector position;
    D3vector velocity;
  
    // atom must be defined with either 3 or 6 arguments
    if ( argc != 6 && argc != 9 ) {
      if ( ui->oncoutpe() )
        cout << "<!-- use: mmatom name species x y z [vx vy vz] -->" << endl;
      return 1;
    }
  
    name = argv[1];
    species = argv[2];
    position.x = atof(argv[3]);
    position.y = atof(argv[4]);
    position.z = atof(argv[5]);
    if ( argc == 9 ) {
      velocity.x = atof(argv[6]);
      velocity.y = atof(argv[7]);
      velocity.z = atof(argv[8]);
    }
  
    Atom *a = new Atom(name,species,position,velocity);
    
    if ( !(s->atoms.addMMAtom( a ) ) ) {
      if ( ui->oncoutpe() )
        cout << "<ERROR> MMAtomCmd: could not add atom " << name << "</ERROR>" << endl;
      delete a;
      return 1;
    }
    
    return 0;
  }
};
#endif
