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
// MMSpeciesCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

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
