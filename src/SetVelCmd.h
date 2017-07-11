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
// SetVelCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SETVELCMD_H
#define SETVELCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"
#include "Species.h"
#include "MMSpecies.h"
#include "Atom.h"

class SetVelCmd : public Cmd {
  public:

  Sample *s;

  SetVelCmd(Sample *sample) : s(sample) { };

  char const*name(void) const { return "set_velocity"; }

  char const*help_msg(void) const {
    return 
    "\n lock\n\n"
    " syntax: set_velocity [atom name] [vx] [vy] [vz]\n\n"
    "   The set_velocity command sets the velocity for a given atom.\n\n";
  }


  int action(int argc, char **argv) {
    string name;

    string vel_unit_name = "atomicvelocity";
 
    if ( argc != 5 && argc != 6) {
      ui->error("<!-- use: set_velocity [atom name] [vx] [vy] [vz] [units] -->");
      return 1;
    }

    if (argc == 5) {
      ui->warning("Units missing for the set_velocity command, assuming 'atomicvelocity'.");
    } else {
      vel_unit_name = string(argv[5]);
    }
    
    Unit vel_unit(Dimensions::velocity, vel_unit_name);
    
    if(!vel_unit.exists()) {
      ui->error("Unknown energy unit '" + vel_unit_name + "'.");
      return 1;
    }
    
    name = argv[1];
    double vx = atof(argv[2]);
    double vy = atof(argv[3]);
    double vz = atof(argv[4]);
    D3vector newvel(vx,vy,vz);

    newvel = vel_unit.to_atomic(newvel);
    
    if (s->atoms.findAtom(name) ) {
      Atom *a = s->atoms.findAtom(name);

      a->set_velocity(newvel);
      if ( ui->oncoutpe() )
        cout << "<!-- Atom " << a->name() << " velocity set to " << newvel << " -->" << endl;
    }
    else if (s->atoms.findMMAtom(name) ) {
      Atom *a = s->atoms.findMMAtom(name);

      a->set_velocity(newvel);
      if ( ui->oncoutpe() )
        cout << "<!-- Atom " << a->name() << " velocity set to " << newvel << " -->" << endl;
    }
    else {
      if ( ui->oncoutpe() )
        cout << "<ERROR>SetVelCmd:  " << name << " not found!</ERROR>" << endl;
      return 1;
    }
    return 0;
  }

};
#endif

// Local Variables:
// mode: c++
// End:
