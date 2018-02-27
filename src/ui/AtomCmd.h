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
// AtomCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef ATOMCMD_H
#define ATOMCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include <ui/UserInterface.h>
#include <qball/Sample.h>
#include "Unit.h"

class AtomCmd : public Cmd
{
  public:

  Sample *s;

  AtomCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "atom"; }
  char const*help_msg(void) const
  {
    return 
    "\n atom\n\n"
    " syntax: atom name species x y z [vx vy vz]\n\n"
    "   The atom command defines a new atom and adds it to the atom list.\n"
    "   The name can be any character string,  the species must be have been\n"
    "   previously defined using the 'species' command.  The position of the\n"
    "   atom is specified by x y and z. Optionally, the atom velocity can be\n"
    "   specified by vx vy and vz.\n\n";
  }

  int action(int argc, char **argv)
  {
    string name;
    string species;
    D3vector position;
    D3vector velocity;
    string pos_unit_name = "bohr";
    string vel_unit_name = "atomicvelocity";

    // atom must be defined with either 3 or 6 arguments
    if (argc != 6 && argc != 9 && argc != 7 && argc != 11) {
      ui->error("<!-- use: atom name species x y z units [vx vy vz units] -->");
      return 1;
    } else if (argc == 6) {
      ui->warning("Units missing for the atom command, assuming 'bohr'.");
    } else if (argc == 9) {
      ui->warning("Units missing for the atom command, assuming 'bohr' and 'atomicvelocity'.");
    }
  
    name = argv[1];
    species = argv[2];

    if(argc == 7 || argc == 11) pos_unit_name = argv[6];

    position.x = atof(argv[3]);
    position.y = atof(argv[4]);
    position.z = atof(argv[5]);

    if(pos_unit_name == "crystal" || pos_unit_name == "direct" || pos_unit_name == "reduced"){
      pos_unit_name = "bohr";
      position = s->atoms.cell().crystal_to_cart(position);
    }
        
    Unit pos_unit(Dimensions::length, pos_unit_name);

    if(!pos_unit.exists()) {
      ui->error("Unknown energy unit '" + pos_unit_name + "'.");
      return 1;
    }

    position = pos_unit.to_atomic(position);
    
    if ( argc == 9 ) {
      velocity.x = atof(argv[6]);
      velocity.y = atof(argv[7]);
      velocity.z = atof(argv[8]);
    }

    if ( argc == 11 ) {
      vel_unit_name = argv[10];
      velocity.x = atof(argv[7]);
      velocity.y = atof(argv[8]);
      velocity.z = atof(argv[9]);
    }

    Unit vel_unit(Dimensions::velocity, vel_unit_name);

    if(!vel_unit.exists()) {
      ui->error("Unknown energy unit '" + vel_unit_name + "'.");
      return 1;
    }
    
    velocity = vel_unit.to_atomic(velocity);
    
    Atom *a = new Atom(name, species, position, velocity);

    const int atoms_nel_before = s->atoms.nel();

    if ( !(s->atoms.addAtom( a ) ) ) {
       if ( ui->oncoutpe() ) cout << " AtomCmd: could not add atom " << name << endl;
       delete a;
       return 1;
    }
    
    const int atoms_nel_after = s->atoms.nel();
    const int delta_nel = (name[0] != '+' ? (atoms_nel_after - atoms_nel_before) : 0);
    const int wf_nel = s->wf.nel();

    if (delta_nel != 0) {
       s->wf.set_nel(wf_nel+delta_nel);
       //ewd s->wf.update_occ(0.0,0);
       if ( s->wfv != 0 )
       {
          s->wfv->set_nel(wf_nel+delta_nel);
          s->wfv->clear();
       }
    } else {
       if ( ui->oncoutpe() )
          cout << " AtomCmd: added ion but no electrons " << endl;
    }

    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
