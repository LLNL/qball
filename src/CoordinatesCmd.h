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

#ifndef COORDINATESCMD_H
#define COORDINATESCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <cstdlib>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"
#include "Unit.h"
#include "Atom.h"

class CoordinatesCmd : public Cmd
{
  public:

  Sample *s;

  CoordinatesCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "coordinates"; }
  char const*help_msg(void) const
  {
    return 
    "\n coordinates\n\n"
    " syntax: coordinates filename\n\n"
    " The coordinates command adds atoms from an xyz coordinate file.\n\n";
  }

  int action(int argc, char **argv){

    string pos_unit_name = "angstrom";

    string filename(argv[1]);
    
    // coordinates must be defined with either 3 or 6 arguments
    if (argc != 2) {
      ui->error("<!-- use: coordinates filename -->");
      return 1;
    }

    std::ifstream file(filename);

    if(!file){
      ui->error("CoordinateCmd: Cannot open coordinates file '" + filename + "'.");
      return 1;
    }

    int natoms;

    file >> natoms;

    cout << "CoordinateCmd: Adding " <<  natoms <<  " atoms from coordinate file '" + filename + "'." << endl;

    string comment_line;
    
    getline(file, comment_line);
    cout << comment_line << endl;
    getline(file, comment_line);
    cout << comment_line << endl;  

    Unit pos_unit(Dimensions::length, pos_unit_name);
    
    for(int iatom = 0; iatom < natoms; iatom++){
      string species;
      D3vector position;
      D3vector velocity(0.0, 0.0, 0.0);
 
      file >> species >> position.x >> position.y >> position.z;

      position = pos_unit.to_atomic(position);
      
      cout << species << ' ' << position << endl;

      stringstream atom_name;
      atom_name << species << iatom;

      const int atoms_nel_before = s->atoms.nel();
 
      Atom *a = new Atom(atom_name.str(), species, position, velocity);

      if ( !(s->atoms.addAtom( a ) ) ) {
	delete a;
	ui->error("CoordinateCmd: could not add atom '" + atom_name.str() + "' from file '" + filename + "'.");
	return 1;
      }
      
      const int atoms_nel_after = s->atoms.nel();
      const int delta_nel = (species[0] != '+' ? (atoms_nel_after - atoms_nel_before) : 0);
      const int wf_nel = s->wf.nel();
      
      if (delta_nel != 0) {
	s->wf.set_nel(wf_nel+delta_nel);

       if ( s->wfv != 0 )
	 {
	   s->wfv->set_nel(wf_nel+delta_nel);
	   s->wfv->clear();
	 }
      } else {
	if ( ui->oncoutpe() )
          cout << " AtomCmd: added ion but no electrons " << endl;
      }
      
    }
    
    return 0;
    
  }
};
#endif

// Local Variables:
// mode: c++
// End:
