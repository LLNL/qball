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

#include <ui/UserInterface.h>
#include <qball/Sample.h>
#include "Unit.h"
#include <qball/Atom.h>

class CoordinatesCmd : public Cmd
{

private:

  enum Format { XYZ, POSCAR };
  
  public:

  Sample *s;

  CoordinatesCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "coordinates"; }
  char const*help_msg(void) const
  {
    return 
    "\n coordinates\n\n"
    " syntax: coordinates filename [units]\n\n"
    " The coordinates command adds atoms from an xyz coordinate file.\n"
    " Optionally,  it is possible to specify the units of the coordi-\n"
    " nates (angstrom by default).  You can specify coordinates rela-\n"
    " tive to the unit cell with the 'crystal' units.\n\n";
  }

  int action(int argc, char **argv){

    Format file_format = XYZ;   
    
    string filename(argv[1]);
    
    // coordinates must be defined with either 2 or 3 arguments
    if (argc != 2 && argc != 3) {
      ui->error("<!-- use: coordinates filename [units] -->");
      return 1;
    }
    
    string pos_unit_name = "angstrom";
    bool crystal_units = false;
      
    if(argc == 3){
      pos_unit_name = argv[2];     
    }

    if(pos_unit_name == "crystal" || pos_unit_name == "direct" || pos_unit_name == "reduced"){
      crystal_units = true;
      pos_unit_name = "bohr";
    }
    
    std::ifstream file(filename);

    if(!file){
      ui->error("CoordinateCmd: Cannot open coordinates file '" + filename + "'.");
      return 1;
    }

    int natoms;
    file >> natoms;

    if(ui->oncoutpe()) {
      cout << "CoordinateCmd: Adding " <<  natoms <<  " atoms from coordinate file '" + filename + "'." << endl;
    }

    string dump_line;
    
    getline(file, dump_line);
    getline(file, dump_line);
    
    Unit pos_unit(Dimensions::length, pos_unit_name);
    
    for(int iatom = 0; iatom < natoms; iatom++){
      string species;
      D3vector position;
      D3vector velocity(0.0, 0.0, 0.0);

      file >> species >> position.x >> position.y >> position.z;
      getline(file, dump_line);
      
      if(crystal_units) position = s->atoms.cell().crystal_to_cart(position);

      position = pos_unit.to_atomic(position);
      
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
