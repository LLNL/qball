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
// SavesysCmd.C:
//
////////////////////////////////////////////////////////////////////////////////


#include "SavesysCmd.h"
#include "fstream"
#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"

#include "UnitCell.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SavesysCmd::action(int argc, char **argv) {

  if ( !(argc == 2) ) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: savesys filename -->" << endl;
    return 1;
  }

  if (s->ctrl.timer_hit) {
    if (s->ctrl.timer_savesyscmd) {
      if ( ui->oncoutpe() )
        cout << " <!-- SavesysCmd: run_timer exceeded: all savesys commands beyond first will be ignored. -->" << endl;
      return 0;
    }
    else
      s->ctrl.timer_savesyscmd = true;
  }

  char* filename = argv[1];
  if (ui->oncoutpe() ){
    
    ofstream os;
    os.open(filename,ofstream::out);

    // cell info
    string cmd("set cell ");
    s->wf.cell().printsys(os,cmd);

    // ref cell info, if necessary
    if ( s->wf.refcell().volume() != 0.0 ) {
      string refcmd("set ref_cell ");
      s->wf.refcell().printsys(os,refcmd);
    }

    // species info
    const int nspqm_ = s->atoms.nsp();
    for (int i=0; i<nspqm_; i++)
      s->atoms.species_list[i]->printsys(os);

    const int nspmm_ = s->atoms.nsp_mm();
    for (int i=0; i<nspmm_; i++)
      s->atoms.mmspecies_list[i]->printsys(os);

    // atom coordinates and info
    s->atoms.printsys(os);
    
    os.close();
  }



  return 0;
}
