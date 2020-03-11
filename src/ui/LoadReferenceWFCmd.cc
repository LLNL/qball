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
// LoadReferenceWFCmd.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "LoadReferenceWFCmd.h"
#include <qball/Sample.h>
#include <qball/Context.h>
#include <qball/ChargeDensity.h>
#include <qball/FourierTransform.h>
#include <qball/Basis.h>
#include "fstream"
#include <iostream>
#include <iomanip>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int LoadReferenceWFCmd::action(int argc, char **argv) {

  if ( !(argc>=2 && argc<=3 ) ) 
  {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: load -states filename -->" << endl;
    return 1;
  }
  
  char* filename = 0;
  string encoding = "states"; 

  // parse arguments
  for ( int i = 1; i < argc; i++ ) 
  {
    string arg(argv[i]);
    
    if ( arg=="-states" )
      encoding = "states";
    else if ( arg[0] != '-' && i == argc-1 )
      filename = argv[i];
    else 
    {
      if ( ui->oncoutpe() )
        cout << "  <!-- use: load -states filename -->" 
             << endl;
      return 1;
    }
  }
  
  if ( filename == 0 ) 
  {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: load -states  filename -->" 
           << endl;
    return 1;
  }

  string filestr(filename);

  if (encoding == "states" ) 
  {
     s->wf.read_states(filestr);
     s->previous_wf = new Wavefunction(s->wf);
     *(s->previous_wf) = s->wf;
     (*(s->previous_wf)).update_occ(0.0,0);
  }
}
