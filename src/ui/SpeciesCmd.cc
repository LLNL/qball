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
// SpeciesCmd.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "SpeciesCmd.h"
#include "SpeciesReader.h"
#include <qball/Species.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SpeciesCmd::action(int argc, char **argv) {
  if (! (argc == 3 || argc == 4)) {
    if ( ui->oncoutpe() )
      cout << "  <!-- use: species name uri [ewald_width] -->" << endl;
    return 1;
  }
  
  if ( ui->oncoutpe() )
    cout << "  <!-- SpeciesCmd: defining species " << argv[1]
         << " as " << argv[2] << " -->" << endl;

  SpeciesReader sp_reader(s->ctxt_);
  
  Species* sp = new Species(s->ctxt_,argv[1]);
  
  try {
    sp_reader.readSpecies(*sp,argv[2]);
    sp_reader.bcastSpecies(*sp);
    if (argc == 4) {
      const double rcpsin = atof(argv[3]);
      s->atoms.addSpecies(sp,argv[1],rcpsin);
    }
    else
      s->atoms.addSpecies(sp,argv[1]);

    if (sp->ultrasoft()) {
      s->ctrl.ultrasoft = true;
      s->wf.set_ultrasoft(true);
    }
    
    if (sp->nlcc()) {
       s->ctrl.nlcc = true;
    }
  }
  catch ( const SpeciesReaderException& e ) {
    cout << " SpeciesReaderException caught in SpeciesCmd" << endl;
    cout << " SpeciesReaderException: cannot define Species" << endl;
  }
  catch (...) {
    cout << " SpeciesCmd: cannot define Species" << endl;
  }
  
  return 0;
}
