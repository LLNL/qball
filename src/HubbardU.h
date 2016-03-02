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
// HubbardU.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef HUBBARDU_H
#define HUBBARDU_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class HubbardU : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "hubbard_u"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> use:  set hubbard_u [species] [U (Ha)] [orbital l local] </ERROR>" << endl;
      return 1;
    }

    s->ctrl.dft_plus_u = true;
    
    string species = argv[1];
    double uval = atof(argv[2]);
    int lval = atoi(argv[3]);

    if ( s->atoms.findSpecies(species)) {
      int isp = s->atoms.isp(species);
      s->atoms.species_list[isp]->set_hubbard_u(uval,lval);
      if ( ui->oncoutpe() )
        cout << "<!-- Hubbard U = " << uval << " Hartree, l = " << lval << " for species " << isp << ", " << species << " -->" << endl;

      // reinitialize to force computation of hubbard terms
      double rcps = s->atoms.species_list[isp]->rcps();
      s->atoms.species_list[isp]->initialize(rcps);
      
    }
    else {
      cout << " <ERROR> set hubbard_u:  species " << species << " not found! </ERROR>" << endl;
      return 1;
    }
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.dft_plus_u;
     return st.str();
  }

  HubbardU(Sample *sample) : s(sample) { s->ctrl.dft_plus_u = false; };
};
#endif
