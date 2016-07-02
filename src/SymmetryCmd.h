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
// SymmetryCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SYMMETRYCMD_H
#define SYMMETRYCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SymmetryCmd : public Cmd {
  public:

  Sample *s;

  SymmetryCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "symmetry"; }
  char *help_msg(void) const {
    return 
    "\n symmetry\n\n"
    " syntax: symmetry s11 s12 s13 s21 s22 s23 s31 s32 s33 [f1 f2 f3]\n\n"
    "   The symmetry command defines a new symmetry operation over which the charge\n"
    "   density will be averaged.  Symmetry operations may include a fractional\n"
    "   translation, (f1,f2,f3).  All s values should be doubles, such that the operation\n"
    "   can act on real-space density grid points.  Fractional translations should be\n"
    "   doubles between 0.0 and 1.0.\n\n";
  }

  int action(int argc, char **argv) {

    // symmetry must be defined with either 10 or 13 arguments
    if ( argc != 10 && argc != 13 ) {
      if ( ui->oncoutpe() )
        cout << "<!-- use: symmetry s11 s12 s13 s21 s22 s23 s31 s32 s33 [f1 f2 f3] -->" << endl;
      return 1;
    }
  
    double ssin[9];
    double ftransin[3];
    for (int i=0; i<3; i++)
      ftransin[i] = 0.;

    for (int i=0; i<9; i++)
      ssin[i] = atof(argv[i+1]);
    if (argc == 13)
      for (int i=0; i<3; i++)
        ftransin[i] = atof(argv[i+10]);

    if (ssin[0] == 1. && ssin[4] == 1. && ssin[8] == 1. && ssin[1] == 0. && ssin[2] == 0. && ssin[3] == 0. && ssin[5] == 0. && ssin[6] == 0. && ssin[7] == 0. && ftransin[0] == 0. && ftransin[1] == 0. && ftransin[2] == 0.) {
      if ( ui->oncoutpe()) {
        cout << "<!-- Symmetry:  unity operator already included by default, skipping. -->" << endl;
      }
      return 0;
    }

    Symmetry *sym = new Symmetry(ssin[0],ssin[1],ssin[2],ssin[3],ssin[4],ssin[5],ssin[6],ssin[7],ssin[8],ftransin[0],ftransin[1],ftransin[2]);
    
    if ( !(s->symmetries.addSymmetry( sym ) ) ) {
      if ( ui->oncoutpe() )
        cout << "<ERROR> SymmetryCmd: could not add symmetry </ERROR>" << endl;
      delete sym;
      return 1;
    }

    return 0;
  }
};
#endif
