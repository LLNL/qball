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
// KpointCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef KPOINTCMD_H
#define KPOINTCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class KpointCmd : public Cmd {
  public:

  Sample *s;

  KpointCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "kpoint"; }
  char const*help_msg(void) const {
    return 
    "\n kpoint\n\n"
    " syntax: kpoint kx ky kz weight\n\n"
    "   The kpoint command defines a new k-point and adds it to the k-point list.\n\n";
  }

  int action(int argc, char **argv) {
    double kx,ky,kz,wt;
  
    // kpoint must be defined with either 3, 4, 5 or 6 arguments
    if ( argc < 4 && argc > 7) {
      if ( ui->oncoutpe() )
        cout << "<!-- use:  kpoint kx ky kz weight [unit keyword] -->" << endl;
      return 1;
    }
  
    kx = atof(argv[1]);
    ky = atof(argv[2]);
    kz = atof(argv[3]);

    wt = 1.0;
    if (argc >= 5)
      wt = atof(argv[4]);

    string units = "crystal";
    if (argc >= 6)
      units = argv[5];

    if ( !(units == "crystal" || units == "CRYSTAL" || units == "cartesian" || units == "CARTESIAN" || units == "Cartesian") ) {
      if (ui->oncoutpe() )
        cout << "<ERROR> kpoint command:  units keyword must be crystal or cartesian </ERROR>" << endl;
      return 1;
    }

    UnitCell wfcell;
    if (s->wf.refcell().volume() > 0.0) 
      wfcell = s->wf.refcell();
    else
      wfcell = s->wf.cell();

    D3vector addkpt;
    if (units == "crystal" || units == "CRYSTAL") {
      addkpt = D3vector(kx,ky,kz);
    }
    else if (units == "cartesian" || units == "CARTESIAN" || units == "Cartesian") {
      double tkx = kx;
      double tky = ky;
      double tkz = kz;
      if (argc == 7) {
        double a = atof(argv[6]);
        //double twopia = 2.0*3.14159265358979/a;
        double twopia = 2.0*M_PI/a;
        tkx *= twopia;
        tky *= twopia;
        tkz *= twopia;
      }
      D3vector a0 = wfcell.a(0);
      D3vector a1 = wfcell.a(1);
      D3vector a2 = wfcell.a(2);
      const double twopiinv = 0.5/M_PI;
      D3vector cartv(tkx,tky,tkz);
      addkpt.x = cartv*a0*twopiinv;
      addkpt.y = cartv*a1*twopiinv;
      addkpt.z = cartv*a2*twopiinv;
    }
    s->wf.add_kpoint(addkpt,wt);

    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
