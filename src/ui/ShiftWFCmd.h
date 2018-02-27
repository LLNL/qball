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
// ShiftWFCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ShiftWFCmd.h,v 1.0 2011-03-25 15:56:17 schleife Exp $

#include <config.h>

#ifndef SHIFT_WF_CMD_H
#define SHIFT_WF_CMD_H

#include <iostream>
#include <ui/UserInterface.h>
#include <qball/Sample.h>
#include <cstdlib>

// AS: this class implements a command to shift state n_state by the vector (shift_x, shift_y, shift_z)
// AS: useful, for instance, to test the time propagation schemes

class ShiftWFCmd : public Cmd
{
  public:

  Sample *s;

  ShiftWFCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "shift_wf"; }
  char const*help_msg(void) const
  {
    return
    "\n shift_wf\n\n"
    " syntax: shift_wf x y z n_state\n\n"
    "   The shift_wf command shifts the wave function of state n_state by the given vector.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 5 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: shift_wf x y z n_state" << endl;
      }
      return 1;
    }

    double shift_x = atof(argv[1]);
    double shift_y = atof(argv[2]);
    double shift_z = atof(argv[3]);
    int n_state = atoi(argv[4]);

    if ( ui->oncoutpe() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout << setprecision(3)
           << " shifting wave function of state " << n_state << " by ("
           << shift_x << " , " << shift_y  << " , " << shift_z << ") " << endl;
    }

    // AS: Before we shift we keep a copy of the unshifted ground-state wave functions and use the pointer
    // AS: hamil_wf to make qbox construct the Hamiltonian from these ground-state wave functions (charge densities)
    if ( ( ( s->ctrl.wf_dyn == "TDEULER" ) ||
           ( s->ctrl.wf_dyn == "SOTD" ) ||
           ( s->ctrl.wf_dyn == "SORKTD" ) ||
           ( s->ctrl.wf_dyn == "ETRS" ) ||
           ( s->ctrl.wf_dyn == "AETRS" ) ||
           ( s->ctrl.wf_dyn == "FORKTD" ) ) && (s->hamil_wf == &(s->wf) ) )
    {
      if ( ui->oncoutpe() )
      {
        cout << "ShiftWFCmd::keeping a copy of the unshifted wave function" << endl;
      }
      s->hamil_wf = new Wavefunction(s->wf);
      (*s->hamil_wf) = s->wf;
      if ( ui->oncoutpe() )
      {
        cout << "ShiftWFCmd::Hamiltonian will be constructed using the unshifted wave function" << endl;
      }

      // AS: set the correct occupations of the copied wave function because this is NOT done above
      // AS: Keep in mind: also the eigenvalues are not set yet for hamil_wf!
      (*s->hamil_wf).update_occ(0.0,0);
    }

    s->wf.shift_wf(shift_x,shift_y,shift_z,n_state);

    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
