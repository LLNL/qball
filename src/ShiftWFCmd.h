////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// ShiftWFCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ShiftWFCmd.h,v 1.0 2011-03-25 15:56:17 schleife Exp $

#ifndef SHIFT_WF_CMD_H
#define SHIFT_WF_CMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

// AS: this class implements a command to shift state n_state by the vector (shift_x, shift_y, shift_z)
// AS: useful, for instance, to test the time propagation schemes

class ShiftWFCmd : public Cmd
{
  public:

  Sample *s;

  ShiftWFCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "shift_wf"; }
  char *help_msg(void) const
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
