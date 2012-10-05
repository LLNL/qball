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
// RandomizeRealWfCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RandomizeRealWfCmd.h,v 1.5 2008-09-08 15:56:19 fgygi Exp $

#ifndef RANDOMIZEREALWFCMD_H
#define RANDOMIZEREALWFCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeRealWfCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeRealWfCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "randomize_real_wf"; }
  char *help_msg(void) const
  {
    return
    "\n randomize_real_wf\n\n"
    " syntax: randomize_real_wf [amplitude]\n\n"
    "   The randomize_real_wf command adds random amplitudes to\n"
    "   the wavefunction Fourier coefficients ensuring that it is real!\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc > 2 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: randomize_real_wf [amplitude]" << endl;
      }
      return 1;
    }
    double amp = 0.02;
    if ( argc == 2 )
      amp = atof(argv[1]);
    s->wf.randomize_real(amp);
    return 0;
  }
};
#endif
