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
// RandomizeWfCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef RANDOMIZEWFCMD_H
#define RANDOMIZEWFCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeWfCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeWfCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "randomize_wf"; }
  char const*help_msg(void) const
  {
    return 
    "\n randomize_wf\n\n"
    " syntax: randomize_wf [amplitude]\n\n"
    "   The randomize_wf command adds random amplitudes to\n"
    "   the wavefunction Fourier coefficients\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc > 2 )
    {
      if ( ui->oncoutpe() )
      {
        cout << "<!-- use: randomize_wf [amplitude] -->" << endl;
      }
      return 1;
    }
    double amp = 0.02;
    if ( argc == 2 )
      amp = atof(argv[1]);

    // use extra memory for SlaterDets if memory variable = normal, large or huge
    bool highmem = false;
    if (s->ctrl.extra_memory >= 3)
      highmem = true;
    if (s->ctrl.ultrasoft)
      s->wf.randomize_us(amp,s->atoms,highmem);
    else
      s->wf.randomize(amp,highmem);
    
    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
