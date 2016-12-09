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
// RandomizeRealWfCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RandomizeRealWfCmd.h,v 1.5 2008-09-08 15:56:19 fgygi Exp $

#include <config.h>

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

  char const*name(void) const { return "randomize_real_wf"; }
  char const*help_msg(void) const
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

// Local Variables:
// mode: c++
// End:
