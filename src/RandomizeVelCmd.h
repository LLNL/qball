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
// RandomizeVelCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef RANDOMIZEVELCMD_H
#define RANDOMIZEVELCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeVelCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeVelCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "randomize_vel"; }
  char const*help_msg(void) const
  {
    return 
    "\n randomize_vel\n\n"
    " syntax: randomize_vel [(optional) temperature]\n\n"
    "   The randomize_vel command randomly sets the ion velocities\n"
    "   to achieve a given temperature, either given as an argument\n"
    "   or with the th_temp variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc > 2 )
    {
      if ( ui->oncoutpe() )
      {
        cout << "<!-- use: randomize_vel [(optional) temperature] -->" << endl;
      }
      return 1;
    }
    double temp = s->ctrl.th_temp;
    if ( argc == 2 )
      temp = atof(argv[1]);

    s->atoms.randomize_velocities(temp);
    
    return 0;
  }
};
#endif
