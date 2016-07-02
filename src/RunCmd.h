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
// RunCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef RUNCMD_H
#define RUNCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"

class RunCmd : public Cmd
{
  private:

  public:

  Sample *s;

  RunCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "run"; }
  char *help_msg(void) const
  {
    return 
    "\n run\n\n"
    " syntax: run n [nite]\n\n"
    "   The run command runs n steps of simulation. Each step\n"
    "   consists of one (optionally nite) electronic steps\n\n";
  }

  int action(int argc, char **argv);

};
#endif
