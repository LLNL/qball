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
// PlotCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PlotCmd.h,v 1.1 2009-06-29 09:59:41 fgygi Exp $

#ifndef PLOTCMD_H
#define PLOTCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class PlotCmd : public Cmd
{
  public:

  Sample *s;

  PlotCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "plot"; }

  char *help_msg(void) const
  {
    return
    "\n plot\n\n"
    " syntax: plot filename\n"
    "         plot -density filename\n"
    "         plot -wf <n> filename\n"
    "         plot -wf <nmin> <nmax> filename\n\n"
    "   The plot command creates a plot file in xyz or cube format.\n\n"
    "   The default format is xyz, used for plotting atoms only.\n"
    "   When using the -density option, the charge density is written\n"
    "   after the atomic positions.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
