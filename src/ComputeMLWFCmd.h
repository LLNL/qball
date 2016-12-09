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
// ComputeMLWFCmd.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef COMPUTEMLWFCMD_H
#define COMPUTEMLWFCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"
#include "MLWFTransform.h"

class ComputeMLWFCmd : public Cmd
{
  private:

  MLWFTransform* mlwft;

  public:
  Sample *s;

  ComputeMLWFCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "compute_mlwf"; }

  char const*help_msg(void) const
  {
    return
    "\n compute_mlwf\n\n"
    " syntax: compute_mlwf\n\n"
    "   The compute_mlwf command computes maximally localized \n"
    " Wannier functions.\n\n";
  }

  int action(int argc, char **argv);

  ComputeMLWFCmd();
  ~ComputeMLWFCmd();
};
#endif

// Local Variables:
// mode: c++
// End:
