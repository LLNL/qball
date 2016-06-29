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
// WFPhaseRealCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WFPhaseRealCmd.h,v 1.0 2011-03-25 15:56:17 schleife Exp $

#ifndef WF_PHASE_REAL_CMD_H
#define WF_PHASE_REAL_CMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

// AS: change phase of the wave function to make it real for Gamma only

class WFPhaseRealCmd : public Cmd
{
  public:

  Sample *s;

  WFPhaseRealCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "wf_phase_real"; }
  char *help_msg(void) const
  {
    return
    "\n wf_phase_real\n\n"
    " syntax: wf_phase_real\n\n"
    "   The wf_phase_real command modifies the phase of the wave function so that it is real for Gamma only.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 1 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: wf_phase_real" << endl;
      }
      return 1;
    }

    if ( ui->oncoutpe() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout << setprecision(3)
           << " modifying the phase of the wave functions " << endl;
    }

    s->wf.phase_wf_real();

    return 0;
  }
};
#endif
