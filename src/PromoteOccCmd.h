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
// PromoteOccCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef PROMOTEOCCCMD_H
#define PROMOTEOCCCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class PromoteOccCmd : public Cmd
{
  public:

  Sample *s;

  PromoteOccCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "promote_occ"; }
  char *help_msg(void) const
  {
    return
    "\n promote_occ\n\n"
    " syntax: promote_occ number_of_electrons origin_orbital destination_orbital\n\n"
    "   The promote_occ command allows you to change the occupation numbers\n"
    "   of specific orbitals. Optionally, one can also specify a spin channel (0 or 1).\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc < 4 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: promote_occ number_of_electrons origin_orbital destination_orbital [spin]" << endl;
      }
      return 1;
    }

    double occ_change;
    int origin_level;
    int destination_level;

    occ_change = atof(argv[1]);
    origin_level = atoi(argv[2]);
    destination_level = atoi(argv[3]);
    int ispin = 0;
    if ( argc == 5 )
       ispin = atoi(argv[4]);
    s->wf.promote_occ(occ_change, origin_level, destination_level, ispin);

    if (argc == 4 && s->wf.nspin() > 1)
    {
       if ( ui->oncoutpe() )
          cout << "PromoteOccCmd:  promote_occ used without spin argument, promoting both spin channels..." << endl;
       ispin = 1;
       s->wf.promote_occ(occ_change, origin_level, destination_level, ispin);
    }
        
    return 0;
  }
};
#endif
