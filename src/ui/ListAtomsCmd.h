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
// ListAtomsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef LISTATOMSCMD_H
#define LISTATOMSCMD_H

#include <iostream>
#include <ui/UserInterface.h>
#include <qball/Sample.h>
#include <cstdlib>

class ListAtomsCmd : public Cmd
{
  public:

  Sample *s;

  ListAtomsCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "list_atoms"; }
  char const*help_msg(void) const
  {
    return 
    "\n list_atoms\n\n"
    " syntax: list_atoms\n\n"
    "   The list_atoms command prints a list of all defined atoms.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 1 )
    {
      if ( ui->oncoutpe() )
      {
        cout << "<!-- use: list_atoms -->" << endl;
      }
      return 1;
    }
    s->atoms.listAtoms();
    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
