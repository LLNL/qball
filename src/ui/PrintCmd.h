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
// PrintCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef PRINTCMD_H
#define PRINTCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include <qball/Sample.h>

class PrintCmd : public Cmd
{
  public:

  Sample *s;

  PrintCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "print"; }

  char const*help_msg(void) const
  {
    return 
    "\n print\n\n"
    " syntax: print variable\n\n"
    "   The print command prints the value of an interface variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( ui->oncoutpe() )
    {
      if ( argc != 2 )
      {
        cout << "<!-- use: print variable -->" << endl;
        return 1;
      }

      Var *varptr = ui->findVar(argv[1]);
      if ( varptr )
      {
        cout << varptr->print() << endl;
      }
      else
      {
        // variable is not in the variable list
        cout << "<WARNING> no such variable: " << argv[1] << " </WARNING>" << endl;
        return 1;
      }
    }
    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
