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
// SetCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SETCMD_H
#define SETCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SetCmd : public Cmd
{
  public:

  Sample *s;

  SetCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "set"; }

  char const*help_msg(void) const
  {
    return 
    "\n set\n\n"
    " syntax: set variable value[s]\n\n"
    "   The set command sets the value of an interface variable.\n\n";
  }

  int action(int argc, char **argv) {
    if ( argc < 3 && ui->oncoutpe() ) {
      cout << "<!-- use: set variable value[s] -->" << endl;
      return 1;
    }

    Var* varptr = ui->findVar(argv[1]);
    
    if ( varptr ) {
      return varptr->set(argc-1,&argv[1]);
    } else {
      // variable is not in the variable list
      if ( ui->oncoutpe() )
        cout << "<WARNING> no such variable: " << argv[1] << " </WARNING>" << endl;
      return 1;
    }
  }
  
};
#endif

// Local Variables:
// mode: c++
// End:
