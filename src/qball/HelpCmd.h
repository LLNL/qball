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
// HelpCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef HELPCMD_H
#define HELPCMD_H

#include <iostream>
#include <iomanip>
#include <string.h>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class HelpCmd : public Cmd
{
  public:

  Sample *s;

  HelpCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "help"; }

  char const*help_msg(void) const
  {
    return 
    "\n help\n\n"
    " syntax: help [command_name]\n\n"
    "   The help command gives a short description of a command. If used\n"
    "   without arguments, help prints a list of valid commands.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( ui->oncoutpe() )
    {
      if ( argc == 1 )  // no arguments
      {
        cout << endl << " valid commands are: " << endl << endl;
        list<Cmd*>::iterator cmd = ui->cmdlist.begin();
        int n = 0;
        while ( cmd != ui->cmdlist.end() )
        {
          n++;
          cout.setf(ios::left,ios::adjustfield);
          cout << " " << setw(16) << (*cmd)->name();
          cout.setf(ios::right,ios::adjustfield);
          if ( n%4 == 0 ) cout << endl;
          cmd++;
        }
        if ( n%4 != 0 ) cout << endl;
        cout << endl;
      }
      else if ( argc == 2 ) // one argument
      {
        // search command list
        Cmd *cmdptr = ui->findCmd(argv[1]);

        if ( cmdptr )
        {
          cout << cmdptr->help_msg();
        }
        else
        {
          cout << " help: " << argv[1] << " is not a valid command" << endl;
        }
      }
      else
      {
        cout << " use: help [command_name]" << endl;
      }
    }
    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
