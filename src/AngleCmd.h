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
// AngleCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef ANGLECMD_H
#define ANGLECMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class AngleCmd : public Cmd
{
  public:

  Sample *s;

  AngleCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "angle"; }
  char const*help_msg(void) const
  {
    return
    "\n angle\n\n"
    " syntax: angle name1 name2 name3\n\n"
    "   The angle command prints the angle defined by three atoms.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 4 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: angle name1 name2 name3" << endl;
      }
      return 1;
    }

    string name1(argv[1]);
    string name2(argv[2]);
    string name3(argv[3]);
    Atom* a1 = s->atoms.findAtom(name1);
    Atom* a2 = s->atoms.findAtom(name2);
    Atom* a3 = s->atoms.findAtom(name3);
    if ( a1 == 0 || a2 == 0 || a3 == 0 )
    {
      if ( ui->oncoutpe() )
      {
        if ( a1 == 0 )
          cout << " AngleCmd: atom " << name1 << " not defined" << endl;
        if ( a2 == 0 )
          cout << " AngleCmd: atom " << name2 << " not defined" << endl;
        if ( a3 == 0 )
          cout << " AngleCmd: atom " << name3 << " not defined" << endl;
      }
      return 1;
    }

    if ( a1 == a2 || a2 == a3 || a3 == a1 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " AngleCmd: replicated atoms in " << name1
             << " " << name2 << " " << name3 << endl;
      }
      return 1;
    }

    D3vector r12(a1->position()-a2->position());
    D3vector r32(a3->position()-a2->position());
    if ( norm(r12) == 0.0 || norm(r32) == 0.0 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " AngleCmd: atoms are too close" << endl;
      }
      return 1;
    }

    const double sp = normalized(r12) * normalized(r32);
    const double c = max(-1.0,min(1.0,sp));
    const double a = (180.0/M_PI)*acos(c);
    if ( ui->oncoutpe() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout << " angle " << name1 << "-" << name2  << "-" << name3
           << ": "
           << setprecision(3)
           << a << " (deg)" << endl;
    }
    return 0;
  }
};
#endif
