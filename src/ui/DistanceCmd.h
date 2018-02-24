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
// DistanceCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef DISTANCECMD_H
#define DISTANCECMD_H

#include <iostream>
#include "UserInterface.h"
#include <qball/Sample.h>
#include <cstdlib>

class DistanceCmd : public Cmd
{
  public:

  Sample *s;

  DistanceCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "distance"; }
  char const*help_msg(void) const
  {
    return
    "\n distance\n\n"
    " syntax: distance name1 name2\n\n"
    "   The distance command prints the distance between two atoms.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 3 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: distance name1 name2" << endl;
      }
      return 1;
    }

    string name1(argv[1]);
    string name2(argv[2]);
    Atom* a1 = s->atoms.findAtom(name1);
    Atom* a2 = s->atoms.findAtom(name2);
    if ( a1 == 0 || a2 == 0 )
    {
      // either a1 or a2 was not found
      if ( ui->oncoutpe() )
      {
        if ( a1 == 0 )
          cout << " DistanceCmd: atom " << name1 << " not defined"
               << endl;
        if ( a2 == 0 )
          cout << " DistanceCmd: atom " << name2 << " not defined"
               << endl;
      }
      return 1;
    }

    if ( ui->oncoutpe() )
    {
      const double d = length(a1->position()-a2->position());
      cout.setf(ios::fixed,ios::floatfield);
      cout << " distance " << name1 << "-" << name2 << ": "
           << setprecision(3)
           << d << " (a.u.) / " << 0.529177*d << " (Ang)" << endl;
    }
    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
