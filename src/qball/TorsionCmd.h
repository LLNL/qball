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
// TorsionCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef TORSIONCMD_H
#define TORSIONCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class TorsionCmd : public Cmd
{
  public:

  Sample *s;

  TorsionCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "torsion"; }
  char const*help_msg(void) const
  {
    return
    "\n torsion\n\n"
    " syntax: torsion name1 name2 name3 name4\n\n"
    "   The torsion command prints the dihedral defined by four atoms.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 5 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: torsion name1 name2 name3 name4" << endl;
      }
      return 1;
    }

    string name1(argv[1]);
    string name2(argv[2]);
    string name3(argv[3]);
    string name4(argv[4]);
    Atom* a1 = s->atoms.findAtom(name1);
    Atom* a2 = s->atoms.findAtom(name2);
    Atom* a3 = s->atoms.findAtom(name3);
    Atom* a4 = s->atoms.findAtom(name4);
    if ( a1 == 0 || a2 == 0 || a3 == 0 || a4 == 0 )
    {
      if ( ui->oncoutpe() )
      {
        if ( a1 == 0 )
          cout << " TorsionCmd: atom " << name1 << " not defined" << endl;
        if ( a2 == 0 )
          cout << " TorsionCmd: atom " << name2 << " not defined" << endl;
        if ( a3 == 0 )
          cout << " TorsionCmd: atom " << name3 << " not defined" << endl;
        if ( a4 == 0 )
          cout << " TorsionCmd: atom " << name4 << " not defined" << endl;
      }
      return 1;
    }

    if ( a1 == a2 || a1 == a3 || a1 == a4 ||
         a2 == a3 || a2 == a4 || a3 == a4 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " TorsionCmd: replicated atoms in "
             << name1 << " " << name2 << " "
             << name3 << " " << name4 << endl;
      }
      return 1;
    }

    D3vector r12(a1->position()-a2->position());
    D3vector r32(a3->position()-a2->position());
    D3vector r43(a4->position()-a3->position());
    if ( norm(r12) == 0.0 || norm(r32) == 0.0 || norm(r43) == 0.0 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " TorsionCmd: atoms are too close" << endl;
      }
      return 1;
    }

    D3vector e12(normalized(r12));
    D3vector e32(normalized(r32));
    D3vector e23(-e32);
    D3vector e43(normalized(r43));

    const double sin123 = length(e12^e32);
    const double sin234 = length(e23^e43);

    double a = 0;
    if ( sin123 != 0.0 && sin234 != 0.0 )
    {
      D3vector e123 = normalized(e12^e32);
      D3vector e234 = normalized(e23^e43);
      double cc = max(min(e123*e234,1.0),-1.0);
      double ss = -max(min((e123^e234)*e32,1.0),-1.0);
      a = (180.0/M_PI) * atan2(ss,cc);
    }

    if ( ui->oncoutpe() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout << " torsion "
           << name1 << "-" << name2 << "-"
           << name3 << "-" << name4 << ": "
           << setprecision(3) << a << " (deg)" << endl;
    }
    return 0;
  }
};
#endif

// Local Variables:
// mode: c++
// End:
