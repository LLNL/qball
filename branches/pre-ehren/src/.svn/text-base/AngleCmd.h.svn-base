////////////////////////////////////////////////////////////////////////////////
//
// AngleCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AngleCmd.h,v 1.2 2009/03/27 00:53:24 draeger1 Exp $

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

  char *name(void) const { return "angle"; }
  char *help_msg(void) const
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
