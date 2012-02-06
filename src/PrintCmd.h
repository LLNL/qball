////////////////////////////////////////////////////////////////////////////////
//
// PrintCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PrintCmd.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef PRINTCMD_H
#define PRINTCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"

class PrintCmd : public Cmd
{
  public:

  Sample *s;

  PrintCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "print"; }

  char *help_msg(void) const
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
