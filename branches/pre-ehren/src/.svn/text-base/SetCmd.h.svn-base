////////////////////////////////////////////////////////////////////////////////
//
// SetCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SetCmd.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

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

  char *name(void) const { return "set"; }

  char *help_msg(void) const
  {
    return 
    "\n set\n\n"
    " syntax: set variable value[s]\n\n"
    "   The set command sets the value of an interface variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc < 3 && ui->oncoutpe() )
    {
      cout << "<!-- use: set variable value[s] -->" << endl;
      return 1;
    }

    Var* varptr = ui->findVar(argv[1]);

    if ( varptr ) 
    {
      varptr->set(argc-1,&argv[1]);
    }
    else
    {
      // variable is not in the variable list
      if ( ui->oncoutpe() )
        cout << "<WARNING> no such variable: " << argv[1] << " </WARNING>" << endl;
      return 1;
    }
    return 0;
  }
};
#endif
