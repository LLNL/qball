////////////////////////////////////////////////////////////////////////////////
//
// ListConstraintsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ListConstraintsCmd.h,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#ifndef LISTCONSTRAINTSCMD_H
#define LISTCONSTRAINTSCMD_H

#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class ListConstraintsCmd : public Cmd
{
  public:

  Sample *s;

  ListConstraintsCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "list_constraints"; }
  char *help_msg(void) const
  {
    return
    "\n list_constraints\n\n"
    " syntax: list_constraints\n\n"
    "   The list_constraints command prints the list"
    " of active constraints.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( s->ctxt_.oncoutpe() ) s->constraints.list_constraints(cout);
    return 0;
  }
};
#endif
