////////////////////////////////////////////////////////////////////////////////
//
// FoldInWsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FoldInWsCmd.h,v 1.2 2009/03/27 00:53:24 draeger1 Exp $

#ifndef FOLDINWSCMD_H
#define FOLDINWSCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class FoldInWsCmd : public Cmd
{
  public:

  Sample *s;

  FoldInWsCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "fold_in_ws"; }
  char *help_msg(void) const
  {
    return
    "\n fold_in_ws\n\n"
    " syntax: fold_in_ws \n\n"
    "   The fold_in_ws command folds all atomic positions back in\n"
    "   the Wigner-Seitz cell of the current unit cell.\n\n";
  }

  int action(int argc, char **argv)
  {
    s->atoms.fold_in_ws();
    return 0;
  }
};
#endif
