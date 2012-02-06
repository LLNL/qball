////////////////////////////////////////////////////////////////////////////////
//
// SaveESPCmd.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveESPCmd.h,v 1.1 2009/07/21 20:31:47 draeger1 Exp $

#ifndef SAVEESPCMD_H
#define SAVEESPCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"

class SaveESPCmd : public Cmd
{
  private:

  public:
  Sample *s;

  SaveESPCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "saveesp"; }

  char *help_msg(void) const
  {
    return
    "\n saveesp\n\n"
    " syntax: saveesp\n\n"
        "   The saveesp command prints the electrostatic potential on a grid \n\n";
  }

  int action(int argc, char **argv);

  SaveESPCmd();
  ~SaveESPCmd();
};
#endif
