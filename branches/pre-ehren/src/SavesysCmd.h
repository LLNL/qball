////////////////////////////////////////////////////////////////////////////////
//
// SavesysCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SavesysCmd.h,v 1.1 2006/07/14 18:45:08 draeger1 Exp $

#ifndef SAVESYSCMD_H
#define SAVESYSCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SavesysCmd : public Cmd
{
  public:

  Sample *s;

  SavesysCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "savesys"; }

  char *help_msg(void) const
  {
    return 
    "\n savesys\n\n"
    " syntax: savesys filename \n\n"
    "   The savesys command saves the coordinates in sys format to the file filename.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
