////////////////////////////////////////////////////////////////////////////////
//
// SaveCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveCmd.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef SAVECMD_H
#define SAVECMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SaveCmd : public Cmd
{
  public:

  Sample *s;

  SaveCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "save"; }

  char *help_msg(void) const
  {
    return 
    "\n save\n\n"
    " syntax: save filename \n\n"
    "   The save command saves the sample to the file filename.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
