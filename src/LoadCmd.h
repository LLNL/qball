////////////////////////////////////////////////////////////////////////////////
//
// LoadCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LoadCmd.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef LOADCMD_H
#define LOADCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class LoadCmd : public Cmd
{
  public:

  Sample *s;

  LoadCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "load"; }

  char *help_msg(void) const
  {
    return 
    "\n load\n\n"
    " syntax: load [-serial] filename \n\n"
    "   The load command loads a sample from the file filename.\n"
    "   The -serial option bypasses the parallel load algorithm.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
