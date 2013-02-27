////////////////////////////////////////////////////////////////////////////////
//
// MDSaveCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDSaveCmd.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef MDSAVECMD_H
#define MDSAVECMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class MDSaveCmd : public Cmd
{
  public:

  Sample *s;

  MDSaveCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "mdsave"; }

  char *help_msg(void) const
  {
    return 
    "\n save\n\n"
    " syntax: mdsave [-dump/-states (optional)] [filebase (optional)] \n\n"
    "   The mdsave command checkpoints the current md iteration to a directory md.ITER#.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
