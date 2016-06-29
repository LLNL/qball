////////////////////////////////////////////////////////////////////////////////
//
// ResetVcmCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ResetVcmCmd.h,v 1.2 2009/03/27 00:53:24 draeger1 Exp $

#ifndef RESETVCMCMD_H
#define RESETVCMCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class ResetVcmCmd : public Cmd
{
  public:

  Sample *s;

  ResetVcmCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "reset_vcm"; }
  char *help_msg(void) const
  {
    return
    "\n reset_vcm\n\n"
    " syntax: reset_vcm \n\n"
    "   The reset_vcm command subtracts the velocity of the center\n"
    "   of mass from the velocity of each atom.\n\n";
  }

  int action(int argc, char **argv)
  {
    s->atoms.reset_vcm();
    return 0;
  }
};
#endif
