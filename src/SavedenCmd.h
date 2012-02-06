////////////////////////////////////////////////////////////////////////////////
//
// SavedenCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SavedenCmd.h,v 1.1 2006/12/22 01:17:11 draeger1 Exp $

#ifndef SAVEDENCMD_H
#define SAVEDENCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SavedenCmd : public Cmd {
  public:

  Sample *s;

  SavedenCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "saveden"; }

  char *help_msg(void) const {
    return 
    "\n saveden\n\n"
    " syntax: saveden filename \n\n"
    "   The saveden command saves the real-space charge density to file.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
