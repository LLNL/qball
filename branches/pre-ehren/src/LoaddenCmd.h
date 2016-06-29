////////////////////////////////////////////////////////////////////////////////
//
// LoaddenCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LoaddenCmd.h,v 1.1 2009/03/25 22:30:34 draeger1 Exp $

#ifndef LOADDENCMD_H
#define LOADDENCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class LoaddenCmd : public Cmd {
  public:

  Sample *s;

  LoaddenCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "loadden"; }

  char *help_msg(void) const {
    return 
    "\n loadden\n\n"
    " syntax: loadden filename \n\n"
    "   The loadden command loads the real-space charge density from file.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
