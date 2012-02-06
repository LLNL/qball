////////////////////////////////////////////////////////////////////////////////
//
// SpeciesCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesCmd.h,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $

#ifndef SPECIESCMD_H
#define SPECIESCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SpeciesCmd : public Cmd
{
  public:

  Sample *s;

  SpeciesCmd(Sample *sample) : s(sample) { s->ctrl.ultrasoft = false; s->ctrl.nlcc = false; };

  char *name(void) const { return "species"; }

  char *help_msg(void) const
  {
    return 
    "\n species\n\n"
    " syntax: species name uri\n\n"
    "   The species command defines a species name.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
