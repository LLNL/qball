////////////////////////////////////////////////////////////////////////////////
//
// ParOptCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ParOptCmd.h,v 1.1 2007/06/05 21:34:10 draeger1 Exp $

#ifndef PAROPTCMD_H
#define PAROPTCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"

class ParOptCmd : public Cmd {
  private:

  public:

  Sample *s;

  ParOptCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "paropt"; }
  char *help_msg(void) const {
    return 
    "\n paropt\n\n"
    " syntax: paropt filename n [nscf] [nnonscf]\n\n"
    "   The paropt command chooses the parallel parameters which will\n"
    " run n steps of simulation as efficiently as possible.\n\n";
  }

  int action(int argc, char **argv);

};
#endif
