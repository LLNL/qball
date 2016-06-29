////////////////////////////////////////////////////////////////////////////////
//
// StatusCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StatusCmd.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef STATUSCMD_H
#define STATUSCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"

class StatusCmd : public Cmd
{
  private:

  int niter, nfi;

  public:

  Sample *s;

  StatusCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "status"; }
  char *help_msg(void) const
  {
    return 
    "\n run\n\n"
    " syntax: status \n\n"
    "   The status command print information about the current\n"
    "   status of the simulation\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( ui->oncoutpe() )
    {
      s->wf.info(cout,"wf");
      if ( s->wfv != 0 )
        s->wfv->info(cout,"wfv");
    }
    return 0;
  }

};
#endif
