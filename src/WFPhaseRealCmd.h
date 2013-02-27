////////////////////////////////////////////////////////////////////////////////
//
// WFPhaseRealCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WFPhaseRealCmd.h,v 1.0 2011-03-25 15:56:17 schleife Exp $

#ifndef WF_PHASE_REAL_CMD_H
#define WF_PHASE_REAL_CMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

// AS: change phase of the wave function to make it real for Gamma only

class WFPhaseRealCmd : public Cmd
{
  public:

  Sample *s;

  WFPhaseRealCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "wf_phase_real"; }
  char *help_msg(void) const
  {
    return
    "\n wf_phase_real\n\n"
    " syntax: wf_phase_real\n\n"
    "   The wf_phase_real command modifies the phase of the wave function so that it is real for Gamma only.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 1 )
    {
      if ( ui->oncoutpe() )
      {
        cout << " use: wf_phase_real" << endl;
      }
      return 1;
    }

    if ( ui->oncoutpe() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout << setprecision(3)
           << " modifying the phase of the wave functions " << endl;
    }

    s->wf.phase_wf_real();

    return 0;
  }
};
#endif
