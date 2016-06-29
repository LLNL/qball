////////////////////////////////////////////////////////////////////////////////
//
// RandomizeWfCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RandomizeWfCmd.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef RANDOMIZEWFCMD_H
#define RANDOMIZEWFCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeWfCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeWfCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "randomize_wf"; }
  char *help_msg(void) const
  {
    return 
    "\n randomize_wf\n\n"
    " syntax: randomize_wf [amplitude]\n\n"
    "   The randomize_wf command adds random amplitudes to\n"
    "   the wavefunction Fourier coefficients\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc > 2 )
    {
      if ( ui->oncoutpe() )
      {
        cout << "<!-- use: randomize_wf [amplitude] -->" << endl;
      }
      return 1;
    }
    double amp = 0.02;
    if ( argc == 2 )
      amp = atof(argv[1]);

    // use extra memory for SlaterDets if memory variable = normal, large or huge
    bool highmem = false;
    if (s->ctrl.extra_memory >= 3)
      highmem = true;
    if (s->ctrl.ultrasoft)
      s->wf.randomize_us(amp,s->atoms,highmem);
    else
      s->wf.randomize(amp,highmem);
    
    return 0;
  }
};
#endif
