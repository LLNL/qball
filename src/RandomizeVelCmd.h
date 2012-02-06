////////////////////////////////////////////////////////////////////////////////
//
// RandomizeVelCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RandomizeVelCmd.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef RANDOMIZEVELCMD_H
#define RANDOMIZEVELCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeVelCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeVelCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "randomize_vel"; }
  char *help_msg(void) const
  {
    return 
    "\n randomize_vel\n\n"
    " syntax: randomize_vel [(optional) temperature]\n\n"
    "   The randomize_vel command randomly sets the ion velocities\n"
    "   to achieve a given temperature, either given as an argument\n"
    "   or with the th_temp variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc > 2 )
    {
      if ( ui->oncoutpe() )
      {
        cout << "<!-- use: randomize_vel [(optional) temperature] -->" << endl;
      }
      return 1;
    }
    double temp = s->ctrl.th_temp;
    if ( argc == 2 )
      temp = atof(argv[1]);

    s->atoms.randomize_velocities(temp);
    
    return 0;
  }
};
#endif
