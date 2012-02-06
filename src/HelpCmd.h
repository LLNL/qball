////////////////////////////////////////////////////////////////////////////////
//
// HelpCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: HelpCmd.h,v 1.2 2008/04/07 22:00:37 draeger1 Exp $

#ifndef HELPCMD_H
#define HELPCMD_H

#include <iostream>
#include <iomanip>
#include <string.h>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class HelpCmd : public Cmd
{
  public:

  Sample *s;

  HelpCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "help"; }

  char *help_msg(void) const
  {
    return 
    "\n help\n\n"
    " syntax: help [command_name]\n\n"
    "   The help command gives a short description of a command. If used\n"
    "   without arguments, help prints a list of valid commands.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( ui->oncoutpe() )
    {
      if ( argc == 1 )  // no arguments
      {
        cout << endl << " valid commands are: " << endl << endl;
        list<Cmd*>::iterator cmd = ui->cmdlist.begin();
        int n = 0;
        while ( cmd != ui->cmdlist.end() )
        {
          n++;
          cout.setf(ios::left,ios::adjustfield);
          cout << " " << setw(16) << (*cmd)->name();
          cout.setf(ios::right,ios::adjustfield);
          if ( n%4 == 0 ) cout << endl;
          cmd++;
        }
        if ( n%4 != 0 ) cout << endl;
        cout << endl;
      }
      else if ( argc == 2 ) // one argument
      {
        // search command list
        Cmd *cmdptr = ui->findCmd(argv[1]);

        if ( cmdptr )
        {
          cout << cmdptr->help_msg();
        }
        else
        {
          cout << " help: " << argv[1] << " is not a valid command" << endl;
        }
      }
      else
      {
        cout << " use: help [command_name]" << endl;
      }
    }
    return 0;
  }
};
#endif
