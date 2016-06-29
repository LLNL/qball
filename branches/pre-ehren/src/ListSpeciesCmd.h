////////////////////////////////////////////////////////////////////////////////
//
// ListSpeciesCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ListSpeciesCmd.h,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#ifndef LISTSPECIESCMD_H
#define LISTSPECIESCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class ListSpeciesCmd : public Cmd
{
  public:

  Sample *s;

  ListSpeciesCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "list_species"; }
  char *help_msg(void) const
  {
    return 
    "\n list_species\n\n"
    " syntax: list_species\n\n"
    "   The list_species command prints a list of all defined species.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 1 )
    {
      if ( ui->oncoutpe() )
      {
        cout << "<!-- use: list_species -->" << endl;
      }
      return 1;
    }
    s->atoms.listSpecies();
    return 0;
  }
};
#endif
