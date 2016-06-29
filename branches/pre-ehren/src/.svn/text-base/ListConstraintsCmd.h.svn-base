////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// ListConstraintsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ListConstraintsCmd.h,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#ifndef LISTCONSTRAINTSCMD_H
#define LISTCONSTRAINTSCMD_H

#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class ListConstraintsCmd : public Cmd
{
  public:

  Sample *s;

  ListConstraintsCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "list_constraints"; }
  char *help_msg(void) const
  {
    return
    "\n list_constraints\n\n"
    " syntax: list_constraints\n\n"
    "   The list_constraints command prints the list"
    " of active constraints.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( s->ctxt_.onpe0() ) s->constraints.list_constraints(cout);
    return 0;
  }
};
#endif
