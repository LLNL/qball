////////////////////////////////////////////////////////////////////////////////
//
// ConstraintCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintCmd.h,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#ifndef CONSTRAINTCMD_H
#define CONSTRAINTCMD_H

#include "UserInterface.h"
#include "Sample.h"

class ConstraintCmd : public Cmd
{
  public:

  Sample *s;

  ConstraintCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "constraint"; }
  char *help_msg(void) const
  {
    return
    "\n constraint\n\n"
    " syntax:\n\n"
    "   constraint define position name atom\n"
    "   constraint define distance name atom1 atom2 distance [velocity]\n"
    "   constraint define angle name atom1 atom2 atom3 angle [velocity]\n"
    "   constraint define torsion name atom1 atom2 atom3 atom4 angle [velocity]\n"
    "   constraint set name value [velocity]\n"
    "   constraint delete name\n"
    "   constraint list\n"
    "   constraint enforce\n\n"
    "   Constraints are enforced at each MD step if ions are allowed to move.\n"
    "   Velocity parameters are optional.\n\n";
  }

  int action(int argc, char **argv)
  {
    const bool oncoutpe = s->ctxt_.oncoutpe();
    if ( argc < 2 )
    {
      if ( oncoutpe )
        cout << help_msg();
      return 1;
    }
    string subcmd(argv[1]);
    if ( subcmd == "define" )
    {
      return s->constraints.define_constraint(s->atoms,argc,argv);
    }
    else if ( subcmd == "set" )
    {
      return s->constraints.set_constraint(argc,argv);
    }
    else if ( subcmd == "delete" )
    {
      return s->constraints.delete_constraint(argc,argv);
    }
    else if ( subcmd == "enforce" )
    {
      s->constraints.enforce(s->atoms);
      // reset velocities to zero to avoid jump in temperature
      s->atoms.reset_velocities();
      // Note: should catch exception here if constraints could not be enforced
    }
    else if ( subcmd == "list" )
    {
      if ( oncoutpe )
        s->constraints.list_constraints(cout);
    }
    else
    {
      if ( oncoutpe )
        cout << help_msg();
    }

    return 0;
  }
};
#endif
