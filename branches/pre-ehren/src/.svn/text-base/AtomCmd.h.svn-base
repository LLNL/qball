////////////////////////////////////////////////////////////////////////////////
//
// AtomCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomCmd.h,v 1.6 2009/12/18 18:40:13 draeger1 Exp $

#ifndef ATOMCMD_H
#define ATOMCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class AtomCmd : public Cmd
{
  public:

  Sample *s;

  AtomCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "atom"; }
  char *help_msg(void) const
  {
    return 
    "\n atom\n\n"
    " syntax: atom name species x y z [vx vy vz]\n\n"
    "   The atom command defines a new atom and adds it to the atom list.\n"
    "   The name can be any character string, the species must be the name\n"
    "   of a file containing the definition of the pseudopotential for this\n"
    "   atomic species. The position of the atom is specified by x y and z.\n"
    "   Optionally, the atom velocity can be specified by vx vy and vz.\n\n";
  }

  int action(int argc, char **argv)
  {
    string name;
    string species;
    D3vector position;
    D3vector velocity;
  
    // atom must be defined with either 3 or 6 arguments
    if ( argc != 6 && argc != 9 )
    {
      if ( ui->oncoutpe() )
        cout << "<!-- use: atom name species x y z [vx vy vz] -->" << endl;
      return 1;
    }
  
    name = argv[1];
    species = argv[2];
    position.x = atof(argv[3]);
    position.y = atof(argv[4]);
    position.z = atof(argv[5]);
    if ( argc == 9 )
    {
      velocity.x = atof(argv[6]);
      velocity.y = atof(argv[7]);
      velocity.z = atof(argv[8]);
    }
  
    Atom *a = new Atom(name,species,position,velocity);
    
    if ( !(s->atoms.addAtom( a ) ) )
    {
      if ( ui->oncoutpe() )
        cout << "<ERROR> AtomCmd: could not add atom " << name << " </ERROR>" << endl;
      delete a;
      return 1;
    }
    
#if DEBUG
cout << " dbg check "<< __FILE__ <<" "<< __LINE__ <<" mype="<< mype << endl;
#endif
 
    s->wf.set_nel(s->atoms.nel());
    if ( s->wfv != 0 )
    {
      s->wfv->set_nel(s->atoms.nel());
      s->wfv->clear();
    }
    
#if DEBUG
cout << " dbg check "<< __FILE__ <<" "<< __LINE__ <<" mype="<< mype << endl;
#endif
    return 0;
  }
};
#endif
