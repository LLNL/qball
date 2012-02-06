////////////////////////////////////////////////////////////////////////////////
//
// Atom.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Atom.h,v 1.7 2007/12/18 19:38:42 draeger1 Exp $

#ifndef ATOM_H
#define ATOM_H

#include "D3vector.h"
#include <string>
using namespace std;

class Atom
{
  private:
  
  string name_;
  string species_;
  D3vector position_;
  D3vector velocity_;
  bool locked_; // if locked_ is true, keep position and velocity fixed at initial values
  bool rescalewhenlocked_; // if atom is locked, whether or not to rescale when cell moves

  public:

  Atom (string name, string species, D3vector position, D3vector velocity);
  string name(void) { return name_; };
  string species(void) { return species_; };
  D3vector position(void) { return position_; };
  D3vector velocity(void) { return velocity_; };
  void set_position(D3vector p);
  void set_position(D3vector p, bool rescale);
  void set_velocity(D3vector v);
  void block(void) { velocity_ = D3vector(0.0,0.0,0.0); };
  void printsys(ostream &os, string atomcmd) const;
  bool islocked(void) { return locked_; };
  void lock_atom(void) { locked_ = true; };
  void unlock_atom(void) { locked_ = false; };
  bool rescalewhenlocked(void) { return rescalewhenlocked_; };
  void set_rescale(bool rescale) { rescalewhenlocked_ = rescale; };
};

ostream& operator << ( ostream &os, Atom &a );
#endif
