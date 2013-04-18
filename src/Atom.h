////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// Atom.h:
//
////////////////////////////////////////////////////////////////////////////////

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
