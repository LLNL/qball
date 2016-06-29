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
// Atom.C:
//
////////////////////////////////////////////////////////////////////////////////

#include "Atom.h"
#include <iomanip>
using namespace std;

Atom::Atom (string newname, string newspecies, D3vector pos, D3vector vel)
{
  name_ = newname;
  species_ = newspecies;
  position_ = pos;
  velocity_ = vel;
  locked_ = false;
  rescalewhenlocked_ = true;
}

void Atom::set_position(D3vector p) { 
  if (!locked_) 
    position_ = p; 
}

void Atom::set_position(D3vector p, bool rescale) { 
  if (!locked_ || (locked_ && rescale && rescalewhenlocked_)) 
    position_ = p; 
}

void Atom::set_velocity(D3vector v) { 
  if (!locked_) 
    velocity_ = v; 
}

void Atom::printsys(ostream& os, string atomcmd) const {
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(8);
  os << atomcmd << " " << name_ << " " << species_ << " " 
     << setw(14) << position_.x << " " << setw(14) << position_.y << " " << setw(14) << position_.z << " " 
     << setw(14) << velocity_.x << " " << setw(14) << velocity_.y << " " << setw(14) << velocity_.z << " " 
     << endl;
  if (locked_)
     os << "lock " << name_ << endl;
}
  
ostream& operator << ( ostream &os, Atom &a )
{
  os.setf(ios::left,ios::adjustfield);
  os << "  <atom name=\"" << a.name() << "\""
     << " species=\"" << a.species() << "\">\n"
     << "    <position> ";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setw(12) << setprecision(8) << a.position().x << " "
     << setw(12) << setprecision(8) << a.position().y << " "
     << setw(12) << setprecision(8) << a.position().z << "  "
     << " </position>\n"
     << "    <velocity> ";
  os.setf(ios::scientific,ios::floatfield);
  os << setw(13) << setprecision(6) << a.velocity().x << " "
     << setw(13) << setprecision(6) << a.velocity().y << " "
     << setw(13) << setprecision(6) << a.velocity().z
     << " </velocity>\n  </atom>" << endl;
  return os;
}
