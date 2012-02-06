////////////////////////////////////////////////////////////////////////////////
//
// Atom.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Atom.C,v 1.8 2007/12/18 19:38:42 draeger1 Exp $

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
