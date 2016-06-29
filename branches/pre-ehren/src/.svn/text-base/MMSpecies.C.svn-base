////////////////////////////////////////////////////////////////////////////////
//
// MMSpecies.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MMSpecies.C,v 1.1 2007/04/10 23:30:12 draeger1 Exp $

#include "MMSpecies.h"
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
MMSpecies::MMSpecies(const Context& ctxt, string name, double mass) : 
  ctxt_(ctxt), name_(name), mass_(mass) {

}
  
ostream& operator << ( ostream &os, MMSpecies &s ) {
  os <<"<mmspecies name=\"" << s.name() 
     << "\" mass=\"" << s.mass() << "\"/>" << endl;
  return os;
}

void MMSpecies::printsys(ostream& os) const {
  os.setf(ios::left,ios::adjustfield);
  os << "mmspecies " << name() << " " << mass() << endl;
}

void MMSpecies::info(ostream &os) {
  os.setf(ios::left,ios::adjustfield);
  os << " name_ = " << name() << endl;
  os << " mass_ = " << mass() << " (amu)" << endl;
}
