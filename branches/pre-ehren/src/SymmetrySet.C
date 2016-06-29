////////////////////////////////////////////////////////////////////////////////
//
// SymmetrySet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SymmetrySet.C,v 1.3 2008/04/07 22:00:37 draeger1 Exp $

#include "SymmetrySet.h"
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
bool SymmetrySet::addSymmetry(Symmetry *sym) {
  symlist.push_back(sym);
  nsym_++;
  assert(symlist.size() == nsym_);
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool SymmetrySet::delSymmetry(int i) {
  delete symlist[i];
  symlist.erase(symlist.begin() + i);
  nsym_--;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool SymmetrySet::reset(void) {
  // delete all symmetry ops from symlist
  for ( int i = 0; i < symlist.size(); i++ )
    delete symlist[i];
  symlist.resize(0);
  nsym_ = 0;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator << ( ostream &os, SymmetrySet &ss ) {
  if ( ss.context().oncoutpe() ) {
    os << "<symmetryset>\n";
    for ( int i = 0; i < ss.symlist.size(); i++ )
      os << *ss.symlist[i];
    os << "</symmetryset>\n";
  }
  return os;
}
