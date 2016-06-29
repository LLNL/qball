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
// SymmetrySet.C
//
////////////////////////////////////////////////////////////////////////////////

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
