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
// SymmetrySet.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SYMMETRYSET_H
#define SYMMETRYSET_H

#include "Context.h"
#include "Symmetry.h"
#include <string>
using namespace std;

class SymmetrySet {
  private:
  
  const Context& ctxt_;
  int nsym_;
  
  public:
  
  SymmetrySet(const Context& ctxt) : ctxt_(ctxt) {
    nsym_ = 0;
  }

  vector<Symmetry *> symlist; // symlist[nsym], symmetry operators

  const Context& context(void) const { return ctxt_; }
  bool addSymmetry(Symmetry *sym);
  bool delSymmetry(int i);
  bool reset(void); // remove all symmetries
  int nsym(void) const { return nsym_; };
  int size(void);
};
ostream& operator << ( ostream &os, SymmetrySet &ss );
#endif
