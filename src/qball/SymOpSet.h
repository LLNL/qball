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
// SymOpSet.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef SYMOPSET_H
#define SYMOPSET_H

#include <math/D3vector.h>
#include "SymOp.h"
#include <qball/UnitCell.h>
#include <vector>
using namespace std;

class SymOpSet {
  private:
  
  int nsym_;
  vector<SymOp*> symset_;

  public:

  SymOpSet();
  ~SymOpSet();

  void generateOps(char const* name);
  void convertOpsToXtal(const UnitCell& uc);
  void printXtal(ostream& os);  
  vector<SymOp*> returnAllOps(void) const {return symset_; }
  SymOp* returnOp(int n) const {return symset_[n]; }
  int nsym(void) { return nsym_; }
  
  void addOp(SymOp* addsym);
  void removeOp(SymOp* delsym);
  void clear(void);
};

#endif

// Local Variables:
// mode: c++
// End:
