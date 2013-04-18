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
// MMSpecies.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MMSPECIES_H
#define MMSPECIES_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;
#include "Context.h"

class MMSpecies {
  private:
  
  const Context& ctxt_;
  string name_;         // name used to refer to species in current application
  double mass_;        // mass in a.m.u (Carbon = 12.0)

  public:

  MMSpecies(const Context& ctxt, string name, double mass);
  
  const Context& context(void) const { return ctxt_; }
  const string& name(void) const { return name_; }
  double mass(void) const { return mass_; }

  void info(ostream& os);
  void printsys(ostream& os) const;
    
};
ostream& operator << ( ostream &os, MMSpecies &a );
#endif
