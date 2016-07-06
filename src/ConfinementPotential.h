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
// ConfinementPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef CONFINEMENTPOTENTIAL_H
#define CONFINEMENTPOTENTIAL_H

class Basis;
#include <valarray>
using namespace std;

class ConfinementPotential
{
  private:
  
  double ecuts_, facs_, sigmas_;
  const Basis& basis_;
  valarray<double> fstress_, dfstress_;

  public:
  
  const valarray<double>& fstress(void) const { return fstress_; }
  const valarray<double>& dfstress(void) const { return dfstress_; }
  
  void update(void);
  
  const Basis& basis() const { return basis_; }

  ConfinementPotential(double ecuts, double facs, double sigmas, 
    const Basis& basis);
  ~ConfinementPotential();
};
#endif
