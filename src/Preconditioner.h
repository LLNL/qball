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
// Preconditioner.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

class Sample;
class EnergyFunctional;
class Wavefunction;

#include <vector>
#include <valarray>
using namespace std;

class Preconditioner
{
  private:
  
  const Sample& s_;
  const EnergyFunctional& ef_;
  const Wavefunction& wf_;
  vector<vector<valarray<double> > > diag_; // diag_[ispin][ikp][ig]

  public:

  void update(void);
  
  const valarray<double>& diag(int ispin, int ikp) const
  { return diag_[ispin][ikp]; }

  Preconditioner(const Sample& s, const Wavefunction& wf, const EnergyFunctional& ef);
  //~Preconditioner();
};
#endif

// Local Variables:
// mode: c++
// End:
