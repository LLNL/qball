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
// HubbardPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef HUBBARDPOTENTIAL_H
#define HUBBARDPOTENTIAL_H

#include "AtomSet.h"
#include "SymmetrySet.h"
#include "Wavefunction.h"
#include "Context.h"
#include "Matrix.h"
#include "D3vector.h"

class HubbardPotential
{
  private:
  
  const Context& ctxt_;
  AtomSet& atoms_;
  const SymmetrySet& symset_;
  const Wavefunction& wf_;
  
  int nsp;   // number of species
  vector<int> na; // number of atoms of each species

  // DFT+U parameters (alpha unused for now, should be zero)
  vector<int> hub_l_;
  vector<double> hub_u_;
  vector<double> hub_alpha_;
  vector<int> lmsize_;
  vector<vector<vector<vector<double> > > > phiylm;
  vector<vector<vector<double> > > ylmsym;
  vector<D3vector> rmrand;
  vector<double> rmrandlen;
  vector<vector<vector<double> > > rmsym;
  vector<vector<double> > rmsymlen;
  
  void init(void);
   
  public:
  
  HubbardPotential(AtomSet& as, const SymmetrySet& symset, const Wavefunction& wf) :  
      wf_(wf), ctxt_(wf.context()), atoms_(as), symset_(symset) { init(); }
  ~HubbardPotential(void);
               
  void update_phiylm(void);
  double energy(bool compute_hpsi, Wavefunction& dwf);

};
#endif
