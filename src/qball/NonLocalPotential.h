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
// NonLocalPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef NONLOCALPOTENTIAL_H
#define NONLOCALPOTENTIAL_H

#include "AtomSet.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "Context.h"
#include "VectorPotential.h"
#include <math/matrix.h>

class StructureFactor;

class NonLocalPotential
{
  private:
  
  const Context& ctxt_;
  AtomSet& atoms_;
  const Basis& basis_;
  Basis* cdbasis_;  // pointer to charge density basis, needed for ultrasoft
  
  int nsp;          // number of species
  int nspnl;        // number of non-local species
  bool ultrasoft_;  // there are ultrasoft potentials
  bool highmem_; 
  
  vector<int>             lmax;     // lmax[is]
  vector<int>             lloc;     // lloc[is]
  vector<int>             na;       // na[is]
  vector<int>             npr;      // npr[is]
  vector<int>             nprna;    // nprna[is]
  vector<vector<int> >    lproj;    // lproj[is][ipr]
  vector<vector<int> >    icproj;   // icproj[is][ipr]
  vector<vector<vector<vector<int> > > > iprojlm; // iprojlm[is][ic][l][m]
  vector<vector<double> > wt;       // wt[is][ipr]
  vector<vector<double> > twnl;     // twnl[is][npr*ngwl]
  vector<vector<double> > dtwnl;    // dtwnl[is][6*npr*ngwl], ij=0,..,5
  
  vector<int>             nquad;    // nquad[is]
  vector<vector<double> > rquad;    // rquad[is][iquad], iquad = 0, nquad[is]-1
  vector<vector<double> > wquad;    // wquad[is][iquad], iquad = 0, nquad[is]-1

  vector<vector<complex<double> > > qnmg_;   // ultrasoft function Q_nm^I(G), wf basis
  vector<vector<complex<double> > > sfactcd_;  // structure factor of local atoms, cd basis
  vector<vector<complex<double> > > sfactwf_;  // structure factor of local atoms, wf basis
  
  
  mutable TimerMap tmap;
  void init(const bool compute_stress);

  VectorPotential * vp;
  
  public:
  
  NonLocalPotential(AtomSet& as, const Context& ctxt, const Basis & basis, VectorPotential * vparg, const bool compute_stress) :  
    ctxt_(ctxt), atoms_(as), basis_(basis), vp(vparg) { init(compute_stress); }
  ~NonLocalPotential(void);
               
  void update_twnl(const bool compute_stress);
  void update_usfns(SlaterDet& sd, Basis* cdbasis);  // update Q_nm^I(G), beta^I(G) when atoms
                                      // move or basis changes
  void use_highmem(void) { highmem_ = true; }  // use extra memory to speed calculation
  double energy(SlaterDet& sd, bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl,
    vector<complex<double> >& veff);

  void print_memory(ostream&os, int kmult, int kmultloc, double& totsum, double& locsum) const;
  void print_timing();
};
#endif

// Local Variables:
// mode: c++
// End:
