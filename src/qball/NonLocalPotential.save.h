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
#include <math/Matrix.h>

class NonLocalPotential
{
  private:
  
  const Context& ctxt_;
  const AtomSet& atoms_;
  const SlaterDet& sd_;
  const Basis& basis_;

  int nsp;   // number of species
  int nspnl; // number of non-local species
  
  vector<int>             lmax;     // lmax[is]
  vector<int>             lloc;     // lloc[is]
  vector<int>             na;       // na[is]
  //vector<int>             naloc;    // naloc[is]
  //vector<int>             nalocmax; // nalocmax[is]
  vector<int>             npr;      // npr[is]
  vector<int>             nprna;    // nprna[is]
  vector<vector<int> >    lproj;    // lproj[is][ipr]
  vector<vector<double> > wt;       // wt[is][ipr]
  vector<vector<double> > twnl;     // twnl[is][npr*ngwl]
  vector<vector<double> > dtwnl;    // dtwnl[is][6*npr*ngwl], ij=0,..,5
  
  //vector<DoubleMatrix*>   anl;      // anl[is][ipr*ia][ig]
  //vector<vector<double> > singr;    // singr[is][naloc*2*ngwloc]
  //vector<vector<double> > cosgr;    // cosgr[is][naloc*2*ngwloc]
  
  vector<int>             nquad;    // nquad[is]
  vector<vector<double> > rquad;    // rquad[is][iquad], iquad = 0, nquad[is]-1
  vector<vector<double> > wquad;    // wquad[is][iquad], iquad = 0, nquad[is]-1
    
  mutable TimerMap tmap;
  void init(void);
   
  public:
  
  NonLocalPotential(const AtomSet& as, const SlaterDet& sd) :  
    ctxt_(sd.context()), atoms_(as), sd_(sd), basis_(sd.basis()) { init(); }
  ~NonLocalPotential(void);
               
  void update_twnl(void);
  //void update_eigr(vector<vector<double> >& tau);
  //void update_anl(void);
  double energy(bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl);
};
#endif

// Local Variables:
// mode: c++
// End:
