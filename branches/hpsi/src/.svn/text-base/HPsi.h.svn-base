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
// HPsi.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef HPSI_H
#define HPSI_H

#include "Wavefunction.h"
#include "AtomSet.h"
#include "FourierTransform.h"
#include "ChargeDensity.h"
#include "Matrix.h"
#include <vector>

class HPsi
{
  private:
  
   const AtomSet& atoms_;
  
   int nsp;          // number of species
   int nspnl;        // number of non-local species
  
   vector<int>             lmax;          // lmax[is]
   vector<int>             lloc;          // lloc[is]
   vector<int>             na;            // na[is]
   vector<int>             npr;           // npr[is]
   vector<int>             nprna;         // nprna[is]
   vector<vector<int> >    lproj;         // lproj[is][ipr]
   vector<vector<double> > wt;            // wt[is][ipr]
   vector<vector<vector<double> > > twnl; // twnl[is][npr*ngwl]
  
   vector<int>             nquad;    // nquad[is]
   vector<vector<double> > rquad;    // rquad[is][iquad], iquad = 0, nquad[is]-1
   vector<vector<double> > wquad;    // wquad[is][iquad], iquad = 0, nquad[is]-1
   vector<vector<double> >& vofr_;
   vector<vector<FourierTransform*> > ft;
   
   void init(const Wavefunction& wfi);
   
  public:
  
   HPsi(const Wavefunction& wfi, const AtomSet& as, const ChargeDensity& cd, vector<vector<double> >& v_r);
   ~HPsi(void);
               
   void cell_moved(const Wavefunction& wfi);
   void compute(const Wavefunction& wf, Wavefunction& dwf);
};
#endif
