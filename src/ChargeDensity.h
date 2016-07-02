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
// ChargeDensity.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef CHARGEDENSITY_H
#define CHARGEDENSITY_H

#include <vector>
#include <valarray>
#include <complex>
#include <string>
#include <map>
#include "Timer.h"
#include "Context.h"

class Sample;
class Wavefunction;
class FourierTransform;
class Basis;
class AtomSet;

typedef map<string,Timer> TimerMap;

class ChargeDensity {
  private:
  
  const Context& ctxt_;
  Context vcontext_;
  Wavefunction& wf_;
  AtomSet& atoms_;
  Basis* vbasis_;
  FourierTransform* vft_;
  int np0v_, np1v_, np2v_;
  vector<vector<FourierTransform*> > ft_; // ft_[ispin][ikp];
  valarray<complex<double> > rhotmp;
  vector<int> symindexloc;
  vector<int> symmultloc;
  int nsym_;
  int nsymgrp_;
  bool ultrasoft_;
  bool highmem_;
  bool nlcc_;
  bool tddft_involved_;
  vector<vector<complex<double> > > qnmg_;
  vector<vector<complex<double> > > sfactloc_;
  vector<double> rhornlcc_; 

  void initialize(const Sample& s);
  void initializeSymmetries(const Sample& s);
  
  public:
  
  mutable TimerMap tmap;

  vector<vector<double> > rhor; // rhor[ispin][i]
  vector<vector<complex<double> > > rhog; // rhog[ispin][ig]
  vector<vector<double> > xcrhor; 
  vector<complex<double> > rhognlcc; 
  vector<vector<complex<double> > > xcrhog; 

  void update_density();
  void update_rhor(void);
  
  Basis* vbasis(void) const { return vbasis_; }
  const Context& vcontext(void) const { return vcontext_; }
  FourierTransform* vft(void) const { return vft_; }
  FourierTransform* ft(int ispin, int ikp) const { return ft_[ispin][ikp]; }
  void reshape_rhor(const Context& oldvctxt,const Context& newvctxt);
  void update_usfns();
  void set_nlcc(bool nlcc) { nlcc_ = nlcc; }
  bool nlcc() { return nlcc_; }
  void update_nlcc();
  void add_nlccden();
  void nlcc_forceden(int is, vector<complex<double> > &rhog);
  void print_memory(ostream&os, double& totsum, double& locsum) const;
  void print_timing();
  
  ChargeDensity(const Sample& s);
  ChargeDensity(const Sample& s, Wavefunction& cdwf);
  ~ChargeDensity();
};
#endif
