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
// XCPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef XCPOTENTIAL_H
#define XCPOTENTIAL_H

#include "ChargeDensity.h"
#include "LDAFunctional.h"
#include "PBEFunctional.h"
#include "PBESolFunctional.h"
#include "PBERevFunctional.h"
#include "BLYPFunctional.h"
#include <vector>
#include <valarray>
#include <complex>
using namespace std;

class Basis;
class FourierTransform;

class XCPotential
{
  private:
  
  const Context& ctxt_;  
  ChargeDensity& cd_;
  ChargeDensity& cd_ecalc_;
  XCFunctional* xcf_;
  
  vector<vector<double> > vxctmp;          // vxctmp[ispin][ir]
  vector<complex<double> > tmpr;           // tmpr[ir]
  vector<complex<double> > tmp1, tmp2;     // tmp1[ig], tmp2[ig]
  
  double exc_, dxc_, dxc0_, dxc1_, dxc2_;
  int nspin_;
  int ngloc_;
  int np012loc_;
  bool tddft_involved_;
  
  FourierTransform& vft_;
  Basis& vbasis_;

  void initialize(string functional_name);
  
  public:

  const XCFunctional* xcf() { return xcf_; }
  XCPotential(ChargeDensity& cd, const string functional_name);
  XCPotential(ChargeDensity& cd, const string functional_name, ChargeDensity& cd_ecalc);
  ~XCPotential();
  void update(vector<vector<double> >& vr);
  void update_exc(vector<vector<double> >& vr);
  void compute_stress(valarray<double>& sigma_exc);
  double exc(void) { return exc_; }
};

class XCPotentialException
{
  public:
  string msg;
  XCPotentialException(string s) : msg(s) {}
};
#endif

// Local Variables:
// mode: c++
// End:
