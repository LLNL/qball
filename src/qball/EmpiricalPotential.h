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
// EmpiricalPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef EMPIRICALPOTENTIAL_H
#define EMPIRICALPOTENTIAL_H

#include "Context.h"
#include <math/D3vector.h>
#include <string>
#include <vector>
using namespace std;

class EmpiricalPotential {
  private:

  const Context& ctxt_;
  string pottype_;
  string spname1_, spname2_;
  string filename_;
  int npts_;
  vector<double> r_;
  vector<double> pot_;
  vector<double> pot_spl_;
  const double param1_;
  const double param2_;
  const double param3_;

  public:

  int is1,is2;
  EmpiricalPotential(const Context& ctxt, string pottype, string spname1, string spname2, double param1, double param2);
  EmpiricalPotential(const Context& ctxt, string pottype, string spname1, string spname2, double param1, double param2, double param3);
  EmpiricalPotential(const Context& ctxt, string spname1, string spname2, string filename);
  ~EmpiricalPotential();

  double r(int i);
  double pot(double rval);
  D3vector force(D3vector r12);  // returns force on r1 (f2 = -f1)
  string spname1(void) { return spname1_; }
  string spname2(void) { return spname2_; }
  double param1(void) { return param1_; }
  double param2(void) { return param2_; }
  double param3(void) { return param3_; }
  string pottype(void) { return pottype_; }
  string filename(void) { return filename_; }
  int npts(void) { return npts_; }
  void printsys(ostream& os) const;
  
};

#endif

// Local Variables:
// mode: c++
// End:
