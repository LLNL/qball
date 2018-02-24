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
//  CoordinateConstraint.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "CoordinateConstraint.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void CoordinateConstraint::setup(const AtomSet& atoms)
{
  // find coordinate in tau array corresponding to atom atom_name
  is1_ = atoms.is(atom_name_);
  ia1_ = atoms.ia(atom_name_);
  assert(is1_>=0);
  assert(ia1_>=0);
}

////////////////////////////////////////////////////////////////////////////////
void CoordinateConstraint::update(double dt)
{
  // nothing to update
}

////////////////////////////////////////////////////////////////////////////////
bool CoordinateConstraint::enforce_r(const vector<vector<double> > &r0,
vector<vector<double> > &rp) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  double* pr1p  = &rp[is1_][3*ia1_];
  D3vector r1p(pr1p);

  double sigma = fabs(r1p[coord_] - r1[coord_]);
  
  if ( sigma < tol_ ) return true;

  pr1p[coord_] = pr1[coord_];

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool CoordinateConstraint::enforce_v(const vector<vector<double> > &r0,
vector<vector<double> > &v0) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  double* pv1 = &v0[is1_][3*ia1_];
  D3vector v1(pv1);

  const double err = fabs(v1[coord_]);
  if ( err < tol_ ) return true;

  pv1[coord_] = 0.0;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void CoordinateConstraint::compute_force(const vector<vector<double> > &r0, const vector<vector<double> > &f) {
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  const double* pf1 = &f[is1_][3*ia1_];
  D3vector f1(pf1);

  force_ = fabs(f1[coord_]);
}

////////////////////////////////////////////////////////////////////////////////
ostream& CoordinateConstraint::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << " <constraint name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"" << atom_name_ << "\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "  value=\"" << setprecision(6) << 0;
  os << "\" velocity=\"" << setprecision(6) << 0 << "\"\n";
  os << "  force=\"" << setprecision(6) << force_;
  os << "\" weight=\"" << setprecision(6) << weight_ << "\"/>";
  return os;
}
