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
//  CoordinateConstraint.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef COORDINATECONSTRAINT_H
#define COORDINATECONSTRAINT_H

#include "Constraint.h"
#include <cassert>
#include <cmath> // fabs

class AtomSet;

class CoordinateConstraint : public Constraint
{
  std::string atom_name_;
  int    ia1_, is1_;
  double force_, weight_, tol_;
  int coord_;

  public:

  CoordinateConstraint(std::string name, std::string atom_name, const int coord, double tolerance):
    atom_name_(atom_name), tol_(tolerance), coord_(coord) {
    name_ = name;
    names_.resize(1);
    names_[0] = atom_name_;
    force_ = 0.0;
    weight_ = 1.0;
  }
  ~CoordinateConstraint(void) {}

  std::string type(void) const { return "coordinate"; }
  double value(void) const { return 0.0; }
  double velocity(void) const { return 0.0; }
  double force(void) const { return force_; }
  double weight(void) const { return weight_; }
  double tolerance(void) const { return tol_; }
  void set_value(double value)
  {
    // no value to set
  }
  void set_velocity(double velocity)
  {
    // no value to set
  }

  void setup(const AtomSet& atoms);
  int dofs(void) const { return 3; }
  void update(double dt);
  bool enforce_r(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &rp) const;
  bool enforce_v(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &v0) const;
  void compute_force(const std::vector<std::vector<double> > &r0,
                     const std::vector<std::vector<double> > &f);
  std::ostream& print( std::ostream& os );

};
#endif

// Local Variables:
// mode: c++
// End:
