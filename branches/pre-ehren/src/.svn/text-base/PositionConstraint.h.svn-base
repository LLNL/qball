////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
//  PositionConstraint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PositionConstraint.h,v 1.1 2010/01/16 01:26:35 draeger1 Exp $

#ifndef POSITIONCONSTRAINT_H
#define POSITIONCONSTRAINT_H

#include "Constraint.h"
#include <cassert>
#include <cmath> // fabs

class AtomSet;

class PositionConstraint : public Constraint
{
  std::string name1_;
  int    ia1_, is1_;
  double force_, weight_, tol_;

  public:

  PositionConstraint(std::string name, std::string name1, double tolerance):
  name1_(name1), tol_(tolerance)
  {
    name_ = name;
    names_.resize(1);
    names_[0] = name1_;
    force_ = 0.0;
    weight_ = 1.0;
  }
  ~PositionConstraint(void) {}

  std::string type(void) const { return "position"; }
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
