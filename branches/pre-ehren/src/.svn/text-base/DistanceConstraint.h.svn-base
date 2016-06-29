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
//  DistanceConstraint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: DistanceConstraint.h,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#ifndef DISTANCECONSTRAINT_H
#define DISTANCECONSTRAINT_H

#include "Constraint.h"
#include <cassert>
#include <cmath> // fabs

class AtomSet;

class DistanceConstraint : public Constraint
{
  std::string name1_, name2_;
  int    ia1_, ia2_, is1_, is2_;
  double m1_, m2_, m1_inv_, m2_inv_;
  double distance_, velocity_, force_, weight_, tol_;

  public:

  DistanceConstraint(std::string name, std::string name1, std::string name2,
                     double distance, double velocity, double tolerance):
  name1_(name1), name2_(name2), distance_(distance),
  velocity_(velocity), tol_(tolerance), m1_(0.0), m2_(0.0)
  {
    name_ = name;
    names_.resize(2);
    names_[0] = name1_;
    names_[1] = name2_;
    force_ = 0.0;
    weight_ = 1.0;
  }
  ~DistanceConstraint(void) {}

  std::string type(void) const { return "distance"; }
  double value(void) const { return distance_; }
  double velocity(void) const { return velocity_; }
  double force(void) const { return force_; }
  double weight(void) const { return weight_; }
  double tolerance(void) const { return tol_; }
  void set_value(double value)
  {
    distance_ = fabs(value);
  }
  void set_velocity(double velocity)
  {
    velocity_ = velocity;
  }

  void setup(const AtomSet& atoms);
  int dofs(void) const { return 1; }
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
