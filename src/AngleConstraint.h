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
//  AngleConstraint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AngleConstraint.h,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#ifndef ANGLECONSTRAINT_H
#define ANGLECONSTRAINT_H

#include "Constraint.h"
#include "D3vector.h"
#include <cassert>

class AtomSet;

class AngleConstraint : public Constraint
{
  std::string name1_, name2_, name3_;
  int    ia1_, ia2_, ia3_, is1_, is2_, is3_;
  double m1_, m2_, m3_, m1_inv_, m2_inv_, m3_inv_;
  double angle_, velocity_, force_, weight_, tol_;
  double sigma(D3vector a, D3vector b, D3vector c) const;
  void grad_sigma(const D3vector &r1, const D3vector &r2, const D3vector &r3,
                  D3vector &g1, D3vector &g2,D3vector &g3) const;
  double bond_angle(D3vector a, D3vector b, D3vector c) const;

  public:

  AngleConstraint(std::string name, std::string name1,
                  std::string name2, std::string name3,
                  double angle, double velocity, double tolerance):
  name1_(name1), name2_(name2), name3_(name3),
  velocity_(velocity),
  tol_(tolerance), m1_(0.0), m2_(0.0), m3_(0.0)
  {
    set_value(angle);
    name_ = name;
    names_.resize(3);
    names_[0] = name1_;
    names_[1] = name2_;
    names_[2] = name3_;
    force_ = 0.0;
    weight_ = 1.0;
  }
  ~AngleConstraint(void) {}

  std::string type(void) const { return "angle"; }
  double value(void) const { return angle_; }
  double velocity(void) const { return velocity_; }
  double force(void) const { return force_; }
  double weight(void) const { return weight_; }
  double tolerance(void) const { return tol_; }
  void set_value(double value)
  {
    angle_ = value;
    if ( angle_ <   0.0 ) angle_ =   0.0;
    if ( angle_ > 180.0 ) angle_ = 180.0;
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
