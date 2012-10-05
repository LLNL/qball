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
// ConstraintSet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintSet.h,v 1.10 2010-02-20 23:13:02 fgygi Exp $

#ifndef CONSTRAINTSET_H
#define CONSTRAINTSET_H

#include <vector>
#include <string>

class Atom;
class AtomSet;
class Constraint;
class Context;

class ConstraintSet
{
  private:

  const Context& ctxt_;
  std::vector<Constraint *> constraint_list;
  // ndofs_: total number of degrees of freedom blocked by the constraints
  int ndofs_;

  public:

  ConstraintSet(const Context& ctxt) : ctxt_(ctxt), ndofs_(0) {}
  ~ConstraintSet();
  bool define_constraint(AtomSet &atoms, int argc, char **argv);
  bool set_constraint(int argc, char **argv);
  bool delete_constraint(int argc, char **argv);
  void list_constraints(std::ostream &os);
  int size(void) const { return constraint_list.size(); }
  int ndofs(void) const { return ndofs_; }
  void enforce(AtomSet& atoms);
  void enforce_r(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &rp);
  void enforce_v(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &v0);
  void compute_forces(const std::vector<std::vector<double> > &r0,
                      const std::vector<std::vector<double> > &f);
  void update_constraints(double dt);
  void setup(AtomSet& atoms);
  void reset(void);
};
#endif
