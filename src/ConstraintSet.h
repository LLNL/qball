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
