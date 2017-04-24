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
// ConstraintSet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintSet.C,v 1.11 2010-02-20 23:13:02 fgygi Exp $

#include <config.h>

#include "ConstraintSet.h"
#include "PositionConstraint.h"
#include "DistanceConstraint.h"
#include "AngleConstraint.h"
#include "TorsionConstraint.h"
//#include "MultiDistanceConstraint.h"
#include "Atom.h"
#include "AtomSet.h"
#include "Context.h"
#include <iostream>
#include <iomanip>
using namespace std;

const int constraints_maxiter = 10;

////////////////////////////////////////////////////////////////////////////////
ConstraintSet::~ConstraintSet(void)
{
  for ( int ic = 0; ic < constraint_list.size(); ic++ )
    delete constraint_list[ic];
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::define_constraint(AtomSet &atoms, int argc, char **argv)
{
  enum constraint_type { unknown, position_type, distance_type,
                         angle_type, torsion_type }
    type = unknown;
  const double position_tolerance = 1.0e-7;
  const double distance_tolerance = 1.0e-7;
  const double angle_tolerance = 1.0e-4;
  const bool oncoutpe = ctxt_.oncoutpe();

  // argv[0] == "constraint"
  // argv[1] == "define"
  // argv[2] == {"position", "distance", "angle", "torsion"}
  // argv[3] == constraint name
  // argv[4-(5,6,7)] == atom names
  // argv[{5,6,7}] == {distance,angle,angle}
  // argv[{6,7,8}] == velocity

  if ( argc < 2 )
  {
    if ( oncoutpe )
    {
      cout << " Use: constraint define position constraint_name atom_name"
           << endl;
      cout << " Use: constraint define distance constraint_name "
           << "atom_name1 atom_name2 distance_value [velocity]"
           << endl;
      cout << "      constraint define angle constraint_name "
           << "name1 name2 name3 angle_value [velocity]"
           << endl;
      cout << "      constraint define torsion constraint_name "
           << "name1 name2 name3 name4 angle_value"
           << " [velocity] "
           << endl;
    }
    return false;
  }
  const string constraint_type = argv[2];
  if ( constraint_type == "position" )
  {
    type = position_type;
  }
  else if ( constraint_type == "distance" )
  {
    type = distance_type;
  }
  else if ( constraint_type == "angle" )
  {
    type = angle_type;
  }
  else if ( constraint_type == "torsion" )
  {
    type = torsion_type;
  }
  else
  {
    if ( oncoutpe )
      cout << " Incorrect constraint type " << constraint_type << endl;
    return false;
  }

  if ( type == position_type )
  {
    // define position name A

    if ( argc != 5 )
    {
      if ( oncoutpe )
        cout << " Incorrect number of arguments for position constraint"
             << endl;
      return false;
    }
    string name = argv[3];
    string atom_name = argv[4];

    Atom *a1 = atoms.findAtom(atom_name);

    if ( a1 == 0 )
    {
      if ( oncoutpe )
      {
        cout << " ConstraintSet: could not find atom " << atom_name << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }

    // check if constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint with same name or with same atom is defined
      if ( pc->type() == "position" )
        found = ( pc->name() == name ) || ( pc->names(0) == atom_name );
    }

    if ( found )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: constraint is already defined:\n"
             << " cannot define constraint" << endl;
      return false;
    }
    else
    {
      PositionConstraint *c =
        new PositionConstraint(name,atom_name,position_tolerance);
      constraint_list.push_back(c);
    }
  }
  else if ( type == distance_type )
  {
    // define distance name A B value
    // define distance name A B value velocity

    if ( argc < 7 || argc > 8 )
    {
      if ( oncoutpe )
        cout << " Incorrect number of arguments for distance constraint"
             << endl;
      return false;
    }
    double distance, velocity=0.0;
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);

    if ( a1 == 0 || a2 == 0 )
    {
      if ( oncoutpe )
      {
        if ( a1 == 0 )
          cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ConstraintSet: could not find atom " << name2 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }
    if ( name1 == name2 )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: cannot define distance constraint between "
             << name1 << " and " << name2 << endl;
      return false;
    }

    distance = atof(argv[6]);
    if ( argc == 8 )
    {
      velocity = atof(argv[7]);
    }

    if ( distance <= 0.0 )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: distance must be positive" << endl
             << " ConstraintSet: could not define constraint" << endl;
      return false;
    }

    // check if constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint (name1,name2) or (name2,name1) is defined
      if ( pc->type() == "distance" )
        found = ( pc->names(0) == name1 && pc->names(1) == name2 ) ||
                ( pc->names(1) == name1 && pc->names(0) == name2 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: constraint is already defined:\n"
             << " cannot define constraint" << endl;
      return false;
    }
    else
    {
      DistanceConstraint *c =
        new DistanceConstraint(name,name1,name2,distance,
                               velocity,distance_tolerance);

      constraint_list.push_back(c);
    }
  }
  else if ( type == angle_type )
  {
    // constraint define angle name A B C value
    // constraint define angle name A B C value velocity

    if ( argc < 8  || argc > 9 )
    {
      if ( oncoutpe )
        cout << " Incorrect number of arguments for angle constraint"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];
    string name3 = argv[6];

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
    Atom *a3 = atoms.findAtom(name3);

    if ( a1 == 0 || a2 == 0 || a3 == 0 )
    {
      if ( oncoutpe )
      {
        if ( a1 == 0 )
          cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ConstraintSet: could not find atom " << name2 << endl;
        if ( a3 == 0 )
          cout << " ConstraintSet: could not find atom " << name3 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }

    if ( name1 == name2 || name1 == name3 || name2 == name3)
    {
      if ( oncoutpe )
        cout << " ConstraintSet: cannot define angle constraint between "
             << name1 << " " << name2 << " and " << name3 << endl;
      return false;
    }

    const double angle = atof(argv[7]);
    double velocity = 0.0;
    if ( argc == 9 )
    {
      velocity = atof(argv[8]);
    }

    if ( angle < 0.0 || angle > 180.0 )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: angle must be in [0,180]" << endl
             << " ConstraintSet: could not define constraint" << endl;
      return false;
    }

    // check if equivalent constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint (name1,name2,name3) or
      // (name3,name2,name1) is defined
      if ( pc->type() == "angle" )
        found = ( pc->names(0) == name1 &&
                  pc->names(1) == name2 &&
                  pc->names(2) == name3 ) ||
                ( pc->names(0) == name3 &&
                  pc->names(1) == name2 &&
                  pc->names(2) == name1 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( oncoutpe )
        cout << " ConstraintSet:set_constraint: an angle constraint "
             << name1 << " " << name2 << " " << name3
             << " was found" << endl
             << " ConstraintSet: cannot define constraint" << endl;
      return false;
    }
    else
    {
      AngleConstraint *c =
      new AngleConstraint(name, name1,name2,name3,angle,
        velocity,angle_tolerance);
      constraint_list.push_back(c);
    }
  }
  else if ( type == torsion_type )
  {
    // constraint define torsion name A B C D angle
    // constraint define torsion name A B C D angle velocity

    if ( argc < 9  || argc > 10 )
    {
      if ( oncoutpe )
        cout << " Incorrect number of arguments for torsion constraint"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];
    string name3 = argv[6];
    string name4 = argv[7];

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
    Atom *a3 = atoms.findAtom(name3);
    Atom *a4 = atoms.findAtom(name4);

    if ( a1 == 0 || a2 == 0 || a3 == 0 || a4 == 0 )
    {
      if ( oncoutpe )
      {
        if ( a1 == 0 )
          cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ConstraintSet: could not find atom " << name2 << endl;
        if ( a3 == 0 )
          cout << " ConstraintSet: could not find atom " << name3 << endl;
        if ( a4 == 0 )
          cout << " ConstraintSet: could not find atom " << name4 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }
    if ( name1 == name2 || name1 == name3 || name1 == name4 ||
         name2 == name3 || name2 == name4 || name3 == name4 )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: cannot define torsion constraint using "
             << name1 << " " << name2 << " " << name3 << " " << name4
             << endl;
      return false;
    }

    double angle = atof(argv[8]);
    if ( angle > 180.0 )
      while ( angle > 180.0 ) angle -= 360.0;
    else if ( angle < -180.0 )
      while ( angle < -180.0 ) angle += 360.0;

    double velocity = 0.0;
    if ( argc == 10 )
    {
      velocity = atof(argv[9]);
    }

    // check if equivalent constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if an equivalent constraint (name1,name2,name3,name4) or
      // (name4,name3,name2,name1) is defined
      if ( pc->type() == "angle" )
        found = ( pc->names(0) == name1 &&
                  pc->names(1) == name2 &&
                  pc->names(2) == name3 &&
                  pc->names(3) == name4 ) ||
                ( pc->names(0) == name4 &&
                  pc->names(1) == name3 &&
                  pc->names(2) == name2 &&
                  pc->names(3) == name1 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( oncoutpe )
        cout << " ConstraintSet: a torsion constraint "
             << name1 << " " << name2 << " " << name3 << " " << name4
             << " is already defined" << endl
             << " ConstraintSet: cannot define constraint" << endl;
      return false;
    }
    else
    {
      TorsionConstraint *c =
      new TorsionConstraint(name,name1,name2,name3,name4,
                            angle,velocity,angle_tolerance);
      constraint_list.push_back(c);
    }
  }
  else
  {
    if ( oncoutpe )
      cout << " ConstraintSet::set_constraint: internal error" << endl;
    return false;
  }

  // update total number of blocked degrees of freedom
  ndofs_ += constraint_list.back()->dofs();

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::set_constraint(int argc, char **argv)
{
  const bool oncoutpe = ctxt_.oncoutpe();
  assert(argc==4||argc==5);
  // argv[0] == "constraint"
  // argv[1] == "set"
  // argv[2] == constraint_name
  // argv[3] == value
  // argv[4] (optional) == velocity
  string name = argv[2];
  const double value = atof(argv[3]);
  double velocity = 0.0;
  const bool set_velocity = ( argc == 5 );
  if ( set_velocity ) velocity = atof(argv[4]);

    // check if constraint is already defined
  bool found = false;
  vector<Constraint*>::iterator i = constraint_list.begin();
  while ( !found && i != constraint_list.end() )
  {
    Constraint *pc = *i;
    assert(pc != 0);

    if ( pc->name() == name )
    {
      found = true;
      pc->set_value(value);
      if ( set_velocity ) pc->set_velocity(velocity);
    }
    i++;
  }

  if ( !found )
  {
    if ( oncoutpe )
      cout << " ConstraintSet: no such constraint" << endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::delete_constraint(int argc, char **argv)
{
  assert(argc==3);
  // argv[0] == "constraint"
  // argv[1] == "delete"
  // argv[2] == constraint_name
  string name = argv[2];
  const bool oncoutpe = ctxt_.oncoutpe();

  bool found = false;
  // note next loop in reverse: avoid use of invalidated iterators
  // after erase operation

  vector<Constraint*>::iterator i = constraint_list.begin();
  while ( !found && i != constraint_list.end() )
  {
    Constraint *pc = *i;
    assert(pc != 0);

    // note structure of if else test to avoid incrementing
    // invalidated iterator after erase (see Meyers STL, p.45)
    if ( pc->name() == name )
    {
      found = true;

      // update total number of blocked degrees of freedom
      ndofs_ -= pc->dofs();

      delete pc;

      // remove constraint pointer from the list
      // note: iterator is incremented before erasing, remains valid
      constraint_list.erase(i++);
    }
    else
    {
      i++;
    }
  }

  if ( !found )
  {
    if ( oncoutpe ) cout << " No such constraint" << endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::list_constraints(ostream &os)
{
  if ( !constraint_list.empty() )
  {
    os << " <constraint_set>" << endl;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      Constraint *c = constraint_list[i];
      os << *c << endl;
    }
    os << " </constraint_set>" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce(AtomSet& atoms)
{
  vector<vector<double> > r0,rp,v0;
  setup(atoms);
  atoms.get_positions(r0);
  rp=r0;
  atoms.get_velocities(v0);
  enforce_r(r0,rp);
  atoms.set_positions(rp);
  enforce_v(r0,v0);
  atoms.set_velocities(v0);
}
////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce_r(const vector<vector<double> > &r0,
                              vector<vector<double> > &rp)
{
  const bool oncoutpe = ctxt_.oncoutpe();
  int iter = 0;
  bool done = false;
  while ( !done && (iter < constraints_maxiter) )
  {
    done = true;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      Constraint *c = constraint_list[i];
      bool b = c->enforce_r(r0,rp);
      done &= b;
    }
    iter++;
  }

  if ( !done )
  {
    if ( oncoutpe )
      cout << " ConstraintSet: could not enforce position constraints in "
           << constraints_maxiter << " iterations" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce_v(const vector<vector<double> > &r0,
                              vector<vector<double> > &v0)
{
  const bool oncoutpe = ctxt_.oncoutpe();
  int iter = 0;
  bool done = false;
  while ( !done && (iter < constraints_maxiter) )
  {
    done = true;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      bool b = constraint_list[i]->enforce_v(r0,v0);
      done &= b;
    }
    iter++;
  }

  if ( !done )
  {
    if ( oncoutpe )
      cout << " ConstraintSet: could not enforce velocity constraints in "
           << constraints_maxiter << " iterations" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::compute_forces(const vector<vector<double> > &r0,
 const vector<vector<double> > &f)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    constraint_list[i]->compute_force(r0,f);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::update_constraints(double dt)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    constraint_list[i]->update(dt);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::setup(AtomSet& atoms)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    constraint_list[i]->setup(atoms);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::reset(void)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
    delete constraint_list[i];
  ndofs_ = 0;
  constraint_list.resize(0);
}
