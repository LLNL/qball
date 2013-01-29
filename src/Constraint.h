////////////////////////////////////////////////////////////////////////////////
//
//  Constraint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Constraint.h,v 1.3 2010/01/16 01:26:35 draeger1 Exp $

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <string>
#include <vector>
#include <cassert>

class AtomSet;

class Constraint
{
  protected:

  std::string name_;          // constraint name
  std::vector<std::string> names_; // names of atoms involved in the constraint

  public:

  virtual ~Constraint(void){}
  virtual std::string type(void) const = 0;
  virtual double value(void) const = 0;
  virtual double velocity(void) const = 0;
  virtual double force(void) const = 0;
  virtual double weight(void) const = 0;
  virtual double tolerance(void) const = 0;
  virtual void set_value(double value) = 0;
  virtual void set_velocity(double velocity) = 0;
  virtual bool enforce_r(const std::vector<std::vector<double> > &r0,
                         std::vector<std::vector<double> > &rp) const = 0;
  virtual bool enforce_v(const std::vector<std::vector<double> > &r0,
                         std::vector<std::vector<double> > &v0) const = 0;
  virtual void compute_force(const std::vector<std::vector<double> > &r0,
                             const std::vector<std::vector<double> > &f) = 0;
  virtual void update(double dt) = 0;
  virtual void setup(const AtomSet& atoms) = 0;
  virtual int dofs(void) const = 0;
  virtual std::ostream& print(std::ostream &os) = 0;
  std::string name(void) const { return name_; }
  std::string names(int i) const
  {
    assert( i >= 0 && i < names_.size() );
    return names_[i];
  }
};
std::ostream& operator << ( std::ostream &os, Constraint &c );
#endif


