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
// IonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef IONICSTEPPER_H
#define IONICSTEPPER_H

#include "Sample.h"
#include "Species.h"
#include <vector>

class IonicStepper
{
  protected:

  Sample& s_;
  AtomSet& atoms_;
  ConstraintSet& constraints_;
  double                    dt_;
  int                       nsp_;
  int                       nspqm_;
  int                       nspmm_;
  int                       ndofs_;
  std::vector<int>               na_;  // number of atoms per species na_[nsp_]
  std::vector<std::vector< double> >  r0_; // r0_[nsp_][3*na_]
  std::vector<std::vector< double> >  rp_; // rp_[nsp_][3*na_]
  std::vector<std::vector< double> >  rm_; // rm_[nsp_][3*na_]
  std::vector<std::vector< double> >  v0_; // v0_[nsp_][3*na_]
  std::vector<double>            pmass_;   // pmass_[nsp_]

  public:

  IonicStepper (Sample& s) : s_(s), atoms_(s.atoms),
    constraints_(s.constraints), dt_(s.ctrl.dt)
  {
    ndofs_ = 3 * atoms_.size() - constraints_.size();
    // if there are more constraints than dofs, set ndofs_ to zero
    if ( ndofs_ < 0 ) ndofs_ = 0;
    nspqm_ = atoms_.nsp();
    nspmm_ = atoms_.nsp_mm();
    nsp_ = nspqm_ + nspmm_;
    na_.resize(nsp_);
    r0_.resize(nsp_);
    rp_.resize(nsp_);
    v0_.resize(nsp_);
    pmass_.resize(nsp_);
    for ( int is = 0; is < nsp_; is++ )
    {
      const int nais = atoms_.na(is);
      na_[is] = nais;
      r0_[is].resize(3*nais);
      rp_[is].resize(3*nais);
      v0_[is].resize(3*nais);
    }
    for ( int is = 0; is < nspqm_; is++ )
      pmass_[is] = atoms_.species_list[is]->mass() * 1822.89;
    for ( int is = 0; is < nspmm_; is++ )
      pmass_[is+nspqm_] = atoms_.mmspecies_list[is]->mass() * 1822.89;

    atoms_.get_positions(r0_,false);
    atoms_.get_velocities(v0_,false);
  }

  double r0(int is, int i) const { return r0_[is][i]; }
  double v0(int is, int i) const { return v0_[is][i]; }
  const std::vector<std::vector<double> >& r0(void) const { return r0_; }
  const std::vector<std::vector<double> >& v0(void) const { return v0_; }
  const std::vector<double>& pmass(void) const { return pmass_; }

  void setup_constraints(void)
  {
    constraints_.setup(atoms_);
  }
  virtual void compute_r(double e0,
    const std::vector<std::vector< double> >& f0) = 0;
  virtual void compute_v(double e0,
    const std::vector<std::vector< double> >& f0) = 0;
  virtual void reset(void) {}
  virtual double ekin(void) const { return 0.0; }
  virtual double temp(void) const { return 0.0; }

  virtual ~IonicStepper() {}

};
#endif
