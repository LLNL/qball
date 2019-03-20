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
// VectorPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>
#include <math/d3vector.h>
#include <qball/Basis.h>
#include "UnitCell.h"

#ifndef VECTORPOTENTIAL_H
#define VECTORPOTENTIAL_H

using namespace std;

class VectorPotential {

public:

  enum class Dynamics {
    NONE,
    POLARIZATION
  };

  VectorPotential(Dynamics dyn, const D3vector & initial_value):
    dynamics_(dyn),
    value_(initial_value),
    value2_(norm(value_))
  {
    velocity_ = D3vector(0.0, 0.0, 0.0);
    accel_ = D3vector(0.0, 0.0, 0.0);
  }

  double * get_kpgpa(const Basis & basis) const {
    const double * kpgpa2 = get_kpgpa2(basis);
    double * kpgpa = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      kpgpa[ig] = sqrt(kpgpa2[ig]);
    }
    delete [] kpgpa2;
    return kpgpa;
  }
  
  double * get_kpgpa2(const Basis & basis) const {
    double * kpgpa2 = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      kpgpa2[ig] = basis.kpg2_ptr()[ig] + value2();
      kpgpa2[ig] -= 2 * value_[0]*basis.kpgx_ptr(0)[ig];
      kpgpa2[ig] -= 2 * value_[1]*basis.kpgx_ptr(1)[ig];
      kpgpa2[ig] -= 2 * value_[2]*basis.kpgx_ptr(2)[ig];
    }
    return kpgpa2;
  }

  double * get_kpgpai(const Basis & basis) const {
    const double * kpgpa2 = get_kpgpa2(basis);
    double * kpgpai = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      if(kpgpa2[ig] > 0.0){
	kpgpai[ig] = 1.0/sqrt(kpgpa2[ig]);
      } else {
	kpgpai[ig] = 0.0;
      }
    }
    delete [] kpgpa2;
    return kpgpai;
  }
  
  double * get_kpgpax(const Basis & basis, int j) const {
    if(j == 0){ 
      double * kpgpax = new double[3*basis.localsize()];
      for(int j2 = 0; j2 < 3; j2++){
	for(int ig = 0; ig < basis.localsize(); ig++){
	  kpgpax[j2*basis.localsize() + ig] = basis.kpgx_ptr(j2)[ig] - value_[j2];
	}
      }
      return kpgpax;
    } else {
      double * kpgpax = new double[basis.localsize()];
      for(int ig = 0; ig < basis.localsize(); ig++){
	kpgpax[ig] = basis.kpgx_ptr(j)[ig] - value_[j];
      }
      return kpgpax;
    }
  }
  
  const D3vector & value() const {
    return value_;
  }

  const double & value2() const {
    return value2_;
  }

  void calculate_acceleration(const double & dt, const D3vector& total_current, const UnitCell & cell){
    //update the velocity to time t - dt/2
    velocity_ += 0.5*dt*accel_;

    if(dynamics_ == Dynamics::POLARIZATION){
      accel_ = -4.0*M_PI*total_current/cell.volume();
    } else {
      accel_ = D3vector(0.0, 0.0, 0.0);
    }

    //update the velocity to time t
    velocity_ += 0.5*dt*accel_;
  }
  
  void propagate(const double & dt){
    value_ += dt*velocity_ + 0.5*dt*dt*accel_;
  }

  
private:
  Dynamics dynamics_;
  D3vector value_;
  D3vector velocity_;
  D3vector accel_;
  double value2_;
  
};
#endif

// Local Variables:
// mode: c++
// End:

