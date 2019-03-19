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

#ifndef VECTORPOTENTIAL_H
#define VECTORPOTENTIAL_H

using namespace std;

class VectorPotential {

public:
  VectorPotential(const D3vector & initial_value):
    value_(initial_value),
    value2_(norm(value_))
  {
  }

  const double * get_kpgpa2(const Basis & basis) const {
    double * kpgpa2 = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      kpgpa2[ig] = basis.kpg2_ptr()[ig] + value2();
      kpgpa2[ig] -= value_[0]*basis.kpgx_ptr(0)[ig];
      kpgpa2[ig] -= value_[1]*basis.kpgx_ptr(1)[ig];
      kpgpa2[ig] -= value_[2]*basis.kpgx_ptr(2)[ig];
    }
    return kpgpa2;
  }

  const double * get_kpgpax(const Basis & basis, int j) const {
    double * kpgpax = new double[basis.localsize()];
    for(int ig = 0; ig < basis.localsize(); ig++){
      kpgpax[ig] = basis.kpgx_ptr(j)[ig] - value_[j];
    }
    return kpgpax;
  }

  const D3vector & value() const {
    return value_;
  }

  const double & value2() const {
    return value2_;
  }
  
private:
  D3vector value_;
  double value2_;
  
};
#endif

// Local Variables:
// mode: c++
// End:

