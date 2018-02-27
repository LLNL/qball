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
// MLWFTransform.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef MLWFTRANSFORM_H
#define MLWFTRANSFORM_H

#include <vector>
#include <complex>
class SlaterDet;
class UnitCell;
class DoubleMatrix;
#include <math/d3vector.h>
#include "BasisMapping.h"

class MLWFTransform
{
  private:

  const SlaterDet& sd_;
  const UnitCell& cell_;
  const Context& ctxt_;

  BasisMapping bm_;
  std::vector<DoubleMatrix*> a_;  // cosine and sine matrices
  DoubleMatrix* u_;               // orthogonal transformation
  std::vector<std::vector<double> > adiag_; // diagonal elements

  public:

  void compute_transform(void);
  void compute_sincos(const int n, const std::complex<double>* f,
    std::complex<double>* fc, std::complex<double>* fs);
  void apply_transform(SlaterDet& sd);

  double spread2(int i, int j);
  double spread2(int i);
  double spread2(void);
  double spread(int i);
  double spread(void);
  D3vector center(int i);
  D3vector dipole(void);

  MLWFTransform(const SlaterDet& sd);
  ~MLWFTransform(void);
};
#endif

// Local Variables:
// mode: c++
// End:
