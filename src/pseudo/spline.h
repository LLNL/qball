////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory.
// Written by Xavier Andrade (xavier@llnl.gov), Erik Draeger
// (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
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

#include <config.h>
/*******************************************************************************
 *
 * spline.h
 *
 ******************************************************************************/
#ifndef SPLINE_H
#define SPLINE_H

#include <vector>

#define SPLINE_FLAT_BC 0.0       /* Flat boundary condition (y'=0) */
#define SPLINE_NATURAL_BC 1.e31  /* Natural boundary condition (Y"=0) */


void spline(const double *x, const double *y, int n, double yp1, double ypn, double *y2);
void splint (const double *xa, const double *ya, const double *y2a, int n, double x, double *y);
void splintd (const double *xa, const double *ya, const double *y2a, int n, double x, double *y, double *dy);

class Spline {

 public:

  Spline(){
  }
  
  void fit(double *x, double *y, int n, double yp1, double ypn){
    x_.resize(n);
    y_.resize(n);
    y2_.resize(n);
    
    for(int ii = 0; ii < n; ii++){
      x_[ii] = x[ii];
      y_[ii] = y[ii];
    }
    spline(x, y, n, yp1, ypn, &y2_[0]);
  }

  double value(const double & x) const {
    double y;
    splint(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y);
    return y;
  }
  
  void derivative(const double & x, double & y, double & dy) const {
    splintd(&x_[0], &y_[0], &y2_[0], x_.size(), x, &y, &dy);
  }
 private :

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> y2_;
  
};

#endif

// Local Variables:
// mode: c++
// End:
