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
// tvalarray.C
// causes error messages related to valarray headers using xlC 7.0
////////////////////////////////////////////////////////////////////////////////
#include <valarray>
using namespace std;
class A
{
  private:
  double v_;
  double amat_inv_t_[9];
  public:
  void smatmult3x3(const double* xs, const double* y, double *z) const;
  void f(const valarray<double>& sigma, valarray<double>& deda) const;
};
////////////////////////////////////////////////////////////////////////////////
void A::f(const valarray<double>& sigma, valarray<double>& deda) const
{
  // local copy of sigma 
  valarray<double> sigma_loc(6);
  sigma_loc = sigma;
  smatmult3x3(&sigma_loc[0],&amat_inv_t_[0],&deda[0]);
  deda *= -v_;
}

////////////////////////////////////////////////////////////////////////////////
void A::smatmult3x3(const double* xs, const double* y, double *z) const
{
  //  | z0 z3 z6 |     | xs0 xs3 xs5 |     | y0 y3 y6 |
  //  | z1 z4 z7 |  =  | xs3 xs1 xs4 |  *  | y1 y4 y7 |
  //  | z2 z5 z8 |     | xs5 xs4 xs2 |     | y2 y5 y8 |
  
  const double z00 = xs[0]*y[0]+xs[3]*y[1]+xs[5]*y[2];
  const double z10 = xs[3]*y[0]+xs[1]*y[1]+xs[4]*y[2];
  const double z20 = xs[5]*y[0]+xs[4]*y[1]+xs[2]*y[2];
  
  const double z01 = xs[0]*y[3]+xs[3]*y[4]+xs[5]*y[5];
  const double z11 = xs[3]*y[3]+xs[1]*y[4]+xs[4]*y[5];
  const double z21 = xs[5]*y[3]+xs[4]*y[4]+xs[2]*y[5];
  
  const double z02 = xs[0]*y[6]+xs[3]*y[7]+xs[5]*y[8];
  const double z12 = xs[3]*y[6]+xs[1]*y[7]+xs[4]*y[8];
  const double z22 = xs[5]*y[6]+xs[4]*y[7]+xs[2]*y[8];
  
  z[0] = z00;
  z[1] = z10;
  z[2] = z20;
  
  z[3] = z01;
  z[4] = z11;
  z[5] = z21;
  
  z[6] = z02;
  z[7] = z12;
  z[8] = z22; 
}
