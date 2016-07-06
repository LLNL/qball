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
// UnitCell.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef UNITCELL_H
#define UNITCELL_H

#include "D3vector.h"
#include <valarray>
using namespace std;

class UnitCell
{
  private:

  D3vector a_[3];
  D3vector b_[3];
  double volume_;
  D3vector an_[13];
  D3vector bn_[13];
  double an2h_[13];
  double bn2h_[13];
  
  // 3x3 matrix forms
  double amat_[9];
  double bmat_[9];
  // 3x3 matrix form of inverse
  double amat_inv_[9];
  // 3x3 matrix form of inverse transpose
  double amat_inv_t_[9];
  
  public:

  const D3vector& a(int i) const { return a_[i]; }
  const D3vector& b(int i) const { return b_[i]; }
  
  UnitCell(void) { set(D3vector(0,0,0),D3vector(0,0,0),D3vector(0,0,0)); }  
  explicit UnitCell(const D3vector& a0, const D3vector& a1, const D3vector& a2)
  { set(a0,a1,a2); }
  
  void set(const D3vector& a0, const D3vector& a1, const D3vector& a2);
  double volume(void) const { return volume_; }
  
  const double* amat(void) const { return &amat_[0]; }
  const double* bmat(void) const { return &bmat_[0]; }
  const double* amat_inv(void) const { return &amat_inv_[0]; }
  double amat(int ij) const { return amat_[ij]; }
  double bmat(int ij) const { return bmat_[ij]; }
  double amat_inv(int ij) const { return amat_inv_[ij]; }
  
  // 3x3 matrix vector multiply Z = X Y where X is a 3x3 matrix, Y,Z 3-vectors
  void vecmult3x3(const double* x, const double* y, double *z) const;

  // 3x3 sym matrix vector multiply Z = X Y where X is a sym 3x3 matrix,
  // Y,Z 3-vectors
  void vecsmult3x3(const double* x, const double* y, double *z) const;
  
  // 3x3 matrix matrix multiply Z = X Y where X, Y are 3x3 matrices
  void matmult3x3(const double* x, const double* y, double *z) const;
  // Z = X Y where X is a symmetric 3x3 matrix and Y a general 3x3 matrix
  // uses only the first 6 elements of array xs
  // where xs[0] = x00, xs[1] = x11, xs[2] = x22,
  // xs[3] = x10, xs[4] = x21, xs[5] = x20
  void smatmult3x3(const double* xs, const double* y, double *z) const;
  void compute_deda(const valarray<double>& sigma, valarray<double>& deda) const;
  
  void cart_to_crystal(const double* scart, double* scryst) const;
  D3vector cart_to_crystal(const D3vector& v) const;
  void crystal_to_cart(const double* scryst, double* scart) const;
  D3vector crystal_to_cart(const D3vector& v) const;
  bool in_ws(const D3vector& v) const;
  double min_wsdist(void) const;
  void fold_in_ws(D3vector& v) const;
  bool in_bz(const D3vector& k) const;
  void fold_in_bz(D3vector& k) const;
  
  bool encloses(const UnitCell& c) const;
  bool contains(D3vector v) const;
  
  void print(ostream& os) const;  
  void printsys(ostream& os, string setcmd) const;  
  bool operator==(const UnitCell& c) const;
  bool operator!=(const UnitCell& c) const;
};
ostream& operator << ( ostream& os, const UnitCell& cell );
#endif
