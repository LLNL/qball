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
// SymOp.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "SymOp.h"
#include "D3vector.h"
#include <iomanip>
using namespace std;

SymOp::SymOp (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33) {
  s11_ = s11;  s12_ = s12;  s13_ = s13;
  s21_ = s21;  s22_ = s22;  s23_ = s23;
  s31_ = s31;  s32_ = s32;  s33_ = s33;
  symcub_.resize(3);
  symxtal_.resize(3);
  for (int i=0; i<3; i++) {
    symcub_[i].resize(3);
    symxtal_[i].resize(3);
  }
  symcub_[0][0] = s11;  symcub_[0][1] = s12;  symcub_[0][2] = s13;  
  symcub_[1][0] = s21;  symcub_[1][1] = s22;  symcub_[1][2] = s23;  
  symcub_[2][0] = s31;  symcub_[2][1] = s32;  symcub_[2][2] = s33;  
  ftrans1_ = 0.;  ftrans2_ = 0.;  ftrans3_ = 0.;
}

SymOp::SymOp (double s11, double s12, double s13, double s21, double s22, double s23, double s31, double s32, double s33, double ftrans1, double ftrans2, double ftrans3) {
  s11_ = s11;  s12_ = s12;  s13_ = s13;
  s21_ = s21;  s22_ = s22;  s23_ = s23;
  s31_ = s31;  s32_ = s32;  s33_ = s33;
  symcub_.resize(3);
  symxtal_.resize(3);
  for (int i=0; i<3; i++) {
    symcub_[i].resize(3);
    symxtal_[i].resize(3);
  }
  symcub_[0][0] = s11;  symcub_[0][1] = s12;  symcub_[0][2] = s13;  
  symcub_[1][0] = s21;  symcub_[1][1] = s22;  symcub_[1][2] = s23;  
  symcub_[2][0] = s31;  symcub_[2][1] = s32;  symcub_[2][2] = s33;  
  ftrans1_ = ftrans1;  ftrans2_ = ftrans2;  ftrans3_ = ftrans3;
}

void SymOp::convertToXtal(const UnitCell& uc) {
  // convert symmetry operation to crystal coordinates:
  //
  //   sym_xtal^T = a^-1 sym a
  //      where a is a 3x3 matrix of cell vectors
  vector<D3vector> a;
  a.resize(3);
  a[0] = uc.a(0);
  a[1] = uc.a(1);
  a[2] = uc.a(2);

  const double a11 = a[0].x; const double a12 = a[1].x; const double a13 = a[2].x;
  const double a21 = a[0].y; const double a22 = a[1].y; const double a23 = a[2].y;
  const double a31 = a[0].z; const double a32 = a[1].z; const double a33 = a[2].z;
  const double deta = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32;
  if (deta == 0.0) {
    cout << "ERROR:  determinant of lattice vector matrix is zero!" << endl;;
    exit(1);
  }

  const double invdeta = 1.0/deta;
  vector< vector<double > > ainv;
  vector< vector<double > > symmulta;
  ainv.resize(3);
  symmulta.resize(3);
  for (int i=0; i<3; i++) {
    ainv[i].resize(3);
    symmulta[i].resize(3);
  }

  // compute inverse of lattice vector matrix
  ainv[0][0] = invdeta*(a22*a33-a23*a32);
  ainv[0][1] = invdeta*(a13*a32-a12*a33);
  ainv[0][2] = invdeta*(a12*a23-a13*a22);
  ainv[1][0] = invdeta*(a23*a31-a21*a33);
  ainv[1][1] = invdeta*(a11*a33-a13*a31);
  ainv[1][2] = invdeta*(a13*a21-a11*a23);
  ainv[2][0] = invdeta*(a21*a32-a31*a22);
  ainv[2][1] = invdeta*(a12*a31-a11*a32);
  ainv[2][2] = invdeta*(a11*a22-a12*a21);

  for (int j=0; j<3; j++) {
    for (int k=0; k<3; k++) {
      symmulta[j][k] = symcub_[j][0]*a[k].x + symcub_[j][1]*a[k].y + symcub_[j][2]*a[k].z;
    }
  }

  for (int j=0; j<3; j++) {
    for (int k=0; k<3; k++) {
      symxtal_[k][j] = 0.0;
      for (int l=0; l<3; l++) {
        symxtal_[k][j] += ainv[j][l]*symmulta[l][k];
      }
    }
  }
}

D3vector SymOp::applyToVector(const D3vector& v, const bool applyfractrans) const {
  D3vector vsym;

  // we use this to identify symmetry-equivalent atoms
  if (applyfractrans) { 
    vsym.x = s11_*v.x + s21_*v.y + s31_*v.z - ftrans1_;
    vsym.y = s12_*v.x + s22_*v.y + s32_*v.z - ftrans2_;
    vsym.z = s13_*v.x + s23_*v.y + s33_*v.z - ftrans3_;
  }
  // we use this to average forces over symmetry equivalent atoms
  else {
    vsym.x = s11_*v.x + s21_*v.y + s31_*v.z;
    vsym.y = s12_*v.x + s22_*v.y + s32_*v.z;
    vsym.z = s13_*v.x + s23_*v.y + s33_*v.z;
  }
  return vsym;
}

void SymOp::applyToTensor(const double* v, double* vsym) {
  // v[0]=v_xx, v[1]=v_yy, v[2]=v_zz, v[3]=v_xy, v[4]=v_yz, v[5]=v_xz, 
  //vsym_ij = sum_kl s(i,k)*v(k,l)*s(j,l) = 
  //  si1_*sj1*v[0] + si1_*sj2*v[3] + si1_*sj3*v[5] + 
  //  si2_*sj1*v[3] + si2_*sj2*v[2] + si2_*sj3*v[4] + 
  //  si3_*sj1*v[5] + si3_*sj2*v[4] + si3_*sj3*v[3];

  vsym[0] = 
    s11_*s11_*v[0] + s11_*s12_*v[3] + s11_*s13_*v[5] + 
    s12_*s11_*v[3] + s12_*s12_*v[1] + s12_*s13_*v[4] + 
    s13_*s11_*v[5] + s13_*s12_*v[4] + s13_*s13_*v[2];
  vsym[1] = 
    s21_*s21_*v[0] + s21_*s22_*v[3] + s21_*s23_*v[5] + 
    s22_*s21_*v[3] + s22_*s22_*v[1] + s22_*s23_*v[4] + 
    s23_*s21_*v[5] + s23_*s22_*v[4] + s23_*s23_*v[2];
  vsym[2] = 
    s31_*s31_*v[0] + s31_*s32_*v[3] + s31_*s33_*v[5] + 
    s32_*s31_*v[3] + s32_*s32_*v[1] + s32_*s33_*v[4] + 
    s33_*s31_*v[5] + s33_*s32_*v[4] + s33_*s33_*v[2];
  vsym[3] = 
    s11_*s21_*v[0] + s11_*s22_*v[3] + s11_*s23_*v[5] + 
    s12_*s21_*v[3] + s12_*s22_*v[1] + s12_*s23_*v[4] + 
    s13_*s21_*v[5] + s13_*s22_*v[4] + s13_*s23_*v[2];
  vsym[4] = 
    s21_*s31_*v[0] + s21_*s32_*v[3] + s21_*s33_*v[5] + 
    s22_*s31_*v[3] + s22_*s32_*v[1] + s22_*s33_*v[4] + 
    s23_*s31_*v[5] + s23_*s32_*v[4] + s23_*s33_*v[2];
  vsym[5] = 
    s11_*s31_*v[0] + s11_*s32_*v[3] + s11_*s33_*v[5] + 
    s12_*s31_*v[3] + s12_*s32_*v[1] + s12_*s33_*v[4] + 
    s13_*s31_*v[5] + s13_*s32_*v[4] + s13_*s33_*v[2];

  return;
}

void SymOp::print(ostream &os) {
  // convertToXtal has to have been called previously
  const double thresh = 1.E-7;
  os.setf(ios::fixed,ios::floatfield);
  os << "symmetry  ";
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      const double stmp = abs(symxtal_[i][j]);
      if (abs(stmp-1.0) < thresh || abs(stmp-0.0) < thresh || abs(stmp-0.25) < thresh || abs(stmp-0.50) < thresh || abs(stmp-0.75) < thresh)
        os << setprecision(2) << symxtal_[i][j] << "  ";
      else
        os << setprecision(10) << symxtal_[i][j] << "  ";
    }
  }
  if (ftrans1_ != 0.0 || ftrans2_ != 0.0 || ftrans3_ != 0.0) 
    os << setprecision(10) << ftrans1_ << "  " << setprecision(10) << ftrans2_ << "  " << setprecision(10) << ftrans3_;

  os << endl;
}

ostream& operator << ( ostream &os, SymOp &sym ) {
  sym.print(os);
  return os;
}
