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
// UnitCell.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <qball/UnitCell.h>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void UnitCell::set(const D3vector& a0, const D3vector& a1, const D3vector& a2)
{
  a_[0] = a0; a_[1] = a1, a_[2] = a2;
  amat_[0] = a0.x;
  amat_[1] = a0.y;
  amat_[2] = a0.z;
  amat_[3] = a1.x;
  amat_[4] = a1.y;
  amat_[5] = a1.z;
  amat_[6] = a2.x;
  amat_[7] = a2.y;
  amat_[8] = a2.z;
  
  // volume = det(A)
  volume_ = a0 * ( a1 ^ a2 );
  if ( volume_ > 0.0 )
  {
    // Compute rows of A-1 (columns of A^-T)
    double fac = 1.0 / volume_;
    D3vector amt0 = fac * a1 ^ a2;
    D3vector amt1 = fac * a2 ^ a0;
    D3vector amt2 = fac * a0 ^ a1;
    
    amat_inv_[0] = amt0.x;
    amat_inv_[1] = amt1.x;
    amat_inv_[2] = amt2.x;
    amat_inv_[3] = amt0.y;
    amat_inv_[4] = amt1.y;
    amat_inv_[5] = amt2.y;
    amat_inv_[6] = amt0.z;
    amat_inv_[7] = amt1.z;
    amat_inv_[8] = amt2.z;
    
    amat_inv_t_[0] = amt0.x;
    amat_inv_t_[1] = amt0.y;
    amat_inv_t_[2] = amt0.z;
    amat_inv_t_[3] = amt1.x;
    amat_inv_t_[4] = amt1.y;
    amat_inv_t_[5] = amt1.z;
    amat_inv_t_[6] = amt2.x;
    amat_inv_t_[7] = amt2.y;
    amat_inv_t_[8] = amt2.z;
    
    // B = 2 pi A^-T
    b_[0] = 2.0 * M_PI * amt0;
    b_[1] = 2.0 * M_PI * amt1;
    b_[2] = 2.0 * M_PI * amt2;

    bmat_[0] = b_[0].x;
    bmat_[1] = b_[0].y;
    bmat_[2] = b_[0].z;
    bmat_[3] = b_[1].x;
    bmat_[4] = b_[1].y;
    bmat_[5] = b_[1].z;
    bmat_[6] = b_[2].x;
    bmat_[7] = b_[2].y;
    bmat_[8] = b_[2].z;
  }
  else
  {
    b_[0] = b_[1] = b_[2] = D3vector(0.0,0.0,0.0);
    amat_inv_[0] =  amat_inv_[1] =  amat_inv_[2] = 
      amat_inv_[3] =  amat_inv_[4] =  amat_inv_[5] = 
      amat_inv_[6] =  amat_inv_[7] =  amat_inv_[8] = 0.0;
    bmat_[0] =  bmat_[1] =  bmat_[2] = 
      bmat_[3] =  bmat_[4] =  bmat_[5] = 
      bmat_[6] =  bmat_[7] =  bmat_[8] = 0.0;
  }
  
  
  
  an_[0]  = a_[0];
  an_[1]  = a_[1];
  an_[2]  = a_[2];
  an_[3]  = a_[0] + a_[1];
  an_[4]  = a_[0] - a_[1];
  an_[5]  = a_[1] + a_[2];
  an_[6]  = a_[1] - a_[2];
  an_[7]  = a_[2] + a_[0];
  an_[8]  = a_[2] - a_[0];
  an_[9]  = a_[0] + a_[1] + a_[2];
  an_[10] = a_[0] - a_[1] - a_[2];
  an_[11] = a_[0] + a_[1] - a_[2];
  an_[12] = a_[0] - a_[1] + a_[2];
 
  bn_[0]  = b_[0];
  bn_[1]  = b_[1];
  bn_[2]  = b_[2];
  bn_[3]  = b_[0] + b_[1];
  bn_[4]  = b_[0] - b_[1];
  bn_[5]  = b_[1] + b_[2];
  bn_[6]  = b_[1] - b_[2];
  bn_[7]  = b_[2] + b_[0];
  bn_[8]  = b_[2] - b_[0];
  bn_[9]  = b_[0] + b_[1] + b_[2];
  bn_[10] = b_[0] - b_[1] - b_[2];
  bn_[11] = b_[0] + b_[1] - b_[2];
  bn_[12] = b_[0] - b_[1] + b_[2];
 
  for ( int i = 0; i < 13; i++ )
  {
    an2h_[i] = 0.5 * norm(an_[i]);
    bn2h_[i] = 0.5 * norm(bn_[i]);
  }
}
 
////////////////////////////////////////////////////////////////////////////////
D3vector UnitCell::cart_to_crystal(const D3vector& v) const {
  D3vector vcryst;
  const double twopiinv = 0.5/M_PI;
  vcryst.x = b_[0]*v*twopiinv;
  vcryst.y = b_[1]*v*twopiinv;
  vcryst.z = b_[2]*v*twopiinv;
  return vcryst;
}
////////////////////////////////////////////////////////////////////////////////
D3vector UnitCell::crystal_to_cart(const D3vector& v) const {
  D3vector vcart = v.x*a_[0] + v.y*a_[1] + v.z*a_[2];
  return vcart;
}
////////////////////////////////////////////////////////////////////////////////
void UnitCell::cart_to_crystal(const double* scart, double* scryst) const {
  // convert symmetric 3x3 matrix scart to crystalline coordinates
  // scart[0] = s_xx, scart[1] = s_yy, scart[2] = s_zz, scart[3] = s_xy, 
  // scart[4] = s_yz, scart[5] = s_xz

  scryst[0] = 
    a_[0].x*scart[0]*a_[0].x + a_[0].x*scart[3]*a_[0].y + a_[0].x*scart[5]*a_[0].z + 
    a_[0].y*scart[3]*a_[0].x + a_[0].y*scart[1]*a_[0].y + a_[0].y*scart[4]*a_[0].z + 
    a_[0].z*scart[5]*a_[0].x + a_[0].z*scart[4]*a_[0].y + a_[0].z*scart[2]*a_[0].z;
  scryst[1] = 
    a_[1].x*scart[0]*a_[1].x + a_[1].x*scart[3]*a_[1].y + a_[1].x*scart[5]*a_[1].z + 
    a_[1].y*scart[3]*a_[1].x + a_[1].y*scart[1]*a_[1].y + a_[1].y*scart[4]*a_[1].z + 
    a_[1].z*scart[5]*a_[1].x + a_[1].z*scart[4]*a_[1].y + a_[1].z*scart[2]*a_[1].z;
  scryst[2] = 
    a_[2].x*scart[0]*a_[2].x + a_[2].x*scart[3]*a_[2].y + a_[2].x*scart[5]*a_[2].z + 
    a_[2].y*scart[3]*a_[2].x + a_[2].y*scart[1]*a_[2].y + a_[2].y*scart[4]*a_[2].z + 
    a_[2].z*scart[5]*a_[2].x + a_[2].z*scart[4]*a_[2].y + a_[2].z*scart[2]*a_[2].z;
  scryst[3] = 
    a_[0].x*scart[0]*a_[1].x + a_[0].x*scart[3]*a_[1].y + a_[0].x*scart[5]*a_[1].z + 
    a_[0].y*scart[3]*a_[1].x + a_[0].y*scart[1]*a_[1].y + a_[0].y*scart[4]*a_[1].z + 
    a_[0].z*scart[5]*a_[1].x + a_[0].z*scart[4]*a_[1].y + a_[0].z*scart[2]*a_[1].z;
  scryst[4] = 
    a_[1].x*scart[0]*a_[2].x + a_[1].x*scart[3]*a_[2].y + a_[1].x*scart[5]*a_[2].z + 
    a_[1].y*scart[3]*a_[2].x + a_[1].y*scart[1]*a_[2].y + a_[1].y*scart[4]*a_[2].z + 
    a_[1].z*scart[5]*a_[2].x + a_[1].z*scart[4]*a_[2].y + a_[1].z*scart[2]*a_[2].z;
  scryst[5] = 
    a_[0].x*scart[0]*a_[2].x + a_[0].x*scart[3]*a_[2].y + a_[0].x*scart[5]*a_[2].z + 
    a_[0].y*scart[3]*a_[2].x + a_[0].y*scart[1]*a_[2].y + a_[0].y*scart[4]*a_[2].z + 
    a_[0].z*scart[5]*a_[2].x + a_[0].z*scart[4]*a_[2].y + a_[0].z*scart[2]*a_[2].z;

  return;
}
////////////////////////////////////////////////////////////////////////////////
void UnitCell::crystal_to_cart(const double* scryst, double* scart) const {
  // convert symmetric 3x3 matrix scryst to Cartesian coordinates
  // scryst[0] = s_xx, scryst[1] = s_yy, scryst[2] = s_zz, scryst[3] = s_xy, 
  // scryst[4] = s_yz, scryst[5] = s_xz
  const double twopiinv = 0.5/M_PI;

  scart[0] = 
    b_[0].x*scryst[0]*b_[0].x + b_[0].x*scryst[3]*b_[1].x + b_[0].x*scryst[5]*b_[2].x + 
    b_[1].x*scryst[3]*b_[0].x + b_[1].x*scryst[1]*b_[1].x + b_[1].x*scryst[4]*b_[2].x + 
    b_[2].x*scryst[5]*b_[0].x + b_[2].x*scryst[4]*b_[1].x + b_[2].x*scryst[2]*b_[2].x;
  scart[1] = 
    b_[0].y*scryst[0]*b_[0].y + b_[0].y*scryst[3]*b_[1].y + b_[0].y*scryst[5]*b_[2].y + 
    b_[1].y*scryst[3]*b_[0].y + b_[1].y*scryst[1]*b_[1].y + b_[1].y*scryst[4]*b_[2].y + 
    b_[2].y*scryst[5]*b_[0].y + b_[2].y*scryst[4]*b_[1].y + b_[2].y*scryst[2]*b_[2].y;
  scart[2] = 
    b_[0].z*scryst[0]*b_[0].z + b_[0].z*scryst[3]*b_[1].z + b_[0].z*scryst[5]*b_[2].z + 
    b_[1].z*scryst[3]*b_[0].z + b_[1].z*scryst[1]*b_[1].z + b_[1].z*scryst[4]*b_[2].z + 
    b_[2].z*scryst[5]*b_[0].z + b_[2].z*scryst[4]*b_[1].z + b_[2].z*scryst[2]*b_[2].z;
  scart[3] = 
    b_[0].x*scryst[0]*b_[0].y + b_[0].x*scryst[3]*b_[1].y + b_[0].x*scryst[5]*b_[2].y + 
    b_[1].x*scryst[3]*b_[0].y + b_[1].x*scryst[1]*b_[1].y + b_[1].x*scryst[4]*b_[2].y + 
    b_[2].x*scryst[5]*b_[0].y + b_[2].x*scryst[4]*b_[1].y + b_[2].x*scryst[2]*b_[2].y;
  scart[4] = 
    b_[0].y*scryst[0]*b_[0].z + b_[0].y*scryst[3]*b_[1].z + b_[0].y*scryst[5]*b_[2].z + 
    b_[1].y*scryst[3]*b_[0].z + b_[1].y*scryst[1]*b_[1].z + b_[1].y*scryst[4]*b_[2].z + 
    b_[2].y*scryst[5]*b_[0].z + b_[2].y*scryst[4]*b_[1].z + b_[2].y*scryst[2]*b_[2].z;
  scart[5] = 
    b_[0].x*scryst[0]*b_[0].z + b_[0].x*scryst[3]*b_[1].z + b_[0].x*scryst[5]*b_[2].z + 
    b_[1].x*scryst[3]*b_[0].z + b_[1].x*scryst[1]*b_[1].z + b_[1].x*scryst[4]*b_[2].z + 
    b_[2].x*scryst[5]*b_[0].z + b_[2].x*scryst[4]*b_[1].z + b_[2].x*scryst[2]*b_[2].z;

  for (int i=0; i<6; i++)
    scart[i] *= twopiinv*twopiinv;

  return;
}
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::in_ws(const D3vector& v) const
{
  bool in = true;
  int i = 0;
  while ( i < 13 && in )
  {
    in = ( abs(v*an_[i]) <= an2h_[i] ) ;
    i++;
  }
  return in;
}

////////////////////////////////////////////////////////////////////////////////
double UnitCell::min_wsdist(void) const {

  double min = sqrt(2.*an2h_[0]);
  for (int i=1; i<13; i++) 
    if (sqrt(2.*an2h_[i]) < min) min = sqrt(2.*an2h_[i]);

  return min;
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::fold_in_ws(D3vector& v) const
{
  const double epsilon = 1.e-10;
  bool done = false;
  const int maxiter = 10;
  int iter = 0;
  while ( !done && iter < maxiter )
  {
    done = true;
    for ( int i = 0; (i < 13) && done; i++ )
    {
      const double sp = v*an_[i];
      if ( sp > an2h_[i] + epsilon )
      {
        done = false;
        do
          v -= an_[i];
        while ( v*an_[i] > an2h_[i] + epsilon );
      }
      else if ( sp < -an2h_[i] - epsilon )
      {
        done = false;
        do
          v += an_[i];
        while ( v*an_[i] < -an2h_[i] - epsilon );
      }
    }
    iter++;
  }
  assert(iter < maxiter);
}
 
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::in_bz(const D3vector& k) const
{
  bool in = true;
  int i = 0;
  while ( i < 13 && in )
  {
    in = ( abs(k*bn_[i]) <= bn2h_[i] ) ;
    i++;
  }
  return in;
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::fold_in_bz(D3vector& k) const
{
  const double epsilon = 1.e-10;
  bool done = false;
  const int maxiter = 10;
  int iter = 0;
  while ( !done && iter < maxiter )
  {
    done = true;
    for ( int i = 0; (i < 13) && done; i++ )
    {
      double sp = k*bn_[i];
      if ( sp > bn2h_[i] + epsilon )
      {
        done = false;
        do
          k -= bn_[i];
        while ( k*bn_[i] > bn2h_[i] + epsilon );
      }
      else if ( sp < -bn2h_[i] - epsilon )
      {
        done = false;
        do
          k += bn_[i];
        while ( k*bn_[i] < -bn2h_[i] - epsilon );
      }
    }
    iter++;
  }
  assert(iter < maxiter);
}
 
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::encloses(const UnitCell& c) const
{
  bool in = true;
  int i = 0;
  while ( i < 13 && in )
  {
    in = ( contains(c.an_[i]) );
    if ( !in )
      cout << "UnitCell::encloses: " << c.an_[i] << " not in cell "
      << c << endl;
    i++;
  }
  return in;
}

////////////////////////////////////////////////////////////////////////////////
bool UnitCell::contains(D3vector v) const
{
  const double fac = 0.5 / ( 2.0 * M_PI );
  const double p0 = fac * v * b_[0];
  const double p1 = fac * v * b_[1];
  const double p2 = fac * v * b_[2];
  return ( (p0 > 0.0) && (p0 <= 1.0) &&
           (p1 > 0.0) && (p1 <= 1.0) &&
           (p2 > 0.0) && (p2 <= 1.0) );
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::print(ostream& os) const
{
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(8);
  os << "<unit_cell " << endl;
  os << "    a=\"" << setw(12) << a_[0].x << " "
     << setw(12) << a_[0].y << " "
     << setw(12) << a_[0].z << "\"" << endl;
  os << "    b=\"" << setw(12) << a_[1].x << " "
     << setw(12) << a_[1].y << " "
     << setw(12) << a_[1].z << "\"" << endl;
  os << "    c=\"" << setw(12) << a_[2].x << " "
     << setw(12) << a_[2].y << " "
     << setw(12) << a_[2].z << "\"" << " />" << endl;
  /*
  os << "    <volume> " << setw(12) << volume_ << " </volume>" << endl;
  os << "    <a0> " << setw(12) << a_[0].x << " " 
                   << setw(12) << a_[0].y << " " 
                   << setw(12) << a_[0].z << " </a0>" << endl;
  os << "    <a1> " << setw(12) << a_[1].x << " " 
                   << setw(12) << a_[1].y << " " 
                   << setw(12) << a_[1].z << " </a1>" << endl;
  os << "    <a2> " << setw(12) << a_[2].x << " " 
                   << setw(12) << a_[2].y << " " 
                   << setw(12) << a_[2].z << " </a2>" << endl;
  os << "    <b0> " << setw(12) << b_[0].x << " " 
                   << setw(12) << b_[0].y << " " 
                   << setw(12) << b_[0].z << " </b0>" << endl;
  os << "    <b1> " << setw(12) << b_[1].x << " " 
                   << setw(12) << b_[1].y << " " 
                   << setw(12) << b_[1].z << " </b1>" << endl;
  os << "    <b2> " << setw(12) << b_[2].x << " " 
                   << setw(12) << b_[2].y << " " 
                   << setw(12) << b_[2].z << " </b2>" << endl;
  */
}
  
////////////////////////////////////////////////////////////////////////////////
void UnitCell::printsys(ostream& os, string setcmd) const {
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(8);
  os << setcmd.c_str()
     << setw(12) << a_[0].x << " " << setw(12) << a_[0].y << " " << setw(12) << a_[0].z << " " 
     << setw(12) << a_[1].x << " " << setw(12) << a_[1].y << " " << setw(12) << a_[1].z << " " 
     << setw(12) << a_[2].x << " " << setw(12) << a_[2].y << " " << setw(12) << a_[2].z
     << " bohr" << endl;
}
  
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::operator==(const UnitCell& c) const
{
  return ( a_[0]==c.a_[0] && a_[1]==c.a_[1] && a_[2]==c.a_[2] );
}
 
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::operator!=(const UnitCell& c) const
{
  return ! ( *this == c );
}

////////////////////////////////////////////////////////////////////////////////
void UnitCell::vecmult3x3(const double* x, const double* y, double *z) const
{
  //  | z0 |     | x0 x3 x6 |     | y0 |
  //  | z1 |  =  | x1 x4 x7 |  *  | y1 |
  //  | z2 |     | x2 x5 x8 |     | y2 |
  
  const double z0 = x[0]*y[0]+x[3]*y[1]+x[6]*y[2];
  const double z1 = x[1]*y[0]+x[4]*y[1]+x[7]*y[2];
  const double z2 = x[2]*y[0]+x[5]*y[1]+x[8]*y[2];
  z[0] = z0;
  z[1] = z1;
  z[2] = z2;
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::vecsmult3x3(const double* xs, const double* y, double *z) const
{
  // multiply a vector by a symmetric 3x3 matrix

  //  | z0 |     | xs0 xs3 xs5 |     | y0 |
  //  | z1 |  =  | xs3 xs1 xs4 |  *  | y1 |
  //  | z2 |     | xs5 xs4 xs2 |     | y2 |

  z[0] = xs[0]*y[0]+xs[3]*y[1]+xs[5]*y[2];
  z[1] = xs[3]*y[0]+xs[1]*y[1]+xs[4]*y[2];
  z[2] = xs[5]*y[0]+xs[4]*y[1]+xs[2]*y[2];

}

////////////////////////////////////////////////////////////////////////////////
void UnitCell::matmult3x3(const double* x, const double* y, double *z) const
{
  //  | z0 z3 z6 |     | x0 x3 x6 |     | y0 y3 y6 |
  //  | z1 z4 z7 |  =  | x1 x4 x7 |  *  | y1 y4 y7 |
  //  | z2 z5 z8 |     | x2 x5 x8 |     | y2 y5 y8 |
  
  const double z00 = x[0]*y[0]+x[3]*y[1]+x[6]*y[2];
  const double z10 = x[1]*y[0]+x[4]*y[1]+x[7]*y[2];
  const double z20 = x[2]*y[0]+x[5]*y[1]+x[8]*y[2];
  
  const double z01 = x[0]*y[3]+x[3]*y[4]+x[6]*y[5];
  const double z11 = x[1]*y[3]+x[4]*y[4]+x[7]*y[5];
  const double z21 = x[2]*y[3]+x[5]*y[4]+x[8]*y[5];
  
  const double z02 = x[0]*y[6]+x[3]*y[7]+x[6]*y[8];
  const double z12 = x[1]*y[6]+x[4]*y[7]+x[7]*y[8];
  const double z22 = x[2]*y[6]+x[5]*y[7]+x[8]*y[8];
  
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
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::smatmult3x3(const double* xs, const double* y, double *z) const
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
////////////////////////////////////////////////////////////////////////////////
void UnitCell::compute_deda(const valarray<double>& sigma,
 valarray<double>& deda) const
{
  // Compute energy derivatives dE/da_ij from a symmetric stress tensor
  assert(sigma.size()==6);
  assert(deda.size()==9);
  
  //!! local copy of sigma to circumvent bug in icc compiler
  valarray<double> sigma_loc(6);
  sigma_loc = sigma;
  
  // deda = - omega * sigma * A^-T
  smatmult3x3(&sigma_loc[0],&amat_inv_t_[0],&deda[0]);
  
  deda *= -volume_;
}
 
////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const UnitCell& cell )
{ 
  cell.print(os); 
  return os;
}
